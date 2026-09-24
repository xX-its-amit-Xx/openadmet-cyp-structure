"""Audit a built submission zip against what the ORGANISERS' scorer actually does.

Written on drop night, 2026-09-24, to settle the §9 "residual risk" of
`docs/DROP_RUN_2026-09-23.md` and — as it turned out — to catch a larger one.

**The organisers' scorer is public**, in `github.com/OpenADMET/blind-challenge-template`
(`backend/evaluate_predictions.py::score_single_structure`, `docs/structure-scoring.md`).
Two facts from reading it drive this script:

1. **Residue numbering does not matter to them.** They run OpenStructure's
   `LDDTPLIScorer`/`SCRMSDScorer` without `resnum_alignments` or `seqres`, and OST's
   default is `resnum_alignments=False`, so `ChainMapper` maps chains and residues by
   **Needleman-Wunsch alignment**, never by residue number. FINDING 021 / trap T1 is a
   hazard of OUR scorer (`cypstruct.pose.lddt_pli` pairs by `(resnum, atom_name)`), not
   of theirs. `numbering` below measures what a numbering error WOULD cost under the
   pessimistic scorer, which is the only honest way to size the risk we were carrying.

2. **`POSEBUSTERS_MAX_FAILURES = 0`.** After OST scoring they run
   `PoseBusters(config="dock")` on the PREDICTED complex, with `mol_cond` built from
   every non-`LIG` fragment combined, and **one** failed check of the 22 overrides
   LDDT-PLI and LDDT-LP to 0.0 and BiSyRMSD to 20 A. Coverage stays 1.0, so a zeroed
   compound does not read as a failure anywhere. `pbgate` replicates that call for call,
   so a hard zero is found here rather than on the board.

    python scripts/submit/submission_audit.py numbering \
        --zip submissions/02_blind_rehearsal_val87b.zip \
        --ligands data/processed/validation_ligands.csv --tag val87
    python scripts/submit/submission_audit.py pbgate \
        --zip submissions/02_drop.zip \
        --ligands data/processed/test_ligands.csv --tag drop
    python scripts/submit/submission_audit.py pbrescue --tag drop \
        --pool C:/cyp_struct/matched_depth/poses/drop_flat \
        --ligands data/processed/test_ligands.csv \
        --xeng data/processed/xeng_drop.csv --ids OCNT-2313038 OCNT-2313429
"""
from __future__ import annotations

import argparse
import json
import sys
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "submit"))

from cypstruct import pose as P  # noqa: E402
from cypstruct.paths import DATA_PROCESSED, REFERENCE, SCRATCH, guard_scratch  # noqa: E402


def _id_col(df: pd.DataFrame) -> str:
    return "structure" if "structure" in df.columns and "id" not in df.columns else "id"


# --------------------------------------------------------------------------
# 1. numbering
# --------------------------------------------------------------------------
def shift(cx: P.Complex, k: int) -> P.Complex:
    """Relabel PROTEIN residue numbers by +k. A pure relabelling.

    Coordinates, atom counts, the ligand, the heme and the axial SG are the same
    objects. Only `prot_key` and `prot_res` change, which is exactly the column a
    renumbering would touch in the PDB.
    """
    return P.Complex(
        name=cx.name, prot_xyz=cx.prot_xyz,
        prot_key=[(c, r + k, a) for (c, r, a) in cx.prot_key],
        prot_res={(c, r + k): v for (c, r), v in cx.prot_res.items()},
        lig_xyz=cx.lig_xyz, lig_elem=list(cx.lig_elem), lig_name=cx.lig_name,
        lig_chain=cx.lig_chain, fe=cx.fe, heme_xyz=cx.heme_xyz,
        heme_atom=list(cx.heme_atom), heme_elem=list(cx.heme_elem),
        axial_sg=cx.axial_sg)


def _load_reference(pdb_id: str, ligand_code: str):
    """The deposited structure, restricted to the ONE chain that holds the ligand."""
    cif = REFERENCE / "rcsb" / f"{pdb_id}.cif"
    if not cif.exists():
        return None, "no_cif"
    import gemmi
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    for chain in st[0]:
        if any(r.name.strip().upper() == ligand_code.upper() for r in chain):
            try:
                return P.load_structure(cif, ligand_code=ligand_code,
                                        assembly_chain=chain.name), "ok"
            except Exception as exc:                     # noqa: BLE001
                return None, f"load:{type(exc).__name__}"
    return None, "ligand_not_in_any_chain"


def numbering(zip_path: Path, ligands_csv: Path, tag: str,
              offsets: list[int]) -> None:
    lig = pd.read_csv(ligands_csv)
    idc = _id_col(lig)
    lig[idc] = lig[idc].astype(str)
    if "pdb" not in lig.columns:
        raise SystemExit(f"{ligands_csv} has no `pdb` column - there is nothing to score "
                         "against. This arm only runs on a set with crystals.")
    z = zipfile.ZipFile(zip_path)
    names = {Path(n).stem: n for n in z.namelist() if n.endswith(".pdb")}
    print(f"{zip_path.name}: {len(names)} pdbs; ligand table {len(lig)} rows", flush=True)

    guard_scratch(0.5)
    tmp = SCRATCH / f"audit_{tag}"
    tmp.mkdir(parents=True, exist_ok=True)
    rows, skipped = [], {}
    for r in lig.itertuples():
        sid = str(getattr(r, idc))
        n = names.get(sid)
        if n is None:
            skipped[sid] = "not_in_zip"
            continue
        p = tmp / f"{sid}.pdb"
        p.write_bytes(z.read(n))
        try:
            model = P.load_structure(p, ligand_code="LIG")
        except Exception as exc:                          # noqa: BLE001
            skipped[sid] = f"model:{type(exc).__name__}"
            continue
        ref, why = _load_reference(r.pdb, sid)
        if ref is None or len(ref.lig_xyz) == 0:
            skipped[sid] = why if ref is None else "ref_no_ligand"
            continue
        try:
            perm = P.best_ligand_mapping(r.smiles, model, ref)
        except Exception:                                 # noqa: BLE001
            perm = None
        row = {"id": sid, "pdb": r.pdb}
        for k in offsets:
            m = model if k == 0 else shift(model, k)
            row[f"shared_{k}"] = len({(rn, a) for _c, rn, a in m.prot_key}
                                     & {(rn, a) for _c, rn, a in ref.prot_key})
            row[f"lddt_{k}"] = P.lddt_pli(m, ref, lig_perm=perm)
        rows.append(row)
        p.unlink(missing_ok=True)

    df = pd.DataFrame(rows)
    df.to_csv(DATA_PROCESSED / f"numbering_penalty_{tag}.csv", index=False)
    out = {"zip": str(zip_path), "n_scored": len(df), "n_skipped": len(skipped),
           "skipped": skipped, "offsets": offsets, "arms": {}}
    print(f"\nscored {len(df)} / {len(lig)}; skipped {len(skipped)}"
          f"{' ' + str(dict(Counter(skipped.values()))) if skipped else ''}")
    for k in offsets:
        c = df[f"lddt_{k}"]
        a = {"mean": float(np.nanmean(c)), "median": float(np.nanmedian(c)),
             "n_nan": int(c.isna().sum()),
             "n_exactly_zero": int((c.fillna(-1) == 0.0).sum()),
             "n_below_0.1": int((c.fillna(-1) < 0.1).sum()),
             "mean_shared_protein_atoms": float(df[f"shared_{k}"].mean())}
        out["arms"][str(k)] = a
        print(f"offset +{k:>2}: mean {a['mean']:.4f}  median {a['median']:.4f}  "
              f"exactly 0.0000 {a['n_exactly_zero']:>3}/{len(df)}  "
              f"<0.1 {a['n_below_0.1']:>3}  NaN {a['n_nan']}  "
              f"shared prot atoms {a['mean_shared_protein_atoms']:.0f}")
    base = out["arms"][str(offsets[0])]["mean"]
    for k in offsets[1:]:
        d = out["arms"][str(k)]["mean"] - base
        out["arms"][str(k)]["delta_vs_first_arm"] = d
        print(f"penalty of a +{k} numbering error: {d:+.4f}")
    (DATA_PROCESSED / f"numbering_penalty_{tag}.json").write_text(json.dumps(out, indent=1))


# --------------------------------------------------------------------------
# 2. the organisers' PoseBusters gate
# --------------------------------------------------------------------------
def _pb_one(pb, pdb_path: Path, smiles: str | None) -> dict:
    """One compound through the backend's exact PoseBusters call."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    pred = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
    if pred is None:
        return {"status": "rdkit_unparseable__backend_KEEPS_ost_scores"}
    frags = Chem.SplitMolByPDBResidues(pred)
    ligmol = frags.pop("LIG", None)
    if ligmol is None:
        return {"status": "no_LIG_fragment__backend_KEEPS_ost_scores"}
    prot = None
    for f in frags.values():
        prot = f if prot is None else Chem.CombineMols(prot, f)
    tmpl = Chem.MolFromSmiles(smiles) if smiles else None
    rec: dict = {}
    if tmpl is not None:
        try:
            ligmol = AllChem.AssignBondOrdersFromTemplate(tmpl, ligmol)
        except Exception:                                 # noqa: BLE001
            rec["bond_orders"] = "FAILED_single_bonds"
            tmpl = None
    try:
        res = pb.bust(mol_pred=ligmol, mol_true=tmpl, mol_cond=prot)
    except Exception as exc:                              # noqa: BLE001
        rec["status"] = f"posebusters_raised_{type(exc).__name__}__NaN_coverage0"
        return rec
    bcols = [c for c in res.columns if res[c].dtype == bool]
    row = res[bcols].iloc[0]
    failed = sorted(row.index[~row].tolist())
    rec.update({"status": "ZEROED" if failed else "pass",
                "n_checks": len(bcols), "n_failed": len(failed),
                "failed": ";".join(failed)})
    return rec


def pbgate(zip_path: Path, ligands_csv: Path, tag: str) -> None:
    from posebusters import PoseBusters
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")

    lig = pd.read_csv(ligands_csv)
    idc = _id_col(lig)
    smi = dict(zip(lig[idc].astype(str), lig["smiles"].astype(str)))
    z = zipfile.ZipFile(zip_path)
    names = sorted(n for n in z.namelist() if n.endswith(".pdb"))
    guard_scratch(0.5)
    tmp = SCRATCH / f"pbgate_{tag}"
    tmp.mkdir(parents=True, exist_ok=True)
    pb = PoseBusters(config="dock")
    rows, fails = [], Counter()
    for i, n in enumerate(names, 1):
        sid = Path(n).stem
        p = tmp / f"{sid}.pdb"
        p.write_bytes(z.read(n))
        rec = {"id": sid}
        rec.update(_pb_one(pb, p, smi.get(sid)))
        for f in filter(None, rec.get("failed", "").split(";")):
            fails[f] += 1
        rows.append(rec)
        print(f"  {i}/{len(names)} {sid}: {rec['status']} "
              f"({rec.get('n_failed', '-')} failed) {rec.get('failed', '')[:90]}", flush=True)
        p.unlink(missing_ok=True)

    df = pd.DataFrame(rows)
    df.to_csv(DATA_PROCESSED / f"posebusters_gate_{tag}.csv", index=False)
    n_zero = int((df.status == "ZEROED").sum())
    n_pass = int((df.status == "pass").sum())
    print(f"\n{tag}: {len(df)} files | ZEROED {n_zero} | pass {n_pass} | "
          f"other {len(df) - n_zero - n_pass}")
    print("failed-check frequency:", dict(fails.most_common()))
    (DATA_PROCESSED / f"posebusters_gate_{tag}.json").write_text(json.dumps(
        {"zip": str(zip_path), "n": len(df), "n_zeroed": n_zero, "n_pass": n_pass,
         "max_pb_failures_upstream": 0,
         "check_frequency": dict(fails.most_common()), "per_file": rows}, indent=1))


# --------------------------------------------------------------------------
# 3. is a passing pose already in the pool?
# --------------------------------------------------------------------------
def pbrescue(pool: Path, ligands_csv: Path, xeng_csv: Path, tag: str,
             ids: list[str]) -> None:
    """For each named ligand, walk its pool poses in xeng order and gate each one.

    Measurement only - it builds no zip and selects nothing. The point is to say
    whether a hard zero is recoverable from the pool already paid for, and at what
    cost in selector rank.
    """
    from posebusters import PoseBusters
    from rdkit import RDLogger

    from build_submission import to_submission_pdb
    RDLogger.DisableLog("rdApp.*")

    lig = pd.read_csv(ligands_csv)
    idc = _id_col(lig)
    smi = dict(zip(lig[idc].astype(str), lig["smiles"].astype(str)))
    xe = pd.read_csv(xeng_csv)
    guard_scratch(0.5)
    tmp = SCRATCH / f"pbrescue_{tag}"
    tmp.mkdir(parents=True, exist_ok=True)
    pb = PoseBusters(config="dock")
    out: dict = {}
    for sid in ids:
        sub = xe[xe.ligand.astype(str) == sid].sort_values("xeng")   # rank 1 == shipped
        recs = []
        for rank, r in enumerate(sub.itertuples(), 1):
            cif = pool / f"{r.sample}.cif"
            if not cif.exists():
                recs.append({"rank": rank, "sample": r.sample, "status": "cif_missing"})
                continue
            p = tmp / f"{sid}_{rank}.pdb"
            to_submission_pdb(cif, p, heme="keep")
            rec = {"rank": rank, "sample": r.sample, "xeng": float(r.xeng)}
            rec.update(_pb_one(pb, p, smi.get(sid)))
            recs.append(rec)
            p.unlink(missing_ok=True)
        n_pass = sum(1 for x in recs if x.get("status") == "pass")
        first = next((x for x in recs if x.get("status") == "pass"), None)
        out[sid] = {"n_poses": len(recs), "n_pass": n_pass,
                    "shipped_rank1_status": recs[0].get("status") if recs else None,
                    "best_passing_rank": first["rank"] if first else None,
                    "best_passing_sample": first["sample"] if first else None,
                    "records": recs}
        print(f"{sid}: {n_pass}/{len(recs)} poses pass; shipped(rank1)="
              f"{out[sid]['shipped_rank1_status']}; first passing at xeng rank "
              f"{out[sid]['best_passing_rank']}", flush=True)
    (DATA_PROCESSED / f"posebusters_pool_rescue_{tag}.json").write_text(
        json.dumps(out, indent=1))


# --------------------------------------------------------------------------
# 4. does gating the SELECTOR on PoseBusters beat the shipped selector?
# --------------------------------------------------------------------------
def _pose_cif(pool: Path, lig: str, sample: str, arm: str, flat: bool):
    cands = [pool / f"{lig}__{arm}__s1" / f"{sample}.cif",
             pool / f"{sample}.cif",
             pool / lig / f"{sample}.cif",
             pool / f"{lig}__{sample}.cif"]
    if flat:
        cands = [c for c in cands if c.parent == pool]
    return next((c for c in cands if c.exists()), None)


def pbselect(pool: Path, ligands_csv: Path, xeng_csv: Path, scored_csv: Path,
             tag: str, arm: str, pool_flat: bool, boot: int, max_walk: int,
             recheck: int) -> None:
    """Rule A (argmin xeng) against rule B (argmin xeng among PB-passing poses).

    Scored the way the ORGANISERS score: a pose that fails one PoseBusters `dock`
    check contributes **0.0**, not its true LDDT-PLI. `lddt_pli` for every pose comes
    from `poses_scored_<tag>.csv`, i.e. from crystals - legitimate here because this is
    the *measurement*, not the rule. Both rules themselves are prediction-side.

    Poses are gated lazily in xeng order and the walk stops at the first pass, so a
    ligand whose shipped pose already passes costs one gate call. A ligand that never
    passes is walked to the end, which is how `no_passing_pose` is counted.
    """
    from posebusters import PoseBusters
    from rdkit import RDLogger

    from build_submission import to_submission_pdb
    RDLogger.DisableLog("rdApp.*")

    lig = pd.read_csv(ligands_csv)
    idc = _id_col(lig)
    smi = dict(zip(lig[idc].astype(str), lig["smiles"].astype(str)))
    xe = pd.read_csv(xeng_csv).dropna(subset=["xeng"])
    xe["ligand"] = xe.ligand.astype(str)
    sc = pd.read_csv(scored_csv)
    sc = sc[sc.arm == arm]
    truth = {(str(r.ligand), str(r.sample)): float(r.lddt_pli) for r in sc.itertuples()}

    guard_scratch(0.5)
    tmp = SCRATCH / f"pbselect_{tag}"
    tmp.mkdir(parents=True, exist_ok=True)
    cache_path = DATA_PROCESSED / f"posebusters_pose_gate_{tag}.csv"
    cache: dict = {}
    if cache_path.exists():
        for r in pd.read_csv(cache_path).itertuples():
            cache[(str(r.ligand), str(r.sample))] = {
                "status": str(r.status),
                "failed": "" if pd.isna(r.failed) else str(r.failed),
                "n_checks": int(r.n_checks) if not pd.isna(r.n_checks) else -1}
        print(f"resuming: {len(cache)} poses already gated", flush=True)
    pb = PoseBusters(config="dock")

    def gate(ligid: str, sample: str) -> dict:
        key = (ligid, sample)
        if key in cache:
            return cache[key]
        cif = _pose_cif(pool, ligid, sample, arm, pool_flat)
        if cif is None:
            rec = {"status": "cif_missing", "failed": "", "n_checks": -1}
        else:
            p = tmp / f"{ligid}_{sample}.pdb"
            to_submission_pdb(cif, p, heme="keep")
            r = _pb_one(pb, p, smi.get(ligid))
            rec = {"status": r.get("status", "?"), "failed": r.get("failed", ""),
                   "n_checks": int(r.get("n_checks", -1))}
            p.unlink(missing_ok=True)
        cache[key] = rec
        return rec

    rows, fails = [], Counter()
    ligs = sorted(xe.ligand.unique())
    for i, ligid in enumerate(ligs, 1):
        sub = xe[xe.ligand == ligid].sort_values("xeng").reset_index(drop=True)
        n_poses = len(sub)
        walk = min(max_walk, n_poses)
        rank1 = str(sub.loc[0, "sample"])
        first_pass, first_rank, n_gated = None, None, 0
        for rank in range(walk):
            s = str(sub.loc[rank, "sample"])
            rec = gate(ligid, s)
            n_gated += 1
            for f in filter(None, rec["failed"].split(";")):
                fails[f] += 1
            if rec["status"] == "pass":
                first_pass, first_rank = s, rank + 1
                break
        r1 = cache[(ligid, rank1)]
        pick_b = first_pass if first_pass is not None else rank1
        t1 = truth.get((ligid, rank1), float("nan"))
        tb = truth.get((ligid, pick_b), float("nan"))
        rows.append({
            "ligand": ligid, "n_poses": n_poses, "n_gated": n_gated,
            "rank1_sample": rank1, "rank1_status": r1["status"],
            "rank1_failed": r1["failed"],
            "b_sample": pick_b, "b_rank": first_rank if first_rank else 1,
            "b_has_passing_pose": first_pass is not None,
            "changed": pick_b != rank1,
            "true_A": t1, "true_B": tb,
            "gated_A": t1 if r1["status"] == "pass" else 0.0,
            "gated_B": tb if first_pass is not None else 0.0,
        })
        print(f"  {i}/{len(ligs)} {ligid}: rank1={r1['status']} gated={n_gated} "
              f"B_rank={first_rank} changed={pick_b != rank1}", flush=True)
        pd.DataFrame([{"ligand": k[0], "sample": k[1], **v}
                      for k, v in cache.items()]).to_csv(cache_path, index=False)

    df = pd.DataFrame(rows)
    df.to_csv(DATA_PROCESSED / f"pbselect_{tag}.csv", index=False)

    a, b = df.gated_A.to_numpy(float), df.gated_B.to_numpy(float)
    d = b - a
    rng = np.random.default_rng(0)
    idx = rng.integers(0, len(d), size=(boot, len(d)))
    bs = d[idx].mean(axis=1)
    lo, hi = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
    p = stat = float("nan")
    try:
        from scipy.stats import wilcoxon
        nz = d[d != 0]
        if len(nz):
            w = wilcoxon(nz)
            p, stat = float(w.pvalue), float(w.statistic)
    except Exception as exc:                              # noqa: BLE001
        print("wilcoxon failed:", exc)

    changed = df[df.changed]
    no_pass = df[~df.b_has_passing_pose]
    out = {
        "tag": tag, "arm": arm, "n_ligands": int(len(df)),
        "pool": str(pool), "xeng": str(xeng_csv), "scored": str(scored_csv),
        "max_walk": max_walk, "n_poses_gated": int(df.n_gated.sum()),
        "board_A_gated": float(a.mean()), "board_B_gated": float(b.mean()),
        "board_A_ungated": float(df.true_A.mean()),
        "board_B_ungated": float(df.true_B.mean()),
        "delta_gated": float(d.mean()),
        "delta_gated_ci95": [lo, hi],
        "delta_ungated": float((df.true_B - df.true_A).mean()),
        "wilcoxon_stat": stat, "wilcoxon_p": p,
        "n_nonzero_pairs": int((d != 0).sum()),
        "n_changed_pick": int(len(changed)),
        "n_rank1_zeroed": int((df.rank1_status != "pass").sum()),
        "n_no_passing_pose": int(len(no_pass)),
        "no_passing_pose": sorted(no_pass.ligand.tolist()),
        "changed": [{"ligand": r.ligand, "b_rank": int(r.b_rank),
                     "rank1_status": r.rank1_status,
                     "true_A": r.true_A, "true_B": r.true_B,
                     "gated_A": r.gated_A, "gated_B": r.gated_B}
                    for r in changed.itertuples()],
        "check_frequency_over_gated_poses": dict(fails.most_common()),
        "b_rank_histogram": {str(k): int(v) for k, v in
                             sorted(Counter(df.b_rank.tolist()).items())},
    }

    # determinism: re-gate a spread of already-gated files and compare verdicts
    if recheck:
        keys = [k for k in cache if cache[k]["status"] in ("pass", "ZEROED")]
        pick = keys if len(keys) <= recheck else [
            keys[i] for i in np.linspace(0, len(keys) - 1, recheck).astype(int)]
        same = 0
        for ligid, s in pick:
            before = cache[(ligid, s)]
            cif = _pose_cif(pool, ligid, s, arm, pool_flat)
            pp = tmp / f"recheck_{ligid}_{s}.pdb"
            to_submission_pdb(cif, pp, heme="keep")
            r = _pb_one(pb, pp, smi.get(ligid))
            pp.unlink(missing_ok=True)
            same += int(r.get("status") == before["status"]
                        and r.get("failed", "") == before["failed"])
        out["determinism_rechecked"] = len(pick)
        out["determinism_identical"] = same
        print(f"\ndeterminism: {same}/{len(pick)} re-runs identical")

    print(f"\n{tag}: {len(df)} ligands, {out['n_poses_gated']} poses gated")
    print(f"  A (shipped argmin(xeng))       gated {out['board_A_gated']:.4f}   "
          f"ungated {out['board_A_ungated']:.4f}")
    print(f"  B (argmin(xeng) | PB-passing)  gated {out['board_B_gated']:.4f}   "
          f"ungated {out['board_B_ungated']:.4f}")
    print(f"  paired delta (gated) {out['delta_gated']:+.4f}  "
          f"95% CI [{lo:+.4f}, {hi:+.4f}]  wilcoxon p={p:.3g} "
          f"(n non-zero {out['n_nonzero_pairs']})")
    print(f"  delta (ungated, cost of moving off rank 1) {out['delta_ungated']:+.4f}")
    print(f"  ligands changing pick {out['n_changed_pick']}; "
          f"rank-1 zeroed {out['n_rank1_zeroed']}; "
          f"no passing pose at all {out['n_no_passing_pose']}")
    print("  checks fired over gated poses:", out["check_frequency_over_gated_poses"])
    print("  B rank histogram:", out["b_rank_histogram"])
    (DATA_PROCESSED / f"pbselect_{tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["numbering", "pbgate", "pbrescue", "pbselect"])
    ap.add_argument("--zip", default=None)
    ap.add_argument("--ligands", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--offsets", type=int, nargs="+", default=[0, 24, 27, 43],
                    help="the first is the as-shipped control and the baseline")
    ap.add_argument("--pool", default=None)
    ap.add_argument("--xeng", default=None)
    ap.add_argument("--ids", nargs="+", default=[])
    ap.add_argument("--scored", default=None)
    ap.add_argument("--arm", default="unsteered")
    ap.add_argument("--pool-flat", action="store_true")
    ap.add_argument("--boot", type=int, default=10000)
    ap.add_argument("--max-walk", type=int, default=20)
    ap.add_argument("--recheck", type=int, default=20)
    a = ap.parse_args()
    if a.cmd == "numbering":
        numbering(Path(a.zip), Path(a.ligands), a.tag, a.offsets)
    elif a.cmd == "pbgate":
        pbgate(Path(a.zip), Path(a.ligands), a.tag)
    elif a.cmd == "pbrescue":
        if not (a.pool and a.xeng and a.ids):
            raise SystemExit("pbrescue needs --pool, --xeng and --ids")
        pbrescue(Path(a.pool), Path(a.ligands), Path(a.xeng), a.tag, a.ids)
    else:
        if not (a.pool and a.xeng and a.scored):
            raise SystemExit("pbselect needs --pool, --xeng and --scored")
        pbselect(Path(a.pool), Path(a.ligands), Path(a.xeng), Path(a.scored),
                 a.tag, a.arm, a.pool_flat, a.boot, a.max_walk, a.recheck)
