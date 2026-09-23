"""Fragment pose transfer: where does a functional group sit, relative to the iron?

Pre-registered in `docs/PREREG_fragment_transfer.md`, committed at `fa6257f` before a
single number existed. Nothing here may be tuned after a score exists. Written up as
FINDING 030.

**The lever.** FINDING 024 refuted every intervention that conditions the PROTEIN on
CYP3A4 - the pocket is already right to 0.73 A and rho(protein error, ligand error) =
+0.03. This conditions the LIGAND instead. A functional group that a P450 binds sits
somewhere specific relative to the iron and the porphyrin plane; that placement is an
empirical crystallographic fact over 185 distinct targets, and the heme frame is the one
landmark that exists in all of them.

Four CPU-only stages, each resumable:

    python scripts/structure/fragment_transfer.py frames    # 1,740 poses -> heme frame
    python scripts/structure/fragment_transfer.py donors    # 1,002 crystal observations
    python scripts/structure/fragment_transfer.py match --workers 12
    python scripts/structure/fragment_transfer.py evaluate

**The trap this script is built around** is leakage, because the feature reads real
crystal structures. Five filters (L1-L5) run per query and every one of them reports a
non-zero count, because a filter that never fires is a bug rather than a pass. The
headline excludes the whole CYP3A subfamily - leave-one-TARGET-out, which 185 targets is
what makes possible - and the with-3A variant is computed only as the optimistic bound it
is.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

DATA = REPO / "data" / "processed"
UNIV = DATA / "p450_universe"
POOL = Path(os.environ.get("FT_POOL", "D:/cyp_scratch/val87b_unsteered"))
SCRATCH = Path(os.environ.get(
    "FT_SCRATCH",
    r"C:\tb\tmp\1\claude\D--Users-ashenoy00000--windsurf-OpenADMET-cyp-structure"
    r"\bd288271-2aa5-4ed4-a343-ea31e5ad8c13\scratchpad\frag_transfer"))

# ---- pre-registered constants; none of these may move ---------------------
MIN_MCS_ATOMS = 6                 # primary arm
MIN_MCS_ATOMS_LOOSE = 5           # sensitivity arm
TANIMOTO_CUT = 0.90               # L3, fixed a priori
FE_CUT = 10.0                     # L5, A
MIN_DONORS = 3                    # coverage bar for the primary subset
MAX_MATCHES = 64                  # per side
MAX_PAIRINGS = 4096
TIE_DRAWS = 64
NULL_DRAWS = 2000
BOOT_DRAWS = 10000
SCRAMBLE_REDRAWS = 8
SEED = 20260922

CYP3A_UNIPROT = {"P08684", "P20815", "P24462", "Q9HB55"}


# ==========================================================================
# stage 1 - query poses into their own heme frames
# ==========================================================================

def stage_frames() -> None:
    """1,740 Boltz-2 poses -> ligand heavy atoms in each pose's OWN heme frame.

    `data/processed/xeng_cache/val87b` looks like it already holds this and does not:
    it stores RAW model coordinates. Reading it as a heme frame would put ligand atoms
    6 A below the porphyrin, on the proximal face, where FINDING 008 measured 0.7% of
    ligand atoms in 116 crystals. Checked, not assumed.
    """
    import pandas as pd

    from cypstruct import pose as P
    from cypstruct import xengine as X

    out = SCRATCH / "frames"
    out.mkdir(parents=True, exist_ok=True)
    val = pd.read_csv(DATA / "validation_ligands.csv")
    done = skipped = 0
    for lig in val["id"]:
        f = out / f"{lig}.npz"
        if f.exists():
            skipped += 1
            continue
        job = POOL / f"{lig}__unsteered__s1"
        cifs = sorted(job.glob("input_model_*.cif"),
                      key=lambda p: int(p.stem.rsplit("_", 1)[1]))
        names, frames, raw, elems = [], [], [], None
        for cif in cifs:
            try:
                m = P.load_structure(cif)
            except Exception:
                continue
            v = X.in_heme_frame(m)
            if v is None:
                continue
            names.append(cif.stem)
            frames.append(v)
            raw.append(np.asarray(m.lig_xyz, float))
            elems = [e.upper() for e in m.lig_elem]
        if not names:
            continue
        np.savez_compressed(f, names=np.array(names),
                            xyz=np.array(frames, np.float32),
                            raw=np.array(raw, np.float32),
                            elems=np.array(elems))
        done += 1
    print(f"frames: {done} written, {skipped} cached", flush=True)


def check_frames() -> dict:
    """C4 - the frame computed here reproduces the SHIPPED `xeng` to 1e-6.

    This is FINDING 027's control C3. If the frame were wrong every number downstream
    would be wrong in a way that no amount of statistics would catch, and the PXR
    campaign lost a whole metric to exactly that (CLAUDE.md, lesson 1).
    """
    import pandas as pd

    from cypstruct import xengine as X

    ref = X.load_reference(DATA / "reference_set_cyp3a4.npz")
    shipped = pd.read_csv(DATA / "xeng_val87b.csv")
    shipped = {(r.ligand, r.sample): r.xeng for r in shipped.itertuples()}
    worst, n = 0.0, 0
    for f in sorted((SCRATCH / "frames").glob("*.npz")):
        lig = f.stem
        refs = ref.get(lig, [])
        if not refs:
            continue
        d = np.load(f, allow_pickle=True)
        for name, v in zip(d["names"], d["xyz"]):
            got = shipped.get((lig, str(name)))
            if got is None:
                continue
            mine = X.xeng_score(np.asarray(v, float), refs)
            worst = max(worst, abs(mine - got))
            n += 1
    return {"rows_checked": n, "max_abs_delta": worst, "pass": bool(worst < 1e-6)}


# ==========================================================================
# stage 2 - donor crystals into their own heme frames
# ==========================================================================

def _mol_from_3d(xyz: np.ndarray, elems: list[str], smiles: str | None):
    """RDKit mol in COORDINATE ORDER, so an MCS atom index is an xyz row index.

    Connectivity is perceived from the geometry, then bond orders are transferred from
    the CCD SMILES. Both halves are needed: perception alone cannot tell an aromatic
    ring from a saturated one, and the SMILES alone carries no atom correspondence to
    the deposited coordinates.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDetermineBonds

    if len(xyz) < 3:
        return None, "too-small"
    mol = Chem.RWMol()
    conf = Chem.Conformer(len(xyz))
    for i, (e, p) in enumerate(zip(elems, xyz)):
        sym = e.capitalize() if len(e) > 1 else e.upper()
        try:
            a = Chem.Atom(sym)
        except Exception:
            return None, "bad-element"
        mol.AddAtom(a)
        conf.SetAtomPosition(i, [float(p[0]), float(p[1]), float(p[2])])
    m = mol.GetMol()
    m.AddConformer(conf)
    try:
        rdDetermineBonds.DetermineConnectivity(m, useHueckel=False, charge=0)
    except Exception:
        return None, "connectivity-failed"
    if smiles is None:
        return None, "no-smiles"
    tmpl = Chem.MolFromSmiles(smiles)
    if tmpl is None:
        return None, "bad-smiles"
    tmpl = Chem.RemoveHs(tmpl)
    if tmpl.GetNumAtoms() != m.GetNumAtoms():
        return None, "atom-count-mismatch"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = AllChem.AssignBondOrdersFromTemplate(tmpl, m)
    except Exception:
        return None, "template-failed"
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None, "sanitise-failed"
    return out, "ok"


def stage_donors() -> None:
    """Every crystal ligand observation, in its own heme frame, with a chemistry graph."""
    import pandas as pd
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")

    atoms = pd.read_parquet(UNIV / "p450_atoms.parquet")
    geom = pd.read_parquet(UNIV / "p450_geometry.parquet")
    cof = pd.read_csv(UNIV / "p450_cofold_set.csv")

    smi = dict(zip(cof["id"], cof["smiles"]))
    uni = dict(zip(cof["pdb"], cof["uniprot"]))
    tkey = dict(zip(cof["pdb"], cof["target_key"]))
    fe_d = {(r.pdb, r.chain, r.lig, r.seqid): float(r.closest_fe)
            for r in geom.itertuples()}

    rows, reasons = [], {}
    for r in atoms.itertuples():
        xyz = np.frombuffer(r.xyz, dtype=np.float32).reshape(-1, 3).astype(float)
        elems = [e for e in r.elems.split(",") if e]
        fe = np.frombuffer(r.fe, dtype=np.float32).astype(float)
        n = np.frombuffer(r.normal, dtype=np.float32).astype(float)
        x = np.frombuffer(r.xaxis, dtype=np.float32).astype(float)
        if len(elems) != len(xyz):
            reasons["elem-xyz-mismatch"] = reasons.get("elem-xyz-mismatch", 0) + 1
            continue
        y = np.cross(n, x)
        d = xyz - fe
        frame = np.column_stack([d @ x, d @ y, d @ n])
        mol, why = _mol_from_3d(xyz, elems, smi.get(r.lig))
        reasons[why] = reasons.get(why, 0) + 1
        if mol is None:
            continue
        rows.append({
            "pdb": r.pdb, "chain": r.chain, "lig": r.lig, "seqid": int(r.seqid),
            "uniprot": uni.get(r.pdb), "target_key": tkey.get(r.pdb),
            "closest_fe": fe_d.get((r.pdb, r.chain, r.lig, int(r.seqid)), np.nan),
            "n_heavy": len(elems),
            "smiles": smi.get(r.lig),
            "molblock": Chem.MolToMolBlock(mol),
            "frame": frame.astype(np.float32).tobytes(),
        })
    df = pd.DataFrame(rows)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    df.to_parquet(SCRATCH / "donors.parquet")
    print(f"donors: {len(df)} of {len(atoms)} observations usable", flush=True)
    print("  reasons:", json.dumps(reasons, sort_keys=True), flush=True)
    (SCRATCH / "donor_build.json").write_text(json.dumps(
        {"observations": int(len(atoms)), "usable": int(len(df)),
         "reasons": reasons,
         "unique_codes": int(df["lig"].nunique()),
         "unique_entries": int(df["pdb"].nunique()),
         "unique_targets": int(df["target_key"].nunique())}, indent=2))


# ==========================================================================
# stage 3 - substructure matching and the transfer distance
# ==========================================================================

def _fp(mol):
    from rdkit.Chem import rdFingerprintGenerator
    g = rdFingerprintGenerator.GetMorganGenerator(radius=2)
    return g.GetCountFingerprint(mol)


def _mcs_smarts(a, b, loose: bool):
    from rdkit.Chem import rdFMCS
    p = rdFMCS.MCSParameters()
    p.AtomTyper = rdFMCS.AtomCompare.CompareElements
    p.BondTyper = (rdFMCS.BondCompare.CompareAny if loose
                   else rdFMCS.BondCompare.CompareOrderExact)
    p.BondCompareParameters.RingMatchesRingOnly = not loose
    p.BondCompareParameters.CompleteRingsOnly = not loose
    p.AtomCompareParameters.RingMatchesRingOnly = not loose
    p.Timeout = 10
    try:
        r = rdFMCS.FindMCS([a, b], p)
    except Exception:
        return None, 0
    if r.canceled and r.numAtoms == 0:
        return None, 0
    return r.smartsString, int(r.numAtoms)


def _matches(mol, patt):
    ms = mol.GetSubstructMatches(patt, uniquify=False, maxMatches=MAX_MATCHES,
                                 useChirality=False)
    return [np.asarray(m, int) for m in ms]


def _transfer_distance(pose_frames, q_idx, donor_frame, d_idx):
    """RMSD between the shared fragment of a pose and of a donor, in their own heme
    frames, with NO superposition - the frames are already common. Minimised over
    symmetry-equivalent atom correspondences only; no truth enters the minimisation.
    """
    nq, nd = len(q_idx), len(d_idx)
    while nq * nd > MAX_PAIRINGS:
        if nq >= nd:
            nq -= 1
        else:
            nd -= 1
    A = pose_frames[:, np.stack(q_idx[:nq]), :]          # (P, nq, k, 3)
    B = donor_frame[np.stack(d_idx[:nd]), :]             # (nd, k, 3)
    sq = ((A[:, :, None, :, :] - B[None, None, :, :, :]) ** 2).sum(-1).mean(-1)
    return np.sqrt(sq.reshape(len(pose_frames), -1).min(axis=1)), nq * nd


def _query_record(lig: str, smiles: str):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    f = SCRATCH / "frames" / f"{lig}.npz"
    if not f.exists():
        return None
    d = np.load(f, allow_pickle=True)
    raw = np.asarray(d["raw"], float)
    elems = [str(e) for e in d["elems"]]
    mol, why = _mol_from_3d(raw[0], elems, smiles)
    if mol is None:
        return {"ligand": lig, "status": f"query-mol-{why}"}
    return {"ligand": lig, "status": "ok", "mol": mol,
            "frames": np.asarray(d["xyz"], float),
            "names": [str(x) for x in d["names"]]}


def _crystal_frame(lig: str, pdb: str, query_mol):
    """The query's OWN crystal ligand, in its own heme frame, in POSE atom order.

    Used only as a diagnostic (D5), never as a feature: it is the answer. Atom order
    is transferred by substructure match rather than by array position, because
    index-for-index is the trap that once halved every LDDT-PLI in this project.
    """
    import pandas as pd
    from rdkit import Chem

    donors = pd.read_parquet(SCRATCH / "donors.parquet")
    sub = donors[(donors["pdb"] == pdb) & (donors["lig"] == lig)]
    if len(sub) == 0:
        return None, "no-crystal-row"
    r = sub.iloc[0]
    frame = np.frombuffer(r["frame"], np.float32).reshape(-1, 3).astype(float)
    cm = Chem.MolFromMolBlock(r["molblock"], sanitize=True)
    if cm is None:
        return None, "crystal-mol-unreadable"
    if cm.GetNumAtoms() != query_mol.GetNumAtoms():
        return None, "crystal-atom-count-mismatch"
    m = cm.GetSubstructMatch(query_mol)
    if not m or len(m) != query_mol.GetNumAtoms():
        m = query_mol.GetSubstructMatch(cm)
        if not m or len(m) != query_mol.GetNumAtoms():
            return None, "no-atom-order-map"
        inv = np.empty(len(m), int)
        inv[list(m)] = np.arange(len(m))
        return frame[inv], "ok"
    return frame[list(m)], "ok"


def _run_query(args):
    lig, smiles, pdb, loose = args[:4]
    want_xtal = len(args) > 4 and bool(args[4])
    import pandas as pd
    from rdkit import Chem, DataStructs, RDLogger
    RDLogger.DisableLog("rdApp.*")

    rec = _query_record(lig, smiles)
    if rec is None or rec["status"] != "ok":
        return {"ligand": lig, "status": (rec or {}).get("status", "no-frames")}

    xtal_status = "not-requested"
    if want_xtal:
        cf, xtal_status = _crystal_frame(lig, pdb, rec["mol"])
        if cf is not None:
            rec["frames"] = np.concatenate([rec["frames"], cf[None, :, :]], axis=0)

    donors = pd.read_parquet(SCRATCH / "donors.parquet")
    qfp = _fp(rec["mol"])
    min_atoms = MIN_MCS_ATOMS_LOOSE if loose else MIN_MCS_ATOMS

    # ---- leakage filters L1-L5, each counted -----------------------------
    n0 = len(donors)
    counts = {"start": n0}
    keep = np.ones(n0, bool)

    l1 = donors["pdb"].values == pdb
    counts["L1_same_entry"] = int(l1.sum())
    keep &= ~l1

    l2 = donors["lig"].values == lig
    counts["L2_same_ligand_code"] = int((l2 & keep).sum())
    keep &= ~l2

    # L3 - count-ECFP4 Tanimoto >= 0.90, one value per CCD code
    sim = {}
    for code, sub in donors.groupby("lig"):
        m = Chem.MolFromMolBlock(sub["molblock"].iloc[0], sanitize=True)
        sim[code] = (0.0 if m is None
                     else DataStructs.TanimotoSimilarity(qfp, _fp(m)))
    l3 = np.array([sim.get(c, 0.0) >= TANIMOTO_CUT for c in donors["lig"].values])
    counts["L3_tanimoto_ge_0.90"] = int((l3 & keep).sum())
    keep &= ~l3

    l5 = ~(donors["closest_fe"].values <= FE_CUT)
    counts["L5_not_active_site"] = int((l5 & keep).sum())
    keep &= ~l5

    unclass = donors["uniprot"].isna().values
    counts["L0_unclassified_target"] = int((unclass & keep).sum())
    keep &= ~unclass

    l4 = np.array([u in CYP3A_UNIPROT for u in donors["uniprot"].values])
    counts["L4_cyp3a_subfamily"] = int((l4 & keep).sum())

    out = {"ligand": lig, "status": "ok", "filters": counts, "names": rec["names"],
           "crystal_row": (len(rec["frames"]) - 1) if xtal_status == "ok" else None,
           "crystal_status": xtal_status}
    for arm, mask in (("no3a", keep & ~l4), ("with3a", keep)):
        idx = np.nonzero(mask)[0]
        per_donor, meta = [], []
        seen: dict[str, tuple] = {}
        for i in idx:
            row = donors.iloc[int(i)]
            code = row["lig"]
            if code not in seen:
                m = Chem.MolFromMolBlock(row["molblock"], sanitize=True)
                if m is None:
                    seen[code] = (None, 0)
                else:
                    seen[code] = _mcs_smarts(rec["mol"], m, loose)
            smarts, natoms = seen[code]
            if smarts is None or natoms < min_atoms:
                continue
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue
            dm = Chem.MolFromMolBlock(row["molblock"], sanitize=True)
            if dm is None:
                continue
            qi = _matches(rec["mol"], patt)
            di = _matches(dm, patt)
            if not qi or not di:
                continue
            frame = np.frombuffer(row["frame"], np.float32).reshape(-1, 3).astype(float)
            dist, npair = _transfer_distance(rec["frames"], qi, frame, di)
            per_donor.append(dist.astype(np.float32))
            meta.append({"pdb": row["pdb"], "chain": row["chain"], "lig": code,
                         "seqid": int(row["seqid"]), "uniprot": row["uniprot"],
                         "target_key": row["target_key"], "k": int(natoms),
                         "n_pairings": int(npair), "n_heavy": int(row["n_heavy"]),
                         "tanimoto": float(sim.get(code, 0.0)),
                         "frag_centroid": [float(v) for v in
                                           frame[di[0]].mean(axis=0)]})
        out[arm] = {"n_donors": len(per_donor), "meta": meta,
                    "d": (np.stack(per_donor, 1).tolist() if per_donor else [])}
    return out


def stage_match(workers: int, loose: bool, xtal: bool = False) -> None:
    import pandas as pd

    val = pd.read_csv(DATA / "validation_ligands.csv")
    tag = ("loose" if loose else "primary") + ("_xtal" if xtal else "")
    out = SCRATCH / f"match_{tag}"
    out.mkdir(parents=True, exist_ok=True)
    jobs = [(r.id, r.smiles, r.pdb, loose, xtal) for r in val.itertuples()
            if not (out / f"{r.id}.json").exists()]
    print(f"match[{tag}]: {len(jobs)} queries to do", flush=True)
    if not jobs:
        return
    if workers <= 1:
        for j in jobs:
            r = _run_query(j)
            (out / f"{j[0]}.json").write_text(json.dumps(r))
            print(f"  {j[0]}: {r.get('status')} "
                  f"no3a={r.get('no3a', {}).get('n_donors')}", flush=True)
        return
    import concurrent.futures as cf
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(_run_query, jobs):
            (out / f"{r['ligand']}.json").write_text(json.dumps(r))
            print(f"  {r['ligand']}: {r.get('status')} "
                  f"no3a={r.get('no3a', {}).get('n_donors')}", flush=True)


# ==========================================================================
# stage 4 - selection, nulls, and the paired test against the incumbent
# ==========================================================================

def _select(score_hi, truth, rng, draws=TIE_DRAWS):
    """Mean selected truth with ties broken at random, averaged over `draws`.

    `score_hi` and `truth` are lists of equal-length per-ligand arrays; HIGHER score
    is better. Returned per ligand so every downstream test can be paired.
    """
    n = len(score_hi)
    acc = np.zeros(n)
    picks = []
    for d in range(draws):
        p = []
        for i, (s, t) in enumerate(zip(score_hi, truth)):
            order = rng.permutation(len(s))
            j = order[np.argmax(s[order])]
            acc[i] += t[j]
            p.append(j)
        picks.append(p)
    return acc / draws, np.array(picks)


def _zwithin(v):
    sd = v.std()
    return (v - v.mean()) / (sd if sd > 1e-12 else 1.0)


def stage_evaluate(loose: bool) -> None:
    import pandas as pd
    from scipy import stats

    tag = "loose" if loose else "primary"
    rng = np.random.default_rng(SEED)

    scored = pd.read_csv(DATA / "poses_scored_val87b.csv")
    scored = scored[scored["arm"] == "unsteered"]
    truth = {(r.ligand, r.sample): float(r.lddt_pli) for r in scored.itertuples()}
    xeng = pd.read_csv(DATA / "xeng_val87b.csv")
    xe = {(r.ligand, r.sample): float(r.xeng) for r in xeng.itertuples()}
    donors_all = pd.read_parquet(SCRATCH / "donors.parquet")

    per = {}
    filt_tot: dict[str, int] = {}
    for f in sorted((SCRATCH / f"match_{tag}").glob("*.json")):
        r = json.loads(f.read_text())
        if r.get("status") != "ok":
            per[r["ligand"]] = r
            continue
        for k, v in r["filters"].items():
            filt_tot[k] = filt_tot.get(k, 0) + int(v)
        per[r["ligand"]] = r

    ligs = sorted(per)
    cov = {arm: np.array([per[l].get(arm, {}).get("n_donors", 0) for l in ligs])
           for arm in ("no3a", "with3a")}

    # ---------------- coverage FIRST ----------------------------------
    coverage = {
        "queries": len(ligs),
        "query_mol_failures": sorted(l for l in ligs if per[l].get("status") != "ok"),
        "filters_summed_over_queries": filt_tot,
        "donor_pool": {
            "observations": int(len(donors_all)),
            "codes": int(donors_all["lig"].nunique()),
            "entries": int(donors_all["pdb"].nunique()),
            "targets": int(donors_all["target_key"].nunique()),
        },
    }
    for arm in ("no3a", "with3a"):
        c = cov[arm]
        coverage[arm] = {
            "with_ge1_donor": int((c >= 1).sum()),
            "with_ge3_donors": int((c >= MIN_DONORS).sum()),
            "median_donors": float(np.median(c)),
            "mean_donors": float(c.mean()),
            "p90_donors": float(np.percentile(c, 90)),
            "max_donors": int(c.max()),
            "zero_donor_queries": int((c == 0).sum()),
        }
    print("=" * 72)
    print("COVERAGE (reported before any selector number)")
    print(json.dumps(coverage, indent=2))
    print("=" * 72, flush=True)

    results = {"arm": tag, "coverage": coverage, "check_frames": check_frames()}
    print("C4 frame check:", results["check_frames"], flush=True)

    # ---------------- assemble per-ligand arrays -----------------------
    def build(arm, subset_ligs):
        T, XE, F = [], [], {"cons": [], "best": [], "med": []}
        dmat, meta = [], []
        for l in subset_ligs:
            r = per[l]
            names = r["names"]
            d = np.asarray(r[arm]["d"], float)          # (poses, donors)
            T.append(np.array([truth[(l, n)] for n in names]))
            XE.append(np.array([xe[(l, n)] for n in names]))
            F["cons"].append(d.mean(1))
            F["best"].append(d.min(1))
            F["med"].append(np.median(d, 1))
            dmat.append(d)
            meta.append(r[arm]["meta"])
        return T, XE, F, dmat, meta

    for arm in ("no3a", "with3a"):
        subset = [l for l, c in zip(ligs, cov[arm]) if c >= MIN_DONORS]
        if not subset:
            results[arm] = {"subset": 0, "note": "no query has >= 3 legal donors"}
            continue
        T, XE, F, dmat, meta = build(arm, subset)
        n = len(subset)
        oracle = float(np.mean([t.max() for t in T]))
        baseline = float(np.mean([t.mean() for t in T]))

        block = {"subset_ligands": n, "oracle": oracle, "random_baseline": baseline}

        # C3 - the feature is not constant within a ligand
        sds = np.array([f.std() for f in F["cons"]])
        block["C3_nonconstant_frac"] = float((sds > 1e-9).mean())
        block["C3_min_within_sd"] = float(sds.min())

        # --- selectors -------------------------------------------------
        sel = {}
        for name, feats in (("frag_cons", F["cons"]), ("frag_best", F["best"]),
                            ("frag_med", F["med"])):
            s = [-f for f in feats]
            v, picks = _select(s, T, np.random.default_rng(SEED))
            sel[name] = (v, picks)
        inc_s = [-x for x in XE]
        inc_v, inc_p = _select(inc_s, T, np.random.default_rng(SEED))
        sel["incumbent"] = (inc_v, inc_p)
        combo = [_zwithin(-x) + _zwithin(-f) for x, f in zip(XE, F["cons"])]
        sel["frag_combo"] = _select(combo, T, np.random.default_rng(SEED))

        block["selection"] = {k: float(v.mean()) for k, (v, _) in sel.items()}
        block["gain_vs_random"] = {k: float(v.mean() - baseline)
                                   for k, (v, _) in sel.items()}

        # --- bar 1: random-feature null on THIS pool -------------------
        rg = np.random.default_rng(SEED + 1)
        null = np.empty(NULL_DRAWS)
        sizes = [len(t) for t in T]
        Tarr = T
        for i in range(NULL_DRAWS):
            tot = 0.0
            for t, m in zip(Tarr, sizes):
                tot += t[int(np.argmax(rg.standard_normal(m)))]
            null[i] = tot / n - baseline
        block["null_random_feature"] = {
            "draws": NULL_DRAWS, "p95": float(np.percentile(null, 95)),
            "p99": float(np.percentile(null, 99)), "max": float(null.max()),
            "sd": float(null.std())}
        block["bar1_pass"] = {
            k: bool(g > block["null_random_feature"]["p99"])
            for k, g in block["gain_vs_random"].items()}
        block["null_exceedance"] = {
            k: float((null >= g).mean())
            for k, g in block["gain_vs_random"].items()}

        # --- within-ligand rho and correct-sign fraction ---------------
        rhos = []
        for f, t in zip(F["cons"], T):
            if f.std() < 1e-12:
                continue
            rhos.append(stats.spearmanr(f, t).statistic)
        rhos = np.array([r for r in rhos if np.isfinite(r)])
        block["within_ligand"] = {
            "n": int(len(rhos)), "mean_rho": float(rhos.mean()),
            "median_rho": float(np.median(rhos)),
            "correct_sign_frac": float((rhos < 0).mean()),
            "sign_binom_p": float(stats.binomtest(int((rhos < 0).sum()),
                                                  len(rhos), 0.5).pvalue)}

        # --- bar 2: PAIRED against the incumbent ------------------------
        paired = {}
        for name in ("frag_cons", "frag_best", "frag_med", "frag_combo"):
            v, picks = sel[name]
            delta = v - inc_v
            tied = int(sum(int((picks[:, i] == inc_p[:, i]).all())
                           for i in range(n)))
            nz = delta[np.abs(delta) > 1e-12]
            bg = np.random.default_rng(SEED + 2)
            boot = np.array([delta[bg.integers(0, n, n)].mean()
                             for _ in range(BOOT_DRAWS)])
            w = (float(stats.wilcoxon(nz).pvalue) if len(nz) >= 1 else 1.0)
            paired[name] = {
                "mean_delta": float(delta.mean()),
                "boot_ci95": [float(np.percentile(boot, 2.5)),
                              float(np.percentile(boot, 97.5))],
                "wilcoxon_p": w, "n_informative": int(len(nz)),
                "ligands_tied_same_pose": tied,
                "better": int((delta > 1e-12).sum()),
                "worse": int((delta < -1e-12).sum())}
        block["paired_vs_incumbent"] = paired

        # --- bar 3: complementarity -------------------------------------
        head = np.array([t.max() for t in T]) - inc_v
        comp = {}
        for name in ("frag_cons", "frag_combo"):
            d = sel[name][0] - inc_v
            r = stats.pearsonr(d, head)
            comp[name] = {"pearson_r": float(r.statistic), "p": float(r.pvalue)}
        block["complementarity"] = comp

        # --- C1: the scrambled-donor control -----------------------------
        block["scrambled"] = _scrambled(arm, subset, per, donors_all, T, baseline,
                                        inc_v, meta)

        # --- donor provenance --------------------------------------------
        tks, ups = set(), {}
        for ms in meta:
            for m in ms:
                tks.add(m["target_key"])
                ups[m["uniprot"]] = ups.get(m["uniprot"], 0) + 1
        block["donor_provenance"] = {
            "distinct_targets_used": len(tks),
            "top_uniprots": sorted(ups.items(), key=lambda kv: -kv[1])[:10],
            "mean_fragment_atoms": float(np.mean([m["k"] for ms in meta for m in ms])),
            "median_fragment_atoms": float(np.median([m["k"] for ms in meta
                                                      for m in ms]))}
        results[arm] = block

    # whole-board number: uncovered queries fall back to the incumbent
    results["whole_board"] = _whole_board(ligs, cov, per, truth, xe)

    outp = DATA / f"frag_transfer_{tag}.json"
    outp.write_text(json.dumps(results, indent=2))
    print(json.dumps({k: v for k, v in results.items()
                      if k not in ("coverage",)}, indent=2)[:12000])
    print(f"\nwrote {outp}", flush=True)


def _scrambled(arm, subset, per, donors_all, T, baseline, inc_v, meta):
    """C1 - the same machinery with the CHEMISTRY removed.

    Each matched donor is replaced by a donor drawn uniformly at random from the same
    legal pool and by a random atom subset of the same size, and the minimum is taken
    over the same number of alternatives, so the scrambled arm keeps every structural
    advantage of the matched arm except knowing which atoms correspond. If matched is
    not better than this, the substructure match is decorative.
    """
    rng = np.random.default_rng(SEED + 7)
    pool = donors_all
    frames = [np.frombuffer(b, np.float32).reshape(-1, 3).astype(float)
              for b in pool["frame"].values]
    uni = pool["uniprot"].values
    legal_global = np.array([isinstance(u, str) and
                             (arm == "with3a" or u not in CYP3A_UNIPROT)
                             for u in uni])
    cand = np.nonzero(legal_global)[0]

    vals = []
    for rep in range(SCRAMBLE_REDRAWS):
        feats = []
        for li, l in enumerate(subset):
            f = np.load(SCRATCH / "frames" / f"{l}.npz", allow_pickle=True)
            pose = np.asarray(f["xyz"], float)
            ms = meta[li]
            acc = np.zeros((len(pose), len(ms)))
            for j, m in enumerate(ms):
                k, npair = m["k"], min(m["n_pairings"], MAX_MATCHES)
                ok = [c for c in rng.choice(cand, size=40, replace=False)
                      if len(frames[c]) >= k]
                c = ok[0] if ok else int(cand[0])
                df = frames[c]
                if len(df) < k:
                    acc[:, j] = np.nan
                    continue
                subs = np.stack([rng.choice(len(df), size=k, replace=False)
                                 for _ in range(npair)])
                B = df[subs]                                  # (npair, k, 3)
                qi = rng.choice(len(pose[0]), size=k, replace=False)
                A = pose[:, qi, :]                            # (P, k, 3)
                sq = ((A[:, None, :, :] - B[None, :, :, :]) ** 2).sum(-1).mean(-1)
                acc[:, j] = np.sqrt(sq.min(axis=1))
            feats.append(np.nanmean(acc, axis=1))
        v, _ = _select([-f for f in feats], T, np.random.default_rng(SEED + rep))
        vals.append(v)
    V = np.stack(vals)
    return {"redraws": SCRAMBLE_REDRAWS,
            "selected": float(V.mean()),
            "gain_vs_random": float(V.mean() - baseline),
            "sd_over_redraws": float(V.mean(axis=1).std()),
            "per_redraw_gain": [float(x - baseline) for x in V.mean(axis=1)]}


def _whole_board(ligs, cov, per, truth, xe):
    """All 87, uncovered queries falling back to the incumbent's pick."""
    out = {}
    for arm in ("no3a", "with3a"):
        rows_t, rows_s = [], []
        for l, c in zip(ligs, cov[arm]):
            r = per[l]
            if r.get("status") != "ok":
                continue
            names = r["names"]
            t = np.array([truth[(l, n)] for n in names])
            x = np.array([xe[(l, n)] for n in names])
            if c >= MIN_DONORS:
                d = np.asarray(r[arm]["d"], float).mean(1)
                s = -d
            else:
                s = -x
            rows_t.append(t)
            rows_s.append(s)
        v, _ = _select(rows_s, rows_t, np.random.default_rng(SEED))
        iv, _ = _select([-np.array([xe[(l, n)] for n in per[l]["names"]])
                         for l in ligs if per[l].get("status") == "ok"],
                        rows_t, np.random.default_rng(SEED))
        base = float(np.mean([t.mean() for t in rows_t]))
        out[arm] = {"ligands": len(rows_t), "random_baseline": base,
                    "oracle": float(np.mean([t.max() for t in rows_t])),
                    "frag_with_incumbent_fallback": float(v.mean()),
                    "incumbent": float(iv.mean()),
                    "mean_delta": float((v - iv).mean())}
    return out


# ==========================================================================
# stage 5 - diagnostics. POST-HOC, and labelled as such everywhere.
# ==========================================================================

def stage_diagnose() -> None:
    """Why a crystallographic prior cannot rank these poses. Five measurements.

    None of this is a candidate feature and none of it may be promoted to one without
    a fresh pre-registration. D5 in particular reads the query's own crystal, which is
    the answer: it is a diagnostic of the PRIOR, not of a predictor.

    | D1 | the donor ORACLE, against a matched best-of-N random-feature null |
    | D2 | gain as a function of keeping only the k most similar legal donors |
    | D3 | what the feature actually ranks by - centroid radius in the heme frame |
    | D4 | is there a consensus at all? the spread of donor fragment centroids |
    | D5 | the decisive one: where the query's own CRYSTAL pose falls on the feature |
    """
    import pandas as pd
    from scipy import stats

    rng = np.random.default_rng(SEED + 11)
    scored = pd.read_csv(DATA / "poses_scored_val87b.csv")
    scored = scored[scored["arm"] == "unsteered"]
    truth = {(r.ligand, r.sample): float(r.lddt_pli) for r in scored.itertuples()}

    prim = {f.stem: json.loads(f.read_text())
            for f in (SCRATCH / "match_primary").glob("*.json")}
    out: dict = {"note": "POST-HOC diagnostics; not pre-registered; not features"}

    rows = []
    for f in sorted((SCRATCH / "match_primary_xtal").glob("*.json")):
        r = json.loads(f.read_text())
        if r.get("status") != "ok" or r.get("crystal_row") is None:
            continue
        d = np.asarray(r["no3a"]["d"], float)
        if d.size == 0 or d.shape[1] < MIN_DONORS:
            continue
        rows.append((r, d))

    # C5 - the crystal re-run reproduces the primary run on the 20 predicted poses
    worst = 0.0
    for r, d in rows:
        p = prim.get(r["ligand"])
        if p is None or not p.get("no3a", {}).get("d"):
            continue
        dp = np.asarray(p["no3a"]["d"], float)
        if dp.shape == d[:dp.shape[0]].shape:
            worst = max(worst, float(np.abs(dp - d[:dp.shape[0]]).max()))
    out["C5_xtal_rerun_reproduces_primary"] = {
        "max_abs_delta": worst, "pass": bool(worst < 1e-5), "ligands": len(rows)}

    # ---- D5: where does the TRUTH fall on the feature? --------------------
    pct, better_than_median, deltas = [], 0, []
    for r, d in rows:
        k = r["crystal_row"]
        cons = d.mean(1)
        c_val, p_vals = cons[k], np.delete(cons, k)
        pct.append(float((p_vals > c_val).mean()))
        better_than_median += int(c_val < np.median(p_vals))
        deltas.append(float(c_val - p_vals.mean()))
    out["D5_crystal_vs_predicted"] = {
        "ligands": len(rows),
        "mean_percentile_of_crystal": float(np.mean(pct)),
        "median_percentile_of_crystal": float(np.median(pct)),
        "crystal_better_than_median_pose": better_than_median,
        "mean_delta_A": float(np.mean(deltas)),
        "sign_binom_p": float(stats.binomtest(better_than_median,
                                              len(rows), 0.5).pvalue),
        "reading": ("1.0 would mean the crystal is the single closest thing to the "
                    "donor consensus; 0.5 means the prior cannot tell the answer "
                    "from a wrong guess")}

    # ---- D4: is there a consensus at all? ---------------------------------
    spreads, own = [], []
    for r, d in rows:
        cents = np.array([m["frag_centroid"] for m in r["no3a"]["meta"]])
        spreads.append(float(np.sqrt(((cents - cents.mean(0)) ** 2).sum(1).mean())))
        fr = np.load(SCRATCH / "frames" / f"{r['ligand']}.npz", allow_pickle=True)
        c = np.asarray(fr["xyz"], float).mean(axis=1)
        own.append(float(np.sqrt(((c - c.mean(0)) ** 2).sum(1).mean())))
    out["D4_consensus_spread"] = {
        "donor_fragment_centroid_rms_spread_A": float(np.mean(spreads)),
        "median_A": float(np.median(spreads)),
        "own_20_pose_centroid_rms_spread_A": float(np.mean(own)),
        "mean_transfer_distance_A": float(np.mean([d.mean() for _, d in rows])),
        "reading": ("a prior whose donors scatter by several A about their own mean "
                    "cannot resolve a 30 deg rotation of a 4.8 A ligand")}

    # ---- D3: what does the feature rank by? -------------------------------
    rr, rg_ = [], []
    for r, d in rows:
        fr = np.load(SCRATCH / "frames" / f"{r['ligand']}.npz", allow_pickle=True)
        X = np.asarray(fr["xyz"], float)[:d.shape[0] - 1]
        cons = d.mean(1)[:len(X)]
        rad = np.linalg.norm(X.mean(1), axis=1)
        gyr = np.sqrt(((X - X.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        if rad.std() > 1e-9:
            rr.append(stats.spearmanr(cons, rad).statistic)
        if gyr.std() > 1e-9:
            rg_.append(stats.spearmanr(cons, gyr).statistic)
    rr = np.array([x for x in rr if np.isfinite(x)])
    rg_ = np.array([x for x in rg_ if np.isfinite(x)])
    out["D3_what_it_ranks_by"] = {
        "mean_rho_vs_centroid_radius": float(rr.mean()),
        "frac_positive_radius": float((rr > 0).mean()),
        "mean_rho_vs_radius_of_gyration": float(rg_.mean()),
        "frac_positive_gyration": float((rg_ > 0).mean())}

    # ---- D1: the donor oracle, against a matched best-of-N null -----------
    T, per_donor_best, ndon = [], [], []
    for r, d in rows:
        names = r["names"]
        t = np.array([truth[(r["ligand"], n)] for n in names])
        dd = d[:len(t)]
        T.append(t)
        per_donor_best.append(np.array([t[int(np.argmin(dd[:, j]))]
                                        for j in range(dd.shape[1])]))
        ndon.append(dd.shape[1])
    base = float(np.mean([t.mean() for t in T]))
    oracle_over_donors = float(np.mean([b.max() for b in per_donor_best]))
    nullbest = []
    for _ in range(200):
        v = []
        for t, m in zip(T, ndon):
            g = rng.standard_normal((m, len(t)))
            v.append(t[g.argmax(axis=1)].max())
        nullbest.append(np.mean(v))
    out["D1_donor_oracle"] = {
        "ligands": len(T), "random_baseline": base,
        "best_single_donor": oracle_over_donors,
        "gain": oracle_over_donors - base,
        "matched_best_of_N_random_null_mean": float(np.mean(nullbest)),
        "matched_null_p95": float(np.percentile(nullbest, 95)),
        "reading": ("with ~500 donors, best-of-500 random rankers is already this "
                    "high; a donor oracle that does not clear it is a counting "
                    "artifact, not chemistry")}

    # ---- D2: keep only the k most similar legal donors ---------------------
    curve = {}
    for k in (1, 3, 10, 30, 100, 300, 10 ** 9):
        feats, TT = [], []
        for r, d in rows:
            tan = np.array([m["tanimoto"] for m in r["no3a"]["meta"]])
            order = np.argsort(-tan)[:k]
            names = r["names"]
            t = np.array([truth[(r["ligand"], n)] for n in names])
            feats.append(d[:len(t)][:, order].mean(1))
            TT.append(t)
        v, _ = _select([-f for f in feats], TT, np.random.default_rng(SEED))
        curve["all" if k > 10 ** 8 else str(k)] = {
            "selected": float(v.mean()), "gain": float(v.mean() - base)}
    out["D2_similarity_ranked_donors"] = curve

    # ---- per-ligand table, so every pooled number above is auditable -------
    xe = pd.read_csv(DATA / "xeng_val87b.csv")
    xe = {(r.ligand, r.sample): float(r.xeng) for r in xe.itertuples()}
    tab = []
    for (r, d), t, pc in zip(rows, T, pct):
        lig = r["ligand"]
        names = r["names"]
        dd = d[:len(t)]
        cons = dd.mean(1)
        x = np.array([xe[(lig, n)] for n in names])
        fv, _ = _select([-cons], [t], np.random.default_rng(SEED))
        iv, _ = _select([-x], [t], np.random.default_rng(SEED))
        rho = (stats.spearmanr(cons, t).statistic if cons.std() > 1e-12 else np.nan)
        tab.append({"ligand": lig, "poses": len(t),
                    "n_donors_no3a": int(dd.shape[1]),
                    "n_donors_with3a": int(len(r["with3a"]["meta"])),
                    "mean_fragment_atoms": float(np.mean([m["k"]
                                                          for m in r["no3a"]["meta"]])),
                    "oracle": float(t.max()), "random_mean": float(t.mean()),
                    "frag_cons_selected": float(fv[0]),
                    "incumbent_selected": float(iv[0]),
                    "delta_frag_minus_incumbent": float(fv[0] - iv[0]),
                    "within_ligand_rho": float(rho),
                    "crystal_percentile_on_feature": pc,
                    "mean_transfer_distance_A": float(cons.mean())})
    pd.DataFrame(tab).to_csv(DATA / "frag_transfer_per_ligand.csv", index=False)

    p = DATA / "frag_transfer_diagnostics.json"
    p.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print(f"\nwrote {p}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=["frames", "donors", "match", "evaluate",
                                      "check", "diagnose"])
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--loose", action="store_true",
                    help="the pre-registered sensitivity arm")
    ap.add_argument("--xtal", action="store_true",
                    help="append the query's own crystal pose - DIAGNOSTIC ONLY")
    a = ap.parse_args()
    SCRATCH.mkdir(parents=True, exist_ok=True)
    if a.stage == "frames":
        stage_frames()
    elif a.stage == "donors":
        stage_donors()
    elif a.stage == "match":
        stage_match(a.workers, a.loose, a.xtal)
    elif a.stage == "check":
        print(json.dumps(check_frames(), indent=2))
    elif a.stage == "diagnose":
        stage_diagnose()
    else:
        stage_evaluate(a.loose)


if __name__ == "__main__":
    main()
