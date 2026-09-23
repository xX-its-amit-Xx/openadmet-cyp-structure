"""Two-ligand co-folding of CYP3A4 on OpenProtein — submit, collect, score.

Pre-registered in `docs/PREREG_two_ligand_cofold.md` (commit f6ec85c) before a single
job was submitted. Nothing in this file may move a threshold that document fixed.

**The question.** The biology map ranks eight candidate co-folding partners and puts a
**second copy of the query ligand** first: the only one with CYP3A4 structural precedent
(6 of 122 entries), the only one that touches the F/G roof, and one extra ligand entity
rather than a second chain. CYB5A and POR are DO NOT SPEND because they bind the
proximal face, 9.7 A from the iron across the porphyrin (FINDING 024).

**Three arms, identical but for the ligand entities.**

    single   protein + HEM + 1 x query
    double   protein + HEM + 2 x query          <- the hypothesis
    decoy    protein + HEM + 1 x query + 1 x ibuprofen   <- "does ANY second entity do this?"

**Design constraints that are measured facts, not guesses.**

* `diffusion_samples` does NOT sample the ligand on OpenProtein, for any engine
  (FINDING 009). So `samples=1`; depth would have to come from replicate jobs.
* protenix_v2 is deterministic across replicate jobs (FINDING 015). Replicates are
  therefore a CONTROL here, not depth — and determinism is an advantage for a paired
  design, because arm A and arm B differ by the intervention and by nothing else.
* An uploaded MSA makes four of the six working engines fail server-side (FINDING 009),
  so anything other than protenix_v2 / protenix / esmfold2 runs with `Protein.NullMSA`.
* Predicted and crystal residue numbering differ and the symptom is an exact 0.0
  LDDT-PLI (FINDING 021). Every pose is renumbered by residue-NAME agreement and the
  offset and identity are recorded per pose.

Usage:

    python scripts/cofold/two_ligand_cofold.py plan
    python scripts/cofold/two_ligand_cofold.py submit  --engine protenix_v2
    python scripts/cofold/two_ligand_cofold.py collect --engine protenix_v2
    python scripts/cofold/two_ligand_cofold.py score   --engine protenix_v2
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

OUT_ROOT = DATA_PROCESSED / "two_ligand"
JOBS = OUT_ROOT / "jobs.json"

USER = os.environ.get("OPENPROTEIN_USER", "shenoy.am@northeastern.edu")
PASS = os.environ.get("OPENPROTEIN_PASS", "Squack123!")

# ---------------------------------------------------------------------------
# everything below this line is fixed by the pre-registration
# ---------------------------------------------------------------------------

# Curated cryoprotectant / polymer component IDs. Same exclusion the biology map used
# over the whole PDB in §7.1. A second copy of a PEG fragment is not the hypothesis.
CRYO = {"PG0", "PG4", "PEG", "P6G", "1PE", "PE4", "EDO", "GOL", "MPD",
        "DMS", "SO4", "PO4", "ACT", "TRS", "IMD"}
N_PILOT = 15

# Ibuprofen. A CYP2C9 substrate, drug-like and lipophilic so it is a fair competitor for
# the cavity, and verified absent from the 87-ligand set by canonical SMILES.
DECOY_NAME = "ibuprofen"
DECOY_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"

ARMS = ("single", "double", "decoy")

# The measured peripheral-groove contact set (biology map §7.2: progesterone in 1W0F plus
# 5A1P/5A1R's Cys239). NOT the residue list the site is usually given - Leu211, Glu218,
# Phe215, Phe241 and Phe304 are REFUTED as peripheral-site residues there.
PERIPHERAL = [212, 213, 214, 217, 219, 220, 238, 239, 240]
FG_SPAN = list(range(210, 217))     # 210-216, FINDING 028's failure span
BACKBONE = ("N", "CA", "C", "O")

# Three-site geometry, biology map §7.1-7.3: active site 0-11 A, channel ~12-14 A,
# peripheral groove ~17-21 A from the iron.
ACTIVE_MAX_FE = 11.0
CHANNEL_MAX_FE = 15.0
PERIPHERAL_MIN_CONTACTS = 3
PERIPHERAL_CONTACT_RADIUS = 4.5


def pilot_ligands() -> tuple[pd.DataFrame, dict]:
    """The 15 pilot ligands, by the pre-registered rule. Prediction-only throughout."""
    from rdkit import Chem

    df = pd.read_csv(DATA_PROCESSED / "validation_ligands.csv")
    n_all = len(df)
    kept = df[~df.id.str.upper().isin(CRYO)].copy()
    n_cryo = n_all - len(kept)
    kept["n_heavy"] = [Chem.MolFromSmiles(s).GetNumHeavyAtoms() for s in kept.smiles]
    kept = kept.sort_values(["n_heavy", "id"]).head(N_PILOT).reset_index(drop=True)
    counts = {"n_validation_ligands": n_all, "n_excluded_cryoprotectant": int(n_cryo),
              "n_eligible": int(n_all - n_cryo), "n_pilot": int(len(kept))}
    return kept, counts


# ---------------------------------------------------------------------------
# submission
# ---------------------------------------------------------------------------

def connect():
    import openprotein
    return openprotein.connect(username=USER, password=PASS)


def build_complex(seq: str, smiles: str, arm: str, msa=None):
    """One protein chain, the heme by CCD code, and one or two ligand entities.

    `msa=None` means explicit single-sequence mode. Note the API gotcha recorded in
    `openprotein_cofold.py`: `Protein.single_sequence_mode` is a CLASS to be PASSED to
    `set_msa`, not a method to call.
    """
    from openprotein.molecules.chains import Ligand
    from openprotein.molecules.complex import Complex
    from openprotein.molecules.protein import Protein

    prot = Protein.from_expr(seq)
    prot.set_msa(msa if msa is not None else Protein.NullMSA)
    cx = Complex()
    cx.set_chain("A", prot)
    cx.set_chain("H", Ligand(ccd="HEM"))
    cx.set_chain("L", Ligand(smiles=smiles))
    if arm == "double":
        cx.set_chain("M", Ligand(smiles=smiles))
    elif arm == "decoy":
        cx.set_chain("M", Ligand(smiles=DECOY_SMILES))
    return cx


def _jobs() -> dict:
    return json.loads(JOBS.read_text()) if JOBS.exists() else {}


def _save_jobs(d: dict) -> None:
    JOBS.parent.mkdir(parents=True, exist_ok=True)
    tmp = JOBS.with_suffix(".tmp")
    tmp.write_text(json.dumps(d, indent=1))
    tmp.replace(JOBS)


def submit(engine: str, replicates: int, batch: int, single_sequence: bool,
           dry: bool = False) -> dict:
    from cypstruct import budget
    from cypstruct.targets import fetch_sequences

    df, counts = pilot_ligands()
    n_complexes = len(df) * len(ARMS) * replicates
    n_jobs = -(-n_complexes // batch)

    ok, why, _est = budget.preflight("openprotein_cofold", n_jobs,
                                     venue="openprotein", cap_key="openprotein_jobs")
    if not ok:
        return {"REFUSED_BY_PREFLIGHT": why, "n_jobs_requested": n_jobs}
    print(f"preflight OK: {n_jobs} jobs, {n_complexes} complexes "
          f"({len(df)} ligands x {len(ARMS)} arms x {replicates} replicates)", flush=True)
    if dry:
        return {"dry_run": True, "n_jobs": n_jobs, "counts": counts,
                "ligands": df.id.tolist()}

    s = connect()
    model = getattr(s.fold, engine)
    seq = fetch_sequences()["cyp3a4"]
    msa = None
    if not single_sequence:
        sys.path.insert(0, str(REPO / "scripts" / "cofold"))
        from openprotein_cofold import get_msa
        msa = get_msa(s)

    jobs = _jobs()
    key = engine
    jobs.setdefault(key, {"engine": engine, "single_sequence": single_sequence,
                          "batches": []})
    claimed = {(b["arm"], b["rep"], sid)
               for b in jobs[key]["batches"] for sid in b["ligands"]}

    n_sub = n_fail = 0
    for rep in range(replicates):
        for arm in ARMS:
            todo = [(r.id, r.smiles) for r in df.itertuples()
                    if (arm, rep, r.id) not in claimed]
            for i in range(0, len(todo), batch):
                chunk = todo[i:i + batch]
                sids = [sid for sid, _ in chunk]
                try:
                    fut = model.fold(
                        sequences=[build_complex(seq, smi, arm, msa)
                                   for _s, smi in chunk],
                        diffusion_samples=1, num_recycles=3)
                    jobs[key]["batches"].append(
                        {"job_id": str(fut.job_id), "ligands": sids, "arm": arm,
                         "rep": rep, "submitted": time.time()})
                    _save_jobs(jobs)
                    n_sub += 1
                    print(f"  {arm} r{rep} {','.join(sids)} -> {fut.job_id}", flush=True)
                except Exception as exc:
                    n_fail += 1
                    print(f"  {arm} r{rep} SUBMIT-FAIL {sids}: "
                          f"{type(exc).__name__}: {exc}", flush=True)
                time.sleep(0.6)

    budget.record(run_id=f"two_ligand_{engine}_{int(time.time())}", venue="openprotein",
                  kind="openprotein_cofold", units=n_sub, est=0.0,
                  note="FINDING 031 two-ligand pilot", n_jobs=n_sub, engine=engine,
                  arms=list(ARMS), replicates=replicates)
    return {"engine": engine, "submitted_jobs": n_sub, "submit_failures": n_fail,
            "counts": counts}


def _split_models(cif_text: str, sid: str, out: Path, arm: str, rep: int) -> list[str]:
    import gemmi
    st = gemmi.read_structure_string(cif_text)
    names = []
    for k in range(len(st)):
        one = st.clone()
        for j in reversed(range(len(one))):
            if j != k:
                del one[j]
        one.setup_entities()
        f = out / f"{sid}__{arm}__r{rep}s{k}.cif"
        f.write_text(one.make_mmcif_document().as_string())
        names.append(f.name)
    return names


def collect(engine: str) -> dict:
    """Fetch whatever finished. Resumable; a batch retires only when every ligand saved."""
    s = connect()
    out = OUT_ROOT / engine
    out.mkdir(parents=True, exist_ok=True)
    jobs = _jobs()
    if engine not in jobs:
        return {"error": f"no submitted jobs for {engine}"}

    n_ok = n_pending = n_fail = 0
    by_status: dict[str, int] = {}
    for b in jobs[engine]["batches"]:
        if b.get("done") is True:
            n_ok += len(b["ligands"])
            continue
        if b.get("done") == "failed":
            n_fail += len(b["ligands"])
            continue
        try:
            fut = s.load_job(b["job_id"])
            status = str(fut.job.status).upper()
            by_status[status] = by_status.get(status, 0) + 1
            if "SUCCESS" not in status:
                if "FAIL" in status or "CANCEL" in status:
                    b["done"] = "failed"
                    # read the failure message rather than guessing (README trap 2)
                    b["failure_message"] = str(getattr(fut.job, "failure_message", ""))[:400]
                    n_fail += len(b["ligands"])
                    print(f"  {b['job_id'][:8]} {status}: {b.get('failure_message','')}",
                          flush=True)
                else:
                    n_pending += len(b["ligands"])
                continue
            results = fut.get()
        except Exception as exc:
            print(f"  batch {b['job_id'][:8]}: {type(exc).__name__}: {exc}", flush=True)
            n_pending += len(b["ligands"])
            continue

        ok_all = True
        for idx, sid in enumerate(b["ligands"]):
            try:
                _split_models(results[idx].to_string(), sid, out, b["arm"], b["rep"])
                n_ok += 1
            except Exception as exc:
                print(f"  {sid} {b['arm']}: save failed ({type(exc).__name__}: {exc})",
                      flush=True)
                n_fail += 1
                ok_all = False
        if ok_all:
            b["done"] = True
    _save_jobs(jobs)
    return {"engine": engine, "collected": n_ok, "pending": n_pending,
            "failed": n_fail, "job_status": by_status}


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------

def renumber_to_reference(pred, ref):
    """FINDING 021. Shift prediction numbering onto the crystal's by residue-NAME
    agreement, which cannot manufacture agreement where there is none."""
    pm = {num: name for (_c, num), name in pred.prot_res.items()}
    rm = {num: name for (_c, num), name in ref.prot_res.items()}
    if len(pm) < 30 or len(rm) < 30:
        return pred, 0, 0.0
    best, best_id = 0, -1.0
    for off in range(int(min(rm) - max(pm)), int(max(rm) - min(pm)) + 1):
        shared = [n for n in pm if n + off in rm]
        if len(shared) < 30:
            continue
        ident = sum(pm[n] == rm[n + off] for n in shared) / len(shared)
        if ident > best_id:
            best, best_id = off, ident
    if best_id < 0:
        return pred, 0, 0.0
    if best != 0:
        pred.prot_key = [(c, n + best, a) for (c, n, a) in pred.prot_key]
        pred.prot_res = {(c, n + best): v for (c, n), v in pred.prot_res.items()}
    return pred, best, best_id


def load_crystal(pdb_id: str, ligand_code: str):
    import gemmi
    from cypstruct import pose as P
    from cypstruct.targets import fetch_cif
    cif = fetch_cif(pdb_id)
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    for chain in st[0]:
        if any(r.name.strip().upper() == ligand_code.upper() for r in chain):
            return P.load_structure(cif, ligand_code=ligand_code,
                                    assembly_chain=chain.name)
    return None


def ligand_entities(path: Path) -> list[dict]:
    """Every non-heme, non-solvent HET group in a prediction, with chain and coords.

    `pose.load_structure` deliberately returns ONE query ligand. A two-copy prediction
    has two, and which one is "the ligand" is a pre-registered decision, so the copies
    are enumerated here rather than left to the loader's largest-HET-group heuristic.
    """
    import gemmi
    from cypstruct.pose import IGNORE_HET, HEME_ALIASES, _looks_like_heme

    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    out = []
    for chain in st[0]:
        for res in chain:
            rn = res.name.strip().upper()
            info = gemmi.find_tabulated_residue(rn)
            if info and info.is_amino_acid():
                continue
            if rn in HEME_ALIASES or _looks_like_heme(res):
                continue
            if rn in IGNORE_HET:
                continue
            xyz = [[a.pos.x, a.pos.y, a.pos.z] for a in res
                   if a.element != gemmi.Element("H")]
            el = [a.element.name.upper() for a in res
                  if a.element != gemmi.Element("H")]
            if len(xyz) < 5:
                continue
            out.append({"chain": chain.name, "resname": rn,
                        "xyz": np.array(xyz, float), "elem": el})
    return out


def _kabsch_R(A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, float]:
    """Rotation taking centred A onto centred B, and the residual RMSD."""
    Ac, Bc = A - A.mean(0), B - B.mean(0)
    U, S, Vt = np.linalg.svd(Ac.T @ Bc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    res = float(np.sqrt(((Ac @ R.T - Bc) ** 2).sum(1).mean()))
    return R, res


def decompose(model_lig: np.ndarray, ref_lig: np.ndarray) -> dict:
    """FINDING 024's decomposition: rotation, translation, internal conformer.

    Both ligands are already in the crystal frame (the model was superposed on the
    binding site by residue number). The ligand is then optimally superposed on itself,
    so what is left of the rotation is the orientation error and what is left of the
    residual is the internal conformer difference.
    """
    R, conf = _kabsch_R(model_lig, ref_lig)
    ang = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0))))
    trans = float(np.linalg.norm(model_lig.mean(0) - ref_lig.mean(0)))
    raw = float(np.sqrt(((model_lig - ref_lig) ** 2).sum(1).mean()))
    return {"rot_deg": ang, "trans_a": trans, "conformer_a": conf, "raw_rmsd": raw}


def backbone_rmsd(a, b, resnums) -> tuple[float, int]:
    ka = {(r, at): a.prot_xyz[i] for i, (_c, r, at) in enumerate(a.prot_key)
          if r in resnums and at in BACKBONE}
    kb = {(r, at): b.prot_xyz[i] for i, (_c, r, at) in enumerate(b.prot_key)
          if r in resnums and at in BACKBONE}
    shared = sorted(set(ka) & set(kb))
    if len(shared) < 4:
        return float("nan"), len(shared)
    A = np.array([ka[k] for k in shared])
    B = np.array([kb[k] for k in shared])
    return float(np.sqrt(((A - B) ** 2).sum(1).mean())), len(shared)


def classify_second_copy(cx, other_xyz: np.ndarray) -> dict:
    """Active site / channel / peripheral groove / elsewhere, by the pre-registered rule."""
    fe_d = (float(np.linalg.norm(other_xyz - cx.fe, axis=1).min())
            if cx.fe is not None else float("nan"))
    contacts = 0
    per_res = []
    for r in PERIPHERAL:
        idx = [i for i, (_c, n, _a) in enumerate(cx.prot_key) if n == r]
        if not idx:
            continue
        d = np.linalg.norm(cx.prot_xyz[idx][:, None, :] - other_xyz[None, :, :],
                           axis=2).min()
        if d <= PERIPHERAL_CONTACT_RADIUS:
            contacts += 1
            per_res.append(r)
    if fe_d < ACTIVE_MAX_FE:
        cls = "active_site"
    elif fe_d <= CHANNEL_MAX_FE:
        cls = "channel"
    elif contacts >= PERIPHERAL_MIN_CONTACTS:
        cls = "peripheral_groove"
    else:
        cls = "elsewhere"
    return {"second_fe_dist": fe_d, "second_peripheral_contacts": contacts,
            "second_peripheral_residues": ",".join(map(str, per_res)), "second_class": cls}


def score(engine: str) -> dict:
    from cypstruct import pose as P
    from cypstruct import xengine as X

    refset = X.load_reference(DATA_PROCESSED / "reference_set_cyp3a4.npz")
    df, counts = pilot_ligands()
    meta = {r.id: r for r in df.itertuples()}
    out = OUT_ROOT / engine
    rows, skipped = [], {"no_crystal": 0, "no_pose_file": 0, "no_ligand": 0,
                         "mapping_failed": 0, "load_failed": 0}

    crystals = {}
    for sid, m in meta.items():
        c = load_crystal(m.pdb, sid)
        if c is None or len(c.lig_xyz) == 0:
            skipped["no_crystal"] += 1
        crystals[sid] = c

    for sid, m in meta.items():
        ref = crystals[sid]
        if ref is None:
            continue
        ref_pocket = P.pocket_residues_from_structure(ref, radius=8.0)
        n_ref = len(ref.lig_xyz)
        for arm in ARMS:
            for f in sorted(out.glob(f"{sid}__{arm}__r*.cif")):
                try:
                    cx = P.load_structure(f)
                except Exception:
                    skipped["load_failed"] += 1
                    continue
                cx, off, ident = renumber_to_reference(cx, ref)
                ents = ligand_entities(f)
                # the heme is filtered out inside ligand_entities; the remaining
                # entities are the query copy/copies plus (arm C) the decoy.
                if not ents:
                    skipped["no_ligand"] += 1
                    continue
                # candidate copies of the QUERY ligand: those matching its heavy count
                cand = [e for e in ents if len(e["xyz"]) == n_ref]
                if not cand:
                    # fall back to the closest heavy-atom count, and record it
                    cand = [min(ents, key=lambda e: abs(len(e["xyz"]) - n_ref))]
                fe = cx.fe
                if fe is None:
                    skipped["no_ligand"] += 1
                    continue
                fed = [float(np.linalg.norm(e["xyz"] - fe, axis=1).min()) for e in cand]

                per_copy = []
                for ci, e in enumerate(cand):
                    trial = P.Complex(
                        name=cx.name, prot_xyz=cx.prot_xyz, prot_key=list(cx.prot_key),
                        prot_res=dict(cx.prot_res), lig_xyz=e["xyz"],
                        lig_elem=list(e["elem"]), lig_name=e["resname"],
                        lig_chain=e["chain"], fe=cx.fe, heme_xyz=cx.heme_xyz,
                        heme_atom=list(cx.heme_atom), heme_elem=list(cx.heme_elem),
                        axial_sg=cx.axial_sg)
                    perm = P.best_ligand_mapping(m.smiles, trial, ref)
                    mapped = perm is not None
                    if perm is None:
                        perm = np.arange(min(len(trial.lig_xyz), n_ref))
                    ld = P.lddt_pli(trial, ref, lig_perm=perm)
                    bs = P.bisy_rmsd(trial, ref, align_resnums=ref_pocket, lig_perm=perm)
                    try:
                        aligned, fit, nfit = P.align_by_residue(trial, ref,
                                                               resnums=ref_pocket)
                        dec = decompose(aligned.lig_xyz[:len(perm)], ref.lig_xyz[perm])
                    except ValueError:
                        aligned, fit, nfit = None, float("nan"), 0
                        dec = {"rot_deg": np.nan, "trans_a": np.nan,
                               "conformer_a": np.nan, "raw_rmsd": np.nan}
                    # the shipped selector's feature: mean Chamfer, in the pose's OWN
                    # heme frame, to the frozen independent-engine references
                    hf = X.in_heme_frame(trial)
                    xe = (X.xeng_score(hf, refset.get(sid, []))
                          if hf is not None else float("nan"))
                    per_copy.append({"ci": ci, "chain": e["chain"], "fe_dist": fed[ci],
                                     "lddt_pli": ld, "bisy_rmsd": bs, "mapped": mapped,
                                     "pocket_ca_rmsd": fit, "n_fit": nfit,
                                     "xeng": xe, "aligned": aligned, **dec})
                    if not mapped:
                        skipped["mapping_failed"] += 1

                i_fe = int(np.argmin([c["fe_dist"] for c in per_copy]))
                i_best = int(np.argmax([c["lddt_pli"] for c in per_copy]))
                sel = per_copy[i_fe]

                # the OTHER entity, for the localisation question
                other = None
                if len(ents) > 1:
                    keep = {id(cand[i_fe])}
                    others = [e for e in ents if id(e) not in keep]
                    if others:
                        other = min(others, key=lambda e: float(
                            np.linalg.norm(e["xyz"] - fe, axis=1).min()))
                loc = (classify_second_copy(cx, other["xyz"]) if other is not None
                       else {"second_fe_dist": np.nan, "second_peripheral_contacts": 0,
                             "second_peripheral_residues": "", "second_class": "none"})
                if other is not None:
                    sel_xyz = cand[i_fe]["xyz"]
                    inter = float(np.linalg.norm(
                        sel_xyz[:, None, :] - other["xyz"][None, :, :], axis=2).min())
                else:
                    inter = float("nan")

                fg_ref, n_fg = backbone_rmsd(
                    per_copy[i_fe]["aligned"] if per_copy[i_fe]["aligned"] is not None
                    else cx, ref, set(FG_SPAN))

                rows.append({
                    "ligand": sid, "arm": arm, "file": f.name,
                    "rep": int(f.stem.split("__r")[1].split("s")[0]),
                    "n_entities": len(ents), "n_query_copies": len(cand),
                    "renumber_offset": off, "renumber_identity": round(ident, 4),
                    "lddt_pli": sel["lddt_pli"], "bisy_rmsd": sel["bisy_rmsd"],
                    "rot_deg": sel["rot_deg"], "trans_a": sel["trans_a"],
                    "conformer_a": sel["conformer_a"], "raw_rmsd": sel["raw_rmsd"],
                    "pocket_ca_rmsd": sel["pocket_ca_rmsd"], "n_fit": sel["n_fit"],
                    "mapped": sel["mapped"], "copy_chain": sel["chain"],
                    "fe_dist": sel["fe_dist"], "xeng": sel["xeng"],
                    "lddt_pli_bestmatch": per_copy[i_best]["lddt_pli"],
                    "bisy_bestmatch": per_copy[i_best]["bisy_rmsd"],
                    "same_copy_fe_and_best": i_fe == i_best,
                    "inter_copy_min_dist": inter,
                    "fg_bb_rmsd_to_crystal": fg_ref, "n_fg_atoms": n_fg,
                    "md5": hashlib.md5(f.read_bytes()).hexdigest(),
                    **loc,
                })
    res = pd.DataFrame(rows)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    res.to_csv(DATA_PROCESSED / f"two_ligand_poses_{engine}.csv", index=False)
    return {"engine": engine, "n_rows": len(res), "counts": counts,
            "skipped": skipped,
            "by_arm": res.groupby("arm").size().to_dict() if len(res) else {}}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["plan", "submit", "collect", "score"])
    ap.add_argument("--engine", default="protenix_v2")
    ap.add_argument("--replicates", type=int, default=2)
    ap.add_argument("--batch", type=int, default=4)
    ap.add_argument("--single-sequence", action="store_true")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()

    if a.cmd == "plan":
        df, counts = pilot_ligands()
        print(json.dumps(counts, indent=1))
        print(df[["id", "n_heavy", "pdb", "cls"]].to_string(index=False))
    elif a.cmd == "submit":
        print(json.dumps(submit(a.engine, a.replicates, a.batch,
                                a.single_sequence, a.dry), indent=1)[:2000])
    elif a.cmd == "collect":
        print(json.dumps(collect(a.engine), indent=2))
    else:
        print(json.dumps(score(a.engine), indent=2, default=str))
