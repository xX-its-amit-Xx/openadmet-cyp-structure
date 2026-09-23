"""Pool DEPTH on the predicted-Type-I stratum — submit, collect, yield, score.

Pre-registered in `docs/PREREG_type_i_depth.md` (committed before the first job).
Nothing in this file may move a threshold that document fixed.

**The question.** FINDING 033 recommended spending remaining budget on pool depth for
the predicted-Type-I ligands, because their oracle is 0.058 below Type II and selection
cannot reach what the pool does not contain. FINDING 004 measured depth at +0.0125 of
selected score per doubling. Four separate pool EXPANSIONS (013, 016, 027, 031) each
raised the oracle and left selection flat or worse. This is depth, and the distinction
is under test.

**Design constraints that are measured facts, not guesses.**

* `diffusion_samples` does NOT sample the ligand on OpenProtein, for any engine
  (FINDING 009) -> `diffusion_samples=1`, depth comes from replicate jobs.
* `boltz2` with an UPLOADED MSA fails server-side on OpenProtein ("internal server
  error", verified on the stored probe job) while single-sequence `boltz2` succeeds.
  Single-sequence mode is therefore forced, and it is confound C1 of the prereg.
* Determinism is a property of the CONFIGURATION, not of the engine (FINDING 031 found
  protenix_v2 bimodal where FINDING 015 found it deterministic), so the distinct-pose
  yield of THIS configuration is measured before any depth is bought.
* Predicted and crystal residue numbering differ and the symptom is an exact 0.0
  LDDT-PLI (FINDING 021). Every pose is renumbered by residue-NAME agreement and the
  offset and identity are recorded per pose.

Usage:

    python scripts/cofold/type_i_depth.py plan
    python scripts/cofold/type_i_depth.py submit --replicates 6 --probe
    python scripts/cofold/type_i_depth.py collect
    python scripts/cofold/type_i_depth.py yield
    python scripts/cofold/type_i_depth.py submit --replicates 24
    python scripts/cofold/type_i_depth.py score
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
sys.path.insert(0, str(REPO / "scripts" / "cofold"))

from cypstruct.paths import DATA_PROCESSED, SCRATCH  # noqa: E402

ENGINE = "boltz2"
OUT_ROOT = SCRATCH / "type_i_depth"          # C:/cyp_struct - never D:
JOBS = DATA_PROCESSED / "type_i_depth_jobs.json"
POSES_CSV = DATA_PROCESSED / "type_i_depth_poses.csv"

USER = os.environ.get("OPENPROTEIN_USER", "shenoy.am@northeastern.edu")
PASS = os.environ.get("OPENPROTEIN_PASS", "REDACTED")

# ---------------------------------------------------------------------------
# fixed by the pre-registration
# ---------------------------------------------------------------------------
COORD_MAX = 2.6          # repo constant, build_reference_set.py. Do not tune.
PROBE_LIGANDS = ["MWS", "YNV"]   # one high-quality, one low-quality Type I; fixed a priori
DEDUPE_TOL = 0.05        # cypstruct.xengine._dedupe default, A


def type_i_ligands() -> pd.DataFrame:
    """The predicted-Type-I stratum, by the pre-registered prediction-side rule.

    median fe_donor_dist over the ligand's OWN 20 existing poses > 2.6 A. Nothing from
    any crystal enters this selection.
    """
    L = pd.read_csv(DATA_PROCESSED / "binding_mode_labels_cyp3a4.csv")
    sel = L[L.pred_fe_donor_median > COORD_MAX].copy()
    sel = sel.sort_values("pred_fe_donor_median").reset_index(drop=True)
    return sel


def counts() -> dict:
    L = pd.read_csv(DATA_PROCESSED / "binding_mode_labels_cyp3a4.csv")
    sel = type_i_ligands()
    return {
        "n_validation_ligands": int(len(L)),
        "filter_pred_fe_donor_median_gt_2.6": int(len(sel)),
        "excluded_predicted_type_II": int((L.pred_fe_donor_median <= COORD_MAX).sum()),
        "label_disagreements_total": int((~L.pred_agrees_crystal).sum()),
        "disagreements_in_stratum": sorted(
            sel.loc[~sel.pred_agrees_crystal, "id"].tolist()),
        "disagreements_excluded": sorted(
            L.loc[~L.pred_agrees_crystal & (L.pred_fe_donor_median <= COORD_MAX),
                  "id"].tolist()),
    }


# ---------------------------------------------------------------------------
# submission  (reuses two_ligand_cofold's job-ledger / resume machinery)
# ---------------------------------------------------------------------------

def connect():
    import openprotein
    return openprotein.connect(username=USER, password=PASS)


def build_complex(seq: str, smiles: str):
    """protein + CCD heme + one ligand, explicit single-sequence mode.

    `Protein.single_sequence_mode` is a CLASS to be PASSED to `set_msa`, not a method
    to call (openprotein_cofold.py records how that one cost a tranche).
    """
    from openprotein.molecules.chains import Ligand
    from openprotein.molecules.complex import Complex
    from openprotein.molecules.protein import Protein

    prot = Protein.from_expr(seq)
    prot.set_msa(Protein.NullMSA)
    cx = Complex()
    cx.set_chain("A", prot)
    cx.set_chain("H", Ligand(ccd="HEM"))
    cx.set_chain("L", Ligand(smiles=smiles))
    return cx


def _jobs() -> dict:
    return json.loads(JOBS.read_text()) if JOBS.exists() else {}


def _save_jobs(d: dict) -> None:
    JOBS.parent.mkdir(parents=True, exist_ok=True)
    tmp = JOBS.with_suffix(".tmp")
    tmp.write_text(json.dumps(d, indent=1))
    tmp.replace(JOBS)


def submit(replicates: int, batch: int, probe: bool = False, dry: bool = False) -> dict:
    from cypstruct import budget
    from cypstruct.targets import fetch_sequences

    df = type_i_ligands()
    if probe:
        df = df[df.id.isin(PROBE_LIGANDS)].reset_index(drop=True)

    jobs = _jobs()
    key = ENGINE
    jobs.setdefault(key, {"engine": ENGINE, "single_sequence": True, "batches": []})
    claimed = {(b["rep"], sid) for b in jobs[key]["batches"] for sid in b["ligands"]}

    todo_pairs = [(rep, r.id, r.smiles) for rep in range(replicates)
                  for r in df.itertuples() if (rep, r.id) not in claimed]
    n_complexes = len(todo_pairs)
    # jobs are batched WITHIN a replicate, so count them that way rather than dividing
    # the total - a preflight that under-counts is not a preflight.
    n_jobs = sum(-(-sum(1 for (r, _s, _m) in todo_pairs if r == rep) // batch)
                 for rep in range(replicates))

    ok, why, _est = budget.preflight("openprotein_cofold", max(n_jobs, 1),
                                     venue="openprotein", cap_key="openprotein_jobs")
    if not ok:
        return {"REFUSED_BY_PREFLIGHT": why, "n_jobs_requested": n_jobs}
    print(f"preflight OK: {n_jobs} jobs, {n_complexes} new complexes "
          f"({len(df)} ligands x {replicates} replicates, {len(claimed)} already claimed)",
          flush=True)
    if dry or not n_complexes:
        return {"dry_run": True, "n_jobs": n_jobs, "n_complexes": n_complexes,
                "ligands": df.id.tolist(), "counts": counts()}

    s = connect()
    model = getattr(s.fold, ENGINE)
    seq = fetch_sequences()["cyp3a4"]

    n_sub = n_fail = 0
    # batch WITHIN a replicate so a job never mixes replicate indices
    for rep in range(replicates):
        chunk_src = [(sid, smi) for (r, sid, smi) in todo_pairs if r == rep]
        for i in range(0, len(chunk_src), batch):
            chunk = chunk_src[i:i + batch]
            sids = [sid for sid, _ in chunk]
            try:
                fut = model.fold(
                    sequences=[build_complex(seq, smi) for _s, smi in chunk],
                    diffusion_samples=1, num_recycles=3)
                jobs[key]["batches"].append(
                    {"job_id": str(fut.job_id), "ligands": sids, "rep": rep,
                     "submitted": time.time(), "probe": bool(probe)})
                _save_jobs(jobs)
                n_sub += 1
                print(f"  r{rep} {','.join(sids)} -> {fut.job_id}", flush=True)
            except Exception as exc:
                n_fail += 1
                print(f"  r{rep} SUBMIT-FAIL {sids}: {type(exc).__name__}: {exc}",
                      flush=True)
            time.sleep(0.6)

    budget.record(run_id=f"type_i_depth_{ENGINE}_{int(time.time())}",
                  venue="openprotein", kind="openprotein_cofold", units=n_sub, est=0.0,
                  note="FINDING 034 Type-I depth", n_jobs=n_sub, engine=ENGINE,
                  replicates=replicates, probe=bool(probe))
    return {"engine": ENGINE, "submitted_jobs": n_sub, "submit_failures": n_fail,
            "n_complexes": n_complexes, "counts": counts()}


def _split_models(cif_text: str, sid: str, out: Path, rep: int) -> list[str]:
    import gemmi
    st = gemmi.read_structure_string(cif_text)
    names = []
    for k in range(len(st)):
        one = st.clone()
        for j in reversed(range(len(one))):
            if j != k:
                del one[j]
        one.setup_entities()
        f = out / f"{sid}__r{rep}s{k}.cif"
        f.write_text(one.make_mmcif_document().as_string())
        names.append(f.name)
    return names


def collect() -> dict:
    """Fetch whatever finished. Resumable; a batch retires only when every ligand saved."""
    from cypstruct.paths import guard_scratch
    guard_scratch()
    s = connect()
    out = OUT_ROOT / ENGINE
    out.mkdir(parents=True, exist_ok=True)
    jobs = _jobs()
    if ENGINE not in jobs:
        return {"error": f"no submitted jobs for {ENGINE}"}

    n_ok = n_pending = n_fail = 0
    by_status: dict[str, int] = {}
    for b in jobs[ENGINE]["batches"]:
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
                    b["failure_message"] = str(
                        getattr(fut.job, "failure_message", ""))[:400]
                    n_fail += len(b["ligands"])
                    print(f"  {b['job_id'][:8]} {status}: "
                          f"{b.get('failure_message','')}", flush=True)
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
                _split_models(results[idx].to_string(), sid, out, b["rep"])
                n_ok += 1
            except Exception as exc:
                print(f"  {sid} r{b['rep']}: save failed "
                      f"({type(exc).__name__}: {exc})", flush=True)
                n_fail += 1
                ok_all = False
        if ok_all:
            b["done"] = True
    _save_jobs(jobs)
    return {"engine": ENGINE, "collected": n_ok, "pending": n_pending,
            "failed": n_fail, "job_status": by_status}


# ---------------------------------------------------------------------------
# distinct-pose yield — measured BEFORE any depth is bought (prereg section 3)
# ---------------------------------------------------------------------------

def pose_yield(verbose: bool = True) -> dict:
    """Distinct poses per replicate for THIS configuration.

    Counting distinct POSES, never distinct jobs or distinct file hashes: FINDING 009's
    20 files all had distinct md5s and shared one ligand conformation.
    """
    from cypstruct import pose as P
    from cypstruct import xengine as X

    out = OUT_ROOT / ENGINE
    rows = []
    per_lig: dict[str, dict[int, np.ndarray]] = {}
    for f in sorted(out.glob("*__r*.cif")):
        lig = f.name.split("__r")[0]
        rep = int(f.name.split("__r")[1].split("s")[0])
        try:
            v = X.in_heme_frame(P.load_structure(f))
        except Exception as exc:
            rows.append({"ligand": lig, "rep": rep, "error": f"{type(exc).__name__}"})
            continue
        if v is None:
            rows.append({"ligand": lig, "rep": rep, "error": "no heme frame"})
            continue
        per_lig.setdefault(lig, {})[rep] = v

    res = {}
    for lig, d in sorted(per_lig.items()):
        poses = [d[r] for r in sorted(d)]
        distinct = X._dedupe(poses, tol=DEDUPE_TOL)
        # pairwise per-atom max deviation, as a spread statistic
        spreads = []
        for i in range(len(poses)):
            for j in range(i + 1, len(poses)):
                if len(poses[i]) == len(poses[j]):
                    spreads.append(float(np.abs(poses[i] - poses[j]).max()))
        res[lig] = {"n_replicates": len(poses), "n_distinct": len(distinct),
                    "yield": len(distinct) / max(1, len(poses)),
                    "median_pairwise_max_dev": (round(float(np.median(spreads)), 4)
                                                if spreads else None)}
    ys = [v["yield"] for v in res.values()]
    summary = {"engine": ENGINE, "per_ligand": res,
               "mean_yield": round(float(np.mean(ys)), 4) if ys else None,
               "min_yield": round(float(np.min(ys)), 4) if ys else None,
               "n_ligands": len(res), "dedupe_tol": DEDUPE_TOL,
               "read_errors": [r for r in rows if "error" in r]}
    if ys:
        y = float(np.mean(ys))
        summary["depth_capable"] = bool(y >= 0.25)
        summary["replicates_for_20_distinct"] = int(np.ceil(20 / max(y, 0.1)))
    if verbose:
        print(json.dumps(summary, indent=2))
    return summary


# ---------------------------------------------------------------------------
# scoring  (renumber + LDDT-PLI + geometry + the shipped xeng column)
# ---------------------------------------------------------------------------

def score() -> dict:
    from cypstruct import pose as P
    from cypstruct import xengine as X
    from cypstruct.qmscore import geometry as G
    from two_ligand_cofold import load_crystal, renumber_to_reference

    refset = X.load_reference(DATA_PROCESSED / "reference_set_cyp3a4.npz")
    df = type_i_ligands()
    vl = pd.read_csv(DATA_PROCESSED / "validation_ligands.csv").set_index("id")
    out = OUT_ROOT / ENGINE
    rows = []
    skipped = {"no_crystal": 0, "no_ligand": 0, "load_failed": 0,
               "mapping_failed": 0, "no_heme_frame": 0, "low_renumber_identity": 0}

    for r in df.itertuples():
        sid = r.id
        smiles = vl.loc[sid, "smiles"]
        pdb = vl.loc[sid, "pdb"]
        ref = load_crystal(pdb, sid)
        if ref is None or len(ref.lig_xyz) == 0:
            skipped["no_crystal"] += 1
            continue
        ref_pocket = P.pocket_residues_from_structure(ref, radius=8.0)
        files = sorted(out.glob(f"{sid}__r*.cif"))
        for f in files:
            try:
                cx = P.load_structure(f)
            except Exception:
                skipped["load_failed"] += 1
                continue
            if len(cx.lig_xyz) == 0:
                skipped["no_ligand"] += 1
                continue
            cx, off, ident = renumber_to_reference(cx, ref)
            if ident < 0.95:
                skipped["low_renumber_identity"] += 1
            perm = P.best_ligand_mapping(smiles, cx, ref)
            mapped = perm is not None
            if perm is None:
                skipped["mapping_failed"] += 1
                perm = np.arange(min(len(cx.lig_xyz), len(ref.lig_xyz)))
            ld = P.lddt_pli(cx, ref, lig_perm=perm)
            bs = P.bisy_rmsd(cx, ref, align_resnums=ref_pocket, lig_perm=perm)
            geo = G.compute(cx.lig_xyz, cx.lig_elem, cx.prot_xyz,
                            cx.heme_xyz, cx.heme_atom, cx.axial_sg)
            hf = X.in_heme_frame(cx)
            if hf is None:
                skipped["no_heme_frame"] += 1
            xe = (X.xeng_score(hf, refset.get(sid, [])) if hf is not None
                  else float("nan"))
            rows.append({
                "ligand": sid, "pdb": pdb, "arm": "op_boltz2_ss",
                "sample": f.stem.split("__")[1],
                "rep": int(f.stem.split("__r")[1].split("s")[0]),
                "lddt_pli": ld, "bisy_rmsd": bs, "xeng": xe,
                "fe_donor_dist": geo.fe_donor_dist,
                "s_fe_donor_angle": geo.s_fe_donor_angle,
                "is_coordinated": geo.is_coordinated,
                "frac_proximal": geo.frac_atoms_proximal,
                "renumber_offset": off, "renumber_identity": round(ident, 4),
                "mapped": mapped,
                "md5": hashlib.md5(f.read_bytes()).hexdigest(),
                "n_ref_poses": len(refset.get(sid, [])),
            })
    res = pd.DataFrame(rows)
    res.to_csv(POSES_CSV, index=False)
    return {"engine": ENGINE, "n_rows": int(len(res)),
            "n_ligands": int(res.ligand.nunique()) if len(res) else 0,
            "skipped": skipped, "counts": counts(),
            "renumber_offsets": (res.renumber_offset.value_counts().to_dict()
                                 if len(res) else {}),
            "min_renumber_identity": (float(res.renumber_identity.min())
                                      if len(res) else None),
            "n_exact_zero_lddt": int((res.lddt_pli == 0).sum()) if len(res) else 0,
            "mean_lddt": float(res.lddt_pli.mean()) if len(res) else None}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["plan", "submit", "collect", "yield", "score"])
    ap.add_argument("--replicates", type=int, default=6)
    ap.add_argument("--batch", type=int, default=7)
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()

    if a.cmd == "plan":
        print(json.dumps(counts(), indent=1))
        print(type_i_ligands()[["id", "pdb", "crystal_mode", "pred_mode",
                                "pred_fe_donor_median", "pool_mean_lddt"]]
              .to_string(index=False))
    elif a.cmd == "submit":
        print(json.dumps(submit(a.replicates, a.batch, a.probe, a.dry),
                         indent=1, default=str)[:4000])
    elif a.cmd == "collect":
        print(json.dumps(collect(), indent=2))
    elif a.cmd == "yield":
        pose_yield()
    else:
        print(json.dumps(score(), indent=2, default=str))
