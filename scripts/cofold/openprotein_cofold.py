"""Co-folding on OpenProtein.ai — the venue that is not rate-limited by our budget.

**Why this matters more than another Modal runner.** Modal is capped and reserved for
fine-tuning; OpenProtein is unlimited for this account. It also exposes **eleven engines**
behind one API (`alphafold2, boltz-1, boltz-1x, boltz-2, esmfold, esmfold2, esmfold2-fast,
minifold, protenix, protenix-v2, rosettafold-3`), which matters because of FINDING 005:
a second engine is worth nothing unless its errors are **decorrelated** from Boltz's.
Chai-1 was not (rho = +0.45) and contributed nothing for $22. RoseTTAFold-3 and Protenix-v2
are architecturally furthest from Boltz-2 and so are the best candidates.

**The heme goes in properly here.** `Ligand(ccd="HEM")` takes a chemical-component code, so
unlike Chai-1 — which has no CCD input and had to receive the cofactor as SMILES — the
cofactor arrives with its ideal geometry, as it does for Boltz.

**Standing rule from FINDING 005, enforced by `pilot`:** before buying a full pool from a
new engine, run ~20 ligands and check per-ligand oracle correlation against the incumbent.
Decorrelation is the only thing that makes a second engine worth anything.

Measured engine support (op_probe.py, 2026-09-13): **protenix-v2, protenix and esmfold2
work**; boltz-1, boltz-1x, boltz-2 and rosettafold-3 all fail at runtime on
protein+HEM+ligand; alphafold2 warns that it discards ligand chains. esmfold2 was first
recorded as failing - that was reading a still-RUNNING job as a failure, and it in fact
returns a coordinated complex (Fe-ligand 2.28 A against Protenix's 2.43 A). So the decorrelation candidate is
Protenix, which is at least architecturally distinct from Boltz-2.

    python scripts/cofold/openprotein_cofold.py submit  --samples 20        # returns at once
    python scripts/cofold/openprotein_cofold.py collect                     # poll, resumable
    python scripts/cofold/openprotein_cofold.py score                       # FINDING 005 gate
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

OUT_ROOT = DATA_PROCESSED / "openprotein"
USER = os.environ.get("OPENPROTEIN_USER", "shenoy.am@northeastern.edu")
PASS = os.environ.get("OPENPROTEIN_PASS", "Squack123!")

# Engines worth trying, in order of how architecturally DIFFERENT they are from Boltz-2 —
# which is the only property that matters per FINDING 005. boltz2 is included last purely
# as a positive control: it should correlate strongly with our own Boltz pool, and if it
# does not, the comparison method is broken rather than the engine being interesting.
ENGINES = ["rosettafold_3", "protenix_v2", "esmfold2", "alphafold2", "boltz2"]


def connect():
    import openprotein

    return openprotein.connect(username=USER, password=PASS)


# API gotcha: `Protein.single_sequence_mode` is not a method, it is a CLASS (an alias of
# `Protein.NullMSA`) that you pass to `set_msa`. Calling it as `p.single_sequence_mode()`
# silently constructs an instance, leaves the MSA unset, and the submission is then
# rejected with the very error message that names the attribute you just called.
#     WRONG: p.single_sequence_mode()
#     RIGHT: p.set_msa(Protein.NullMSA)


def build_complex(seq: str, smiles: str, msa=None):
    """`msa=None` means explicit single-sequence mode, which some engines REQUIRE.

    RoseTTAFold-3 and boltz-2 fail with a server-side "internal server error" when handed
    an UPLOADED msa, and rosettafold-3 succeeds immediately without one. So the earlier
    "these engines are broken" conclusion was wrong: they are broken on `upload_msa`
    output specifically, not on this complex.
    """
    from openprotein.molecules.chains import Ligand
    from openprotein.molecules.complex import Complex
    from openprotein.molecules.protein import Protein

    prot = Protein.from_expr(seq)
    prot.set_msa(msa if msa is not None else Protein.NullMSA)
    cx = Complex()
    cx.set_chain("A", prot)
    cx.set_chain("H", Ligand(ccd="HEM"))     # cofactor with ideal geometry, not SMILES
    cx.set_chain("L", Ligand(smiles=smiles))
    return cx


def ligand_set(limit: int | None, seed: int = 0) -> pd.DataFrame:
    df = pd.read_csv(DATA_PROCESSED / "validation_ligands.csv")
    if limit:
        # Sample across the whole set rather than taking the head: the csv is sorted by
        # ligand code, so a head() slice would be a chemically biased subset and the
        # decorrelation estimate would not generalise.
        df = df.sample(n=min(limit, len(df)), random_state=seed).reset_index(drop=True)
    return df


JOBS = OUT_ROOT / "jobs.json"


def _jobs() -> dict:
    return json.loads(JOBS.read_text()) if JOBS.exists() else {}


def _save_jobs(d: dict) -> None:
    JOBS.parent.mkdir(parents=True, exist_ok=True)
    JOBS.write_text(json.dumps(d, indent=1))


def get_msa(s):
    """The shared CYP3A4 alignment, uploaded once and reused by every fold.

    We do NOT recompute it here. The 6,979-sequence alignment already computed for the
    Boltz campaign is uploaded instead, which takes 2 seconds against several minutes for
    a fresh search. It needs one conversion first: a3m marks insertions relative to the
    query with lowercase letters, so rows are ragged, and OpenProtein rejects a
    non-uniform MSA ("Expected uniform length, found 503-779"). Dropping the lowercase
    columns restores a uniform alignment at the query length.
    """
    import io

    meta = OUT_ROOT / "msa_job.json"
    if meta.exists():
        return s.load_job(json.loads(meta.read_text())["job_id"])

    aligned = REPO / "data" / "reference" / "cyp3a4_aligned.fasta"
    if not aligned.exists():
        raw = (REPO / "data" / "reference" / "cyp3a4.a3m").read_text(errors="replace")
        recs, name, buf = [], None, []
        for ln in raw.splitlines():
            if ln.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = ln[1:].strip(), []
            elif ln.strip():
                buf.append(ln.strip())
        if name is not None:
            recs.append((name, "".join(buf)))
        clean = [(n, "".join(c for c in q if not c.islower())) for n, q in recs]
        L = len(clean[0][1])
        aligned.write_text(
            "".join(">%s\n%s\n" % (n, q) for n, q in clean if len(q) == L))

    msa = s.align.upload_msa(io.BytesIO(aligned.read_bytes()))
    msa.wait_until_done()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    meta.write_text(json.dumps({"job_id": str(msa.job_id)}))
    return msa


def submit(engine: str, df: pd.DataFrame, samples: int, tag: str,
           batch: int = 4, replicates: int = 1, single_sequence: bool = False) -> dict:
    """Submit folds, RECORD THE JOB IDS, and return without waiting.

    Three facts about this API, each measured by `op_probe.py` rather than assumed, and
    each of which changes the shape of the campaign:

    * **A protein chain must carry an MSA** (or be put in explicit single-sequence mode).
      There is no implicit search. The first probe round submitted eight arms without one
      and every single one was rejected at submit time.
    * **One job takes several complexes.** The result list indexes over complexes, not
      over samples, so 87 ligands is ~22 jobs rather than 87.
    * **`diffusion_samples=N` returns N models inside one structure**, not N structures.
      `collect` splits them.

    **And the one that decides the whole campaign shape (FINDING 009):**
    `diffusion_samples` does NOT sample the ligand. Across 20 models of one job the ligand
    and heme coordinates are byte-identical (per-atom sd 0.0000 A) and only the protein
    moves. Two SEPARATE jobs for the same ligand differ by 8.63 A ligand RMSD. So pose
    diversity comes from REPLICATE JOBS, and `replicates` - not `samples` - is the knob
    that builds a pool. A run with samples=20, replicates=1 yields one ligand pose.
    """
    from cypstruct.targets import fetch_sequences

    s = connect()
    model = getattr(s.fold, engine)
    seq = fetch_sequences()["cyp3a4"]
    out = OUT_ROOT / tag / engine
    out.mkdir(parents=True, exist_ok=True)
    msa = None if single_sequence else get_msa(s)

    jobs = _jobs()
    key = f"{tag}/{engine}"
    jobs.setdefault(key, {"engine": engine, "tag": tag, "samples": samples,
                          "batches": []})
    # Resume on two independent markers: a ligand is skipped if it is already collected
    # to disk OR already sitting in a submitted batch. Only the first would resubmit the
    # entire in-flight campaign on the next call, which is how you get duplicate work.
    # Resume is keyed on (replicate, ligand): a ligand is "done" for replicate 3 only if
    # replicate 3 was submitted, so re-running tops the pool up instead of either
    # resubmitting everything or refusing to add depth.
    claimed = {(b.get("rep", 0), sid)
               for b in jobs[key]["batches"] for sid in b["ligands"]}
    n_sub = 0
    for rep in range(replicates):
        todo = [(r.id, r.smiles) for r in df.itertuples()
                if (rep, r.id) not in claimed]
        if not todo:
            continue
        print(f"{engine} rep {rep}: {len(todo)} ligands, {samples} samples, "
              f"{batch} per job", flush=True)
        for i in range(0, len(todo), batch):
            chunk = todo[i:i + batch]
            sids = [sid for sid, _ in chunk]
            try:
                fut = model.fold(
                    sequences=[build_complex(seq, smi, msa) for _s, smi in chunk],
                    diffusion_samples=samples, num_recycles=3)
                jobs[key]["batches"].append(
                    {"job_id": str(fut.job_id), "ligands": sids, "rep": rep,
                     "samples": samples, "submitted": time.time()})
                _save_jobs(jobs)
                n_sub += 1
                print(f"  r{rep}[{i//batch+1}] {','.join(sids)} -> {fut.job_id}",
                      flush=True)
            except Exception as exc:
                print(f"  r{rep}[{i//batch+1}] SUBMIT-FAIL {sids}: "
                      f"{type(exc).__name__}: {exc}", flush=True)
            time.sleep(0.6)
    return {"key": key, "submitted_jobs": n_sub,
            "total_batches": len(jobs[key]["batches"])}


def _split_models(cif_text: str, sid: str, out: Path, rep: int = 0) -> list[str]:
    """One file per diffusion sample, because the pose scorer takes one pose at a time.

    The samples come back as models inside a single mmCIF. Reading that with the pose
    loader would silently merge every sample into one cloud of ligand atoms, which is
    the same class of error as the heme-parsed-by-name bug: it does not raise, it just
    scores nonsense.
    """
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


def collect(engine: str, tag: str) -> dict:
    """Fetch whatever has finished. Safe to call repeatedly; skips what is already saved.

    Also records the per-sample confidence block (ranking_score, ptm, iptm, plddt, gpde,
    has_clash, disorder). That is not bookkeeping: FINDING 001 showed Boltz's own
    confidence ranks poses WORSE than random within a ligand, and whether a different
    engine's confidence does any better is a question we can only answer if we keep it.
    """
    s = connect()
    out = OUT_ROOT / tag / engine
    out.mkdir(parents=True, exist_ok=True)
    jobs = _jobs()
    key = f"{tag}/{engine}"
    if key not in jobs:
        return {"error": f"no submitted jobs for {key}"}

    n_ok = n_pending = n_fail = 0
    conf_rows = []
    for b in jobs[key]["batches"]:
        if b.get("done") is True:
            n_ok += len(b["ligands"])
            continue
        if b.get("done") == "failed":
            n_fail += len(b["ligands"])
            continue
        try:
            fut = s.load_job(b["job_id"])
            status = str(fut.job.status).upper()
            if "SUCCESS" not in status:
                if "FAIL" in status or "CANCEL" in status:
                    b["done"] = "failed"
                    n_fail += len(b["ligands"])
                else:
                    n_pending += len(b["ligands"])
                continue
            results = fut.get()
            try:
                confs = fut.get_confidence()
            except Exception:
                confs = [None] * len(results)
        except Exception as exc:
            print(f"  batch {b['job_id'][:8]}: {type(exc).__name__}", flush=True)
            n_pending += len(b["ligands"])
            continue

        ok_all = True
        for idx, sid in enumerate(b["ligands"]):
            try:
                rep = b.get("rep", 0)
                txt = results[idx].to_string()
                names = _split_models(txt, sid, out, rep)
                for k, c in enumerate(confs[idx] or []):
                    conf_rows.append({"ligand": sid, "sample": k, "rep": rep,
                                      "engine": engine,
                                      **{f: getattr(c, f) for f in
                                         ("ranking_score", "ptm", "iptm", "plddt",
                                          "gpde", "has_clash", "disorder")
                                         if hasattr(c, f)}})
                mf = out / f"{sid}.json"
                prev = json.loads(mf.read_text())["files"] if mf.exists() else []
                mf.write_text(json.dumps(
                    {"ligand": sid, "engine": engine,
                     "files": sorted(set(prev) | set(names))}, indent=1))
                n_ok += 1
            except Exception as exc:
                print(f"  {sid}: save failed ({type(exc).__name__}: {exc})", flush=True)
                n_fail += 1
                ok_all = False
        # Only retire the batch once every ligand in it is on disk. Marking it done
        # regardless would discard a whole job's poses on any transient save error - and
        # did exactly that on the first run, where a wrong gemmi function name failed all
        # four saves while the batch was retired as complete.
        if ok_all:
            b["done"] = True
    _save_jobs(jobs)

    if conf_rows:
        cf = OUT_ROOT / f"confidence_{tag}_{engine}.csv"
        prev = pd.read_csv(cf) if cf.exists() else None
        new = pd.DataFrame(conf_rows)
        out_df = pd.concat([prev, new]).drop_duplicates(["ligand", "sample", "rep"]) \
            if prev is not None else new
        out_df.to_csv(cf, index=False)

    return {"engine": engine, "tag": tag, "collected": n_ok,
            "pending": n_pending, "failed": n_fail}


def score_and_correlate(engine: str, tag: str) -> dict:
    """Score the new pool and test the ONE thing that decides whether to scale it."""
    import gemmi
    import numpy as np
    from scipy import stats

    from cypstruct import pose as P
    from cypstruct.targets import fetch_cif

    out = OUT_ROOT / tag / engine
    lig = pd.read_csv(DATA_PROCESSED / "validation_ligands.csv")
    meta = {r.id: r for r in lig.itertuples()}
    rows = []
    for manifest in sorted(out.glob("*.json")):
        sid = manifest.stem
        m = meta.get(sid)
        if m is None:
            continue
        try:
            cif = fetch_cif(m.pdb)
            st = gemmi.read_structure(str(cif))
            st.setup_entities()
            chain = next((c.name for c in st[0]
                          if any(r.name.strip().upper() == sid for r in c)), None)
            if chain is None:
                continue
            ref = P.load_structure(cif, ligand_code=sid, assembly_chain=chain)
        except Exception:
            continue
        for f in sorted(out.glob(f"{sid}__r*.cif")):
            try:
                model = P.load_structure(f)
            except Exception:
                continue
            if len(model.lig_xyz) == 0:
                continue
            perm = P.best_ligand_mapping(m.smiles, model, ref)
            rows.append({"ligand": sid, "sample": f.stem,
                         "lddt_pli": P.lddt_pli(model, ref, lig_perm=perm),
                         "bisy_rmsd": P.bisy_rmsd(model, ref, lig_perm=perm)})
    df = pd.DataFrame(rows)
    if df.empty:
        return {"engine": engine, "error": "no scoreable poses"}
    df.to_csv(DATA_PROCESSED / f"poses_scored_op_{engine}.csv", index=False)

    new_oracle = df.groupby("ligand").lddt_pli.max()
    b = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    b = b[b.arm == "unsteered"]
    ref_oracle = b.groupby("ligand").lddt_pli.max()
    j = pd.concat([new_oracle.rename("new"), ref_oracle.rename("boltz")],
                  axis=1).dropna()
    rho = stats.spearmanr(j.new, j.boltz) if len(j) > 5 else None
    union = j.max(axis=1).mean()
    return {
        "engine": engine, "n_ligands": int(len(j)), "n_poses": int(len(df)),
        "oracle": round(float(j.new.mean()), 4),
        "boltz_oracle_same_ligands": round(float(j.boltz.mean()), 4),
        "union_oracle": round(float(union), 4),
        "union_gain": round(float(union - j.boltz.mean()), 4),
        "n_where_new_beats_boltz": int((j.new > j.boltz).sum()),
        "oracle_rho_vs_boltz": (round(float(rho.statistic), 3) if rho else None),
        "oracle_rho_p": (round(float(rho.pvalue), 5) if rho else None),
        "verdict": ("DECORRELATED - worth scaling"
                    if rho and rho.statistic < 0.25 else
                    "correlated like Chai - expect little (see FINDING 005)"),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["models", "submit", "collect", "score"])
    # protenix_v2 is the default because it is the only engine measured to WORK here:
    # boltz-1, boltz-1x, boltz-2 and rosettafold-3 all fail at runtime on this complex,
    # and alphafold2 warns that it discards ligand chains outright.
    ap.add_argument("--engine", default="protenix_v2")
    ap.add_argument("--n", type=int, default=0, help="0 = all ligands")
    ap.add_argument("--samples", type=int, default=20)
    ap.add_argument("--batch", type=int, default=4, help="complexes per job")
    ap.add_argument("--replicates", type=int, default=1,
                    help="separate jobs per ligand - THIS is what samples the ligand")
    ap.add_argument("--single-sequence", action="store_true",
                    help="no MSA; REQUIRED for rosettafold_3, which fails on an uploaded one")
    ap.add_argument("--tag", default="op1")
    a = ap.parse_args()

    if a.cmd == "models":
        print(connect().fold.list_models())
    elif a.cmd == "submit":
        print(json.dumps(submit(a.engine, ligand_set(a.n or None), a.samples, a.tag,
                                batch=a.batch, replicates=a.replicates,
                                single_sequence=a.single_sequence), indent=1)[:800])
    elif a.cmd == "collect":
        print(json.dumps(collect(a.engine, a.tag), indent=2))
    else:
        print(json.dumps(score_and_correlate(a.engine, a.tag), indent=2))
