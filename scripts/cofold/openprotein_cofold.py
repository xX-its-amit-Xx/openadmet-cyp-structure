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
PASS = os.environ.get("OPENPROTEIN_PASS", "REDACTED")

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


def ligand_set(limit: int | None, seed: int = 0, csv: str | None = None) -> pd.DataFrame:
    df = pd.read_csv(csv or (DATA_PROCESSED / "validation_ligands.csv"))
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


def parse_sweep(spec: str) -> list[tuple[int, int]]:
    """`'3x200,10x200,3x50,3x400'` -> [(3,200),(10,200),(3,50),(3,400)].

    Refuses `num_recycles < 2`: FINDING 016 measured `num_recycles=1` at **-0.0596**
    against the default over 100 pairs - the one setting that degrades rather than
    diversifies. A reference pose that is simply worse is not an independent opinion.
    """
    out: list[tuple[int, int]] = []
    for tok in spec.split(","):
        tok = tok.strip().lower()
        if not tok:
            continue
        rec, steps = tok.split("x")
        rec, steps = int(rec), int(steps)
        if rec < 2:
            raise SystemExit(
                f"refusing sweep setting {tok!r}: num_recycles < 2 measured -0.0596 "
                "against the default; it degrades rather than diversifies (FINDING 016)")
        out.append((rec, steps))
    if len(out) != len(set(out)):
        raise SystemExit(f"duplicate settings in --sweep {spec!r}: each wave must differ")
    return out


def submit(engine: str, df: pd.DataFrame, samples: int, tag: str,
           batch: int = 4, replicates: int = 1, single_sequence: bool = False,
           sweep: list[tuple[int, int]] | None = None) -> dict:
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

    **`--replicates` is now DEAD and `--sweep` replaces it (FINDINGs 015/016/034).** Both
    OpenProtein engines became deterministic per configuration: a 12->24 replicate doubling
    moved the oracle on 0 of 489 pairs and replicates 16-29 returned 196 of 196
    byte-identical. What still varies the output is the SAMPLER, and distinct poses track
    submission WAVES rather than jobs. `sweep=[(3,200),(10,200),(3,50),(3,400)]` submits
    one wave per setting, each as its own set of jobs, and records the setting alongside
    the wave index so `refset` can tell four opinions from four copies. Four settings gave
    4 distinct placements in FINDING 016 - which is exactly the depth-4 minimum the
    cross-engine feature refuses below.
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
    # One wave per sampler setting (the live depth lever), else `replicates` waves at the
    # default setting (kept so existing call sites and reruns behave exactly as before).
    waves = ([(i, rec, steps) for i, (rec, steps) in enumerate(sweep)] if sweep
             else [(rep, 3, 200) for rep in range(replicates)])
    # A wave index already used at a DIFFERENT setting would be indistinguishable from a
    # replicate of the first one downstream, since the reference is keyed on (engine, rep).
    prior = {b["rep"]: b.get("setting") for b in jobs[key]["batches"] if "rep" in b}
    n_sub = 0
    for rep, rec, steps in waves:
        setting = f"{rec}x{steps}"
        if prior.get(rep) not in (None, setting):
            raise SystemExit(
                f"wave {rep} of {key} was already submitted at setting {prior[rep]!r}, "
                f"not {setting!r}. Waves are the unit of independence - pick a free "
                f"index (used: {sorted(k for k in prior if prior[k])}) rather than "
                "reusing one, or the reference set will count two settings as one.")
        todo = [(r.id, r.smiles) for r in df.itertuples()
                if (rep, r.id) not in claimed]
        if not todo:
            continue
        print(f"{engine} wave {rep} ({setting}): {len(todo)} ligands, {samples} samples, "
              f"{batch} per job", flush=True)
        for i in range(0, len(todo), batch):
            chunk = todo[i:i + batch]
            sids = [sid for sid, _ in chunk]
            try:
                fut = model.fold(
                    sequences=[build_complex(seq, smi, msa) for _s, smi in chunk],
                    diffusion_samples=samples, num_recycles=rec, num_steps=steps)
                jobs[key]["batches"].append(
                    {"job_id": str(fut.job_id), "ligands": sids, "rep": rep,
                     "setting": setting, "samples": samples, "submitted": time.time()})
                _save_jobs(jobs)
                n_sub += 1
                print(f"  w{rep}[{i//batch+1}] {','.join(sids)} -> {fut.job_id}",
                      flush=True)
            except Exception as exc:
                print(f"  w{rep}[{i//batch+1}] SUBMIT-FAIL {sids}: "
                      f"{type(exc).__name__}: {exc}", flush=True)
            time.sleep(0.6)
    return {"key": key, "submitted_jobs": n_sub, "waves": [w[0] for w in waves],
            "settings": {str(w[0]): f"{w[1]}x{w[2]}" for w in waves},
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


def refset(tag: str, engines: list, out_npz: Path, min_depth: int = 4,
           ligands_csv: str | None = None, tol: float = 0.05,
           allow_thin: bool = False) -> dict:
    """Turn collected OpenProtein waves into the frozen reference set the selector needs.

    **This is the drop-day half of G6, and the loud failure is the point of it.** The
    shipped selector scores a Boltz pose by its mean Chamfer distance to independent-engine
    poses of the SAME ligand; `reference_set_cyp3a4.npz` covers only the 87 validation
    ligands, so a blind set has NO reference until this runs. Below 4 independent poses the
    feature measured **-0.0055** - it makes selection worse, not weaker - and the only
    symptom downstream is `build_submission` quietly choosing the FINDING 003 rule at
    +0.0265 instead of +0.0395. So this REFUSES, by default, rather than writing a thin set.

    Three filters, each counted in the report because a count is the only thing that shows
    a filter fired (T7/T8):

      1. **one file per (engine, wave)** - within a single job every diffusion sample shares
         one ligand conformation (FINDING 009), so N models are 1 opinion, not N;
      2. **md5 of the file bytes** - replicates of a deterministic engine come back
         byte-identical (FINDING 015: 196 of 196), and historical pools are ~85% duplicates;
      3. **coordinates in the heme frame at `tol`** - two files can differ in a header and
         be the same pose. Bytes prove identity; they do not prove difference.
    """
    import hashlib
    import re as _re

    sys.path.insert(0, str(REPO / "src"))
    from cypstruct import pose as P
    from cypstruct import xengine as X

    want = None
    if ligands_csv:
        df = pd.read_csv(ligands_csv)
        idc = "structure" if "structure" in df.columns else "id"
        want = set(df[idc].astype(str))

    fname = _re.compile(r"(.+)__r(\d+)s(\d+)\.cif$")
    counts = {"files_seen": 0, "unparseable_name": 0, "not_in_csv": 0,
              "dropped_same_wave": 0, "dropped_same_md5": 0,
              "dropped_same_coords": 0, "unreadable": 0, "no_heme_frame": 0}
    picked: dict = {}
    for eng in engines:
        d = OUT_ROOT / tag / eng
        for f in sorted(d.glob("*__r*.cif")):
            counts["files_seen"] += 1
            m = fname.match(f.name)
            if m is None:
                counts["unparseable_name"] += 1
                continue
            lig, wave = m.group(1), int(m.group(2))
            if want is not None and lig not in want:
                counts["not_in_csv"] += 1
                continue
            key = (eng, wave)
            if key in picked.get(lig, {}):
                counts["dropped_same_wave"] += 1     # FINDING 009: samples are not poses
                continue
            picked.setdefault(lig, {})[key] = f

    ref: dict = {}
    per_lig: dict = {}
    for lig, by in sorted(picked.items()):
        seen_md5 = set()
        vs = []
        for key, f in sorted(by.items()):
            h = hashlib.md5(f.read_bytes()).hexdigest()
            if h in seen_md5:
                counts["dropped_same_md5"] += 1      # FINDING 015 / 034
                continue
            seen_md5.add(h)
            try:
                v = X.in_heme_frame(P.load_structure(f))
            except Exception:
                counts["unreadable"] += 1
                continue
            if v is None:
                counts["no_heme_frame"] += 1
                continue
            vs.append(v)
        n_before = len(vs)
        vs = X._dedupe(vs, tol=tol)
        counts["dropped_same_coords"] += n_before - len(vs)
        if vs:
            ref[lig] = vs
        per_lig[lig] = {"waves": len(by), "after_md5": n_before, "depth": len(vs)}

    thin = sorted(k for k, v in per_lig.items() if v["depth"] < min_depth)
    missing = sorted(want - set(per_lig)) if want else []
    depth = X.reference_depth(ref) if ref else {"ligands": 0, "usable": False}
    report = {"tag": tag, "engines": list(engines), "filters": counts,
              "ligands_with_poses": len(per_lig), "reference_depth": depth,
              "min_depth": min_depth, "below_min_depth": thin,
              "no_poses_at_all": missing, "per_ligand": per_lig}

    if (thin or missing) and not allow_thin:
        print(json.dumps({k: v for k, v in report.items() if k != "per_ligand"},
                         indent=2, default=str))
        print("")
        print("REFUSING to write %s." % out_npz)
        print("  %d ligand(s) below depth %d: %s" % (len(thin), min_depth, thin[:10]))
        if missing:
            print("  %d ligand(s) with no reference pose at all: %s"
                  % (len(missing), missing[:10]))
        print("At one reference pose this feature measured -0.0055: a thin reference set")
        print("does not give a weaker selector, it gives a HARMFUL one, and the only")
        print("downstream symptom is a silent fall back to FINDING 003 (+0.0265).")
        print("Buy depth with SAMPLER WAVES - `submit --sweep 3x200,10x200,3x50,3x400`,")
        print("one wave per setting (FINDING 016). NOT --replicates, which is")
        print("byte-identical on both engines (FINDING 015).")
        print("Pass --allow-thin only if you will also pass --skip-thin downstream and")
        print("REPORT which ligands fell back.")
        raise SystemExit(2)

    out_npz = Path(out_npz)
    report["saved"] = X.save_reference(ref, out_npz)
    report["out"] = str(out_npz)
    return report


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["models", "submit", "collect", "score", "refset"])
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
    ap.add_argument("--csv", default=None, help="ligand csv; defaults to validation set")
    ap.add_argument("--sweep", default=None,
                    help="sampler settings as '<recycles>x<steps>,...', one submission "
                         "WAVE each, e.g. '3x200,10x200,3x50,3x400'. This is the only "
                         "live depth lever: --replicates is byte-identical on both "
                         "engines (FINDING 015/016/034).")
    ap.add_argument("--engines", nargs="+", default=["protenix_v2", "protenix"],
                    help="refset: engine pools under data/processed/openprotein/<tag>/")
    ap.add_argument("--out", default=None, help="refset: output .npz")
    ap.add_argument("--min-depth", type=int, default=4,
                    help="refset: independent poses per ligand below which it REFUSES")
    ap.add_argument("--allow-thin", action="store_true",
                    help="refset: write anyway (then you MUST pass --skip-thin to "
                         "build_xeng_feature and report the fallbacks)")
    a = ap.parse_args()

    if a.cmd == "models":
        print(connect().fold.list_models())
    elif a.cmd == "submit":
        print(json.dumps(submit(a.engine, ligand_set(a.n or None, csv=a.csv), a.samples, a.tag,
                                batch=a.batch, replicates=a.replicates,
                                single_sequence=a.single_sequence,
                                sweep=parse_sweep(a.sweep) if a.sweep else None),
                         indent=1)[:800])
    elif a.cmd == "refset":
        out = a.out or str(DATA_PROCESSED / ("reference_set_%s.npz" % a.tag))
        r = refset(a.tag, a.engines, Path(out), min_depth=a.min_depth,
                   ligands_csv=a.csv, allow_thin=a.allow_thin)
        print(json.dumps({k: v for k, v in r.items() if k != "per_ligand"},
                         indent=2, default=str))
    elif a.cmd == "collect":
        print(json.dumps(collect(a.engine, a.tag), indent=2))
    else:
        print(json.dumps(score_and_correlate(a.engine, a.tag), indent=2))
