"""Co-fold the whole P450 universe: 493 pairs, 367 ligands, 185 targets.

This is the campaign that lifts n. Every selection experiment so far has run on 87 CYP3A4
ligands, where FINDING 007 puts the noise floor at +0.0138 and therefore makes a genuine
+0.015 effect unprovable. 493 pairs with crystal ground truth is 5.7x that, and - more
valuable than the raw count - it spans 185 different P450s, which converts
leave-one-ligand-cluster-out into leave-one-TARGET-out. A selection term that survives
being tested on a P450 it has never seen is a different class of evidence from one tuned
within a single pocket.

Three stages, each resumable, because an MSA search takes >10 minutes and a fold takes
minutes more. Nothing here blocks:

    python scripts/cofold/p450_campaign.py msa       # submit 185 homology searches
    python scripts/cofold/p450_campaign.py msa-status
    python scripts/cofold/p450_campaign.py submit    # fold every pair whose MSA is ready
    python scripts/cofold/p450_campaign.py collect

**Why one MSA per target and not one per pair.** 493 pairs share 185 sequences; the
alignment depends only on the sequence. Searching per pair would be 2.7x the work for
identical results.

**Why the deposited construct sequence and not the UniProt canonical.** Crystallised P450s
are truncated, mutated and tagged. Folding the wild-type sequence against a pose solved
with an engineered construct compares two different proteins, and the residue numbering
that `align_by_residue` depends on would not line up.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "cofold"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402
from openprotein_cofold import _split_models, connect  # noqa: E402

SET = DATA_PROCESSED / "p450_universe" / "p450_cofold_set.csv"
OUT = DATA_PROCESSED / "p450_universe" / "poses"
STATE = DATA_PROCESSED / "p450_universe" / "campaign.json"


def _state() -> dict:
    return json.loads(STATE.read_text()) if STATE.exists() else {"msa": {}, "folds": {}}


def _save(d: dict) -> None:
    """Write the campaign state atomically, keeping the previous copy.

    This file is the ONLY record of in-flight OpenProtein jobs - at the time this was
    written, 2,300 of them - and a submitter and two collectors all read-modify-write it
    concurrently. A bare `write_text` truncates the file before it refills, so any crash
    or interleave in that window leaves unparseable JSON and every queued job becomes
    unrecoverable: the job ids exist only here.

    Temp-file-and-replace makes the swap atomic, so a reader sees either the old state or
    the new one and never a half-written one. `.prev` is the recovery copy for the case
    atomicity does not cover - a concurrent writer clobbering another's new records.
    """
    STATE.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(d, indent=1)
    if STATE.exists():
        try:
            shutil.copyfile(STATE, STATE.with_suffix(".json.prev"))
        except OSError:
            pass
    tmp = STATE.with_suffix(f".json.tmp{os.getpid()}")
    tmp.write_text(payload)
    os.replace(tmp, STATE)


def _save_merged(d: dict) -> None:
    """Save by overlaying our changes onto whatever is on disk NOW.

    Atomicity was not enough, and that showed up as 349 vanished job records. A collector
    reads the state, spends many minutes downloading, and saves - writing back a dict that
    predates every job the submitter queued in the meantime. The newest records are the
    ones destroyed, which is the worst possible selection: those jobs are running on
    OpenProtein with nothing left that knows their ids.

    Re-reading at save time and overlaying makes concurrent writers additive. Disk is the
    base, so records this process never saw survive; our own fields win per record, so the
    `done` flags we just earned are applied. Two writers can still race inside the
    read-write window here, but it is microseconds wide instead of minutes.
    """
    try:
        base = json.loads(STATE.read_text()) if STATE.exists() else {}
    except (OSError, json.JSONDecodeError):
        _save(d)
        return
    for section in ("folds", "msa"):
        merged = dict(base.get(section, {}))
        for k, v in d.get(section, {}).items():
            if isinstance(v, dict) and isinstance(merged.get(k), dict):
                merged[k] = {**merged[k], **v}
            else:
                merged[k] = v
        d[section] = merged
    for k, v in base.items():
        d.setdefault(k, v)
    _save(d)


def targets() -> pd.DataFrame:
    """One row per distinct sequence: what an MSA is actually keyed on."""
    df = pd.read_csv(SET)
    return df.drop_duplicates("target_key")[["target_key", "sequence", "uniprot"]]


def cmd_msa(limit: int | None) -> dict:
    s = connect()
    st = _state()
    tg = targets()
    todo = [r for r in tg.itertuples() if r.target_key not in st["msa"]]
    if limit:
        todo = todo[:limit]
    print(f"{len(todo)} MSA searches to submit ({len(st['msa'])} already)", flush=True)
    for i, r in enumerate(todo):
        try:
            f = s.align.create_msa(r.sequence.encode())
            st["msa"][r.target_key] = {"job_id": str(f.job_id), "uniprot": r.uniprot,
                                       "submitted": time.time()}
            _save(st)
            print(f"  [{i+1}/{len(todo)}] {r.target_key} ({r.uniprot}) -> {f.job_id}",
                  flush=True)
        except Exception as exc:
            print(f"  [{i+1}/{len(todo)}] {r.target_key} FAIL "
                  f"{type(exc).__name__}: {exc}", flush=True)
        time.sleep(0.5)
    return {"submitted": len(todo), "total": len(st["msa"])}


def cmd_msa_status() -> dict:
    s = connect()
    st = _state()
    counts: dict[str, int] = {}
    for key, rec in st["msa"].items():
        if rec.get("status") == "SUCCESS":
            counts["SUCCESS"] = counts.get("SUCCESS", 0) + 1
            continue
        try:
            f = s.load_job(rec["job_id"])
            status = str(f.job.status).split(".")[-1]
        except Exception as exc:
            status = f"ERR_{type(exc).__name__}"
        rec["status"] = status
        counts[status] = counts.get(status, 0) + 1
    _save(st)
    return counts


def cmd_submit(samples: int, batch: int, limit: int | None,
               replicates: int = 1, engine: str = "protenix_v2",
               single_sequence: bool = False) -> dict:
    """Fold every pair whose target MSA is ready, grouped so a job shares one target.

    `replicates`, not `samples`, is what builds a pool - see FINDING 009. Within one job
    every model carries an identical ligand and only the protein moves, so a run with
    samples=10 and replicates=1 yields exactly one pose per pair.
    """
    from openprotein_cofold import build_complex

    s = connect()
    st = _state()
    df = pd.read_csv(SET)
    if single_sequence:
        # No MSA needed, so EVERY target is available immediately. The alignment queue is
        # serial at roughly one per 20 minutes, so waiting for all 185 costs about two
        # days; single-sequence mode trades some prediction quality for all of it back.
        ready = set(df.target_key)
        print(f"single-sequence mode: all {len(ready)} targets available", flush=True)
    else:
        ready = {k for k, v in st["msa"].items() if v.get("status") == "SUCCESS"}
        print(f"{len(ready)} of {df.target_key.nunique()} targets have an MSA", flush=True)

    df["pair"] = df.pdb + "_" + df.id
    # Some ligands make an engine fail server-side every single time - esmfold2 cannot
    # handle the Ir/Ru organometallics, and each replicate re-failed on the same pairs,
    # burning queue slots for nothing. Skip what is already known not to work here.
    skip_path = (DATA_PROCESSED / "p450_universe" / f"{engine}_unsupported.json")
    if skip_path.exists():
        skip = set(json.loads(skip_path.read_text()).get("pairs", []))
        before = len(df)
        df = df[~df.pair.isin(skip)]
        if before != len(df):
            print(f"skipping {before - len(df)} pairs {engine} cannot fold", flush=True)
    n = 0
    # Group by target: every complex in one job must carry the same MSA future, and
    # grouping also means one MSA object is reused rather than reloaded per complex.
    for rep in range(replicates):
      claimed = {p for f in st["folds"].values()
                 if f.get("rep", 0) == rep
                 and f.get("engine", "protenix_v2") == engine
                 and f.get("single_sequence", False) == single_sequence
                 for p in f["pairs"]}
      todo = df[df.target_key.isin(ready) & ~df.pair.isin(claimed)]
      print(f"rep {rep}: {len(todo)} pairs to fold", flush=True)
      for key, grp in todo.groupby("target_key"):
        if limit and n >= limit:
            break
        msa = None if single_sequence else s.load_job(st["msa"][key]["job_id"])
        seq = grp.sequence.iloc[0]
        rows = list(grp.itertuples())
        for i in range(0, len(rows), batch):
            if limit and n >= limit:
                break
            chunk = rows[i:i + batch]
            try:
                fut = getattr(s.fold, engine).fold(
                    sequences=[build_complex(seq, r.smiles, msa) for r in chunk],
                    diffusion_samples=samples, num_recycles=3)
                st["folds"][str(fut.job_id)] = {
                    "target_key": key, "pairs": [r.pair for r in chunk],
                    "ligands": [r.id for r in chunk], "samples": samples,
                    "rep": rep, "engine": engine, "single_sequence": single_sequence,
                    "submitted": time.time()}
                _save(st)
                n += 1
                print(f"  {key} [{','.join(r.id for r in chunk)}] -> {fut.job_id}",
                      flush=True)
            except Exception as exc:
                print(f"  {key} FAIL {type(exc).__name__}: {str(exc)[:160]}", flush=True)
            time.sleep(0.5)
    return {"jobs_submitted": n, "total_jobs": len(st["folds"])}


def cmd_collect() -> dict:
    s = connect()
    st = _state()
    OUT.mkdir(parents=True, exist_ok=True)
    n_ok = n_pend = n_fail = 0
    for job_id, rec in st["folds"].items():
        if rec.get("done") is True:
            n_ok += len(rec["pairs"])
            continue
        try:
            f = s.load_job(job_id)
            status = str(f.job.status).upper()
            if "SUCCESS" not in status:
                if "FAIL" in status or "CANCEL" in status:
                    rec["done"] = "failed"
                    n_fail += len(rec["pairs"])
                else:
                    n_pend += len(rec["pairs"])
                continue
            res = f.get()
        except Exception:
            n_pend += len(rec["pairs"])
            continue

        ok_all = True
        for idx, pair in enumerate(rec["pairs"]):
            rep = rec.get("rep", 0)
            eng = rec.get("engine", "protenix_v2")
            if rec.get("single_sequence"):
                eng = eng + "_ss"      # kept apart: different inputs, different pool
            # per-engine subdirectory: the cross-engine feature needs to know WHICH
            # engine produced a pose, and a flat directory silently merges them
            d = OUT / pair / eng
            mf = d / "manifest.json"
            if any(d.glob(f"{pair}__r{rep}s*.cif")):
                n_ok += 1
                continue
            try:
                d.mkdir(parents=True, exist_ok=True)
                names = _split_models(res[idx].to_string(), pair, d, rep)
                prev = json.loads(mf.read_text())["files"] if mf.exists() else []
                mf.write_text(json.dumps(
                    {"pair": pair, "files": sorted(set(prev) | set(names)),
                     "engine": eng}, indent=1))
                n_ok += 1
            except Exception as exc:
                print(f"  {pair}: {type(exc).__name__}: {exc}", flush=True)
                n_fail += 1
                ok_all = False
        if ok_all:
            rec["done"] = True
    # a collect pass takes minutes; the submitter appends records throughout it, so
    # writing our own stale dict back would delete them (it did - 349 of them)
    _save_merged(st)
    return {"collected": n_ok, "pending": n_pend, "failed": n_fail}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["msa", "msa-status", "submit", "collect", "status"])
    ap.add_argument("--samples", type=int, default=10)
    ap.add_argument("--batch", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--replicates", type=int, default=1,
                    help="separate jobs per pair - THIS is what samples the ligand")
    ap.add_argument("--single-sequence", action="store_true",
                    help="no MSA; unblocks every target at once and works on all engines")
    ap.add_argument("--engine", default="protenix_v2",
                    help="protenix_v2 | protenix | esmfold2 (all verified to run here)")
    a = ap.parse_args()
    lim = a.limit or None

    if a.cmd == "msa":
        print(json.dumps(cmd_msa(lim), indent=2))
    elif a.cmd == "msa-status":
        print(json.dumps(cmd_msa_status(), indent=2))
    elif a.cmd == "submit":
        print(json.dumps(cmd_submit(a.samples, a.batch, lim, a.replicates, a.engine,
                                    a.single_sequence), indent=2))
    elif a.cmd == "collect":
        print(json.dumps(cmd_collect(), indent=2))
    else:
        st = _state()
        print(json.dumps({"msa": len(st["msa"]), "fold_jobs": len(st["folds"])}, indent=2))
