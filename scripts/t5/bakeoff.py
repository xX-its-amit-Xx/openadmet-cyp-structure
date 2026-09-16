"""T5 acceptance gate: the 200-compound co-folding bake-off, run on OpenProtein.

The BU team's T5 spec (tasks/T5-structure-generation.md) makes this the non-negotiable
gate before any production run: 200 compounds spanning the pIC50 range on CYP3A4 and
CYP2D6, measuring pose validity, completion rate, and whether derived features reduce
downstream RMSE. It also says GPU access requires an SCC allocation.

That last part is what this script disproves. OpenProtein folds protein+HEM+ligand on six
engines, unmetered, with no local GPU — so the co-folding arm of the gate can run now.
See https://github.com/buaiml/openadmet-cyp/issues/4.

**Completion is the measurement, not an obstacle.** T5 is explicit that co-folders "fail
entirely on flexible pockets rather than producing poor poses", which is why its output
schema carries a `completed` flag. So a failed job here is a data point and is recorded
as one; it is not retried into silence.

Sampling spans the pIC50 range by decile rather than taking the strongest binders,
because the question is whether structure helps across the range - a bake-off on potent
compounds only would answer a different and easier question.

    python scripts/t5/bakeoff.py select          # choose + freeze the compound list
    python scripts/t5/bakeoff.py msa             # build the two isoform MSAs
    python scripts/t5/bakeoff.py submit          # fold them
    python scripts/t5/bakeoff.py collect         # emit the T5 schema csv
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

SIBLING = Path(r"D:\Users\ashenoy00000\.windsurf\OpenADMET-cyp-challenge")
TRAIN = SIBLING / "data" / "raw" / "cyp-challenge-TRAIN_inhibition.csv"

OUT = DATA_PROCESSED / "t5_bakeoff"
STATE = OUT / "state.json"
POSES = OUT / "poses"

# T5 names CYP3A4 and CYP2D6. Axial cysteines from cypstruct.targets.
ISOFORMS = {
    "CYP3A4": {"uniprot": "P08684", "col": "CYP3A4_pIC50_direct_inhibition"},
    "CYP2D6": {"uniprot": "P10635", "col": "CYP2D6_pIC50_direct_inhibition"},
}
ENGINE = "protenix_v2"


def _state() -> dict:
    if STATE.exists():
        return json.loads(STATE.read_text())
    return {"msa": {}, "folds": {}, "compounds": None}


def _save(d: dict) -> None:
    """Atomic + merge, for the reason p450_campaign learned the hard way.

    A collector that reads state, spends minutes downloading, then writes its own stale
    dict back deletes every job the submitter queued meanwhile. That destroyed 349 job
    records once; the ids existed nowhere else.
    """
    STATE.parent.mkdir(parents=True, exist_ok=True)
    try:
        base = json.loads(STATE.read_text()) if STATE.exists() else {}
    except (OSError, json.JSONDecodeError):
        base = {}
    for section in ("msa", "folds"):
        merged = dict(base.get(section, {}))
        for k, v in d.get(section, {}).items():
            merged[k] = {**merged[k], **v} if isinstance(v, dict) and isinstance(
                merged.get(k), dict) else v
        d[section] = merged
    if STATE.exists():
        try:
            shutil.copyfile(STATE, STATE.with_suffix(".json.prev"))
        except OSError:
            pass
    tmp = STATE.with_suffix(f".json.tmp{os.getpid()}")
    tmp.write_text(json.dumps(d, indent=1))
    os.replace(tmp, STATE)


def inchikey_block(smiles: str) -> str:
    """T5 joins on `inchikey_block` - the first 14 chars, i.e. the skeleton layer."""
    from rdkit import Chem
    from rdkit.Chem import inchi
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return ""
    try:
        return inchi.MolToInchiKey(m).split("-")[0]
    except Exception:
        return ""


def cmd_select(n: int, seed: int) -> dict:
    """Pick n compounds spanning the pIC50 range on BOTH isoforms."""
    df = pd.read_csv(TRAIN)
    c3, c2 = ISOFORMS["CYP3A4"]["col"], ISOFORMS["CYP2D6"]["col"]
    both = df.dropna(subset=[c3, c2]).copy()
    print(f"{len(df)} train rows; {len(both)} measured on BOTH 3A4 and 2D6")

    # Stratify on the 3A4 decile so the sample spans potency rather than concentrating
    # where the data is dense. T5 asks for "spanning the pIC50 range" explicitly.
    both["bin"] = pd.qcut(both[c3], q=10, labels=False, duplicates="drop")
    per = max(1, n // both["bin"].nunique())
    picked = (both.groupby("bin", group_keys=False)
                  .apply(lambda g: g.sample(min(per, len(g)), random_state=seed)))
    if len(picked) < n:                      # top up from the remainder, still at random
        rest = both.drop(picked.index)
        picked = pd.concat([picked, rest.sample(min(n - len(picked), len(rest)),
                                                random_state=seed)])
    picked = picked.head(n).reset_index(drop=True)
    picked["inchikey_block"] = [inchikey_block(s) for s in picked.SMILES]
    bad = (picked.inchikey_block == "").sum()
    if bad:
        print(f"  WARNING {bad} compounds gave no InChIKey - dropped")
        picked = picked[picked.inchikey_block != ""].reset_index(drop=True)

    OUT.mkdir(parents=True, exist_ok=True)
    keep = ["Molecule_Name", "SMILES", "inchikey_block", c3, c2]
    picked[keep].to_csv(OUT / "compounds.csv", index=False)
    st = _state()
    st["compounds"] = str(OUT / "compounds.csv")
    _save(st)
    return {
        "compounds": len(picked),
        "pIC50_3A4": [round(float(picked[c3].min()), 2), round(float(picked[c3].max()), 2)],
        "pIC50_2D6": [round(float(picked[c2].min()), 2), round(float(picked[c2].max()), 2)],
        "deciles_covered": int(picked["bin"].nunique()) if "bin" in picked else None,
    }


def cmd_msa() -> dict:
    """One MSA per isoform. This is the expensive shared input; build it once."""
    from openprotein_cofold import connect

    from cypstruct import targets as T
    s = connect()
    seqs = T.fetch_sequences()
    st = _state()
    for iso in ISOFORMS:
        key = iso.lower()
        if st["msa"].get(iso, {}).get("job_id"):
            continue
        seq = seqs.get(key)
        if not seq:
            print(f"  {iso}: no sequence")
            continue
        fut = s.align.create_msa(seq.encode())
        st["msa"][iso] = {"job_id": str(fut.job_id), "len": len(seq)}
        print(f"  {iso}: MSA {fut.job_id} ({len(seq)} aa)", flush=True)
    _save(st)
    return {iso: st["msa"].get(iso, {}) for iso in ISOFORMS}


def cmd_msa_status() -> dict:
    from openprotein_cofold import connect
    s = connect()
    st = _state()
    out = {}
    for iso, rec in st["msa"].items():
        try:
            j = s.load_job(rec["job_id"])
            rec["status"] = str(j.job.status).upper()
        except Exception as exc:
            rec["status"] = f"ERR {type(exc).__name__}"
        out[iso] = rec["status"]
    _save(st)
    return out


def cmd_submit(limit: int | None) -> dict:
    from openprotein_cofold import build_complex, connect
    s = connect()
    st = _state()
    df = pd.read_csv(st["compounds"])
    n = 0
    for iso in ISOFORMS:
        rec = st["msa"].get(iso, {})
        if str(rec.get("status", "")).upper() != "SUCCESS":
            print(f"  {iso}: MSA not ready ({rec.get('status')}) - skipping")
            continue
        msa = s.load_job(rec["job_id"])
        claimed = {f["key"] for f in st["folds"].values()}
        for r in df.itertuples():
            key = f"{iso}:{r.inchikey_block}"
            if key in claimed:
                continue
            if limit and n >= limit:
                break
            try:
                fut = getattr(s.fold, ENGINE).fold(
                    sequences=[build_complex(seqs_for(iso), r.SMILES, msa)],
                    diffusion_samples=1, num_recycles=3)
            except Exception as exc:
                # A submission that never leaves the client is a completion failure too,
                # and T5's schema wants it recorded rather than dropped.
                st["folds"][f"unsent:{key}:{n}"] = {
                    "key": key, "isoform": iso, "inchikey_block": r.inchikey_block,
                    "smiles": r.SMILES, "done": "failed",
                    "error": f"{type(exc).__name__}: {exc}"[:200]}
                _save(st)
                n += 1
                continue
            st["folds"][str(fut.job_id)] = {
                "key": key, "isoform": iso, "inchikey_block": r.inchikey_block,
                "smiles": r.SMILES, "submitted": time.time()}
            _save(st)
            n += 1
            if n % 25 == 0:
                print(f"  submitted {n}", flush=True)
    return {"submitted": n, "total_records": len(st["folds"])}


_SEQ_CACHE: dict[str, str] = {}


def seqs_for(iso: str) -> str:
    if not _SEQ_CACHE:
        from cypstruct import targets as T
        _SEQ_CACHE.update(T.fetch_sequences())
    return _SEQ_CACHE[iso.lower()]


def cmd_collect() -> dict:
    from openprotein_cofold import _split_models, connect
    s = connect()
    st = _state()
    POSES.mkdir(parents=True, exist_ok=True)
    ok = pend = fail = 0
    for job_id, rec in st["folds"].items():
        if rec.get("done") is True:
            ok += 1
            continue
        if rec.get("done") == "failed":
            fail += 1
            continue
        try:
            f = s.load_job(job_id)
            status = str(f.job.status).upper()
            if "FAIL" in status or "CANCEL" in status:
                rec["done"] = "failed"
                rec["error"] = f"job status {status}"
                fail += 1
                continue
            if "SUCCESS" not in status:
                pend += 1
                continue
            res = f.get()
        except Exception:
            pend += 1
            continue
        d = POSES / rec["isoform"]
        d.mkdir(parents=True, exist_ok=True)
        try:
            names = _split_models(res[0].to_string(), rec["inchikey_block"], d, "t5")
            rec["done"] = True
            rec["files"] = names
            ok += 1
        except Exception as exc:
            rec["done"] = "failed"
            rec["error"] = f"parse: {type(exc).__name__}: {exc}"[:200]
            fail += 1
    _save(st)
    return {"collected": ok, "pending": pend, "failed": fail}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["select", "msa", "msa-status", "submit", "collect"])
    ap.add_argument("-n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.cmd == "select":
        print(json.dumps(cmd_select(a.n, a.seed), indent=2))
    elif a.cmd == "msa":
        print(json.dumps(cmd_msa(), indent=2))
    elif a.cmd == "msa-status":
        print(json.dumps(cmd_msa_status(), indent=2))
    elif a.cmd == "submit":
        print(json.dumps(cmd_submit(a.limit or None), indent=2))
    else:
        print(json.dumps(cmd_collect(), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
