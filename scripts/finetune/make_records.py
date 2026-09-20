"""Convert our crystal structures into Boltz training records.

Runs ON EXPLORER, against the staged inputs. Uses `boltz.data.parse.mmcif.parse_mmcif` -
Boltz's own parser, not a reimplementation - because training data parsed differently from
how inference parses it degrades a fine-tune silently instead of failing. That class of
mismatch is the most expensive kind of bug available here.

Emits, per arm:
    <arm>/structures/<pdb>.npz   the parsed structure
    <arm>/manifest.json          records, with our train/test split carried through

The split comes from build_arms.py, which holds out whole Murcko-scaffold clusters and
checks the held-out novelty against the challenge's measured 0.587 median (FINDING 017).
Nothing here re-splits; it only carries that decision forward.

    ./env/bin/python make_records.py --arm arm4_mix
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

ROOT = Path("/scratch/shenoy.am/cyp-finetune")
CCD = ROOT / "boltz_cache" / "ccd.pkl"
# mols.tar unpacks to mols/mols/*.pkl - one level deeper than the obvious path.
# Pointing at the outer dir yields moldir=None and reproduces the very failure
# this fixes: 79 structures whose ligands are absent from ccd.pkl died with
# "expected str, bytes or os.PathLike object, not NoneType".
MOLDIR = ROOT / "boltz_cache" / "mols" / "mols"


def load_mols():
    """The CCD component dictionary the parser needs to recognise HEM and the ligands."""
    import pickle
    if CCD.exists():
        with CCD.open("rb") as fh:
            return pickle.load(fh)  # noqa: S301
    return None


def main() -> int:
    import pandas as pd
    from boltz.data.parse.mmcif import parse_mmcif

    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="arm4_mix")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    arm_csv = ROOT / "finetune_arms" / f"{a.arm}.csv"
    df = pd.read_csv(arm_csv)
    if a.limit:
        df = df.head(a.limit)
    out = ROOT / "records" / a.arm
    (out / "structures").mkdir(parents=True, exist_ok=True)

    mols = load_mols()
    moldir = str(MOLDIR) if MOLDIR.exists() else None
    print(f"{a.arm}: {len(df)} pairs | mols={'yes' if mols else 'no'} moldir={moldir}")

    records, ok, fail = [], 0, 0
    reasons: dict[str, int] = {}
    for r in df.itertuples():
        cif = ROOT / "rcsb" / f"{r.pdb}.cif"
        if not cif.exists():
            reasons["missing cif"] = reasons.get("missing cif", 0) + 1
            fail += 1
            continue
        try:
            # compute_interfaces=False is a WORKAROUND for an upstream bug, not a
            # preference. boltz/data/parse/mmcif.py defines compute_interfaces() as a
            # module function at line 319 AND takes a parameter of the same name; at
            # line 1206 the parameter shadows the function, so `compute_interfaces(
            # atoms, chains)` calls a bool and raises TypeError. The default is True,
            # so every caller using defaults hits it. False takes the else-branch and
            # returns an empty interface array.
            #
            # What that costs: `interfaces` is chain-chain contact metadata. Our
            # training signal is the LIGAND pose, carried by diffusion_loss_weight 4.0
            # on coordinates, so an empty interface list does not remove the gradient
            # we are training on. Worth revisiting if interface-weighted loss is ever
            # turned on.
            parsed = parse_mmcif(str(cif), mols=mols, moldir=moldir,
                                 use_assembly=False, compute_interfaces=False)
            npz = out / "structures" / f"{r.pdb}.npz"
            parsed.data.dump(npz) if hasattr(parsed, "data") else None
            rec = parsed.info if hasattr(parsed, "info") else None
            records.append({
                "pair": r.pair, "pdb": r.pdb, "ligand": r.id,
                "split": r.split, "npz": str(npz),
                "record": rec.to_dict() if hasattr(rec, "to_dict") else None,
            })
            ok += 1
        except Exception as exc:  # noqa: BLE001
            fail += 1
            k = f"{type(exc).__name__}: {str(exc)[:60]}"
            reasons[k] = reasons.get(k, 0) + 1
            if fail <= 3:
                traceback.print_exc(limit=2)
        if (ok + fail) % 50 == 0:
            print(f"  {ok+fail}/{len(df)}  ok={ok} fail={fail}", flush=True)

    (out / "manifest.json").write_text(json.dumps(records, indent=1, default=str))
    print(json.dumps({
        "arm": a.arm, "parsed": ok, "failed": fail,
        "train": sum(1 for x in records if x["split"] == "train"),
        "test": sum(1 for x in records if x["split"] == "test"),
        "failure_reasons": dict(sorted(reasons.items(), key=lambda kv: -kv[1])[:5]),
        "out": str(out),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
