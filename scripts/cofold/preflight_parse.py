"""Dry-parse every planned co-folding job on CPU before any GPU is allocated.

This is a standing gate, not a debugging aid. It has already paid for itself twice:

- It caught a `contact` constraint with three list elements where Boltz wants two. On the
  GPU path that failure appears as **return code 0 with zero structures in 25 seconds**,
  because `boltz predict` catches a schema error per input and moves on. Return code 0
  with no output is the shape of failure most easily mistaken for success.
- It caught **12 organometallic ligands** (iridium and ruthenium complexes with dative-bond
  SMILES) that Boltz cannot parse at all, before 21 GPU-hours were committed to a batch
  that would have dropped them.

Cost: seconds of CPU on Modal. Run it before every submit.

    python scripts/cofold/preflight_parse.py --csv <ligands.csv> --tag <tag> --samples 20
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "cofold"))

import modal_boltz as MB  # noqa: E402

from cypstruct.paths import DATA_PROCESSED  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--tag", default="preflight")
    ap.add_argument("--samples", type=int, default=20)
    ap.add_argument("--seeds", default="1")
    a = ap.parse_args()

    seeds = [int(s) for s in a.seeds.split(",")]
    jobs = MB.plan(a.csv, a.tag, samples=a.samples, seeds=seeds)
    variants = {j["job_id"]: j["yaml"] for j in jobs}
    print(f"dry-parsing {len(variants)} planned jobs on CPU ...", flush=True)

    with MB.app.run():
        res = MB.probe_yaml.remote(variants)

    sig = res.pop("_signature", None)
    bad = {k: v for k, v in res.items() if not v.get("ok")}
    ok_n = len(res) - len(bad)
    print(f"{ok_n}/{len(res)} parse cleanly   (parser signature: {sig})")
    for k, v in sorted(bad.items())[:30]:
        print(f"  FAIL {k:32s} {v.get('error', '')[:110]}")

    out = DATA_PROCESSED / f"preflight_{a.tag}.json"
    out.write_text(json.dumps(
        {"tag": a.tag, "n_jobs": len(res), "n_ok": ok_n,
         "failures": {k: v.get("error") for k, v in bad.items()}}, indent=1))
    print(f"-> {out}")

    # Non-zero exit so a caller can gate a submit on this without parsing stdout.
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
