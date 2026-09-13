"""Produce the cross-engine feature the submission builder selects on.

This is the drop-day bridge between a pool of poses and `build_submission.py`. It writes
`data/processed/xeng_<tag>.csv` with {ligand, sample, xeng}, which `choose_poses` picks up
automatically and prefers over the older selectors.

**It refuses rather than degrades.** The feature is worth +0.0381 at >= 4 genuinely
independent reference poses and **-0.0055 at one** - so a silently thin reference set does
not give a weaker selector, it gives a harmful one. `reference_depth` is checked before
anything is written and the run aborts with the numbers if the pool is too shallow.

Three traps, all of them already paid for once:

* `diffusion_samples` does not diversify the ligand on OpenProtein - within one job every
  model shares a single ligand conformation, for Protenix and esmfold2 alike (FINDING 009).
  Only replicate jobs count, which is why references are keyed on replicate.
* Replicates are not automatically distinct either: 39% of nominal Protenix replicates
  were duplicates, and rosettafold-3 in single-sequence mode is deterministic for about
  half of ligands. `reference_poses` dedupes.
* More engines is NOT better. The best measured reference set is the two Protenix
  checkpoints; adding esmfold2 at matched depth drops +0.0380 to +0.0178 (FINDING 011).

    python scripts/structure/build_xeng_feature.py --tag val87b \\
        --pool D:/cyp_scratch/val87b_unsteered \\
        --refs protenix_v2 protenix
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

OP_ROOT = DATA_PROCESSED / "openprotein" / "op1"


def main() -> int:
    from cypstruct import pose as P
    from cypstruct import xengine as X

    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--pool", required=True,
                    help="directory of <LIG>__*/ dirs holding the poses to be SELECTED")
    ap.add_argument("--refs", nargs="+", default=["protenix_v2", "protenix"],
                    help="engine pool names under data/processed/openprotein/op1/")
    ap.add_argument("--pattern", default="*.cif")
    ap.add_argument("--min-depth", type=int, default=4)
    ap.add_argument("--force", action="store_true",
                    help="write even if the reference set is too shallow (do not)")
    a = ap.parse_args()

    ref = X.reference_poses([OP_ROOT / e for e in a.refs], P.load_structure)
    depth = X.reference_depth(ref)
    print(f"reference engines: {a.refs}")
    print(f"reference depth: {  {k: v for k, v in depth.items() if k != 'note'} }")
    thin = [k for k, v in ref.items() if len(v) < a.min_depth]
    if thin and not a.force:
        print(f"\nREFUSING: {len(thin)} of {len(ref)} ligands have fewer than "
              f"{a.min_depth} independent reference poses, e.g. {thin[:5]}.")
        print("At one reference pose this feature measured -0.0055 - it would make the")
        print("selection WORSE, not weaker. Add replicate jobs (not samples) and re-run.")
        return 2

    rows = []
    for dirp in sorted(Path(a.pool).glob("*__*")):
        lig = dirp.name.split("__")[0]
        if lig not in ref or not ref[lig]:
            continue
        for f in sorted(dirp.glob(a.pattern)):
            try:
                v = X.in_heme_frame(P.load_structure(f))
            except Exception:
                continue
            if v is None:
                continue
            rows.append({"ligand": lig, "sample": f.stem,
                         "xeng": X.xeng_score(v, ref[lig])})
    if not rows:
        print("no poses read from the pool - check --pool and --pattern")
        return 1

    df = pd.DataFrame(rows)
    out = DATA_PROCESSED / f"xeng_{a.tag}.csv"
    df.to_csv(out, index=False)
    print(f"\n{len(df)} poses over {df.ligand.nunique()} ligands -> {out.name}")
    print("build_submission.py will now select on this automatically.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
