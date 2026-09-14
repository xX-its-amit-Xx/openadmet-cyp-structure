"""Before submitting: is this pool one where cross-engine selection is worth much?

FINDING 012 makes the payoff predictable from the pool alone. The gain tracks how often
the pool contains catastrophic poses, measured across four independent pools (the
catastrophe->gain relation rests on four points; the SPREAD PROXY below rests on two):

    catastrophic 19.3% -> +0.3006      catastrophic 12.0% -> +0.1448
    catastrophic 14.5% -> +0.2045      catastrophic  0.2% -> +0.0381

The catastrophe rate itself needs ground truth, so it cannot be computed on a blind
release. Its **proxy can**: within-ligand pose spread. Ligands whose poses disagree are
the ones where a catastrophe is hiding, and FINDING 012's addendum measured the split
directly - wide-spread CYP3A4 ligands gained +0.0704, narrow-spread ones +0.0063 (null).

So this reports what the pool looks like and which regime to expect, **without any
reference structures**, and it reports its own calibration alongside the answer,
because that calibration rests on exactly TWO measured pools.

**The first version of this tool was wrong**, and instructively so: its thresholds were
guessed before either anchor was measured, and it called the CYP3A4 pool "moderate,
+0.05 to +0.15" when the true gain is +0.0381. The proxy itself was fine - CYP3A4 sits at
0.51 and the P450 pool at 0.83, in the right order - the numbers bolted onto it were not.

    python scripts/ops/pool_diagnostics.py --xeng data/processed/xeng_val87b.csv
    python scripts/ops/pool_diagnostics.py --pool-dir <dir> --refs protenix_v2 protenix
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

# CALIBRATION: exactly two anchor pools, and the tool says so rather than implying more.
# The first threshold table here was wrong - it called the CYP3A4 pool "moderate,
# +0.05 to +0.15" when its true gain is +0.0381, because the thresholds were guessed
# before either anchor was measured. Measured anchors:
#
#   median relative spread 0.51  ->  gain +0.0381   (CYP3A4 Boltz pool, 0.1% catastrophic)
#   median relative spread 0.83  ->  gain +0.1448   (P450 protenix pool, 12% catastrophic)
#
# Two points define a line and nothing more, so this reports the measurement next to the
# anchors and interpolates coarsely. It does NOT extrapolate far beyond them.
ANCHORS = [(0.51, 0.0381, "CYP3A4 Boltz, 0.1% catastrophic"),
           (0.83, 0.1448, "P450 protenix, 12% catastrophic")]


def regime(med: float) -> tuple[str, str, str]:
    lo_s, lo_g, _ = ANCHORS[0]
    hi_s, hi_g, _ = ANCHORS[1]
    if med < lo_s - 0.10:
        return ("tighter than any pool measured", "~+0.02 or less",
                "the pool already agrees; generation effort likely pays more")
    if med > hi_s + 0.15:
        return ("wider than any pool measured", ">= +0.15 (extrapolated - unreliable)",
                "the feature should be worth a lot, but this is beyond calibration")
    t = (med - lo_s) / (hi_s - lo_s)
    est = lo_g + t * (hi_g - lo_g)
    band = "tight" if t < 0.33 else ("moderate" if t < 0.67 else "wide")
    return (band, f"about +{est:.3f} (interpolated between two anchors)",
            "worth applying; confirm the reference set has >= 4 independent poses")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xeng", default=None,
                    help="an xeng_<tag>.csv; spread is measured on the xeng values")
    ap.add_argument("--scored", default=None,
                    help="optional poses_scored csv - if given, reports the TRUE "
                         "catastrophe rate as well (needs ground truth)")
    a = ap.parse_args()

    if not a.xeng:
        print("give --xeng (build it with scripts/structure/build_xeng_feature.py)")
        return 1
    df = pd.read_csv(a.xeng)
    if "ligand" not in df or "xeng" not in df:
        print("expected columns ligand,sample,xeng")
        return 1

    # Disagreement proxy: spread of the cross-engine distance WITHIN each ligand,
    # normalised by its own median so ligand size does not dominate.
    g = df.groupby("ligand").xeng
    rel = ((g.transform("max") - g.transform("min")) /
           g.transform("median").replace(0, np.nan))
    per_lig = rel.groupby(df.ligand).first().dropna()
    med = float(per_lig.median())

    print(f"pool: {len(df)} poses over {df.ligand.nunique()} ligands")
    print(f"within-ligand cross-engine spread (relative): "
          f"p25 {per_lig.quantile(0.25):.2f}  median {med:.2f}  "
          f"p75 {per_lig.quantile(0.75):.2f}")
    # NB: the fraction of ligands above a fixed cut is NOT usable - it runs the wrong way
    # (CYP3A4 77% vs the P450 pool's 65%, while the P450 pool has the larger gain). Only
    # the median relative spread tracks the regime.
    label, rng, advice = regime(med)
    print("\nanchors (measured):")
    for sp, gn, what in ANCHORS:
        print(f"   spread {sp:.2f} -> gain +{gn:.4f}   {what}")
    print(f"\nREGIME: {label}")
    print(f"  expected cross-engine gain: {rng}")
    print(f"  {advice}")
    print("\n  TWO anchor pools only. The DIRECTION is solid - more within-ligand")
    print("  disagreement means more for the feature to resolve - but any number here")
    print("  is an interpolation between two points. Treat it as an order of magnitude.")

    if a.scored:
        sc = pd.read_csv(a.scored)
        col = "lddt_pli"
        if col in sc:
            cat = float((sc[col] < 0.1).mean())
            print(f"\n  TRUE catastrophe rate (needs ground truth): {cat:.1%}")
            print("  On a blind release this line is unavailable - that is the whole")
            print("  reason the spread proxy exists.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
