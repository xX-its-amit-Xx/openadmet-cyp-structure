# Finding 004 — I was wrong: more sampling still pays

**Date:** 2026-09-11 · **Method:** subsample the existing 20-sample pools, 30 draws per
ligand per size · **n = 87 ligands, 1,740 poses**

FINDING 002 ended with "**Do not spend more GPU on sampling.** The oracle is 0.6975 and
the PXR winning entry was 0.564. There is nothing wrong with the poses we have."

That was an assertion, not a measurement, and it was wrong. It confused the pool's
*ceiling* with what a selector can actually reach.

## The curve

| samples | oracle | selector (FINDING 003) | random |
|---|---|---|---|
| 2 | 0.6114 | 0.5798 | 0.5768 |
| 3 | 0.6285 | 0.5836 | 0.5769 |
| 5 | 0.6464 | 0.5863 | 0.5769 |
| 8 | 0.6637 | 0.5886 | 0.5765 |
| 10 | 0.6739 | 0.5923 | 0.5779 |
| 15 | 0.6860 | 0.5984 | 0.5770 |
| **20** | **0.6975** | **0.6048** | 0.5769 |

**The oracle is still climbing at 20 samples** — about **+0.024 per doubling** — with no
sign of saturation.

**And the selector tracks it: +0.0125 per doubling.** Random is flat at 0.577 throughout,
as it must be, which is a useful check that the subsampling is unbiased.

## Why the original claim was wrong

"The oracle already beats the PXR winner" is true and irrelevant. The oracle is what a
*perfect* selector would score. What we submit is what our *actual* selector picks, and
that number rises with pool size because a bigger pool contains better poses for the
selector to find — even a weak one.

## What this is worth, against the alternative

Two levers, priced the same way:

| lever | gain | cost |
|---|---|---|
| doubling samples 20 → 40 | **+0.0125** | ~9.6 GPU-h |
| the FINDING 003 selector, over the incumbent | +0.028 | ~0 (analysis only) |

Selector work remains the better buy per unit of effort, and it compounds — every future
sampling gain is converted at the selector's rate. But sampling is no longer "don't", it
is "cheap, additive, and worth doing in parallel with the analysis".

## Action

Run a second seed on the 87-ligand set at 20 samples (≈9.6 GPU-h), taking each pool to 40
samples. Then re-measure both curves. If the oracle is still climbing at 40, run a third.

Budget check: 19.6 of 60 GPU-h used, plus 7.4 in flight for Chai.
