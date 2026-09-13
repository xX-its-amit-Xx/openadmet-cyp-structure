# FINDING 012 — cross-engine agreement generalises to 17 unseen P450s, and it is a catastrophe detector rather than a fine-ranker

**Date:** 2026-09-13
**Pool:** 732 Protenix-v2 poses over **141 pairs across 17 non-CYP3A4 targets**
**Reference:** esmfold2 poses of the same pair, one per replicate, deduplicated
**Script:** `scripts/structure/score_p450_pool.py`, `cypstruct.xengine`

---

## The result the harvest was built for

Every selector number in this repo, including FINDING 011's +0.0381, rests on **one
protein**. The whole point of harvesting the P450 superfamily was to find out whether that
survives contact with a pocket the method has never seen. It does:

| | CYP3A4 (FINDING 011) | **17 held-out P450 targets** |
|---|---|---|
| random | 0.5784 | 0.4241 |
| oracle | 0.6975 | 0.7360 |
| **selected** | 0.6164 | **0.7247** |
| **gain** | **+0.0381** | **+0.3006** |
| within-pair rho | −0.258 | **−0.763** |
| correct direction | 75.9% | **96.7%** |
| p | 0.0000 | 0.0000 (null 99th +0.0600) |

No parameters are fitted, so every target is held out by construction. The feature was
developed entirely on CYP3A4 and applied unchanged.

## Why the effect is 8x larger, and why that is NOT a better result

An eightfold jump is the kind of number that should be distrusted before it is reported.
It is not leakage — both sides are predictions, neither sees ground truth. It is pool
composition:

| | P450 pool | CYP3A4 pool |
|---|---|---|
| poses below 0.1 LDDT-PLI (catastrophic) | **19.2%** | 0.2% |
| poses above 0.6 | 36.1% | 51.9% |
| within-pair range (median max−min) | **0.631** | 0.234 |

The P450 pool is **bimodal**: on an unfamiliar target Protenix either places the ligand
essentially right or misses entirely. Separating those two modes is easy, and cross-engine
agreement does it almost perfectly (96.7%). The CYP3A4 Boltz pool contains virtually no
catastrophes — with the heme bonded to Cys442 the model is reliable — so there the same
feature has to *fine-rank* among mostly-decent poses, which is a much harder problem and
yields ρ = −0.26 rather than −0.76.

**So the mechanism is: cross-engine agreement detects catastrophic failure, and its value
scales with how often the pool contains any.** That reframes it from "a better selector"
to "a safety net whose payoff depends on how bad the pool's worst poses are".

## What this predicts for the challenge

The released CYP3A4 structures are the trained-on, well-behaved case, so the honest
expectation is the **+0.038 regime, not +0.30**. The +0.30 is what the feature is worth
when co-folding actually fails, which on CYP3A4 with a bonded heme it rarely does.

Two things would move the release toward the P450 regime, and both are plausible given
OpenADMET's own analysis: the F/G loop remodelling they name as why co-folding does poorly
here, and ligands with multiple mutually exclusive conformations. If the release behaves
like an unfamiliar target rather than like our validation set, this feature is worth far
more than +0.038.

## Caveats

- P450 pools are shallow (3-7 poses per pair) against CYP3A4's 20, so the P450 oracle is
  reached more easily and `selected` sits closer to it.
- esmfold2 reference depth is mostly 2, below the 4 the CYP3A4 dose-response wants. The
  effect is large enough that this is unlikely to reverse it, but the number will move.
- 17 targets is not 185. The remaining MSAs are still computing.
