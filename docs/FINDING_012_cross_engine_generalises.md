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

---

## Addendum — the mechanism is confirmed on CYP3A4, and it makes the feature self-limiting

The catastrophe-detector explanation predicts that the gain should track how uncertain the
pool is, not which protein it came from. Testing that on the CYP3A4 set, split at the
median within-ligand spread:

| CYP3A4 subset | ligands | random | selected | gain | rho | p |
|---|---|---|---|---|---|---|
| **wide pool spread (> 0.23)** | 43 | 0.587 | 0.657 | **+0.0704** | −0.387 | 0.0000 |
| **narrow pool spread (<= 0.23)** | 44 | 0.570 | 0.577 | **+0.0063** | −0.132 | 0.17 |

An **11x difference** on the same protein, the same engines and the same feature, driven
only by how much disagreement there is to resolve. The narrow-spread half is a clean null.

This is the mechanism holding up under a test that could have refuted it, and the P450
half could not supply that test: only 13 of 141 P450 pairs have a worst pose above 0.2, so
there is no clean subset there to compare against. The converse split does work - P450
pairs that contain a catastrophe score +0.3343 against +0.3099 overall.

### Why this is a good property rather than a caveat

The feature **pays where the pool is uncertain and costs nothing where it is not**:
+0.0704 when there is disagreement to resolve, +0.0063 (indistinguishable from zero) when
there is not. It is a safety net, not a gamble. Adding it cannot meaningfully hurt a
well-behaved target, which is exactly the property that makes it safe to apply blind to a
release whose difficulty we cannot know in advance.

It also explains the headline +0.0381 as an average over two regimes rather than a single
effect: roughly half the CYP3A4 ligands are getting +0.07 and the other half nothing.

### One number that does not fit, and is not being explained away

Splitting instead on whether the pool contains a genuinely bad pose (< 0.25) gives
**gain +0.0767 with rho +0.005** on 15 ligands - a large mean gain with *zero* rank
correlation. The two statistics disagree completely. The plausible reading is that
avoiding one catastrophic pose lifts the mean while leaving the ordering of the remaining
good poses uninformative, but n=15 is too small to call it, and it is recorded as
unresolved rather than folded into the story.
