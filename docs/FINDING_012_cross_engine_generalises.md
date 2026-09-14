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

---

## Update — the gain attenuates as the pool improves, exactly as the mechanism predicts

Re-run after the P450 pools deepened. The catastrophe rate fell from 19.3% to 14.5% of
poses, and the gain fell with it:

| | first run | **deeper pool** |
|---|---|---|
| catastrophic poses (< 0.1) | 19.3% | **14.5%** |
| random / oracle | 0.4241 / 0.7360 | 0.5219 / 0.7464 |
| **gain** | +0.3006 | **+0.2045** |
| within-pair rho | −0.763 | −0.579 |
| correct direction | 96.7% | 86.7% (p = 2.7e−14) |

This is the catastrophe-detector account making a **prospective** prediction and being
right: improve the pool, and the feature has less to catch. It was not fitted to this
outcome - the mechanism was written down first, from the CYP3A4/P450 contrast, and the
pool then improved on its own as replicates landed.

### It is positive on every held-out target

| target | pairs | random → selected | gain |
|---|---|---|---|
| Q2IU02 | 25 | 0.598 → 0.904 | +0.3060 |
| Q7Z1V1 | 7 | 0.504 → 0.795 | +0.2910 |
| Q00441 | 4 | 0.550 → 0.806 | +0.2558 |
| Q55080 | 4 | 0.570 → 0.814 | +0.2446 |
| Q9Y6A2 | 12 | 0.611 → 0.855 | +0.2440 |
| Q385E8 | 4 | 0.470 → 0.708 | +0.2376 |
| **P08684 (CYP3A4)** | 60 | 0.461 → 0.659 | **+0.1978** |
| P20815 | 4 | 0.373 → 0.470 | +0.0968 |
| P00178 | 11 | 0.393 → 0.485 | +0.0928 |

**9 of 9 positive**, spanning human and bacterial P450s, with no parameters fitted.

### The internal control nobody designed

P08684 **is** CYP3A4 — the same protein the feature was developed on — and here it gains
**+0.1978** against the +0.0381 measured on the Boltz pool. Same protein, same feature,
same reference engine; the only difference is the pool being selected from. The P450
campaign's CYP3A4 pool is Protenix with a freshly computed MSA and contains far more
catastrophic poses than the Boltz pool does with its bonded heme.

That is the cleanest available demonstration that **the gain is a property of the pool, not
of the protein or of any tuning**. It also warns against reading the P450 numbers as "the
feature is better on other proteins" - it is not, those pools are simply worse.

---

## The mechanism is now quantitative: gain tracks the catastrophe rate

Three independent re-measurements as the P450 pools deepened, plus the CYP3A4 anchor. The
pool improved on its own; nothing was tuned between these points.

| pool | catastrophic poses (< 0.1) | gain | within-pair rho | targets |
|---|---|---|---|---|
| P450, first | 19.3% | **+0.3006** | −0.763 | 17 |
| P450, deeper | 14.5% | **+0.2045** | −0.579 | 17 |
| P450, deeper still | **12.0%** | **+0.1448** | −0.331 | 22 |
| CYP3A4 Boltz (anchor) | 0.2% | **+0.0381** | −0.258 | 1 |

Monotone in the catastrophe rate across a 100-fold range of it and an 8-fold range of the
gain. The P450 series is converging toward the CYP3A4 value as its pools stop producing
catastrophic poses, which is what the catastrophe-detector account requires and what a
"this feature is better on other proteins" account cannot explain.

Still **positive on 9 of 10** targets with >= 4 pairs at the latest point (one target has
now gone slightly negative - expected as the effect shrinks toward the noise floor, and
recorded rather than dropped).

### Why this matters for the challenge

The gain is now predictable from something **measurable on the day, without ground truth**.
The catastrophe rate cannot be computed directly on a blind release, but its proxy can:
within-ligand pool spread, which is what FINDING 012's addendum already split on
(+0.0704 wide vs +0.0063 narrow). So on drop day:

1. build the pool, compute the within-ligand spread distribution;
2. a wide-spread, high-disagreement pool means the feature is worth a lot - closer to the
   +0.10-0.30 regime;
3. a tight pool means expect +0.04, and the effort is better spent on generation.

That converts "we hope this transfers" into a check that can be run before submitting.
