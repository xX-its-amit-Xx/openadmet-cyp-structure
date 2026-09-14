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


---

## Fifth measurement, more proteins — and a terminology correction

**Correction first, because it changes how every count above reads.** This repo's scoring
output says "targets" for **distinct construct sequences**, while the per-target tables
group by **UniProt accession**, i.e. distinct proteins. They are not the same: the current
pool is **44 construct sequences over 27 proteins**. The honest denominator for "does this
generalise across proteins" is the protein count, and earlier sections that say "17
targets" or "22 targets" mean *proteins* where they appear in per-target tables and
*sequences* where they come from the scorer. Numbers below state which.

| pool | catastrophic | gain | rho | proteins | pairs |
|---|---|---|---|---|---|
| P450, first | 19.3% | +0.3006 | −0.763 | 17 | 141 |
| P450, deeper | 14.5% | +0.2045 | −0.579 | 17 | 141 |
| P450, deeper still | 12.0% | +0.1410 | −0.352 | 22 | 173 |
| **P450, current** | **10.3%** | **+0.1098** | **−0.296** | **27** | **227** |
| CYP3A4 Boltz (anchor) | 0.2% | +0.0381 | −0.258 | 1 | 87 |

Five points, monotone in the catastrophe rate, converging on the CYP3A4 anchor as the pools
improve. **Positive on 13 of 16 proteins** with >= 4 pairs.

The per-protein hit rate falls as coverage grows (9/10 → 13/16) and that is expected rather
than worrying: each new protein arrives with a better pool than the early ones had, so
there is less for a catastrophe detector to catch. The clearest case is **P11509**, whose
pool is already excellent - random 0.909, selected 0.909, **delta −0.0006**. Nothing to
catch, and nothing lost by looking.

---

## The "shallow pools" caveat was wrong, and testing it retires it

Every version of this finding has carried a caveat: *the P450 pools are shallow (3-7 poses
per pair) against CYP3A4's 20, so the comparison may be unfair.* The pools have now grown
enough to test that directly rather than keep repeating it.

| subset | pairs | proteins | catastrophic | oracle | selected | gain |
|---|---|---|---|---|---|---|
| all (>= 3 poses) | 234 | 28 | 8.9% | 0.761 | 0.738 | **+0.0908** |
| deeper (>= 6 poses) | 234 | 28 | 8.9% | 0.761 | 0.738 | **+0.0908** |
| **deepest (>= 8 poses)** | 147 | 16 | **11.7%** | 0.739 | 0.712 | **+0.1402** |

**Depth is not the driver.** The deepest subset gives the *larger* gain, which rules out
"shallow pools inflate the effect" outright. It gives the larger gain because it happens to
contain harder pairs - 11.7% catastrophic against 8.9% - which is the mechanism again, now
confirmed on a split chosen for depth rather than for difficulty.

Two more points for the series, both consistent:

| catastrophic | gain |
|---|---|
| 11.7% | +0.1402 |
| 8.9% | +0.0908 |

So the caveat is retired. What remains true, and is a different claim, is that these pools
are **built differently** from the CYP3A4 Boltz pool - Protenix with a fresh MSA rather than
Boltz with a bonded heme - and that difference is exactly why their catastrophe rates
differ. The gain follows the catastrophe rate across seven measurements now, whether the
pool is shallow or deep, CYP3A4 or bacterial.

## The reference-depth caveat: partly real, and it needs reading carefully

The other standing caveat was that esmfold2 supplies only ~2 reference poses per pair,
below the >= 4 the CYP3A4 dose-response wants. Tested:

| esmfold2 reference depth | pairs | catastrophic | gain | p |
|---|---|---|---|---|
| >= 2 | 234 | 8.9% | +0.0908 | 0.0000 |
| **>= 3** | 170 | 9.7% | **+0.1028** | 0.0000 |
| >= 4 | **21** | **3.3%** | +0.0162 | 0.29 |

**The >= 4 row is not a depth failure and must not be read as one.** It is 21 pairs at a
3.3% catastrophe rate - an unusually *easy* subset, where the mechanism predicts almost no
gain and delivers almost none. Reading it as "deep references stop working" would invert
the actual cause, which is that those pairs have nothing to catch.

What is genuinely true: **two to three independent references are enough here**, against
the >= 4 the CYP3A4 dose-response required. The likely reason is that the two settings ask
different questions - catching a catastrophe is a coarse judgement that a couple of
independent opinions settle, while fine-ranking mostly-decent poses needs a better estimate
of where the other engine thinks the ligand goes. The threshold is a property of the task,
not a constant of the feature.

**A measurement worth carrying separately:** esmfold2's raw replicates dedupe heavily -
median 5 raw collapses to 2-4 distinct - so its replicate count overstates independence
roughly as much as Protenix's 39% did. Deduplication is doing real work on every engine
tested so far.


---

## Reproducible in one command, and the eighth point

`python scripts/structure/score_p450_pool.py validate` now *is* this test. It replaced an
earlier `validate` that fitted three selector weights leave-one-target-out and measured
+0.0149 - inside the noise, because fitting weights on this much data overfits, exactly as
FINDING 002 found for a fitted ranker. The unweighted single term is what ships, so it is
what gets tested.

Latest run - **30 proteins, 50 construct sequences, 247 pairs, 2,325 poses**:

| | value |
|---|---|
| catastrophic poses | **7.74%** |
| random / oracle | 0.6685 / 0.7619 |
| **selected** | **0.7382** |
| **gain** | **+0.0697** (null 99th +0.0084, empirical p = 0.0) |
| within-pair rho | −0.2314 |
| **proteins positive** | **17 of 19** with >= 4 pairs |

The eighth measurement, and the series still runs monotone with the catastrophe rate as
the pools improve: 19.3% → +0.3006, 14.5% → +0.2045, 12.0% → +0.1410, 10.3% → +0.1098,
8.9% → +0.0908, **7.74% → +0.0697**, against the CYP3A4 anchor at 0.2% → +0.0381.

It is converging on the anchor from above, which is what the catastrophe-detector account
requires and what no "better on other proteins" account would produce. **17 of 19 proteins
positive** is the widest protein coverage yet, and the gain is now within a factor of two
of the CYP3A4 number - as the P450 pools become as good as the CYP3A4 one, the two
measurements are becoming the same measurement.

## Ninth measurement, and an honest look at where it fails

**34 proteins, 56 construct sequences, 268 pairs, 2,781 poses.** Catastrophic 7.05%,
random 0.6691, oracle 0.7584, selected 0.7347, **gain +0.0655** (null 99th +0.0138,
empirical p = 0.0), rho −0.2068, **positive on 17 of 20** proteins with >= 4 pairs.

Ninth point, still monotone in the catastrophe rate as the pools improve, still converging
on the CYP3A4 anchor (0.2% → +0.0381).

### The three proteins where it loses

All three losses are near zero (−0.0, −0.0056, −0.0106), but "near zero" is not an
explanation, so:

| protein | pairs | pool mean | catastrophic | oracle | explained? |
|---|---|---|---|---|---|
| P11509 | 5 | 0.895 | **0.000** | 0.914 | yes - an excellent pool with nothing to catch |
| Q16696 | 4 | 0.589 | **0.000** | 0.591 | yes - pool mean ≈ oracle, every pose is equally good |
| **P20815** | 4 | 0.461 | 4.5% | 0.591 | **no - a genuine miss** |

The negative set averages **1.5% catastrophic against 4.9% elsewhere**, so the pattern
holds in aggregate: the feature loses where there is nothing to win. Two of the three are
textbook - Q16696's pool mean and oracle differ by 0.002, meaning its poses are
interchangeable and any selector is picking between equivalents.

**P20815 is not explained by the mechanism.** It has real spread (0.461 mean against a
0.591 oracle) and 4.5% catastrophic poses, so there was something to catch, and the feature
caught none of it. At 4 pairs this is well inside sampling noise and should not be
over-read - but it is recorded as an unexplained miss rather than folded into the
"nothing to catch" story it does not fit.

## Tenth measurement — and the failures are the SAME proteins every time

**38 proteins, 63 construct sequences, 283 pairs, 3,412 poses.** Catastrophic 6.07%,
selected 0.7338 against random 0.6773, **gain +0.0566** (null 99th +0.0126, p = 0.0),
rho −0.2032, **positive on 20 of 23** proteins with >= 4 pairs.

The full series, ten measurements, nothing tuned between them:

| catastrophic | 19.3% | 14.5% | 12.0% | 10.3% | 8.9% | 7.74% | 7.05% | **6.07%** | 0.2% |
|---|---|---|---|---|---|---|---|---|---|
| **gain** | +.3006 | +.2045 | +.1410 | +.1098 | +.0908 | +.0697 | +.0655 | **+.0566** | +.0381 |

Monotone throughout, converging on the CYP3A4 anchor as the P450 pools stop producing
catastrophic poses.

### The negatives are stable, which is the useful part

Across the last three measurements the losing proteins are **the same three every time** —
Q16696, P11509, P20815 — not a rotating cast. If these were sampling noise, different
proteins would drop below zero on each re-run as the pools changed. They do not.

Two are explained: P11509 and Q16696 have **zero** catastrophic poses, and Q16696's pool
mean and oracle differ by 0.002, so its poses are interchangeable and any selector is
choosing between equivalents. **P20815 remains unexplained** and has drifted slightly
further negative (−0.0 → −0.0144) as its pool grew. At 4 pairs that is still inside
sampling noise, but it is now a *persistent* unexplained miss rather than a one-off, which
is worth more attention than its magnitude suggests. If it survives to 10+ pairs it is a
real counterexample to the mechanism and should be treated as one.
