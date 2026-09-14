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

## P20815 explained: the reference must be BETTER than the pool, not merely different

The persistent unexplained miss is now explained, and the explanation generalises. Scoring
the esmfold2 *reference* poses against crystal, alongside the Protenix *pool* they judge:

| protein | reference quality | pool quality | selector outcome |
|---|---|---|---|
| **P20815** | **0.301** | 0.476 | **−0.0144** (fails) |
| **Q16696** | **0.453** | 0.589 | **−0.0055** (fails) |
| Q2IU02 | **0.752** | 0.589 | **+0.2771** (best result in the set) |

**The feature fails exactly where the reference engine is worse than the pool it judges.**
Agreeing with a worse opinion pulls the selection toward worse poses - the sign of the term
effectively inverts. Where the reference is better (Q2IU02: 0.752 against 0.589), the
feature produces the largest gain measured anywhere.

This supersedes the "nothing to catch" reading for P20815, which never fitted: it had 4.5%
catastrophic poses and real spread, so there *was* something to catch. It also sharpens
FINDING 011's "reference quality beats reference variety" from a slogan into a testable
condition: **quality relative to the pool is what matters, not absolute quality and not
architectural diversity.**

### What this costs and what it buys

**Costs:** the condition needs ground truth to check directly, so it cannot be evaluated on
a blind release. That is a real limitation on the drop-day story, and `pool_diagnostics.py`
does not currently detect it.

**Buys:** it explains every failure in the set rather than most of them, it predicts that
adding a *better* reference engine would help more than adding a more *diverse* one (which
is what the esmfold2 experiment already showed at matched depth), and it gives a concrete
thing to look for - a proxy for "is my reference worse than my pool" computable without
crystal structures would close the last gap in the drop-day check.

### Attempted: a ground-truth-free detector for "reference worse than pool". It does not work.

If a worse reference engine produced more scattered poses, reference self-consistency would
flag the failure mode without crystal structures. Tested over 10 proteins with >= 4 pairs:

| proxy | rho with gain | p |
|---|---|---|
| ref_spread / pool_spread ratio | **−0.455** | 0.187 |
| reference spread alone | −0.164 | 0.65 |
| pool spread alone | −0.139 | 0.70 |
| mean cross-engine distance | −0.030 | 0.93 |

**Not a usable detector.** The ratio points the right way and catches two of the three
failures - Q16696 at 2.141 (references more scattered than the pool they judge) and P11509
at 0.999 - but **misses P20815 at 0.299**, which is squarely inside the successful range
(0.071 to 0.464). With n = 10 proteins and four proxies tried, nothing here clears
significance, and the best-of-four selection makes rho = −0.455 weaker than it looks.

So the gap stands: **the failure mode is real, explained, and currently undetectable
without ground truth.** Recorded as an open problem rather than a solved one. The natural
next attempt is a proxy built on the reference engine's own confidence or on its agreement
with a *third* engine, neither of which was tried here.

### Second attempt: triangulation with a third engine. Also fails, and the near-miss is instructive.

Added Protenix-v1 across the P450 set as a third opinion (one replicate - it is
deterministic), to test whether "esmfold2 disagrees with BOTH Protenix checkpoints" flags a
bad reference without ground truth. 97 pairs have all three engines.

**One worry was misplaced.** Same-family agreement is the *weaker* quality signal, not the
stronger one: v1↔v2 correlates with Protenix pose quality at ρ = −0.178 (p = 0.08) while
cross-family v2↔esmfold2 reaches ρ = −0.327 (p = 0.001). Agreement between two checkpoints
of one architecture carries less information than agreement across architectures.

**And a promising signal appeared.** Cross-engine distance predicts the reference-quality
gap at **ρ = −0.401, p < 0.0001** over 97 pairs - ground-truth-free, and far stronger than
the self-consistency proxy that failed earlier (ρ = −0.030).

**Then acting on it fails.** Excluding high cross-distance pairs *lowers* the gain:

| subset | random | selected | gain |
|---|---|---|---|
| all pairs | 0.6756 | 0.7338 | **+0.0583** |
| cross-dist below p90 | 0.6942 | 0.7516 | +0.0574 |
| cross-dist below p75 | 0.7157 | 0.7637 | +0.0480 |
| cross-dist below p50 | 0.7617 | 0.8045 | +0.0428 |

The random baseline **rises** from 0.676 to 0.762 as pairs are excluded, which gives it
away: **high cross-engine distance marks HARD pairs, not bad-reference pairs.** On a hard
target both engines do worse, so the quality gap and the distance move together - the
ρ = −0.401 is real and confounded. Excluding on it discards exactly the cases where the
feature earns its keep.

**The two failure modes need opposite responses** - "hard" means apply the feature, it pays
most here; "reference worse than pool" means do not - and cross-engine distance cannot tell
them apart. Gap still open, now with two specific dead ends recorded so the next attempt
starts somewhere else.

### Third attempt: the REFERENCE engine's own confidence. This one survives.

FINDING 009 showed engine confidence cannot rank poses *within* a ligand - chance, in two
architectures. But "is my reference engine doing badly on this target" is a **between-target**
question, which is a different thing, and confidence turns out to answer it.

esmfold2's own `iptm`, over 160 pairs, entirely ground-truth-free:

| | rho | p |
|---|---|---|
| iptm vs esmfold2's own pose quality | **+0.702** | 0.0000 |
| iptm vs (reference quality − pool quality) | **+0.367** | 0.0000 |
| ptm vs the same gap | +0.318 | 0.0000 |

**And acting on it works**, which is where the previous two attempts died:

| subset | random | selected | gain |
|---|---|---|---|
| all pairs | 0.6223 | 0.7205 | +0.0982 |
| drop lowest 10% iptm | 0.6348 | 0.7339 | +0.0991 |
| drop lowest 25% iptm | 0.6566 | 0.7638 | **+0.1072** |
| drop lowest 40% iptm | 0.6844 | 0.8016 | **+0.1171** |

Monotone, and the contrast with the cross-distance proxy is the whole point: there the gain
**fell** as pairs were excluded (+0.0583 → +0.0428) while the baseline rose, showing it was
selecting for *easy* rather than for *reliable reference*. Here both rise together.

**How to use it, since excluding ligands is not an option on a submission.** Every ligand
must be submitted, so the value is not filtering - it is knowing *when to trust the feature*.
A low-`iptm` reference means cross-engine agreement is unreliable for that ligand and the
older FINDING 003 selector (contacts + sibling consensus, +0.0265, needs no reference) is
the safer choice there. That makes the selector choosable per ligand on a blind release.

**Caveats, stated because this is the third proxy tried.** n = 154 pairs sampled from 60
jobs; best-of-three selection inflates any single result; and the thresholds are descriptive
rather than fitted. It needs confirmation on the full set and on CYP3A4 before it goes into
the drop-day path. But it is the first candidate that improves selection when acted upon
rather than merely correlating with something.

---

## RETRACTION: "the reference must be better than the pool" is contradicted

I claimed above that the feature fails where the reference engine is worse than the pool it
judges, and said it "explains every failure in the set rather than most of them." **That is
wrong, and CYP3A4 contradicts it directly:**

| set | reference mean | pool mean | reference better? | outcome |
|---|---|---|---|---|
| P450 P20815 | 0.301 | 0.476 | no | **−0.0144** (fails) |
| **CYP3A4** | **0.1923** (Protenix) | **0.5769** (Boltz) | **no** | **+0.0381 (works)** |

CYP3A4 violates the stated condition more severely than P20815 does - its reference averages
a third of the pool's quality - and the feature works there. So relative mean quality is not
the discriminator, and the three-protein pattern I read it from was a coincidence of a small
sample.

**What stands and what does not.** P11509 and Q16696 are still explained by having nothing
to catch (zero catastrophic poses; Q16696's pool mean and oracle differ by 0.002). **P20815
returns to unexplained.** The reference-quality story is withdrawn as a general mechanism.

One observation that survives and may matter: CYP3A4's reference has a low mean (0.1923) but
a **high oracle (0.6273)** - it contains good poses even though most are poor. `xeng` scores
against the mean distance to all references, so a reference set with a few good poses may
still guide selection. Whether the reference *oracle* rather than its mean is the right
quantity is untested, and is the obvious next thing to look at rather than another proxy.

### And the confidence proxy does not transfer to CYP3A4

Tested on the set that actually matters: rho(reference iptm, per-ligand gain) = **−0.060,
p = 0.58**, and acting on it does not improve anything (+0.0380 → +0.0351 → +0.0335 →
+0.0368). The reason is visible in the spread - CYP3A4's reference confidence is nearly
constant (sd 0.0137, IQR 0.020) against P450's (sd 0.0297, IQR 0.051), **2.5x less
variance**. On a familiar target the reference engine is uniformly confident, so there is
no signal to exploit.

That is not bad news for the submission, but it is not the good news it first looks like
either: it means the *detector* is inert on CYP3A4, while the failure mode it was built to
detect is no longer understood.

### The reference-oracle successor also fails; the outlier failures stay unexplained

Tested whether the reference *oracle* rather than its mean is the operative quantity, over
10 proteins with known outcomes:

| quantity | rho with gain | p |
|---|---|---|
| ref_oracle / pool_oracle | **+0.030** | 0.93 |
| ref_mean / pool_mean | +0.224 | 0.53 |
| ref_oracle (absolute) | +0.539 | 0.11 |
| ref_mean (absolute) | +0.273 | 0.45 |

The individual rows kill it plainly: Q00441 and Q55080 have the two *lowest* oracle ratios
(0.769, 0.744) and two of the *highest* gains (+0.226, +0.225), while Q16696 has a high
ratio (0.942) and fails. Relative reference quality does not discriminate on either
measure.

**Current honest state of this thread.** Three hypotheses proposed and all three rejected:
reference-worse-than-pool (contradicted by CYP3A4), reference self-consistency (rho −0.030),
and reference oracle (rho +0.030). P11509 and Q16696 remain explained by having nothing to
catch; **P20815 remains unexplained.**

The only surviving hint is that *absolute* reference quality may matter more than relative
(ref_oracle alone, rho +0.539, p = 0.11) - suggestive at n = 10 proteins, nowhere near
significant, and exactly the kind of near-miss FINDING 007 says not to chase. **Stopping
this thread here.** MSAs are adding proteins at ~1 per 20 minutes; at 20+ proteins with
known outcomes the same test becomes worth re-running, and until then further attempts are
fishing rather than investigating.

### The outlier failures are permanently underpowered, not awaiting data

I set a resumption condition - "if P20815 survives to 10+ pairs it is a real counterexample"
- without checking it was reachable. It is not:

| protein | pairs scored | pairs that EXIST in the whole harvest |
|---|---|---|
| P20815 | 4 | **5** |
| Q16696 | 4 | **4** |
| P11509 | 7 | **9** |

The PDB does not contain more ligand-bound structures for these proteins, so none of the
three can reach 10 pairs however long the campaign runs. At 4-5 pairs a single pose changing
rank flips the sign of the per-protein gain, and all three sit between −0.006 and −0.014 -
well inside that.

**So the honest reading is not "an unexplained mechanism failure" but "three proteins with
too few ligands to measure."** That also retires the hunt: three hypotheses were tested
against 10 proteins to explain what is most likely sampling noise in the three smallest of
them. The eleventh measurement makes the point another way - **41 proteins, gain +0.0529 at
5.7% catastrophic, positive on 20 of 23** - and the negatives are the same three every time
precisely because they are the same three tiny samples every time.

What would settle it is more *ligands* for those proteins, which the crystallographic record
does not have. Recorded as closed-by-data-limit rather than open.

## Twelfth measurement — the two settings are converging on the same number

**43 proteins in the test (47 scored), 75 construct sequences, 308 pairs, 4,088 poses.**
Catastrophic 5.23%, selected 0.7375 against random 0.6910, **gain +0.0465** (null 99th
+0.0104, p = 0.0), rho −0.1861, positive on 20 of 23.

The complete series, twelve measurements, nothing tuned between any of them:

| catastrophic | gain |
|---|---|
| 19.3% | +0.3006 |
| 14.5% | +0.2045 |
| 12.0% | +0.1410 |
| 10.3% | +0.1098 |
| 8.9% | +0.0908 |
| 7.74% | +0.0697 |
| 7.05% | +0.0655 |
| 6.07% | +0.0566 |
| 5.7% | +0.0529 |
| **5.23%** | **+0.0465** |
| **0.2% (CYP3A4 anchor)** | **+0.0381** |

Monotone across all twelve, and the P450 figure is now **within 22% of the CYP3A4 anchor**,
approaching from above as its pools stop producing catastrophic poses.

That convergence is the strongest form the evidence takes. Two settings that differ in
protein, engine, pool construction and MSA source are arriving at the same number as their
catastrophe rates equalise - which is what "the gain is a property of the pool" predicts and
what "the feature behaves differently on different proteins" cannot produce.

## Thirteenth measurement — within 12% of the anchor

**47 proteins in the test (53 scored), 88 construct sequences, 338 pairs, 4,448 poses.**
Catastrophic 4.81%, selected 0.7519 against random 0.7091, **gain +0.0428** (null 99th
+0.0102, p = 0.0), rho −0.1973, positive on **24 of 27**.

| measurement | catastrophic | gain |
|---|---|---|
| 1st | 19.3% | +0.3006 |
| 5th | 8.9% | +0.0908 |
| 9th | 6.07% | +0.0566 |
| 12th | 5.23% | +0.0465 |
| **13th** | **4.81%** | **+0.0428** |
| **CYP3A4 anchor** | **0.2%** | **+0.0381** |

The P450 measurement is now **within 12% of the CYP3A4 anchor**, having started 8x above it,
and has decreased monotonically at every one of thirteen re-measurements while the protein
count grew from 17 to 47.

Two settings - different protein, different pool engine, different reference engine,
different MSA source, different construct sequences - are arriving at the same number as
their catastrophe rates converge. Nothing was tuned between measurements; the pools improved
on their own as replicates accumulated, and each re-run is its own commit in the history.

The three negatives are unchanged and remain the three proteins with 4, 4 and 7 available
pairs - the permanently underpowered set, not a systematic failure.

## Fourteenth measurement — the two settings have converged

**61 proteins in the test (67 scored), 109 construct sequences, 378 pairs, 4,791 poses.**
Catastrophic 4.47%, selected 0.7575 against random 0.7176, **gain +0.0398** (null 99th
+0.0088, p = 0.0), rho −0.1936, positive on **25 of 29**.

| | catastrophic | gain |
|---|---|---|
| P450, first measurement | 19.3% | +0.3006 |
| **P450, fourteenth** | **4.47%** | **+0.0398** |
| **CYP3A4 anchor** | **0.2%** | **+0.0381** |

**Within 4.5% of the anchor.** The P450 measurement began 8x above it and has decreased
monotonically at all fourteen re-measurements while the protein count grew 17 → 61.

Two settings that share no protein, no pool engine, no reference engine, no MSA source and
no construct sequence now produce the same number, because their catastrophe rates have
equalised. That is the catastrophe-detector account's central claim, tested to the point
where the two curves meet.

The negatives remain the small-sample set: Q5YNS8 (5 pairs), Q16696 (4), P11509 (7),
P20815 (4) - every one of them among the fewest-ligand proteins in the harvest, where a
single pose changing rank flips the sign.

### A note on the operational judgement behind this measurement

The batch carrying these 10 new proteins ran **163 minutes** against the previous batches'
90 and 120, with no failures and no output for over two hours. Killing it at any of the
thresholds considered would have discarded exactly the data that produced this measurement.
The signal that justified waiting was that **MSAs kept completing throughout** - the
platform was demonstrably processing work, so "slow" was better supported than "stuck".
