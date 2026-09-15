# FINDING 015 — protenix_v2 is deterministic, and the pool is ~70% duplicates

**Status: stands.** Measured 2026-09-15 on 12,230 protenix_v2 poses over 491 P450 pairs.

## The measurement

The pool was doubled from 12 to 24 replicates per pair — 5,892 new poses, every
replicate complete at 491 pairs, nothing failed. Then the oracle was compared at the two
depths over identical files and references (`validate --max-rep`):

| | poses | poses/pair | pool mean | oracle |
|---|---|---|---|---|
| depth 12 | 6,314 | 12.9 | 0.6835 | **0.753818** |
| depth 24 | 12,182 | 24.9 | 0.7024 | **0.753818** |

**Pairs where the oracle improved: 0 of 489. Delta +0.000000.**

An oracle is a per-pair maximum. It can only rise when depth rises. Rising on *zero* of
489 pairs is not saturation — saturation still moves some pairs. It means the new poses
are not new.

## What is actually on disk

md5 over the pose files, per pair:

```
1EA1_TPF   reps 0-11: 4 distinct / 12    reps 12-23: 1 distinct / 12   (already present earlier)
1EGY_9AP   reps 0-11: 4 distinct / 12    reps 12-23: 1 distinct / 12
1EUP_ASD   reps 0-11: 4 distinct / 12    reps 12-23: 1 distinct / 12
1IZO_PAM   reps 0-11: 3 distinct / 12    reps 12-23: 1 distinct / 12
```

Replicates 12–23 are **one file written twelve times**, byte-identical to a file that
already existed. `lddt_pli` collapses to a single repeated value on **489 of 489 pairs**.

And the older half is barely better: 3–4 distinct files out of 12. **The pool's nominal
depth of 12 is a real depth of about 4.**

Submission parameters were identical across every replicate (`samples=1`,
`single_sequence=False`, `num_recycles=3`), so the diversity in replicates 0–3 did not
come from how we called the engine. The engine's behaviour changed under us.

## Why the pool mean rose while the oracle did not

This is the part that nearly sold the result as progress. Pool mean went 0.6835 → 0.7024
and would have been reported as a gain. It is an artifact: duplicating whichever pose the
engine deterministically returns reweights the mean toward that pose. Nothing improved.
**A mean that moves while its own oracle is frozen to six decimals is a duplication
signature, not an improvement.**

## What this corrects

- RUNBOOK said "esmfold2 and protenix_v2 are properly diverse." **Wrong for
  protenix_v2**, at least as the API behaves now. It belongs on the deterministic list
  next to protenix-v1 and rosettafold-3.
- FINDING 004's "+0.0125 per doubling" **cannot be bought from this engine**. Not
  refuted — untestable here, because depth cannot be purchased at all.
- Every "poses per pair" figure over this set overstates the real choice set ~3×.
  `validate` reporting `poses_per_pair 12.91` is counting duplicates.

## What it does NOT invalidate

The +0.0358 generalisation gain stands. Selection and its random baseline are both drawn
from the same pool, so the comparison is internally consistent — the selector really does
beat random *on the pool we have*. What changes is the interpretation: it is choosing
among ~4 genuine options per pair, not 13. The follow-up worth doing is re-measuring on
**deduplicated** pools, to state the gain per distinct opinion rather than per file.

## The lever, and what is left to test

`protenix_v2.fold` exposes no seed. Its knobs are `diffusion_samples` (FINDING 009: does
not diversify the ligand), `num_recycles`, `num_steps`, `templates`. So replicate count is
not a lever and never was; **`num_recycles` / `num_steps` are the untested ones**, and
they change the trajectory rather than resample it.

## The cost, stated plainly

About 3,800 jobs and 5,892 poses (~1.8 GB) bought nothing. The check that would have
caught it costs one command — hash the files, or count distinct values per pair — and it
is the same check FINDING 009 and FINDING 011 already demanded for the *reference* set
("39% of replicates were duplicates", "count distinct poses, never jobs"). The rule was
written down and applied to the reference while the pool went unchecked.

**Dedupe the pool before believing any depth number.** Not the reference only. The pool.

---

## Addendum, same day — what duplicates do to the measured gain

Matched comparison: the same 342 pairs, the same references, only duplicates removed.

| arm | poses/pair | random | selected | oracle | gain |
|---|---|---|---|---|---|
| with duplicates | 25.13 | 0.7129 | **0.7436** | 0.7677 | +0.0306 |
| deduplicated | 3.62 | 0.7022 | **0.7436** | 0.7673 | **+0.0414** |

`selected` is **identical to four decimals** and the oracle is unchanged. The whole
+0.0108 comes from the *null* moving, because duplicating the pose the engine
deterministically returns drags a random draw toward it. Duplicate fraction: 85.6%.

**Do not read this as "the selector got better".** It did not move at all. The two
numbers answer different questions, and both are defensible:

- **+0.0306** is the gain over *pick a random file from what you generated*. On drop day
  the duplicates are real — that is genuinely what a naive baseline would hand you.
- **+0.0414** is the gain *per distinct opinion*, which is the honest measure of the
  selector's discriminative power.

The thing that matters for the leaderboard is unchanged either way: **selected = 0.7436**.
No submission gets better because we chose a different denominator. So `--dedupe-pool`
stays opt-in rather than becoming the default — silently switching it would make the
seventeen prior measurements incomparable while improving nothing.

What it does change is the read on FINDING 011/012: their +0.0381 and +0.0357 were
measured against nulls that included duplicates, so the selector's *discriminative power*
is stronger than those figures imply, even though the poses it picks are the same ones.

---

## Second addendum, same day — esmfold2 is deterministic too. Replicates are dead.

The reference-depth question needed a within-pair test, so esmfold2 was deepened 6 → 14
replicates. Applying this finding's own rule *while the batch ran* rather than after:

```
pairs with new esmfold2 replicates on disk: 469
  reps 0-5 : median 3 distinct of 6
  reps 6+  : median 1 distinct of 2
  NEW poses not already present: median 0.0, MEAN 0.00
  pairs where the new reps added NOTHING: 468 of 469
```

Zero. Killed at 652 jobs instead of 1,200.

**Both engines now return a deterministic pose per input.** The ~50% duplication in
esmfold2 reps 0-5 and the 3-4-distinct-of-12 in protenix_v2 reps 0-11 are *historical* —
from when the service behaved differently. Nothing submitted today diversifies at all.

### What this means, and it is not small

**Pose depth on OpenProtein is frozen at what we already have** — about 4 distinct
protenix_v2 poses and 3 distinct esmfold2 poses per pair. It cannot be increased by
buying replicates, from either engine, at any count. The `--replicates` lever that this
whole campaign was built on has stopped working.

Remaining levers, none of them yet measured:

1. **`num_recycles` / `num_steps`** — change the trajectory rather than resampling it.
   Cheapest to test and the obvious next probe.
2. **Other engines** — rosettafold_3, boltz2, boltz_1x run here. RF3 was already
   half-deterministic; the Boltz pair is untested for diversity on this venue.
3. **MSA variation** — subsample the 6,979-sequence MSA to different depths. Changes the
   input, so determinism of the engine does not prevent diversity of the output.
4. **Modal-hosted Boltz** — genuinely stochastic, but over its spend cap.

Until one of those is shown to work, **stop buying replicates**. They are pure cost.

### The pattern worth naming

Three times now the same check has decided the outcome, and each time it was applied one
step later than it should have been: FINDING 009 (samples do not diversify), FINDING 011
(39% of reference replicates were duplicates), and this one twice over — the pool, then
the reference. The check is always the same and always cheap: **hash the outputs and
count distinct ones before believing a depth number.** The cost of skipping it here was
~3,800 protenix_v2 jobs plus ~650 esmfold2 jobs. The cost of running it is one command.

---

## Third addendum — `num_recycles` / `num_steps` DO diversify. The lever is not dead.

Replicates are dead; the sampler knobs are not. Four pairs, five settings
(recycles × steps: 3×200 default, 1×200, 10×200, 3×50, 3×400):

```
4BJK_18I   5 distinct, max atom displacement 18.570 A
4FDH_0T3   5 distinct, max atom displacement  2.056 A
8EWN_X1O   5 distinct, max atom displacement 10.068 A
8EXB_X4E   5 distinct, max atom displacement 11.368 A
```

Five distinct ligand placements from five settings, on every pair, nothing unretrievable.

### And the diversity is not degradation

This is the check that disqualified single-sequence mode, which produced plenty of
"diversity" at 75% catastrophic poses and an oracle of 0.0505. Scored against crystal
ground truth:

| setting | mean LDDT-PLI | min | max |
|---|---|---|---|
| 3×200 (default) | 0.5997 | 0.4002 | 0.9875 |
| 3×50 | 0.6047 | 0.4018 | 0.9882 |
| 10×200 | 0.5972 | 0.3969 | 0.9875 |
| 3×400 | 0.5898 | 0.3952 | 0.9882 |
| 1×200 | 0.5648 | 0.3449 | 0.9770 |

**Catastrophic poses: 0.0%.** Every setting but `num_recycles=1` matches the default's
quality. Five-setting oracle 0.6095 against the default's 0.5997, **+0.0098**.

**n = 4 pairs and 20 poses. The +0.0098 is not significant at that size** and is quoted
only to say the direction is right; a 100-pair run is queued to measure it properly.

One thing to watch: 4BJK_18I shows 18.6 Å of atom displacement while its LDDT-PLI moves
only 0.345 → 0.402. Large coordinate spread with small score spread — either the
differences sit where LDDT-PLI is insensitive, or the ligand has a symmetry the heme-frame
comparison is treating as motion. Worth understanding before trusting displacement as a
diversity proxy.

### Corrects this finding's own earlier conclusion

The second addendum said "no replicate-style lever remains on this venue" as an open
question; the probe's first run then *printed* that as a verdict from one pair with one
setting, after 19 of 20 jobs failed retrieval against a loaded queue. That was wrong, and
it was the dangerous direction to be wrong in — it argued for abandoning the last
untested lever. Re-run on a quiet queue: 4 of 4 pairs, 0 failures, unambiguous.

**Depth is purchasable again — through the sampler settings, not through replicates.**
