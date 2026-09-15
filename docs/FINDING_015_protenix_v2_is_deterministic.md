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
