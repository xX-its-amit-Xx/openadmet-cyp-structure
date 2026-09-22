# PRE-REGISTRATION — what a 15-ligand structure test set can and cannot tell us

**Written 2026-09-21, BEFORE the test set is released.** `readiness.py` reports
`STRUCTURE_DATASET_SIZE = 15` (the 184 was a PXR placeholder) with
`STRUCTURE_TRACK_LIVE = False`. Interim deadline 2026-09-24.

This exists so that the number we get back is interpreted against an expectation set in
advance, rather than explained after the fact.

## The measurement

Per-ligand gain of the shipped selector over random selection, on our 87-ligand pool:
**mean +0.0383, standard deviation across ligands 0.0944.** The per-ligand spread is
**2.5x the mean effect** — most of the selector's value comes from a minority of ligands
where it avoids a bad pose, which is what a catastrophe detector looks like.

Bootstrapping test sets of each size from those 87 ligands, 20,000 draws:

| test set n | mean gain | sd | P(gain > 0) | 5th pct | 95th pct |
|---|---|---|---|---|---|
| **15** | **+0.0385** | **0.0243** | **0.959** | **+0.0017** | **+0.0814** |
| 30 | +0.0384 | 0.0171 | 0.994 | +0.0118 | +0.0677 |
| 50 | +0.0381 | 0.0133 | 0.999 | +0.0172 | +0.0609 |
| 87 | +0.0382 | 0.0100 | 1.000 | +0.0224 | +0.0553 |

## What this means for the interim result

**The selector will almost certainly help: P(gain > 0) = 0.96 at n = 15.**

**The magnitude is close to uninformative.** The 90% interval at n = 15 runs from
**+0.0017 to +0.0814** — from "indistinguishable from random" to "more than double the
measured effect". A single 15-ligand score cannot separate those.

So, pre-registered:

1. **A weak interim score is not evidence the selector failed.** One draw in twenty lands
   at essentially zero gain purely by ligand sampling.
2. **A strong interim score is not evidence it is better than +0.038.** The top of the
   interval is reachable by luck alone.
3. **Leaderboard position at n = 15 is substantially noise** for every entrant, not only
   us. Differences between neighbouring entries should not be read as method differences.
4. **Nothing gets retuned on the basis of the interim result.** The selector has replicated
   at +0.0395 / +0.0357 / +0.0383 across three pools of 87, 81 and 87 units; a 15-ligand
   sample is the weakest evidence we hold, and it would be the worst thing to fit to.

## The one thing worth watching instead

`pool_diagnostics.py` reports the pool's catastrophe rate before submission. The selector's
gain tracks that rate (FINDING 011: +0.30 / +0.07 / +0.006 across regimes), so the
diagnostic predicts the regime *before* we see a score. If the released pool is
well-behaved, expect the low end — and expect it for a reason we can state in advance
rather than discover afterwards.


---

## Update 2026-09-22 — the test set grew from 15 to 20

`readiness.py` now reports `STRUCTURE_DATASET_SIZE = 20`, confirmed on two consecutive
reads. `STRUCTURE_TRACK_LIVE` is still `False` and `dataset last modified` is unchanged at
2026-08-27, so the size is being read from configuration rather than from file mtime —
the entry count moved without the files being touched. **First upstream change in 21
hourly polls.**

Re-running the same bootstrap at the new size:

| test set n | mean gain | sd | P(gain > 0) | 5th pct | 95th pct |
|---|---|---|---|---|---|
| 15 (previous) | +0.0385 | 0.0243 | 0.959 | +0.0017 | +0.0814 |
| **20 (current)** | **+0.0386** | **0.0211** | **0.978** | **+0.0065** | **+0.0751** |
| 30 | +0.0380 | 0.0172 | 0.993 | +0.0112 | +0.0675 |

Five extra ligands move P(gain > 0) from 0.959 to **0.978** and lift the 5th percentile
off zero, from +0.0017 to +0.0065. Better, and **not enough to change any of the four
pre-registered conclusions** — the 90% interval is still [+0.0065, +0.0751], roughly an
eleven-fold spread. A single 20-ligand score still cannot distinguish our measured
+0.038 from half or double it.

The pre-registration stands unchanged.
