# FINDING 038 — reference depth is deterministic in the sweep size; the only thing that removes it is a failed job

**Status: measured.** 2026-09-23. 9 validation ligands × 4 sampler settings × 2 Protenix
checkpoints on OpenProtein, through the shipped `submit --sweep` / `collect` / `refset`
commands. **16 jobs, 25.7 minutes wall clock, $0.**

This closes the volume half of G6. It is a measurement, not a feature: no code path
changed, no selector was touched, and nothing here was tuned after the numbers arrived.

---

## 1. What was asked, and what it was asked against

The shipped selector (`cypstruct.xengine.select`, FINDING 011) scores a Boltz pose by its
mean Chamfer distance, in the heme frame, to poses of the **same** ligand from independent
engines. A blind test set has no such reference poses; they must be generated. Below **4**
independent poses the feature measured **−0.0055** — it makes selection *worse*, not
weaker — so `refset` refuses below depth 4 and names the ligands.

Before today the only evidence that a sweep clears depth 4 was the `g6probe` run:
**one ligand, two settings**. So "which ligands fall below depth 4" was unknowable, and the
release schedule could not be planned on it. The test ligands do not exist yet, so the
answer we need is a **rate**, measured on ligands we do have.

---

## 2. The ligand rule, fixed before any job was submitted

Stated first, then applied. The script is reproduced in
`data/processed/refdepth_ligands.csv`; no LDDT-PLI, oracle, gain or selection column is
read anywhere in it.

1. For each of the 87 validation ligands compute, **from the SMILES alone**,
   `n_heavy` (RDKit `GetNumHeavyAtoms`) and `n_rot` (`CalcNumRotatableBonds`).
2. Binding mode is the **prediction-side** label (FINDING 033: median `fe_donor_dist` over
   that ligand's own pool poses, cut at 2.6 Å), stored as `pred_mode`. It needs no crystal
   and is not a score.
3. Split the 87 into heavy-atom **terciles**: small `[8, 35]` n=30, medium `[36, 40]` n=32,
   large `[41, 54]` n=25.
4. **a — size × flexibility arm.** From each tercile take the Type II ligand with the
   *fewest* rotatable bonds and the one with the *most*. (6 ligands)
   **b — binding-mode arm.** From the 13 Type I (non-coordinating) ligands take the
   smallest, the median and the largest by `n_heavy`. (3 ligands)
   Every tie breaks on ligand id ascending.

**Amendment, recorded before a single job was submitted.** The first draft of 4b asked for
one Type I *per tercile*. The medium tercile contains **zero** Type I ligands (small 9 /
medium 0 / large 4 against 21 / 31 / 21 Type II, plus one peripheral), so that rule is
unsatisfiable on this data. 4b above replaces it and spans Type I over its own size range.
No score was read in making the change.

| ligand | slot | pred_mode | n_heavy | n_rot | MW |
|---|---|---|---|---|---|
| A1ASV | small / Type II / min rot | type_II | 22 | 2 | 311.3 |
| 8AW | small / Type II / max rot | type_II | 35 | 13 | 493.7 |
| KLN | medium / Type II / min rot | type_II | 36 | 7 | 531.4 |
| T0K | medium / Type II / max rot | type_II | 40 | 15 | 552.7 |
| QDJ | large / Type II / min rot | type_II | 44 | 11 | 605.8 |
| A1APA | large / Type II / max rot | type_II | 54 | 19 | 776.0 |
| PG0 | Type I / min heavy | **type_I** | 8 | 5 | 120.1 |
| MWS | Type I / median heavy | **type_I** | 27 | 7 | 372.4 |
| ERY | Type I / max heavy | **type_I** | 51 | 7 | 733.9 |

Span: heavy atoms **8 → 54**, rotatable bonds **2 → 19**, 6 Type II / 3 Type I.

---

## 3. What was run

```bash
python scripts/cofold/openprotein_cofold.py submit --engine protenix_v2 \
    --csv data/processed/refdepth_ligands.csv --tag refdepth \
    --sweep 3x200,10x200,3x50,3x400 --batch 5 --samples 5
python scripts/cofold/openprotein_cofold.py submit --engine protenix    ... (identical)
python scripts/cofold/openprotein_cofold.py collect --engine protenix_v2 --tag refdepth
python scripts/cofold/openprotein_cofold.py collect --engine protenix    --tag refdepth
python scripts/cofold/openprotein_cofold.py refset  --tag refdepth \
    --engines protenix_v2 protenix --csv data/processed/refdepth_ligands.csv \
    --out data/processed/refdepth_reference_set.npz
```

**One wave per setting**, per the FINDING 034 trap — distinct poses track submission waves,
not jobs, so `--replicates` is dead and was not used. `num_recycles < 2` is refused by
`parse_sweep` and was not attempted.

**Budget.** `cypstruct.budget.preflight('openprotein_cofold', 16, venue='openprotein',
cap_key='openprotein_jobs')` returned `ok=True`, est 0.32, and the launch proceeded. It was
not routed around. **Job count: 16 used.** `data/processed/openprotein/jobs.json` now holds
**610** OpenProtein jobs for this project against the **2,000/month** cap — 30.5% used.
*(Caveat, unchanged from the Modal incident: the ledger's `spent()` for OpenProtein reads
0.32 GPU-hours against a cap whose units are jobs. The honest count is the 610 above, read
from `jobs.json`, not from the ledger.)*

---

## 4. Wall clock — the driver is the slowest job, not the number of waves

All 8 waves were submitted within **1.3 min** of each other and ran **concurrently**: the
15 successful jobs sum to ~160 min of individual duration inside a 25.7 min wall clock,
i.e. **≈6.2× parallelism, with no queueing penalty observed at 16 jobs.**

Minutes from first submit (t₀ = 08:25:01Z):

| engine | setting | jobs | submitted | last terminal | wave duration | status |
|---|---|---|---|---|---|---|
| protenix_v2 | 3×200 | 2 | 0.0 | 10.1 | 10.1 | SUCCESS |
| protenix_v2 | 10×200 | 2 | 0.2 | 15.9 | 15.7 | SUCCESS |
| protenix_v2 | 3×50 | 2 | 0.3 | 17.4 | 17.1 | SUCCESS |
| protenix_v2 | **3×400** | 2 | 0.5 | **25.7** | 25.1 | **FAILURE + SUCCESS** |
| protenix | 3×200 | 2 | 0.8 | 7.3 | 6.5 | SUCCESS |
| protenix | 10×200 | 2 | 0.9 | 9.6 | 8.6 | SUCCESS |
| protenix | 3×50 | 2 | 1.1 | 10.7 | 9.6 | SUCCESS |
| protenix | 3×400 | 2 | 1.3 | 14.9 | 13.6 | SUCCESS |

* **Total: 25.7 min.** Last SUCCESS at **21.1 min**.
* Per-job duration over the 15 successes: **min 5.9 / median 10.1 / max 20.6 min**.
* `3×400` is the slowest setting on both engines and the only one that failed.
* `protenix` is uniformly faster than `protenix_v2` (6.5–13.6 vs 10.1–25.1 min per wave).
* `collect` + `refset` are negligible: `refset` reads all 340 files, dedupes and writes the
  `.npz` in **2.9 s**.

### The one job failure

`6e3f95dc…`, `protenix_v2` wave 3 (`3×400`), complexes `A1APA, PG0, MWS, ERY`:
`JobStatus.FAILURE`, `failure_message='internal server error'`, after **24.7 min** with
`progress_counter=75`. **1 of 16 jobs = 6.25%.** It is reported, not hidden, and it is the
single most consequential number in this document — see §7.

---

## 5. The deliverable: the depth distribution

### Filter counts — every filter, including the ones that fired zero times

| filter | count | what it removes |
|---|---|---|
| `files_seen` | **340** | 68 complexes × 5 diffusion samples |
| `unparseable_name` | **0** | — |
| `not_in_csv` | **0** | — |
| **`dropped_same_wave`** | **272** | FINDING 009: within one job every sample shares one ligand conformation — 4 of every 5 files were discarded unread |
| **`dropped_same_md5`** | **0** | FINDING 015/034: byte-identical replicates. **Fired zero times.** |
| **`dropped_same_coords`** | **0** | byte-different, geometrically identical poses at 0.05 Å in the heme frame. **Fired zero times.** |
| `unreadable` | **0** | — |
| `no_heme_frame` | **0** | — |

**A filter that never fires beats any passing check (T8), so both zero counts were proved
to be real zeros rather than dead code.** Positive controls, all passing:

| control | result |
|---|---|
| md5 filter on an exact byte copy | fires |
| coordinate filter on an identical vector | fires |
| coordinate filter on a geometrically identical pose (+1e-4 Å, byte-different) | fires |
| coordinate filter **keeps** a genuinely distinct pose (+5.0 Å) | keeps 2 |
| same-wave filter on duplicate samples | fires (272×) |

The zeros are the interesting result, not an absence of one: **every (engine, wave) pair
returned a genuinely distinct pose.** Because each wave is a *different sampler setting*,
FINDING 015's determinism — which is per-configuration — never gets a chance to bite.

### Per-ligand depth

| ligand | poses returned | after same-wave | after md5 | after coords | **depth** | min pairwise Chamfer | median | max |
|---|---|---|---|---|---|---|---|---|
| 8AW | 40 | 8 | 8 | 8 | **8** | 0.351 | 0.605 | 0.838 |
| A1APA | 35 | 7 | 7 | 7 | **7** | 0.395 | 0.871 | 1.038 |
| A1ASV | 40 | 8 | 8 | 8 | **8** | 0.288 | 0.528 | 0.957 |
| ERY | 35 | 7 | 7 | 7 | **7** | 0.486 | 0.954 | 1.191 |
| KLN | 40 | 8 | 8 | 8 | **8** | 0.672 | 1.032 | 1.383 |
| MWS | 35 | 7 | 7 | 7 | **7** | 0.342 | 0.659 | 0.876 |
| PG0 | 35 | 7 | 7 | 7 | **7** | 0.516 | 0.955 | **15.667** |
| QDJ | 40 | 8 | 8 | 8 | **8** | 0.440 | 0.933 | 1.646 |
| T0K | 40 | 8 | 8 | 8 | **8** | 0.403 | 1.016 | 1.627 |

### The rate

* **Clear depth 4: 9 of 9 = 100%.**
* **Clear depth 6: 9 of 9 = 100%.**
* **Median depth 8. Minimum 7. Maximum 8.**
* `refset` wrote the `.npz` — `below_min_depth: []`, `no_poses_at_all: []`, 68 poses,
  2,396 atoms, 26 kB.

### The mechanism, stated plainly

`poses_returned / 5 = after_same_wave = after_md5 = after_coords = depth` on **every one of
the nine ligands**. Therefore:

> **Depth = (number of engines) × (number of sampler settings), minus one for every failed
> job that contained the ligand.** There is no attrition from duplication at all.

Five ligands got 8 = 2 × 4. Four ligands got 7 = 2 × 4 − 1, and those four are *exactly*
`A1APA, PG0, MWS, ERY` — the contents of the one failed job. Nothing else varies.

The distinctness has margin: the **minimum** pairwise Chamfer anywhere in the reference set
is **0.288 Å**, against a dedupe tolerance of 0.05 Å — these poses are ~6× clear of the
threshold, not marginal passes. Median pairwise Chamfer across all 224 pairs is **0.831 Å**.
PG0 (8 heavy atoms, MW 120) carries one placement **15.7 Å** from its siblings — a genuine
catastrophic outlier in the reference set, which is what a Chamfer *mean* is designed to
absorb; it is recorded, not acted on.

---

## 6. Does depth depend on chemistry? No.

| relation | statistic | verdict |
|---|---|---|
| depth ~ n_heavy | ρ = **−0.087**, p = 0.82 | null |
| depth ~ n_rot | ρ = **+0.088**, p = 0.82 | null |
| depth ~ binding mode | Type I 7.00 vs Type II 7.83, MWU **p = 0.0369** | **artifact — see below** |
| pose spread ~ n_heavy | ρ = +0.217, p = 0.58 | null |
| pose spread ~ n_rot | ρ = +0.102, p = 0.79 | null |
| pose spread ~ binding mode | Type I 0.856 vs Type II 0.831 Å, MWU p = 0.90 | null |

**The p = 0.0369 is not chemistry. It is batch packing, and it is a trap worth naming.**

`--batch 5` packs the CSV **in file order**. The CSV is written by the rule in §2, which
emits the six Type II size/flexibility picks first and the three Type I picks last. So the
two chunks were:

```
chunk 1: A1ASV, 8AW, KLN, T0K, QDJ        -> 5 Type II
chunk 2: A1APA, PG0, MWS, ERY             -> 1 Type II + ALL THREE Type I
```

The one job that failed was chunk 2. **All three Type I ligands were in it**, so all three
lost a wave, and the Mann–Whitney picks that up as a mode effect at p = 0.037 on n = 3 vs
n = 6. Remove the failed job and Type I and Type II depths are *identical by construction*
(8 and 8).

So the answer to the question as posed — "if Type I ligands yield systematically less
depth, that compounds with FINDING 033" — is **they do not**. Type I ligands are not harder
to reference. But there is a real operational hazard here, and it is new:

> **An ordered CSV plus `--batch` packing aligns job boundaries with strata, so a single
> job failure can take out an entire stratum.** This was a 1-in-16 event that happened to
> land on 3 of 3 Type I ligands — the stratum FINDING 033 already identifies as having the
> worse pool. **Shuffle or interleave the ligand CSV before submitting**, so a failure
> spreads across chemistry instead of destroying one class of it.

---

## 7. The failure is the binding constraint, and `submit` cannot currently recover from it

Two facts, measured:

1. **1 of 16 jobs failed** with an opaque `internal server error` after 24.7 minutes.
2. **`submit` will not resubmit it.** Its resume set is
   `claimed = {(rep, ligand) for every batch}` with **no check of `b["done"]`**. All 4
   lost ligand-waves are still marked claimed, verified directly: `4 of 4`. Re-running the
   identical `submit --sweep` command prints "0 ligands" and buys nothing. The only
   recovery is a **fresh wave index at a new setting**, because `submit` also (correctly)
   refuses to reuse a wave index at a different setting.

This is the reason the recommended sweep in §9 carries margin rather than sitting on the
depth-4 line. A 6.25% per-job failure rate against a hard refusal at depth 4 is the whole
risk of drop day's longest pole.

---

## 8. The minimum viable sweep

Every sub-configuration of the poses already on disk, re-deduped from scratch
(`data/processed/refdepth_subsweep.csv`). Jobs are per 9 ligands at `--batch 5`.

**Two engines** (protenix_v2 + protenix):

| settings | n | jobs | min depth | median | frac ≥4 | frac ≥6 |
|---|---|---|---|---|---|---|
| any one setting | 1 | 4 | 1–2 | 2.0 | **0.00** | 0.00 |
| 3×200 + 10×200 | 2 | 8 | 4 | 4.0 | **1.00** | 0.00 |
| 3×200 + 3×50 | 2 | 8 | 4 | 4.0 | **1.00** | 0.00 |
| 10×200 + 3×50 | 2 | 8 | 4 | 4.0 | **1.00** | 0.00 |
| any 2 incl. 3×400 | 2 | 8 | 3 | 4.0 | 0.56 | 0.00 |
| **3×200 + 10×200 + 3×50** | **3** | **12** | **6** | **6.0** | **1.00** | **1.00** |
| any other 3 (incl. 3×400) | 3 | 12 | 5 | 6.0 | **1.00** | 0.56 |
| all four | 4 | 16 | 7 | 8.0 | **1.00** | 1.00 |

**One engine:**

| engine | settings | jobs | min depth | frac ≥4 | frac ≥6 |
|---|---|---|---|---|---|
| protenix | 3 | 6 | 3 | 0.00 | 0.00 |
| protenix | 4 | 8 | 4 | **1.00** | 0.00 |
| protenix_v2 | 3 | 6 | 2 | 0.00 | 0.00 |
| protenix_v2 | 4 | 8 | 3 | 0.56 | 0.00 |

Read it straight off the mechanism in §5: depth is engines × settings. Every `frac ≥4 =
0.56` above is the four ligands from the one failed job, and every `min depth` one below
median is the same four. **The 3×400 column looks weak because it is the column that
failed, not because the setting is weak** — its 5 surviving ligands are indistinguishable
from the other settings.

**Conclusions.**

* **The arithmetic minimum is 2 engines × 2 settings = 4 poses = 8 jobs per 9 ligands.**
  That is exactly the refusal threshold with **zero margin**: one failed job drops those
  ligands to 3 and `refset` refuses them.
* **The minimum viable sweep is 2 engines × 3 settings — `3x200,10x200,3x50`.** Depth 6 on
  every ligand, which survives one job failure per ligand (6 → 5 ≥ 4) and costs 12 jobs per
  9 ligands. It also **drops `3×400`, which is the slowest setting on both engines and the
  only one that failed** — pure wall clock removed along with the risk.
* **The fourth setting adds nothing but insurance.** It raises depth 6 → 8 at +33% jobs and
  +8 min of wall clock. Take it only if the release schedule has slack, or if the first
  `refset` names anyone.
* **A single engine cannot do the job.** Four settings on protenix alone lands *exactly* on
  depth 4 with no margin; on protenix_v2 alone it does not even reach the rate. This is
  consistent with FINDING 011's "best reference is the two Protenix checkpoints".

---

## 9. Cost model for the release schedule

At the recommended sweep (**2 engines × 3 settings = 6 waves**) and `--batch 5`:

`jobs = 6 × ceil(N / 5)`

| test set | jobs (3 settings) | jobs (4 settings) | disk at `--samples 5` | disk at `--samples 1` |
|---|---|---|---|---|
| 20 ligands | **24** | 32 | ~250 MB | ~50 MB |
| 50 ligands | **60** | 80 | ~620 MB | ~125 MB |
| 100 ligands | **120** | 160 | ~1.2 GB | ~250 MB |

All of these fit the 2,000/month cap with 610 already spent (100 ligands × 4 settings =
160 jobs → 770 of 2,000).

**Wall clock.** Concurrency was measured at **≥16 jobs with no queueing penalty**; above
that it is **not measured**, so the schedule is quoted as a band. The lower bound is the
slowest single wave (17.4 min at 3 settings, 25.7 min at 4); the upper bound assumes the
observed 16 is a hard ceiling and serialises in blocks of 16 at ~18 min each.

| test set | jobs | optimistic (full concurrency) | conservative (16 at a time) |
|---|---|---|---|
| 20 ligands | 24 | **~20 min** | ~36 min |
| 50 ligands | 60 | **~20 min** | ~72 min |
| 100 ligands | 120 | **~20 min** | ~2 h 25 m |

Add ~5 min for `collect` and **under 1 minute** for `refset` at any of these sizes
(2.9 s for 340 files, and it is linear).

**Even the conservative column is well inside a drop day.** G6 is the longest pole of the
reference-set step, and at 100 ligands it is under two and a half hours — far shorter than
Step 2's pool generation on Explorer (~6.6 min/ligand × 100 = 11 h). **Reference depth is
no longer the schedule risk.**

**Use `--samples 1`, not the default 20.** 272 of 340 downloaded files — **80%** — were
discarded unread by the same-wave filter, because `diffusion_samples` does not sample the
ligand (FINDING 009) and `refset` keeps exactly one file per (engine, wave). At
`--samples 20` on 100 ligands that is ~5 GB downloaded to keep ~250 MB, onto a D: drive
with 12 GB free. The confidence block is the only thing lost, and Protenix confidence ranks
poses at chance (FINDING 009).

---

## 10. Recommendation for the release schedule

1. **Shuffle the test CSV before submitting.** One line. It prevents a single job failure
   from destroying a whole stratum (§6), which is what happened here.
2. **Submit `--sweep 3x200,10x200,3x50` on both Protenix checkpoints, `--batch 5`,
   `--samples 1`.** 6 waves, `6 × ceil(N/5)` jobs, depth 6 on every ligand, ~20 min to
   ~2 h 25 m wall clock for N = 20 to 100.
3. **Run `refset` and read `below_min_depth` before committing to anything.** At depth 6 it
   should be empty. If it is not, the cause will be a failed job, not chemistry — check
   `jobs.json` for `"done": "failed"` first.
4. **Recover a failed job by adding a FRESH wave at a new setting** (`--sweep 3x400`, or
   `5x200`), not by re-running the same command. `submit`'s resume set marks failed batches
   as claimed and will skip them (§7). Fixing that is a small change to
   `openprotein_cofold.submit` and is the obvious follow-up to this finding; it was not
   made here because this is a measurement.
5. **Budget the fourth setting as slack, not as baseline.** +33% jobs, +8 min, depth 6 → 8.
6. **Do not add a third engine.** FINDING 011: two Protenix checkpoints measure +0.0380;
   adding esmfold2 at matched depth drops it to +0.0178. Depth is not the constraint here —
   it is now free.

---

## 11. What this does and does not establish

**Establishes.** On 9 CYP3A4 ligands spanning 8–54 heavy atoms, 2–19 rotatable bonds and
both binding modes, a 4-setting × 2-engine sweep yields **8 distinct reference poses per
ligand with zero duplicate attrition**, and depth is **deterministic in the sweep size**.
100% clear depth 4; 100% clear depth 6. Chemistry does not predict depth (all p > 0.5 once
the failed job is accounted for). The whole thing cost 16 jobs and 25.7 minutes.

**Does not establish.** n = 9 ligands, all CYP3A4, all already known to fold. The blind set
may contain chemistry that fails to *parse* (organometallics — Step 1's gate) or that the
engine rejects; that is a different failure than thin depth and this says nothing about it.
Job concurrency above 16 is unmeasured. The per-job failure rate is **1 observation** —
6.25% is a point estimate with a 95% interval of roughly 0.2–30%, which is precisely why
the recommendation carries margin rather than sitting on the threshold.

**Explicitly not claimed.** That these 8 reference poses *select better* than 4. This
finding measures supply, not value. FINDING 011 measured the value of depth ≥ 4 and found
+0.0381; whether 8 beats 4 is untested, and the §8 recommendation is framed as
failure-tolerance, not as a score claim.

---

## Artifacts

* `data/processed/refdepth_ligands.csv` — the 9 picks with the slot each one fills
* `data/processed/refdepth_chemistry_all87.csv` — n_heavy / n_rot / mode for all 87
* `data/processed/refdepth_per_ligand.csv` — the depth table of §5
* `data/processed/refdepth_subsweep.csv` — every sub-configuration of §8
* `data/processed/refdepth_summary.json` — filter counts, positive controls, rate, chemistry
* `data/processed/refdepth_job_timeline.json` — per-job submit/terminal timestamps
* `data/processed/refdepth_reference_set.npz` — the frozen reference set (9 ligands, 68 poses)
* `data/processed/openprotein/refdepth/` — 340 mmCIFs, **gitignored** (115 MB)
