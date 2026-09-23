# FINDING 041 — reference depth 4 → 8 is worth **+0.0076** on CYP3A4, and it costs **jobs, not wall clock**

**Status: measured.** 2026-09-23. 37 CYP3A4 validation ligands, **370 new reference
poses** from 80 OpenProtein jobs (**0 failures**), **38.3 minutes** wall clock, $0.
Pre-registered in `docs/PREREG_reference_depth_value.md`, committed at `40e0f55` **before
a single job was submitted** and not edited afterwards.

The depth curve below is **exact** — every one of the `C(10, d)` reference subsets is
enumerated. The only sampled numbers are the noise floor (2 × 10⁶ draws, SE reported) and
the across-draw percentiles (2 × 10⁵ draws), and both say so where they appear.

---

## What was asked

Two findings gave opposite release-day advice about reference depth, and both cannot be
right for the same hour of wall clock:

* **FINDING 016** — *"Reference depth saturates at 4. Buy exactly 4 settings. Eight costs
  twice the jobs for a difference indistinguishable from zero."* Measured within-pair on
  **487 P450 pairs**: **+0.0315** at depth 4 against **+0.0345** at depth 8, Δ **+0.0030**.
* **FINDING 040 §7** — *"Spend the wall clock on reference depth instead. Depth in the pool
  is worth +0.0027; depth in the reference is what makes the +0.0395 gain exist at all."*

**FINDING 038** measured that depth 8 is cheaply *achievable* on CYP3A4 and said in its own
§11 that it **explicitly does not claim** 8 reference poses select better than 4 — *"this
finding measures supply, not value."*

CYP3A4 is a known family outlier (FINDING 022: the P450 family is solved at 0.86 / 89%
sub-2 Å, CYP3A4 alone sits at 0.555 / 33%), so FINDING 016's family answer is not
automatically CYP3A4's. **This measures the value on CYP3A4.**

---

## 1. What was run

**Ligands — the rule, fixed in §2 of the pre-registration and applied without reading a
score.** All **13 prediction-side Type I** ligands (the scarce stratum, FINDING 033's named
exposure) plus **8 Type II per heavy-atom tercile** drawn with
`default_rng(41).choice(sorted(ids), 8, replace=False)`; the single `peripheral` ligand
dropped. **n = 37**, heavy atoms **8 → 51**, rotatable bonds **0 → 17**, 13 Type I / 24
Type II across all three size terciles. `data/processed/refvalue_ligands.csv`.

**Reference — ten waves, five settings, two checkpoints.**

```bash
python scripts/cofold/openprotein_cofold.py submit --engine protenix_v2 \
    --csv data/processed/refvalue_ligands.csv --tag refvalue \
    --sweep 3x200,10x200,3x50,5x200,3x150 --batch 5 --samples 1
python scripts/cofold/openprotein_cofold.py submit --engine protenix     ... (identical)
#   collect (poll until terminal)
python scripts/structure/refvalue_refset.py 8        # refset, report saved
```

**Five settings, not four, and the reason is a methods trap, not a budget.** FINDING 035
established that *a subsample curve's last rung is a deterministic maximum, so its final
slope over-prices new samples.* Had available depth been exactly 8 there would be exactly
one subset of size 8, and the depth-8 point would be a **maximum** where depths 4–7 are
**averages** — a bias pointing straight at the conclusion under test. At available depth
**10** every rung from 4 to 8 is a genuine average over `C(10, d) ≥ 45` subsets.
`3x400` was dropped (slowest setting on both engines in FINDING 038, and the only one that
failed); `5x200` and `3x150` come from FINDING 016's second four.

**Pool — held completely fixed.** The existing **20-pose unsteered Boltz-2** pool, already
scored (`data/processed/poses_scored_val87b.csv`). 740 poses. Nothing was added, removed
or re-scored. The only thing that varies across the whole primary endpoint is **which
reference poses the feature is allowed to see.**

**Budget.** `preflight('openprotein_cofold', 80, venue='openprotein',
cap_key='openprotein_jobs')` → `ok=True`, est 1.60. It was not routed around. **80 jobs
used**; `jobs.json` goes **610 → 690** of the **2,000/month** cap (**34.5%**). The honest
count is read from `jobs.json`, not from the ledger, whose `spent()` reports GPU-hours
against a cap whose units are jobs (the caveat FINDING 038 recorded, unchanged).

---

## 2. Supply — depth 10 on 37 of 37, with every filter count including the zeros

| filter | count | what it removes |
|---|---|---|
| `files_seen` | **370** | 37 ligands × 10 waves × 1 sample |
| `unparseable_name` | 0 | — |
| `not_in_csv` | 0 | — |
| `dropped_same_wave` | **0** | FINDING 009. **Zero here because `--samples 1` is now the default** — FINDING 038 discarded 272 of 340 files (80%) unread at `--samples 5`. The fix works. |
| `dropped_same_md5` | **0** | FINDING 015/034 byte-identical replicates |
| `dropped_same_coords` | **0** | geometrically identical at 0.05 Å in the heme frame |
| `unreadable` | 0 | — |
| `no_heme_frame` | 0 | — |

**Reference depth: min 10, median 10, max 10, on all 37 ligands.** `below_min_depth: []`,
`no_poses_at_all: []`, 370 poses, 12,530 atoms, 136 kB.
**Depth = engines × settings** (FINDING 038's arithmetic) with **no subtraction, because
nothing failed: 0 of 80 jobs**, against FINDING 038's 1 of 16.

**A filter that never fires beats any passing check**, so all five zeros were proved real
(`refvalue_filter_controls.json`, all pass):

| positive control | result |
|---|---|
| md5 on an exact byte copy | fires |
| coordinate filter on an identical vector | fires |
| coordinate filter on a geometrically identical pose (+1e-4 Å) | fires |
| coordinate filter **keeps** a pose 5.0 Å away | keeps 2 |
| coordinate filter **keeps** two real reference poses of the same ligand | keeps 2 |

**Distinctness margin.** Median pairwise Chamfer inside a ligand's reference set is
**0.998 Å**; the **tightest pair anywhere is 0.091 Å** (CFF), against a 0.05 Å tolerance.
That is **1.8×** clear, where FINDING 038 was 5.8× clear at 0.288 Å — a real narrowing,
worth naming: at five settings the two closest sampler settings do start to converge on
some ligands, and a sixth setting would be the one to check rather than assume.

**Waves, not jobs (C6).** 80 jobs produced **370 distinct poses** = 37 × 10 waves. The two
numbers are never conflated anywhere below.

---

## 3. Controls — all pass, before any curve is believed

| # | control | result |
|---|---|---|
| **C1** | numbering (FINDING 021) | **0 of 37** ligands read LDDT-PLI identically 0.0. Four individual poses are exactly 0.0 — PG0 m19 and PG4 m17/18/19 — and all four carry `bisy_rmsd` **23.1–26.9 Å** with `valid=False`: genuine ejections from the site, not a numbering offset. Range 0.0 – 0.908. |
| **C2** | reproduce the shipped `xeng` column | 740 of 740 poses recomputed from `reference_set_cyp3a4.npz`, **max abs diff 2.45e-07** (the npz stores float32), mean 1.44e-08. |
| **C3** | `argmin(xeng)` **is** `xengine.select()` | **37 of 37** identical picks. |
| **C4/C5** | `refset` filters + distinctness | §2 above |
| **C6** | waves ≠ jobs | §2 above |

**The answer-recognition gate was NOT run**, per FINDING 036: R1 applies to *priors*, not
to within-ligand comparators, and the shipped selector fails it by construction while
working.

---

## 4. The deliverable — selection as a function of reference depth

Pool, reported first and unchanged at every depth: **oracle 0.70182**, **random 0.56094**,
740 poses over 37 ligands.

**Noise floor, recomputed on THIS population** (37 ligands × 20 poses, 2 × 10⁶ draws):
**p95 = +0.02272** (SE **±0.000023**), p99 = +0.03209, sd of gain 0.01412. Nothing is
quoted from another finding — FINDING 039 established that a floor is a number *about a
population*, and this one is 37 ligands, so it sits well above FINDING 040's n=73 value of
+0.01383 and well below its own n=14 values.

| reference depth | subsets/ligand | **selected (exact)** | **gain over random** | clears floor | **sd across draws** | p5 | p95 | mean within-ligand ρ | correct sign |
|---|---|---|---|---|---|---|---|---|---|
| **4** | 210 | 0.59993 | **+0.03899** | yes | **0.01215** | 0.57953 | 0.61961 | −0.2613 | 83.8% |
| **5** | 252 | 0.60189 | **+0.04095** | yes | 0.01147 | 0.58232 | 0.62017 | −0.2676 | 83.8% |
| **6** | 210 | 0.60420 | **+0.04326** | yes | 0.01031 | 0.58637 | 0.62025 | −0.2661 | 83.8% |
| **7** | 120 | 0.60698 | **+0.04604** | yes | 0.00892 | 0.59115 | 0.62066 | −0.2749 | 83.8% |
| **8** | 45 | 0.60756 | **+0.04662** | yes | **0.00812** | 0.59317 | 0.62006 | −0.2782 | 83.8% |
| *10 (all available)* | 1 | *0.61271* | *+0.05177* | yes | *0* | — | — | −0.2713 | 83.8% |

The depth-10 row is the **deterministic maximum** the pre-registration was designed around
and is printed in italics for the **A0 gate only** — it is a single subset and its slope
must not be read as a rate. **A0 passes**: the new sweep reference at full depth scores
**+0.0518** against this population's floor of **+0.0227**, so the instrument works before
any depth claim is made.

Increments of the exact gain: **+0.00196 / +0.00231 / +0.00277 / +0.00058**
(4→5, 5→6, 6→7, 7→8). **4 → 6 is +0.00427; 6 → 8 is +0.00336; 4 → 8 is +0.00763.**

### The across-draw spread — this is the part that is not in the mean

A "draw" is one independent reference subset per ligand, which is exactly what a blind set
gets: **one**. The pooled mean's distribution over 200,000 such draws:

| depth | sd | min | p5 | p95 | max |
|---|---|---|---|---|---|
| 4 | **0.01215** | 0.55129 | **0.57953** | 0.61961 | 0.64778 |
| 6 | 0.01031 | 0.55767 | 0.58637 | 0.62025 | 0.64113 |
| 8 | **0.00812** | 0.56540 | **0.59317** | 0.62006 | 0.62927 |

**Which four you pick matters more than whether you pick four or eight.** The depth-4
p5-to-p95 range is **0.0401 LDDT-PLI wide — 5.3× the whole 4 → 8 effect of 0.0076.** And
the shape of what depth buys is asymmetric and unambiguous:

* **p5 rises +0.0136** (0.5795 → 0.5932),
* **p95 rises +0.0004** (0.6196 → 0.6201),
* **sd falls 33%** (0.01215 → 0.00812), **worst draw rises +0.0141**, **best draw falls
  0.0185**.

> **Reference depth does not make the selector better. It makes a bad draw less likely.**
> That is the same character FINDING 012 gave the selector itself — a catastrophe
> detector, not a ranker — reappearing one level up, in the reference set.

### Ranking does not improve, and the clean number is real

Mean within-ligand ρ moves **−0.2613 → −0.2782** across the whole doubling, and the
correct-sign fraction is **83.8% (31 of 37) at every single depth**. A figure that
identical at five depths is exactly the "too-clean number" this repo treats as a tell, so
it was checked independently with a different seed and 200 subsets per ligand: it is the
**same six ligands** (`A1CIW, MF8, MWY, NI6, TCI, YNV`) that carry ρ ≥ 0 at depth 4 and at
depth 8. The number is real, and what it says is that **depth does not convert any ligand
from "the feature ranks it backwards" to "the feature ranks it right."**

---

## 5. The verdict against the pre-registered clauses

Primary quantity **Δ = selected(8) − selected(4) = +0.00763**, exact, same 37 ligands,
same 740 poses.

| clause | required | measured | fires? |
|---|---|---|---|
| **A0** (gate) | full-depth gain clears this population's floor | **+0.0518 vs +0.0227** | **passed** |
| **T1 — magnitude** | Δ ≥ **+0.0125** (FINDING 004's per-doubling bar) | **+0.00763** | **NO** |
| **T2 — inference** | bootstrap CI excludes 0 | **[−0.0021, +0.0190]**, 93.1% of 10,000 draws positive, Wilcoxon **p = 0.261** | **NO** |
| **T3 — not one ligand** | sign survives leave-one-out | drop-1 range **[+0.0042, +0.0093]**, **sign stable on all 37** | **yes** |
| **T4 — monotone** | selected(d) non-decreasing, d = 4…8 | **strictly increasing at all four steps** | **yes** |

> ## **VERDICT: "detectable, not worth the wall clock" — the pre-registration's middle band.**
>
> `0 < Δ < +0.0125`. Reference depth beyond 4 **does** pay on CYP3A4, monotonically and
> with a sign that survives dropping any ligand, and it pays **+0.0076 per doubling**
> against the **+0.0125** that would make wall clock worth spending on it.

**And then the cost turned out not to be wall clock at all — see §7.** That is a fact
measured in this run, not a change to the acceptance rule, and it is reported separately
because the rule was fixed before the cost was known.

### Leverage — reported whichever way it lands

| | value |
|---|---|
| Δ, all 37 | **+0.00763** |
| drop-1 minimum (drops **CL6**) | **+0.00419** |
| drop-1 maximum | +0.00934 |
| ligands whose **modal pick** changes between depth 4 and 8 | **2 of 37** (`80K`, `A1AST`) |
| ligands that **cannot move at all** (Δ exactly 0, disagreement 0) | **3 of 37** (`ERY`, `X6V`, `X7M`) |
| mean fraction of (depth-4 subset, depth-8 subset) pairs that **disagree on the pick** | **33.2%** (median 37.1%, max 73.4%) |

**This is the honest split, and it is the mirror image of FINDING 040's.** There, the pool
depth point estimate **flipped sign** when one ligand was removed (+0.0027 → −0.0012).
Here the estimate is **sign-stable across all 37 leave-one-outs** — the strongest thing
that can be said for it — while being built out of enormous cancelling per-ligand swings:
`CL6` **+0.132**, `D0R` **+0.089**, `08J` **+0.061** against `CFF` **−0.054**, `PG0`
**−0.044**, `X6J` **−0.025**. Only **two ligands** change their modal pick, yet **a third
of all subset pairs disagree**: the churn is huge and almost all of it cancels.

**The oracle-to-selection gap narrows**, which no pool-depth purchase in this repo has ever
managed: **0.1019 at depth 4 → 0.0943 at depth 8 → 0.0891 at depth 10.** FINDING 040's
pool doubling *widened* it 0.0770 → 0.0876. Reference depth spends nothing on the ceiling
and takes all of its (small) gain out of the gap.

### A second reference source, descriptive, not pre-registered

The same 37 ligands and the same 740 poses, scored against the **shipped** reference set
(`reference_set_cyp3a4.npz` — engine *replicates*, not a sampler sweep; depth 6–11 on these
ligands, median 7):

| depth | 4 | 5 | 6 | full (6–11) |
|---|---|---|---|---|
| gain | +0.0396 | +0.0420 | **+0.0391** | +0.0428 |

**Not monotone, total range 0.003 over 4→6** — and the sweep arm over the same 4→6 range is
**+0.0043**. So the two reference sources agree that **nothing much happens below depth 6**,
and the sweep arm's whole effect is loaded into **6 → 7**. This is a second, independent
reason not to over-read +0.0076.

---

## 6. Reconciling FINDING 016 and FINDING 040 — which advice survives

**They do not actually disagree about the measurement. They disagree about the price, and
the price changed.**

| | FINDING 016 | **FINDING 041 (this)** |
|---|---|---|
| population | 487 P450 pairs, 81 proteins | **37 CYP3A4 ligands** |
| reference | 4 vs 8 sampler settings | 4 vs 8, subsampled from 10 |
| gain at depth 4 | +0.0315 | **+0.0390** |
| gain at depth 8 | +0.0345 | **+0.0466** |
| **Δ** | **+0.0030** | **+0.0076** |
| verdict as pre-registered | "not measurably better" | "detectable, below the bar" |

**Both are right on their own populations, and they are the same answer to two decimal
places of importance.** CYP3A4's Δ is 2.5× the family's, which is what FINDING 022 would
predict — the outlier target has more headroom to recover — but **both are small, and
neither clears +0.0125.** The disagreement was never about the size of the effect.

**What FINDING 016 got right, and what does not survive.** Its *measurement* replicates:
reference depth is flat-to-slightly-positive past 4, and the depth-4 refusal threshold is
a knee rather than a lower bound. Its *recommendation* — "buy exactly 4; eight costs twice
the jobs" — rested on jobs and wall clock being the same currency. **They are not.**
Measured here on 80 concurrent jobs: **37.6× parallelism, zero queueing penalty, zero
failures**, and the wall clock is set by the **slowest single wave**, not by the number of
waves. Going from 2 settings to 4 costs **+32 jobs and a median +2.9 minutes**. FINDING 016
priced eight settings at 2× the budget; the real price on this venue is 2× the **jobs** and
about **3 minutes**.

**What FINDING 040 got right, and what needs correcting.** Its *ranking* survives and is
now quantified: reference depth **+0.0076** (sign-stable under drop-1) beats pool depth
**+0.0027** (sign flips under drop-1) — reference depth is the better of two small buys,
by a factor of ~2.8 in the point estimate and by a wide margin in robustness. But its
*reasoning* — *"depth in the reference is what makes the +0.0395 gain exist at all"* —
**conflates two different quantities and should not be quoted again in that form**:

* the **first four** reference poses are worth the whole **+0.039** (below 4 the feature
  measures **−0.0055**, FINDING 011 — that is the existence claim, and it is true);
* the **next four** are worth **+0.0076** (this finding — that is the marginal claim, and
  it is the one that governs the spare hour).

**Surviving release-day rule.** Submit **2 engines × 4 settings** in one concurrent wave
set. Not because depth 8 selects meaningfully better than depth 4 — it does not, by
+0.0076 with a CI touching zero — but because it is **nearly free in the currency that is
actually scarce on release day**, and because it cuts the across-draw sd by a third on a
set where you get exactly one draw.

---

## 7. The release-day question, answered directly

> **Given one spare hour, does it go to pool depth or reference depth?**

**Neither, and that is the real answer: reference depth does not need the hour.**

| | **pool depth** (20 → 40 samples) | **reference depth** (4 → 8 poses) |
|---|---|---|
| measured gain | **+0.0027** (FINDING 040, Type II, n=73) | **+0.0076** (this, n=37) |
| p | 0.870 | 0.261 |
| **drop-1** | **−0.0012 — sign flips** | **+0.0042 — sign stable, all 37** |
| oracle-to-selection gap | **widens** 0.0770 → 0.0876 | **narrows** 0.1019 → 0.0943 |
| across-draw sd | n/a (deterministic given the pool) | **0.01215 → 0.00812, −33%** |
| what it costs | GPU: ~1.6 min/ligand on an H200 → **~2.7 h for 100 ligands** | **+32 OpenProtein jobs per 37 ligands; +2.9 min median wall clock** |
| venue | Explorer, serial, blocks the submission | OpenProtein, concurrent, runs alongside everything |

**The measured wall clock, from the 80-job timeline** (`refvalue_job_timeline.json`, all
80 SUCCESS): first create to last end **38.3 min**; submit span 7.4 min; per-job duration
**4.8 / 17.0 / 37.3 min** (min/median/max); 1,443 job-minutes inside a 38.3-minute wall
clock = **37.6× parallelism with no queueing penalty at 80 jobs** — 5× the concurrency
FINDING 038 could observe at 16.

| sweep | depth | jobs (37 ligands) | wall clock, min / median / max over setting choices |
|---|---|---|---|
| 2 settings | 4 | 32 | 26.2 / **35.4** / 38.3 min |
| 3 settings | 6 | 48 | 31.1 / **38.3** / 38.3 min |
| **4 settings** | **8** | **64** | 35.4 / **38.3** / 38.3 min |
| 5 settings | 10 | 80 | **38.3** min |

Because the waves run concurrently, **the marginal wall clock of depth 8 over depth 4 is
the difference between two slowest waves — a median of 2.9 minutes, at most 12.1.** So:

> ## **Submit four settings on both Protenix checkpoints at t = 0, in one go, and spend the spare hour on nothing.**
>
> Reference depth 8 costs **+32 jobs and ~3 minutes**, buys **+0.0076** with a
> **sign-stable** drop-1, **narrows** the oracle gap, and **cuts the across-draw sd by a
> third**. Pool depth 20 → 40 costs **hours of GPU** and buys **+0.0027** whose sign
> **flips** when one ligand is removed. If the hour truly cannot be given back, it goes to
> reference depth — but the correct action is to make the hour unnecessary by submitting
> all four settings in the first wave set.

**The conditions on that recommendation, stated so they can be checked on the day:**

1. **Concurrency.** Measured to **80 jobs**. At 100 ligands × 4 settings that is 160 jobs,
   which is **unmeasured**. If queueing appears, the wall clock rises toward FINDING 038's
   conservative block-of-16 model and the third and fourth settings become a real cost.
   `refset` is the check: run it, read `below_min_depth`, and drop a setting if the clock
   is binding.
2. **Failures.** 0 of 80 here, 1 of 16 in FINDING 038 — a pooled **1 of 96 ≈ 1.0%**, with
   a wide interval. Depth 8 survives two failures per ligand; depth 4 survives none. That
   margin is worth more than the +0.0076.
3. **Distinctness.** At five settings the tightest pair reached **0.091 Å** against a
   0.05 Å tolerance. A sixth setting must be **checked**, not assumed, before it is
   counted as a sixth opinion.
4. **The bar is unchanged.** +0.0076 does **not** clear +0.0125 and its CI touches zero.
   Nothing here licenses calling reference depth a lever. It is a cheap insurance policy
   with a small positive expected value.

---

## 8. Deviations from the pre-registration

**None to the experiment.** The ligand rule, the depths, the exhaustive-subset scheme, the
draw counts, the floor construction, the acceptance clauses and the verdict map are all as
committed at `40e0f55`, and the verdict falls in a band the pre-registration named in
advance.

Two additions, both labelled where they appear and **neither entering the verdict**:

* **§5's shipped-reference arm** — the same pool scored against the *engine-replicate*
  reference set at depths 4–6. Descriptive; it is a second reference source, not a second
  test of the same one.
* **§7's job timeline** — pulled from the OpenProtein API after the fact. It is a cost
  measurement, and it is reported *after* the pre-registered verdict precisely because the
  verdict's bands were written in units of wall clock before the wall clock was known.

One methods note that is *not* a deviation: the within-ligand ρ in §4 is averaged over 50
random subsets per depth (the pre-registration did not fix the number). The identical
83.8% across all five depths was verified independently at a different seed and 200
subsets; it is the same six ligands at every depth.

---

## Artefacts

| what | where |
|---|---|
| pre-registration | `docs/PREREG_reference_depth_value.md` (`40e0f55`) |
| ligand list + the rule's slots | `data/processed/refvalue_ligands.csv` |
| analysis (frames / controls / curve / filter controls) | `scripts/structure/refvalue_analysis.py` |
| `refset` wrapper that saves the full report | `scripts/structure/refvalue_refset.py` |
| table renderer | `scripts/structure/refvalue_report.py` |
| frozen reference set — 37 ligands, **370 poses**, 136 kB | `data/processed/refvalue_reference_set.npz` |
| `refset` report, all filter counts + per-ligand depth | `data/processed/refvalue_refset_report.json` |
| filter positive controls (C4) | `data/processed/refvalue_filter_controls.json` |
| C1 / C2 / C3 | `data/processed/refvalue_controls.json`, `refvalue_xeng_shipped_repro.csv` |
| the exact curve, spread, floor, primary, leverage | `data/processed/refvalue_curve.json` |
| per-ligand table | `data/processed/refvalue_per_ligand.csv` |
| per-ligand ρ at depth 4 and 8, independent seed | `data/processed/refvalue_rho_per_ligand.json` |
| shipped-reference arm | `data/processed/refvalue_shipped_ref_arm.json` |
| 80-job timeline, all SUCCESS | `data/processed/refvalue_job_timeline.json` |
| the 370 mmCIFs (125 MB, **gitignored**) | **archived** to `onedrive:rclone-offload/cyp-structure/pool/refvalue` via `cypstruct.storage.push(move=True)`; the local `data/processed/openprotein/refvalue/` is now empty. The frozen `.npz` is the reproducible artefact — pull the mmCIFs back with `storage.pull` only if the dedupe has to be re-run |
