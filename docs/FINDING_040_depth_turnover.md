# FINDING 040 — the selected depth curve does NOT turn over. FINDING 039's pool-33 peak was 08J, and it does not survive n = 73

**Date:** 2026-09-23 · **Status:** measured · **4 Slurm jobs, 1,460 new complexes,
0 failed examples, ~7.9 GPU-hours, $0** · **Verdict: NO TURNOVER (monotone)** — the exact
union curve is **strictly increasing at every one of 39 steps**, minimum increment
**+0.00036**, argmax at depth **40**.

Pre-registered in `docs/PREREG_depth_turnover.md`, committed at **`c4a18b8`** before a
single Boltz job was submitted. No ligand, rung, threshold or acceptance clause was moved
afterwards.

Generation path and its traps: `docs/RUNBOOK_explorer_boltz.md`. Nothing was rebuilt —
this is `scripts/explorer/boltz_depth.py` with four new tags.

---

## What was asked

`FINDING_035` doubled a 14-ligand predicted-Type-I stratum from 20 poses to 40 at matched
conditioning: Δ selected **+0.0116**, Wilcoxon **p = 0.715**, **10 of 14 tied**.
`FINDING_039` §3b then recomputed that curve in **closed form**, with zero sampling error,
and found something the body of 035 had not remarked on: the **selected** score peaks at
pool **33** (0.5424) and **falls** to 0.5386 by pool 40, while the oracle stays monotone —
and the entire effect is **one ligand** (08J, +0.1721 between pool 30 and 40, against CFF's
−0.1195), with ten of fourteen unable to move at all. 039 closed by naming its successor:

> What would settle it: the same exact closed-form curve on a **second matched-depth
> stratum of ≥ 30 ligands**.

This is not a question about the record. On release day we choose a generation depth. If
depth keeps paying past 20 poses we buy more; if the selected curve genuinely turns over,
buying past ~30 is wasted wall clock and may actively cost score.

---

## 1. The stratum, and the rule — fixed before any pose existed

> **Every validation ligand in the prediction-side Type II majority: every one of the 87
> with `pred_fe_donor_median <= 2.6 Å`. By construction this is the exact complement of
> `FINDING_035`'s 14-ligand predicted-Type-I stratum — overlap is zero by definition, not
> by inspection.**

| | |
|---|---|
| validation ligands | **87** |
| pass `pred_fe_donor_median <= 2.6` | **73** ← the stratum |
| excluded as `FINDING_035`'s Type-I stratum | **14** |
| **overlap with 035's stratum** | **[] — asserted, not assumed** |
| reference depth (the `FINDING_011` precondition) | min **6**, median **7** (6:18, 7:43, 8:5, 9:2, 10:1, 11:1, 12:3); **0** below 4 |
| ligands dropped for any other reason | **0** |

Three reasons this rule and not a convenient one: it **cannot** overlap 035; it is where
the test set will sit (`FINDING_033` measured the validation set at 83% Type II and named
the mix as the one unmeasured exposure, and 035 has **no** coverage of the 83%); and it is
**prediction-side**, so it is computable on a blind set at inference time.

**Batching.** 73 × 20 samples exceeds the 7-hour walltime the `gpu` fallback partition
allows, so generation was split into four jobs by **round-robin over the alphabetically
sorted ids** (`i % 4`) — *not* contiguous blocks. That is a direct response to
`FINDING_038`, where an **ordered** csv under `--batch` packing put all three Type I picks
into the one job that failed and manufactured a chemistry signal out of one job's death.

---

## 2. The conditioning match, reported before any curve — **A0 = +0.0005**

| | Modal pool | **Explorer arm** | matched | evidence |
|---|---|---|---|---|
| `build_yaml` | — | **AST-lifted from `modal_boltz.py`**, sha1 `aefda513…` | ✔ | the same sha1 `FINDING_035` recorded |
| boltz / torch | 2.2.1 / 2.5.1 | **2.2.1 / 2.5.1** | ✔ | job `10539911` |
| checkpoint | auto-downloaded | `boltz2_conf.ckpt` | ✔ | md5 `2f0a1775bf8fc366a1a85e2019eca288`, hashed **in a job** |
| **MSA** | `cyp3a4.a3m`, 6,979 seqs | **the same file** | ✔ | md5 `6de0ee1330ea9aa1aef1efb6ec6a2e3d` on both sides, 0 NUL bytes, re-hashed far-side |
| heme bond | `[A,442,SG] ↔ [H,1,FE]` | identical | ✔ | Boltz's own parser, **73/73** |
| chains / steering | A,H,L / unsteered | identical | ✔ | no pocket, contact or template block |
| samples / recycling / sampling steps | 20 / 3 / 200 | 20 / 3 / 200 | ✔ | |
| kernels · GPU · CUDA · seed · batching | — | `--no_kernels` · V100/H200 · cu121 · 202 · 18–19/job | ✘ **K1–K6** | pre-registered §4 |

**A0, the gate:**

| | |
|---|---|
| new arm pool mean, 1,460 poses | **0.5999** |
| existing Modal pool mean, same 73 ligands | **0.5994** |
| **Δ** | **+0.0005** — inside the 0.020 bar by a factor of 40 |
| `FINDING_035`'s arm, for reference | −0.0019 |
| `FINDING_034`'s unmatched OpenProtein arm | **−0.057** |
| fraction coordinating | **1.000** new vs **0.9993** existing |
| median Fe–donor | **2.219 Å** new vs **2.219 Å** existing |

This is the closest conditioning match the repo has achieved, and the Type II character —
which is what defines this stratum — is preserved to three decimals on both statistics.

**K6, the new confound, measured rather than argued.** Four jobs, four nodes, and the
per-batch gap against the pool each batch joins is **−0.0013 / −0.0011 / +0.0019 /
+0.0026**. No job-level effect is hiding inside the arm mean.

| batch | node | ligands | wall | new mean | existing mean, same ligands |
|---|---|---|---|---|---|
| `dt2a` job `10539925` | d4079 (**H200**) | 19 | **0:30:40** | 0.5525 | 0.5538 |
| `dt2b` job `10539927` | c2205 (V100) | 18 | 2:37:02 | 0.6042 | 0.6053 |
| `dt2c` job `10539928` | c2207 (V100) | 18 | 2:39:45 | 0.6327 | 0.6308 |
| `dt2d` job `10539929` | c2206 (V100) | 18 | 2:07:49 | 0.6128 | 0.6102 |

All four **COMPLETED** (the `|| true` fix from `FINDING_035` held), all four logged
`Number of failed examples: 0`, **1,460 / 1,460** mmCIFs, **20 per ligand on all 73**.
A venue number worth recording: **the H200 is 4.6× the V100** — 19 ligands in 30 minutes
against 18 in 2 h 38 m.

---

## 3. Controls, with counts, including the zeros

| id | control | result |
|---|---|---|
| **N2** | the shipped board reproduces on all 87 | selected **0.6164**, oracle **0.6975**, random **0.5769**, gain **+0.0395**, ρ **−0.2582**, correct sign **75.86%** — `FINDING_011`/`033` to four decimals |
| **N2b** | `argmin(xeng)` **is** `cypstruct.xengine.select()` | max abs diff **0.0**, equal to 1e-9 on all 87. Asserted |
| **N1** | numbering (`FINDING_021`) on 1,460 new poses | offset **0 on 1,460 of 1,460**, min residue-name identity **0.9956**, **0** below the 0.95 bar, `mapped` True on **1,460 of 1,460** |
| **C7** | exact-zero LDDT-PLI rows | **0** (035 had 2, both PEG ejections; this stratum has none) |
| **N4** | distinct poses, not jobs | 1,460 rows → **1,460 distinct md5** → **1,460 distinct at 0.05 Å in the heme frame** = **20.0 per ligand**, min 20 |
| **N5** | filters | 87 → **73** pass, **14** excluded. Scoring skips: no_crystal 0, no_ligand 0, load_failed 0, **mapping_failed 0**, no_heme_frame 0, low_renumber_identity 0, no_files 0 |
| **REF** | reference depth | min **6**, all ≥ the 4 `FINDING_011` requires |
| | existing pool diversity in the stratum | **20 unique `xeng` on every one of 73**, median within-ligand LDDT sd 0.0577, zero ligands with sd 0 |

### X1 — the closed forms reproduce `FINDING_039` §3b exactly, before they are used here

New estimators on a new stratum are not believed until they replay a published exact curve.
Run on 035's 14-ligand Type-I stratum, **maximum absolute deviation 0.00000** across all
eight published rungs:

| pool | 20 | 24 | 28 | 30 | **33** | 36 | 38 | **40** |
|---|---|---|---|---|---|---|---|---|
| `FINDING_039` | 0.5270 | 0.5349 | 0.5400 | 0.5415 | **0.5424** | 0.5418 | 0.5406 | **0.5386** |
| **here** | 0.5270 | 0.5349 | 0.5400 | 0.5415 | **0.5424** | 0.5418 | 0.5406 | **0.5386** |

and the union curve reproduces 039's argmax at depth **34** (peak 0.5411 vs its 0.5410),
and 035's oracle endpoints **0.6499 / 0.6595** exactly. A brute-force Monte-Carlo
cross-check agrees to 0.0004 at 3,000 draws.

### The noise floor, recomputed on **this** population

`FINDING_039` established that a floor is a number *about a population* and that a
4,000-draw p95 carries ±0.0009. Nothing is quoted from another finding; both floors below
are **2,000,000 draws**.

| population | **p95** | SE of the p95 | **p99** | sd of the gain |
|---|---|---|---|---|
| **prediction-side Type II, n = 73, depth 20** | **+0.01383** | **±0.000013** | +0.01935 | 0.00858 |
| **prediction-side Type II, n = 73, depth 40** | **+0.01367** | **±0.000012** | +0.01909 | 0.00852 |

For orientation: 039's three authoritative floors are +0.04306 / +0.04490 / +0.04400, all
at **n = 14**. Five times the ligands buys a floor three times lower — which is exactly why
035 could not settle this and this can.

**The shipped selector clears its floor on this stratum, comfortably, at both depths:**
gain **+0.03415** at depth 20 and **+0.03659** at depth 40, against +0.01383 / +0.01367.
Within-ligand ρ **−0.2625 → −0.2698**, correct sign **76.71% → 75.34%**.

---

## 4. The closed-form curve — the primary endpoint

Exact, not subsampled. `argmin` over a uniformly random subset has a combinatorial
distribution, so every rung below is an expectation with **zero sampling error**.

**Union curve** — `d` poses drawn uniformly from each ligand's 40-pose union (20 Modal + 20
Explorer), n = 73:

| depth | **oracle** | **selected** | random |
|---|---|---|---|
| 1 | 0.59962 | 0.59962 | 0.59962 |
| 2 | 0.63155 | 0.61263 | 0.59962 |
| **5** | 0.66186 | **0.61965** | 0.59962 |
| **10** | 0.68256 | **0.62353** | 0.59962 |
| **15** | 0.69443 | **0.62621** | 0.59962 |
| **20** | 0.70279 | **0.62852** | 0.59962 |
| **25** | 0.70933 | **0.63062** | 0.59962 |
| **28** | 0.71270 | **0.63180** | 0.59962 |
| **30** | 0.71478 | **0.63256** | 0.59962 |
| **33** | 0.71770 | **0.63367** | 0.59962 |
| **35** | 0.71953 | **0.63440** | 0.59962 |
| **38** | 0.72215 | **0.63549** | 0.59962 |
| **40** | **0.72382** | **0.63622** | 0.59962 |

> **The selected column is STRICTLY INCREASING at all 39 steps. Minimum increment
> +0.00036. argmax = 40. The oracle is strictly increasing too.**

**Augment curve** — the *purchase*: all 20 existing poses kept, `k` new added. Also
monotone, also argmax at pool 40.

| pool | 20 | 24 | 28 | 30 | **33** | 35 | 38 | **40** |
|---|---|---|---|---|---|---|---|---|
| oracle | 0.70663 | 0.71145 | 0.71519 | 0.71682 | 0.71910 | 0.72051 | 0.72253 | **0.72382** |
| **selected** | 0.63353 | 0.63370 | 0.63405 | 0.63430 | **0.63477** | 0.63514 | 0.63576 | **0.63622** |

| the 20 → 40 purchase | value | 95% CI | |
|---|---|---|---|
| **Δ selected** | **+0.0027** | [−0.0130, +0.0182] | **19 better / 16 worse / 38 unchanged** |
| **Δ oracle** | **+0.0172** | **[+0.0102, +0.0251]** | excludes zero |
| Wilcoxon signed-rank | | | **p = 0.870** |
| `FINDING_035`, same purchase on Type I | +0.0116 | [−0.0687, +0.1034] | p = 0.715, 10 of 14 tied |

**Rates per doubling, on the unbiased union curve:**

| octave | oracle | **selected** |
|---|---|---|
| 2 → 5 | +0.0303 | +0.0070 |
| 5 → 10 | +0.0207 | +0.0039 |
| 10 → 20 | +0.0202 | **+0.0050** |
| **20 → 40** | +0.0210 | **+0.0077** |

**The last-rung caveat, stated rather than glossed.** `FINDING_035`'s methods correction
holds here too: depth 40 is the whole pool, so *both* the oracle and the selected value at
that rung are deterministic — `max` and `argmin` over everything, with zero variance —
while every earlier rung is an expectation over random subsets. The closed form removes
Monte-Carlo error; it does **not** remove that structural fact. So the 20 → 40 slope is
still a slope measured *into* a fixed endpoint and should not be extrapolated past 40.
Nothing in the recommendation below extrapolates past 40.

---

## 5. Per-ligand decomposition — so a one-ligand artefact cannot hide again

### How many ligands can move at all

| | this stratum | `FINDING_035` / `039` |
|---|---|---|
| ligands with **≥ 1 new pose beating their incumbent on `xeng`** | **36 of 73 — 49.3%** | 4 of 14 — 28.6% |
| ligands that **cannot** move at any depth | **37** | 10 |
| ligands whose selected score actually changed at k = 20 | **35** (19 up, 16 down) | 4 |

**Half the stratum is unreachable by depth**, and that is the hard ceiling on any depth
purchase: 37 of 73 ligands have no new pose the selector would ever take, so the other
half carries the whole effect. Better than 035's 10-of-14, and still a bound worth saying
plainly. (One ligand is movable but tied — its new winner has an identical LDDT-PLI.)

### The 30 → 40 comparison — the one 039 found inverted

| | this stratum, n = 73 | `FINDING_039`, n = 14 |
|---|---|---|
| **mean sel(30) − sel(40)** | **−0.00366** (40 is better) | **+0.0029** (30 was better) |
| ligands with a positive contribution | **31 of 73** | 4 of 14 non-zero |
| **largest single-ligand contribution** | **+0.0604** (X71) | **+0.1721** (08J) |

No single ligand can invert it: the largest positive contribution is +0.060 against a total
of −0.267. The adversarial-pose mechanism 039 identified is **still present** — X71's one
new `xeng`-winning pose costs it 0.252 LDDT-PLI, D7J's costs 0.150 — it is simply
**outvoted** at n = 73 by A1A4V (+0.281), A1ASQ (+0.218) and X5Y (+0.154).

### Leverage — the curve is robust, the *point estimate* of the purchase is not

| | n | Δ selected, 20 → 40 (purchase) | union sel(20) → sel(40) |
|---|---|---|---|
| all 73 | 73 | **+0.0027** | 0.62852 → 0.63622 (**+0.0077**) |
| drop A1A4V | 72 | **−0.0012** | +0.0059 |
| drop top 2 | 71 | +0.0024 | +0.0076 |
| drop top 3 | 70 | +0.0056 | +0.0081 |
| drop top 5 | 68 | +0.0003 | +0.0059 |

**This is the honest split.** The *shape* of the curve is robust — the union delta stays
between +0.0059 and +0.0081 under every deletion — while the *purchase* point estimate
swings sign when one ligand is removed. Quote **+0.0027 with `drop-1 = −0.0012`
attached**, exactly as 035 must be quoted with its `drop-both = +0.0020`.

### The combined 87 — not pre-registered, reported as a descriptive aggregation

Stitching this stratum to 035's gives every validation ligand at matched depth 40:

| depth | 5 | 10 | 20 | 30 | **33** | 38 | **40** |
|---|---|---|---|---|---|---|---|
| oracle | 0.64405 | 0.66743 | 0.69062 | 0.70398 | 0.70711 | 0.71176 | **0.71347** |
| **selected** | 0.59919 | 0.60480 | 0.61256 | 0.61769 | **0.61877** | 0.62012 | **0.62051** |

**Monotone, argmax 40.** 039's pool-33 peak does not survive at n = 87 either.

---

## 6. The verdict against the pre-registered clauses

| clause | required | measured | fires? |
|---|---|---|---|
| **A0** (gate) | pool-mean gap ≤ 0.020 | **+0.0005** | **passed** |
| yield | ≥ 30 ligands with 20 distinct new poses | **73 of 73** | **passed** |
| **T1 — existence** | `d* < 40` | **`d* = 40`** | **NO** |
| **T2 — magnitude** | Δ_peak ≥ +0.0050 | 0.0 | no |
| **T3 — not one ligand** | sign survives leave-one-ligand-out | n/a | no |
| **T4 — inference** | CI excludes 0 and p < 0.05 | n/a | no |
| | | | |
| **NO TURNOVER (monotone)** | ¬T1 | | **YES** |

> ## **NO TURNOVER.**
>
> **The selected depth curve does not turn over on this stratum, at any depth from 1 to
> 40. It rises at every single step. `FINDING_039`'s pool-33 peak was 08J, it is a
> property of n = 14, and it does not generalise — not to 73 Type II ligands, and not to
> the combined 87.**

### What this settles, and what it does not

1. **The turning point does not exist at this depth.** It is not "too small to detect" —
   the curve is **strictly** increasing, with zero sampling error, at every one of 39
   steps. There is no depth in 1–40 at which generating more poses costs selected score
   on the Type II majority.
2. **The mechanism 039 described is real and is a variance term, not a trend.** Adversarial
   poses that win on `xeng` and lose on LDDT-PLI exist here too, and cost up to 0.252 on
   a single ligand. At n = 14 one of them steers the mean; at n = 73 they cancel. This is
   `FINDING_035`'s own thesis — the near-tie is a coin flip with an enormous payoff — with
   the sample size that shows it is a coin flip.
3. **Depth still does not convert, and the reason is unchanged.** The oracle climbs
   +0.0172 [+0.0102, +0.0251] over the doubling and selection takes **+0.0027** of it,
   p = 0.870, with 38 of 73 tied. The **oracle-to-selection gap at depth 40 is 0.0876** on
   this stratum, larger than the 0.0770 it was at depth 20: **buying depth widens the gap
   the selector fails to close.** `FINDING_001` is unmoved.
4. **`FINDING_035`'s "the conversion rate decays with depth" is Type-I-specific.** 035
   measured +0.0227 at 10→20 falling to +0.0095 at 20→40. Type II does the opposite:
   **+0.0050 at 10→20 and +0.0077 at 20→40**, rising. Both are small; the *decay* is not a
   law and should not be quoted as one.
5. **Beyond 40 is still unmeasured**, and the last rung's slope must not be extrapolated
   (§4). This says nothing about 40 → 80.
6. **Half the stratum is beyond reach of depth at all** — 37 of 73 have no new pose the
   selector would take. That bounds every future depth purchase and is measurable before
   buying one.

---

## 7. The decision — how deep to generate on release day

**Twenty diffusion samples per ligand, in one Boltz job. Buy the second twenty only if
wall clock is free; never delay a submission for it.**

The reasoning, and both halves matter:

* **Depth is never harmful.** The curve is monotone on the Type II majority (this
  finding), on the combined 87, and the Type I stratum's apparent turnover was one ligand.
  There is no "buy too much" failure mode to avoid. This retires a live worry.
* **Depth is also barely worth anything.** A full 20 → 40 doubling at matched conditioning
  buys **+0.0077** selected on the union curve and **+0.0027** as an actual purchase
  (CI [−0.0130, +0.0182], p = 0.870, drop-1 = −0.0012) — on the 83% of the test set that
  will be Type II. On Type I it is +0.0116 (`FINDING_035`). Neither clears `FINDING_004`'s
  +0.0125-per-doubling benchmark, and **no octave of this curve does**: +0.0070, +0.0039,
  +0.0050, +0.0077.
* **The cost is real.** 73 ligands × 20 samples was 7.9 GPU-hours. Doubling it doubles
  that, and on the V100 fallback a single 18-ligand batch is 2 h 40 m.
* **Spend the wall clock on reference depth instead**, which `FINDING_038` priced at
  25.7 minutes for 9 ligands and which `FINDING_011` makes a hard precondition of the
  selector. Depth in the *pool* is worth +0.0027; depth in the *reference* is what makes
  the +0.0395 gain exist at all.

A one-line recommendation for the playbook:

> **`--samples 20`, one job per ~18 ligands, and stop. The selected-score curve is
> monotone to depth 40 — more never hurts, and a full doubling buys +0.0027 on the Type II
> majority (p = 0.87) against +0.0125 needed to be worth the wall clock.**

---

## 8. Deviations from the pre-registration

**One, and it is to a check rather than to the experiment.**
`scripts/explorer/verify_inputs.sbatch` had the list `("smoke", "stratum")` hardcoded into
its YAML-parser check, so for any later tag it would have verified **zero files and
reported success** — the repo's standing trap, sitting inside the tool built to avoid it.
It now enumerates every staged subdirectory and prints `NOTHING VERIFIED` for an empty one.
Changed *before* `verify` ran; the run then exercised all six subdirectories
(`dt2a` 19, `dt2b` 18, `dt2c` 18, `dt2d` 18, `smoke` 1, `stratum` 14 — 0 bad).

Nothing in §2, §6 or §7 of the pre-registration was edited. The two additions marked
**not pre-registered** in this document — the combined-87 curve (§5) and the leverage
table (§5) — are descriptive, are labelled where they appear, and neither enters the
verdict.

---

## Artefacts

| what | where |
|---|---|
| pre-registration | `docs/PREREG_depth_turnover.md` (`c4a18b8`) |
| generation path + traps | `docs/RUNBOOK_explorer_boltz.md` |
| runner (unchanged) | `scripts/explorer/boltz_depth.py` |
| analysis | `scripts/structure/depth_turnover_analysis.py` (`x1`, `controls`, `distinct`, `curves`, `decomp`) |
| ligand list and batches | `data/processed/depth_turnover_stratum.csv`, `depth_turnover_dt2[a-d].csv`, `depth_turnover_batches.json` |
| the 1,460 new poses, scored | `data/processed/depth_turnover_poses.csv` (per batch: `matched_depth_dt2[a-d]_poses.csv`) |
| controls incl. A0, N1, N4, N5, C7, K6 | `data/processed/depth_turnover_controls.json` |
| X1, the 039 replay | `data/processed/depth_turnover_x1.json` |
| distinct-pose dedupe | `data/processed/depth_turnover_distinct.json` |
| exact curves, both, + floors + acceptance | `data/processed/depth_turnover_curves.json` |
| per-ligand decomposition | `data/processed/depth_turnover_per_ligand.csv`, `depth_turnover_decomp.csv`, `depth_turnover_decomp.json` |
| job ledger | `data/processed/matched_depth_jobs.json` (jobs `10539911`, `10539925`, `10539927`, `10539928`, `10539929`) |
| the 1,460 mmCIFs | `C:\cyp_struct\matched_depth\poses\dt2*_flat\` and `/scratch/shenoy.am/cyp-depth/out/dt2*_s202/` — never `D:` |
