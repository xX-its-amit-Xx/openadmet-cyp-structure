# FINDING 033 — the shipped selector does not care about binding mode; the POOL does

**Date:** 2026-09-22 · **Status:** measured, zero new inference, zero new downloads, CPU-only ·
**Verdict:** **ROBUST** — on the family-wide set, where the power is. The CYP3A4-only test
is **underpowered** and is reported as such rather than dressed up in either direction.

**Script:** `scripts/structure/binding_mode_robustness.py` (`labels` / `cyp3a4` / `p450`).
Thresholds are the repo's existing constants, unchanged and fixed before any number was
computed: `COORD_MAX = 2.6 Å`, `OVERHEAD_MAX = 6.5 Å`.

---

## What was asked

`docs/worldmodel/CYP3A4_BIOLOGY_MAP.md` §4 and §5 measured that the 87-ligand validation
set is **83% Type II** — the ligand coordinates the heme iron directly — and that Reactome
places CYP3A4 in exactly one of six CYP-by-substrate-type pathways (Xenobiotics) and in
none of the sterol, fatty-acid, eicosanoid or vitamin ones, although UniProt documents all
of those substrates. Both frames we sample from are biased toward heme-coordinating
chemistry.

Everything we ship rests on that set. `cypstruct.xengine.select()` has been validated at
+0.0395 / +0.0357 / +0.0383 across three pools. It scores a pose by its mean symmetric
Chamfer distance **in the heme frame** to independent-engine poses of the same ligand — and
that frame is built from the heme. **A Type I ligand is not anchored to the landmark that
defines the frame.** So there is a specific, mechanistic reason the gain could be
mode-dependent, and if it is, and the blind test set is Type I-rich, every selector number
in this repo was measured on the wrong population.

> **The question.** Does the shipped selector's gain over random survive on ligands that do
> not touch the iron — and if not, is that the selector failing or the pool being worse?

The test set is not released and the interim deadline is on 2026-09-24, so this is the last
cheap chance to know.

---

## What was done

| | |
|---|---|
| CYP3A4 pool | `poses_scored_val87b.csv`, unsteered arm — 87 ligands × 20 Boltz-2 poses = **1,740** |
| CYP3A4 feature | `xeng_val87b.csv`, the shipped column, Protenix ×2 deduped reference |
| family pool | `p450_universe/` — protenix_v2, deduped, **342 pairs / 1,239 poses / 66 proteins (65 non-CYP3A4)** |
| family reference | esmfold2 of the same pair, one per replicate, deduped — FINDING 012's setup exactly |
| labels | crystal Fe-to-closest-ligand-atom, re-derived; prediction-side label computed separately |
| random baseline | **exact per-ligand expectation** (mean of that ligand's pool) *and* a 300-draw empirical average; they agree to 0.0014 |
| null | a random feature, **4,000 draws inside each stratum**, because the FINDING 007 floor is n-dependent |
| CIs | 10,000-resample percentile bootstrap over ligands/pairs |

### Two controls first

**C1 — the shipped column reproduces exactly.** Random 0.5783, oracle 0.6975, selected
**0.6164**, gain **+0.0381**, within-ligand ρ **−0.2582**, correct sign on **75.9%**.
FINDING 011's full-depth deduped row is 0.6164 / +0.0380 / −0.258 / 75.9%. Nothing has
drifted. `per_unit_table()` is asserted against `cypstruct.xengine.select()` to 1e-9 on
every run, so the stratified arithmetic and the shipped function are the same selector.

**C2 — the four exact-zero LDDT-PLI rows are ejections, not FINDING 021.** A numbering
failure is a property of the **pair**, not of the pose: if the offset were wrong, all 20
poses of that ligand against that crystal would read 0.000.

| ligand | pdb | zero poses | BiSyRMSD of the zeros | median LDDT-PLI of the ligand's OTHER poses |
|---|---|---|---|---|
| PG0 | 9BV6 | 1 / 20 | 23.1 Å | 0.441 |
| PG4 | 9YK4 | 3 / 20 | 25.9, 26.9, 26.1 Å | 0.577 |

`mapped` is `True` on every row. Both ligands are PEG fragments, both are Type I, and the
zeros carry 23–27 Å ligand RMSD. These are genuine ejections. Confirmed, not assumed.

---

## 1. The labels — one is wrong, and the classes separate perfectly

Re-derived from the crystal Fe-to-closest-ligand-atom distance in
`validation_ligands.csv`, against the stored `cls` string.

| | stored | **re-derived** |
|---|---|---|
| Type II (≤ 2.6 Å) | 72 | **73** |
| Type I (2.6–6.5 Å) | 15 | **14** |

**`PK9` (4D6Z) is mislabelled.** It is stored `type_I_active_site` and its crystal Fe
distance is **2.379 Å** — inside the coordination window by 0.22 Å, and 1.5 Å below the
nearest genuine Type I. The biology map's "72 of 87 / 15" should read **73 / 14**. This is
the only disagreement in 87.

**The separation is complete — AUC = 1.00, and there is nobody in the gap.**

| | n | min | p5 | median | p95 | max |
|---|---|---|---|---|---|---|
| Type II | 73 | 1.935 | 2.026 | **2.200** | 2.399 | **2.581** |
| Type I | 14 | **2.824** | 2.949 | 4.109 | 5.217 | 5.303 |

The closest Type I sits **0.244 Å** beyond the furthest Type II (Mann–Whitney
p = 3.7e−09). The 2.6 Å cut could be moved anywhere in [2.59, 2.82] without relabelling a
single ligand, so no result below turns on where it was put. That is the check the task
asked for, and it passes.

**The label is computable from the PREDICTION on 84 of 87 (96.6%).** Rule fixed a priori:
median `fe_donor_dist` over the ligand's own 20 poses, same 2.6 Å cut, nothing from the
crystal. Median fraction of poses coordinating is **1.00** for crystal Type II and **0.00**
for crystal Type I — the co-folder already knows.

The three it gets wrong are worth naming, because they are the ones that would be
mislabelled on a blind test ligand:

| ligand | pdb | crystal Fe | predicted median | poses coordinating | verdict |
|---|---|---|---|---|---|
| D0R | 3TJS | 2.107 Å (II) | 5.430 Å | 25% | model fails to coordinate a genuine Type II |
| **QDY** | 6UNJ | 2.824 Å (I) | 2.231 Å | **100%** | model **forces** coordination on a non-coordinator |
| MWV | 6OO9 | 3.763 Å (I) | 7.695 Å | 0% | right mode, ejected past the active site |

QDY is the only ligand in the set whose label "looks wrong" from the prediction side, and
the crystal is unambiguous: 2.824 Å is not a bond. It is the cleanest single example of
the bias this finding was commissioned to look for — the model's prior is Type II, and it
applies it to a ligand that does not obey it.

**All numbers below use the re-derived 73 / 14 split.**

---

## 2. CYP3A4, stratified — the point estimates are identical and the CI is enormous

| | Type II | Type I | all |
|---|---|---|---|
| n ligands | **73** | **14** | 87 |
| pool mean (= random, exact) | 0.5985 | **0.4644** | 0.5769 |
| oracle | 0.7068 | **0.6488** | 0.6975 |
| selected | 0.6376 | **0.5055** | 0.6164 |
| **gain** | **+0.0392** | **+0.0411** | +0.0395 |
| **gain 95% CI** | **[+0.0228, +0.0569]** | **[−0.0073, +0.1012]** | [+0.0232, +0.0571] |
| headroom (oracle − random) | 0.1084 | **0.1843** | 0.1206 |
| fraction of headroom captured | 36.2% | 22.3% | 32.8% |
| within-ligand ρ | −0.2632 | −0.2319 | −0.2582 |
| ρ correct sign | 76.7% | 71.4% | 75.9% |
| ligands with positive gain | 71.2% (52/73) | **50.0% (7/14)** | 67.8% |
| its own null, p95 / p99 | +0.0144 / +0.0198 | **+0.0433 / +0.0622** | +0.0139 / +0.0193 |
| empirical p vs its own null | **0.0000** | **0.0595** | 0.0000 |

**Read the last two rows before the gain row.** At n = 14 a *random* feature scores
**+0.0433** at the 95th percentile. Type I's +0.0411 does not clear its own noise floor
(p = 0.0595). So on CYP3A4 Type I ligands **we cannot demonstrate that the selector works
at all** — not that it fails, that the measurement has no resolution. This is FINDING 007's
arithmetic, applied to a stratum instead of a feature.

The leverage confirms it. Jackknifing one ligand out:

| | gain | drop the best 1 | drop the best 2 |
|---|---|---|---|
| Type II (n=73) | +0.0392 | +0.0350 | +0.0324 |
| **Type I (n=14)** | **+0.0411** | **+0.0196** | **+0.0048** |

Two ligands — YNV (+0.320) and 08J (+0.198) — carry the entire Type I mean. Seven of
fourteen have a negative gain.

### The difference, tested properly

| | |
|---|---|
| difference (Type II − Type I) | **−0.0019** |
| 95% CI | **[−0.0634, +0.0497]**, width **0.1131** |
| Welch p / Mann–Whitney p / permutation p | 0.95 / 0.56 / 0.94 |
| pooled per-ligand SD | 0.0807 |
| **minimum detectable difference, 80% power, α=0.05** | **0.0660** |
| observed / MDD | **0.03** |

**This test cannot see anything smaller than 0.066 LDDT-PLI**, which is 1.7× the whole
effect being protected. The observed difference is 3% of that. The honest statement is not
"the strata are the same" — it is that a difference as large as **±0.05** in either
direction is entirely consistent with these data, and the CI width is 0.113. **At n = 14,
CYP3A4 alone cannot answer this question.** The headroom-matched version (§4) only reaches
two usable bins and is equally uninformative: +0.0234, CI [−0.0331, +0.0918].

### What IS resolved at n = 87, and matters more

The absolute numbers, not the gain. Type I ligands have a **pool that is 0.134 LDDT-PLI
worse** (0.4644 vs 0.5985) and an **oracle 0.058 worse** (0.6488 vs 0.7068). After
selection they land at **0.5055 against 0.6376** — a 0.132 gap that has nothing to do with
selection and everything to do with generation. That gap is 3.5× the selector's entire
contribution, and it is measured on 14 ligands against 73 with a CI on the pool means that
does not come close to overlapping.

---

## 3. The family-wide set — where the n is, and where the answer is

342 pairs, 1,239 poses, 66 proteins (65 of them not CYP3A4), reproducing FINDING 012's
configuration. Selected 0.7436 against an empirical random of 0.7022 — **+0.0414**, the
published number to four decimals. The `all` row below uses the exact expectation
(+0.0402); both are reported.

Every pair in this set has its ligand inside the active site — `closest_fe` runs
1.62–6.26 Å — so the family splits cleanly into Type II and Type I with **no peripheral
stratum**, and the two strata are nearly balanced.

| | Type II | Type I | all |
|---|---|---|---|
| n pairs | **173** | **169** | 342 |
| pool mean (random, exact) | 0.6617 | **0.7461** | 0.7034 |
| oracle | 0.7475 | 0.7874 | 0.7673 |
| selected | 0.7176 | 0.7702 | 0.7436 |
| **gain** | **+0.0559** | **+0.0241** | +0.0402 |
| **gain 95% CI** | **[+0.0440, +0.0685]** | **[+0.0159, +0.0329]** | [+0.0326, +0.0478] |
| headroom | 0.0858 | 0.0413 | 0.0638 |
| fraction of headroom captured | 65.1% | 58.3% | 62.9% |
| within-pair ρ | −0.2305 | −0.2560 | −0.2431 |
| ρ correct sign | 67.3% | 65.9% | 66.6% |
| pairs with positive gain | 70.5% | 62.7% | 66.7% |
| its own null p95 / p99 | +0.0197 / +0.0261 | +0.0120 / +0.0159 | +0.0121 / +0.0163 |
| empirical p | 0.0000 | 0.0000 | 0.0000 |
| drop the best 1 / best 2 pairs | +0.0548 / +0.0538 | +0.0227 / +0.0215 | +0.0396 / +0.0390 |

**On the family set the selector clears its own null on BOTH strata, by a factor of two or
more, with no leverage problem.** That is the result CYP3A4 alone could not produce: the
shipped selector demonstrably works on ligands that do not touch the iron, on 169 pairs
across 65 unseen proteins.

### The raw difference looks large — and it is entirely pool composition

| | |
|---|---|
| raw difference (Type II − Type I) | **+0.0318** |
| 95% CI | [+0.0172, +0.0469] |
| permutation p / Welch p / Mann–Whitney p | 0.0000 / <1e−4 / 0.0004 |
| MDD at this n | 0.0213 — **observed is 1.5× MDD, so this is a real difference in the raw statistic** |

At face value that is MODE-DEPENDENT. It is not, and the repo already knows why: FINDING
012's eightfold family-wide gain was **pool composition, not a better selector**, and the
CLAUDE.md rule says to report gain-relative-to-random *within* stratum and to check whether
a between-stratum difference is just different oracles. Here the two strata have
**twice the headroom on one side** (0.0858 vs 0.0413). Binning both strata into five
headroom quintiles and differencing inside each:

| headroom quintile | mean headroom | n Type II | n Type I | gain II | gain I | difference |
|---|---|---|---|---|---|---|
| 1 | 0.002 | 19 | 50 | −0.0009 | −0.0001 | −0.0008 |
| 2 | 0.009 | 30 | 38 | −0.0000 | −0.0008 | +0.0007 |
| 3 | 0.025 | 33 | 35 | −0.0018 | +0.0075 | −0.0094 |
| 4 | 0.088 | 37 | 31 | +0.0594 | +0.0482 | +0.0112 |
| 5 | 0.194 | 54 | 15 | +0.1398 | +0.1563 | −0.0164 |
| **weighted** | | **173** | **169** | | | **−0.0030** |

**Headroom-matched difference −0.0030, 95% CI [−0.0132, +0.0075].** The +0.0318 raw
difference vanishes, and this CI is **tight** — at n = 342 we can exclude a mode effect
larger than about ±0.013 once pool quality is held fixed. The fraction of headroom each
stratum converts is 65.1% vs 58.3%, difference CI **[−0.084, +0.234]**: indistinguishable.

The quintile table also re-states FINDING 012 in a new place: **the gain lives entirely in
the top two headroom bins** (+0.06 and +0.14) and is exactly zero in the bottom three, in
*both* strata. The selector is a catastrophe detector, and what it detects has nothing to
do with the iron.

### The two alternative explanations, both measured, both dead

| alternative | measurement | runs the right way? |
|---|---|---|
| Type I pairs got shallower reference sets | reference depth **3.55 (I) vs 2.91 (II)** | **no — Type I has MORE references** |
| Type I pairs got shallower pools | poses/pair 3.51 (I) vs 3.73 (II) | marginal, and the matched test already controls the consequence |
| the mode effect is really a target effect | within the 5 proteins holding ≥3 of each mode, the Type II − Type I difference is **−0.020, +0.028, −0.007, +0.030, +0.070** | inconsistent in sign — no per-protein pattern |

---

## 4. The frame explanation — tested directly, and it does not hold

The mechanism proposed in the brief: the heme frame is built from a landmark a Type I
ligand is not anchored to, so the gain should **fall as the ligand sits further from the
iron** — within Type II as well as between strata.

| test | CYP3A4 (n=87) | family (n=342) |
|---|---|---|
| gain vs crystal Fe distance, all | ρ = **−0.008**, p = 0.94 | ρ = **−0.214**, p = 6.4e−05 |
| gain vs crystal Fe distance, **Type II only** | ρ = **+0.049**, p = 0.68 | ρ = −0.127, p = 0.095 |
| gain vs **predicted** Fe distance, all | ρ = +0.058, p = 0.59 | — |
| within-ligand ρ vs Fe distance | ρ = +0.013, p = 0.90 | — |
| **the confound:** gain vs pool headroom | ρ = **+0.440**, p < 1e−4 | ρ = **+0.650**, p = 2e−42 |
| **partial:** gain vs Fe distance, headroom controlled, all | −0.090, CI [−0.305, +0.140] | **−0.033, CI [−0.132, +0.075]** |
| partial, **Type II only** | +0.033, CI [−0.199, +0.255] | −0.138, CI [−0.273, **+0.011**] |

The one number that looks like the hypothesis — family-wide ρ = −0.214 at p = 6e−05 —
**collapses to −0.033 when headroom is partialled out**, with a CI that straddles zero. The
within-Type-II partial of −0.138 is the single surviving hint in the predicted direction,
and its CI touches zero (approximate permutation p = 0.079); at n = 173 that is a
non-result, and it is contradicted by the decisive test, the headroom-matched between-strata
difference of −0.003 ± 0.010.

**The frame explanation fails, and there is a reason it should have.** The Chamfer feature
is z-scored *within* the ligand, so an absolute scale shift cannot affect the ranking — and
Type I ligands do not have a degenerate frame, they have a **wider** one:

| | mean Chamfer | within-unit SD of Chamfer |
|---|---|---|
| CYP3A4 Type II | 1.084 Å | 0.153 Å |
| CYP3A4 **Type I** | **1.937 Å** | **1.000 Å** |
| family Type II | 1.839 Å | 1.767 Å |
| family Type I | 1.284 Å | 0.696 Å |

On CYP3A4 the reference engines disagree **6.5× more** about where a Type I ligand goes,
within a single ligand. That is more dynamic range for the feature to rank with, not less.
The frame is not the problem. The heme is a landmark for the *coordinate system*, not an
anchor the ligand has to touch, and z-scoring makes the distinction moot.

---

## Verdict — **ROBUST**

Stated with its scope, because the two halves of the evidence are not equally strong:

1. **Family-wide (n = 342 pairs, 169 of them Type I, 65 unseen proteins): robust, and
   properly powered.** Both strata clear their own nulls at p = 0.0000. The raw +0.0318
   difference is pool composition; matched on headroom it is **−0.0030, CI [−0.0132,
   +0.0075]** — a tight null that excludes any mode effect above ~0.013.
2. **CYP3A4 alone (n = 14 Type I): underpowered, and reported as such.** Difference
   −0.0019, **CI [−0.0634, +0.0497]**, MDD 0.0660. Type I's own +0.0411 does not clear its
   own n=14 noise floor of +0.0433 (p = 0.0595) and collapses to +0.0048 if two ligands are
   dropped. **Nothing here is evidence either way**; the wide CI is the deliverable.

The overall verdict is ROBUST because (1) is a direct, adequately powered test of exactly
the question and (2) is silent — not because (2) agreed.

**What replaces the worry.** The exposure to a Type I-rich test set is real, but it is not
in the selector. It is in the pool:

| | CYP3A4 Type II | CYP3A4 Type I | cost |
|---|---|---|---|
| random | 0.5985 | 0.4644 | −0.134 |
| **oracle** | **0.7068** | **0.6488** | **−0.058** |
| selected (shipped) | 0.6376 | 0.5055 | **−0.132** |

A fully Type I test set would score about **0.51 instead of 0.64**, and **selection cannot
recover it**: the oracle itself is 0.058 lower, so more than half the gap is unreachable
with this pool no matter how well we choose. That is a *generation* problem, and it is the
thing to act on.

---

## What to do at submission time if the test set turns out to be Type I-rich

The test ligands are not released, so this has to be a decision rule, not a plan.

1. **Label the test ligands on the prediction side, before anything else.** The rule is
   measured at 96.6% here and needs nothing from a crystal: fraction of a ligand's own
   pool poses with `fe_donor_dist ≤ 2.6 Å`. Median 1.00 for Type II, 0.00 for Type I. One
   line, no new inference. Then report the composition next to the submission.
2. **Do not change the selector.** `-z(xeng)` is unchanged on Type I by the only adequately
   powered test available (−0.003 ± 0.010), and every alternative in this repo is worse on
   both strata. Changing a selector on the strength of an n=14 stratum is exactly the
   post-hoc pick FINDING 011 retracted `-zm - zx` for.
3. **Spend the remaining budget on POOL DEPTH for the Type I ligands, not on scoring.**
   This is the actionable consequence. Type I oracle is 0.058 below Type II, and FINDING
   004 measured the oracle still climbing with sampling at +0.0125 of selected score per
   doubling. Doubling the pool on the predicted-Type-I subset is a small number of jobs and
   attacks the part of the gap that selection cannot.
4. **Expect, and pre-announce, a lower absolute score if the set is Type I-rich.** ~0.51
   against ~0.64 on CYP3A4-like chemistry. Do not read that as the selector failing — the
   fraction of available headroom captured is statistically indistinguishable between the
   strata on both sets (CYP3A4 36% vs 22%, difference CI [−0.17, +0.44]; family 65% vs
   58%, difference CI [−0.08, +0.23]).
5. **Watch QDY's failure mode specifically.** On 1 of 14 Type I ligands the co-folder forced
   coordination on a ligand that does not coordinate — all 20 poses at 2.23 Å against a
   crystal 2.82 Å. If the test set is Type I-rich this error becomes systematic, and it is
   *detectable without the answer*: a ligand whose poses all coordinate but whose scaffold
   has no accessible sp²/sp³ nitrogen donor (`cypstruct.chem.coordinating_atoms`, 98.8%
   recall) is a flagged case. **Not yet measured as a selector — logged as the one concrete
   Type I-specific follow-up this finding produced.**
6. **Fix the label in the biology map.** 72/15 → **73/14**; PK9 (4D6Z) coordinates at
   2.379 Å. Three places in `CYP3A4_BIOLOGY_MAP.md` cite the old split. The argument the map
   builds on it — that both sampling frames are biased toward coordinating chemistry — is
   unaffected: 84% instead of 83%.

**What this does NOT license.** It does not say the selector is mode-independent on
*CYP3A4*; n=14 cannot support that. It does not extend to peripheral binders — the family
set contains none (max `closest_fe` 6.26 Å) and CYP3A4 none, so ligands sitting outside the
active site are untested by anything here. And it does not address the multi-conformer
ligands the challenge announcement describes; this is one pose per ligand throughout.

---

## Artefacts

| what | where |
|---|---|
| script, three stages | `scripts/structure/binding_mode_robustness.py` |
| verified labels, both sides, 87 rows | `data/processed/binding_mode_labels_cyp3a4.csv` |
| CYP3A4 strata, comparison, frame tests | `data/processed/binding_mode_strata_cyp3a4.json` |
| CYP3A4 per-ligand gains | `data/processed/binding_mode_per_ligand_cyp3a4.csv` |
| family xeng cache (lets the poses go back to cold storage) | `data/processed/binding_mode_p450_xeng.csv` |
| family per-pair gains | `data/processed/binding_mode_per_pair_p450.csv` |
| family strata, matched comparison, controls | `data/processed/binding_mode_strata_p450.json` |

No downloads, no inference, no GPU. Intermediates went to `C:\Temp`, never to `D:`.

---

# Correction, 2026-09-23 — the n = 14 floor quoted here is right, and it is **not** the one 034/035/036 quote

Appended by `FINDING_039`. The body is unchanged; this note records what was re-measured.

**The number in this finding is correct.** The `+0.0433 / +0.0622` floor in the Type I column
is computed on the **crystal-side** Type I stratum (`crystal_mode == type_I`, random baseline
**0.4644**). Recomputed at 2,000,000 draws, that population's floor is
**p95 +0.04306 · p99 +0.06160** (SE ±0.00004). The published +0.0433 is **+0.27 sd** of a
4,000-draw estimate away from it — a good draw.

**It is a different quantity from the one FINDINGs 034, 035 and 036 call "the n = 14 floor".**
Those use the **prediction-side** stratum, which swaps `D0R` in and `QDY` out and moves the
random baseline to 0.4597. Its floor is **+0.04490** at depth 20 and **+0.04400** at depth 40.
Do not cross-quote. `FINDING_039` §1 has the full mapping.

**The verdict is unaffected.** Type I's +0.0411 against the authoritative crystal-side floor
gives an exact one-sided **p = 0.0581** (2 × 10⁶ draws), against the +0.0595 published here.
The stratum still does not clear its own floor and the conclusion — *the measurement has no
resolution at n = 14*, not *the selector fails* — stands exactly as written.
