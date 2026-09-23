# FINDING 036 — the near-tie is a variance mechanism, not a prize

**Date:** 2026-09-23 · **Status:** measured, zero new inference, zero downloads, CPU-only,
one sitting · **Verdict:** **MEASURED** at the pre-registered stage-1 stop — **a perfect
top-2 tie-break is worth +0.0129 against a +0.0134 noise floor recomputed on this pool.**
The five candidates were run past the stop anyway and are **REFUTED**: best paired
+0.0011, p = 0.68, and two of the five are significantly *negative*.

Pre-registered in `docs/PREREG_tiebreak.md`, committed at `a8ae652` before a single
endpoint existed. One deviation, declared in §5.

---

## What was asked

`FINDING_035` doubled the pool at genuinely matched conditioning and found Δ selected
**+0.0116**, Wilcoxon **p = 0.715**, **10 of 14 ligands tied**. It named the mechanism:

> A deeper pool does not give the selector a better pose to find, it gives it more
> near-ties — and a near-tie is a coin flip with an enormous payoff.

Both swing ligands turned on a **sub-0.06 Å Chamfer margin** producing a **> 0.3 LDDT-PLI**
swing (CFF +0.4826 on 0.057 Å, 08J −0.3442 on 0.045 Å), and 035 item 5 named the successor:
*"A tie-break for the top-2 `xeng` poses is a better-evidenced next experiment than more
depth — and it is measurable on this very pool, with no new inference."*

> **The question, in the pre-registered order.** First: **how much score is actually at
> stake inside the top-k**, i.e. what would a *perfect* tie-break be worth? Only then:
> can any selector-native rule claim it?

---

## 1. Controls — all pass, with counts including the zeros

| id | control | measured |
|---|---|---|
| **C-XENG** | recompute `xeng` from the 1,740 pool mmCIFs + the frozen reference | max abs diff **2.45e-07**, 1,740 / 1,740 rows joined, 0 left-only, 0 right-only |
| **C-XENG-B** | same on the 280 `FINDING_035` Explorer poses | max abs diff **1.92e-07**, 280 rows |
| **C-SEL** | `argmin(xeng)` **is** `cypstruct.xengine.select()` | max abs diff **0.0**, equal to 1e-9 on all **87**. Re-verified, not inherited |
| **C-N2** | the shipped board reproduces | selected **0.6164**, oracle **0.6975**, random **0.5769**, gain **+0.0395**, ρ **−0.2582**, correct sign **75.86%** — `FINDING_011`/`033`/`035` to four decimals |
| **C-NUM** | `FINDING_021` numbering | `mapped` True on **1,740 of 1,740**; **4** exact-zero LDDT-PLI rows (PG4 ×3, PG0 ×1) at BiSyRMSD **23.1–26.9 Å** → **ejections, not numbering**, the same four `FINDING_032` found |
| **C-REF** | reference depth | min **6**, median **7**, max **12**; all ≥ the 4 `FINDING_011` requires |
| **C-XTAL** | crystals for the R1 gate | **87 of 87** loaded, one chain each, heme frame constructible, **0 missing** |
| **C-FILT** | every filter fires | 87 validation ligands → **14** predicted-Type-I, **73** excluded; 1,740 pool-A poses, **1,740** with features; 280 pool-B poses; **0** dropped for no reference, **0** for no truth |

One control read **0 rows** on its first run — pool B's files are named `<LIG>__<sample>.cif`
while the truth table keys on `<sample>` — and was fixed rather than accepted. A control
that silently measures nothing is `too-clean-numbers-are-the-tell` in its cheapest form.

---

## 2. The prize, measured before any tie-breaker existed

Pool A: **87 ligands × 20 unsteered Boltz-2 poses**. Poses sorted ascending by `xeng`
(lower is better). `stake_k` = best LDDT-PLI in the top-k minus the incumbent's pick.

**Pool oracle 0.6975 · incumbent 0.6164 · random 0.5769 · oracle gap 0.0811.**

### 2a. The near-tie is not a regime. It is the whole pool.

| `margin_12` = `xeng(2) − xeng(1)` | p5 | p25 | **p50** | p75 | p95 |
|---|---|---|---|---|---|
| Å of Chamfer | 0.0013 | 0.0097 | **0.0194** | 0.0493 | 0.1473 |
| in within-ligand sd | 0.006 | 0.067 | **0.146** | 0.426 | 1.075 |

**70 of 87 ligands** have a top-2 margin below **0.06 Å** — the number `FINDING_035`
flagged as razor-thin on its two swing ligands. **77 of 87** are below 0.10 Å. The median
top-1/top-2 decision is made on **0.15 within-ligand standard deviations**.

So 035's observation is correct and is *more* general than it looked: the selector is
almost always choosing between two poses it cannot meaningfully distinguish. That framing
predicts a large prize. It is wrong.

### 2b. The oracle tie-break — the entire budget for this line of work

| k | **ORACLE tie-break** | median | worst-case tie-break | ligands with stake > 0 | > 0.05 | share of the 0.0811 oracle gap |
|---|---|---|---|---|---|---|
| **2** | **+0.0129** | 0.0000 | −0.0451 | 33 / 87 | 7 | **15.9%** |
| **3** | **+0.0210** | 0.0072 | −0.0552 | 53 / 87 | 11 | 25.8% |
| **5** | **+0.0318** | 0.0113 | −0.0670 | 60 / 87 | 19 | 39.2% |

**A perfect top-2 tie-break is worth +0.0129 LDDT-PLI per ligand.** The random-feature
noise floor, recomputed on this pool over 4,000 draws, is **p95 +0.0134 / p99 +0.0193**
(`FINDING_007` measured +0.0137 / +0.0193 — reproduced).

**The ceiling of the entire family sits below its own noise floor.** The pre-registered
stop rule fires.

### 2c. And the incumbent already takes most of what is there

The number that makes the ceiling small is not that the top-2 poses are similar — they are
not, the median stake among the 33 ligands that have one is large. It is that
**`argmin(xeng)` is already the right one of the two far more often than a coin flip.**

| | pool A, k = 2 |
|---|---|
| random pick among the top 2 | **0.6003** |
| **incumbent (`argmin`)** | **0.6164** |
| perfect tie-break | 0.6293 |
| worst-case tie-break | 0.5713 |

The incumbent sits **78% of the way** from the worst possible top-2 choice to the best. It
beats the coin flip by **+0.0161**, which is *more than the entire remaining prize*
(+0.0129). The matched null says the same thing without an assumption: over **4,000**
random tie-breaks of the same top-2, the mean is **−0.0161** and even the **95th
percentile is −0.0071** — not one draw in twenty beats simply keeping the argmin.

**A 0.0194 Å median margin is not noise.** The shipped Chamfer resolves the choice at a
scale far below where it looked meaningful.

### 2d. The stake does not live in the tight margins

If near-ties were where the score is decided, `stake_2` would rise as the margin falls. It
does not.

| margin quartile | mean margin (Å) | mean `stake_2` | mean downside |
|---|---|---|---|
| Q1 tightest | 0.0045 | **0.0065** | 0.0243 |
| Q2 | 0.0151 | 0.0227 | 0.0125 |
| Q3 | 0.0322 | 0.0076 | 0.0515 |
| Q4 widest | 0.1052 | 0.0146 | 0.0925 |

ρ(margin, `stake_2`) = **−0.1436, p = 0.18**. No trend, in either direction. The
**tightest** quartile has the **smallest** stake — when two poses are that close in the
heme frame they are usually close in LDDT-PLI too, which is the feature behaving exactly
as designed.

### 2e. The `FINDING_035` stratum, at depth 40, says the same thing with n = 14

| | pool B, 14 predicted-Type-I ligands × **40** | same ligands × **20** |
|---|---|---|
| oracle / incumbent / random | 0.6595 / 0.5386 / 0.4587 | 0.6499 / 0.5270 / 0.4597 |
| median `margin_12` | 0.0363 Å | 0.0705 Å |
| **ORACLE top-2 tie-break** | **+0.0447** | +0.0164 |
| noise floor **inside n = 14**, 4,000 draws | **p95 +0.0435 / p99 +0.0623** | |
| incumbent − random-in-top-2 | **−0.0153** | **+0.0334** |

Two things worth separating. **The ceiling rises with depth** (+0.0164 → +0.0447) — 035's
"depth buys near-ties" is confirmed quantitatively, the margin halves and the stake almost
triples. And **the deepened pool is the one place where `argmin(xeng)` is worse than a coin
flip** (−0.0153): at depth 40 on this stratum the top-1/top-2 decision genuinely is a coin
flip, which is exactly the mechanism 035 described. But the ceiling still lands **on** the
n = 14 floor (+0.0447 against +0.0435), so even there nothing is provable.

---

## 3. R1 — the answer-recognition gate, and a correction it forces on the repo

The pre-registration carried an a-priori calibration branch: the shipped `xeng` is scored
on the identical gate, because every candidate is a function of the same pose–reference
distance matrix `xeng` averages. **That branch fired, and it is the most consequential
thing in this finding.**

Percentile of the query's own **crystal ligand** among its 20 predicted poses (1.0 = the
crystal is the single best-scoring object in the set), 87 ligands, all 87 crystals loaded:

| feature | mean percentile | median | crystal beats the median pose | binomial p | passes ≥ 0.65? |
|---|---|---|---|---|---|
| **`xeng` — THE SHIPPED SELECTOR** | **0.3437** | 0.20 | **28 / 87** | **0.0012** | **no** |
| `xeng_2nd` | 0.3259 | 0.15 | 23 / 87 | 1.3e-05 | no |
| `xeng_borda` | 0.3471 | 0.25 | 30 / 87 | 0.0050 | no |
| `xeng_sd` | 0.5615 | 0.65 | 48 / 87 | 0.39 | no |
| `xeng_votes` | 0.2822 | 0.00 | 27 / 87 | 5.2e-04 | no |

**The shipped selector — +0.0395 on this pool, +0.0357 across 81 held-out P450 proteins —
puts the true answer at the 34th percentile of its own predictions, significantly below
chance.** That is the same signature `FINDING_032` used to close the internal conformer
(31st percentile, p = 2.6e-04) and `FINDING_030` reported at chance (56th, p = 0.44).

**So R1 is not a valid gate for this family, and the pre-registration said so in advance.**
A gate that kills a feature which demonstrably selects, on the very pool it selects on, is
not measuring selection value.

### Why, and what the gate is actually good for

The mechanism is not mysterious and it is already in the file. Cross-engine agreement
measures distance to the *consensus of independent co-folders*, and `FINDING_005`/`024`
established that those engines **share CYP3A4's error** — Chai ρ = +0.45, Protenix
ρ = +0.47–0.60 against Boltz, and the error is a 30° rotation every engine makes. The
consensus is therefore **displaced from the truth**, so the truth is an outlier to it. What
`xeng` exploits is that *among the predictions*, moving toward the displaced consensus is
still moving toward the answer. A **biased but informative** estimator.

**The carry-forward, which supersedes the RUNBOOK's current wording:** R1 tests whether a
feature can **recognise the truth as an object**. That is the right test for a feature
whose claim is *"the crystallographic record constrains where this goes"* — a **prior**,
which is what 030 and 032 both were, and both verdicts stand unchanged. It is the **wrong**
test for a feature whose claim is *"these predictions can be ordered relative to one
another"* — a **within-ligand comparator**, where the reference is allowed to be biased as
long as the bias is common to the pool. Apply R1 to priors. Do not apply it to comparators.
`FINDING_032`'s **term oracle** rung is unaffected and remains the stronger instrument:
it is what closed this experiment, four hours before any candidate was scored.

---

## 4. The candidates, and the three bars

### 4a. The deviation, declared

The stop rule fired at k = 2 and the pre-registration says *"no tie-breaker is built."* The
candidates were nevertheless computed, because the k = 3 (+0.0210) and k = 5 (+0.0318)
ceilings **do** clear the +0.0134 floor, and because the computation is free on a cached
distance matrix. **Everything in §4 is therefore post-stop and is not eligible to ship
whatever it shows.** It is reported as confirmation, not as a test.

### 4b. Candidates, all native to the selector, no new information, signs fixed a priori

| id | rule | sign | oriented statistic |
|---|---|---|---|
| **T1** | medoid of the top-k (mean pose-to-pose Chamfer inside the top-k) | lower better | `med_k`, undefined at k = 2 |
| **T2** | distance to the **2nd-nearest** reference, not the mean | lower better | `xeng_2nd` |
| **T3** | **Borda** mean rank across references | lower better | `xeng_borda` |
| **T4** | within-pose **sd** of the Chamfer across references | lower better | `xeng_sd` |
| **T5** | **plurality vote** — references whose single nearest pose this is | higher better | `−xeng_votes` |

### 4c. BAR 1 — versus random, nulls recomputed for these features on this pool

| null, 4,000 draws | pool A (n = 87) | pool B (n = 14) |
|---|---|---|
| random-feature gain, p95 / p99 | **+0.0134 / +0.0193** | **+0.0435 / +0.0623** |
| random tie-break vs incumbent, k = 2, mean / p95 | **−0.0161 / −0.0071** | +0.0151 / +0.0365 |
| random tie-break vs incumbent, k = 3 | −0.0192 / −0.0098 | −0.0159 / +0.0276 |
| random tie-break vs incumbent, k = 5 | −0.0188 / −0.0095 | −0.0334 / +0.0134 |

The pooled floor reproduces `FINDING_007` (+0.0137) and the n = 14 floor reproduces
`FINDING_035` (+0.0453) — both recomputed here rather than quoted, and each is used only
for its own set.

### 4d. BAR 2 — paired against `cypstruct.xengine.select()`, identical poses, n = 87

Pool oracle **0.6975**, incumbent **0.6164**, random **0.5769**.

| candidate | k | selected | gain vs random | **paired Δ** | 95% CI | Wilcoxon p | **picks CHANGED** | **TIED** |
|---|---|---|---|---|---|---|---|---|
| **T3 Borda** | **2** | 0.6175 | +0.0406 | **+0.0011** | [−0.0036, +0.0059] | 0.679 | **16** | **71** |
| T3 Borda | 3 | 0.6160 | +0.0391 | −0.0004 | [−0.0064, +0.0052] | 0.884 | 22 | 65 |
| T3 Borda | 5 | 0.6147 | +0.0378 | −0.0017 | [−0.0099, +0.0053] | 0.951 | 23 | 64 |
| T2 2nd-nearest | 3 | 0.6134 | +0.0365 | −0.0030 | [−0.0121, +0.0049] | 0.485 | 26 | 61 |
| T2 2nd-nearest | 2 | 0.6104 | +0.0335 | −0.0060 | [−0.0143, +0.0010] | 0.122 | 21 | 66 |
| T5 plurality | 2 | 0.6102 | +0.0333 | −0.0062 | [−0.0130, −0.0008] | 0.149 | 27 | 60 |
| T5 plurality | 5 | 0.6089 | +0.0320 | −0.0075 | [−0.0156, −0.0010] | 0.122 | 33 | 54 |
| T2 2nd-nearest | 5 | 0.6088 | +0.0319 | −0.0076 | [−0.0187, +0.0025] | 0.187 | 28 | 59 |
| T4 reference sd | 2 | 0.6077 | +0.0308 | −0.0087 | [−0.0191, +0.0007] | 0.093 | 37 | 50 |
| T1 medoid | 5 | 0.6013 | +0.0244 | −0.0151 | [−0.0332, +0.0013] | 0.390 | 62 | 25 |
| T1 medoid | 3 | 0.5988 | +0.0219 | −0.0176 | [−0.0344, −0.0032] | 0.215 | 48 | 39 |
| **T4 reference sd** | 3 | 0.5986 | +0.0217 | **−0.0178** | [−0.0338, −0.0031] | **0.025** | 59 | 28 |
| **T4 reference sd** | 5 | 0.5961 | +0.0192 | **−0.0203** | [−0.0385, −0.0044] | **0.043** | 70 | 17 |
| **T4 reference sd** | *standalone, all 20* | 0.5619 | **−0.0150** | **−0.0545** | [−0.0854, −0.0260] | **0.001** | 80 | 7 |

**Fifteen of sixteen rows are negative.** The single positive is **+0.0011 at p = 0.679**
with a CI that straddles zero and **71 of 87 picks unchanged** — a null with 16 degrees of
freedom. No row approaches the +0.0134 floor, let alone the SHIPS bar.

### 4e. On exactly the ligands whose pick CHANGES — reported separately, as required

A rule that moves 16 of 87 picks cannot be judged by a pooled mean.

| candidate | k | changed | better | worse | mean Δ on the changed | incumbent → candidate, those ligands |
|---|---|---|---|---|---|---|
| T3 Borda | 2 | 16 | **8** | **8** | +0.0061 | 0.5484 → 0.5544 |
| T2 2nd-nearest | 2 | 21 | 8 | 13 | −0.0247 | 0.6239 → 0.5992 |
| T5 plurality | 2 | 27 | 12 | 15 | −0.0200 | 0.6272 → 0.6072 |
| T4 reference sd | 2 | 37 | 13 | 24 | −0.0205 | 0.6253 → 0.6048 |
| T1 medoid | 3 | 48 | 21 | 27 | −0.0319 | 0.6107 → 0.5788 |
| T4 reference sd | *standalone* | 80 | 29 | 51 | −0.0593 | 0.6196 → 0.5603 |

The best candidate moves 16 picks and wins **8 of them**. That is a coin flip, stated
without a p-value.

### 4f. BAR 3 — ranking, beside every pooled mean, and `FINDING_011`'s trap again

| feature | within-ligand ρ | correct sign | binomial p | best paired Δ |
|---|---|---|---|---|
| **`xeng_borda`** | **−0.2671** | **77.01%** | 4.3e-07 | **+0.0011** |
| `xeng` (incumbent) | −0.2582 | 75.86% | 1.4e-06 | — |
| `xeng_2nd` | −0.2303 | 71.26% | 9.1e-05 | −0.0030 |
| `xeng_votes` | −0.1083 | 65.52% | 0.0050 | −0.0062 |
| `xeng_sd` | −0.0305 | 51.72% | 0.83 | −0.0087 |

**Borda ranks the pool better than the shipped feature does** — ρ −0.2671 against −0.2582,
correct sign on 77.0% against 75.9% — **and buys +0.0011.** That is the third time in this
repo that ρ and top-1 selection have separated (`FINDING_011`, twice). Judging on ρ would
have shipped a change worth nothing.

`xeng_sd`, the one candidate with an independent rationale rather than a re-weighting of
`xeng`, ranks at **chance** (ρ −0.031, 51.7%) and is the most harmful of the five.

### 4g. The n = 14 stratum, and what a 14-ligand board is made of

| candidate | k | selected | paired Δ | 95% CI | p | changed | tied |
|---|---|---|---|---|---|---|---|
| **T4 reference sd** | 2 | 0.5778 | **+0.0392** | [−0.0015, +0.0961] | 0.263 | 8 | 6 |
| T5 plurality | 3 | 0.5540 | +0.0154 | [−0.0026, +0.0469] | 0.273 | 4 | 10 |
| T3 Borda | 2 | 0.5443 | +0.0057 | [−0.0036, +0.0207] | 0.655 | 2 | 12 |
| T4 reference sd | 5 | 0.5010 | −0.0376 | [−0.1301, +0.0356] | 0.790 | 11 | 3 |
| T4 reference sd | *standalone* | 0.3950 | −0.1436 | [−0.2767, −0.0148] | 0.135 | 14 | 0 |

`+0.0392` is below the n = 14 floor (**+0.0435**), its CI includes zero, p = 0.263 — and it
is **one ligand**. The post-hoc trace, on 035's own two swing cases, top-3 by `xeng`:

| | incumbent | oracle in top-3 | T2 | T3 | **T4 sd** | T5 | T1 medoid |
|---|---|---|---|---|---|---|---|
| **CFF** | **0.9547** | 0.9547 | 0.9349 | 0.4721 | 0.4721 | 0.9349 | 0.9547 |
| **08J** | **0.1364** | 0.4806 | 0.1364 | 0.1364 | **0.4806** | 0.1364 | 0.1579 |

At k = 2, `xeng_sd` keeps CFF and rescues 08J — **+0.344 on one ligand out of fourteen**,
which is the whole of its +0.0392. The identical rule on 87 ligands scores **−0.0087**
(k = 2), **−0.0178** (k = 3, p = 0.025) and **−0.0203** (k = 5, p = 0.043). This is
`within-ligand-variance-is-what-matters` and `cyp-selector-noise-floor` in one table: the
n = 14 stratum cannot distinguish a rule from a ligand.

---

## 5. Verdict — **MEASURED**, and the family is **REFUTED**

Against the pre-registered table:

| clause | required | measured | fires? |
|---|---|---|---|
| **stage-1 stop** | O2 < pooled p95 floor | **+0.0129** vs **+0.0134** | **YES** |
| SHIPS | paired ≥ +0.0137 **and** CI excludes 0 **and** p < 0.05 | best **+0.0011**, CI [−0.0036, +0.0059], p = 0.679 | no |
| REFUTED (as written) | best paired ≤ 0 **and** CI upper < +0.0137 | best paired is **+0.0011** (> 0) | no, on a technicality |
| **MEASURED** | anything else, **including the stage-1 stop** | | **YES** |

**MEASURED** is the pre-registered verdict and the stage-1 stop is the clause that fires.
The candidate half is a **refutation in substance** — five rules, sixteen configurations,
fifteen negative, the one positive a p = 0.68 coin flip on 16 of 87 picks, and two rules
significantly harmful — and it misses the literal REFUTED clause only because the best of
sixteen post-stop numbers happens to be +0.0011 rather than −0.0011.

### What it establishes

1. **A perfect top-2 tie-break is worth +0.0129 and a perfect top-3 +0.0210**, against a
   noise floor of +0.0134 / p99 +0.0193 recomputed on this pool. **The ceiling of the whole
   family is at or below its own floor.** This closes tie-breaking the way `FINDING_032`
   closed the internal conformer — by the **oracle**, not by a gate. No tie-breaker of any
   shape, from any information source, can clear the floor on this pool at k = 2.
2. **`FINDING_035`'s near-tie observation is correct, more general than it looked, and does
   not imply a prize.** The median top-2 margin is **0.0194 Å** and **70 of 87** ligands sit
   under the 0.06 Å that looked razor-thin — the near-tie is the normal case, not a regime.
   But the stake does **not** concentrate there (ρ = −0.14, p = 0.18; the tightest quartile
   has the *smallest* stake), so the near-ties are a **variance mechanism** — they explain
   why depth adds symmetric coin flips with large payoffs and therefore does not convert —
   and **not** a recoverable prize. That is the branch the pre-registration wrote down in
   advance for exactly this outcome.
3. **The 0.0194 Å margin is signal, not noise.** `argmin(xeng)` beats a coin flip over the
   top 2 by **+0.0161** — more than the entire remaining prize — and beats **4,000 of 4,000**
   random tie-breaks with the 95th percentile still at −0.0071. The shipped Chamfer resolves
   the choice at a scale two orders of magnitude below the ligand error it is ranking.
4. **The answer-recognition gate does not apply to within-ligand comparators.** The shipped
   selector puts the crystal at the **34th percentile** of its own predictions, p = 0.0012 —
   `FINDING_032`'s exact signature, on the one feature in this repo that works. The gate
   tests whether a **prior** recognises the truth, and 030's and 032's verdicts are
   untouched; it is invalid for a feature whose claim is *ordering predictions against each
   other*, because the co-folders' shared 30° error (`FINDING_024`) displaces the consensus
   from the answer while leaving it informative about the ordering. Pre-declared as a
   branch; reported because it fired.
5. **Where the 0.0811 oracle gap actually lives: outside the top-k.** Only **15.9%** of it
   is reachable inside the top 2, 25.8% inside the top 3, 39.2% inside the top 5. **Over
   60% of the gap requires promoting a pose the incumbent ranks 6th or worse** — a
   different problem from tie-breaking and the only one left with room in it.

### What it forecloses

- **Tie-breaking the top-2 by `xeng`, by any rule.** Closed by the term oracle.
- **The pose–reference distance matrix as a source of anything beyond its mean.** Five
  different reductions of the same `D[pose, ref]` — second order statistic, Borda rank,
  within-pose sd, plurality vote, pool-internal medoid — and the plain mean beats all five.
  `FINDING_011` retired the *combination* `-zm - zx`; this retires the *decomposition*.
- **Reading a small margin as an absence of information.** It cost this experiment nothing
  to check and it reversed the premise.
- **Using n = 14 to grade a selection rule.** The best number on that board, +0.0392, is one
  ligand, and the same rule is significantly negative at n = 87.

### What it does NOT license

- It does not say the near-ties are harmless: the **worst-case** top-2 tie-break is
  **−0.0451**, three times the best case. A selector change that randomised the top-2 would
  cost real score. The asymmetry is why the incumbent must be left alone.
- It tests **one pool, one reference set**. A pool with a genuinely bimodal top-2 — the
  `FINDING_012` catastrophe regime — would have a larger stake, and the blind release may be
  such a pool. **There is no drop-day proxy for the ceiling and none is claimed:** `stake_k`
  needs truth, and the one quantity that is computable without it — the margin distribution
  — does not predict the stake here (ρ = −0.14, p = 0.18).
- The n = 14 stratum's ceiling (+0.0447 at depth 40) sits **on** its floor (+0.0435), so
  that set neither confirms nor refutes; it is quoted, not concluded from.

### What to do instead

`FINDING_035` item 5 said the near-tie was better evidenced than more depth. It was worth
measuring and the measurement says no. The remaining **0.081** per ligand is 60% outside the
top 5, which means the question is not *"which of these two"* but *"why is the right pose
ranked 8th"* — and every candidate that has failed since `FINDING_025` either scored a pose
from its own geometry (025, 027, 029, 030, 031, 032) or, as here, re-read the same
consensus. `docs/README.md`'s standing observation holds and is now one finding deeper.

---

## 6. Reproduce

```
python scripts/structure/tiebreak_analysis.py frames     # 2,020 poses + 87 crystals, ~3 min
python scripts/structure/tiebreak_analysis.py controls   # C-XENG C-SEL C-N2 C-NUM C-REF C-XTAL
python scripts/structure/tiebreak_analysis.py regime     # STAGE 1 -- the prize, and the stop
python scripts/structure/tiebreak_analysis.py gate       # R1, incl. the xeng calibration
python scripts/structure/tiebreak_analysis.py select     # candidates + the three bars
python scripts/structure/tiebreak_analysis.py swing      # post-hoc, CFF and 08J
```

| what | where |
|---|---|
| pre-registration | `docs/PREREG_tiebreak.md` (`a8ae652`) |
| script, six stages | `scripts/structure/tiebreak_analysis.py` |
| controls | `data/processed/tiebreak_controls.json` |
| **stage 1, the prize** | `data/processed/tiebreak_regime.json` |
| R1 answer recognition | `data/processed/tiebreak_answer_recognition.json` |
| candidates, nulls, bars | `data/processed/tiebreak_selectors.json` |
| per-ligand margins and stakes | `data/processed/tiebreak_per_ligand.csv` |
| post-hoc swing trace | `data/processed/tiebreak_swing.json` |

Heme-frame caches and the pose–reference distance matrices are intermediates and live in
`C:\Temp\cyp_tiebreak` (1.1 MB), never on `D:`. No inference was run and nothing was
downloaded.
