# FINDING 037 — the right pose IS a minority report, and the majority vote is still the right bet

**Date:** 2026-09-23 · **Status:** measured, zero new inference, zero downloads, CPU-only,
one sitting · **Verdict:** **MEASURED** at the pre-registered clause. The **mechanism is
CONFIRMED on all four pre-registered predictions** — `xeng` ranks by orientation-cluster
membership (ρ = −0.25 to −0.34, correct on 69–85% of ligands), and where the stake is real
the best pose sits in a **minority** cluster **twice as often as chance** (46.7% vs 23.1%,
p < 1e-4). The **rescue is REFUTED in substance**: a **perfect** *"the consensus is wrong
here, take the minority cluster"* rule is worth **+0.0427** at a cut that makes a genuine
partition, of which only **+0.0135** survives a size-matched shuffle — and the same
clustering shows that picking the **majority** cluster is worth **+0.0546**, *more than the
entire perfect rescue*, at p = 1.3e-06.

Pre-registered in `docs/PREREG_minority_report.md`, committed at `6f547ef` before a single
endpoint existed. One methods correction and four post-hoc diagnostics, all declared in §6.

---

## What was asked

`FINDING_036` closed tie-breaking by its own term oracle — a **perfect** top-2 tie-break is
worth **+0.0129** against a **+0.0134** floor — and named the successor in its own words:

> **60% of the 0.081 oracle gap sits outside the top 5.** The question is not *"which of
> these two"* but *"why is the right pose ranked 8th?"*

036 also supplied the mechanism to test first. The shipped selector fails the
answer-recognition gate (the crystal sits at the **34th percentile** of its own
predictions, p = 0.0012) while still selecting at **+0.0395**, because the co-folders
**share** CYP3A4's 30° orientation error (`FINDING_024`). If that is right, then when every
engine makes the same mistake the best pose in the pool is a **minority report** — an
outlier the consensus actively down-ranks — and `xeng` is majority-voting for a shared
error.

Four measurements, in the pre-registered order: **where the oracle pose sits**; **whether
it is a minority report**; **when the consensus is confidently wrong**; and **what a
perfect minority-rescue rule would be worth**.

---

## 1. Controls — all pass, with counts including the zeros

| id | control | measured |
|---|---|---|
| **C-XENG** | recompute `xeng` from the cached pose–reference Chamfer matrix | max abs diff **2.45e-07**, **1,740 / 1,740** joined, 0 left-only, 0 right-only |
| **C-SEL** | `argmin(xeng)` **is** `cypstruct.xengine.select()` | max abs diff **0.0**, equal to 1e-9 on all **87**. Re-verified, not inherited |
| **C-N2** | the shipped board reproduces | selected **0.6164**, oracle **0.6975**, random **0.5769**, gain **+0.0395**, ρ **−0.2582**, correct sign **75.86%** — `FINDING_011`/`033`/`035`/`036` to four decimals |
| **C-NUM** | `FINDING_021` numbering | `mapped` True on **1,740 of 1,740**; **4** exact-zero LDDT-PLI rows (PG4 ×3, PG0 ×1) at BiSyRMSD **23.1–26.9 Å** → ejections, not numbering. The same four 032 and 036 found |
| **C-REF** | reference depth | min **6**, median **7**, max **12**; all ≥ the 4 `FINDING_011` requires |
| **C-ATOM** | every pose of a ligand has one atom count — required before the rotation metric may use identity correspondence | **87 of 87** homogeneous, 0 heterogeneous |
| **C-FILT** | every filter fires | 87 validation ligands → 1,740 pool poses → 1,740 with features → **87 ligands assembled, 20 poses each**; **0** dropped for no reference, **0** for no truth; stake split **45 HIGH / 42 LOW** |

`FINDING_036`'s top-k ceilings also reproduce to the digit as a by-product of step 1
(**k2 +0.0129 / 15.92%**, **k3 +0.0210 / 25.84%**, **k5 +0.0318 / 39.23%**), which is the
cheapest available confirmation that this is the same pool and the same feature.

---

## 2. Step 1 — where the oracle pose sits. The answer is: **rank 8**

Poses sorted ascending by `xeng`. `rank_oracle` is the 1-indexed position of the
best-LDDT-PLI pose; ties at the maximum take the **smallest** rank, which is conservative
against this hypothesis.

Pool oracle **0.6975** · incumbent **0.6164** · random **0.5769** · **oracle gap 0.0811**.

### 2a. Pooled, the distribution is almost exactly uniform

| | mean rank | median | rank 1 | top 3 | top 5 | bottom half (>10) | dead last (20) |
|---|---|---|---|---|---|---|---|
| **all 87** | **9.43** | **8.0** | 14 (16.1%) | 23 (26.4%) | 30 (34.5%) | 34 (39.1%) | 5 (5.7%) |
| uniform null | 10.5 | 10.5 | 5% | 15% | 25% | 50% | 5% |

Wilcoxon against the uniform median, **p = 0.11**. So the literal question in 036's last
paragraph has a literal answer: **the median oracle pose is ranked 8th**, and pooled, the
ranking of the best pose is barely distinguishable from a lottery.

### 2b. Split by what is at stake, it is two completely different distributions

`stake_all` = per-ligand (oracle − incumbent's pick). HIGH = > 0.05, fixed in advance.

| stratum | n | mean stake | mean rank | median | rank 1 | top 3 | top 5 | bottom half | dead last | p vs uniform |
|---|---|---|---|---|---|---|---|---|---|---|
| **LOW** (≤ 0.05) | 42 | 0.0174 | **6.62** | **3.5** | 14 (33.3%) | 21 (50.0%) | 23 (54.8%) | 9 (21.4%) | 1 | **4.4e-04** |
| **HIGH** (> 0.05) | 45 | 0.1406 | **12.04** | **11.0** | 0 (0%)\* | **2 (4.4%)** | **7 (15.6%)** | 25 (55.6%) | 4 (8.9%) | 0.063 |

\* **rank 1 forces stake 0 by construction**, so the 0 of 45 is definitional and is
reported as such, not as evidence. Everything else in that row is not: **2 of 45** in the
top 3 against a uniform 15%, **7 of 45** in the top 5 against 25%, and a mean rank of
**12.04 against 10.5** — a *trend toward worse than a lottery*, p = 0.063.

Terciles of the stake say the same thing monotonically:

| tercile | n | mean stake | mean rank | median | top 3 | dead last |
|---|---|---|---|---|---|---|
| T1 least | 29 | 0.0027 | **6.00** | 1.0 | 19 (65.5%) | 1 |
| T2 | 29 | 0.0545 | 10.35 | 9.0 | 3 (10.3%) | 0 |
| **T3 most** | 29 | **0.1805** | **13.03** | **14.0** | **1 (3.4%)** | **4 (13.8%)** |

**ρ(rank, stake) = +0.601, p = 7.4e-10.** The stake-weighted mean rank is **12.15**, worse
than the 10.5 a lottery would give. Post-hoc de-confounding (§6, A1): restricted to the 73
ligands with any stake at all, ρ = **+0.327, p = 0.0047** — smaller, and still there.

**The one-line reading.** The incumbent finds the best pose reliably *when finding it is
worth nothing* (median rank 3.5, 33% at rank 1) and buries it *when it is worth 0.14
LDDT-PLI* (median rank 11, 2 of 45 in the top 3). That is not a fine-ranking failure. It is
the ranking inverting exactly where the score lives, and it is the strongest form of 036's
"60% outside the top 5" yet measured.

---

## 3. Step 2 — is the oracle pose a minority report? **Yes, on all four predictions**

Poses clustered per ligand in their own heme frame, average linkage. Primary distance
**pose-to-pose Chamfer** (permutation-invariant, the same metric the selector uses);
secondary **`FINDING_024`'s centroid-aligned rotation angle** (Kabsch, identity
correspondence, legal by C-ATOM). Eight Chamfer cuts and six rotation cuts, **all fixed in
advance and all reported**. Nothing is tuned.

### 3a. P2 — `xeng` ranks by cluster membership. Confirmed at every cut.

Within-ligand ρ(cluster size, `xeng` rank) — negative means a bigger cluster gets a better
rank:

| cut | mean clusters | **ρ(size, rank)** | correct sign | ligands defined | p |
|---|---|---|---|---|---|
| Chamfer 0.50 Å | 12.71 | **−0.246** | 75.3% | 77 | 9.8e-06 |
| Chamfer 0.75 Å | 6.59 | **−0.306** | 76.0% | 79 | 4.2e-06 |
| **Chamfer 1.00 Å** | **3.55** | **−0.341** | **79.7%** | 74 | **2.6e-07** |
| Chamfer 1.25 Å | 2.10 | −0.250 | 68.9% | 45 | 0.016 |
| Rotation 15° | 4.13 | **−0.309** | **81.1%** | 74 | **6.2e-08** |
| Rotation 30° | 2.61 | −0.291 | 78.3% | 60 | 1.2e-05 |
| Rotation 45° | 2.08 | −0.275 | 81.3% | 48 | 1.5e-05 |
| Rotation 90° | 1.60 | −0.276 | 85.0% | 40 | 8.4e-06 |

**Every cut, both metrics, the predicted sign, 69–85% of ligands.** The majority-vote
mechanism is not an inference from 024 — it is measured directly, and it is about the same
magnitude as the selector's own within-ligand ρ of −0.258. `xeng` and "how many other poses
agree with this one" are close to the same quantity.

### 3b. P1 and P3 — and a double dissociation

**A methods correction, declared.** The pre-registration set the chance rate for E2a as
`mean(size_largest / n)`. That **undercounts whenever several clusters tie for largest**,
which is the whole fine-cut regime. The correct null replaces the oracle pose by a **random
pose of the same ligand** and recomputes the rate (4,000 draws, §6 A2). Both are in the
output; **the simulation is what is reported**, and it is the stricter of the two at every
informative cut.

| cut | **oracle in a majority cluster** | chance | **p** | incumbent's pick | chance |
|---|---|---|---|---|---|
| Chamfer 0.50 Å | **0.299** | 0.423 | **0.0025** | 0.425 | 0.423 |
| Chamfer 0.75 Å | **0.448** | 0.606 | **0.0000** | 0.563 | 0.606 |
| **Chamfer 1.00 Å** | **0.644** | 0.773 | **0.0008** | **0.770** | 0.773 |
| Chamfer 1.25 Å | 0.816 | 0.884 | 0.018 | 0.897 | 0.884 |
| **Rotation 15°** | **0.586** | 0.735 | **0.0003** | **0.736** | 0.735 |
| Rotation 30° | 0.747 | 0.835 | 0.010 | 0.862 | 0.835 |

**The dissociation is the result.** The **oracle** pose is in a majority cluster
*significantly below* chance at every cut that makes a real partition. The **incumbent's**
pick is at or *above* chance at **14 of 14** cuts. The two statistics are computed from the
same clustering, on the same ligands, in the same run.

P3 agrees and is signed the same way at **14 of 14** cuts — the oracle pose's cluster is
always smaller than the incumbent's:

| cut | size(oracle) | size(incumbent) | Δ | Wilcoxon p |
|---|---|---|---|---|
| Chamfer 0.50 Å | 2.98 | 3.64 | **−0.67** | 0.034 |
| Chamfer 0.75 Å | 7.55 | 9.08 | **−1.53** | 0.020 |
| Chamfer 1.00 Å | 11.56 | 13.02 | −1.46 | 0.057 |
| Rotation 15° | 10.68 | 12.54 | **−1.86** | **0.005** |
| Rotation 30° | 13.69 | 15.31 | **−1.62** | 0.030 |
| Rotation 120° | 17.62 | 18.69 | **−1.07** | 0.018 |

### 3c. P4 — conditioning on the stake **doubles** the effect

Required by the pre-registration, because a ligand whose poses are all equivalent must not
dilute one that matters. On the 45 HIGH-stake ligands, **minority rate** = 1 − (in a
majority cluster):

| cut | **minority rate, observed** | chance | p |
|---|---|---|---|
| Chamfer 0.50 Å | **66.7%** | 52.1% | 0.0055 |
| Chamfer 0.75 Å | **64.4%** | 41.8% | 0.0010 |
| **Chamfer 1.00 Å** | **46.7%** | **23.1%** | **0.0000** |
| Chamfer 1.25 Å | **31.1%** | 15.0% | 0.0003 |
| **Rotation 15°** | **53.3%** | **29.4%** | **0.0003** |
| Rotation 30° | 35.6% | 20.3% | 0.0043 |

**Where the score is, the best pose is in a minority orientation cluster about twice as
often as chance.** Six cuts, two independent distance definitions, every one in the
predicted direction, four of six at p < 0.005. The cluster-size gap widens with it
(Chamfer 1.00 Å: 10.22 vs 12.76 on HIGH stake, against 11.56 vs 13.02 pooled).

**P1, P2, P3, P4 — all four pre-registered predictions hold.** The minority-report account
of `FINDING_036` §3 is confirmed, and this is the first direct measurement of the mechanism
that 024 and 036 inferred.

---

## 4. Step 3 — when is the consensus confidently wrong? Rarely, and it is the wrong failure

Prediction-side only. TIGHT = bottom tercile of mean pairwise pose-to-pose Chamfer
(≤ 0.768 Å); WRONG = incumbent-selected LDDT-PLI < 0.50. Both fixed before looking.

| | WRONG | RIGHT | rate wrong |
|---|---|---|---|
| **TIGHT** (n = 29) | **6** | 23 | **20.7%** |
| **WIDE** (n = 58) | 17 | 41 | **29.3%** |

`FINDING_012`'s catastrophe-detector account predicts the inverse of the naive worry, and
that is what happens: **a tight pool is *less* often wrong, not more.** Confidently wrong is
**6 of 87 ligands (6.9%)** — A1A4T, A1BNX, CFF, CL6, NJ0, T0K.

And those six are not a rescue opportunity. Their **oracle is 0.526** against a selected
0.364: when the pool agrees and is wrong, *the whole pool is wrong*, so there is no minority
cluster to promote. That is the same object `FINDING_024` called "dispersed inside a wrong
basin", here at its tight extreme.

Nothing prediction-side separates them. Seven covariates fixed in advance —
`spread`, `n_clusters`, `margin_12`, mean `xeng`, `n_heavy`, `pred_fe_donor_median`,
`pred_frac_coordinated` — Mann-Whitney TIGHT∧WRONG (n = 6) vs TIGHT∧RIGHT (n = 23):
**every raw p ≥ 0.23, every Holm-corrected p = 1.00, 0 of 7 survive.** At n = 6 the test is
underpowered and is reported as such; it is not evidence of absence, it is an absence of
evidence, and it is why no selector is proposed from it.

---

## 5. Step 4 — the term oracle. The number that decides everything

The perfect rule gets to choose the **cluster** knowing the truth, then must take the pose a
**prediction-side** rule names inside it (`argmin xeng` within the cluster). It does **not**
get to choose the pose — that would be the pool oracle and would test nothing.

**Two bars, both recomputed here.** **F87**, the `FINDING_007` random-feature floor over
4,000 draws: **p95 +0.0134 / p99 +0.0193** (007 measured +0.0137 / +0.0193 — reproduced).
**N-A**, a matched null that shuffles cluster labels within each ligand **preserving the
size multiset** and recomputes O-A identically, 4,000 draws — because **O-A is a
best-of-`k` oracle and rises with `k` by construction**.

| cut | mean clusters | **O-A** | 95% CI | switched | share of the 0.0811 gap | N-A mean | **excess over N-A** |
|---|---|---|---|---|---|---|---|
| Chamfer 0.50 Å | 12.71 | **+0.0745** | [+0.057, +0.094] | 64 | 91.9% | 0.0631 | +0.0114 |
| Chamfer 0.75 Å | 6.59 | +0.0595 | [+0.042, +0.079] | 50 | 73.4% | 0.0437 | **+0.0158** |
| **Chamfer 1.00 Å** | **3.55** | **+0.0427** | [+0.028, +0.059] | 35 | 52.7% | 0.0292 | **+0.0135** |
| Chamfer 1.25 Å | 2.10 | +0.0241 | [+0.013, +0.037] | 20 | 29.7% | 0.0170 | +0.0071 |
| Chamfer 1.50 Å | 1.41 | +0.0040 | [0.000, +0.010] | 3 | 4.9% | 0.0081 | −0.0041 |
| Chamfer 2.00 Å | 1.07 | +0.0032 | [0.000, +0.009] | 2 | 3.9% | 0.0021 | +0.0011 |
| Chamfer 2.50–3.00 Å | 1.03 | **0.0000** | — | 0 | 0% | 0.0014 | −0.0014 |
| Rotation 15° | 4.13 | +0.0537 | [+0.036, +0.074] | 38 | 66.2% | 0.0346 | **+0.0191** |
| Rotation 30° | 2.61 | +0.0391 | [+0.022, +0.058] | 23 | 48.2% | 0.0224 | +0.0167 |
| Rotation 45° | 2.08 | +0.0328 | [+0.017, +0.052] | 16 | 40.4% | 0.0167 | +0.0161 |
| Rotation 90° | 1.60 | +0.0225 | [+0.008, +0.040] | 9 | 27.7% | 0.0101 | +0.0124 |
| Rotation 120° | 1.38 | +0.0186 | [+0.005, +0.035] | 7 | 22.9% | 0.0058 | +0.0128 |

### 5a. The answer depends entirely on the cut, and that dependence is a degeneracy

`O-A` ranges from **0.0000 to +0.0745** across a grid nobody is allowed to tune, and it is
**monotone in the number of clusters**. The reason is not subtle: as the cut tightens toward
singletons, `O-A` converges on the **full pose oracle, +0.0811** — at which point *"take the
minority cluster"* means *"take any pose you like"* and carries **no information at all**.
At the loose end every pose is one cluster and the rule cannot fire. **The pre-registration
required this to be stated plainly, and it is: a single headline number for this ceiling
does not exist.**

The honest quantity is the **excess over the size-matched shuffle** — the part attributable
to clustering by **orientation** rather than into same-shaped groups:

| | median over grid | max over grid | **F87 p95** |
|---|---|---|---|
| Chamfer | **+0.0041** | +0.0158 | **+0.0134** |
| Rotation | **+0.0150** | +0.0191 | **+0.0134** |

**Straddling the floor.** The best single number anywhere on the grid, +0.0191 at
rotation 15°, is **1.4× the floor** and **below the p99 of +0.0193** — for a rule that is
handed the answer.

### 5b. And the majority vote pays more than the minority rescue is worth

The decisive comparison, declared post-hoc (§6, A4), asks the inverse question on the same
clustering: is the incumbent's pick better than the representative of a **randomly chosen**
orientation mode?

| cut | mean clusters | **incumbent** | random cluster | worst cluster | **incumbent − random** | p | perfect rescue (O-A) |
|---|---|---|---|---|---|---|---|
| Chamfer 0.75 Å | 6.59 | 0.6164 | 0.5642 | 0.4497 | **+0.0522** | 2.9e-06 | +0.0595 |
| **Chamfer 1.00 Å** | 3.55 | **0.6164** | **0.5617** | **0.4655** | **+0.0546** | **1.3e-06** | **+0.0427** |
| Chamfer 1.25 Å | 2.10 | 0.6164 | 0.5821 | 0.5279 | +0.0343 | 9.2e-04 | +0.0241 |
| **Rotation 15°** | 4.13 | 0.6164 | 0.5619 | 0.4573 | **+0.0545** | 1.0e-06 | +0.0537 |
| Rotation 30° | 2.61 | 0.6164 | 0.5644 | 0.4805 | **+0.0520** | 7.6e-06 | +0.0391 |
| Rotation 45° | 2.08 | 0.6164 | 0.5675 | 0.4925 | +0.0489 | 1.8e-05 | +0.0328 |

At Chamfer 1.00 Å — the cut with the strongest mechanism signal in §3 — **choosing the
majority mode is worth +0.0546 and perfectly rescuing the minority mode is worth +0.0427.**
The downside of getting it wrong is **−0.151** (worst cluster 0.4655 against 0.6164). Five
of six rows are the same shape.

This is `FINDING_036` §2c one level up, and it is the same arithmetic. There,
`argmin(xeng)` beat a coin flip over its own top 2 by **+0.0161**, *more than the entire
remaining top-2 prize of +0.0129*. Here, `argmin(xeng)` beats a random orientation mode by
**+0.0546**, *more than the entire perfect minority rescue of +0.0427*, with an asymmetry
three times larger.

**So both halves of the hypothesis are true and they do not add up to a prize.** The truth
*is* a minority report. `xeng` *is* majority-voting. And the majority vote is still the
right bet, because on the ligands where the consensus is right — which is most of them —
the minority clusters are catastrophically worse. Without an oracle, the expected value of
"switch to the minority" is strongly negative.

---

## 6. Deviations and post-hoc work, declared

| id | what | why it is not a pre-registered endpoint |
|---|---|---|
| **M1** | E2a's chance rate was pre-registered as `mean(size_largest / n)`; it undercounts when several clusters tie for largest. Replaced by a 4,000-draw simulation | a **methods correction**, reported alongside the pre-registered version; the simulation is stricter at every informative cut, so it cannot flatter the hypothesis |
| **A1** | ρ(rank, stake) recomputed on the 73 ligands with stake > 0 | `rank_oracle == 1` *forces* `stake == 0`, so the pooled ρ = +0.601 is partly tautological. De-confounded: **+0.327, p = 0.0047** |
| **A2** | the E2a simulation null, pooled and on HIGH stake | see M1 |
| **A3** | O-A minus the matched-shuffle mean, per cut | the pre-registration named the N-A p95 bar; the *mean* excess is the readable version of the same quantity |
| **A4** | incumbent vs a random cluster representative | **not pre-registered.** It asks the inverse question and it is what turns a bounded ceiling into a verdict. It is reported because it fired against the hypothesis, not for it |

No fitted parameters. No tuned cut. No conclusion drawn from n = 14 — pool B was not used.
No answer-recognition gate was run on any comparator, per `FINDING_036` §3.

---

## 7. Verdict — **MEASURED**; mechanism **CONFIRMED**, rescue **REFUTED** in substance

Against the pre-registered table:

| clause | required | measured | fires? |
|---|---|---|---|
| **CLOSED / REFUTED** | **max** over the Chamfer grid < F87 p95 | **+0.0745** vs +0.0134 | no |
| **LICENSED** | median-over-grid ≥ 2 × F87 **and** beats N-A p95 at ≥ half the cuts | median **+0.0140** vs the required +0.0268; 4 of 8 cuts beat N-A | **no** (fails the median clause) |
| **MEASURED** | anything else | | **YES** |

**MEASURED** is the pre-registered verdict, and no selector was built, as the
pre-registration required under every outcome.

### What it establishes

1. **The oracle pose's median rank is 8** — 036's rhetorical question has a literal answer —
   **and the ranking inverts where the money is.** LOW stake: median rank 3.5, 33% at rank 1.
   HIGH stake: median rank **11**, **2 of 45** in the top 3, **7 of 45** in the top 5, mean
   rank **12.04 against a lottery's 10.5**. ρ(rank, stake) = **+0.601** pooled and **+0.327**
   (p = 0.0047) once the rank-1 tautology is removed. The incumbent solves the ligands that
   are already solved.
2. **The minority-report mechanism is confirmed, directly, on all four pre-registered
   predictions.** `xeng` ranks by orientation-cluster membership (ρ = −0.25 to −0.34,
   correct on 69–85% of ligands, p to 2.6e-07, at **every** cut of **both** metrics); the
   oracle pose is in a majority cluster **significantly below chance** while the incumbent's
   pick is at or **above** chance (**14 of 14** cuts); the oracle's cluster is smaller at
   **14 of 14** cuts; and **conditioning on the stake doubles it** — a **46.7%** minority
   rate against a **23.1%** chance rate at Chamfer 1.00 Å. This is the first direct
   measurement of what `FINDING_024` and `FINDING_036` §3 inferred.
3. **And it is not a prize.** The perfect rescue is worth **+0.0427** at a cut that makes a
   real partition, **but the ceiling is a monotone function of the cut** and converges on the
   full pose oracle as clusters become singletons, where the rule means nothing. The
   orientation-specific excess over a size-matched shuffle is **+0.0041 median / +0.0158 max**
   (Chamfer) and **+0.0150 median / +0.0191 max** (rotation) — **straddling the +0.0134
   floor, with the single best number below p99.**
4. **The majority vote pays more than the rescue is worth.** Choosing the majority mode
   beats choosing a random one by **+0.0546** (p = 1.3e-06) where the perfect rescue is worth
   **+0.0427**, and the worst mode costs **−0.151**. `FINDING_036` §2c, one level up, with a
   three-times-larger asymmetry. **A biased-but-informative consensus is worth defending.**
5. **"Confidently wrong" is rare and is a different failure.** TIGHT∧WRONG is **6 of 87**,
   and a tight pool is **less** often wrong (20.7% vs 29.3%) — `FINDING_012`'s account
   holding. Those six have an **oracle of 0.526**: when the pool agrees and is wrong, the
   whole pool is wrong and there is no minority to promote. **0 of 7** prediction-side
   covariates separate them (all Holm p = 1.00), at n = 6, which is underpowered and said so.

### What it forecloses

- **Promoting a minority orientation cluster, by any rule, on this pool.** Closed by the
  term oracle once the best-of-`k` component is subtracted, and closed again by A4 — the
  rule it would have to beat is the one it was proposed to replace.
- **Reading `FINDING_036` §3's mechanism as an opportunity.** It was the right diagnosis and
  it is now measured; it is an explanation of why `xeng` is biased, not a route to fixing it.
  A displaced-but-informative estimator is the *working* case, not the broken one.
- **Cluster-structure statistics of the pool as a source of selection value.** Cluster size,
  membership and count are all functions of the same pose-to-pose geometry whose *mean*
  `xeng` already reads. `FINDING_036` retired the decomposition of `D[pose, ref]`; this
  retires the decomposition of `D[pose, pose]`.
- **Using "the consensus is tight" as a warning flag.** It points the wrong way (20.7% vs
  29.3%), which is `FINDING_012` making a correct prospective prediction for the fourth time.

### What it does NOT license

- It does not say the clustering is uninformative. The orientation-specific excess is
  **positive at 11 of 14 cuts** and reaches **+0.0191**; it is *at* the floor, not below it.
  A pool with genuinely bimodal orientations — the `FINDING_012` catastrophe regime, which
  the blind cryoEM release may be, and which OpenADMET's "multiple mutually exclusive
  conformations" explicitly describes — would have more clusters and a larger stake. **There
  is no drop-day proxy for it and none is claimed.**
- It does not explain **why** the oracle is ranked 11th where the stake is. It measures that
  it is, and that the pool's own orientation modes cannot be used to find it.
- The step-3 separator test is **n = 6** and proves nothing about those ligands.
- One pool, one reference set, one protein.

### What to do instead

Three findings have now converged on the same wall from different sides. `FINDING_036`: the
top-2 prize does not exist. `FINDING_032`: the torsion oracle is +0.0213 and loses to the
incumbent. Here: the best pose is a genuine minority report, and rescuing it perfectly is
worth less than defending the majority. **Every route that reads the pool's own geometry —
its mean (011), its higher moments (036), its modal structure (037) — is now measured and
spent.** The remaining **0.081** is not hiding in a statistic of the 20 poses; §2b says it
is concentrated on 45 ligands whose best pose the pool ranks 11th on average, and §3 says
the pool's own agreement structure actively points away from it. Anything that moves it must
bring information the pool does not contain.

---

## 8. Reproduce

```
python scripts/structure/minority_report.py controls   # C-XENG C-SEL C-N2 C-NUM C-REF C-ATOM C-FILT
python scripts/structure/minority_report.py rank       # STEP 1 -- where the oracle pose sits
python scripts/structure/minority_report.py cluster    # STEPS 2 + 4 -- minority + the term oracle
python scripts/structure/minority_report.py confident  # STEP 3 -- confidently wrong
python scripts/structure/minority_report.py addendum   # A1 A2 A3 A4, all post-hoc
```

| what | where |
|---|---|
| pre-registration | `docs/PREREG_minority_report.md` (`6f547ef`) |
| script, five stages | `scripts/structure/minority_report.py` |
| controls | `data/processed/minority_controls.json` |
| **step 1, the ranking** | `data/processed/minority_rank.json`, `minority_per_ligand.csv` |
| **steps 2 + 4, every cut** | `data/processed/minority_clusters.json`, `minority_cluster_per_ligand.csv` |
| step 3, confidently wrong | `data/processed/minority_confident.json`, `minority_confident_per_ligand.csv` |
| post-hoc addendum | `data/processed/minority_addendum.json` |

Heme-frame and pose–reference caches are reused from `C:\Temp\cyp_tiebreak`; the
pose-to-pose distance matrices are this experiment's own intermediate and live in
`C:\Temp\cyp_minority` (224 KB), never on `D:`. No inference was run and nothing was
downloaded.
