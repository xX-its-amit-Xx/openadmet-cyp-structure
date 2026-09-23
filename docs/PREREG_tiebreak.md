# PRE-REGISTRATION — breaking the near-tie

**Date:** 2026-09-23 · **Committed before any endpoint was computed.** Nothing below may be
edited afterwards; anything added later is labelled post-hoc in the finding.

Successor to `docs/FINDING_035_matched_depth.md`, which named this experiment in its
item 5: *"The near-tie is the target, not the pose."*

---

## Where this comes from

`FINDING_035` doubled the pool at genuinely matched conditioning (20 → 40 on the 14
predicted-Type-I ligands) and found Δ selected **+0.0116**, Wilcoxon **p = 0.715**, with
**10 of 14 ligands unchanged**. It identified the mechanism precisely: depth does not hand
the selector a better pose, it hands it more **near-ties**. Both swing ligands turned on a
**sub-0.06 Å Chamfer margin** that produced a **> 0.3 LDDT-PLI swing** — CFF **+0.4826**
on a 0.057 Å margin, 08J **−0.3442** on a 0.045 Å margin.

So a large amount of score is being decided by margins far below the meaningful resolution
of the statistic deciding them. That makes tie-breaking the best-evidenced remaining lever
— *if* there is anything systematically recoverable there, which is exactly what stage 1
is for and why it runs first.

---

## Data — everything already on disk, zero new inference, CPU only

| | |
|---|---|
| **pool A (primary)** | `val87b_unsteered` — 87 CYP3A4 ligands × 20 unsteered Boltz-2 poses = **1,740**, mmCIFs at `D:/cyp_scratch/val87b_unsteered/` |
| truth | `data/processed/poses_scored_val87b.csv`, `arm == "unsteered"` |
| shipped feature | `data/processed/xeng_val87b.csv` |
| **pool B (secondary)** | the `FINDING_035` stratum — **14** predicted-Type-I ligands at depth **40** (20 Modal + 20 Explorer), `data/processed/matched_depth_poses.csv` plus pool A's rows for those ligands |
| reference set | `data/processed/reference_set_cyp3a4.npz` — frozen, two Protenix checkpoints, deduplicated, depth 6–11 |
| crystals | `data/reference/rcsb/<PDB>.cif`, one chain, chosen by `score_pool.load_reference`'s rule |
| selector under test | `cypstruct.xengine.select()` — `argmin(xeng)`, no fitted parameters |

---

## Controls — all must pass and be reported before any endpoint is read

| id | control | rule |
|---|---|---|
| **C-XENG** | recompute `xeng` from the pool mmCIFs + the frozen reference | max abs diff vs the shipped column **< 1e-6** on all 1,740 rows |
| **C-SEL** | `argmin(xeng)` **is** `cypstruct.xengine.select()` | equal to **1e-9** on all 87 (`FINDING_035` N2b — re-verified, not assumed) |
| **C-N2** | the shipped board reproduces | selected **0.6164**, oracle **0.6975**, random **0.5769**, gain **+0.0395**, ρ **−0.2582**, correct sign **75.86%**, each to 5e-4 |
| **C-NUM** | `FINDING_021` numbering | `mapped` True on every scored row; every exact-zero LDDT-PLI row must be an **ejection** (BiSyRMSD > 10 Å), not a numbering failure |
| **C-REF** | reference depth | ≥ **4** independent poses on every ligand (`FINDING_011`'s condition) |
| **C-XTAL** | crystal load, for the R1 gate | ligand and heme present in one chain, heme frame constructible; count of ligands that fail is reported, not absorbed |
| **C-FILT** | every filter fires | each stage reports a count **including the zeros** |

---

## Stage 1 — characterise the near-tie regime. This runs FIRST and can stop the experiment.

Per ligand, poses sorted **ascending** by `xeng` (lower is better; `x(1) ≤ x(2) ≤ …`):

| quantity | definition |
|---|---|
| `margin_12` | `x(2) − x(1)`, in Å of Chamfer |
| `zmargin_12` | `margin_12 / sd(xeng)` within that ligand |
| `stake_k` | `max(lddt_pli over the k lowest-xeng poses) − lddt_pli(argmin xeng)`, **≥ 0** by construction, for **k = 2, 3, 5** |
| **oracle tie-break at k** | `mean(stake_k)` over ligands — **the entire budget for this line of work** |
| `downside_k` | `lddt_pli(argmin) − min(lddt_pli over the k lowest-xeng)` — what a *worst-case* tie-break costs |
| top-1 vs random-in-top-k | selected(argmin) − mean(lddt_pli over the top k) — does the `xeng` ordering carry information at sub-margin resolution at all? |

Reported: the full distribution of `margin_12` and `zmargin_12`; mean/median `stake_k`;
mean `stake_2` within **margin quartiles**, to show whether the stake concentrates in the
near-ties as `FINDING_035` implies it should.

### The stop rule, fixed now

Let `O2 = mean(stake_2)` on pool A (n = 87).

* **O2 < +0.0137** (the `FINDING_007` pooled p95 random-feature floor on this pool) →
  the line is **CLOSED**. No tie-breaker is built. Verdict **MEASURED**; the finding
  reports the ceiling and says so plainly.
* **+0.0137 ≤ O2 < +0.0274** → candidates are tested, but the expectation recorded *here
  in advance* is that the paired bar will not be cleared, because a rule capturing even
  half of a ceiling that small cannot beat the floor.
* **O2 ≥ +0.0274** → a genuine prospect; proceed.

---

## Stage 2 — the two diagnostic rungs, BEFORE any selector is built

Both are leaky by construction and **never enter a score**.

### R1 — the answer-recognition gate (`FINDING_030`, sharpened by `032`)

Score the query's own **crystal ligand** on each candidate feature — its heavy atoms in its
own heme frame, compared to the same reference poses. Chamfer is permutation-invariant, so
no atom-order mapping is required and the `FINDING_011` ordering trap cannot apply. Record
the crystal's **percentile** among that ligand's 20 predicted poses (1.0 = the crystal is
the single best-scoring object in the set).

**Pass rule, fixed a priori: mean percentile ≥ 0.65 and binomial p < 0.05 in the correct
direction.**

**The calibrator, fixed a priori.** The shipped `xeng` is scored on the identical gate as a
positive control, because every candidate here is a function of the same pose–reference
distance matrix that `xeng` is.

* If **`xeng` passes R1**, the gate is informative for this family and any candidate that
  fails it is **dead** and is not carried to the paired bar.
* If **`xeng` fails R1**, the gate is **not valid for this family** — a feature that
  demonstrably selects cannot be killed by a gate it also fails — and R1 is reported as a
  diagnostic only, with every candidate still carried forward. This branch is declared now
  so that the outcome cannot be read either way after the fact.

### R2 — the term oracle (`FINDING_032`)

For a tie-break at depth `k`, the term oracle **is** `mean(stake_k)` from stage 1: it is
exactly what perfect knowledge of the tie-break is worth. It is reported before any
candidate number appears, and it bounds every prior of this shape.

---

## Stage 3 — the candidates. Signs fixed here, from geometry, never per fold.

All are functions of

* `D[p, r]` — Chamfer between pool pose `p` (own heme frame) and reference pose `r` (own
  heme frame), `m` ≥ 4 deduplicated references, the **same matrix `xeng` averages**, and
* `C[p, p']` — Chamfer between two pool poses of the same ligand.

**No new information. No new inference. No fitted parameters.** Every candidate operates
**only inside the top-`k` by `xeng`**: the incumbent's ordering is kept and only the top `k`
is reordered. `k ∈ {2, 3, 5}`. Ties *within* a candidate are broken at random over **64
draws**.

| id | rule | statistic | **sign, fixed a priori** | why that sign |
|---|---|---|---|---|
| **T1** | **medoid of the top-k**: mean Chamfer from this pose to the **other top-k poses** | `med_k` | **lower is better** | `FINDING_003` — a consensus medoid was the first selector that ever worked here. Among near-ties, the pose the other near-ties cluster around is the modal opinion. **Undefined at k = 2** (the mutual distance is symmetric); reported N/A, not computed |
| **T2** | distance to the **2nd-nearest** reference, not the mean | `xeng_2nd` | **lower is better** | same direction as `xeng`. The min is the noisiest order statistic and the mean is dragged by a reference that is far from everything; the 2nd is the robust middle |
| **T3** | **Borda rank**: rank all 20 poses by `D[·, r]` for each reference `r`, take the mean rank | `xeng_borda` | **lower is better** | removes per-reference scale, so one badly placed reference cannot dominate the average distance |
| **T4** | **within-pose sd** of `D[p, ·]` across references | `xeng_sd` | **lower is better** | a pose that *all* references agree is close to is a firmer consensus than one a single reference likes; `FINDING_011`'s "quality of the opinions" argument applied within a pose |
| **T5** | **plurality vote**: how many references have this pose as their single nearest among the ligand's 20 | `xeng_votes` | **higher is better** | the argmin of a mean can be nobody's favourite. This is the pose the largest number of independent opinions individually prefer |

Each candidate is additionally reported as a **standalone argmin over all 20** as a
descriptive control — that is not a tie-break and is not eligible to ship.

---

## Stage 4 — the bars

**BAR 1 — versus random.** Gain = selected − the per-ligand mean LDDT-PLI (the exact
expectation of a random pick, not a simulated one). Two nulls, both recomputed **for these
features on this pool**:

1. the **matched tie-break null** — 4,000 draws of a *random* reorder of the same top-`k`,
   which has exactly the degrees of freedom a tie-break has;
2. the `FINDING_007` **random-feature null** — 4,000 draws of a random within-ligand
   feature put through the identical pipeline.

Quote p95 and p99. **The floor is n-dependent** (`FINDING_035`: p95 **+0.0453** inside
n = 14 against **+0.0137** pooled), so each set is quoted against its own floor.

**BAR 2 — paired against the shipped incumbent, on identical poses.** This is the only
comparison that decides anything; nine successive candidates have cleared bar 1 and died
here.

* per-ligand difference vs `cypstruct.xengine.select()`, mean and **10,000-sample bootstrap
  95% CI**;
* **Wilcoxon signed-rank p**;
* the **TIE COUNT** — ligands where the pick is unchanged;
* **how many picks the rule CHANGES**, and the outcome on **exactly those ligands**,
  reported separately from the pooled mean. A rule that moves 2 of 87 picks cannot be
  judged by a pooled average.

**BAR 3 — ranking, beside every pooled mean.** Within-ligand Spearman ρ against LDDT-PLI
and the correct-sign fraction, for every candidate. `FINDING_011` records that ρ and top-1
selection have moved in **opposite** directions twice, so ρ is reported and is never the
criterion.

**The pool oracle is reported before every selection number.**

---

## Acceptance rule

| verdict | requires |
|---|---|
| **SHIPS** | all four: (i) paired mean vs the incumbent ≥ **+0.0137** on pool A; (ii) the 95% bootstrap CI of the paired difference **excludes 0**; (iii) Wilcoxon **p < 0.05**; (iv) the sign is the one fixed above, not one chosen after seeing the data |
| **REFUTED** | the best candidate's paired mean is **≤ 0** and its CI upper bound is below **+0.0137** |
| **MEASURED** | anything else — **including** the stage-1 stop, where the oracle ceiling is below the floor and no candidate is built |

A candidate that clears bar 1 and fails bar 2 is a **negative** and is logged as one. A
candidate selected as the best of five after the fact is labelled post-hoc and is not
eligible to ship, per `FINDING_011`'s retraction of `-zm - zx`.

---

## What I will conclude if the oracle ceiling is small

Stated now, so the conclusion is not written to fit the number.

If a perfect top-2 tie-break is worth less than **+0.0137**, then **no tie-breaker of any
shape can clear the noise floor on this pool**, and this line closes the way `FINDING_032`
closed the internal conformer — by the **oracle**, not by the gate. The correct reading of
`FINDING_035`'s two swing ligands would then be that the near-tie is a **variance
mechanism** — it explains why depth does not convert, because it adds coin flips with large
payoffs in both directions — and **not a prize**, because the flips are symmetric and
nothing systematic is recoverable from them. The near-ties would be recorded as an
irreducible coin flip, and the recommendation would be to spend the remaining runway on the
**0.081 LDDT-PLI of oracle gap that lies OUTSIDE the top-k**, not on the sliver inside it.

Equally, if the ceiling is large and every candidate still fails, that is the more
interesting negative: it would say the pose–reference distance matrix has been fully
exploited by the mean, and that the information needed to break a tie is **not in the
reference set at all**.

---

## Traps, explicitly handled

1. **Reproduce before trusting** — C-XENG and C-SEL run first; `FINDING_035` verified
   `argmin(xeng) == select()` to 0.0 and it is **re-verified here**, not assumed.
2. **Numbering** (`FINDING_021`) — C-NUM, with the exact-zero rows classified as ejections
   or flagged.
3. **A rule that changes few picks** — the changed-pick count and the outcome on exactly
   those ligands are reported separately from the pooled mean.
4. **Signs** — fixed above from geometry. `loo-sign-selection-fakes-negatives`: picking a
   sign per fold manufactures p = 9e-09 from noise.
5. **Random tie-breaking** — ≥ 64 draws wherever a tie is broken at random, 4,000 for nulls.
6. **Too-clean numbers** (`too-clean-numbers-are-the-tell`) — a filter that never fires or
   an sd of exactly 0 is reported as a failure, not as a pass.
7. **Disk** — intermediates to `C:\Temp` or the session scratchpad. Never `D:`.

## Artefacts this will produce

| what | where |
|---|---|
| script | `scripts/structure/tiebreak_analysis.py` |
| controls | `data/processed/tiebreak_controls.json` |
| stage 1, the prize | `data/processed/tiebreak_regime.json` |
| R1 answer-recognition | `data/processed/tiebreak_answer_recognition.json` |
| candidates + bars | `data/processed/tiebreak_selectors.json` |
| per-ligand table | `data/processed/tiebreak_per_ligand.csv` |
| finding | `docs/FINDING_036_tiebreak.md` |
