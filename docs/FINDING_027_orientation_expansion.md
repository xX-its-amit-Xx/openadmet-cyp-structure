# FINDING 027 — the missing poses are not one rotation away, because the model's own pocket does not admit them

**Date:** 2026-09-22 · **Status:** measured, zero new inference, CPU-only ·
**Verdict:** **MEASURED** — and by 0.0008. The oracle rises +0.0108 against a
pre-registered bar of +0.0100, selection **falls** −0.0145, and the mechanism behind both
numbers refutes the lever outright.

Pre-registered in `docs/PREREG_orientation_expansion.md`, committed at `c23fb76` before a
single candidate existed. Nothing below was tuned afterwards.

---

## What was asked

FINDING 024 decomposed CYP3A4's failure: the pocket is right (pocket-lining CA **0.73 Å**),
the ligand centroid is right (**1.07 Å**, isotropic, no direction), and only the
**orientation** is wrong (**30.1°** against **6.1°** for the rest of the family), with
ρ(protein error, ligand error) = **+0.030, p = 0.80** inside CYP3A4 against **+0.510** 
everywhere else. Every protein-conditioning remedy is refuted by that mechanism, and
FINDING 026 closed the selector-side proposal that followed from it.

One lever was left standing, and it is the one this tests:

> If the receptor is correct and only the rotation is wrong, the missing poses are
> reachable by searching **orientation** space inside the model's **own** protein, with no
> new network inference at all.

Test that ceiling honestly. **44 of 87 pairs** (41 in this pool) never produce a sub-2 Å
pose in 20 diffusion samples.

## What was done

`scripts/structure/orientation_expansion.py`, three CPU stages, ~2.5 minutes on 12 cores.

| | |
|---|---|
| pool | `val87b_unsteered` — 87 ligands × 20 Boltz-2 samples = **1,740 poses** |
| grid | **1,024** super-Fibonacci orientations of SO(3) (Alexa, CVPR 2022), deterministic, identical for every pose. Realised nearest-neighbour geodesic angle **22.43° mean / 22.53° median / 23.96° max** |
| transform | `R·(L − c) + c` about the ligand's **own centroid**, which is held fixed because 024 measured it as already correct. Protein, heme and internal conformer untouched |
| candidates | 1,740 × 1,024 = **1,781,760** |
| filters (prediction-only) | **(a)** no ligand heavy atom within **2.2 Å** of any protein or heme heavy atom, Fe excluded; **(b)** if the unrotated pose coordinates the iron (`is_coordinated`, Fe–donor ∈ [1.85, 2.55] Å and S–Fe–donor ≥ 150°), the candidate must too |
| survivors | **23,688 — a pass rate of 1.33%** (clash alone: 2.06%) |
| selector | `cypstruct.xengine.select()`, unchanged, on `xeng` = mean symmetric Chamfer in the pose's own heme frame to the frozen independent-engine references |

84.2% of the original poses coordinate the iron, so filter (b) binds on most of the pool.

### The four controls, all passed before any result was read

| control | required | measured |
|---|---|---|
| **C1** recomputed LDDT-PLI = shipped `poses_scored_val87b.csv` | \|Δ\| < 1e-9 on ≥ 99% of 1,740 rows | **max \|Δ\| = 1.1e-16, 100% of rows** |
| **C2** FINDING 021 numbering | no exact 0.0 except a genuine ejection | offset **0 for all 1,740**, min residue-name identity **0.996**; the 4 exact zeros all have BiSyRMSD > 10 Å |
| **C3** recomputed `xeng` = shipped `xeng_val87b.csv` | \|Δ\| < 1e-6 | **max \|Δ\| = 2.8e-8** |
| **C4** atom counts / mapping | agree on every row | **1,740 / 1,740**, all symmetry-mapped |

C1 and C3 together say the thing under test genuinely *is* the shipped metric and the
shipped selector, not a lookalike.

---

## The headline

Random-pose baseline **0.5769**. Ties broken over 64 draws; the selector has no fitted
parameters, so leave-one-ligand-out is **vacuous here** and is said to be vacuous rather
than performed decoratively — no cross-ligand quantity is fitted anywhere in this work.

| pool | poses per ligand | **oracle** | **selected** |
|---|---|---|---|
| original | 20 | **0.6975** | **0.6164** |
| **+ filtered rotations** | 20 + 272 | **0.7083** (+0.0108) | **0.6019** (−0.0145) |
| + size-matched unfiltered null | 20 + 272 | 0.7028 (+0.0053) | 0.6132 (−0.0032) |
| whole 1,024-grid, no filter at all | 20 + 20,480 | 0.7277 (+0.0302) | — |

**The oracle rises and selection falls.** That is CLAIM F of `CLAUDE.md` reproduced on a
pool 14× larger, and FINDING 013 again in a new costume.

### The paired tests, all three, with Holm over all three

| comparison | mean Δ | bootstrap 95% CI | Wilcoxon | Holm | better/worse/tied |
|---|---|---|---|---|---|
| filtered expansion vs original | **−0.0145** | **[−0.0314, −0.0020]** | 0.080 | 0.239 | 1 / 4 / **82** |
| unfiltered null vs original | −0.0032 | [−0.0103, +0.0006] | 0.655 | 0.655 | 1 / 1 / **85** |
| filtered vs unfiltered null | −0.0113 | [−0.0273, +0.0007] | 0.144 | 0.288 | 1 / 3 / **83** |

**Two pre-registered tests disagree on the same data and the disagreement is instructive.**
The bootstrap CI on the top row excludes zero; the Wilcoxon does not. The tie column is
why: 82 of 87 ligands select the *identical* pose before and after, so the Wilcoxon sees
**five** informative pairs and its smallest attainable two-sided p at n = 5 is 0.0625. The
bootstrap resamples the full 87-vector including its 82 zeros and is driven by a handful of
large negatives. A bootstrap CI on a mostly-zero paired vector is not a hypothesis test
with the same operating characteristics, and where they disagree the Wilcoxon is the one
this project's gates are written against. The honest reading is *the drop is real in
direction and unproven in significance*, which is the same shape as FINDING 026's +0.0041.

### The rescue statistic

| | n |
|---|---|
| pairs reaching sub-2 Å BiSyRMSD, original pool | **46 of 87** |
| ... with the filtered expansion | **52 of 87** (**6 rescued**) |
| ... with the size-matched unfiltered null | 47 of 87 |
| ... with the whole unfiltered 1,024-grid | **56 of 87** (10 of the 41 failures) |

So the expansion does rescue pairs — and **31 of the 41 failures cannot be rescued by any
rotation on the grid, filtered or not.** Median minimum BiSyRMSD over the 41 failures:
3.73 Å original → 3.25 Å filtered → **2.86 Å with the entire grid and no filter at all**.

The grid is not the reason. At 22.43° mean nearest-neighbour spacing the worst-case
residual is 11.2°, which for the median radius of gyration of 4.78 Å is a ligand RMSD of
**0.93 Å** — subtracted in quadrature from 2.86 Å it leaves 2.70 Å. **A rigid rotation
about the predicted centroid, holding the predicted internal conformer, does not reach the
crystal pose for three quarters of the failures.** FINDING 024's own decomposition already
implied this and it was not read: CYP3A4's error is 74% rigid-body and **26% internal
conformer (1.60 Å)**, plus 0.94 Å of translation. Rotation is the largest term, not the
only one, and "only the orientation is wrong" is a statement about a *mean*, not a recipe.

---

## Why selection falls: the mechanism, measured

| | |
|---|---|
| surviving orientations per pose | **mean 13.6, median 0.3** of 1,024 |
| ligands with **zero** surviving rotations anywhere | **4 of 87** |
| mean LDDT-PLI of a surviving candidate | **0.5037** against 0.5769 for an original pose |
| survivors whose `xeng` beats that ligand's best original pose | **0.17%** |
| ligands where the selector picks a **rotated** candidate at all | **5 of 87** — and it loses on 4 of them |
| poses whose own best surviving rotation beats the pose itself | 16.3% |
| oracle gain, failed ligands (n = 41) | **+0.0215** |
| oracle gain, already-solved ligands (n = 46) | +0.0012 |
| ligands with any oracle gain | 15 of 87 |

Two things follow, and the second was **not** what the pre-registration expected.

**1. The gain is where it should be.** The oracle gain concentrates 18-fold on the ligands
that were failing. That is the correct shape for a real effect and it is why the verdict is
not a flat refutation.

**2. The Chamfer is not "given more ways to be wrong" — it is given five, and takes four
of them.** The pre-registered expectation said a 14× larger pool of plausible-but-
unsupported orientations would give the consensus term far more opportunities to err.
It does not: 99.83% of the injected candidates are rejected outright because a rotated
ligand disagrees more with six independent engines than the model's own pose does. The
selector is **far more robust to pool contamination than predicted**. It changes its mind
on five ligands and is wrong on four. That is a small number of large mistakes, not a
diffuse degradation, and the pre-registration is recorded here as wrong about the mechanism
while right about the outcome.

**3. Ranking improved while selection got worse — for the third time in this project.**

| | original pool | expanded pool |
|---|---|---|
| within-ligand ρ(`xeng`, LDDT-PLI) | −0.258 | **−0.390** |
| ligands with the correct sign | 75.9% | **86.2%** |
| within-ligand sd of LDDT-PLI | 0.067 | 0.083 |
| **top-1 selection** | **0.6164** | **0.6019** |

The Chamfer ranks the expanded pool *better* by every correlational measure and *selects*
worse. Trap 3 of `docs/README.md` — "reading ρ as selection value" — now has a third
instance, and this is the cleanest one, because ρ improved by 0.13 while the metric that
will actually be used fell by 0.015. **Judge on the metric that will be used.**

### The filter works on the ceiling and hurts the selection

Filtered oracle 0.7083 against the size-matched unfiltered 0.7028: the physics filter is
worth **+0.0055 of oracle** over drawing the same number of candidates at random, so it is
not inert. On *selection* it is the other way round — 0.6019 filtered against 0.6132
unfiltered. The filter concentrates the injected candidates in exactly the snug,
sterically-plausible region where a consensus distance is most easily fooled. **A filter
that improves the ceiling can make the selector worse, and reporting only the ceiling would
have inverted the conclusion.**

### Sensitivity, declared in advance and used for nothing

| variant | survivors | oracle | selected |
|---|---|---|---|
| **pre-registered** (clash 2.2 Å, window 1.85–2.55 Å) | 23,688 | 0.7083 | 0.6019 |
| clash 3.0 Å | 2,085 | 0.6981 | 0.6110 |
| coordination window 1.90–2.45 Å | 23,548 | **0.7083** | **0.6019** |

The third row is identical to the first to four decimals, which is the shape that
`too-clean-numbers-are-the-tell` says to check rather than report. It was checked: the
narrow window **does** fire — it removes 155 candidates and admits 15 others — it simply
never changes any ligand's per-ligand maximum or its argmin-`xeng`. That is the inertness
FINDING 008 predicted for the iron anchor (claim I: "the iron anchor is saturated"), now
confirmed a fifth time. The 3.0 Å row is the consistent story in miniature: 11× fewer
injected candidates, almost no oracle gain, and less damage to selection.

---

## The measurement that actually decides this — and it is not a selection number

If the receptor were correct, the true orientation would **fit** inside the model's own
pocket. It does not. Take the crystal ligand, superpose the binding site (matched by
residue number, FINDING 021 / PXR lesson 1), drop it into the **model's** protein, and
measure its closest heavy-atom contact.

| ligand, in which protein | median min contact | fraction below the 2.2 Å cutoff |
|---|---|---|
| **crystal ligand in its own crystal protein** (the control) | **2.70 Å** | **0.0%** (minimum over 87 = 2.22 Å) |
| predicted pose in its own predicted protein | 2.57 Å | 3.8% |
| **crystal ligand, in place, in the MODEL's protein** | **1.42 Å** | **71.2%** |
| **crystal ligand, at the predicted centroid, in the MODEL's protein** | **1.55 Å** | **75.5%** |

The control is the point of that first row. The 2.2 Å cutoff was pre-registered from
chemistry, before any of this was computed, and **not one of the 87 deposited complexes has
a ligand–protein contact below it** — the minimum observed anywhere is 2.22 Å. So the
filter admits 100% of true poses in their native protein, and rejects the same true poses
in the co-folded protein **71% of the time**.

**The model builds a pocket that its own answer fits and the real answer does not.** The
0.73 Å pocket accuracy of FINDING 024 is a **CA** measurement; the side chains are moulded
around the orientation the model chose. Within the model's own protein, orientation is not
under-determined by sterics at all — it is over-determined, to 2.06%, around the wrong
answer.

This also reconciles the apparent contradiction with FINDING 024's "the cavity is large
enough that the correct orientation is not determined by sterics". 024 measured the
**ligand-centre-accessible** volume — grid points more than 3.1 Å from any heavy atom,
i.e. where a *centre* may sit. That is a different question from how a 35-atom rigid body
may be *turned*, and the two answers are 881 Å³ and 2%.

**The honest limit of this claim:** a 1.42 Å contact is not proof that the pocket is
unrepairable — a rotamer flip could clear it, and this experiment holds the protein rigid
by construction. What is proven is narrower and is exactly what was asked: **a search over
ligand orientations inside a frozen co-folded protein cannot find the crystal pose, because
that protein excludes it.** Any successor must move side chains, and that puts it back on
the protein side, where ρ = +0.03 says there is no measured route.

---

## Verdict

**MEASURED**, by the pre-registered rules and by 0.0008.

- **SHIPS** — fails. Selection *falls* by 0.0145; the CI lower bound is negative.
- **MEASURED** — holds: the oracle rises **+0.0108 > +0.0100** while SHIPS fails.
- **REFUTED** — does not formally trigger: the oracle rise is above the bar, and the
  selection drop carries Holm p = 0.239, not < 0.05.

That margin is 0.0008 of LDDT-PLI. On a grid one step coarser, or with the 3.0 Å clash
cutoff, the same experiment reads **REFUTED**. I am recording the rule's output, not
rounding it into a story: the formal verdict is MEASURED and the substantive verdict is
that the lever is dead. Both belong in the record, which is why the threshold was written
down first.

## What this changes

1. **The last unrefuted consequence of FINDING 024 is now refuted, and by its own
   decomposition.** "Only the orientation is wrong" is true of the *mean* and false as a
   recipe: 31 of 41 failures are unreachable by any rigid rotation about the predicted
   centroid, because 26% of the error is internal conformer and 0.94 Å is translation.
   **Stop reading a decomposition's largest term as its only term.**
2. **"The protein is right" must be qualified: right at CA, wrong at the side chains, and
   wrong in the one way that matters.** FINDING 024's ρ(pocket CA, ligand RMSD) = +0.03
   stands, and it now means something sharper than "the protein is not the problem": the
   *backbone* is not the problem and is uncorrelated with the ligand error, while the
   side-chain packing is moulded around the wrong pose and excludes the right one 71% of
   the time. This is induced fit running the wrong way, and it explains why sampling is
   deterministic per ligand (024: 96.2% of variance between ligands, ρ = +0.94 seed to
   seed): each seed re-derives the same self-consistent wrong pocket.
3. **Pool expansion is closed as a strategy on this target.** Three independent attempts
   now — a second engine (013), the sampler sweep (016), and orientation search (this one)
   — have added oracle that selection cannot reach. The prize in CLAIM F is not reachable
   by making the pool bigger.
4. **The incumbent is stronger than it was given credit for.** It rejected 99.83% of
   1.78 million injected orientations and changed its mind on five ligands out of 87. That
   is the sixth consecutive line where the honest answer is that
   `cypstruct.xengine.select()` still wins, and the first where it was actively attacked.
5. **For the world model:** a physics filter that raises the ceiling can lower the score.
   Report both, always, and in that order.

## Method notes worth keeping

- **C1 at 1.1e-16 over all 1,740 rows** is what makes every comparison here a comparison
  with the shipped pipeline rather than with a reimplementation of it. It cost twenty lines
  and it is the reason the negative is believable.
- **All 1,024 rotations were scored for every pose**, not only the survivors. That is what
  made the size-matched unfiltered null free, and what allowed the "whole grid, no filter"
  ceiling — the number that killed the lever — to be read off without a second run.
- **The symmetry mapping was frozen per pose**, not re-maximised per candidate, so the
  metric is exactly `cypstruct.pose.lddt_pli` as shipped. Re-maximising would have raised
  the rotated candidates' scores and only the rotated candidates' scores.
- **The 2.2 Å cutoff was set from chemistry before the experiment and validated after**:
  0 of 87 crystal complexes violate it. Setting it afterwards from the same crystals would
  have been leakage; setting it first and checking it afterwards is the control.
- **Three comparisons were pre-registered, three were made, three are reported**, with
  Holm over all three.
