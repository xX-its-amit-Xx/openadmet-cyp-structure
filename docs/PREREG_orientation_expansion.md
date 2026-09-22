# PRE-REGISTRATION — orientation expansion of the CYP3A4 pool

**Written:** 2026-09-22, **before any candidate was generated or scored.**
**Frozen:** nothing below may be changed after the first number is produced. If a result
contradicts this document, the contradiction is reported, not the document edited.

---

## Why this experiment exists

FINDING 024 decomposed CYP3A4's failure and found that

- the pocket the model builds is right (pocket-lining CA **0.73 Å**),
- the ligand centroid is right (**1.07 Å**, isotropic scatter, no preferred direction),
- and only the **orientation** is wrong (**30.1°** against **6.1°** for the rest of the
  P450 family),
- with ρ(protein error, ligand error) = **+0.030, p = 0.80** inside CYP3A4 against
  **+0.510, p = 1e-24** everywhere else.

Every protein-conditioning remedy is therefore refuted by mechanism, and FINDING 026
closed the one selector-side proposal that followed from the decomposition.

**The one lever the decomposition leaves standing:** if the receptor is correct and only
the rotation is wrong, the missing poses are already reachable *inside the model's own
protein* by searching orientation space — with no new network inference at all. This
experiment measures that ceiling and, more importantly, measures whether selection can
reach any of it.

---

## Hypotheses

**H1 (oracle).** Rigidly rotating each pose's ligand about its own centroid, inside that
pose's own protein, and keeping only candidates that pass prediction-only physics filters,
raises the pool oracle LDDT-PLI.

**H2 (rescue).** The number of the 87 pairs that never reach BiSyRMSD < 2 Å falls.

**H3 (selection — the hypothesis that decides the verdict).** The **shipped** selector
`cypstruct.xengine.select()`, run unchanged on the expanded pool, selects a higher mean
LDDT-PLI than the same selector run on the original pool, paired per ligand.

H1 and H2 are *ceilings*. An oracle rises mechanically with pool size, so H1 alone proves
nothing (CLAUDE.md claim F, FINDING 013: a union pool added +0.0375 of oracle that
selection could not reach). **H3 is the experiment.**

---

## Data, fixed

| | |
|---|---|
| pool | `D:/cyp_scratch/val87b_unsteered` — 87 ligands × 20 Boltz-2 samples = **1,740 poses** |
| truth | deposited crystals in `data/reference/rcsb`, ligand pinned by CCD code, **one chain** |
| ligand table | `data/processed/validation_ligands.csv` (id, smiles, pdb) |
| references for the selector | `data/processed/reference_set_cyp3a4.npz` (frozen, independent engines) |
| incumbent numbers to reproduce | `data/processed/poses_scored_val87b.csv`, `data/processed/xeng_val87b.csv` |

## Metrics, fixed

- **LDDT-PLI** — `cypstruct.pose.lddt_pli`, inclusion radius 6.0 Å, thresholds
  (0.5, 1.0, 2.0, 4.0), using **the pose's own** `best_ligand_mapping` permutation, held
  fixed across every rotation of that pose. This is exactly the shipped metric; the
  alternative of re-maximising over graph automorphisms per candidate is reported as a
  secondary column only, and is not used to decide anything.
- **BiSyRMSD** — `cypstruct.pose.bisy_rmsd`, binding-site superposition computed once per
  pose (the protein is unchanged by a ligand rotation) and reused for that pose's
  candidates, same permutation.

## Controls that must pass before any result is believed

- **C1.** Recomputed original-pose LDDT-PLI must reproduce the unsteered rows of
  `poses_scored_val87b.csv` to |Δ| < 1e-9 on at least 99% of the 1,740 rows.
- **C2 (FINDING 021).** `renumber_to_reference` is run on every pose; the chosen offset and
  the residue-name identity are reported per ligand. **No pair may score exactly 0.0**
  except where BiSyRMSD > 10 Å, i.e. a genuine ejection rather than a numbering artefact.
  Exact zeros are the known symptom and are checked first.
- **C3.** Recomputed `xeng` on the original poses must reproduce `xeng_val87b.csv` to
  |Δ| < 1e-6, so the selector under test genuinely is the incumbent and not a lookalike.
- **C4.** Pose and reference heavy-atom counts must agree on every scored row.

## Method, fixed

1. **Grid.** `N = 1024` orientations from **super-Fibonacci sampling of SO(3)**
   (Alexa, CVPR 2022) — deterministic, no RNG, identical for every pose. The realised mean
   and median nearest-neighbour geodesic angle are reported.
2. **Application.** Each candidate is `R · (L − c) + c` with `c` the ligand's own centroid.
   **The centroid is held fixed**, because FINDING 024 measured it as already correct.
   The protein, the heme and the internal ligand conformer are untouched.
3. **Filters — prediction-only, no crystal, fixed a priori from chemistry:**
   - **(a) Clash.** Reject a candidate if any rotated ligand heavy atom lies within
     **2.2 Å** of any protein or heme heavy atom, **excluding the heme Fe**. 2.2 Å is
     below any non-bonded heavy-atom contact observed in organic crystal structures; the
     Fe is excluded because a dative Fe–donor bond is 1.85–2.55 Å *by construction* and
     would otherwise be scored as a clash.
   - **(b) Coordination.** If the **unrotated** pose satisfies
     `cypstruct.qmscore.geometry`'s `is_coordinated` (Fe–donor ∈ [`COORD_LO`, `COORD_HI`]
     = [1.85, 2.55] Å and S(Cys442)–Fe–donor ≥ 150°), then the candidate must satisfy it
     too. If the unrotated pose does **not** coordinate, no coordination constraint is
     applied to its candidates. The window is the shipped p1–p99 acceptance window of
     FINDING 008, **not** the p5–p95 window 1.90–2.45 Å quoted in CLAUDE.md claim B —
     FINDING 008 measured that the latter rejects 10% of genuinely coordinated crystal
     poses and is a category error as an acceptance test.
   - The overall filter **pass rate** is reported.
   - **Sensitivity, declared now and reported, never used to decide:** clash at 3.0 Å;
     coordination window 1.90–2.45 Å.
4. **Expanded pool** = the original 20 poses of a ligand + every surviving rotation of
   every one of those 20 poses.
5. **Size-matched null pool** = the original 20 poses + **the same number** of candidates
   per pose, drawn uniformly at random from the same 1,024-point grid with **no filter at
   all** (fixed seed 20260922). Matched in size, unmatched in physics. If the filtered
   expansion is not better than this, the filter is doing nothing and must be said to be
   doing nothing.
6. **Selector.** `cypstruct.xengine.select()`, unchanged, on
   `xeng = mean symmetric Chamfer in the pose's own heme frame to that ligand's frozen
   reference poses`. The heme frame is unchanged by a ligand rotation, so the same frame
   serves every candidate. The selector has **no fitted parameters**, so there is nothing
   to hold out: leave-one-ligand-out is *vacuous here* and is stated as such rather than
   performed decoratively. No cross-ligand quantity is fitted anywhere in this experiment.
7. **Random baseline and ties.** Random-pose baselines and any tie-breaking are averaged
   over **64 draws** with a fixed seed.
8. **Statistics.** Per-ligand paired differences over the 87 ligands; 20,000-draw
   bootstrap 95% CI on the paired mean; two-sided Wilcoxon signed-rank; **Holm** over the
   three primary comparisons (filtered vs original, unfiltered vs original, filtered vs
   unfiltered). Within-ligand statistics are reported alongside every pooled mean
   (FINDING 026 / the within-ligand-variance lesson), never a pooled mean alone.

## Verdict rules, fixed

- **SHIPS** iff **all** of: paired mean Δ(selected LDDT-PLI, filtered-expanded − original)
  **> 0**; bootstrap 95% CI **lower bound > 0**; Wilcoxon **Holm-adjusted p < 0.05**; and
  the filtered expansion is **≥** the size-matched unfiltered expansion on the same
  statistic.
- **MEASURED** iff the pool oracle rises by **> +0.010** LDDT-PLI while the SHIPS
  criterion fails. This is the outcome CLAUDE.md claim F predicts and it is a publishable
  negative: *the poses exist in the model's own protein and the selector cannot find them.*
- **REFUTED** iff the oracle rise is **≤ +0.010**, or selection **falls** with Holm-adjusted
  Wilcoxon p < 0.05.
- Reported regardless of verdict: the change in the count of the 87 pairs with no sub-2 Å
  pose, the filter pass rate, and the oracle/selection of the unfiltered null.

## What is *not* allowed

- No threshold, sign, cutoff or grid size may be chosen after seeing a scored result.
  Every number in "Method, fixed" is set above, from chemistry, before the first candidate
  existed.
- No truth-guided refinement of the grid. A local search around the best-scoring
  orientation would inflate the oracle and would not be a pool a selector could ever see.
- No per-ligand or per-fold sign selection for any feature (the LOO-sign-selection trap).
- The expanded pool used for the oracle must be **the identical pool** used for selection.

## Pre-registered expectation

The author's prior, recorded so it can be wrong: **MEASURED**. FINDING 013 (a union pool
adds oracle that selection cannot reach), FINDING 004 (the selector tracks the oracle at
only +0.0125 per doubling) and claim F together predict a large oracle rise and a
selection change indistinguishable from zero or slightly negative, because a ~50× larger
pool of physically plausible but consensus-unsupported orientations gives the Chamfer far
more ways to be wrong.
