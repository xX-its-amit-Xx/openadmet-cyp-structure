# FINDING 029 — induced-fit DEMAND: the probe works, the feature does not select

**Date:** 2026-09-22 · **Status:** measured, zero new inference, CPU-only, ~2 hours on
12 shared cores · **Pre-registered:** `docs/PREREG_induced_fit_demand.md`, committed
before any score was computed · **Verdict:** **MEASURED** by the pre-registered rule —
and **nothing ships**.

The repack machinery does exactly what it was built to do: a 90° rotation of a ligand
inside its own protein takes the demand from **0 residues to 4** and the closest contact
from **2.52 Å to 0.79 Å**, at p = 3 × 10⁻¹⁴². That detector is worth **nothing** against
the errors the pool actually makes. Not one of 31 testable demand terms clears the
+0.0140 noise floor; the best single term is **+0.0012**, and the largest single effect in
the finding, −0.0632, points the *wrong way*. The 44-column fitted ensemble reaches
+0.0228 — the 99th percentile of its own matched-dimensionality noise null, whose 200
draws reached **+0.0264** — and is beaten by the zero-parameter incumbent by **0.0167
LDDT-PLI** with only **5 of 87 ligands tied**. Its per-ligand outcome correlates **+0.858**
with the incumbent's and its gain correlates **−0.078** with where the incumbent errs.

This is `rigid-perturbation-detectors` reproducing on hand-built physics rather than on a
neural denoiser, and it is one more candidate that clears the first bar and dies
at the paired one.

---

## What was asked

FINDING 028 measured that the co-folded CYP3A4 pocket is **rigid**: chi1 sd **0.63°**
across 87 different ligands where the crystals give **20.63°**, side chains moving
**0.077 Å** from ligand to ligand against **0.723 Å**. Every scoring experiment in this
repo, FINDING 025's 34-term physics scorer included, was therefore computed against a
receptor with no give — and a clash term against a frozen pocket mostly measures the
pocket.

So: **how much induced fit does each pose DEMAND?** Repack the pocket around each pose
independently and make the repack cost the feature. A physically correct pose should be
accommodated by a small, plausible, low-strain rotamer adjustment; a wrong one should
demand an implausible adjustment, or none at all. Nothing in the quantity touches a
crystal, so unlike FINDING 028 it is admissible as a selection feature.

## What was done

`scripts/structure/induced_fit_demand.py`, seven stages, importing FINDING 028's chemistry
tables and stage-4 acceptance criteria directly from `side_chain_diagnosis.py` so that a
divergence is impossible.

| | |
|---|---|
| pool | `val87b_unsteered` — 87 ligands × 20 Boltz-2 samples = **1,740 poses** |
| truth | `poses_scored_val87b.csv`, arm = unsteered |
| clash cut / self-clash cut | **2.2 Å** / **2.6 Å** — FINDING 028's, verbatim |
| chi grid | chi1 × chi2 at **10°**, escalated to chi1 × chi2 × chi3 at 10/10/**20°** only when chi1 × chi2 fails — FINDING 028's, verbatim |
| accepted repack | the **minimum max\|Δchi\|** state clearing the ligand and staying self-consistent — FINDING 028's rule |
| receptor sets | **BLOCKERS** = the six residues 028 named (Arg212, Phe215, Phe304, Phe241, Phe213, Ser119); **POCKET** = every rotameric residue within 6.0 Å of the *predicted* ligand (median 21 residues) |
| selection | top-1 per ligand, ties broken at random over **64 draws**; gain = selected − per-ligand mean |

### The degeneracy that forced the design, and it is quantified

A predicted pose sits inside the protein the same model built around it. FINDING 028
measured that at a median closest contact of **2.571 Å** with only **3.79%** of poses
below the cut, so repacking a pose against its own receptor finds nothing to do. The
receptor is therefore **de-moulded** first: every rotameric chi1 snapped to the nearest
staggered well, {−60°, +60°, 180°} — one textbook rule, no tuning.

**The pre-registered de-moulding is too weak, and the reason is a number.** The model's
own chi1 angles already sit a median **6.7°** (POCKET) and **7.1°** (BLOCKERS) from the
nearest staggered well, so snapping them there barely moves anything. The moulding is in
*which well* each residue chose, not in how far off a well it sits.

| variant | poses with any demand (`n_demand` > 0) | columns constant within > 50% of ligands |
|---|---|---|
| SELF · BLOCKERS | **4.1%** | 8 of 11 |
| SELF · POCKET | **8.2%** | 3 of 11 |
| CROSS · XBLOCKERS | 100% of ligands vary; 83.7% of poses non-zero | 1 of 11 |
| CROSS · XPOCKET | **83.7%** | 1 of 11 |

**The amendment, declared in code before any LDDT-PLI was read** (commit `9353100`, and
the degeneracy was found on a 120-pose probe that never touched a truth value): the
**CROSS** variant asks the identical question of a receptor that was *not* built around
the pose — place pose *i* in the protein of sibling sample *j* of the same ligand,
superposed on `CYP3A4_RIGID_CORE` (the ligand-free frame declared in `targets.py`),
de-mould, repack, and average the demand over all **19** siblings. Same features, same
a-priori signs, same bars, no new thresholds, and still computable at inference.

### One correctness fix to the inherited code, declared in advance

FINDING 028's `scan()` applies chi1 before chi2 while taking chi2's CB→CG axis from the
**pre-chi1** coordinates. chi1 moves CG, so the composite stretches the CG–CD bond — by
up to 2.6 Å at a 180° chi1. Here the torsions are applied **distal-first** (chi3, then
chi2, then chi1), which is exact because every later torsion's axis atoms are proximal and
have not moved. This is a correctness fix, not a knob, and it was written down before any
number was read. It does not affect FINDING 028's published conclusions, which rest on
*whether* a residue could be cleared far more than on the exact angle.

---

## Controls — all three fired

| control | required | measured |
|---|---|---|
| **FINDING 021 numbering**, resolved against the **canonical CYP3A4 sequence** and never against a crystal | offset 0, high identity | offset **0 on 1,740 / 1,740**, minimum residue-name identity **1.000**; all six named blockers match the sequence on **1,740 / 1,740** |
| **no silent zeros** | no pair at exactly 0.0 for a numbering reason | **4** rows at exactly 0.0, all genuine — PG0/PG4 PEG fragments with BiSyRMSD **23–27 Å**, ligand on the surface; **0** unmapped, **0** null |
| **incumbent reproduction** | `xengine` recomputed from the frozen reference must match the shipped column | **1,740 / 1,740**, max \|Δ\| **2.8 × 10⁻⁸**, median 5.8 × 10⁻⁹; reference depth min 6, median 7, max 12 |

### Every acceptance criterion fires, with counts

Across the feature, cross and rotation stages:

| event | count |
|---|---|
| residue cleared by the chi scan | **21,679** |
| the 2.6 Å self-clash test rejected a candidate block | **13,036** |
| residue blocked *only* by the self-clash test | **1,379** |
| no clearing rotamer at any grid point | **471** |
| chi3 escalation reached | **207** |

No filter is inert and none is saturated.

### A duplication the "too-clean" check caught

Two pairs of single-term gains came back **identical to four decimals**, which is the tell
for one column measured twice. They are: `skeleton_clear` does not depend on the receptor
set (the immovable backbone + CB shell is global), so `BLOCKERS_skeleton_clear` ≡
`POCKET_skeleton_clear` and `XBLOCKERS_skeleton_clear` ≡ `XPOCKET_skeleton_clear`; and in
the BLOCKERS set at most one residue is ever cleared, so `sum_pert` ≡ `max_pert` and
`n_uncleared` ≡ `unsolved` there. Of the 44 modelled columns, **39 are distinct**. Reported
rather than quietly deduplicated, because the identical numbers are what exposed it.

---

## The wrong-pose control — it FIRES, and that is the least interesting result here

Reported before any selector number, as pre-registered. For each ligand, sample 0's pose
is rigidly rotated about its own centroid by 15° and by 90°, five random axes each, and
the POCKET demand is recomputed — in the pose's own receptor (`self`) and in a sibling
sample's receptor (`sibling`), which is the receptor the CROSS features actually score in.
435 rotated poses per magnitude per frame.

| frame · magnitude | `n_demand` | `clash_before` (Å) | `n_uncleared` | `skeleton_clear` (Å) | criteria passing |
|---|---|---|---|---|---|
| self · **90°** | 0 → **4**, p = 2.5e-142, 95.2% | 2.52 → **0.79**, p = 1.1e-131, 99.8% | 0 → **1**, p = 1.8e-94, 67.6% | 2.83 → **1.30**, p = 1.1e-135, 98.4% | **4 / 4** |
| self · 15° | 0 → 1, p = 4.8e-50, 51.7% ✗ | 2.52 → 2.12, p = 1.4e-44, 84.6% | 0 → 0, p = 3.0e-14, 12.2% ✗ | 2.83 → 2.33, p = 5.2e-59, 88.0% | 2 / 4 |
| sibling · **90°** | 0 → **4**, p = 8.4e-116, 91.3% | 2.41 → **0.78**, p = 2.0e-110, 97.5% | 0 → **1**, p = 4.3e-100, 71.3% | 2.79 → **1.17**, p = 1.0e-136, 98.4% | **4 / 4** |
| sibling · 15° | 0 → 1, p = 3.2e-27, 50.8% ✗ | 2.41 → 2.00, p = 6.4e-25, 79.1% | 0 → 0, p = 6.8e-14, 14.3% ✗ | 2.79 → 2.39, p = 1.2e-48, 81.4% | 2 / 4 |

**The control fires** — ≥ 2 criteria at both magnitudes and all four at 90°, in both
frames. The two ✗ cells are honest: at 15° the *counts* are integers and tie far more
often than they differ, so the paired-sign fraction sits near 0.5 even where the
Mann–Whitney is overwhelming. A tie is not a correct sign and was not counted as one.

**And it is worth almost nothing.** `rigid-perturbation-detectors` says exactly this:
detectors ace synthetic rotations (AUC 0.96) and score real errors below chance. The
pre-registration said in advance that passing this control is necessary and not
sufficient. Everything below is the sufficient half, and it fails.

---

## The pool, the incumbent, and the three nulls re-derived here

| | value |
|---|---|
| poses / ligands | 1,740 / 87 |
| **pool oracle** | **0.6975** |
| random selection | **0.5769** |
| **incumbent `cypstruct.xengine.select()`** (no fitted parameters) | **0.6164 — +0.0395**, within-ligand ρ **+0.258**, correct sign on **75.9%** of ligands |

| null | p95 | p99 | max |
|---|---|---|---|
| **N1** random selection, 4,000 draws | **+0.0140** | +0.0195 | +0.0299 |
| **N2** within-ligand feature shuffle, 200 reps, k = 11 | +0.0171 | +0.0198 | +0.0216 |
| **N3** matched-dimensionality Gaussian, 200 reps, k = 11 | +0.0153 | +0.0192 | +0.0214 |

N1's +0.0140 is re-derived from scratch on these poses and lands on FINDING 007's +0.0138.
The noise floor has not moved in seven findings.

---

## BAR 1 — not one single term clears the floor

Sign fixed a priori by chemistry (less demand, less clash depth, less strain, more
clearance = better pose). No fitting, not even a sign. **13 of 44 columns were excluded
before testing** for being constant within more than half the ligands — the pre-registered
guard, and it fires hardest exactly where the SELF variant is degenerate.

| term | gain | within-ρ | frac ρ correct | n | p |
|---|---|---|---|---|---|
| `POCKET_n_demand` | **+0.0012** | +0.041 | 0.58 | 53 | 0.24 |
| `POCKET_sum_pert` | +0.0006 | +0.059 | 0.61 | 51 | 0.99 |
| `POCKET_max_pert` | +0.0005 | +0.060 | 0.61 | 51 | 0.74 |
| `POCKET_rot_strain` | +0.0001 | +0.029 | 0.57 | 47 | 0.99 |
| `POCKET_d_clash` | −0.0005 | +0.044 | 0.60 | 50 | 0.69 |
| `skeleton_clear` (set-independent) | −0.0046 | +0.057 | 0.56 | 87 | 0.47 |
| `XPOCKET_n_uncleared` | −0.0094 | −0.063 | 0.48 | 58 | 0.07 |
| `XPOCKET_unsolved` | −0.0095 | −0.051 | 0.49 | 63 | **0.004** |
| `X*_skeleton_clear` | −0.0132 | +0.008 | 0.54 | 87 | 0.70 |
| `XPOCKET_rot_strain` | −0.0165 | +0.018 | 0.52 | 87 | 0.79 |
| `XPOCKET_n_demand` | −0.0236 | +0.047 | 0.57 | 87 | 0.26 |
| `XPOCKET_clash_before` | −0.0246 | +0.034 | 0.55 | 87 | 0.98 |
| `XBLOCKERS_n_demand` | −0.0268 | +0.007 | 0.55 | 87 | **0.023** |
| `BLOCKERS_clash_before` | −0.0549 | −0.066 | 0.40 | 87 | **0.003** |
| `POCKET_clash_before` | **−0.0632** | −0.051 | 0.36 | 87 | **3e-4** |
| `POCKET_clash_after` | **−0.0634** | −0.060 | 0.36 | 87 | **2e-4** |

(31 terms tested, full table in `induced_fit_selector_val87b.json`.)

**The best single term is +0.0012 against a +0.0140 floor.** The CROSS terms, which are
the ones with real dynamic range, are *worse* than random — `POCKET_clash_before` selects
at 0.5137 where random gives 0.5769, with only **33%** of ligands beating random.

That negative is not a hidden positive. "More clearance is better" is the physics sign;
the data prefer the snugger pose. But FINDING 025's arithmetic applies unchanged — argmin
and argmax of one feature can both sit below the pool mean — and a sign chosen after
seeing this table is exactly the artifact `loo-sign-selection-fakes-negatives` is about.
It is logged, not harvested.

---

## The fitted ensembles, and BAR 2 — which is where it actually dies

Leave-one-**ligand**-out `HistGradientBoostingRegressor`, FINDING 025's fit unchanged.

| model | selected | gain | ρ | frac ρ > 0 | beats random | p |
|---|---|---|---|---|---|---|
| ENSEMBLE BLOCKERS (11) | 0.5774 | +0.0005 | +0.034 | 0.51 | 0.53 | 0.68 |
| ENSEMBLE POCKET (11) | 0.5834 | +0.0065 | +0.020 | 0.45 | 0.58 | 0.18 |
| ENSEMBLE XBLOCKERS (11) | 0.5867 | +0.0099 | +0.095 | 0.63 | 0.56 | 0.27 |
| ENSEMBLE XPOCKET (11) | 0.5879 | +0.0110 | +0.079 | 0.64 | 0.60 | 0.030 |
| **ENSEMBLE ALL (44 columns, 39 distinct)** | **0.5997** | **+0.0228** | +0.153 | 0.71 | 0.66 | 4e-4 |
| **INCUMBENT** | **0.6164** | **+0.0395** | +0.258 | 0.76 | 0.68 | 4e-5 |

The ALL ensemble must be judged against a null of **its own dimensionality**. The first
pass compared it against a k = 11 null — which is the mistake FINDING 025 exists to warn
about — so both nulls were re-run at k = 44.

| null at k = 44 | mean | sd | p95 | p99 | max | p(null ≥ +0.0228) |
|---|---|---|---|---|---|---|
| N2 within-ligand shuffle, 200 reps | +0.0008 | 0.0089 | +0.0147 | +0.0218 | **+0.0266** | **0.010** |
| N3 matched-dimensionality Gaussian, 200 reps | −0.0000 | 0.0096 | +0.0142 | +0.0217 | **+0.0264** | **0.010** |

**BAR 1 passes on the ensemble and fails on every single term** — and it passes the way
FINDING 025's did: 44 noise features fitted identically beat the physics twice in 200
tries, and went *above* it. Two draws in two hundred is the number to quote, not the
Wilcoxon 4 × 10⁻⁴.

### The paired test, on identical poses

| ensemble | mean difference vs incumbent | bootstrap 95% CI | Wilcoxon p | ligands tied | better / worse |
|---|---|---|---|---|---|
| BLOCKERS | −0.0390 | [−0.0637, −0.0158] | 0.013 | 2 | 37 / 48 |
| POCKET | −0.0329 | [−0.0542, −0.0133] | 0.0044 | 4 | 31 / 52 |
| XBLOCKERS | −0.0296 | [−0.0482, −0.0119] | 0.019 | 9 | 37 / 41 |
| XPOCKET | −0.0284 | [−0.0526, −0.0071] | 0.071 | 5 | 35 / 47 |
| **ALL** | **−0.0167** | **[−0.0370, +0.0019]** | 0.264 | **5** | 36 / 46 |

**Every ensemble loses, and the ties do not explain it.** At most **9 of 87** ligands are
tied and at most 9 pick the same pose as the incumbent, so this is not the failure mode
the pre-registration warned about — the two selectors genuinely disagree on 78–86 ligands
and the demand scorer is the one that is wrong more often.

### BAR 3 — complementarity, failed the same way FINDING 025 failed

| ensemble | r(gain, incumbent shortfall) | p | corr of per-ligand outcomes |
|---|---|---|---|
| BLOCKERS | +0.030 | 0.78 | +0.807 |
| POCKET | −0.054 | 0.62 | +0.845 |
| XBLOCKERS | −0.152 | 0.16 | +0.881 |
| XPOCKET | −0.147 | 0.17 | +0.819 |
| ALL | −0.078 | 0.47 | +0.858 |

FINDING 025 died here at r = −0.047 with outcome correlation +0.834. This is the same
number to two decimals, from an entirely different feature vocabulary. **The demand
scorer is not seeing a different failure mode; it is a weaker read of the same one.**

The raw features are, by contrast, close to orthogonal to the incumbent — the mean
within-ligand correlation between `xeng` and any CROSS demand feature is between −0.12 and
+0.24. So the failure is not that the feature duplicates the incumbent. It is that
whatever it measures does not rank poses.

---

## Verdict

**REFUTED as a selector.** The probe is real and the detector works on synthetic
wrongness; against the errors the pool makes it is a null at the single-term level and a
loss at the ensemble level, with no complementarity.

Against the pre-registered decision table this is **MEASURED**, not REFUTED, on a
technicality worth stating plainly: the single-term half of BAR 1 failed outright and only
the **fitted** 44-column ensemble cleared, at p = 0.01 against a null whose own draws went
higher than the observed value. BAR 2 and BAR 3 then failed with no ambiguity. Nothing
here enters the submission pipeline.

Mixing it in does not help either. Rank-averaging the demand ensemble with the incumbent
scores **0.5969**, *below* the incumbent's 0.6164 — mean difference **−0.0195**, bootstrap
CI [−0.0358, −0.0065], Wilcoxon p = 0.012, 24 ligands tied. Every way of combining it
makes the board worse, which is FINDING 025's result repeating verbatim.

### What this licenses, and what it forecloses

**Forecloses.**

1. **Receptor flexibility is not the missing ingredient in CYP3A4 scoring.** The
   standing explanation for FINDING 025's null — "every term was computed against a frozen
   pocket" — is now tested and is not the explanation. Give the pocket a 10°-resolution
   chi1 × chi2 × chi3 repack with physical acceptance criteria and the answer does not
   change. That closes the flexible-receptor branch of `QM_SCORER_DESIGN.md` at the same
   place the rigid one closed.
2. **Clash-based terms are anti-selective on this pool, at any receptor rigidity.**
   `POCKET_clash_before` at −0.0632 is the largest single-term effect in the whole finding
   and it points the wrong way. FINDING 025's `max_clash` and `n_clash` were nulls; with a
   repacked receptor they become significantly harmful. There is nothing left to try on
   the steric axis.
3. **A synthetic-rotation control is not evidence about a selector.** Four criteria at
   p ≤ 1e-94 preceded a −0.0167 paired loss. `rigid-perturbation-detectors` was recorded
   about denoising models; it is now a property of hand-built physics too, and any future
   candidate that offers a rotation-recovery AUC in place of a paired number should be
   read as offering nothing.

**Licenses.**

1. **The repack machinery is reusable and now correct.** `induced_fit_demand.demand()`
   repacks a pocket around an arbitrary ligand placement at ~7 ms, with exact torsion
   composition, and is the fastest way in this repo to ask "does this placement fit".
   FINDING 027's orientation expansion rejected 99.83% of injected poses on a *rigid*
   receptor; the same expansion can now be re-scored on a repacked one, which is the one
   question 027 could not answer.
2. **The incumbent's margin is unchanged and the oracle gap is unchanged.** Selection is
   0.6164 against an oracle of 0.6975; **0.081 LDDT-PLI per ligand is still unclaimed**,
   and six candidates have now failed to claim it. Every one of them scored a pose from
   its own geometry. The only thing that has ever worked here scores a pose against *other
   predictions of the same molecule*. That contrast is now seven findings deep and should
   be the prior for what is tried next.

---

## Reproduce

```
python scripts/structure/induced_fit_demand.py controls         # numbering, zeros, xeng
python scripts/structure/induced_fit_demand.py features  --workers 12   # SELF variant
python scripts/structure/induced_fit_demand.py cross     --workers 12   # CROSS variant
python scripts/structure/induced_fit_demand.py rotate    --workers 12   # wrong-pose control
python scripts/structure/induced_fit_demand.py rotate_report
python scripts/structure/induced_fit_demand.py evaluate  --null-reps 200
python scripts/structure/induced_fit_demand.py addendum  --null-reps 200
```

Artifacts: `data/processed/induced_fit_controls.json`,
`induced_fit_features_val87b.csv` (1,740 × 61),
`induced_fit_rotation_control.json`, `induced_fit_selector_val87b.json`,
`induced_fit_addendum.json`. Per-ligand intermediates stay in the session scratchpad.
