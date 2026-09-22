# PRE-REGISTRATION — induced-fit DEMAND as a selection feature

**Written:** 2026-09-22, **before any score was computed.** Committed before the feature
script was run. Not to be edited afterwards: if the result contradicts anything below, the
contradiction is reported in FINDING 029 and this file stays as written.

---

## 1. The question, and why it could not be asked before

FINDING 028 measured that the co-folded CYP3A4 pocket is effectively **rigid**: chi1 sd
**0.63°** across 87 different ligands against **20.63°** in the corresponding crystals,
**0.077 Å** of ligand-to-ligand side-chain motion against **0.723 Å**. Every scoring
experiment in this repo — including FINDING 025's 34-term physics scorer — was computed
against a receptor with no give. A clash term against a frozen pocket mostly measures the
pocket.

The question asked here is the inverse: **how much induced fit does each pose DEMAND?**
A physically correct pose should be accommodated by a small, plausible, low-strain rotamer
adjustment. A wrong pose should demand an implausible one, or none at all. The demand is
the feature, and it is computable from a prediction alone — no crystal, no truth.

### The degeneracy that forces the design, stated up front

A predicted pose sits inside the protein the same model built around it, so it does not
clash with it: FINDING 028 measured the predicted pose against its own model protein at a
median closest contact of **2.571 Å**, with only **3.79%** of poses below the 2.2 Å cut.
Repacking a pose against its own receptor would therefore find nothing to do on 96% of
poses and the feature would be constant. That is not a result, it is a definition.

So the receptor is **de-moulded** first, and the demand is measured against the de-moulded
pocket. De-moulding is one fixed, universal rule with no tuning:

> **Every rotameric side chain in the receptor set has its chi1 snapped to the nearest
> staggered well, {−60°, +60°, 180°}.**

sp3–sp3 torsions are staggered; this is textbook and is not a parameter chosen here. The
de-moulded pocket is what a pocket looks like before a model moulds it around a chosen
orientation, and the repack cost from that state is the induced fit the pose demands.

## 2. Receptor sets — two, fixed in advance, reported side by side

| set | definition |
|---|---|
| **BLOCKERS** | the six residues FINDING 028 named — **Arg212, Phe215, Phe304, Phe241, Phe213, Ser119** — which between them do 97% of the exclusion |
| **POCKET** | every rotameric residue with a heavy atom within **6.0 Å** of the *predicted* ligand (FINDING 028's `contact_radius_A`) |

Both are computed for every pose, and FINDING 029 reports whether the answer depends on
the choice. Residue numbering is resolved by **residue-name agreement against the canonical
CYP3A4 sequence** (UniProt P08684), never against a crystal, so the method is non-leaky and
transfers to a blind target.

## 3. Acceptance criteria and grid — FINDING 028's, taken verbatim, not tuned

| | value | provenance |
|---|---|---|
| ligand clash cut | **2.2 Å** | `PREREG["clash_cut_A"]`, validated on 87 crystals (0 violations) |
| self-clash cut | **2.6 Å** | `PREREG["self_clash_cut_A"]` |
| chi1 × chi2 grid | **10°** | `PREREG["chi_grid_deg"]` |
| chi3 escalation | **20°**, only when chi1 × chi2 fails | `PREREG["chi3_grid_deg"]` |
| contact radius | **6.0 Å** | `PREREG["contact_radius_A"]` |
| accepted solution | the **minimum max\|delta\|** state that clears the ligand and stays self-consistent | FINDING 028 stage 4 |

Self-clash is tested against the fixed environment (everything except the movable side
chains, excluding own residue and i ± 1) plus the heme with Fe excluded — 028's rule.

**One declared correction to the inherited code.** 028's `scan()` applies chi1 before
chi2 while using the *original* CB→CG axis for chi2; because chi1 moves CG, the composite
distorts the CG–CD bond. The torsions are applied distal-first here (chi2, then chi3, then
chi1) so the composite is exact. This is a correctness fix, not a tuning knob, and it is
declared before any number was read.

## 4. Features and their signs — fixed a priori by chemistry

Sign is the direction that should indicate a **better** pose. No sign is chosen by looking
at a fold; `loo-sign-selection-fakes-negatives` is why this paragraph exists.

| feature | definition | sign |
|---|---|---|
| `n_demand` | residues of the set clashing with the ligand (< 2.2 Å) in the de-moulded pocket | **−** |
| `sum_pert` | Σ over cleared residues of the accepted max\|chi delta\| (deg) | **−** |
| `max_pert` | max over cleared residues of the accepted max\|chi delta\| (deg) | **−** |
| `n_uncleared` | residues no grid point can accommodate | **−** |
| `unsolved` | 1 if the concerted repacked pocket still contacts the ligand below 2.2 Å | **−** |
| `clash_before` | min ligand↔receptor distance in the de-moulded pocket (Å) | **+** |
| `clash_after` | min ligand↔receptor distance after the concerted repack (Å) | **+** |
| `d_clash` | `clash_after − clash_before`: how far the pocket had to open | **−** |
| `rot_strain` | Σ over cleared residues of \|chi1_final − nearest staggered well\| (deg, ≤ 60) — landing *between* wells is strained, a 120° flip is not | **−** |
| `n_skel_block` | residues whose ligand clash involves a backbone atom or CB — no torsion can fix it | **−** |
| `skeleton_clear` | min ligand distance to the immovable skeleton (all backbone + CB + heme, Fe excluded) (Å) | **+** |

`rot_implaus_model` (Σ \|chi1_model − nearest well\| over the set, on the un-de-moulded
model) is computed as a descriptive column only. It is a property of the receptor, so it is
expected to be near-constant within a ligand; its within-ligand CV is reported and it is
excluded from single-term selection if it is constant within > 50% of ligands.

An **ENSEMBLE** per receptor set: leave-one-**ligand**-out `HistGradientBoostingRegressor`
over all features of that set, identical fitting to FINDING 025's so the two are comparable.

## 5. Evaluation protocol

- Pool: `D:/cyp_scratch/val87b_unsteered`, 87 ligands × 20 Boltz-2 samples = 1,740 poses,
  truth from `data/processed/poses_scored_val87b.csv` (arm = unsteered).
- Selection: top-1 per ligand, **ties broken at random, averaged over 64 draws**.
- Random baseline: the per-ligand mean LDDT-PLI. **Gain = selected − random.**
- Reported with every pooled number: the **pool oracle**, mean **within-ligand Spearman ρ**
  against LDDT-PLI, and the **fraction of ligands with the correct sign**.

### Controls that must be run first, and must fire

1. **FINDING 021 numbering.** Offset against the canonical CYP3A4 sequence must be 0 on
   all 1,740 poses with high residue-name identity; **no (ligand, sample) may score
   exactly 0.0** LDDT-PLI for a numbering reason. Reported with counts.
2. **Incumbent reproduction.** `cypstruct.xengine` recomputed from the frozen reference
   set must reproduce `data/processed/xeng_val87b.csv` to < 1e-6 on all 1,740 rows.
3. **Every criterion must be shown to fire**, with counts: how often the 2.2 Å cut binds,
   how often the 2.6 Å self-clash cut binds, how often the chi3 escalation is reached, how
   often no rotamer clears. A filter that never fires is reported as such.

### The wrong-pose sanity control — reported BEFORE any selector number

For each ligand's sample 0 pose, generate rigid rotations of the ligand about its own
centroid inside its own protein at two fixed magnitudes, **15°** and **90°**, 5 random axes
each, and recompute the POCKET demand features.

**Required to proceed:** at the 90° magnitude, the rotated poses must demand *more* induced
fit than the originals on **at least 2 of** {`n_demand`, `clash_before`, `n_uncleared`,
`skeleton_clear`}, each with Mann–Whitney **p < 0.01** in the pre-declared direction, and
the paired per-pose sign correct on **≥ 60%** of poses.

**If this fails the feature is measuring nothing and the work stops there and says so.**

Noted in advance: a synthetic rigid rotation is an easy target — `rigid-perturbation-detectors`
records detectors that score AUC 0.96 on synthetic rotations and below chance on real
errors. Passing this control is therefore **necessary and not sufficient**, and FINDING 029
will say so however it lands.

### The three bars, in order

**BAR 1 — beat random, against a null computed for THESE features.** Three nulls, all
re-derived on this pool rather than quoted:

- **N1 random selection**, 4,000 draws → p95, p99.
- **N2 within-ligand feature shuffle**, 200 reps: permute feature rows inside each ligand,
  leaving the labels alone, and rerun the identical pipeline. Preserves every marginal and
  the exact dimensionality; destroys the pose-to-feature pairing.
- **N3 matched-dimensionality Gaussian**, 200 reps, same feature count, same LOO fit.

An **unfitted single term** must clear N1's p95. The **fitted ensemble** must clear the p95
of **both** N2 and N3.

**BAR 2 — the paired test against the shipped incumbent, on identical poses.** The
incumbent is `cypstruct.xengine.select()` (`−z(xeng)` within ligand). Reported: mean
per-ligand difference, **bootstrap 95% CI** (10,000 resamples), **Wilcoxon** p, and the
**number of ligands where the two selectors pick the same pose** (ties). A large raw gain
with most ligands tied is not a result.

**Pass requires:** mean difference > 0, bootstrap 95% CI excluding 0, Wilcoxon p < 0.05 on
the untied ligands.

**BAR 3 — complementarity.** Pearson r between the per-ligand gain and the incumbent's
per-ligand shortfall (oracle − incumbent) must be **> 0 with p < 0.05**. FINDING 025 died
here at r = −0.047.

### Verdicts, decided by the bars and not by the size of a number

| outcome | verdict |
|---|---|
| wrong-pose control fails | **REFUTED** — the feature measures nothing; stop |
| BAR 1 fails | **REFUTED** as a selector |
| BAR 1 passes, BAR 2 fails | **MEASURED** — real but absorbed; logged, does not ship |
| BARS 1, 2, 3 all pass | **SHIPS** |

### Traps this pre-registration is guarding against

- No cutoff is chosen by looking at a scored result. Every threshold above comes from
  FINDING 028 or is a textbook constant.
- Within-ligand statistics are reported next to every pooled number; a between-ligand
  effect has repeatedly been worthless here (FINDING 024, `within-ligand-variance`).
- Every feature's within-ligand constancy is reported before it is tested
  (`mode-level-features-cannot-select`).
- An sd of exactly 0, a filter that never fires or a perfect correlation is reported as a
  defect, not as a clean result (`too-clean-numbers-are-the-tell`).
