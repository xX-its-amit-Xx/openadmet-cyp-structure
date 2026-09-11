# Finding 001 — steering is a null, and Boltz's confidence is worse than random

**Date:** 2026-09-11 · **Run:** `val87b` · **n = 87 ligands, 3,360 poses, all scored**
**Data:** `data/processed/poses_scored_val87b.csv`, `pool_scorecard.json`

This supersedes the two-ligand smoke result recorded earlier. That result suggested
steering lifted the pool ceiling by +0.12 LDDT-PLI. **At n=87 it does not.**

---

## The numbers

| | steered | unsteered |
|---|---|---|
| pool oracle (LDDT-PLI) | 0.6952 | **0.6975** |
| selected by Boltz confidence | 0.5428 | 0.5706 |
| **mean of a random pose** | **0.5585** | **0.5769** |
| worst pose in pool | 0.4205 | 0.4338 |
| selection gap to oracle | 0.1524 | 0.1269 |
| mean BiSyRMSD (Å) | 4.029 | 3.831 |
| ligands with a sub-2 Å pose | 48.1% | 52.9% |
| poses inside the Fe window | 90.2% | 84.1% |

---

## 1. Steering the donor contact is a null, slightly negative

Oracle −0.0023, selected −0.0278, sub-2 Å ligand rate −4.7 points. Every sign is
against it.

**Why it was redundant:** the unsteered arm already coordinates. Its median Fe–donor
distance is **2.229 Å**, against a crystallographic median of **2.20 Å**, with 84% of
poses inside the 1.90–2.45 Å window. Forcing the contact moved that to 90% and bought
nothing, because the geometry was already right.

**What actually fixed coordination is in both arms:** the explicit `bond` constraint from
Cys442 SG to the heme iron, plus a 6,979-sequence MSA. Anchor the cofactor correctly and
Boltz places the donor nitrogen essentially on the crystallographic distance by itself.

**So a claim in this repo's own CLAUDE.md was wrong.** The sibling repo's pilot put zero
of 120 poses inside the coordination window, closest approach 2.80 Å, and that was
generalised to "co-folders cannot reach Fe-coordination geometry." They can. That pilot
did not bond the heme to the axial cysteine and used no deep MSA. The failure was in the
setup, not in the model.

## 2. Boltz's confidence does not rank poses — it anti-ranks them

Within-ligand Spearman correlation between `complex_ipde` and true LDDT-PLI:

| arm | mean ρ | median ρ | ligands with ρ > 0 |
|---|---|---|---|
| steered | **−0.092** | −0.086 | 36 of 81 |
| unsteered | **−0.033** | −0.027 | 43 of 87 |

Positive in about half the ligands, which is chance. And the consequence is blunt:
**picking the highest-confidence pose is worse than picking one at random**, in both
arms (0.5428 vs 0.5585; 0.5706 vs 0.5769).

The pooled correlation looks mildly positive (+0.05, +0.08) because it is confounded by
ligand difficulty — easy ligands get both better poses and higher confidence. Within a
ligand, which is the only comparison selection actually makes, the signal is gone.

**This means the cross-model z-hybrid, the selector that won the PXR structure track, has
no foundation on this target.** It ranks within an engine by that engine's native
confidence. Here that confidence is noise.

## 3. The prize is large, and it is entirely in selection

Oracle **0.6975** against selection **0.5706**: **0.127 LDDT-PLI per ligand** sitting
unclaimed in a pool that already exists and is already paid for. For scale, the winning
entry in the analogous PXR structure track scored **0.564** — below our pool's *current
selection*, and far below its ceiling.

Generation is not the bottleneck. This is the scoring-bottleneck thesis, confirmed at a
sample size that can carry it.

## 4. The coordination term survives — as a *scorer*, not a steerer

| arm | coordinated poses | non-coordinated |
|---|---|---|
| steered | 0.5704 (n=1462) | 0.4482 (n=158) |
| unsteered | 0.5992 (n=1464) | 0.4588 (n=276) |

A pose that coordinates the iron is worth about **+0.14 LDDT-PLI** over one that does
not. That is real, CYP-specific signal — and it is exactly the kind of thing the physics
scorer is built from. The mistake was using it as a constraint on generation, where it
was redundant, instead of as a feature for selection, where the incumbent is empty.

---

## What this changes

1. **Drop the steered arm.** It costs half the GPU budget and returns nothing. Future
   batches run one arm: heme-bonded, deep MSA, more samples.
2. **The physics scorer is now the whole project, not a component of it.** It no longer
   has to beat a strong incumbent; it has to beat random, which the incumbent does not.
3. **Re-examine the tail-rescue rule before porting it.** It ranks by pool confidence,
   and pool confidence here is noise.
4. **Spend the freed budget on samples, not arms.** The oracle rises with pool size and
   the oracle is what a working scorer converts into score.
