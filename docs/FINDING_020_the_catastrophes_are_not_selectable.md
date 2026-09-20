# FINDING 020 — The held-out catastrophes are not selectable, and three explanations for them are dead

**Date:** 2026-09-20 · **Status:** measured on 84–85 held-out P450 pairs · **Verdict:** a hard capability floor, not a selection gap

## 1. There is nothing to select

Three independent prediction sets now exist for the same pairs — unbonded base, bonded
base, and the 350-step fine-tune — at 5 diffusion samples each. Pooling them:

| pool | poses/pair | sample 0 | **oracle** | headroom | oracle catastrophes |
|---|---|---|---|---|---|
| base | 5 | 0.3430 | 0.3608 | +0.0178 | 38 |
| bonded | 5 | 0.3349 | 0.3516 | +0.0167 | 38 |
| fine-tuned | 5 | 0.3301 | 0.3370 | +0.0069 | 38 |
| base + bonded | 10 | 0.3430 | 0.3631 | +0.0201 | 38 |
| **all three** | **15** | 0.3430 | **0.3634** | **+0.0204** | **38** |

Of the 39 pairs whose base pose is a catastrophe, **exactly 1** has any pose among all 15
exceeding 0.3, and the best-in-pool mean across those 39 is **0.0141**.

Fifteen poses from three different setups, and the catastrophe count moves from 39 to 38.
The failure is *deterministic per target*, not stochastic. **This corrects the closing
line of FINDING 019**, which pointed at selection as the lever. That holds on the CYP3A4
pool, where the oracle is 0.6975 against 0.5706 selected — a 0.127 gap. Here the entire
gap is 0.020, and no selector can recover what no sample contains.

## 2. The split is absolute, not a continuum

| group | n | median LDDT-PLI | median BiSyRMSD | sub-2 Å | >10 Å |
|---|---|---|---|---|---|
| catastrophic | 39 | 0.0000 | 10.93 Å | 0/39 | 24/39 |
| ok | 46 | 0.5737 | 0.74 Å | 44/46 | 0/46 |

Cross-tabulated, `lddt < 0.1 & rmsd < 2 Å` is **0** and `lddt ≥ 0.1 & rmsd ≥ 2 Å` is **2**.
28 of the 39 score *exactly* 0.0 — not one preserved contact. The model either places the
ligand at 0.74 Å or throws it out of the site entirely. There is no middle.

## 3. Three explanations, all dead

**Face flip — dead.** The proximal face is where the Cys thiolate binds and where
deposited structures put 0.7% of ligand atoms. Projecting each ligand centroid onto the
Fe→SG axis: **0 of 85 poses are proximal**, catastrophic and correct alike (median
projection −6.47 Å vs −5.34 Å). Boltz never flips the heme face.

**Wrong site — dead.** If the crystal ligand sat somewhere other than the heme, a model
that always docks at the heme would fail exactly there. But **84 of 85 crystal ligands
are within 6 Å of the iron.** The one that is not does score 0.0 — a single case, not an
explanation for 39.

**Broken iron coordination — dead, and inverted.** Per-pair, the catastrophic poses
coordinate the iron *more* than the good ones: median Fe–ligand 3.415 Å with 41.0% inside
2.45 Å, against 3.959 Å and 17.4% for the correct poses. The failing poses are anchored
at the iron, on the correct face, and still land 10.9 Å from the truth.

## 4. What does separate them

| | catastrophic | ok |
|---|---|---|
| ligand heavy atoms (median) | **26** | 21 |
| sequence length (median) | **480** | 411 |
| crystal ligand coordinates Fe | 0.436 | 0.217 |

Bigger ligands, longer proteins, and — counter-intuitively — ligands that genuinely
coordinate the iron in the crystal. These are cheap, non-leaky priors available before
any pose is generated, so they can flag which targets are likely unrecoverable. They do
not fix anything.

## 5. What this means for the campaign

The held-out set contains a hard floor: ~46% of these targets cannot be solved by this
model at any dose, any constraint, or any number of samples. Four interventions
(fine-tune ×3, heme bond) and 15 poses per pair have not moved a single one.

Two honest options remain, and neither is "try another generation knob":

1. **Accept the floor and work the other 54%**, where median RMSD is already 0.74 Å.
   Selection there is worth little (+0.02), so the win would have to come from ranking
   the *pairs*, not the poses — i.e. knowing when to trust a prediction.
2. **Check whether the floor is CYP3A4's problem at all.** This set is 85 diverse P450
   targets. The challenge is one target with known structure. A 46% catastrophe rate
   across the superfamily may simply not describe CYP3A4, where the repo's own pool
   selects at 0.5706. **That is the cheapest next measurement and it should come first.**
