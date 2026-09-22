# W002 — The multi-conformer risk is foreclosed by the submission format

**Date:** 2026-09-22 · **Status:** resolved · **Effect:** downgrades the risk I ranked highest

## The worry

OpenADMET's announcement says the cryoEM density for several ligands fits **multiple
mutually exclusive conformations, not one pose**. Our pipeline predicts one pose per
ligand, and `cypstruct.pose.load_structure` calls `st.remove_alternative_conformations()`
(pose.py:117), which silently keeps altloc A. I ranked this as the single item most likely
to invalidate a submission rather than merely underperform.

## Two measurements settle it

**1. The submission format mandates one pose.** From the challenge Space's own
`submission.py`, read directly rather than inferred:

> "Submit a .zip archive containing exactly {STRUCTURE_DATASET_SIZE} .pdb files, **one per
> compound**, named after the compound identifier. Each file must be a full protein-ligand
> complex with the ligand residue named LIG."

There is no multi-conformer submission format. We could not express an ensemble if we
wanted to, and neither can anyone else. Whatever the references contain, reducing them to
a comparable single answer is the organisers' problem and it is identical for every
entrant.

**2. Our altloc handling is a non-issue on the proxy set.** Counting alternate
conformations on the query ligand across all 87 CYP3A4 references:

| | count |
|---|---|
| ligands checked | 87 |
| **with ligand altlocs** | **2** (6BDI_DEJ A/B; 9YK4_PG4, a PEG cryoprotectant) |
| partial occupancy without altlocs | 3 (min occ 0.67–0.92) |

2 of 87. `remove_alternative_conformations()` changes essentially nothing about our proxy
scores.

## What remains unknown, and why it is not actionable

These are **X-ray** crystals. There is no cryoEM CYP3A4 structure in the PDB (all 122
entries are X-ray), so this measurement says nothing about how the released cryoEM
references will encode an ensemble, or how the official metric will score a single pose
against one.

But that uncertainty is now **unactionable in the right way**: the format forecloses any
response we could make to it. There is no build decision waiting on the answer.

## Correction to my own priority call

Yesterday I told the user this was the top risk because it "could invalidate a submission
rather than merely underperform". That was wrong in its premise — a format that permits
only one pose cannot be violated by submitting one pose. The risk is not invalidation; at
worst it is a scoring convention we cannot influence.

The genuine top risk reverts to the one below it: **all our validation is X-ray and the
test set is cryoEM**, with zero evidence the proxy transfers.
