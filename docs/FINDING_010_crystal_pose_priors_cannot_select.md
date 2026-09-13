# FINDING 010 — a prior over where crystal ligands sit cannot rank predicted poses

**Date:** 2026-09-13
**Prior fitted on:** 747 in-pocket crystal poses from **non-CYP3A4** P450s
**Tested on:** 1,196 unsteered Boltz poses, 60 CYP3A4 ligands
**Scripts:** `scripts/structure/occupancy_prior.py`, `scripts/structure/orientation_prior.py`

---

## The idea, and why it was worth trying

Findings G and I together say something specific: coordination is worth +0.14 LDDT-PLI as
a feature, but the iron anchor is **saturated** — 84% of poses already coordinate, so
nothing measured at the metal discriminates. Whatever separates a good pose from a bad one
must live in where the *rest* of the molecule goes.

The P450 harvest made it possible to ask that empirically rather than from theory: put
every crystal ligand into its own heme frame, build a density, and score a predicted pose
by how typical its atoms are under it.

The design also gives an unusually clean leakage guarantee. The prior is fitted on
**non-CYP3A4** P450s and tested on CYP3A4 predictions, so nothing about a CYP3A4 ligand —
not the ligand, not the entry, not even the protein — can enter its own score. A positive
result would have been transferable P450 physics, exactly the claim the scorer needs.

---

## Both versions are null

| variant | within-ligand ρ | frac ρ > 0 | gain vs random | empirical p |
|---|---|---|---|---|
| cylindrical `(r, z)` atom density | +0.044 | 0.60 | **+0.0018** | 0.36 |
| full frame: centroid + azimuth + tilt + R_gyr, KDE | +0.022 | 0.53 | **+0.0008** | 0.40 |

Null measured on this same pool: sd 0.0102, 95th pct **+0.0142**, 99th pct **+0.0209**
(random baseline 0.5855, oracle 0.7089). Both gains are an order of magnitude inside it.

The second version exists because the first had an obvious defect: `(r, z)` averages over
the azimuth, and the azimuth is the one dimension in which poses sharing a coordinated
anchor actually differ. Completing the frame — an in-plane x-axis from the heme's
propionate oxygens, the only O atoms in a heme, so the direction is fixed by chemistry
rather than atom names — and scoring placement *and* orientation under a KDE made it
**worse**, not better. The defect was not what limited it.

---

## Why it fails, which is the part worth keeping

Both variants describe a pose's **gross placement and orientation** in the heme frame. The
predicted poses are all grossly plausible: they coordinate the iron, they sit in the
pocket, they lie on the distal face. A superfamily-averaged density cannot separate them
because the thing it measures is the thing they already agree on — the same saturation
that made recalibrating the coordination window worth −0.0002.

The discriminating error is finer and **ligand-specific**: which substituent points into
which subpocket, for *this* molecule. A density averaged over 369 chemically unrelated
ligands has no way to express that, because the answer depends on the ligand's own
topology rather than on where P450 ligands sit in general.

This retro-explains why the one selector that does work — pocket contacts plus consensus
RMSD, +0.0220 held-out — is built entirely from **within-ligand** comparisons. Consensus
asks whether independent samples *of the same molecule* agree. It never needs a
population-level notion of a correct pose, and FINDING 002's addendum already warned that
between-ligand effects do not transfer to within-ligand selection. A crystal-pose prior is
a between-ligand statistic wearing a physics costume.

**Retired.** No third variant. Adding bins, elements, or ligand-size conditioning all vary
the same saturated quantity, and FINDING 007's arithmetic is unforgiving: with ~20 tries
the best random feature scores +0.0138 to +0.0160, so a family that produces +0.0018 and
+0.0008 is not a near miss to be tuned.

**What this does NOT rule out:** a term that compares a pose against *its own ligand's*
alternatives (consensus, strain, per-ligand QM), or against that ligand's specific
chemistry. That is where the remaining upside is, and it is what
`docs/QM_SCORER_DESIGN.md` tier 1 — which is per-ligand by construction — was already
scoped to do.
