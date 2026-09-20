# FINDING 021 — A residue-numbering offset in MY scorer invented a 46% catastrophe rate

**Date:** 2026-09-20 · **Status:** bug found, all held-out numbers rescored · **Severity:** retracts FINDING 020 outright and revises 018 and 019

## The bug

`lddt_pli` pairs protein atoms by **residue number**. Boltz numbers its output `1..N` from
the input sequence. A crystal uses its own auth numbering — which for CYP3A4 starts near
29, for 3NA0 is offset by 43. Where those disagree, every ligand–protein contact is
compared against the *wrong residue* and the score collapses to exactly 0.0.

It fails silently, and only for the targets whose offset happens to be large. Targets
whose crystal numbering already started near 1 scored correctly, which is exactly why the
result looked like a clean bimodal split rather than a bug.

The tell was there and I read it the wrong way round: **28 of 85 pairs scored *exactly*
0.0000**, and all 15 CYP3A4 pairs failed while aromatase (P11511) scored 0.863 and
P10614 scored 0.909. A model that solves aromatase at 0.9 and fails every CYP3A4 at 0.00
is not describing chemistry.

## The fix and its control

Scan integer offsets, take the one maximising residue-**name** agreement between model
and reference, renumber, rescore. Matching on names keeps it superposition-free and
cannot manufacture agreement — a wrong offset scores near-zero identity.

| pair | old | offset | identity | **new** |
|---|---|---|---|---|
| 3NA0_2DC | 0.000 | +43 | 1.00 | **0.986** |
| 3DSJ_243 | 0.000 | +27 | 0.94 | **0.941** |
| 3B6H_MXD | 0.000 | +10 | 0.94 | **0.834** |
| 3NXU_RIT | 0.003 | +22 | 0.94 | **0.755** |
| 2J0D_ERY | 0.010 | +22 | 0.92 | **0.515** |
| 1IZO_PAM | 0.740 | **0** | 0.99 | 0.740 |
| 2CIB_CM6 | 0.906 | **0** | 0.92 | 0.906 |
| 2JJO_EY5 | 0.855 | **0** | 0.96 | 0.855 |

The already-aligned pairs come back at offset 0, unchanged. That is the control: the fix
is a no-op exactly where nothing was broken.

## What the held-out set actually looks like

| | reported before | **true** |
|---|---|---|
| base mean LDDT-PLI | 0.3477 | **0.8110** |
| median BiSyRMSD | 5.85 Å | **0.716 Å** |
| sub-2 Å | 44/84 | **67/84** |
| catastrophes (<0.1) | 39/84 | **1/84** |
| 15-pose oracle headroom | +0.0204 | **+0.0632** |

Boltz-2 with a real MSA is *strong* on held-out P450 targets — 0.81 mean LDDT-PLI, 0.72 Å
median ligand RMSD. There was never a bimodal failure mode, a capability floor, or a
catastrophe population.

## Consequences for the earlier findings

**FINDING 018 (fine-tuning is negative) — STANDS, effect is ~3× larger.**

| dose | old Δ | **true Δ** | wins | p |
|---|---|---|---|---|
| 87 steps | −0.0038 | **−0.0135** | 19/84 | 1.0e-03 |
| 174 steps | −0.0138 | **−0.0369** | 15/84 | 6.6e-08 |
| 350 steps | −0.0125 | **−0.0301** | 20/84 | 6.3e-05 |

Still monotonic to 174 steps, still plateauing, and now clearly significant. The decision
not to train the other three arms is unaffected.

**FINDING 019 (the heme bond is a null) — STANDS, and is now a cleaner null.**
Δ −0.0069 at **p = 0.43** (was −0.0081 at p = 0.086). The coordination measurement is
unaffected — it never used residue numbering — so "unbonded predictions already reproduce
crystal Fe geometry" holds. What must be struck is the closing claim that *catastrophes
are rotations about a correct anchor*: there are no catastrophes.

**FINDING 020 (catastrophes are not selectable) — RETRACTED IN FULL.** Every premise was
an artifact: the bimodal split, the 39 catastrophes, the absolute 10.9 Å/0.74 Å partition,
the "deterministic per-target failure", the ligand-size and sequence-length predictors
(they predicted *which crystals are numbered oddly*), and the conclusion that selection
had no headroom. Real headroom is +0.0632 over 15 poses.

## The lesson, which this repo already wrote down

`CLAUDE.md` lesson 1 from the PXR campaign: *"Validate the reference frame before trusting
any cross-model metric. Some engine exports landed ~20 Å off the crystal and inflated
every agreement number about twofold before anyone noticed."* This is the same failure in
a new coordinate: not a spatial frame, a **numbering** frame.

`cypstruct.pose.align_by_residue` matches by residue number "never by array position" —
which is correct and is why it was trusted. But matching by residue number is only safe
when both structures share a numbering origin, and a co-folder's output never does.

**Control added to the scorer:** every pair now records `resnum_offset` and
`seq_identity_at_offset`, and the run reports how many fall below 0.80 identity. All 85
pass at 0.80+. A bad alignment is now visible in the output instead of being absorbed
into a low score.
