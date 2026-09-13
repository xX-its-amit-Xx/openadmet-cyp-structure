# FINDING 008 — the whole P450 superfamily replicates the coordination thesis, and shows the scorer's windows were rejecting one true pose in ten

**Date:** 2026-09-13
**Data:** PF00067, resolution ≤ 3.2 Å → 585 entries → 874 in-pocket (chain, ligand, copy)
pairs across 369 unique drug-like ligands, of which **463 coordinate the iron**.
For comparison, the CYP3A4-only calibration behind `qmscore/geometry.py` used 116 entries
and 104 coordinated chain observations.
**Scripts:** `scripts/structure/harvest_p450_universe.py`,
`scripts/structure/p450_geometry_pass2.py`

---

## Why go wide

Every negative result in this repo bottoms out in the same place: n = 87. FINDING 007
measured the noise floor on that set at +0.0138 (95th pct for a random feature), which
means a real effect of +0.015 is indistinguishable from luck. Cleverness does not fix
that; more ligands do.

The P450 fold is one superfamily with one cofactor and one axial cysteine. A term that is
real CYP physics should hold across it. A term that only works on the 87 is probably
fitting the 87.

---

## 1. The coordination thesis replicates on an independent set

Claim B of the repo thesis was measured on CYP3A4 alone. It survives a 4.5× larger,
chemically and phylogenetically far more diverse set essentially unchanged:

| term | CYP3A4 (n = 104) | all P450 (n = 463) |
|---|---|---|
| Fe–donor, p5 / p50 / p95 | 1.94 / **2.20** / 2.38 Å | 1.93 / **2.18** / 2.48 Å |
| S(Cys)–Fe–donor, p5 / p50 / p95 | 159.5 / **171.2** / 177.6° | 162.0 / **172.6** / 177.9° |
| donor element | N in 103 of 104 | N in 448, O in 15 |
| Cys SG–Fe, p50 | — | 2.33 Å |

Bacterial P450s (BM3, cam, eryF) are included deliberately: same chemistry, different
pocket, and therefore the strongest available test of whether these numbers are physics
or CYP3A4 idiosyncrasy. They are physics.

---

## 2. The CYP3A4 validation set is biased toward coordinators

**83%** of CYP3A4's deposited ligands coordinate the iron. Across the superfamily it is
**53%** (463 of 874).

This is a caution, not a defect. The CYP3A4 crystallographic record is enriched in
azole-type inhibitors, which is precisely the chemistry that coordinates. Any term
calibrated only on CYP3A4 has seen roughly one type I complex for every five type II, and
FINDING 001-G — coordinating poses score +0.14 LDDT-PLI over non-coordinating ones — was
measured under that enrichment. If the challenge release turns out to be substrate-heavy,
that +0.14 is the number most likely to shrink.

---

## 3. The bug this exposed: a p5–p95 window is not an acceptance test

`COORD_LO, COORD_HI = 1.90, 2.45` was set from the CYP3A4 p5–p95 "with a small margin."
That is a category error. **A p5–p95 interval excludes 10% of the distribution by
construction**, so using one as an acceptance window guarantees rejecting about one true
positive in ten before any pose is ever scored.

Measured against the 463 genuinely coordinated crystal poses:

| window | admits |
|---|---|
| 1.90 – 2.45 (previous) | **89.8%** |
| 1.85 – 2.55 (now) | 97.4% |
| 1.80 – 2.60 | 99.6% |

One coordinated crystal pose in ten was being scored as non-coordinating. Since
coordination is worth +0.14 LDDT-PLI as a selection feature, that is a real loss and not
a rounding detail. The angle window was the same mistake in the other direction —
`>= 155°` was never binding, since p1 is 157.0° and even `>= 145°` admits 100%.

Acceptance windows are now p1–p99. Applied in `qmscore/geometry.py`:
`COORD 1.85–2.55`, `TRANS_ANGLE_MIN 150`, `SOM 2.8–5.8` (type I standoff, n = 411).

---

## 4. The distal-face constraint is harder than reported

The old note read "only 0.7% of ligand atoms fall on the proximal face." The true figure
across 874 in-pocket pairs is that **not one ligand heavy atom does**: the minimum signed
height above the porphyrin plane is +1.72 Å at the 1st percentile, and only 2 of 1,000
measured copies have any atom below the plane — both of them 9.4 Å and 26 Å from the iron,
i.e. surface sites rather than pocket ligands.

So this is not a soft term to be weighted into a score. It is a hard validity filter that
real structures never violate, and any pose that violates it can be discarded outright.

---

## 5. A methodological trap worth recording

Verifying §4 produced an apparent contradiction: the cached atoms said two in-pocket pairs
had atoms below the plane, the summary table said zero. The summary was right.

**`(pdb, chain, lig)` is not a unique key.** A single chain routinely holds two copies of
the same compound — one in the pocket, one at the peripheral access-channel site that
CYP3A4 is known for. There were 125 such duplicates. Joining on that key cross-joins the
copies and pairs one copy's atoms with the other copy's distances, which is what
manufactured the contradiction.

The key is now `(pdb, chain, lig, seqid)`. Note the shape of the error: **the check was
wrong, not the thing being checked**. That is the same lesson as the metformin validator
false alarm, and it is the second time a verification step has been the buggy half. A
disagreement between two computations is not evidence about which one to trust.

---

## What this does and does not buy

**Does:** a 4.5× larger calibration set; an independent replication of thesis claim B;
corrected acceptance windows that stop discarding 10% of true coordination; a hard
validity filter; and a cached atom store (`p450_atoms.parquet`) so no future geometric
question needs another 585 downloads.

**Does not:** improve pose *selection* by itself. None of this is yet a selection
experiment — it recalibrates the physics terms that selection will use. The n = 87
constraint on *validating a selector* is only lifted if these 496 pairs are actually
co-folded, which needs per-target sequences and MSAs and is the obvious next step.
