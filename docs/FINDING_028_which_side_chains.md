# FINDING 028 — which side chains are wrong, which ones exclude the truth, and how far they would have to move

**Date:** 2026-09-22 · **Status:** measured, zero new inference, CPU-only, ~5 minutes on
12 cores · **Verdict:** **DIAGNOSED.** The exclusion is not diffuse and it is not a
rotamer problem in general. **Three residues do 83% of the blocking — Phe215, Arg212,
Phe304 — and six do 97%.** The model's pocket is **rigid, not wrongly adaptive**: its
side chains move **0.077 Å** from one ligand to the next where the crystals move
**0.723 Å**. Where a repack does work it is **one residue turning 30°**; where it fails,
it fails on an atom no torsion can move — **CB or backbone, at 211–216, inside the F/G
loop**.

Nothing here is or may become a selection feature: every quantity is computed against the
crystal. This is diagnosis.

---

## What was asked

FINDING 024 said CYP3A4 fails purely on ligand orientation because the pocket is built to
0.73 Å at CA. FINDING 027 qualified that decisively: the crystal ligand, dropped unmoved
into the model's protein, makes a median closest heavy-atom contact of **1.42 Å** and
clashes below 2.2 Å on **71.2%** of poses, while **0 of 87** deposited complexes clash in
their own protein (minimum 2.22 Å). So the backbone is right and the side chains are
wrong — and 31 of 41 failures are unreachable by any rigid ligand rotation.

Four questions, four numbers. Which residues are wrong (Q1). Which ones do the excluding
(Q2) — a different question. Is the pocket rigid or adaptive-in-the-wrong-direction (Q3).
And how much motion would be enough (Q4).

## What was done

`scripts/structure/side_chain_diagnosis.py`, four stages, reusing FINDING 027's frame,
its reference loading and its `renumber_to_reference`, so the numbers are comparable to
027 by construction rather than by assertion.

| | |
|---|---|
| pool | `val87b_unsteered` — 87 ligands × 20 Boltz-2 samples = **1,740 poses** |
| frame (Q2) | pocket-CA superposition, residues within 8 Å of the crystal ligand, matched by residue **number** — identical to FINDING 027 |
| frame (Q0/Q1/Q3) | `CYP3A4_RIGID_CORE` from `targets.py`, a **ligand-free** set declared long before this experiment, so the fit cannot be dragged by whichever pocket residue moved |
| core pocket | residues contacting the crystal ligand within 6.0 Å in **≥ 20 of 87** pairs → **32 residues** |
| chi grid (Q4) | chi1 × chi2 at **10°** (1,296 points/residue), escalated to chi1 × chi2 × chi3 at 10°/10°/**20°** (23,328) only when chi1 × chi2 fails |
| acceptance (Q4) | ligand contact ≥ **2.2 Å** *and* repacked atoms ≥ **2.6 Å** from the rest of the protein (excluding i ± 1) and from the heme, Fe excluded |

Every threshold above, and the expected **sign** of every correlation in Q3, is in the
`PREREG` dict at the top of the script and is echoed into
`data/processed/side_chain_diagnosis.json`. Signs were fixed before any correlation was
computed. `loo-sign-selection-fakes-negatives` is the reason that sentence is here.

### Chi symmetry, done properly

- **chi2 wrapped modulo 180** for Phe, Tyr, Asp. The brief said "Phe/Tyr/Asp/Glu chi2";
  that is right for three of the four and **wrong for Glu**, whose degenerate torsion is
  chi3 (CB–CG–CD–OE1), not chi2 (CA–CB–CG–CD). Implemented correctly; Glu chi2 is treated
  as non-degenerate. Arg chi5 is also in the degenerate set and is never computed here.
- **Val CG1/CG2 and Leu CD1/CD2 are canonicalised geometrically**, by a fixed sign of
  `(CG1−CB) × (CG2−CB) · (CA−CB)`, applied identically to the crystal and the prediction.
  This removes a spurious ~120° chi1 error from a naming-convention swap without
  collapsing the three genuine rotamers, which is what "minimise over the swap" would do.
  It **fires**: 560 swaps across the 87 crystals (6.4 per structure) and 11,240 across the
  1,740 predictions (6.5 per pose) — the same rate, which is the point.
- Asn OD1/ND2 and Gln OE1/NE2 are **not** swapped. They are chemically distinct, however
  often they are misassigned in deposition.
- Heavy-atom RMSDs minimise over Asp/Glu/Phe/Tyr/Arg/Val/Leu label swaps.

### Controls — every one of them fired

| control | required | measured |
|---|---|---|
| **FINDING 021 numbering** | offset 0, residue-name identity high | offset **0 on all 1,740**, min identity **0.9956** |
| **ligand atom mapping** | no silent index-for-index | 1,740 / 1,740 symmetry-mapped, **0 poses dropped** |
| **reproduces 027** | same measurement, same number | median **1.4096 Å** vs 027's 1.42; **71.38%** below cut vs 71.2%; control median **2.7025 Å**, min **2.2228 Å**, **0.0%** below cut |
| **frame-independence** | 027's pocket frame is not doing the work | rigid-core frame: median **1.396 Å**, **71.09%** below cut |
| **crystal side chains are real** | occupancy, completeness, B | 1,965 contact-residue observations: **2 incomplete** (0.10%), **0 zero-occupancy**, 18 below occupancy 1.0, 18 with altlocs |
| **chi rows** | the completeness filter fires | 92,940 rows, **140 dropped** for an incomplete crystal side chain, 92,800 used; 15,540 have no chi1 (Ala/Gly) |
| **self-clash cutoff not over-strict** | the model's own side chains satisfy it | median **3.25 Å** to the rest of the protein, p5 2.44 Å, **8.4%** below 2.6 Å — and it bound on only **3** residue scans out of ~1,400 |

**The crystal-quality caveat that matters.** Median resolution **2.65 Å** (worst 3.4 Å)
and median pocket side-chain **B = 81.8 Å²** (p90 129 Å²). The side chains are modelled —
they are not modelled *well*. **Arg212 carries a median B of 137 Å²** and is incomplete in
1 of the 29 crystals that contact it. It is the second-largest blocker below, and it is
the conclusion in this document with the least crystallographic support. Phe215, Phe304
and Phe213 sit at 60–107 Å² and are better determined.

---

## Q0 — the backbone, split. This was not asked and it changes the reading

FINDING 024's "pocket-lining CA 0.73 Å" is an average over a region that is two different
regions. Superposed on the **ligand-free** rigid core:

| | median over 1,740 poses |
|---|---|
| rigid-core CA (the fit itself) | **0.47 Å** |
| whole-chain CA, RMS | 0.71 Å (per-residue median 0.34 Å) |
| pocket CA, RMS | 0.65 Å |
| **pocket CA, residues NOT in the F/G span** | **0.44 Å** — and **0.0%** of poses exceed 2 Å |
| **pocket CA, residues IN the F/G span (202–260)** | **1.19 Å**, p90 **4.61 Å** — **39.0%** of poses exceed 2 Å |

A median pocket has **54** residues of which **8** are in the F/G span. **The backbone is
built to 0.44 Å everywhere in the pocket except the F/G loop, where two poses in five are
more than 2 Å out.** OpenADMET named that span as the remodelling hotspot; FINDING 024
measured it as complete in the crystals and concluded it was not what limits us. Both are
true. It is modelled, and it is modelled in the wrong place.

---

## Q1 — which residues are wrong, and by how much

Pooled over the 32 core-pocket residues, per (pair, residue), after collapsing the 20
samples with a circular mean:

| | |
|---|---|
| median \|Δchi1\| | **8.0°** (IQR 3.3–19.0) |
| fraction of (pair, residue) with \|Δchi1\| > 40° | **12.6%** |
| median \|Δchi2\| | 9.1°; 10.2% over 40° |
| median side-chain RMSD, own-backbone frame | **0.39 Å** |
| median side-chain RMSD, rigid-core frame | 0.53 Å |
| median CA error | 0.34 Å |
| per pair: contact residues / of which chi1 off > 40° | median **22** / median **2** (11.1%) |

**The error is concentrated, not spread.** The five worst residues carry **56%** of the
summed `frac_chi1_off` across all 32; the median residue is at **5.0%** and the worst at
**83.3%**. Four residues are at exactly 0.0% (Met114, Arg105, Ile300, Cys442).

| rank | residue | n pairs | median \|Δchi1\| | IQR | frac > 40° | sc-RMSD (local) | CA err | F/G |
|---|---|---|---|---|---|---|---|---|
| 1 | **Leu210** | 24 | **99.8°** | 62–113 | **0.833** | 2.10 Å | 2.61 Å | yes |
| 2 | **Leu211** | 22 | 61.9° | 17–106 | 0.636 | 2.02 Å | **4.25 Å** | yes |
| 3 | **Phe108** | 79 | 38.2° | 23–59 | 0.456 | 1.04 Å | 0.81 Å | no |
| 4 | Leu482 | 47 | 21.7° | 7–88 | 0.404 | 0.97 Å | 0.69 Å | no |
| 5 | Phe213 | 44 | 14.2° | 6–90 | 0.386 | 0.71 Å | 0.65 Å | yes |
| 6 | Ile301 | 86 | 12.7° | 5–48 | 0.314 | 0.40 Å | 0.31 Å | no |
| 7 | **Arg212** | 28 | 17.4° | 9–54 | 0.286 | **2.58 Å** | 0.58 Å | yes |
| 8 | Phe304 | 86 | 9.1° | 4–41 | 0.256 | 0.61 Å | 0.59 Å | no |
| … | Thr224, Arg106, Phe57, Pro107 | | 7–19° | | 0.13–0.24 | | | |
| ↓ | Ile369, Thr309, Phe241, Ser119, Met371, Glu374 | | 5–8° | | 0.04–0.09 | | | |
| = 0 | Met114, Arg105, Ile300, Cys442 | | 3–5° | | **0.000** | | | |

Split by region, the same gradient as Q0:

| | median \|Δchi1\| | frac > 40° | sc-RMSD local | n |
|---|---|---|---|---|
| F/G span | **11.9°** | **26.3%** | **0.72 Å** | 331 |
| everywhere else in the pocket | 7.3° | 9.8% | 0.35 Å | 1,632 |

Mann-Whitney **p = 2.9e-12**.

**Answer to Q1: concentrated, and concentrated on the F/G loop.** Name them: Leu210,
Leu211, Arg212, Phe213 — consecutive, all inside 202–260 — plus Phe108 and Leu482 outside
it. Two of those four carry a CA error of 2.6 and 4.2 Å, so their "rotamer" error is
partly a backbone error wearing a rotamer's clothes.

---

## Q2 — which residues cause the exclusion. A **different** list

Crystal ligand placed in the model's protein, worst contact attributed to the atom that
makes it:

| | |
|---|---|
| worst contact is a **side-chain** atom | **89.8%** of poses |
| worst contact is the heme | 4.9% |
| worst contact lies in the F/G span | **66.1%** |

Among the 1,242 poses that clash below 2.2 Å, the residue making the worst contact:

| residue | poses worst | ligands | median worst contact | share | cumulative |
|---|---|---|---|---|---|
| **Phe215** | 445 | 31 | 0.78 Å | **35.8%** | 35.8% |
| **Arg212** | 412 | 33 | 1.13 Å | **33.2%** | **69.0%** |
| **Phe304** | 169 | 32 | 0.91 Å | 13.6% | **82.6%** |
| Phe241 | 71 | 6 | 0.50 Å | 5.7% | 88.3% |
| Phe213 | 67 | 9 | 0.74 Å | 5.4% | 93.7% |
| Ser119 | 44 | 5 | 1.94 Å | 3.5% | **97.3%** |
| Leu211, Asp214, Arg106, Phe108, Leu216, … | ≤10 | | | ≤0.8% | 100% |

**Three residues do 83% of it. Six do 97%.**

### The cross-tabulation, which is the point of asking twice

Across the 32 core-pocket residues, ρ(frac chi1 off > 40°, frac blocking) = **+0.337,
p = 0.059** — a weak positive that does not clear significance at n = 32. Wrongness and
blocking are **largely dissociated**, and the dissociation is legible residue by residue:

| | frac chi1 off > 40° | median \|Δchi1\| | frac of poses it blocks |
|---|---|---|---|
| **wrong, harmless** — Leu210 | **0.833** | 99.8° | **0.000** |
| Leu211 | 0.636 | 61.9° | 0.001 |
| Leu482 | 0.404 | 21.7° | 0.000 |
| Phe108 | 0.456 | 38.2° | 0.023 |
| **blocking, barely moves** — **Phe215** | **0.024** | **8.7°** | **0.195** |
| Ser119 | 0.058 | 7.6° | 0.081 |
| Phe241 | 0.083 | 7.8° | 0.087 |
| Phe304 | 0.256 | 9.1° | 0.222 |
| **both** — **Arg212** | 0.286 | 17.4° (sc-RMSD **2.58 Å**) | **0.361** |

**Phe215 is the single largest blocker in the pool and its chi1 is off by a median of
8.7°.** It is not in the wrong rotamer. It is in the right rotamer at the wrong place —
its CB is carried there by the F/G backbone. Leu210 is turned 100° the wrong way and
never touches the ligand.

**This kills "get the Phe rotamers right" as a slogan.** The eight pocket phenylalanines
do not behave alike: Phe108 is the most rotamerically wrong Phe and blocks 2.3% of poses;
Phe215 is among the least wrong and blocks 19.5%. Any restraint or scorer built on the Phe
roof has to name *which* Phe and *for which purpose*.

---

## Q3 — rigid, or wrongly adaptive? **Rigid.** Not close

Circular sd of chi1 **across the 87 different ligands**, per residue, model against
crystal:

| | median over 53 residues |
|---|---|
| **model** | **0.63°** |
| **crystal** | **20.63°** |
| ratio | **0.047** |
| residues where the model is more variable | **0 of 53** |
| Wilcoxon | **p = 2.4e-10** |
| model, *within* a ligand (across its 20 samples) | 0.75° |

Residue by residue, the pattern is uniform, and the biggest blockers are among the most
frozen:

| residue | crystal sd chi1 | model sd chi1 |
|---|---|---|
| Leu210 | 74.3° | 12.4° |
| Ile223 | 57.3° | **1.2°** |
| Arg212 | 57.1° | 20.0° |
| Leu482 | 55.3° | **1.4°** |
| Arg106 | 43.2° | **1.1°** |
| Val240 | 41.7° | **0.3°** |
| Phe304 | 38.9° | 12.2° |
| **Phe215** | 19.4° | **2.0°** |
| Phe57 | 36.4° | 1.3° |

### The same answer with no torsions in it

Because a chi-angle result can always be doubted on symmetry handling, the whole thing was
repeated in Cartesian space: superpose residue *r* of structure *i* onto residue *r* of
structure *j* by **that residue's own N/CA/C**, then take the side-chain heavy-atom RMSD,
over all 3,741 pairs of different ligands. Same operation for the 87 predictions and for
the 87 crystals.

| | median over 27 residues | p90 |
|---|---|---|
| **model** | **0.077 Å** | 0.55 Å |
| **crystal** | **0.723 Å** | 1.82 Å |
| crystal more variable | **25 of 27** | |
| model more variable | **2 of 27** — Phe108 (3.71 vs 3.10 Å) and Pro107 (0.67 vs 0.38 Å) | |
| Wilcoxon | **p = 5.5e-6** | |

Two residues go the other way, which is what a real comparison looks like rather than a
filter that never fires.

**Diagnosis (a): the model collapses to one generic rotamer set and does not induce fit at
all.** It is not adapting in the wrong direction — it is not adapting. CYP3A4's pocket is
famously plastic; the crystals move their pocket side chains by 0.72 Å from ligand to
ligand and the predictions move them by 0.08 Å. That is a **9.4× deficit in plasticity**,
and it is the mechanism behind 027's "the model builds a pocket its own answer fits": the
pocket it builds is very nearly the *same* pocket every time, and the ligand is then
placed to suit it.

### Does pocket side-chain error predict pose error? Only within a ligand

Signs fixed a priori: negative against LDDT-PLI, positive against BiSyRMSD.

| feature | between-ligand ρ (LDDT-PLI) | within-ligand mean ρ (LDDT-PLI) | correct sign | Wilcoxon |
|---|---|---|---|---|
| mean \|Δchi1\| over contacts | −0.211 | −0.168 | 72.4% | 2.8e-05 |
| mean side-chain RMSD, local frame | −0.270 | −0.184 | 72.4% | 2.6e-06 |
| **mean side-chain RMSD, core frame** | **−0.082** | **−0.308** | **87.4%** | **4.1e-12** |
| mean CA error over contacts | −0.098 | −0.245 | 81.6% | 1.5e-10 |

Against BiSyRMSD the strongest is the same feature: between +0.131 within-ligand (70.1%,
p = 2.1e-04) against −0.006 between-ligand.

**This is `within-ligand-variance-is-what-matters` in its cleanest form to date.** The
between-ligand column reproduces FINDING 024's ρ(pocket CA, ligand RMSD) = +0.03 null: how
good a ligand's protein is tells you almost nothing about how good its pose is. The
within-ligand column is strong and consistent: among a ligand's own 20 samples, the ones
whose pocket side chains sit closer to the crystal score better, on 87% of ligands.

**And it is not a signal.** It needs the crystal, and more importantly it is a *consistency
relation between two things that were moulded together* — the co-folder places the side
chains around the ligand it chose, so a pose nearer the truth drags its side chains with
it. It is evidence about the mechanism, not a feature. It is logged here so that it is
never proposed as one.

---

## Q4 — how much side-chain motion would be enough

Scope: the **31** rotation-unreachable failures from FINDING 027 × 20 poses = **620**
poses. The 10 rotation-rescuable failures and the 46 already-solved pairs were scanned
too, as strata, so the machinery can be seen to behave differently on easy and hard cases
rather than failing everywhere.

**What was searched:** for every residue with any atom under 2.2 Å of the placed crystal
ligand, a Cartesian grid over chi1 × chi2 at 10° resolution (36 × 36 = 1,296 rotamers),
escalated to chi1 × chi2 × chi3 at 10°/10°/20° (23,328) only when the two-torsion scan
failed. A rotamer is accepted when the residue's own side chain clears the ligand to
≥ 2.2 Å **and** stays ≥ 2.6 Å from the rest of the protein (excluding i ± 1, with the other
clashing residues treated as absent) and from the heme. Every clashing residue is then set
to its minimum-perturbation accepted rotamer simultaneously and the whole-pocket contact is
recomputed. Perturbation magnitude is `max |Δchi|` over the scanned torsions.

### Before the repack

| stratum | ligands | poses | median start contact | frac clashing | median clashing residues | frac with a backbone clash |
|---|---|---|---|---|---|---|
| **rotation-unreachable failures** | 31 | 620 | **0.87 Å** | **80.6%** | **2** (IQR 1–3, max 7) | **21.6%** |
| rotation-rescuable failures | 10 | 200 | 1.22 Å | 75.5% | 1 | 6.0% |
| already-solved control | 46 | 920 | 1.84 Å | 64.2% | 1 | 6.7% |

The gradient is monotone in all four columns and in the right direction, which is the
internal control on the whole stage.

### After the repack

| stratum | frac of poses cleared to ≥ 2.2 Å | median max \|Δchi\| on cleared poses |
|---|---|---|
| rotation-unreachable failures | **58.2%** (361 / 620) | **30°** |
| rotation-rescuable failures | 93.5% | 20° |
| already-solved control | 82.1% | 10° |

On the 31 hard ligands: **24 of 31** have at least one pose the scan can open, **13 of 31**
have all twenty, and **7 of 31 have none**.

### The magnitude — it is one residue, turning 30°

Over the 361 cleared poses:

| | |
|---|---|
| median residues moved | **1** (IQR 0–2) |
| need ≤ 1 residue | **63.7%** |
| need ≤ 2 residues | 78.7% |
| needed **no** move (they were not clashing) | 33.2% |
| median max \|Δchi\| | **30°** (IQR 0–70°) |
| max \|Δchi\| ≤ 30° | 51.5% |
| max \|Δchi\| ≤ 60° | 72.3% |
| median summed \|Δchi\| | 40° |

Restricted to the 241 poses that actually needed a move: **45.6% need exactly one
residue**, and the max-\|Δchi\| deciles are 20° / 30° / 60° / 80° / 90° (p10/p25/p50/p75/p90).

Self-consistency of the concerted solution: only **5 of 361** put two repacked side chains
within 2.6 Å of each other. Dropping those leaves 356 and changes the per-ligand count not
at all.

**It is not a concerted repack. It is one Phe or one Arg turning 30–60°.**

### Where the scan fails, it fails on an atom no torsion can move

Tally over every residue the scan could not clear, across all 620 poses:

| reason | count |
|---|---|
| **the residue's own CB is inside the ligand** | **303** |
| the rest of the protein blocks every clearing rotamer | 3 |
| no clearing rotamer for another reason | **0** |
| no rotatable axis | **0** |

**Every single failure is a CB or backbone collision.** CB's position is fixed by the
backbone, so this is a backbone error, not a rotamer error. The residues:

| residue | poses | ligands | median CB-to-ligand distance | in F/G span |
|---|---|---|---|---|
| **Phe215** | 77 | 4 | **1.38 Å** | yes |
| **Arg212** | 77 | 5 | **1.42 Å** | yes |
| Phe304 | 62 | 6 | 2.08 Å | no |
| Phe213 | 48 | 3 | 1.10 Å | yes |
| Leu211 | 21 | 3 | 1.11 Å | yes |
| Leu216 | 18 | 1 | 2.06 Å | yes |

### The 7 pairs no side-chain scan can rescue

| ligand | pocket CA RMSD | F/G-span CA RMSD | min backbone-to-ligand contact | worst backbone atom |
|---|---|---|---|---|
| MWY | 0.78 Å | 1.48 Å | **0.25 Å** | 212:N |
| A1A4T | 0.94 Å | 1.49 Å | **0.67 Å** | 214:C |
| ERY | 1.62 Å | **2.04 Å** | **0.82 Å** | 213:CA |
| 1RD | 1.53 Å | 1.94 Å | **0.80 Å** | 216:N |
| 5AW | 1.41 Å | 1.77 Å | **1.09 Å** | 213:CA |
| X7P | 0.59 Å | 0.83 Å | 2.53 Å | 215:CA — **CB of Phe215 at 1.44 Å** |
| QEP | 0.29 Å | 0.52 Å | 2.98 Å | 212:CA — **CB of Arg212 at 1.64 Å** |

Five of the seven have the crystal ligand passing straight through **main-chain** atoms of
residues **212–216**. The other two have a clean backbone and a CB sitting 1.4–1.6 Å inside
the ligand, which is the same statement one atom further out. **All seven are the same
lesion: the 210–216 stretch of the F/G loop is in the wrong place.** Note that QEP has the
*best* pocket CA RMSD in the whole table (0.29 Å) and is still unrescuable — a
whole-pocket CA average hides a five-residue displacement completely.

---

## Verdict

**DIAGNOSED.** Four answers, all with numbers:

1. **Which residues are wrong:** concentrated, not spread. Leu210 (83% of pairs off by
   > 40°), Leu211 (64%), Phe108 (46%), Leu482 (40%), Phe213 (39%), Ile301 (31%),
   Arg212 (29%), Phe304 (26%). Top 5 carry 56% of the total; four residues are at 0%.
   F/G-span residues are 2.7× more likely to be off by > 40° (26.3% vs 9.8%, p = 2.9e-12).
2. **Which residues exclude the truth:** a *different* list. Phe215 (36%), Arg212 (33%),
   Phe304 (14%) — 83% between them, 97% with Phe241, Phe213 and Ser119 added.
   ρ(wrongness, blocking) = +0.34, p = 0.06: the two are largely independent.
   **Phe215 blocks most and barely moves; Leu210 moves 100° and blocks nothing.**
3. **Rigid or wrongly adaptive:** **rigid**. 0.077 Å of ligand-to-ligand side-chain motion
   against the crystals' 0.723 Å; chi1 sd 0.63° against 20.63°; 0 of 53 residues more
   variable than reality. No induced fit at all.
4. **How much motion is enough:** for 58% of the hardest poses, **one residue turning 30°**
   (median; 64% need ≤ 1 residue, 52% need ≤ 30°). For the other 42% — and for 7 of the 31
   ligands outright — **no amount**, because the obstruction is CB or main chain at
   211–216.

---

## What this licenses, and what it forecloses

### LICENSED — repack six named residues, then re-search orientation. The mechanism is measured

FINDING 027 froze the protein and found that only **2.06%** of orientations survive a
steric filter, around the wrong answer. This finding says why and says what to unfreeze:
**97.3% of the blocking is done by Arg212, Phe215, Phe304, Phe241, Phe213 and Ser119**,
and the median repair is **one of them turning 30°**. The concrete experiment is 027's
pipeline with one line changed — before rotating the ligand, enumerate rotamers of those
six residues on a 10–30° chi grid, accept on the same pre-registered 2.2 Å contact and
2.6 Å self-consistency rules, and re-run the orientation scan inside each repacked pocket.
It is CPU-only and needs no new inference, exactly as 027 was.

**With the prior attached, because it is not optimistic.** 027 is the direct precedent and
it went: oracle **+0.0108**, selection **−0.0145**. FINDING 013 and 016 went the same way.
The honest expectation is that this raises the ceiling and that
`cypstruct.xengine.select()` does not reach it. Run it as a **ceiling measurement**, report
the oracle *and* the selection, in that order, and pre-register both — the value is in
learning whether the true pose is reachable at all without new inference, not in a score.

### LICENSED, and it reverses a call in FINDING 024 — a holo template is aimed at exactly the F/G lesion

024 dismissed holo templating partly because "a template carries no ligand and **no side
chains** — only backbone frames and CB", and argued that the Phe rotamers are what would be
needed. **This finding inverts that for the cases that matter most.** Of the residues whose
blocking cannot be repacked away, every one fails on **CB or main chain** (303 of 306
uncleared residues are "CB inside the ligand", 0 are "no clearing rotamer"), and 7 of the
31 hardest ligands are unrescuable for exactly that reason. **Backbone frames and CB are
precisely what a template carries, and precisely what those 7 pairs need.**

024's two other objections still stand and must be respected: the proxy set is circular
unless the template is forbidden from the query's scaffold cluster, and CYP3A4 crystals
differ from each other by 1.00 Å of pocket CA against the model's 0.73 Å. But that 1.00 Å
is a **whole-pocket** number, and QEP shows a whole-pocket number hiding a fatal
five-residue displacement. **The measurement that decides it is the crystal-to-crystal
spread of CA over residues 210–216 alone, against the model's 1.19 Å median / 4.61 Å p90 in
that span.** If the crystals agree with each other there better than the model agrees with
any of them, an F/G-restricted template — or an F/G backbone ensemble — is the first
protein-side intervention on this target with a measured mechanism behind it. That is one
cheap CPU measurement and it is the next thing to run.

### FORECLOSED — treating the co-folded protein as a receptor

Rigid ligand search inside it (027), classical rigid docking into it, and any physics
scorer that assumes it. The pocket does not move: 0.077 Å from ligand to ligand. It is a
single generic conformation that the model then fits its own ligand into, and FINDING 025's
34 physics terms were all scored against that conformation. **A receptor with no induced
fit cannot be the frame for a physics score on a protein whose defining property is induced
fit.**

### FORECLOSED — side-chain accuracy as a selection feature

Its between-ligand ρ with LDDT-PLI is −0.08. Its within-ligand ρ is −0.31 at 87% correct
sign, and that is a consistency relation between a ligand and the side chains that were
moulded around it, not information about the truth. It also requires the crystal. Logged so
it is not proposed.

### FORECLOSED — "fix the Phe roof" as a single intervention

The eight pocket phenylalanines are not one object. Phe108 is the most rotamerically wrong
and blocks 2.3% of poses; Phe215 is among the least wrong and blocks 19.5%. Any restraint,
loss term or scorer naming "the Phe cluster" is averaging over residues that fail in
opposite ways.

---

## Method notes worth keeping

- **Reproducing 027 to three decimals (1.4096 vs 1.42; 71.38% vs 71.2%; control min 2.2228)
  before measuring anything new** is what makes every attribution here an attribution of
  *that* measurement rather than of a lookalike. It cost one function call, because the
  frame and the loaders were imported from 027's script rather than rewritten.
- **Asking "which residue is wrong" and "which residue blocks" as two questions was the
  whole value of the exercise.** One question would have produced the Leu210 list and aimed
  the next experiment at a residue that never touches the ligand.
- **A ligand-free alignment frame changes the story.** A pocket-CA fit is dragged by
  whichever pocket residue moved; on 08J it reported 3.19 Å of pocket CA RMSD where a
  rigid-core fit reports 0.63 Å core / 1.08 Å non-F/G pocket / 8.88 Å F/G pocket. The
  decomposition only exists in the ligand-free frame.
- **A whole-pocket average hid a five-residue lesion.** QEP: pocket CA 0.29 Å, F/G CA
  0.52 Å, and an Arg212 CB 1.64 Å inside the crystal ligand. Report the span, not the mean.
- **The prochiral canonicalisation fired at the same rate in both structures** (6.4 vs 6.5
  per structure). Had it fired on one side only, every Val chi1 and Leu chi2 in this
  document would have been ~120° of pure convention.
- **The Glu chi2 correction.** The brief asked for Glu chi2 modulo 180; Glu's degenerate
  torsion is chi3. Implementing what was asked would have wrapped a non-degenerate angle
  and silently halved every Glu chi2 deviation.

---

# Addendum, 2026-09-22 — the template ceiling over the lesion. **Refuted there too**

The body of this finding said holo templates come back on the table only if one deposited
CYP3A4 crystal predicts another's **210–216** better than the model does, and named that as
the deciding measurement. It has now been run. **It does not.** FINDING 024's dismissal of
holo templating stands unqualified, and the partial reversal proposed above is withdrawn.

**Decision rule, fixed before the numbers were read** (recorded in the `templates` stage's
docstring and echoed into the output JSON): *licensed only if the holo–holo
different-ligand spread over 210–216 is clearly below the model's **median and p90** error
over the same span; comparable or worse refutes it.*

Run: `python scripts/structure/side_chain_diagnosis.py templates` →
`data/processed/template_ceiling_lesion.json`, `template_ceiling_pairs.csv`. 107 deposited
CYP3A4 (P08684) entries, all already in `data/reference/rcsb`, no new downloads. One chain
per entry — two chains of one crystal are not two opinions.

## First, as instructed: is the span even modelled?

| | |
|---|---|
| deposited CYP3A4 entries used | **107** |
| **all seven of 210–216 modelled** | **54 (50.5%)** |
| 4–6 of 7 modelled | 10 |
| **none of the span modelled** | **17** |
| median CA B-factor inside the span | **104.7 Å²** (against 81.8 Å² for pocket side chains generally) |
| minimum CA occupancy inside the span | 1.0 |

Per residue: 210 **72.9%**, 211 58.9%, **212 53.3%**, 213 57.0%, 214 60.7%, 215 63.6%,
216 68.2%. **Arg212 and Phe213 — the second- and fifth-largest blockers — are missing from
nearly half of all deposited CYP3A4 structures.**

That is already most of the answer. **A template cannot carry what the crystal does not
contain**, and on a coin-flip of entries it contains nothing here. It also means every
spread below is computed on the *better half* of the archive and is therefore optimistic.

The same audit on the 87 validation crystals: **48 of 87 model all seven**, 58 model at
least four. The F/G placement of 29 of 87 pairs cannot be scored at all.

## Independence of the entries

| | |
|---|---|
| distinct space groups | 6 — but **92 of 107 are `I 2 2 2`** |
| others | `C 1 2 1` 11, and one each of `C 2 2 21`, `P 31`, `P 21 21 21`, `I 1 2 1` |
| resolution | median **2.55 Å**, range 1.78–3.10 Å |
| distinct ligand-code sets | 98 of 107 |
| **apo entries in this harvest** | **0** |

Two things follow. The archive is **one crystal form with a long tail**, so the spread
below is *not* inflated by crystal-packing diversity — if anything it is suppressed by the
lack of it. And **there are no apo entries here at all**: 1TQN and 1W0E are not in the
P450-universe harvest, so FINDING 024's "the only always-legal template is apo" could not
be tested and the **apo–apo** and **apo–holo** strata below are empty rather than small.

## The three stratified spreads

CA deviation over residues 210–216 after superposing two entries on `CYP3A4_RIGID_CORE`
(ligand-free; median core fit **0.56 Å**), over all 1,431 pairs of the 54 fully-modelled
entries:

| stratum | pairs | median | IQR | **p90** | max | > 1 Å | > 2 Å |
|---|---|---|---|---|---|---|---|
| apo–apo | **0** | — | — | — | — | — | — |
| apo–holo | **0** | — | — | — | — | — | — |
| holo–holo, **same** ligand | 7 | 0.78 Å | 0.48–2.63 | 4.51 Å | 4.79 Å | 28.6% | 28.6% |
| **holo–holo, DIFFERENT ligand** | **1,424** | **4.25 Å** | 0.70–4.74 | **7.20 Å** | 9.76 Å | **64.7%** | **58.4%** |
| … of those, cross-space-group only | 532 | 4.42 Å | 4.24–4.74 | 7.40 Å | 9.76 Å | 94.2% | 88.2% |
| all pairs | 1,431 | 4.25 Å | 0.70–4.74 | 7.20 Å | 9.76 Å | 64.6% | 58.2% |

**Against the model, measured over the same seven residues in the same frame:**

| | median | p90 | > 1 Å | > 2 Å |
|---|---|---|---|---|
| **model vs its own query crystal** (1,160 poses, 58 pairs) | **1.03 Å** | **4.98 Å** | 50.3% | 45.2% |
| **crystal vs a different-ligand crystal** (1,424 pairs) | **4.25 Å** | **7.20 Å** | 64.7% | 58.4% |

**A blind deposited template is 4.1× worse at the median and 1.4× worse at the p90 than
what Boltz-2 already produces, at exactly the residues that do the blocking.** This is
FINDING 024's whole-pocket argument (model 0.73 Å against crystal-to-crystal 1.00 Å)
repeated in the lesion, where it is far more extreme, not less.

### The spread is real, not a numbering artifact

The distribution is strongly bimodal, so it was checked rather than reported:

- **Numbering is identical across entries** — every entry examined reads
  Leu210-Leu211-Arg212-Phe213-Asp214-Phe215-Leu216. This is FINDING 021's trap and it is
  not what is happening here.
- **In the same superposition, a control span moves normally.** For pairs showing 3–12 Å
  in the lesion, the I-helix (300–320) agrees to **0.4–1.3 Å**: 5TE8 vs 1W0G, lesion
  per-residue 3.5 / 3.1 / 6.9 / 5.5 / 5.3 / 8.9 / **11.8 Å**, I-helix **0.60 Å**; 7UFA vs
  9BV5, lesion 7.1 / 8.5 / 8.5 / 7.2 / 6.3 / 5.9 / 4.9 Å, I-helix **0.60 Å**. The rest of
  the protein superposes; this loop does not.
- The bimodality is **discrete loop states**, not a few outliers. Scoring each entry by its
  median deviation to every other entry: **27 of 54 form a tight majority cluster** (within
  it, median **0.53 Å**, p90 **0.78 Å**), **22 are minority conformers ≥ 3 Å** from the
  median other entry (3NXU, 4D78, 4D7D, 4I4G, 4I4H, 4K9T, 4K9W, 5TE8, 5VC0, 6UNJ, 7KVH,
  7KVK, 7KVN, 7KVO, 7KVQ, 7KVS, 7UFA, 7UFE, 7UFF, 8EWS, 8EXB, 8SPD), and 5 are intermediate.

**The one number that could have rescued templates is 0.53 Å — and it is unavailable.**
Inside the majority cluster a crystal predicts another crystal's F/G loop to 0.53 Å, twice
as well as the model. But that cluster is **half** the fully-modelled entries and **a
quarter of the archive**, and choosing it requires knowing which state the query is in,
which is the prediction. Blind, you draw from the whole set, and the whole set is 4.25 Å.

## The seven unrescuable ligands specifically

For each, every fully-modelled entry was superposed onto the query crystal by rigid core
and asked whether **its** 210–216 backbone + CB would clear the query's crystal ligand to
2.2 Å:

| ligand | its PDB | own span fully modelled? | donors clearing | best donor | best contact |
|---|---|---|---|---|---|
| **5AW** | 4K9U | no | **0 / 54** | 5VCE | 0.99 Å |
| **ERY** | 2J0D | no | **0 / 54** | 8EXB | 2.12 Å |
| A1A4T | 9COY | no | 3 / 54 (5.6%) | 8EWS | 4.72 Å |
| **1RD** | 4K9T | **yes — a minority conformer, 4.51 Å from the median entry** | 5 / 53 (9.4%) | 8EXB | 4.53 Å |
| MWY | 6OOA | no | 6 / 54 (11.1%) | 7UFA | 2.75 Å |
| X7P | 7KVI | no | 21 / 54 (38.9%) | 8EXB | 7.17 Å |
| QEP | 6UNG | no | 22 / 54 (40.7%) | 8EWS | 5.88 Å |

For two of the seven, **not one deposited CYP3A4 structure** has an F/G loop that would
admit their ligand. For three more, fewer than 12% would. And **six of the seven do not
model 210–216 in their own deposition**, so even a hypothetical oracle template has no
ground truth to be graded against on those cases; the one that does is in the minority
state that only 22 of 54 entries share.

## Verdict

**REFUTED in the lesion too.** A blind holo template over 210–216 is worse than the model
at the median (4.25 Å against 1.03 Å) and at the p90 (7.20 Å against 4.98 Å), is absent
from half the archive, and would fail to admit the true ligand for the specific pairs that
need it. The pre-registered rule required "clearly below" on both statistics; the
measurement came back 4× above on one and 1.4× above on the other.

**FINDING 024's dismissal of holo templates stands unqualified**, and the partial reversal
this document proposed is withdrawn. The reversal was reasoned from the right mechanism —
the lesion needs backbone and CB, which is exactly what a template carries — and was wrong
on the empirical premise that any crystal knows where this loop goes. **None of them do.**

**What survives.** The F/G loop is a discrete multi-state element with at least two
populated conformations 4–8 Å apart, half the archive in one of them, and the model
predicting it to 1.03 Å median. That is not a template problem, it is an **ensemble**
problem: the useful object is not "the right crystal" but the *set* of observed 210–216
conformations, used as alternative receptors to be scored rather than as a prior to
condition on. That is a different experiment with a different failure mode, and it inherits
FINDING 027's warning in full — 22 minority conformers × a repack of the six named blockers
is a far larger pool, and every pool expansion on this target so far has raised the oracle
and lowered selection. It should be run, if at all, as a **ceiling** measurement with both
numbers pre-registered.

**Method note worth keeping.** The bimodality was the tell. A median of 4.25 Å with an IQR
of 0.70–4.74 Å is not a spread, it is two populations, and reporting the median alone would
have described a distribution that no pair of crystals actually occupies.

---

# Addendum 2, 2026-09-22 — the APO template, the caveat the first addendum left open. **Refuted like the holo ones**

The first addendum refuted **holo** templates over 210–216 but could not test apo ones:
the 107-entry harvest contained **zero** apo CYP3A4 structures, so FINDING 024's separate
claim — that the only always-legal blind template is an **apo** structure — was recorded
as untested rather than refuted. That caveat was published on the dashboard and is now
closed. **It closes the same way.**

**Decision rule, fixed before any number was read** (in the script's docstring and echoed
into the output JSON as `prereg`): *an apo template is LICENSED only if its CA deviation
over 210–216 against the 87 validation crystals is below **the model** — median < 1.03 Å
**and** p90 < 4.98 Å. The bar is the model, not the holo templates. Beating the 4.25 Å
holo–holo spread is **not** sufficient and is not a win.*

Run: `python scripts/structure/apo_template_ceiling.py` →
`data/processed/apo_template_ceiling.json`, `apo_entry_classification.json`,
`apo_template_pairs.csv`, `apo_clearance.csv`, `apo_model_span_error.csv`. The span, the
rigid-core frame, the per-entry chain reader, the crystal loader, FINDING 021's
renumbering and the 2.2 Å clearance cutoff are **imported** from
`side_chain_diagnosis.py`, not rewritten, so these numbers are commensurable with the
first addendum by construction.

## Finding the apo entries

A live RCSB query for UniProt **P08684** returns **122** entries — 15 more than the
107-entry harvest, which is why the first addendum saw none. All 122 mmCIFs were
**already on disk** in `data/reference/rcsb`; **0 new downloads** were needed and nothing
was deleted.

**Rule used:** an entry is apo when its only non-water heteroatom components are **HEM**
and members of `cypstruct.targets.IGNORE_HET` (buffers and cryoprotectants). To stop a
cryoprotectant hiding in the pocket, the closest non-HEM het atom to the heme Fe is
reported for every entry alongside the call.

| | |
|---|---|
| deposited P08684 entries | **122** |
| **apo** | **6** — 1TQN, 1W0E, 4I3Q, 5VCC, 5VCD, 9GK1 |
| holo | 116 |
| all six outside the 107-entry harvest | yes — which confirms the first addendum's "0 apo" |
| space group | **all six `I 2 2 2`**, the archive's dominant form |
| resolution | 1.70 / 1.95 / 2.05 / 2.60 / 2.80 / 2.95 Å |

**The two calls worth seeing.** **9YK4** (1.78 Å) looks ligand-free by chemical-component
name but carries **imidazole 2.15 Å from the Fe** — an axial ligand, catalogued in
`IGNORE_HET` — plus tetraethylene glycol at 5.07 Å. Its pocket is occupied; it is **not**
apo. **5VCC and 5VCD** are counted apo and carry **glycerol 4.41 / 4.58 Å from the Fe**;
the distance is reported rather than argued away, so a reader can discount those two. The
conclusion below does not depend on either call — all six apo entries give the same answer
to within 0.37 Å.

## First, as instructed: do the apo entries model the span?

**All six model all seven of 210–216.** This is the opposite of the holo archive, where
only 54 of 107 do.

| entry | res. | 210–216 modelled | CA B-factors 210→216 (Å²) | median CA B | min CA occ | names |
|---|---|---|---|---|---|---|
| **5VCC** | 1.70 | **7 / 7** | 35.6 38.2 36.5 33.0 30.8 28.0 30.7 | **33.0** | 1.0 | ✓ |
| **1TQN** | 2.05 | **7 / 7** | 41.9 41.4 38.7 34.1 32.9 31.8 29.8 | **34.1** | 1.0 | ✓ |
| **1W0E** | 2.80 | **7 / 7** | 47.5 45.9 45.9 37.5 35.9 37.1 42.1 | **42.1** | 1.0 | ✓ |
| **5VCD** | 1.95 | **7 / 7** | 50.4 48.7 55.0 49.7 47.2 51.4 42.6 | **49.7** | 1.0 | ✓ |
| **4I3Q** | 2.60 | **7 / 7** | 96.1 99.2 92.8 74.0 72.9 75.3 71.9 | **75.3** | 1.0 | ✓ |
| **9GK1** | 2.95 | **7 / 7** | 126.9 122.6 117.2 108.6 106.5 102.7 105.7 | **108.6** | 1.0 | ✓ |

Against the holo archive's **50.5% fully modelled, median CA B 104.7 Å²**, the apo
entries are *better* determined in the lesion: 6 of 6 complete, four of six at median
B 33–50 Å². **So "the crystal does not contain it" — which was most of the answer for
holo templates — is not available here.** The apo entries do contain the loop. The
question is forced onto the geometry, which is where it should be.

**Numbering trap, checked before any geometry was trusted.** Every one of the six reads
**Leu210-Leu211-Arg212-Phe213-Asp214-Phe215-Leu216**, and so do all **87** validation
crystals: **0 mismatches in 93 structures**. FINDING 021's trap is not what is happening
here.

## Q1 — how well does an apo structure predict the holo F/G span?

Each donor superposed onto each of the 87 validation crystals on `CYP3A4_RIGID_CORE`; CA
deviation over 210–216 alone, in that superposition. **Median core fit 0.556 Å**, matching
the first addendum's 0.56 Å.

**Every filter, with counts:**

| | |
|---|---|
| donors | **6 apo** + **61 holo** with the span fully modelled |
| pairs attempted | 5,782 |
| dropped, query models < 4 of the span | **1,943** (29 of 87 queries) |
| dropped, < 50 shared rigid-core CAs | **0** |
| pairs used | **3,839** — 348 apo, 3,491 holo |
| validation crystals modelling all seven | **48 of 87** |

| donor set | n | median | IQR | **p90** | max | > 1 Å | > 2 Å |
|---|---|---|---|---|---|---|---|
| **apo → holo query** | **348** | **1.08 Å** | 0.53–4.62 | **5.15 Å** | 7.42 Å | 51.1% | 46.6% |
| holo → holo query, *same pairing* | 3,491 | 4.25 Å | 0.70–4.69 | 7.07 Å | 10.33 Å | 66.0% | 59.7% |
| **the MODEL, recomputed here** | **1,160** | **1.028 Å** | 0.47–4.35 | **4.977 Å** | 7.43 Å | **50.3%** | **45.2%** |

The model row is not quoted from the first addendum; it is **recomputed in this script**,
on the same 58 scorable queries, and lands on **1.0277 Å / 4.977 Å / 50.3% / 45.2%**
against that addendum's **1.03 / 4.98 / 50.3% / 45.2%**. That reproduction is what makes
the comparison below a comparison rather than a lookalike.

**Against the rule fixed in advance, apo fails both legs.** 1.084 Å > 1.03 Å on the
median; 5.154 Å > 4.98 Å on the p90. It is a near-tie, not a win, and the rule was written
as "below the model on both" precisely so a near-tie could not be talked into one.

### It is not one unlucky donor, and no donor rescues it

| donor | median | p90 | > 2 Å |
|---|---|---|---|
| 4I3Q | **0.963 Å** ✓ | 5.010 Å ✗ | 46.55% |
| 1W0E | 1.043 Å ✗ | 5.113 Å ✗ | 46.55% |
| 5VCD | 1.120 Å ✗ | 4.9805 Å ✗ | 46.55% |
| 5VCC | 1.164 Å ✗ | 4.995 Å ✗ | 46.55% |
| 9GK1 | 1.163 Å ✗ | 5.042 Å ✗ | 46.55% |
| 1TQN | 1.199 Å ✗ | 5.113 Å ✗ | 46.55% |
| **oracle: best apo donor per query** | **0.901 Å** | **4.966 Å** | 46.55% |

**Not one of the six is below the model on both statistics.** 4I3Q clears the median and
misses the p90; 5VCD ties the p90 to four decimals and misses the median. Picking 4I3Q
*after* seeing this table is `loo-sign-selection-fakes-negatives` run forwards — it is
recorded, and it is not a licence. Even the **oracle** — the best of the six chosen per
query, which is not available blind — lands at 0.901 Å / 4.966 Å, i.e. **level with the
model, not above it.**

### The head-to-head, per query

Paired, on the 58 queries both can be scored on:

| | |
|---|---|
| queries where **the model** is closer to the truth | **46 of 58** |
| queries where **apo** is closer | **12 of 58** |
| median (apo − model) | **+0.131 Å**, in the model's favour |
| binomial p | **8.2e-06** |
| Wilcoxon p | **4.9e-08** |

The model is better than an apo template on this span, and significantly so.

### Why — the mechanism: **the six apo entries are one conformation**

Every donor returned `frac_over_2A` of exactly **0.4655**. A number that clean is the tell,
so it was checked rather than reported.

| | |
|---|---|
| apo–apo pairs (15) over 210–216 | median **0.330 Å**, max **0.477 Å**, **0%** over 1 Å |
| I-helix control on those pairs | 0.28 Å |
| per query, range across the six donors | median **0.21 Å**, max **0.37 Å** |
| queries where the six agree within 0.5 Å | **58 of 58** |

And the 58 scorable queries then split in two, not into a spread:

| | |
|---|---|
| queries where **every** apo donor is under 1 Å | **27** |
| queries where **every** apo donor is over 2 Å | **27** (4.3–7.3 Å) |
| intermediate | 4 |

0.4655 is 27/58. The identical fractions are **explained** — the apo entries are a single
loop conformer, so they succeed and fail on the same queries — not a filter that never
fires.

**This is the answer in one sentence.** The F/G loop has at least two populated states
4–8 Å apart, roughly half the archive in each, and **all six apo structures sit in the
same one**. On the 27 queries that share that state an apo template is excellent
(0.3–0.9 Å, better than the model); on the 27 that do not it is 4.3–7.3 Å wrong, and
nothing about being apo tells you which you are facing. **A single conformation cannot
cover a two-state loop**, and crystallising without a ligand does not give you the other
state — it picks one.

## Controls — every one fired, with counts

| control | required | measured |
|---|---|---|
| **numbering (FINDING 021)** | L210-L211-R212-F213-D214-F215-L216 | **0 mismatches in 93 structures** (6 apo + 87 validation) |
| **I-helix 300–320, same superposition** | 0.4–1.3 Å even where the lesion is 3–12 Å | all 3,839 pairs: median **0.605 Å**, max 1.63 Å, **0%** over 2 Å |
| … on the 2,148 pairs with lesion > 3 Å | " | median **0.884 Å**, **86.6%** inside 0.4–1.3 Å |
| … on the **150 apo** pairs with lesion > 3 Å | " | median **1.08 Å**, max **1.40 Å**, **89.3%** inside |
| … on the model's 1,160 rows | " | median **0.40 Å**, max 1.06 Å |
| **reproduces the model bar** | 1.03 / 4.98 from the first addendum | **1.0277 / 4.977**, 50.3% / 45.2% |
| **reproduces the holo–holo spread** in a new pairing | ≈ 4.25 Å | **4.25 Å** median, 3,491 pairs |
| **core-CA sufficiency** | no pair scored on a thin frame | **0** pairs dropped; median fit 0.556 Å |
| **own-crystal clearance** | the truth must clear its own ligand | 48 queries, median **4.11 Å**, min **2.55 Å**, **0%** below 2.2 Å |

The named examples, in the form the first addendum used: **1W0E vs 8SPD** — lesion
**7.42 Å**, I-helix **0.52 Å**, core fit 0.50 Å. **1W0E vs 5TE8** — lesion **7.31 Å**,
I-helix **0.54 Å**. The rest of the protein superposes; this loop does not. The
superposition is right and the biology is the discrepancy.

The control is not vacuous in the other direction either: **13.4%** of the > 3 Å pairs put
the I-helix slightly above 1.3 Å (max 1.63 Å), all of them holo pairs with a poorer core
fit. A window nothing ever violated would be the thing to distrust.

## Q2 — would an apo 210–216 backbone + CB admit the crystal ligand?

Same a-priori test as the first addendum: superpose the donor on the rigid core, graft its
210–216 **backbone + CB** into the query frame, and require ≥ **2.2 Å** to the query's
crystal ligand.

| | |
|---|---|
| grafts tested (6 donors × 87 queries) | **522** |
| grafts clearing 2.2 Å | **406 (77.8%)**, median min contact 3.27 Å |
| ligands with ≥ 1 clearing apo donor | **70 of 87** |
| **ligands with NO clearing apo donor** | **17 of 87** |
| ligands where all six clear | 64 |
| **control — the query's OWN span vs its own ligand** (48 fully modelled) | median **4.11 Å**, min **2.55 Å**, **0%** below the cut |

The control is the point: a *correct* 210–216 always admits the ligand it was solved with,
so 2.2 Å is calibrated rather than arbitrary. Split by which state the query is in:

| query's loop state | ligands | grafts clearing | ligands with no clearing donor |
|---|---|---|---|
| in the apo state (apo span dev < 1 Å) | 27 | 92.6% | 2 |
| in a different state (> 2 Å) | 27 | 83.3% | 4 |
| span not scorable (< 4 residues) | 29 | **55.7%** | **11** |

### The seven unrescuable ligands: **0 of 42**

| ligand | its PDB | own span modelled | apo donors clearing | best donor | best contact | own crystal |
|---|---|---|---|---|---|---|
| **QEP** | 6UNG | 4 / 7 | **0 / 6** | 4I3Q | **2.18 Å** | 3.71 Å |
| **X7P** | 7KVI | 1 / 7 | **0 / 6** | 5VCC | 1.04 Å | 8.85 Å |
| **ERY** | 2J0D | 4 / 7 | **0 / 6** | 1TQN | 0.95 Å | 4.59 Å |
| **5AW** | 4K9U | 1 / 7 | **0 / 6** | 9GK1 | 0.87 Å | 8.52 Å |
| **1RD** | 4K9T | **7 / 7** | **0 / 6** | 5VCC | 0.70 Å | 2.55 Å |
| **A1A4T** | 9COY | 0 / 7 | **0 / 6** | 9GK1 | 0.62 Å | — |
| **MWY** | 6OOA | 1 / 7 | **0 / 6** | 5VCC | 0.42 Å | 3.32 Å |

**Not one of the 42 apo grafts admits the ligand of any of the seven.** This is *worse*
than the holo donors on exactly the cases that need it: the first addendum found holo
donors clearing for five of the seven (up to 40.7% of donors, best contacts 4.5–7.2 Å),
while apo clears for **none**, with a best contact of 2.18 Å that misses the cut.

Two rows carry the mechanism. **1RD** models all seven of its own span and is a minority
conformer — the apo donors are **4.58 Å** away from it and clear its ligand by 0.70 Å.
**QEP** is the opposite, and more instructive: the apo span reproduces QEP's crystal
backbone to **0.64 Å**, the best agreement among the seven — and still fails, at 2.18 Å
where QEP's own crystal gives 3.71 Å. **A sub-ångström backbone difference is enough to
fail the clearance test**, which is the quantitative content of "the lesion is a CB and
main-chain problem": it has no tolerance to spend.

## Verdict

**REFUTED, like the holo templates, and for a different reason.** Holo templates failed
because half the archive does not model 210–216 and those that do disagree by 4.25 Å.
**Apo templates fail although all six model the span, at better B-factors than the holo
archive.** They fail because:

1. Against the bar fixed in advance they are **1.084 Å / 5.154 Å** against the model's
   **1.028 Å / 4.977 Å** — worse on both legs, and worse on **46 of 58** queries
   individually (p = 4.9e-08).
2. All six occupy **one** loop conformation (apo–apo 0.33 Å median, 0.48 Å max). They are
   excellent on the 27 queries in that state and 4.3–7.3 Å wrong on the 27 that are not,
   and being apo carries no signal about which is which.
3. On the seven ligands the repack scan cannot rescue, **0 of 42 apo grafts** admit the
   crystal ligand — worse than holo donors on the same seven.

**FINDING 024's "the only always-legal blind template is an apo structure" is now tested
and refuted.** It is always *legal*; it is not *informative*. The dashboard caveat is
closed: templating over the F/G lesion is dead in both flavours, and this document's
withdrawn partial reversal of 024 stays withdrawn, now without the untested qualifier.

**What survives is unchanged, and slightly sharpened.** The useful object is still the
**ensemble** of observed 210–216 conformations used as alternative receptors to be scored,
not a prior to condition on — and the apo entries are worth exactly **one** member of that
ensemble, not six, because they are one conformation. The first addendum's warning applies
in full: every pool expansion on this target so far has raised the oracle and lowered
selection, so it is a ceiling measurement with both numbers pre-registered, or nothing.

**Method note worth keeping.** *Six independent donors returning `frac_over_2A = 0.4655`
to four decimals was the tell, and chasing it produced the finding.* The suspicion was a
stuck filter; the reality was that the six apo structures are one conformer and therefore
fail on the identical 27 queries. `too-clean-numbers-are-the-tell` held — but what it
pointed at was the mechanism, not a bug, and the only way to know which was to measure the
apo–apo spread and the per-query range instead of reporting the table.
