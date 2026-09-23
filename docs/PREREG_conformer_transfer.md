# PRE-REGISTRATION — the ligand's INTERNAL CONFORMER as a scoring prior

**Written:** 2026-09-22, before a single score exists. Committed before
`scripts/structure/conformer_transfer.py` computed anything. **Nothing below may be edited
after the first number.** Post-hoc measurements, if any, are labelled post-hoc in the
write-up and may not be promoted to features without a fresh pre-registration.

---

## Why this, and why now

FINDING 030 refuted fragment pose **transfer** and named what it left standing. It tested
**absolute placement only** — where a shared fragment sits in the heme frame — and it failed
with a mechanism: the crystallographic consensus is **4.48 Å RMS wide** while the error to
fix is about **2.5 Å**. The prior was twice as wide as the thing it had to resolve.

FINDING 024's error decomposition leaves **26% of CYP3A4's ligand error in the ligand's own
internal conformer (1.60 Å against the family's 0.34 Å)**. That term has never been touched.
It is the largest untested term in the decomposition and it is **frame-free**: a torsion does
not care where the pocket is or how the ligand is turned, so the placement scatter that
killed 030 does not apply to it.

**The question.** Does the crystallographic record constrain a CYP3A4 ligand's internal
conformation well enough to rank poses?

Two arms, deliberately different in kind:

- **Arm A — donor torsion transfer.** Does the P450 superfamily's crystallographic record
  say what torsion a shared rotatable bond should take?
- **Arm B — prior-free internal plausibility.** Forget P450s: is this ligand conformation
  chemically reasonable at all, judged against a general small-molecule knowledge base?

Arm B is cheaper, simpler and prior-free, and is run as a full arm rather than a footnote
because it may be the better feature.

**Explicitly out of scope.** FINDING 029 measured clash and strain terms against the
**co-folded pocket** at **−0.0632** — anti-selective. Everything here is **intramolecular
and protein-free**. If any quantity in this work requires a protein atom, it is out of
scope and the experiment stops; that is a closed experiment being re-run.

---

## Data, fixed

| | |
|---|---|
| queries | 87 CYP3A4 ligands, `data/processed/validation_ligands.csv` |
| pool | `val87b_unsteered` — 87 × 20 = 1,740 Boltz-2 poses, unchanged from 030 |
| truth | `data/processed/poses_scored_val87b.csv`, `arm == "unsteered"`, LDDT-PLI |
| incumbent | the shipped `xeng` column, `data/processed/xeng_val87b.csv` |
| donors | `donors.parquet` as built by `scripts/structure/fragment_transfer.py donors` — 936 usable crystal ligand-chain observations, 354 CCD codes, 183 target keys. **Rebuilt by the same code, not edited** |
| seed | 20260922 throughout |

Torsions are computed from **raw Cartesian coordinates** — the pose's model coordinates and
the donor's deposited coordinates — never from heme-frame coordinates, because a frame with
a left-handed basis would silently flip every dihedral's sign. This is the point of the
experiment: the feature is frame-free.

---

## Arm A — donor torsion transfer

### A.1 Which donors are legal

Identical to FINDING 030, filter for filter, with the counts reported per filter because a
filter that never fires is a bug, not a pass.

| filter | rule |
|---|---|
| L0 | donor target not classifiable to a UniProt |
| L1 | the query's own PDB entry |
| L2 | the same CCD ligand code, any entry |
| L3 | count-ECFP4 (Morgan r=2) Tanimoto ≥ **0.90** to the query |
| L5 | `closest_fe` > **10 Å** — not an active-site copy |
| **L4** | **the entire CYP3A subfamily** — P08684, P20815, P24462, Q9HB55 |

**Primary analysis excludes L4** (leave-one-TARGET-out). The with-3A variant lifts only L4
and is reported **separately, as the optimistic bound it is**, never as a headline.

### A.2 Which torsions are shared

MCS parameters identical to 030's primary arm: `CompareElements`, `CompareOrderExact`,
`RingMatchesRingOnly`, `CompleteRingsOnly`, timeout 10 s, **≥ 6 heavy atoms**, one MCS per
(query, donor CCD code).

A **shared torsion** is an ordered quadruple (a, b, c, d) such that:

1. the bond b–c is **single, acyclic**, and matches the standard rotatable-bond pattern
   `[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]` (amide bonds are **kept**; terminal-atom bonds and
   triple bonds are excluded by the pattern itself);
2. **all four** of a, b, c, d lie inside the MCS match, so the quadruple corresponds
   one-to-one between query and donor;
3. a is a neighbour of b other than c, d is a neighbour of c other than b;
4. a and d are chosen **canonically** as the in-MCS neighbour with the smallest index in
   the MCS **pattern** numbering — a numbering shared by both molecules, so the same
   physical torsion is measured on both sides.

### A.3 Symmetry — how it is handled

Torsions are angles, and terminal phenyls, carboxylates, nitros and t-butyls make a raw
angular difference meaningless. Handling, fixed a priori:

- Symmetry classes come from `Chem.CanonicalRankAtoms(mol, breakTies=False,
  includeChirality=False)`, which assigns equal rank exactly to graph-automorphic atoms —
  the two ortho carbons of a phenyl, the two oxygens of a carboxylate, the three methyls of
  a t-butyl.
- For a torsion (a,b,c,d) in a molecule: `s_b` = the number of neighbours of b other than c
  that share a's symmetry class (including a itself); `s_c` likewise for d about c.
- The **period** of that torsion in that molecule is `360 / lcm(s_b, s_c)` degrees.
- The period used for a (query, donor) pair is the **coarser of the two molecules'** periods,
  i.e. `360 / lcm(s_b^q, s_c^q, s_b^D, s_c^D)`. If either side is symmetric, the fold applies.
- Folded deviation: `δ = ((Δθ mod P) + P) mod P`, then `δ = min(δ, P − δ)`, so
  `δ ∈ [0, P/2]`. A phenyl torsion therefore lives in [0°, 90°], a t-butyl in [0°, 60°],
  an ordinary torsion in [0°, 180°].

Deviations are reported and aggregated in **degrees**, unnormalised, so that a torsion whose
period is small cannot dominate by having a wide range. (Normalising by P/2 is the obvious
alternative; it is *not* used, and this sentence is here so that choosing it later would be
visibly post-hoc.)

**Control C-SYM.** Synthetic verification that the machinery returns P = 180° for a
monosubstituted phenyl and for a carboxylate, P = 120° for a t-butyl, and P = 360° for an
ordinary sp3–sp3 torsion; and that relabelling a donor's symmetry-equivalent atoms leaves δ
unchanged. Reported in the write-up whether it passes or fails.

### A.4 The features, with signs fixed from chemistry

Per donor *i*: `d_i` = **mean over that pair's shared torsions of δ, in degrees**.

| name | definition | sign |
|---|---|---|
| **`tor_cons`** — **PRIMARY** | mean of `d_i` over legal donors | **lower is better** |
| `tor_best` | min of `d_i` | lower is better |
| `tor_med` | median of `d_i` | lower is better |
| `tor_cmean` | per shared torsion, δ between the pose's angle and the **circular mean** of that torsion's donor values (circular mean taken on the P-folded angle, i.e. mean of `exp(i·2π·θ/P)`, argument back-transformed); averaged over torsions | lower is better |
| `tor_wgt` | `tor_cmean` with each torsion weighted by the donors' resultant length R (concentration), so a torsion the record agrees on counts more | lower is better |

The sign is fixed from chemistry in every case and identically: **a pose whose torsions sit
closer to what crystallography puts there is the more plausible pose.** No sign is chosen
per fold, per arm, or after seeing a number. (FINDING 031's lesson, and the standing memory
note: per-fold sign selection manufactures p = 9e-09 from noise.)

**Inclusion.** A query enters the primary subset if it has **≥ 3 legal donors each carrying
≥ 1 shared torsion**. Queries with zero rotatable torsions (rigid ligands — caffeine and its
kind) are **structurally uncoverable by arm A** and are reported as such in coverage, not
quietly dropped.

---

## Arm B — prior-free internal plausibility

No P450 data of any kind. No protein atom of any kind.

**Ensemble.** For each query, from its `validation_ligands.csv` SMILES: `EmbedMultipleConfs`
with ETKDGv3, **K = 300**, `randomSeed = 20260922`, `pruneRmsThresh = 0.5`, then MMFF94s
minimisation (500 steps, failures kept as embedded). ETKDG's torsion terms are derived from
small-molecule crystallography (CSD torsion library), so this ensemble is precisely the
"populated rotamer well" knowledge base the question asks about, and it is available for any
molecule with no lookup.

| name | definition | sign |
|---|---|---|
| **`etkdg_tfd_min`** — **PRIMARY of arm B** | min over the ensemble of the **Torsion Fingerprint Deviation** between the pose and the conformer (`rdkit.Chem.TorsionFingerprints.CalculateTFDBetweenConformers`, which performs its own symmetry handling) | **lower is better** |
| `etkdg_tfd_mean` | mean over the ensemble | lower is better |
| `etkdg_well` | mean over the pose's rotatable torsions of the folded circular distance to the **nearest** value that torsion takes anywhere in the ensemble — "how far is this torsion from a populated well" | lower is better |
| `mmff_strain` | MMFF94s energy of the pose (hydrogens added and relaxed, heavy atoms fixed) minus the energy after unconstrained minimisation | lower is better |

`mmff_strain` is declared **in advance as a replication**, not a candidate: FINDING 025
already measured an MMFF94s local-strain term at **+0.0105 against a +0.0138 floor with a
negative within-ligand ρ**. It is computed so the new arm has a calibrated neighbour, and a
positive result on it would be a red flag about this pipeline rather than a discovery.

Atom-order mapping between the pose molecule and the ensemble molecule is by
**substructure match**, never index-for-index. (`cyp-pose-atom-mapping-trap`: index-for-index
once halved every LDDT-PLI in this project.)

---

## THE ANSWER-RECOGNITION TEST — run first, and it can stop the experiment

This is FINDING 030's methodological carry-forward, now standard. Before any selector is
built:

> Score the query's **own crystal ligand** with the feature, as if it were one more pose,
> and report where it ranks among that ligand's 20 predictions.

Crystal ligand coordinates come from `donors.parquet` (the query's own entry — which L1
excludes from every scored quantity), in pose atom order **transferred by substructure
match**. Arm A's donor set for the diagnostic is the primary, CYP3A-free set, exactly as in
the primary run; control **C5** requires the diagnostic run to reproduce the primary run's
`d` matrix on the 20 predicted poses to < 1e-5.

**Percentile** = the fraction of that ligand's 20 predicted poses whose feature value is
**worse** (higher) than the crystal's. 1.0 means the crystal is the single most plausible
thing in the set; 0.5 means the feature cannot tell the answer from a wrong guess.

**Pass rule, fixed now:**

> An arm proceeds to bar 1 **only if** its mean percentile is **≥ 0.65** **and** the count of
> ligands where the crystal beats the median predicted pose is significant at **binomial
> p < 0.05** in the correct direction.

If an arm fails, it is **not** carried to a selector, on any variant, weighting or
calibration — that is what 030 established the test for. If **both** arms fail, the
experiment stops there and is written up as a closed negative.

The crystal is a **diagnostic only**. It never enters a donor set, a null, a selector, a
z-score or any reported selection number.

---

## The bars, in order, and only after the answer-recognition test passes

**Bar 1 — beat random.** Gain over the random-pose baseline must exceed the **p99 of a
2,000-draw random-feature null recomputed for these features on this pool and this subset**.
No figure from FINDING 007 or 030 is quoted as the floor.

**Bar 2 — beat the incumbent, paired.** Against `cypstruct.xengine.select()` on identical
poses and identical ligands: mean Δ, **bootstrap 95% CI (10,000 draws)**, **Wilcoxon**, and
the **tie count** (ligands where both selectors pick the same pose in every tie-break draw).
Eight successive candidates have cleared bar 1 and died here; bar 2 is the bar.

**Bar 3 — complementarity.** Pearson r between per-ligand Δ and the incumbent's headroom
(oracle − incumbent selection), reported with the standing caveat that Δ and headroom share
the term −(incumbent selection) by construction, so the number is comparable across findings
and is not evidence on its own.

**Reported with every pooled mean, without exception:** the **pool oracle**, the random
baseline, the **within-ligand Spearman ρ**, and the **correct-sign fraction** with its
binomial p. Ties broken at random over **64 draws**.

---

## Controls, each of which must report a number

| id | requirement |
|---|---|
| **C-NUM** | FINDING 021. Recompute `lddt_pli` from scratch for 15 randomly drawn unsteered poses and match `poses_scored_val87b.csv` to < 1e-6; and **zero** unsteered poses score exactly 0.0 |
| **C4** | the heme frame recomputed here reproduces the shipped `xeng` column to < 1e-6 over all available rows |
| **C2** | every leakage filter L0–L5 reports a **non-zero** count, summed over queries |
| **C3** | the feature is non-constant within a ligand for > 95% of queries; the smallest within-ligand sd is reported |
| **C-SYM** | the symmetry-folding machinery passes the synthetic cases in §A.3 |
| **C1** | scrambled donors: each matched donor's torsion values replaced by values drawn from a random legal donor's torsions, same count, same aggregation, 8 redraws. Matched must beat scrambled, or the substructure match is decorative |

**Too-clean numbers are the tell.** An sd of exactly 0, a filter that never fires, or a
control that passes to more digits than the arithmetic supports is treated as a bug and
chased before anything is reported.

---

## What each outcome licenses

| outcome | reading |
|---|---|
| answer-recognition fails in both arms | **REFUTED.** The internal conformer is not a rankable term at this pool's resolution, by either a P450 prior or a general one. Closes 024's 26% residual to scoring, and with 030 closes the ligand-side prior programme |
| answer-recognition passes, bar 1 fails | **MEASURED.** The feature knows the truth when handed it but cannot find it among 20 guesses — a discrimination failure rather than a knowledge failure, and a different diagnosis from 030's |
| bar 1 passes, bar 2 fails | the ninth candidate to die at bar 2. Logged, with the tie count, as an incumbent-absorption result |
| both bars pass | a candidate for shipping, subject to bar 3 and to leave-one-TARGET-out |

**Artefacts, fixed now:** script `scripts/structure/conformer_transfer.py`; results
`data/processed/conformer_*.json` / `.csv`; write-up
`docs/FINDING_032_conformer_transfer.md`. Intermediates go to the session scratchpad on
`C:`, never to `D:`. No downloads, no inference, CPU only.
