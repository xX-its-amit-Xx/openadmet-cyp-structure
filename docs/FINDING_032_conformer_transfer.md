# FINDING 032 — the crystallographic torsion record is tight, and it points away from the answer

**Date:** 2026-09-22 · **Status:** measured, zero new inference, zero new downloads, CPU-only ·
**Verdict:** **REFUTED**, at the pre-registered answer-recognition gate, and refuted with a
mechanism that is the *opposite* of FINDING 030's.

Pre-registered in `docs/PREREG_conformer_transfer.md`, committed at `f9d951e` before a
single number existed. Nothing below was tuned afterwards. Everything post-hoc is labelled
post-hoc.

---

## What was asked

FINDING 030 refuted fragment pose transfer and named exactly what it left standing. It
tested **absolute placement** — where a shared fragment sits in the heme frame — and it
failed because the crystallographic consensus is **4.48 Å RMS wide** while the error to fix
is about 2.5 Å. The prior was twice as wide as the thing it had to resolve.

FINDING 024's decomposition leaves **26% of CYP3A4's ligand error (1.60 Å) in the ligand's
INTERNAL CONFORMER** — its own torsions, independent of where it sits or how it is turned.
That term had never been touched. It is the largest untested term in the decomposition and
it is **frame-free**: a torsion does not care about the pocket, so the 4.48 Å placement
scatter that killed 030 does not apply to it.

> **The question.** Does the crystallographic record constrain a CYP3A4 ligand's internal
> conformation well enough to rank poses?

Two arms:

- **A — donor torsion transfer.** For each query ligand, find P450 crystal donors sharing a
  substructure and compare the **torsion angles** of the shared rotatable bonds. Not
  coordinates, not placement.
- **B — prior-free internal plausibility.** Forget P450s: is the conformation chemically
  reasonable at all, against an ETKDG/CSD-derived conformer ensemble?

Both arms are **intramolecular and protein-free**. FINDING 029 measured clash and strain
against the co-folded pocket at **−0.0632**, anti-selective; nothing here touches a protein
atom, which is the only reason arm B was worth running.

## What was done

`scripts/structure/conformer_transfer.py`. The donor harvest, the six leakage filters, the
MCS and the heme frame are **imported from `fragment_transfer.py`**, not rewritten, so 030
and 032 share one implementation and one donor table.

| | |
|---|---|
| queries | 87 CYP3A4 ligands, `data/processed/validation_ligands.csv` |
| pool | `val87b_unsteered` — 87 × 20 = **1,740** Boltz-2 poses |
| donors | `donors.parquet` — **936 usable** crystal ligand-chain observations, 354 CCD codes, 183 target keys |
| geometry | **raw Cartesian coordinates** on both sides, never heme-frame coordinates: a left-handed basis would silently flip every dihedral's sign. The feature is frame-free by construction |
| torsion | (a,b,c,d) with b–c acyclic and rotatable **in both molecules**, all four atoms inside the MCS, a and d the in-MCS neighbours with the smallest **pattern** index — a numbering both molecules share |
| symmetry | period `360 / lcm(s_b, s_c)` per molecule, **coarser of the two** per pair; `s` counted from `CanonicalRankAtoms(breakTies=False)` symmetry classes. Folded deviation `δ ∈ [0, P/2]`, in degrees |
| features | `tor_cons` (mean over donors — PRIMARY), `tor_best`, `tor_med`, `tor_cmean` (deviation from the donors' **circular mean**), `tor_wgt` (same, weighted by resultant length R) |

### How symmetry was handled, and where the pre-registered rule bites — C-SYM

Torsions are angles and a raw angular difference is meaningless for a terminal phenyl or a
t-butyl, so every deviation is folded to the torsion's own period before anything is
averaged. The period is computed from graph-automorphism classes, on **both** molecules,
and the coarser one wins.

| synthetic case | period expected | measured | |
|---|---|---|---|
| monosubstituted phenyl | 180° | **180°** | pass |
| t-butyl | 120° | **120°** | pass |
| ordinary sp3–sp3 | 360° | **360°** | pass |
| carboxylate `C(=O)[O-]` | 180° | **360°** | **fail** |
| carboxylic acid `C(=O)OH` | 180° | **360°** | **fail** |
| nitro `[N+](=O)[O-]` | 180° | **360°** | **fail** |

`CanonicalRankAtoms` is a graph invariant and is **not resonance aware**: it sees a
carboxylate's two oxygens as `=O` and `[O-]`, which are not automorphic, so those torsions
are folded at 360° where 180° is the physical period. This was pre-registered, wrongly, as
covered. It is reported as a failure rather than patched, because the pre-registration is
fixed once committed. **Exposure is 8 of 1,800 shared torsion keys — 0.44%** (see D4) — and
the direction of the error is to *inflate* a deviation, symmetrically for pose and crystal
alike, so it cannot manufacture the result below.

The relabelling control passes: swapping a phenyl's two ortho carbons as the reference atom
moves the folded deviation by **0.0010°**, which is embedding noise on a real 3D geometry,
not a machinery failure. (The coded flag compares against 1e-6 and therefore reads `false`;
the number, not the flag, is the control.)

---

## Coverage first, as pre-registered

**Coverage is not the problem for arm A either, but it is much thinner than 030's.**

| | ≥ 1 legal donor | ≥ 3 legal donors | median donors | p90 | max | zero |
|---|---|---|---|---|---|---|
| **primary (no CYP3A donors)** | **83 / 87** | **81 / 87** | **143** | 179 | 257 | 4 |
| with-3A (optimistic bound) | 83 / 87 | 82 / 87 | 229 | 265 | 348 | 4 |

86 of 87 queries build a molecule (**MF8** fails atom-count assignment, as in 030). Donors
come from **136 distinct P450 targets**. The four zero-donor queries are **CFF** (caffeine,
**0 rotatable bonds** — structurally uncoverable by any torsion feature), **MF8**, and
**08J**/**08Y**, whose shared fragments are ring systems with no rotatable bond inside them.

The funnel is the interesting part, and it is where 030's 522 donors become 143:

| stage | query–donor cells |
|---|---|
| legal after L0–L5 | **62,764** |
| MCS ≥ 6 heavy atoms | **39,486** |
| MCS matched on both molecules | 39,481 |
| **carries a shared rotatable torsion** | **10,388** (26%) |

A maximum common substructure under `CompleteRingsOnly` is usually a *ring system*, and a
ring system contains no acyclic rotatable bond. Three quarters of the chemistry that 030
could use is invisible to a torsion feature. What survives is thin per donor: **median 1
shared torsion, mean 1.78**, over a median MCS of 8 heavy atoms.

### Every leakage filter fires — C2

Summed over 87 × 936 = **80,496** query–donor cells. These reproduce FINDING 030's table
exactly, which is the cross-check that the imported filters are the same filters.

| filter | rule | cells removed |
|---|---|---|
| L1 | the query's own PDB entry | **108** |
| L2 | the same CCD ligand code, any entry | **61** |
| L3 | count-ECFP4 Tanimoto ≥ 0.90 (fixed a priori) | **153** |
| L5 | `closest_fe` > 10 Å — not an active-site copy | **4,887** |
| L0 | target not classifiable to a UniProt | **258** |
| **L4** | **the whole CYP3A subfamily** (P08684/P20815/P24462/Q9HB55) | **12,265** |

L4 is what makes this leave-one-**TARGET**-out. The with-3A variant lifts only L4.

---

## THE ANSWER-RECOGNITION TEST — and the experiment stops here

FINDING 030's methodological carry-forward, now standard, with a **numeric pass rule fixed
in advance**: mean percentile ≥ 0.65 and binomial p < 0.05 in the correct direction.
Percentile = the fraction of that ligand's 20 predicted poses whose feature value is
**worse** than the crystal's. 1.0 means the crystal is the single most plausible pose in the
set; 0.5 means the feature cannot tell the answer from a wrong guess.

### Arm A — the crystal ranks at the 31st percentile, significantly **below** chance

79 ligands with ≥ 3 legal CYP3A-free donors and a readable crystal ligand.

| feature | mean percentile | median | crystal beats the median pose | binomial p | mean Δ |
|---|---|---|---|---|---|
| **`tor_cons`** (PRIMARY) | **0.308** | 0.15 | **23 / 79** | **0.00026** | +4.17° |
| `tor_best` | 0.389 | 0.35 | 26 / 79 | 0.0032 | −0.30° |
| `tor_med` | 0.347 | 0.20 | 28 / 79 | 0.013 | +4.87° |
| `tor_cmean` | 0.378 | 0.20 | 25 / 79 | 0.0015 | +3.67° |
| `tor_wgt` | 0.386 | 0.25 | 28 / 79 | 0.013 | +4.17° |

**Every one of the five fails the gate, and four of them fail it significantly in the wrong
direction.** This is not 030's chance result (0.558, p = 0.44). The P450 crystallographic
torsion record does not merely fail to recognise the true CYP3A4-bound conformer — it
**prefers the model's wrong poses to it**, reliably, across 79 ligands.

### C1 — scrambled donors, applied to the diagnostic

The pre-registration attaches C1 to a selector. No selector was built, so C1 was applied
where the measurement actually lives: replace every donor's torsion value with a value drawn
at random from the **global pool of donor torsion values** — same machinery, same periods,
same donor count, chemistry removed — and re-measure the crystal's percentile. 8 redraws,
seed 20260922.

| | mean percentile of the crystal | crystal beats median |
|---|---|---|
| **matched donors** | **0.308** | 23 / 79 |
| scrambled donors | **0.410** (sd over redraws 0.026) | 32.6 / 79 |
| chance | 0.50 | 39.5 / 79 |

Two things, and both matter. A *generic* crystallographic torsion angle already disfavours
the crystal (0.410 against 0.50), because the pool of deposited torsion values is
concentrated near the wells the co-folder also prefers. And **matching the substructure
makes it worse, by a further 0.10** — so the substructure correspondence is not decorative
(C1 fires, as it did in 030), it is actively pointing the wrong way.

It is not a near miss at the tail, either. The crystal is the **single worst of the 21** on
**23** ligands and the single best on **6**; **44 of 79** sit below the 25th percentile.

Three post-hoc robustness checks, none of which rescues it:

| variant | ligands | mean percentile | binomial p |
|---|---|---|---|
| primary, no CYP3A donors | 79 | **0.308** | 2.6e-04 |
| **with-3A** (the optimistic bound — the query's own subfamily allowed to donate) | 79 | **0.253** | **9.4e-08** |
| only the tightest-consensus torsions (donor resultant R ≥ 0.8, ≥ 3 donors) | 78 | 0.360 | 9.0e-04 |

**Letting CYP3A4's own subfamily donate makes it worse.** That is the reverse of 030, where
the with-3A bound bought +0.011, and it is the single most leakage-flavoured thing this
experiment can legally do. And it is not dilution by loose torsions: keeping only the
torsions the record agrees on most tightly leaves the crystal at the 36th percentile.

It is also flat across the pool: ρ(percentile, pool oracle) = **−0.06**, ρ with ligand heavy
atoms = −0.07, ρ with rotatable-bond count = −0.08. This is not a catastrophe artefact and
not a big-ligand artefact.

### Arm B — a prior-free small-molecule term says the same thing, harder

82 ligands (ETKDG ensembles, mean 243 conformers after pruning; no P450 data of any kind).

| feature | mean percentile | median | crystal beats the median pose | binomial p | mean Δ |
|---|---|---|---|---|---|
| **`etkdg_tfd_min`** (PRIMARY of arm B) | **0.303** | 0.12 | **24 / 82** | **2.2e-04** | +0.026 TFD |
| `etkdg_tfd_mean` | 0.294 | 0.15 | 23 / 82 | 8.7e-05 | +0.022 TFD |
| `etkdg_well` | 0.258 | 0.05 | 20 / 81 | 5.7e-06 | +4.77° |
| `mmff_strain` (declared **replication**) | 0.162 | 0.05 | 8 / 82 | 1.7e-14 | +21.2 kcal/mol |

In absolute terms: the crystal ligand's torsion-fingerprint distance to its **nearest** of
~243 ETKDG conformers is **0.179**, against **0.154** for the predicted poses, and the
crystal is the worse of the two on **58 of 82**.

**`mmff_strain` must be read with a caveat and is not counted as evidence.** A crystal
ligand's bond lengths and angles come from experimental refinement and a model's come from
an idealised generator, so an MMFF energy is not comparable between them; the +21.8 kcal/mol
is partly that artefact. `etkdg_tfd_min`, `etkdg_tfd_mean` and `etkdg_well` are
**torsion-only** and immune to it, and they say the same thing. The replication itself
stands: FINDING 025 measured MMFF strain as a null selector and nothing here contradicts it.

### The gate

**Neither arm passes. Nine of nine features fail, and all nine fail in the same direction**
— the crystal ranks *below* the middle of its own 20 predictions, six of them at
p < 0.005.

**Within-ligand ranking of the predicted poses is, separately, at chance** — which is the
statistic the prereg requires beside every pooled number, and it is worth seeing next to
the gate:

| feature | within-ligand ρ vs LDDT-PLI | correct sign | binomial p |
|---|---|---|---|
| `tor_cons` (n = 81) | −0.013 | 55.6% | 0.37 |
| `etkdg_tfd_min` (n = 82) | −0.032 | 56.1% | 0.32 |
| `etkdg_well` (n = 81) | −0.084 | 53.1% | 0.66 |
| `mmff_strain` (n = 82) | −0.049 | 59.8% | 0.097 |

So the features are *at chance among the predictions* and *significantly wrong about the
truth*. The truth is out-of-distribution with respect to the pool, which is the whole
content of this finding.

**Arm B coverage, stated plainly:** 86 of 87 ETKDG ensembles completed within the sitting.
**ERY** — erythromycin, a 14-membered macrolide, whose 300-conformer embedding ran for over
an hour without finishing — is excluded and named rather than absorbed. Of the 86, 82 are
usable: MF8 fails the same molecule build as in 030, and three have no readable crystal
ligand for the diagnostic.

**No selector was built**, on either arm, in any variant, weighting or calibration. That is
the pre-registered consequence of failing the gate and it is what FINDING 030 established
the gate for: a prior that ranks the truth below the wrong answers cannot be rescued by the
thing that consumes it. Bars 1, 2 and 3 were therefore not run, and no gain-versus-random or
paired-versus-incumbent number appears in this finding.

---

## Why it fails — five post-hoc measurements

Labelled post-hoc; none is a candidate feature and none may become one without a fresh
pre-registration. `data/processed/conformer_diagnostics.json`.

### D1. **There IS a consensus. That is the surprise.**

This is the direct analogue of FINDING 030's D4, and it comes out the other way round.

| quantity | value |
|---|---|
| shared torsion keys measured | **1,800** (median 4 donors per key) |
| **circular sd of the donors' torsion values, per shared torsion** | **25.0°** (per-ligand median 24.5°) |
| circular sd of one ligand's **own 20 predicted poses** on the same torsions | **35.9°** |
| **folded deviation of a predicted pose from the CRYSTAL, same torsions** | **59.3°** |

In 030 the prior scattered **eight times more widely** than the error it had to resolve, and
that was the whole explanation. Here the prior is **2.4× tighter** than the error to fix and
tighter even than the pool's own spread. **The superfamily does agree about what torsion a
shared rotatable bond takes.** The prior exists, it is sharp, and it is still useless —
because it is sharp about the wrong value.

### D2. **The torsion ORACLE — the ceiling of every possible internal-conformer feature**

The decisive measurement. Score each pose by the folded deviation of its torsions from the
**crystal's own torsions** — the answer itself, handed over completely — and select on it.
Leaky by construction, a diagnostic of the *term*, never a predictor.

| | |
|---|---|
| ligands | 79 |
| random-pose baseline | **0.5863** |
| pool oracle | **0.6993** |
| **selected by the TRUE torsions** | **0.6075** |
| **gain vs random** | **+0.0213** |
| the shipped incumbent `cypstruct.xengine.select()`, same ligands | **0.6258** |
| within-ligand ρ | −0.141, correct sign on **65.8%** |

**Perfect knowledge of the ligand's true internal conformation is worth +0.021 LDDT-PLI, and
loses to the incumbent by 0.018.** For scale — and this is a comparison, not a bar 1 test,
which was not run — FINDING 030 recomputed the random-feature null on this same pool at
**p95 +0.0140 / p99 +0.0193**. The oracle of this entire family of features sits *at* the
p99 of that null.

This is what closes the term, and it closes it far beyond arm A. It does not matter how good
a torsion prior is, where it comes from, or how it is weighted: the internal conformer does
not carry enough of the within-ligand variance in LDDT-PLI to rank these poses. FINDING
024's 26% residual is real as an *error* and is not available as a *signal*.

### D3. Arm A and arm B are not the same feature, and both fail

Within-ligand Spearman between `tor_cons` and each prior-free term, over the ligands where
both exist:

| | n | mean within-ligand ρ | fraction positive |
|---|---|---|---|
| `etkdg_tfd_min` | 78 | **+0.142** | 68% |
| `etkdg_well` | 78 | **+0.108** | 68% |
| `mmff_strain` | 78 | **+0.132** | 69% |

Weakly positive and no more. The P450 torsion record and a generic small-molecule torsion
library are **not** measuring the same thing here, and they nevertheless fail the same gate
in the same direction — which is why the conclusion is about the target rather than about
either prior.

### D4. The C-SYM resonance gap touches almost nothing

| | |
|---|---|
| shared torsion keys | **1,800** |
| keys touching a carboxyl, nitro, sulfonyl or phosphonate group | **8** |
| fraction | **0.44%** |

The pre-registered symmetry rule's one real defect is confined to eight torsions out of
1,800. It cannot account for the result, and a resonance-aware fold would not change it.

### D5. The mechanism, stated plainly

Put D1, the answer-recognition test and C1 together:

- the superfamily's torsion consensus is **sharp** (25.0°);
- the co-folder's poses sit **near** it — near enough that a generic deposited torsion value
  already favours them over the crystal;
- the true CYP3A4-bound conformer sits **59.3° away** from where the poses are, and *further*
  from the consensus than the poses are.

**CYP3A4 binds its ligands in torsionally unusual conformations, and Boltz-2 already gives
them ordinary ones.** A prior that rewards ordinariness therefore rewards exactly the poses
that are wrong. That is a different failure from 030's — there the prior was empty; here the
prior is sharp, correct about the superfamily, and anti-correlated with the target.

It is also consistent with the rest of the file: FINDING 027 found the model's pocket
**excludes** the true pose on 71% of poses, and FINDING 028 found that pocket **rigid**
(0.077 Å of ligand-to-ligand side-chain motion against the crystals' 0.723 Å). A ligand
forced into a pocket that has not opened for it will be pushed to a strained, unusual
conformer. The unusual torsions are not noise in the reference; they are the induced fit.

---

## Controls

| control | required | measured |
|---|---|---|
| **C-NUM** (FINDING 021) | 15 poses rescored from scratch match the shipped truth | **15 / 15, max \|Δ\| = 5.6e-17** |
| **C-NUM**, exact zeros | zero silent numbering failures | **4** unsteered poses score exactly 0.0 — **all four with BiSyRMSD 23–27 Å** (PG0 ×1, PG4 ×3), i.e. the ligand left the protein. Not the 021 signature, which was 0.0 at *small* RMSD |
| **C4** frame identity | this frame reproduces the shipped `xeng` | **1,740 / 1,740 rows, max \|Δ\| = 2.45e-07** |
| **C5** crystal re-run | the diagnostic run reproduces the primary `d` matrix | **max \|Δ\| = 0.0 exactly**, 79 ligands |
| **C2** every filter fires | non-zero counts | **all six**, table above |
| **C3** non-constant feature | > 95% of queries | **100%** for all four headline features; smallest within-ligand sd **0.85°** (`tor_cons`), 0.0022 (`etkdg_tfd_min`) |
| **C-SYM** | synthetic periods | **3 of 6 pass**; the three resonance cases fail and are reported, not patched |
| **C1** scrambled donors | matched must differ from scrambled | **fires**: 0.308 matched vs 0.410 scrambled, 3.9 sd |

---

## Verdict

**REFUTED**, at the pre-registered answer-recognition gate.

Arm A's primary feature puts the query's own crystal ligand at the **31st percentile** of its
own 20 predicted poses (23 of 79 better than the median pose, binomial p = **0.00026**) —
not at chance, as in FINDING 030, but **significantly below it**. Arm B, which uses **no P450 data at all**, puts it at the
**30th** percentile (24 of 82, p = 2.2e-04), so the result does not depend on the
crystallographic record being a P450 record — a general small-molecule torsion library
reaches the same verdict. The gate
was pre-registered with a numeric pass rule and no arm cleared it, so no selector was built
and no bar was run.

The post-hoc mechanism is stronger than the gate. The prior is **not** too wide this time —
donor torsions agree to 25.0° where the error to fix is 59.3°. The prior is sharp and points
the wrong way, and the **torsion oracle** says that even a perfect internal-conformer feature
is worth **+0.0213** against a +0.0193 p99 noise floor and **loses to the incumbent**.

## What it licenses, and what it closes

**Closed — the internal-conformer term, as a source of selection signal, by any prior.**
FINDING 024 left 26% of CYP3A4's ligand error (1.60 Å) in the ligand's own torsions and
called it the largest untouched term in the decomposition. It is now touched, and it is not
a signal: the ceiling of the whole family is +0.021, below the incumbent and at the noise
floor's p99. This is a stronger closure than a failed feature, because it is the *oracle*
that fails. Do not propose another torsion prior, another conformer library, another
strain term, or a fine-tune aimed at conformer quality, for **selection** purposes.

**Closed — "match the crystallographic record" as a family.** Two independent tests, on the
same donor set and with the same filters, now say that agreement with deposited P450 geometry
is worth nothing on this target: **placement** (030, prior 8× too wide) and **torsion** (032,
prior sharp but anti-correlated). The remaining ligand-side priors all share the assumption
that CYP3A4 binds like the superfamily, and 032 is the first direct evidence that it does
not.

**NOT closed.** (a) The internal conformer remains a real component of the *error*; a
generation-side intervention that produced better bound conformers would still improve
LDDT-PLI directly. What is refuted is using conformer quality to **choose among poses**.
(b) Arm B's terms are intramolecular; FINDING 029's −0.0632 for pocket-referenced clash and
strain is untouched and still stands separately. (c) The 59.3° pose-to-crystal torsion
deviation is a measurement worth having on its own: it is the first number this project has
put on *how* wrong the bound conformer is, per torsion, and it is large.

**One methodological carry-forward.** 030 added "hand the feature the answer before you build
a selector". 032 adds the next rung: **hand the answer to the feature's own scoring function
and select with it.** The torsion oracle cost one extra stage and closed a term that the
gate alone would only have wounded — the gate says *this prior* fails, the oracle says *no
prior can succeed*. Both belong in the RUNBOOK.

---

## Artefacts

| what | where |
|---|---|
| script, five stages + diagnostics | `scripts/structure/conformer_transfer.py` |
| pre-registration | `docs/PREREG_conformer_transfer.md` (committed `f9d951e`) |
| answer-recognition gate | `data/processed/conformer_answer_recognition.json` |
| coverage + controls (no selector: gate failed) | `data/processed/conformer_transfer_primary.json` |
| C1 scrambled-donor control | `data/processed/conformer_scrambled_control.json` |
| post-hoc diagnostics D1–D4 | `data/processed/conformer_diagnostics.json` |
| per-ligand table | `data/processed/conformer_per_ligand.csv` |

Intermediates (per-query torsion JSONs, ETKDG ensembles) are in the session scratchpad on
`C:`, never on the full `D:` drive. Nothing was downloaded and no inference was run.
