# FINDING 025 — The CYP3A4 physics scorer, built above the iron, is a null

**Date:** 2026-09-22 · **Check:** 1,740 Boltz-2 poses, 87 CYP3A4 ligands, crystal truth,
leave-one-ligand-out · **Verdict:** every term null; the 34-feature ensemble lands on its
own null's observed maximum and loses to a parameter-free incumbent. **Does not ship.**

---

## What was built, and what was deliberately not

`QM_SCORER_DESIGN.md` sections A, B and C are all first-coordination-sphere terms. They are
implemented in `qmscore/geometry.py` and they are inert — five independent confirmations
(FINDINGs 008, 010, 014, 019, 021) say the iron anchor is already reproduced by unbonded
co-folding. They were not rebuilt.

What was built, in `src/cypstruct/qmscore/pocket.py`, is everything the design assigns to
the region *above* the heme, plus what the evolutionary decomposition supplies that the
design predates:

| block | terms | provenance |
|---|---|---|
| **polar** (15) | Ser119, Arg106, Arg212, Asp214 wells + distances; **Thr224**; backbone amide/carbonyl; buried-unsatisfied polar count | design §E; Thr224 and the Gly481-backbone motif from `worldmodel/CYP3A4_EVOLUTION.md` §4.4 |
| **phe** (9) | ligand aromatic ring vs **Phe57/108/213/215/219/220/241/304** — centroid distance, interplanar angle, lateral offset, parallel-displaced vs T-shaped counts, Phe304 gatekeeper standoff | absent from the design entirely; `CYP3A4_EVOLUTION.md` §4.3 |
| **fg** (6) | F/G roof (202–260) contact count, side-chain-only count, share of all contacts, min distance, residues engaged, ligand centroid height over the porphyrin | "discrimination lives under the F/G roof" |
| **strain** (4) | MMFF94s local strain (pose with H relaxed − nearest minimum), protein clash depth and count, intramolecular clash | design §D, included to log the PXR replication |

Ligand rings are perceived from the pose's own coordinates (planarity test), not read off
the input SMILES, so the term is computable from a prediction alone.

---

## The numbering control, which had to pass first

Predictions number 1..N from the construct; every residue named above is UniProt P08684
numbering. Getting this wrong is silent — FINDING 021 is a whole finding about it.

The offset is resolved by **residue-name agreement against the canonical CYP3A4 sequence**,
never against the crystal, so the method works on a blind target. It agrees with the
crystal-derived offset on **1,740 of 1,740 poses** (offsets +20/+21/+22/+23 across four
construct variants; minimum sequence identity at the chosen offset 0.984). All five polar
anchors and all eight cluster phenylalanines are present in every model.

A planar ring was found in 94.5% of poses; the misses are macrolides, which have none.

---

## The pool, the incumbent, and two nulls re-derived here

| | value |
|---|---|
| poses / ligands | 1,740 / 87 |
| **pool oracle** | **0.6931** |
| random selection | **0.5757** |
| **incumbent (XENG, frozen independent references, no fitted parameters)** | **0.6140 — +0.0383**, ρ +0.231, 71.3% of ligands positive |

**Null 1, random selection** (4,000 draws): p95 **+0.0140**, p99 +0.0196.

**Null 2, matched dimensionality** — 34 Gaussian features through the identical
leave-one-ligand-out HistGradientBoosting fit and identical randomised tie-breaking, **200
repetitions**: mean +0.0008, sd 0.0083, p95 **+0.0138**, p99 +0.0196, **max +0.0252**.

That p95 of +0.0138 is re-derived from scratch and lands on FINDING 007's +0.0138. The
noise floor was not borrowed and it has not moved.

---

## Every individual term is a null

Sign fixed **a priori by chemistry** (more H-bonds better, shorter better, more stacking
better, less strain and clash better) — no fitting, not even a sign. Ties broken at random
over 64 draws.

| term | gain | within-ρ | frac ρ>0 | n | p |
|---|---|---|---|---|---|
| `strain` (MMFF94s) | **+0.0105** | −0.037 | 0.43 | 86 | 0.76 |
| `hb_arg212` | +0.0043 | +0.036 | 0.53 | 34 | 0.66 |
| `hb_backbone` | +0.0006 | −0.009 | 0.49 | 57 | 0.91 |
| `phe_n_res` | +0.0003 | +0.073 | 0.54 | 83 | 0.64 |
| `n_hb_backbone` | −0.0005 | +0.000 | 0.51 | 57 | 0.98 |
| `phe_stack_score` | −0.0015 | +0.053 | 0.56 | 81 | 0.94 |
| `anchor_score` | −0.0018 | +0.016 | 0.46 | 85 | 0.53 |
| `phe_min_dist` | −0.0021 | +0.056 | 0.55 | 82 | 0.76 |
| `d_ser119` | −0.0048 | −0.011 | 0.52 | 87 | 0.51 |
| `hb_ser119` | −0.0095 | −0.017 | 0.43 | 82 | 0.67 |
| *CONTROL* generic contact count | −0.0125 | +0.119 | 0.64 | 87 | 0.48 |
| `phe_heavy_contacts` | −0.0159 | +0.113 | 0.59 | 87 | 0.47 |
| `fg_contacts` | −0.0214 | +0.024 | 0.51 | 87 | 0.06 |
| `fg_frac` | −0.0268 | −0.039 | 0.46 | 87 | 0.006 |
| `f304_dist` | −0.0279 | −0.027 | 0.42 | 87 | 0.01 |
| `fg_depth` | −0.0376 | −0.106 | 0.39 | 87 | 0.006 |

(24 further terms between −0.02 and +0.005; full table in `pocket_singles_apriori.json`.)

**Not one term clears +0.0138.** The best is strain, at +0.0105 with a *negative* within-
ligand ρ — which is the PXR campaign's MMFF strain verdict replicating exactly, on a
different protein, a different pool and a different implementation. Design §D predicted
this and it is now logged rather than quietly dropped.

The four terms at the bottom are significantly *negative* with the physics sign. That is
not a hidden positive waiting to be harvested: reversing `fg_depth` gives **−0.0049**, not
+0.0376, because argmin and argmax of the same feature can both sit below the pool mean.
Top-1 selection is not a monotone function of ρ — the third trap in `docs/README.md`.

### Two anchors from the design do not exist in this pocket

The FINDING 014 pre-check, applied first: `hb_asp214` is **constant within 98.9% of
ligands** (median ligand-to-Asp214 distance 11.4 Å), `hb_thr224` within 89.7%, `hb_arg106`
within 74.7%, `hb_arg212` within 60.9%. Only **Ser119** genuinely varies (constant within
5.7%, mean nearest-polar distance 3.10 Å) — and it is a null anyway. Design §E named four
residues; three of them cannot rank a CYP3A4 pose because they are nowhere near the ligand.

---

## The blocks are null too, and the ensemble is exactly at its noise ceiling

Leave-one-ligand-out gradient boosting, per block and over all 34:

| model | selected | gain | ρ | frac ρ>0 | beats random | p |
|---|---|---|---|---|---|---|
| GROUP polar (15) | 0.5769 | +0.0012 | +0.116 | 0.64 | 0.53 | 0.69 |
| GROUP phe (9) | 0.5808 | +0.0051 | +0.025 | 0.53 | 0.54 | 0.63 |
| GROUP fg (6) | 0.5757 | −0.0001 | +0.068 | 0.59 | 0.51 | 0.93 |
| GROUP strain (4) | 0.5854 | +0.0096 | −0.040 | 0.46 | 0.56 | 0.37 |
| **ENSEMBLE (34)** | **0.6010** | **+0.0252** | +0.093 | 0.64 | **0.61** | 0.0022 |
| ABLATE −polar | 0.5818 | +0.0060 | +0.057 | 0.55 | 0.49 | 0.50 |
| ABLATE −phe | 0.5962 | +0.0204 | +0.099 | 0.62 | 0.62 | 0.0065 |
| ABLATE −fg | 0.5786 | +0.0029 | +0.081 | 0.60 | 0.48 | 0.87 |
| ABLATE −strain | 0.5703 | −0.0054 | +0.069 | 0.54 | 0.48 | 0.97 |
| ENSEMBLE + incumbent, rank-avg | 0.6019 | +0.0262 | +0.196 | 0.68 | 0.58 | 0.013 |
| *CONTROL* Boltz sample rank | 0.5869 | +0.0112 | +0.063 | 0.60 | 0.60 | 0.037 |

Four blocks that are each individually null combine to +0.0252, and removing **any** of
them collapses it. That is the signature of a fit assembling a weak multivariate read, not
of chemistry — and the null says how weak:

> **p(matched-dim null ≥ +0.0252) = 0.0050, one draw in 200 — and that draw reached
> +0.0252, the observed value exactly.** Thirty-four Gaussian noise features, fitted the
> same way, matched the physics scorer's best result once in two hundred tries.

Significant at p = 0.005, and sitting precisely on the maximum its own null reached. This
is the number to quote, not the Wilcoxon 0.0022.

---

## Controls

| control | result | reads as |
|---|---|---|
| within-ligand feature shuffle, 5 reps | mean **+0.0012** (−0.016…+0.013) | the pose-to-feature pairing is what carries the +0.0252; marginals and dimensionality alone give nothing |
| numbering offset vs crystal | **1,740/1,740 agree**, min identity 0.984 | non-leaky numbering, transfers to a blind target |
| generic contact count | −0.0125 | the FINDING 003 baseline is also null on *this* pool, so nothing here is that result in disguise |
| Boltz sample rank | +0.0112 | tie-breaking matters: a constant score would have reported this |
| seed resampling | sd **0.0000** | **vacuous, not clean** — `HistGradientBoostingRegressor`'s `random_state` only drives the binning subsample, which does not trigger below 10,000 rows. The deploy rule's "survives resampling" clause is untested here, and reporting sd 0 as stability would have been the "too-clean numbers" trap |

---

## Why it does not ship, which is decided by the gate rather than by the gain

`CLAUDE.md`: *a term ships only if it beats the incumbent on held-out folds and its gain
correlates with where the incumbent errs.* Both halves fail.

| | |
|---|---|
| physics ensemble vs incumbent | **0.6010 vs 0.6140** — loses by 0.0130, while fitting 34 features against the incumbent's zero |
| per-ligand, physics beats incumbent | 48.3% (Wilcoxon p = 0.71) |
| correlation of physics per-ligand gain with incumbent shortfall | **r = −0.047, p = 0.66** |
| correlation of the two selectors' per-ligand outcomes | **+0.834** |
| combined, rank-average | +0.0262 — below the incumbent alone |
| tail rescue on the 10 / 20 / 30 widest-spread ligands | 0.6030 / 0.5991 / 0.6022 — all below 0.6140 |

The physics scorer is not seeing a different failure mode. It is a weaker read of the same
one, and every way of mixing it in makes the board worse.

---

## One post-hoc direction, tested where it was not found

`fg_depth` was the most negative single term, i.e. the data said "lower ligand is better"
against the physics sign. Reversing a sign after seeing the answer is a hypothesis, so it
was taken to the **arm4 P450 holdout — 85 pairs of other P450 proteins, other ligands,
5 samples each**, scored during the fine-tuning campaign and never touched here.

| feature, "lower is better" | CYP3A4 pool (87 lig, 20 poses) | P450 holdout (85 pairs, 5 poses) |
|---|---|---|
| `fg_depth` (centroid height over porphyrin) | −0.0049 | −0.0051 |
| `fe_centroid_dist` | +0.0013 | **+0.0408**, ρ +0.138, p_perm **0.000** |
| `fe_min_dist` | +0.0040 | +0.0291, ρ +0.133, p_perm 0.006 |

`fg_depth` is a null in both directions and in both pools — there was nothing to harvest.
"Closer to the iron" *does* select strongly on the holdout, where 39 of 85 pairs are
catastrophes, and is worth nothing on the CYP3A4 pool, where they are rare. That is
FINDING 012's catastrophe-detector scaling reproducing on a new feature, and it is a
statement about pool quality, not about CYP3A4 chemistry.

---

## Verdict

**The CYP3A4-specific physics scorer is a null.** Fifteen residue-specific polar-anchor
terms, nine phenylalanine-cluster stacking terms, six F/G-roof engagement terms and four
strain/clash terms, measured over 1,740 poses with crystal truth: not one of the 34 clears
a +0.0138 noise floor, not one of the four blocks clears it either, and the fitted ensemble
reaches the exact maximum that 200 draws of matched-dimensionality noise reached.

Three specific things this kills, beyond the aggregate:

1. **Design §E is wrong about three of its four residues.** Arg106, Arg212 and Asp214 are
   constant within 61–99% of ligands — they are not in contact. Ser119 is real and is a
   null.
2. **The phenylalanine cluster, the most CYP3A4-specific recognition feature available and
   the one the design omits, does not select** — +0.0051 as a block, and removing it from
   the ensemble costs only 0.0048.
3. **MMFF strain replicates the PXR negative**, on a different protein with a different
   implementation. Design §D's own warning was correct and the term can now be closed.

This is the fifth negative on the same axis, and it says the same thing the other four do:
the discriminating error is ligand-specific substituent placement, and a fixed vocabulary
of pocket features — however CYP3A4-specific, however correct as chemistry — averages over
exactly the thing that varies. The incumbent works because it compares a pose to *other
predictions of the same molecule*, which is the only representation so far that is
conditioned on the ligand.

**Nothing from this finding enters the submission pipeline.** `pocket.py` stays in the
tree as a measured negative and as a place to hang a term that is conditioned on the
ligand, should one be found.

---

## Reproduce

```
scripts/qmscore/extract_pocket_terms.py   # 1,740 poses -> pocket_terms_run5_s4*.csv
scripts/qmscore/singles_apriori.py        # per-term, sign fixed by chemistry
scripts/qmscore/eval_pocket_scorer.py     # blocks, ensemble, ablations, LOO
scripts/qmscore/null_matched_dim.py       # 200-rep matched-dimensionality null (array)
scripts/qmscore/controls_pocket.py        # shuffle, complementarity, tail rescue
scripts/qmscore/depth_replication.py      # the post-hoc direction, on the P450 holdout
```

Artifacts on Explorer under `/scratch/shenoy.am/zexp/`: `pocket_terms_run5_s4*.csv`,
`pocket_eval.json`, `pocket_singles_apriori.json`, `pocket_controls.json`,
`pocket_within_cv.csv`, `nullgains_all.npy`, `depth_replication_*.json`.

### A method note worth keeping

The first pass oriented each single feature by its **out-of-fold mean within-ligand ρ**.
That rule degenerates on precisely the features it is meant to test: when a term's true
mean ρ is zero, the leave-one-out mean is dominated by the removal of the held-out
ligand's own ρ, so the chosen sign is anti-correlated with it. `n_hb_backbone` came back at
ρ **−0.316 with 0 of 57 ligands positive and p = 9 × 10⁻⁹**, a result that would have been
extremely publishable and is entirely an artifact — its raw mean ρ is **+0.0005 on 51% of
ligands**. Fixing the sign a priori removes it. **Leave-one-out sign selection manufactures
confident negatives from nulls**, and nothing else in this repo has used it before.
