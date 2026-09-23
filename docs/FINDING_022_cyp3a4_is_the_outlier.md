# FINDING 022 — Co-folding has essentially solved the P450 superfamily. CYP3A4 is the exception.

**Date:** 2026-09-20 · **Status:** measured, 85 held-out pairs, numbering-corrected · **Confidence:** high, cross-validated against an independent pool

## The measurement

Same checkpoint, same MSAs, same 5 samples, same scorer, numbering aligned (FINDING 021):

| | n | mean LDDT-PLI | median | sub-2 Å | median BiSyRMSD |
|---|---|---|---|---|---|
| **CYP3A4 (P08684)** | 15 | **0.5550** | 0.5313 | **5/15 (33%)** | **2.842 Å** |
| every other P450 | 70 | **0.8648** | 0.9480 | **62/70 (89%)** | **0.615 Å** |

Per target, for those with n ≥ 4:

| uniprot | n | mean | sub-2 Å |
|---|---|---|---|
| **P08684 (CYP3A4)** | 15 | **0.555** | 0.333 |
| P11511 (aromatase) | 6 | 0.863 | 0.833 |
| P10614 (CYP51A1) | 5 | 0.909 | 1.000 |
| Q2IU02 | 27 | 0.931 | 0.963 |

A gap of **0.31 LDDT-PLI** and **4.6× in ligand RMSD**, inside one enzyme family, from one
model on one day.

## It is not a pipeline artifact

This repo's CYP3A4 pool — built on Modal, different MSAs, different sampling, scored
months earlier by the same `cypstruct.pose` — reports a per-pose mean LDDT-PLI of
**0.5680** across 3,360 poses. The held-out CYP3A4 pairs here score **0.5550**.

Two independent pipelines landing within 0.013 of each other on CYP3A4, while the same
code puts every other P450 at 0.86, is about as clean a cross-validation as this project
is going to get. The CYP3A4 pool also carries essentially no exact-zero scores (4 of
3,360), so it was never touched by the numbering bug that contaminated the held-out set.

## What it changes

**1. "Co-folding is bad at ligand orientation" is too broad.** It is near-solved for most
P450s: 89% of held-out poses under 2 Å, median 0.615 Å. The failure is *specific to
CYP3A4*, which is exactly what OpenADMET said when they attributed it to F/G loop
remodelling. Their claim now has an independent number attached: 0.555 against a
family baseline of 0.865.

**2. It retro-explains FINDING 018.** arm4_mix trained on a set where ~82% of targets are
ones the model already solves at 0.86+. Gradient from an already-correct pose is at best
noise and at worst a pull away from whatever CYP3A4-specific behaviour the pretrained
weights encode. The fine-tune being monotonically negative is consistent with having
spent 350 steps mostly on solved problems. *This is a plausible mechanism, not a
demonstrated one* — the dose-response says the lever is wrong, and this says one reason
it might be.

**3. It makes superfamily benchmarking nearly worthless.** The validation ladder this
project planned — optimise on an analogous target, hop to a second, then to CYP3A4 —
needs analogues that are HARD. Picking P450s by family similarity selects targets with
**0.14 of headroom and 89% already solved**. Any method tuned there is tuned on noise,
and a gain of +0.02 on a 0.86 baseline says nothing about a 0.55 one.

**4. It sharpens where the remaining work is.** On CYP3A4 the model is at 0.555 with
10/15 poses over 2 Å. That is the regime this repo's selection results were measured in,
and where the 15-pose oracle headroom of +0.0632 mostly lives. The target is not "make
co-folding better"; it is "make co-folding better *on this one enzyme*".

## What to do with it

- **Reselect the analogue set by difficulty, not by family.** The criterion should be a
  measured sub-2 Å rate near CYP3A4's 33%, whatever family that comes from. There are
  185 targets in the P450 universe harvest and most are the wrong kind of easy.
- **Re-examine arm3_cyp3a4_only.** It was set aside for having held-out novelty 0.758
  (too close to its training set) — but it is the only arm whose training signal is
  concentrated where the model actually fails. That is a different objection than the one
  that killed arms 1 and 2, and worth one run *if* a difficulty-matched holdout can be
  built for it.
- **Report CYP3A4 numbers separately from family numbers, always.** A mean over 85 P450
  pairs is 0.81 and says nothing about the challenge; the 15 CYP3A4 pairs say 0.555.

---

# Correction, 2026-09-23 — the per-target table reports **sample 0**, not the sample mean

Appended by `FINDING_039`. No number in the body is wrong; the column label is ambiguous and
it is what made this table look as though it disagreed with `FINDING_024`.

Every figure in the "mean LDDT-PLI" column — 0.5550 for CYP3A4, 0.8648 for every other P450,
0.863 / 0.909 / 0.931 per target — is the mean of **`lddt_sample0`**, the engine's own rank-0
pose. It is **not** the mean over the five samples. Both reproduce exactly:

| target | n | **mean of sample 0 (this table)** | mean of all 5 (what `FINDING_024` prints) |
|---|---|---|---|
| P08684 CYP3A4 | 15 | **0.5550** | 0.5526 |
| **P11511 aromatase** | 6 | **0.8629** | **0.6313** |
| P10614 CYP51A1 | 5 | **0.9093** | 0.9076 |
| Q2IU02 | 27 | **0.9306** | 0.9241 |
| every other P450 | 70 | **0.8648** | 0.8466 |

**Aromatase is the one target where the two differ materially** — by **0.232**, against
≤ 0.007 everywhere else. Boltz-2 puts a very good aromatase pose at rank 0 and much worse
ones behind it, which is why 024's 0.631 and this table's 0.863 are both right.

**Prefer the all-5 mean for cross-target comparison.** Rank 0 is the engine's own confidence
ordering, and `FINDING_001` is that this ordering does not rank poses. The headline gap
survives either way: CYP3A4 0.555 / 0.553 against every other P450 0.865 / 0.847.
