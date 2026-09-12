# Finding 007 — the noise floor, and what it says about everything tested so far

**Date:** 2026-09-12 · **Method:** 2,000 random features drawn per-pose, each used to select
one pose per ligand, on the real `val87b` pool (1,740 poses, 87 ligands)

I had tested roughly twenty features and variants, describing several as "promising but not
significant" at +0.015 to +0.017. Before trusting that vocabulary — or the selector — it is
worth knowing what **pure noise** scores on this exact test.

---

## The null

A random value per pose, argmax per ligand, delta against the per-ligand mean. Pool sizes
and per-ligand structure preserved.

| statistic | delta |
|---|---|
| mean | +0.0001 |
| sd | 0.0087 |
| 90th percentile | +0.0113 |
| 95th percentile | +0.0138 |
| **99th percentile** | **+0.0195** |

## 1. The selector is real

**+0.0279, empirical p = 0.0010** against this null — above the 99th percentile of random
features. That is an independent confirmation of FINDING 003 by a completely different
route from the Wilcoxon test, and it survives the leave-one-scaffold-cluster-out value too:
**+0.0220 is still above the 99th percentile.**

## 2. Every "promising" candidate was noise, and now demonstrably so

| candidate | delta | beaten by this share of *single* random features |
|---|---|---|
| `frac_access` (sub-pocket) | +0.0174 | 1.7% |
| `missing_common` | +0.0168 | 2.1% |
| `is_coordinated` | +0.0099 | 13.8% |
| `xeng_support` (cross-engine) | +0.0066 | 23.8% |
| `contact_consensus` | +0.0028 | 38.3% |
| `prior_recall` | −0.0160 | 96.5% |

At face value 1.7% and 2.1% look like near-misses. **They are not, because I did not draw
one feature — I drew about twenty.** The expected maximum of 20 random draws lands near the
95th–97.5th percentile, i.e. **+0.0138 to +0.0160**. So +0.0168 and +0.0174 are almost
exactly what the best of twenty random tries produces. They were never signal.

This explains cleanly why every one of them degraded the selector when combined: there was
nothing to combine.

## 3. The decision rule I should have had from the start

> **A single-feature delta below about +0.020 on this dataset is indistinguishable from the
> best of ~20 random tries. Do not call it promising, do not combine it, do not spend a
> tick on it.**

Two ways past that bar, and they are the only two:

- **Clear it outright** — above +0.020, with the multiple-comparison count stated.
- **Compose with something orthogonal** — which is how the selector itself works: two terms
  at +0.0108 and +0.0100, both individually inside the noise, correlated at only ρ = +0.030,
  giving +0.0279 together. Orthogonality is what converts two noise-level signals into one
  real one, and it is measurable in advance.

## 4. Consequence for the sample size

The null's sd of 0.0087 over 87 ligands is the binding constraint on discovery here. To
resolve a +0.010 feature reliably would need roughly 4× the ligands, which we do not have
and cannot buy — the validation set is bounded by how many CYP3A4 holo structures exist.

So the realistic path is **not** a search for more single features. It is either a term
with a large effect, or orthogonal composition, and both should be screened with the
within-ligand CV check from FINDING 006 first and this noise floor second.
