# FINDING 023 — Boltz-2's trunk pair representation cannot rank poses, because it does not depend on the pose

**Date:** 2026-09-20 · **Status:** measured, 87 ligands x 20 poses, all controls clean · **Verdict:** the premise was wrong

## The premise, and why it was wrong

`docs/worldmodel/README.md` Revision 2 argued that the literature's case against learned
embeddings — that they lose to count-ECFP4 — does not apply to Boltz-2's pair
representation `z`, because **"`z` is a function of protein x ligand x pose"** and is
therefore joint by construction rather than assembled after the fact.

That sentence is false. Boltz-2's trunk runs **once per input, before diffusion
sampling.** Its pair representation is a function of protein x ligand and nothing else,
so it is *identical across every pose of the same ligand*. A feature that is constant
within a ligand cannot choose between that ligand's poses, however rich it is.

The measurement says exactly that:

| feature | selected | gain | within-rho | frac ligands positive |
|---|---|---|---|---|
| `z_trunk_static` (ridge) | 0.5710 | **-0.0047** | **0.011** | 0.506 |
| `z_trunk_static` (medoid) | 0.5705 | -0.0052 | -0.042 | 0.391 |
| ECFP4 ligand-only *control* | 0.5737 | -0.0020 | **0.000** | 0.000 |

The trunk representation lands on the ligand-only control, because it **is** a
ligand-only feature.

## What does carry signal, and what it actually is

| selector | selected | gain | rho | frac pos | p |
|---|---|---|---|---|---|
| `gbm_LOO_z_conf_all` | 0.6074 | **+0.0316** | 0.209 | 0.759 | <1e-4 |
| `gbm_LOO_z_trunk_contactgated` | 0.6071 | +0.0314 | 0.241 | 0.713 | <1e-4 |
| `ridge_LOO_z_conf_all_wcenter` | 0.6072 | +0.0315 | 0.191 | 0.667 | <1e-4 |
| `boltz_confidence_score` alone | 0.5899 | +0.0141 | 0.088 | 0.644 | 0.046 |
| **incumbent cross-engine (reference)** | **0.6164** | **+0.0395** | -0.258 | 0.759 | — |
| best **unsupervised** z medoid | 0.5892 | +0.0134 | 0.153 | 0.713 | 0.055 |

Two things follow.

**The working features are pose-dependent reads, not the trunk.** `z_conf_*` comes from
the confidence module, which *does* see sampled coordinates. `z_trunk_contactgated` is
the trunk representation gated by the pose's actual contacts — the gating is what
reintroduces pose dependence. Remove the pose dependence (`z_trunk_static`) and the
signal vanishes entirely.

**So this is a richer read of confidence, not a new signal.** The scalar
`boltz_confidence_score` gives +0.0141 on this pool; a 128-dimensional learned read of
the same module gives +0.0316. A real improvement over the scalar — and still **below
the parameter-free incumbent at +0.0395.**

## The comparison is worse than the numbers suggest

Every `z` selector above is **fitted** — ridge or gradient boosting, leave-one-ligand-out
over these 87 ligands. The incumbent has **zero fitted parameters**. A fitted selector
that loses to an unfitted one has lost twice.

The like-for-like comparison is the unsupervised medoid, and it does not clear
significance: best +0.0134 at **p = 0.055**.

## Controls, all clean

| control | gain | rho | p | reads as |
|---|---|---|---|---|
| ECFP4 ligand-only | -0.0020 | 0.000 | 0.60 | no within-ligand information, by construction |
| `z_conf` shuffled within ligand | +0.0107 | -0.011 | 0.10 | destroying the pose-to-feature pairing kills it |
| random noise, matched dim | +0.0091 / -0.0096 | ~0 | 0.14 / 0.87 | floor is where it should be |
| `boltz_iptm`, `boltz_ligand_iptm` | -0.0000 | -0.03 | 0.52 | exactly chance |
| `boltz_complex_ipde` | -0.0100 | -0.047 | 0.88 | worse than random, as previously measured |

The ECFP4 control is what validates the split: its within-ligand rho is *exactly* 0.000
because it is constant within a ligand. On a random-row split it would have looked
predictive by memorising ligand identity; leave-one-ligand-out reduces it to nothing,
which is the correct answer.

A methodological note worth keeping: **ties were broken at random over 64 draws.** A
constant score would otherwise always select `_model_0`, the highest-confidence sample,
and report Boltz confidence's number under another name. `z_trunk_static` is precisely
such a constant score, so without this the null would have looked like a small positive.

## What survives

1. **The shared-embedding premise, as stated, is dead for pose selection.** The
   representation argued to be joint-by-construction is not, and the thing that does work
   is a confidence read that loses to a parameter-free geometric consensus.
2. **One open question, cheap to answer:** the fitted z-confidence selector (+0.0316,
   rho +0.209) and the incumbent (+0.0395, rho -0.258) have *opposite-signed* within-ligand
   correlations while agreeing on the same 75.9% of ligands. Whether they are
   complementary or redundant is a combination test on data already on disk.
3. **For the world model more broadly:** a genuinely pose-aware joint representation from
   a co-folder has to come from the diffusion or confidence path. The trunk is the wrong
   tap, and "it is a big tensor from a structure model" is not an argument that it knows
   about the structure.


---

## Addendum — the comparison in this finding was not like-for-like

The verdict above said the fitted z-confidence selector (+0.0316) "still loses to the
incumbent at +0.0395". That +0.0395 was measured on a **different pool**, and the
comparison does not hold once both are computed on the same one.

Geometric cross-engine consensus, run with `cypstruct.xengine`'s `in_heme_frame` and
`chamfer` unchanged, on **this** 20-pose pool (1,740 poses, 87 ligands, all merged):

| | value |
|---|---|
| pool oracle | 0.6931 |
| random | 0.5757 |
| **geometric consensus** | **0.5927 (+0.0170)** |
| within-ligand rho | −0.201, correct direction on **72.4%** of ligands |
| permutation null (2,000 draws) | p95 +0.0139, p99 +0.0199 |
| p(null ≥ observed) | **0.0195** |

So on the same pool, **+0.0170 for consensus against +0.0316 for the z-confidence read.**
The z feature wins, and the consensus barely clears its own 95th-percentile null.

### Why consensus collapses here, and why that is not a refutation of FINDING 011

This pool is **four seeds of one engine**, not four engines. FINDING 011 is explicit that
the selector needs *"≥ 4 GENUINELY independent reference poses"*, and same-engine
replicates are not independent — FINDING 015 went further and found both OpenProtein
engines deterministic, so replicate depth was already known to be thin.

That is exactly what +0.0170 with 72.4% directional accuracy looks like: the mechanism
still points the right way, but with correlated references it has little to work with.

### The corrected statement

- Against a **same-engine** pool, the fitted z-confidence read is the better selector
  (+0.0316 vs +0.0170).
- The **deployed** incumbent at +0.0395 uses genuinely independent engine references, and
  remains ahead of both — but it buys that with extra generation compute the z feature
  does not need, since `z` falls out of a run we already pay for.
- Everything else in this finding stands: the **trunk** representation is pose-independent
  and is a null (−0.0047, rho 0.011), and the working feature is a confidence read.

### A sign error, caught by a control rather than by inspection

The first computation of the consensus number used `max()` over the consensus score and
returned **−0.0598** — a large negative. The score is a *distance*, so lower is better,
and `max()` selects the worst pose in every pool. `CLAUDE.md` lesson 2 warns about this
exact failure ("PDE and PAE are *errors*, so a raw `max()` picks the worst pose") and it
still happened.

What caught it was not reading the code again. It was that the within-ligand rho came back
**−0.201 with only 27.6% of ligands positive** — a feature that anti-correlates that
consistently is a correctly-signed feature being read backwards. Reporting the selected
mean alone would have produced a confident, wrong, and very publishable-looking result.


---

## Addendum 2 — the first addendum overcorrected. The original verdict was right.

Addendum 1 claimed the like-for-like comparison was +0.0170 for consensus against +0.0316
for the z read, and concluded "the z feature wins". **That substituted an invalid
reference set for a valid one and called the result like-for-like.**

The +0.0170 was cross-engine consensus computed with the *four Boltz seeds as each
other's references* — which FINDING 011 explicitly says is not the selector. The deployed
selector scores against a **frozen reference set built from independent engines**. Running
that, on exactly these poses:

| | frozen independent references |
|---|---|
| reference depth | 87 of 87 ligands, min 6 / median 7 / max 12, **all ≥ 4** |
| poses scored | 1,740 over 87 ligands |
| pool oracle | 0.6931 |
| random | 0.5757 |
| **XENG (deployed, no fitted parameters)** | **0.6140 — gain +0.0383** |
| within-ligand rho | −0.231, correct direction on 71.3% |
| null (2,000 draws) | p95 +0.0139, p99 +0.0199, **p = 0.0000** |

### The complete ranking, one pool, one truth table, identical tie-breaking

| selector | gain | fitted? |
|---|---|---|
| **XENG, frozen independent references** | **+0.0383** | no |
| `gbm_LOO_z_conf_all` | +0.0316 | yes (LOO ridge/GBM) |
| geometric consensus, same-engine seeds | +0.0170 | no |
| scalar `boltz_confidence_score` | +0.0141 | no |
| null floor | +0.0139 (p95) | — |

**The original verdict stands.** The incumbent beats the fitted z-confidence read on the
same poses, while fitting nothing — and it replicates almost exactly on a pool it was
never tuned on (+0.0383 here against +0.0395 before, which is the strongest evidence yet
that FINDING 011 is real rather than a fit to its own validation set).

### What I got wrong, and the general form of it

I corrected a *valid* comparison by substituting a *degraded* version of the incumbent —
same-engine references instead of independent ones — and labelled the substitute
"like-for-like" because it ran on the same poses. Running on the same poses is not the
same as being the same method. The reference set is part of the selector, not part of the
data.

The tell was available and I did not act on it: the same-engine consensus scored +0.0170
against a p95 null of +0.0139, i.e. *barely above its own noise floor*, while the method
it was standing in for had been measured at +0.0395 across 81 proteins. A stand-in that
lands on the noise floor is not a stand-in.

### Still open

The combination question is unchanged and now sharper: XENG (+0.0383, rho −0.231) and
the z-confidence read (+0.0316, rho +0.209) are both real, and the second is *free* —
`z` falls out of a run already being paid for, while XENG needs poses from other engines.
If they are complementary, the combination is worth more than either. Testing it needs
per-pose z scores, which `analyse.py` currently computes internally and does not emit.


---

## Addendum 3 — the combination is a null. The line closes with the incumbent unchanged.

XENG and the z-confidence read have opposite-signed within-ligand correlations and only
+0.182 rho between them, so combining them looked promising. On the headline numbers it
was: rank-average +0.0431 against XENG's +0.0383, with a broad weight plateau rather than
a tuned spike.

**Paired over the 87 ligands, it is not a result.**

| | |
|---|---|
| XENG alone | 0.6140 |
| rank-average combination | 0.6187 |
| paired difference | **+0.0047** |
| bootstrap 95% CI (20,000 draws) | **[−0.0079, +0.0167]** |
| P(difference ≤ 0) | 0.227 |
| Wilcoxon p | **0.38** |
| ligands improved / worse / tied | **33 / 28 / 26** |

33 against 28 is a coin flip. The confidence interval crosses zero comfortably.

### Why this needed its own test

The +0.0139 figure this project uses as a noise floor is the 95th percentile of a gain
**versus random selection**. It is the right reference for "does this feature beat
picking a pose at random" and the *wrong* reference for "does selector A beat selector
B". A +0.0048 difference sails past the first bar and dies at the second, and reporting
it against the first would have been a real-looking result built on the wrong null.

### What the whole z line now says

| step | outcome |
|---|---|
| trunk pair representation `z` | **null** — pose-to-pose spread 1e-8, it is the same tensor for every pose |
| z-confidence read, fitted LOO | +0.0271 to +0.0316 — real, but a richer read of confidence |
| combination with XENG | **null** — +0.0047, p = 0.38 |
| **incumbent** | **unchanged at +0.0383 on this pool**, nothing fitted |

Three interventions, one mechanism each, and the shipped selector survives all of them.
That is the fourth consecutive line (fine-tuning, heme bond, ATOMICA, now embeddings)
where the honest answer was "the incumbent still wins", and the incumbent has now
replicated at +0.0395, +0.0357 and +0.0383 on three different pools without being
retuned once.
