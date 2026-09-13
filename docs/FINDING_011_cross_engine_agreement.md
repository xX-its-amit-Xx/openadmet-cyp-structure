# FINDING 011 — cross-engine agreement is the strongest single ranking signal yet, and it still does not add to the selector

**Date:** 2026-09-13
**Pool:** 1,500 unsteered Boltz poses, 75 CYP3A4 ligands, against Protenix-v2 predictions
of the same ligands
**Script:** `scripts/structure/cross_engine_agreement.py`

---

## A different question from FINDING 005

FINDING 005 asked whether a second engine can **replace** the first on the ligands it gets
wrong. No: Chai's errors correlate with Boltz's at rho +0.45, Protenix at +0.474, and tail
rescue does not transfer.

This asks whether a Boltz pose that an independent engine **reproduces** is more likely to
be right. Consensus is the one thing that has ever worked here, and FINDING 010 explains
why - it is a within-ligand, ligand-specific comparison. But `mean_rmsd_to_others` only
measures whether Boltz agrees with *itself*. A different architecture agreeing is strictly
more information than the same model agreeing twice.

Agreement is a symmetric **Chamfer distance** between ligand atom clouds, each pose placed
in its own heme frame. Permutation-invariant, so the atom-ordering trap that once halved
every LDDT-PLI cannot apply, and superposition-free, so there is no reference frame to get
wrong.

---

## The ranking signal is real and strong

| statistic | cross-engine `xeng_min` | Boltz's own confidence (FINDING 001-E) |
|---|---|---|
| within-ligand Spearman vs LDDT-PLI | **-0.202** | -0.033 |
| ligands in the correct direction | **72%** | ~50% |
| Wilcoxon p | **0.00053** | 0.25 |

Six times the magnitude of the engine's own confidence, in the right direction on nearly
three ligands in four. This is the strongest within-ligand ranking statistic measured in
this project.

## The selection gain does not clear the bar, and the combination is worse

Random baseline 0.5769, oracle 0.6967, null 95th pct +0.0124 / 99th pct +0.0178.

| selector | selected | gain | p |
|---|---|---|---|
| pocket contacts alone | 0.5805 | +0.0036 | 0.28 |
| consensus alone | 0.5899 | +0.0131 | 0.045 |
| **cross-engine alone** | 0.5923 | **+0.0154** | 0.024 |
| **incumbent `0.5*zc - zm`** | 0.6009 | **+0.0240** | 0.0003 |
| incumbent - 0.5*zx | 0.5968 | +0.0199 | 0.0037 |
| incumbent - 1.0*zx | 0.5954 | +0.0186 | 0.0073 |
| incumbent - 1.5*zx | 0.5941 | +0.0173 | 0.013 |

Alone it sits between the 95th and 99th percentile of the null. Added to the incumbent it
**costs -0.004 at every weight tried**, monotonically worse as the weight grows.

And it is not redundancy in the usual sense. Within-ligand correlations are
`zx` vs `zc` **-0.036**, `zx` vs `zm` **+0.037** - as orthogonal as the incumbent's own two
terms are to each other (-0.035), which is the property that made *those* combine. So this
is a third pattern, distinct from the two FINDING 006 recorded: not "redundant" (it is
not correlated with what is there) and not "orthogonal but empty" (it carries real
signal). It is orthogonal, individually positive, and still subtractive.

---

## The caveat that decides what to do next

The Protenix side of this comparison is **one ligand pose per ligand**. Because of FINDING
009, replicate 0's twenty models share a single ligand conformation, so `xeng_min` is
currently the distance to *one* independent guess rather than to a distribution. That is
the noisiest possible version of this feature, and it still produces rho = -0.20 at
p = 0.0005.

So this is a **wait for data**, not a dead end. The replicate campaign is building 6
independent Protenix poses per ligand; re-run this then. If the signal scales with the
number of independent poses to agree with, the combination may stop being subtractive. If
it does not, the feature is retired.

Logged as a negative for now, per the repo rule: anything absorbed or degraded by the
incumbent is a negative, and "promising" is not a result.
