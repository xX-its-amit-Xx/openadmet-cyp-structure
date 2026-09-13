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

---

# UPDATE — the pre-registered re-test passed, and cross-engine agreement now beats the incumbent

**Date:** 2026-09-13, same day. **Pool:** 1,260 Boltz poses, **63 ligands with >= 4
independent Protenix poses each** (median 8).

The prediction above was: *if the signal scales with the number of independent poses to
agree with, the combination may stop being subtractive; if not, retire the feature.*
It scaled.

## The dose-response is clean

| feature | 1 pose/ligand | **>= 4 poses/ligand** |
|---|---|---|
| `xeng_min` | +0.0154, p = 0.024 | **+0.0225, p = 0.0035** |
| `xeng_mean` | **-0.0055**, p = 0.65 | **+0.0284, p = 0.0000** |

`xeng_mean` inverted from useless to best, which is the mechanism working exactly as
stated: averaged over one pose it *was* that pose, and averaged over eight independent
ones it is an estimator of where the other engine actually thinks the ligand goes.

**A bug caught before it could fake this result.** The collection code took the first N
*files* per ligand, and replicate 0 alone contributes 20 models sharing a single ligand
conformation (FINDING 009). That would have handed back eight copies of one pose and
called them eight independent opinions - the FINDING 009 trap, reintroduced by the code
written to exploit FINDING 009. Poses are now keyed on replicate index.

## It beats the incumbent, with no fitted parameters

Random 0.5736, oracle 0.6963, null 95th +0.0134 / 99th **+0.0204** (3,000 draws).

| selector | selected | gain | p | fitted params |
|---|---|---|---|---|
| incumbent `0.5*zc - zm` | 0.5943 | +0.0207 | 0.0097 | none |
| **`-zx` alone** | 0.6019 | **+0.0284** | **0.0000** | **none** |
| `-zm - zx`, equal weights | 0.6041 | **+0.0305** | 0.0000 | none |
| `0.5*zc - zm - 2*zx` | 0.6012 | +0.0276 | 0.0007 | none |
| **LOCO, weights fit out-of-fold** | 0.5885 | **+0.0149** | 0.036 | 3 |

Two things to read carefully here.

**The headline is honest.** `-zx` has **no fitted parameters**, so there is no train/test
split to get wrong and no held-out fold to demand - a fixed feature cannot overfit a
weight it does not have. It was pre-registered before the data existed, and 0 of 3,000
random draws beat it. Both variants the script computes cleared the 99th percentile, so
this does not rest on picking the better of two after the fact.

**The fitted combination does NOT validate.** Fitting three weights leave-one-scaffold-
cluster-out drops it to +0.0149, inside the noise floor. That replicates FINDING 002's
result that a fitted ranker (+0.0048) loses to simple fixed selectors, and it is the
reason the shipping recommendation is a fixed unweighted term, not a tuned blend.
`-zm - zx` at +0.0305 is the best number on the board but was chosen post-hoc from about
eight tried, so treat it as provisional against the fixed `-zx`.

## Why this one works when thirty others did not

Both surviving terms are **within-ligand consensus**: does this pose agree with other
opinions about *this molecule*. `zc` (pocket contacts), the one geometric term in the
incumbent, is the part that drops out - `-zm - zx` beats `0.5*zc - zm - zx`.

That is FINDING 010's thesis holding up under a real test rather than a negative one.
Population-level priors cannot select because predicted poses already agree on everything
population-level. Agreement between independent opinions about one specific molecule can.
And an architecturally different engine is a better second opinion than the same engine
sampled twice, **even though its errors correlate with Boltz's at rho = +0.596** - which
is the genuinely surprising part, and the thing FINDING 005's framing would not have
predicted.

## Caveats

- **n = 63**, not 87; the incumbent scores +0.0207 on this subset against +0.0279
  reported on the full set, so the subset is not a neutral slice.
- Requires a second engine's pool at inference. That is affordable here (OpenProtein is
  unmetered) but it is a real dependency, not a free feature.
- Re-run at 12 replicates and on all 87 ligands before this goes into a submission.
