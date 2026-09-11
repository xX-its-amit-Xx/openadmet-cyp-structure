# Finding 002 — nothing beats random yet, and every pool is multimodal

**Date:** 2026-09-11 · **Run:** `val87b`, unsteered arm · **n = 87 ligands, 1,740 poses**
**Data:** `physics_selector_val87b_unsteered.json`, `consensus_selector_val87b_unsteered.json`

FINDING 001 established that Boltz's own confidence selects *worse* than random, leaving
**0.127 LDDT-PLI per ligand** unclaimed between the pool oracle (0.6975) and the selected
pose (0.5706). This is the first attempt to claim it. **It did not work.**

---

## Every selector tested, against a random pose from the same pool

Ranked by gain. Oracle is 0.6975 throughout; random is 0.5769.

| selector | selected | Δ vs random | beats random on | p (Wilcoxon) | oracle captured |
|---|---|---|---|---|---|
| consensus: medoid | 0.5868 | **+0.0100** | 60.9% | **0.087** | 8.3% |
| single: `is_coordinated` | 0.5868 | +0.0099 | 50.6% | 0.57 | 8.2% |
| single: `frac_proximal` (low) | 0.5854 | +0.0085 | 49.4% | 0.74 | 7.0% |
| fitted ranker, leave-one-cluster-out | 0.5817 | +0.0048 | 55.2% | 0.52 | 3.9% |
| single: `fe_donor_dist` (low) | 0.5780 | +0.0011 | 52.9% | 0.69 | 0.9% |
| consensus: largest-cluster medoid | 0.5793 | +0.0024 | 50.6% | 0.92 | 2.0% |
| **Boltz confidence (incumbent)** | 0.5706 | **−0.0063** | 46.0% | 0.25 | −5.2% |
| single: `n_contacts` (high) | 0.5647 | −0.0122 | 49.4% | 0.50 | −10.1% |
| single: `s_fe_donor_angle` (high) | 0.5411 | **−0.0358** | 34.5% | **0.0007** | −29.7% |
| single: `max_clash` (low) | 0.5410 | **−0.0359** | 36.8% | **0.0024** | −29.8% |

**Nothing clears significance.** The best is the consensus medoid at +0.0100, p=0.087,
winning on 61% of ligands. That is suggestive and no more. The fitted ranker, trained on
all tier-0 terms and validated leave-one-scaffold-cluster-out over 33 clusters, does
worse than the single best raw term — the honest reading is that there is not enough
signal in these features for a model to find.

## The two significant results are both *negative*, and they are the informative ones

Selecting the pose with the most textbook trans Fe–donor angle, or the least steric
clash, is **significantly worse than random** (p=0.0007 and p=0.0024, both around
−0.036 LDDT-PLI).

That is not noise and it is not a sign error. It says the coordination anchor is
**necessary but not sufficient**. Once the donor nitrogen is pinned at 2.2 Å — which
FINDING 001 showed happens by itself, in 84% of poses — the molecule can still rotate
about the Fe–donor axis and fold its substituents anywhere. LDDT-PLI scores the contacts
of the *whole* ligand. So optimising the anchor geometry selects for an idealised local
arrangement while the rest of the molecule is free to be wrong, and the most
geometrically perfect anchor is if anything a mild marker of an over-collapsed pose.

**The discriminating information is in where the rest of the molecule sits, not at the
iron.** Every term tested so far is an anchor-local term. That is why they fail.

## Every single pool is multimodal

At a 2 Å clustering cutoff, **87 of 87 ligands** split into more than one cluster. Not
a handful of hard cases — all of them.

This matches what OpenADMET reported from the cryoEM density: for several ligands the map
is consistent with multiple mutually exclusive conformations rather than one pose. Our
pools are reproducing that ambiguity rather than resolving it.

## One thing that does work: pool agreement predicts difficulty

Largest-cluster fraction versus achievable quality: **ρ = +0.287, p = 0.0070.**
Mean pairwise spread versus quality: ρ = −0.149, p = 0.17 (not significant).

So how concentrated a pool is tells us *which ligands we can hope to get right*, at real
significance. That is an uncertainty estimate, not a selector — it cannot pick the pose,
but it can say when the pool is worth trusting, and it is the natural gate for a
tail-rescue or abstention strategy.

---

## What to do next, in order

1. **Build orientation-aware features.** The failure is diagnosed: everything tested is
   anchor-local. Needed instead: where the ligand's *substituents* sit relative to the
   pocket. Candidates — per-residue contact fingerprint against the CYP3A4 pocket prior,
   rotation angle about the Fe–donor axis, buried fraction of each ring system, and the
   position of the predicted site of metabolism relative to the iron.
2. **Add a second engine and check cross-engine agreement.** Consensus within one engine
   is weak (p=0.087). Cross-engine agreement is a different and historically stronger
   signal, and Chai-1 is already wired up.
3. **Then tier-1 QM.** Its per-atom donor and reactivity terms are also anchor-local, so
   on this evidence they should be expected to behave like the others. Worth running,
   but the expectation should be set low in advance rather than after the fact.
4. **Do not spend more GPU on sampling.** The oracle is 0.6975 and the PXR winning entry
   was 0.564. There is nothing wrong with the poses we have.
