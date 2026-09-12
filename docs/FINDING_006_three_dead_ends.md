# Finding 006 — three more features rejected, and a cheap pre-screen for the next one

**Date:** 2026-09-12 · **Data:** `val87b` unsteered, 1,740 poses, 87 ligands
Written so none of these is tried again, and so the next candidate is screened in seconds
rather than an hour.

---

## The three, and why each failed

| feature | Δ vs random | p | effect on the FINDING 003 selector | why |
|---|---|---|---|---|
| contact-fingerprint consensus (Jaccard) | +0.0028 | 0.58 | degrades, +0.0279 → +0.0145 | **redundant**: ρ = +0.177 with the contact count already in the selector |
| azimuth about the Fe–donor axis, as circular consensus | −0.0078 | 0.57 | degrades → +0.0167 | **orthogonal but empty**: ρ = −0.019 with the selector, yet carries no signal |
| crystallographic contact prior (leave-one-ligand-out) | −0.0160 | 0.61 | degrades → −0.0073 | **too little within-ligand contrast** |

The azimuth one deserves a note: it was tested first as a raw scalar, which is meaningless
for a circular quantity where −179° and +179° are neighbours. Retested properly as angular
distance to the per-ligand circular mode, it still does nothing. The failure is real, not
an artefact of the first mistake.

The contact prior is the most interesting failure. It was built from 107 ligand-contact
sets across 106 deposited CYP3A4 entries, with a strict leave-one-ligand-out split so no
ligand votes on its own pose. The prior itself is *correct* — it recovers Ala305 (99%),
Thr309 (98%), Ser119 (97%), Phe304 (94%), Ile301 and Phe108, i.e. exactly the I-helix and
Phe-cluster residues the literature names. It is right and it is useless, because **almost
every predicted pose already touches those residues**.

## The cheap pre-screen

Within-ligand coefficient of variation, and the fraction of ligands where a feature is
effectively flat (CV < 0.05):

| feature | within-ligand CV | ~constant on | selects? |
|---|---|---|---|
| `n_pocket_residues_touched` | 0.079 | 9% | **yes**, +0.0108 |
| `mean_rmsd_to_others` | 0.072 | 13% | **yes**, +0.0100 |
| `prior_recall` | 0.070 | 33% | no |
| `contact_consensus` | 0.065 | 37% | no |
| `is_coordinated` | 0.023 | **98%** | no |

**Flatness is a necessary condition, cheap to check, and it disqualifies decisively at the
extreme.** `is_coordinated` is constant on 98% of ligands — every sample already
coordinates — so it cannot rank anything, which is why its +0.330 between-class effect at
p = 1e-9 never converted (see the FINDING 002 addendum).

**But it is not sufficient**, and saying otherwise would overclaim. `prior_recall` has
almost the same spread as the two features that work, and still fails: it varies, but not
in a direction that tracks quality.

So the screen is: compute within-ligand CV first, discard anything below ~0.03 immediately,
and accept that everything above it still needs the full test.

## The pattern across all failures so far

Every rejected feature answers **"is this pose in the right binding mode?"** — does it
coordinate, does it touch the usual residues, does it sit where the other samples sit.
On this target that question is already settled before selection begins: the engines place
the ligand in the right pocket, in the right general mode, almost every time.

What remains undiscriminated is *which near-miss is nearest*, and none of the
mode-level descriptors speak to it.

## What is left to try

1. **Per-residue contact *deviation***, not recall — which residues does this pose touch
   that its sibling samples do not. The contrast is between poses, not against a prior.
2. **Sub-pocket occupancy**: CYP3A4's cavity has distinguishable lobes; which one the bulk
   of the ligand occupies may vary between near-misses where total contact count does not.
3. **Physics with real within-ligand contrast** — the tier-2 xTB interaction energy on a
   pocket cutout, which was specified in `docs/QM_SCORER_DESIGN.md` and never run. It is
   the one remaining candidate that is not a mode-level descriptor.

---

## Addendum, same day: contact *deviation* also fails, and the pre-screen is vindicated

Candidate 1 from the list below — per-residue contact deviation between sibling poses —
was the next tried. Three variants, all contrasts between poses of the same ligand rather
than against any prior:

| feature | within-ligand CV | Δ vs random | p |
|---|---|---|---|
| `missing_common` (common contacts this pose lacks) | **1.69** | +0.0168 | 0.15 |
| `rare_contact_count` (idiosyncratic contacts) | **1.49** | +0.0136 | 0.28 |
| `weighted_agreement` (entropy-weighted sibling agreement) | 0.70 | +0.0056 | 0.55 |

**The pre-screen worked exactly as intended.** These have within-ligand CV of 1.5–1.7,
twenty times the spread of anything previously tested, so they sail past the flatness
filter — confirming it is a genuine necessary-condition check and not a proxy for
"features I happened to like".

**And they still fail.** `missing_common` at +0.0168 is the largest single-feature effect
measured so far, bigger than either term inside the working selector, yet p = 0.15 and it
*degrades* the selector at every weight from +0.0279 down to +0.0125.

The reason is redundancy, and it is measurable: **ρ = +0.357 between `missing_common` and
the contact count** already in the selector. A pose that lacks the commonly-touched
residues necessarily touches fewer residues overall, so the two terms are two views of one
quantity. Compare the +0.030 between the selector's own two terms — that is what
independence looks like, and it is why they compose.

Running tally: **eight features tested, eight rejected.** The selector is unchanged at
+0.0279 on all data, +0.0220 on held-out clusters.

Remaining untried: sub-pocket lobe occupancy, and the tier-2 xTB interaction energy — the
only candidate left that is neither mode-level nor a re-expression of contact count.
