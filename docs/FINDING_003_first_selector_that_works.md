# Finding 003 — the first selector that beats random, and why it works

**Date:** 2026-09-11 · **Run:** `val87b`, unsteered arm · **n = 87 ligands, 1,740 poses**
**Data:** `orientation_selector_val87b_unsteered.json`, `consensus_features_*.csv`

FINDING 002 left 0.127 LDDT-PLI per ligand unclaimed and diagnosed why nothing could
claim it: every feature tried was **anchor-local**, describing the ligand atom nearest the
iron, which barely varies within a ligand. The prescription was to build features
describing **where the rest of the molecule sits**. That worked.

---

## The selector

```
score = 0.5 · z(n_pocket_residues_touched)  −  z(mean_rmsd_to_other_samples)
```

Both terms are z-scored **within the ligand**, which is the only comparison selection
makes. Neither term uses the reference structure, so nothing here leaks.

| | Δ vs random | beats random on | p |
|---|---|---|---|
| Boltz confidence (incumbent) | −0.0063 | 46.0% | 0.25 |
| contacts alone | +0.0108 | 59.8% | 0.045 |
| coordinate medoid alone | +0.0100 | 60.9% | 0.087 |
| **combined, fitted on all data** | **+0.0279** | **70.1%** | **0.0002** |
| **combined, leave-one-cluster-out** | **+0.0220** | **66.7%** | **0.0038** |

The honest number is the last row: **+0.0220**, with the weight chosen on the other
scaffold clusters and evaluated on held-out ligands. It recovers about **17% of the
oracle gap** and beats the incumbent by **+0.028**.

## Why it works: the two signals are orthogonal

Mean within-ligand correlation between the two terms is **+0.030** — essentially
independent. That is the whole mechanism. Each alone sits just at or below significance;
together they clear it comfortably, which is what genuinely orthogonal signals do.

- **`n_pocket_residues_touched`** is orientation-aware by construction. It counts how many
  CYP3A4 pocket residues the ligand contacts, so a pose rotated about the iron anchor
  scores differently even though its anchor atom has barely moved. This is precisely the
  axis FINDING 002 identified as carrying the information.
- **coordinate medoid** asks a different question entirely — not "is this pose good" but
  "do independent samples agree on it".

## It is not a knife-edge

| weight on contacts | 0.0 | 0.2 | 0.3 | **0.5** | 0.7 | 1.0 | 2.0 |
|---|---|---|---|---|---|---|---|
| Δ | +0.0100 | +0.0190 | +0.0211 | **+0.0279** | +0.0234 | +0.0188 | +0.0137 |
| p | 0.087 | 0.013 | 0.006 | **0.0002** | 0.0008 | 0.014 | 0.057 |

A smooth plateau, p < 0.01 across 0.3–0.8. Per-fold weights chosen independently ranged
only 0.4–0.6, median 0.5.

## What it is not

- **It is not large.** +0.022 against an oracle gap of 0.127 leaves 83% still on the table.
- **It does not rescue the anchor terms.** Coordination, trans angle, clash and the
  per-atom QM descriptors remain useless for ranking, for the reason in the FINDING 002
  addendum: they are near-constant within a ligand.
- **87 ligands is 33 scaffold clusters.** The effective sample size is the cluster count,
  and p = 0.0038 should be read against that, not against 87.

## Next

1. **More orientation features in the same family.** Contact fingerprints per residue
   rather than a bare count; contact-pattern consensus (the Jaccard version tested here
   was weak at +0.0028, worth revisiting with a better similarity); azimuth clustering
   about the Fe–donor axis.
2. **Cross-engine agreement.** Within-engine consensus contributes half of this result.
   Chai-1 now runs correctly, and cross-engine agreement is a different, historically
   stronger signal.
3. **Do not spend GPU on more sampling.** The oracle is already 0.6975 against a PXR
   winning entry of 0.564.
