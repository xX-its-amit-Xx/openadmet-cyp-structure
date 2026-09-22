# FINDING 024 — CYP3A4 fails by ROTATION in a pocket the model builds correctly

**Date:** 2026-09-22 · **Status:** measured, no new inference · **Confidence:** high on the
decomposition (n = 425 held-out poses + 1,660 CYP3A4 pool poses + 406 crystals), medium on
the cavity claim (the family comparison is n = 13 targets)

FINDING 022 established the gap. This asks why, tests six candidate explanations against
data already on disk, and kills three of them along with two of the three experiments that
were proposed to follow.

**The one-sentence answer.** Boltz-2 builds CYP3A4's active site correctly to 0.73 Å at
the pocket lining, drops the ligand within 1 Å of its true centre, and then turns it
**30° the wrong way** — and the protein's accuracy does not predict the ligand's
(ρ = +0.03, p = 0.8) where in every other P450 it does (ρ = +0.51, p = 1e-24). Everything
that conditions the *protein* — templates, a partner chain, orthologs — is therefore
aimed at a part of the problem that is already solved.

## Where the measurements came from

| what | source | n |
|---|---|---|
| held-out scores, 5 samples/pair, all targets | `data/processed/finetune/scores_arm4_mix_base.json` | 85 pairs / 425 poses |
| CYP3A4 deep pool, 4 seeds × 5 samples | `/scratch/shenoy.am/zexp/truth_run5_s4*.csv` | 87 pairs / 1,740 poses |
| cavity volume, ligand geometry | `why3a4.py` → `zexp/why3a4/cavity.csv` | 406 crystals |
| error decomposition + protein accuracy | `why3a4c.py` → `zexp/why3a4/decomp.csv` | 425 poses |
| crystal-to-crystal variation, b5 reach | `why3a4d.py` → `zexp/why3a4/crystal_variation.csv`, `b5_reach.csv` | 406 crystals / 25 CYP3A4 |
| MSA depth | `grep -c "^>"` over `cyp-finetune/msa/*.a3m` | 185 targets |

Scripts are in `/scratch/shenoy.am/zexp/why3a4{,b,c,d}.py`; all ran under `sbatch`.
`why3a4b.py` is **superseded** — it chose the residue offset by maximising overlap count
and produced 16.8 Å CA RMSDs. `why3a4c.py` uses the scorer's own residue-**name**-agreement
offset and is the one to trust. That is FINDING 021's trap, hit again, in one session.

---

## Part 1 — what died

### Ligand size: real for the family, **null inside CYP3A4**

Family-wide it is a genuine driver: ρ(n_heavy, LDDT-PLI) = **−0.405, p = 5e-4** across the
70 non-CYP3A4 pairs. Inside CYP3A4 it is absent, and if anything reversed:

| set | n | ρ(n_heavy, score) | p |
|---|---|---|---|
| non-CYP3A4 held-out | 70 | **−0.405** | 0.0005 |
| CYP3A4 held-out | 15 | +0.061 | 0.83 |
| CYP3A4 deep pool | 83 | **+0.243** | 0.027 |

**Bigger CYP3A4 ligands score better.** Caffeine (14 heavy atoms, 8SO1/8SO2) scores
0.487/0.539, while 14-atom ligands in other P450s average 0.911. Ritonavir at 50 heavy
atoms scores 0.751/0.764.

Size-matched against the 11 non-CYP3A4 pairs with ≥ 31 heavy atoms (CYP3A5 excluded, it is
the same phenomenon — see below), CYP3A4 still fails:

| | poses | raw RMSD | conformer | **rotation** | translation | pocket CA |
|---|---|---|---|---|---|---|
| size-matched others | 55 | 0.96 Å | 0.81 Å | **6.1°** | 0.30 Å | 0.51 Å |
| **CYP3A4** | 75 | **2.85 Å** | 1.60 Å | **30.1°** | 0.94 Å | 0.73 Å |
| CYP3A5 (P20815) | 10 | 3.05 Å | 2.40 Å | 23.5° | 1.03 Å | 1.56 Å |

Mann-Whitney CYP3A4 vs size-matched others: rotation **p = 4.2e-09**, conformer p = 5.5e-04,
pocket CA p = 1.9e-03. With size held fixed the rotation error is **5× worse** and the
protein error only 1.4× worse. **Size is not the cause.**

### MSA depth: dead, and in the wrong direction

CYP3A4's MSA is **14,571 sequences** — the second deepest of any test target. Within the
non-CYP3A4 pairs, depth vs score is ρ = **+0.031, p = 0.8**. In a joint regression it
carries p = 0.28. CYP3A4 is not alignment-starved.

### Template availability / training exposure: dead, and in the wrong direction

CYP3A4 has **107 harvested PDB entries**, more than twice any other target, and is the
worst performer. Per-target ρ(n_entries, score) = **−0.138, p = 0.65**. Nor is there a
memorisation rescue for the older entries: within CYP3A4, ρ(PDB-ID era, score) = **−0.43**
(p = 0.11) — the 2-, 3- and 4-series entries are no better than the 8-series. 4D78 (2014)
scores 0.307.

### The decoy sub-pocket: does not exist

Across 1,660 CYP3A4 pool poses, predicted minus crystal ligand centroid in the heme frame:

- median |offset| **1.07 Å**
- mean signed offset **(+0.38, +0.15, −0.16) Å**, sd (1.20, 1.03, 1.06) — |mean|/sd =
  0.32 / 0.14 / 0.15 per axis, i.e. **isotropic scatter about the true centre, no direction**
- predicted centroid height above the heme plane 7.16 Å vs crystal 7.08 Å (Δ = **−0.05 Å**)
- failed pairs only (median BiSyRMSD > 2 Å): mean Δz = −0.19 Å, mean radial offset 1.38 Å

The ligand lands in the right cavity, at the right height, on the right point, and is
**turned the wrong way**. Clustering the failures would cluster noise.

### "F/G remodelling" as the proximate cause of *our* numbers: not supported

OpenADMET's explanation may be right about the physics, but it is not what limits this
pipeline. Three measurements:

1. **The protein is right.** Global CA 0.86 Å; pocket-lining CA **0.73 Å** (p90 1.93 Å;
   only 15% of poses exceed 1.5 Å).
2. **The F/G span is not missing from the crystals.** Residues 202–260 are modelled in a
   median of **56 of 59** positions across 25 CYP3A4 entries; 17 of 25 are ≥ 90% complete.
3. **Decisive:** within CYP3A4, protein accuracy does not predict ligand accuracy.

| | ρ(pocket CA, ligand RMSD) | p |
|---|---|---|
| every other P450 | **+0.510** | 1.4e-24 |
| **CYP3A4** | **+0.030** | 0.80 |
| CYP3A4, global CA instead | +0.025 | 0.83 |

Elsewhere, a wrong ligand comes with a wrong pocket — one coupled failure. On CYP3A4 the
model builds the pocket correctly and *still* cannot orient the ligand in it. **Fixing the
protein cannot fix the ligand here, and this is the measurement that says so.**

---

## Part 2 — what survived

### Pocket volume: survives as the family-level driver, insufficient on its own

Cavity volume measured in the heme frame with the ligand removed, as the
**ligand-centre-accessible** volume: grid points > 3.1 Å from every protein/heme heavy
atom, buried along ≥ 80% of 26 rays, flood-filled from the column above the iron. These are
*not* comparable to published solvent-excluded Å³ figures (which are 3–5× larger); only the
ranking is meaningful.

| | entries | median cavity |
|---|---|---|
| **CYP3A4** | 107 | **881 Å³** (p10 513, p90 1130) |
| every other P450 | 299 | **262 Å³** (p90 657) |
| Q2IU02 (best target, 0.924) | 48 | **50 Å³** |
| P11511 aromatase (0.631) | 7 | 117 Å³ |
| P20815 CYP3A5 (0.506) | 5 | 733 Å³ |

CYP3A4 ranks **3rd of 50** targets by median cavity. Volume predicts score family-wide:
ρ = **−0.514, p = 5e-07** over 85 pairs, and **−0.385, p = 0.001** within the non-CYP3A4
pairs alone, so it is not simply the CYP3A4 rows driving it.

**But it is not sufficient.** Five non-CYP3A4 pairs with cavities of 635–885 Å³ — as large
as CYP3A4's — average **0.851** (3R9B 0.836, 7CL9 0.838, 2JJO 0.859, 6Q2T 0.816,
3ZKP 0.906). And inside CYP3A4, fill fraction does not predict the score (ρ = +0.075,
p = 0.5), though the sub-2 Å *rate* does climb across fill terciles (0.14 → 0.30 → 0.39,
n = 83). A big cavity is necessary-ish, not sufficient.

Joint OLS on the 85 held-out pairs, LDDT-PLI ~ log(cavity) + n_heavy + MSA depth + CYP3A4:

| model | R² | CYP3A4 coefficient |
|---|---|---|
| CYP3A4 alone | 0.245 | **−0.294** (p = 1.5e-06) |
| + log cavity + n_heavy | 0.273 | **−0.246** (p = 1.9e-04) |
| + MSA depth | 0.283 | **−0.203** (p = 0.008) |

Every measurable covariate together absorbs **31%** of the raw gap. **0.20 LDDT-PLI is not
explained by cavity, size or alignment depth.**

### It is the CYP3A subfamily, not CYP3A4 alone

CYP3A5 (P20815, 84% identity to CYP3A4, 733 Å³ cavity) scores **0.506** with raw 3.05 Å,
rotation 23.5° and translation 1.03 Å — the same signature. n = 2 pairs, so this is a
pointer, not a result. But it means the analogue-set problem of FINDING 022 has a candidate
answer that is leak-free and already on disk.

### Pose dispersion: dispersed, but dispersed inside a wrong basin

Two facts that look contradictory and are not.

**Within a single 5-sample job**, CYP3A4's poses are 3.5× more dispersed than the family:
median pairwise ligand RMSD in the heme frame **2.02 Å vs 0.58 Å** (Mann-Whitney
p = 0.0011). But dispersion scales with error everywhere — disp/err is **0.58** for CYP3A4
and **0.62** for the other targets' failures. CYP3A4 is not distinctively "confidently
wrong" at the pose level; it has the ordinary failure profile, larger.

**Across independent seeds**, which ligands are hard is deterministic:

- seed-to-seed ρ of the per-pair mean score: **+0.941 to +0.966** (all six pairs of seeds)
- **96.2%** of the variance is between ligands, 3.8% between seeds
- **44 of 87 pairs never produce a sub-2 Å pose in 20 tries**; only 6 of 87 are uniformly
  correct
- 20-pose oracle **0.693** against a per-pose mean of 0.576 — and a family mean of 0.865
  at only 5 samples

So the model scatters within a ligand-specific basin that is, for half the ligands, the
wrong basin. **Sampling harder does not reach the answer**; four seeds bought +0.117 of
oracle and left 0.17 to the family baseline.

### The error decomposition — the headline

Each pose superposed on the binding site, then the ligand optimally superposed on itself
under the scorer's own symmetry-aware atom mapping (`best_ligand_mapping`, 100% success):

| | raw | internal conformer | rotation | translation |
|---|---|---|---|---|
| every other P450 | 0.66 Å | 0.34 Å | 6.1° | 0.24 Å |
| **CYP3A4** | **2.85 Å** | **1.60 Å** | **30.1°** | **0.94 Å** |

74% of CYP3A4's mean-square error is rigid-body, 26% internal conformer. The rotation-error
distribution is a **continuum, not a flip**: 31% under 15°, 33% at 15–45°, 16% at 45–90°,
15% at 90–150°, 5% over 150° (the others: 74% under 15°). There is no 180° mode, so this is
not a head-to-tail ambiguity that a symmetry fix would catch. Within CYP3A4, what predicts
the pose error is rotation (ρ = +0.795), translation (ρ = +0.849) and conformer
(ρ = +0.595) — and *not* the protein (ρ = +0.03).

**Reading:** the cavity is large enough that the correct orientation is not determined by
sterics, and the model has no other term that determines it.

---

## Part 3 — the experiments that follow, and the two that do not

### ✗ Cytochrome b5 as a co-folded partner — DO NOT SPEND

Three independent reasons, in order of force.

1. **It targets a part of the problem that is already solved.** b5 could only act by
   stabilising the protein. The protein is correct to 0.73 Å at the pocket and its accuracy
   is uncorrelated with the ligand error (ρ = +0.03, p = 0.80). There is no measured route
   from a better protein to a better pose on this target.
2. **The interface is on the far side of the heme.** Measured on 25 CYP3A4 crystals, the
   12 published proximal-face interface positions (K96, K115, K121, K127, K130, K141, K421,
   R422, K424, K428, K440, R446) sit **9.7 Å from the iron, 8.7 Å from the nearest ligand
   heavy atom, and 10.7 Å from the nearest modelled F/G residue**, across the porphyrin.
   Ligands are 100% distal (`frac_proximal` = 0.0 across the harvest). A CA–CA gap of
   ~11 Å to the F/G span is not so large that allostery is *impossible*, so I will not
   claim the mechanism is physically excluded — only that it is indirect and unevidenced,
   while the direct route is measured to be closed.
3. **There is no ground truth to fold against.** No CYP3A4–b5 complex is deposited. The
   co-folded interface would itself be an unvalidated prediction, and FINDING 009's lesson
   applies: an extra chain adds degrees of freedom, not opinions.

### ✗ Holo templates — DO NOT SPEND, and here is the killer

Mechanically it works. Boltz-2 **2.2.1** (the version in
`/scratch/shenoy.am/cyp-finetune/env`, unmodified) accepts

```yaml
templates:
  - cif: /path/to/1TQN.cif
    chain_id: A          # optional; a QUERY protein chain, never the ligand chain
    template_id: A       # optional; a chain of the template
    force: false         # true requires `threshold`
```

Two facts from reading the code, both fatal:

1. **A template carries no ligand and no side chains.** `schema.py` restricts template
   chains to `mol_type == PROTEIN` and raises if a `chain_id` is not a query protein chain,
   so the ligand chain can never be templated. `featurizerv2.compute_template_features`
   emits only `template_restype`, `template_ca`, `template_cb` and
   `template_frame_rot/frame_t` — **backbone frames and CB**. The Phe-roof rotamers
   (Phe57/108/213/215/219/220/241/304, with Phe304 the gatekeeper) that could break the
   rotational degeneracy are exactly what a template cannot carry.
2. **There is no right template to pick.** Pairwise over 30 CYP3A4 crystals (435 pairs),
   the pocket lining varies by a median of **1.00 Å CA RMSD (p90 1.46 Å)** — the largest of
   any target measured, against 0.26 Å for Q2IU02 and 0.12 Å for aromatase:

   | uniprot | structures | global CA | **pocket CA** |
   |---|---|---|---|
   | **P08684 (CYP3A4)** | 30 | 1.19 Å | **1.00 Å** |
   | Q7Z1V1 | 4 | 1.72 Å | 0.90 Å |
   | Q16850 | 5 | 0.65 Å | 0.86 Å |
   | Q2IU02 | 30 | 0.28 Å | 0.26 Å |
   | P11511 (aromatase) | 7 | 0.18 Å | 0.12 Å |

   **The model's prediction (0.73 Å from the query crystal) is already closer than a
   CYP3A4 crystal chosen without knowing the answer (1.00 Å).** A template selected
   blind is, on average, a *worse* prior than what Boltz already produces. This is §4.2 of
   `docs/worldmodel/CYP3A4_EVOLUTION.md` — "the pocket that would orient the ligand is
   created by the ligand" — turned into a number.

**If it is run anyway, the non-circular protocol is:** for the *challenge*, the references
are new cryoEM structures, so any deposited X-ray template is legitimate and unrestricted.
For the **87-ligand proxy it is circular** unless the template is forbidden by construction
from being the query's own entry *or any entry in the query's scaffold cluster* — the
`scaffold` column in `arm4_mix.csv` is the handle, and ritonavir alone spans 3NXU/5VC0/5VCE.
The only always-legal template is apo (1TQN, 1W0E). Report cluster-held-out and challenge
numbers separately, never pooled. Given the 1.00 Å crystal-to-crystal spread, my prediction
is that apo templating is **negative** and cluster-held-out holo templating is a null.

### ✗ Animal orthologs — DEAD, no data exists

The entire CYP3A subfamily in the 493-pair harvest is **human**: CYP3A4 107 entries,
CYP3A5 5, CYP3A7 1, CYP3A43 0. **Zero** ligand-bound structures for rat 3A1/3A2, dog 3A12,
or any macaque CYP3A. So orthologs cannot supply templates or training pairs, and the only
remaining channel — MSA content — is the one axis already measured to be irrelevant
(14,571 sequences, ρ = +0.03 with score). There is nothing to run.

---

## Part 4 — what to run instead, ranked

### R1. Split the cross-engine score into translation and rotation. Zero GPU.

> **RUN, AND CLOSED — see FINDING 026.** The orientation-only term is the strongest single
> unfitted term measured in this project (+0.0424 vs random, ρ = −0.241, correct on 75.9%
> of ligands) and it does **not** beat the incumbent paired: +0.0041, Wilcoxon p = 0.61,
> Holm p = 1.00, same pose selected on 49 of 87 ligands. The mechanism claim below is
> refuted: within-ligand ρ(mix, ori) = **+0.879** against ρ(mix, cen) = +0.667, so the
> Chamfer was already ranking by orientation, and the centroid term is not a passenger —
> `cen` alone scores **+0.0309**, 81% of the incumbent's whole gain. Consensus agreement
> and error-against-truth are different quantities.

**Mechanism, measured:** the incumbent `cypstruct.xengine.xeng_score` is a Chamfer distance
in the heme frame, so it is dominated by *where* the ligand is. On CYP3A4 the centroid is
already right to 1.07 Å with no systematic direction, while the rotation is wrong by 30°.
The incumbent is spending most of its dynamic range on a coordinate that is not the one that
is wrong.

**Do:** for each pose and each reference of the same ligand, recentre both on their own
centroids before the Chamfer, giving an orientation-only agreement term; keep the centroid
offset as a separate feature. Test the orientation term alone, the pair, and the incumbent,
on the 1,660 poses in `zexp/out/run5_s4*` against the frozen reference set
(`cyp-finetune/reference_set_cyp3a4.npz`) — everything is already on `/scratch`.

**Gate:** leave-one-ligand-out, beat the incumbent's **+0.0383** on the same pool. Noise
floor is +0.0139 (p95) / +0.0196 (p99); measure the null on shuffled references first.
**Cost: one CPU `sbatch`, no inference.** This is the only proposal here whose mechanism is
a direct consequence of the decomposition.

**Honest caveat:** FINDING 006 already killed "azimuth consensus". That measured the
azimuth about the Fe axis; this measures the full 3D rotation about the ligand's own
centroid, which is a different quantity for a ligand whose centre sits 7.1 Å above the iron
and whose long axis is not radial. If it reproduces 006's null, orientation consensus is
closed and should be written up as closed.

### R2. Rebuild the analogue set on cavity volume. One small inference job.

**Mechanism, measured:** FINDING 022 asked for analogues chosen by difficulty rather than
family and had no criterion. There is now one that is leak-free (computable from any
predicted structure), family-independent, and validated at ρ = −0.385 within the
non-CYP3A4 pairs: **cavity volume ≥ 600 Å³**. That selects CYP3A5 (733), Q16850 (760),
Q82GL5 (720), P48635 (994) and the large tail of Q7Z1V1 (to 1054) out of the 493-pair
universe — roughly 40 pairs, and it includes the one target that already reproduces
CYP3A4's exact failure signature.

**Gate:** the selected set must show a sub-2 Å rate near CYP3A4's 33%, not the family's 89%.
If it does not, difficulty is not cavity-shaped and FINDING 022's programme needs a
different criterion. Either outcome is worth the job.

### R3. A second copy of the same ligand, small scope, weak evidence — declare it a probe

**Mechanism:** §4.4 of the worldmodel — one substrate in a cavity that big rattles and is
catalytically uncoupled; 8SO1/8SO2 show **two** caffeines and 8GK3 four DHEA-S. Boltz-2
accepts two ligand chains with the same SMILES; score only the copy nearest the iron.

**Why it is ranked last and not recommended outright:** the supporting correlation is weak
and partly null. Sub-2 Å rate rises across fill terciles (0.14 → 0.30 → 0.39, n = 83) but
mean LDDT-PLI does not (ρ(fill, score) = +0.075, p = 0.5), and multiple occupancy is rare
in the harvest itself — **2 of 126** CYP3A4 entries carry more than one copy of the query
ligand in a chain. Run it, if at all, on the 8 smallest-ligand CYP3A4 pairs only, as a
probe with a pre-registered null, not as a pool change.

---

## What this changes

- **FINDING 022's "CYP3A4 is the outlier" is now "the CYP3A subfamily's cavity is the
  outlier, and the failure is rotational."** Report the decomposition, not just the score.
- **Stop proposing protein-side fixes for this target.** ρ(pocket CA, ligand RMSD) = +0.03
  is a general veto: templates, partner chains, orthologs, more MSA — all condition the
  protein, and the protein is not what is wrong.
- **The 0.865 family baseline is inflated.** 27 of the 70 non-CYP3A4 pairs are one bacterial
  target (Q2IU02) with 12-heavy-atom ligands in a 50 Å³ pocket. The unweighted per-target
  mean excluding CYP3A4 is **0.798** (12 labelled targets), not 0.865.
- **`why3a4b.py` is a worked example of FINDING 021's trap** and is kept in `/scratch`
  unfixed, next to `why3a4c.py`, for exactly that reason.
