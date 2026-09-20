# EXPERIMENT — ATOMICA as a CYP3A4 pose selector

**Date:** 2026-09-20 · **Verdict: it loses.** ATOMICA cannot rank CYP3A4 poses. Its best
feature scores **+0.0148** against random, below this repo's **+0.0138 noise floor**
(FINDING 007) and far below the shipped cross-engine selector's **+0.0395**. The
pre-registered ablation came back exactly as feared: on the *same* poses, ATOMICA separates
a deliberately rotated ligand from an unrotated one at **AUC 0.95–0.97**, while ranking real
prediction errors at chance. **It is a rigid-perturbation detector, not a pose-quality
scorer.**

---

## What was asked

`docs/worldmodel/MULTIMODAL_JEPA.md` §1.1 marked ATOMICA "✅ directly usable as a pose
scorer" — it embeds an already-bound interface, it ships an ATOMICA-Ligand **HEM**
checkpoint, and CYP3A4 is a heme protein. It also pre-registered the ablation, because
ATOMICA's pretraining objective *is* denoising rigid rotations and translations, and 86% of
its corpus is CSD small-molecule crystal packing. Both halves of that note are now measured.

## Setup

| item | value |
|---|---|
| pool | `data/processed/poses_scored_val87b.csv`, **unsteered arm**: 1,740 Boltz-2 poses, 87 CYP3A4 ligands, 20 samples each |
| pose files | `D:/cyp_scratch/val87b_unsteered/<LIG>__unsteered__s1/input_model_*.cif` (606 MB), staged to `/scratch/shenoy.am/atomica_exp/val87b_unsteered` |
| truth | per-pose `lddt_pli` already in the CSV; nothing was re-scored |
| random baseline | **0.5769** (pool mean) · oracle **0.6975** · incumbent `xeng` **0.6164** (+0.0395) · Boltz `confidence` **0.5706** (−0.0063) |
| null | 4,000 random-pick draws → 95th pct **+0.0141**, 99th pct **+0.0196** |
| venue | Explorer, `gyorilab` partition, CPU only. ~1 h wall for 6,786 structures. |
| env | `/scratch/shenoy.am/atomica_exp/env` — a **separate** uv venv (python 3.12, torch 2.11+cu128). `cyp-finetune/env` was not touched. |
| code | ATOMICA @ `mims-harvard/ATOMICA` HEAD, MIT. Weights `ada-f/ATOMICA`, CC-BY-4.0 → `/scratch/shenoy.am/atomica_exp/ckpt` |
| checkpoints | `pretrain_model_weights.pt` (34 MB, ~8M params) and `ligand/small_molecules/HEM/HEM_v1.pt` |

Two receptor definitions were processed for every pose, because this project's thesis is
that pose quality lives in the heme frame:

* **`prot`** — receptor = protein chain A, ligand = chain L / `LIG1`. ATOMICA's canonical
  protein–small-molecule use. 8 Å interface, `PS_300` fragmentation (matches the
  checkpoint's `fragmentation_method`).
* **`protheme`** — same, but the HEM cofactor's 43 atoms are appended to the receptor
  blocks before the interface cut, so the iron is inside the graph.

0 failures out of 1,740 in each mode. Typical graph: 35 blocks / 249 atoms (`prot`).

## Features scored (all label-free)

1. **`tr`, `ro` — the denoising heads read as an energy.** ATOMICA was pretrained by score
   matching on rigid translation / rotation noise, so on a clean input the predicted score
   should be ~0 for a native-looking complex. `score.py` reproduces
   `DenoisePretrainModel.forward` up to the scaled global heads (`top_translation_scale_ffn`,
   `top_rotation_scale_ffn`) with `perturbed_segment=1` (the ligand), and reports their
   norms. Torsion denoising is switched off (`encoder.encoder.remove_torsion_denoiser()`),
   which also removes the `torch_cluster` dependency.
2. **`medoid` — embedding-space consensus.** `z_graph` with
   `pool='mean_component_normalized'`, then per ligand the mean cosine distance to that
   ligand's other 19 poses. This is the ATOMICA analogue of the incumbent's Chamfer
   consensus, and the fair way to use a frozen embedding for selection.
3. **`hem` — the ATOMICA-Ligand HEM checkpoint.** Run on the protein-only side (segment 0)
   of each pose's interface, which is what that head expects.

## Result 1 — selection

Mean LDDT-PLI of the selected pose, one pose per ligand, n = 87.

| selector | selected | gain vs random | p vs null | within-ligand ρ | ρ sign correct |
|---|---|---|---|---|---|
| **oracle (best in pool)** | **0.6975** | +0.1206 | — | — | — |
| **incumbent `cypstruct.xengine`** | **0.6164** | **+0.0395** | <0.00025 | **−0.258** | 75.9% |
| `protheme_medoid` (best ATOMICA) | 0.5917 | +0.0148 | 0.040 | −0.165 | 71.3% |
| `prot_medoid` | 0.5884 | +0.0115 | 0.089 | −0.154 | 73.6% |
| `prot_ro` (max) | 0.5868 | +0.0099 | 0.123 | +0.089 | 57.5% |
| `protheme_tr` (max) | 0.5858 | +0.0089 | 0.152 | +0.078 | 57.5% |
| **random** | **0.5769** | 0 | — | — | — |
| Boltz `confidence` | 0.5706 | −0.0063 | 0.775 | −0.033 | 49.4% |
| `prot_tr` (min, the objective's direction) | 0.5549 | −0.0220 | 0.995 | +0.036 | 39.1% |
| `protheme_tr` (min) | 0.5436 | −0.0333 | 1.000 | +0.078 | 42.5% |

**The denoising heads are at chance or worse.** In the direction the pretraining objective
implies — low predicted score = native-like — they select *worse than random*, by as much as
−0.033. Flipping the sign to whatever the data prefers buys +0.009, under the noise floor,
and sign-flipping after seeing the answer is not a result.

**The embedding consensus is the only positive**, and it is the feature this repo already
has: plain geometric consensus measured **+0.0131** in FINDING 011. ATOMICA's embedding
reproduces that, does not improve on it, and does not reach +0.020.

**It does not combine.** z-scoring within ligand and adding to the incumbent:

| weight on `z(protheme_medoid)` | 0.25 | 0.5 | 1.0 |
|---|---|---|---|
| selected | 0.6169 | 0.6156 | 0.6054 |

+0.0005 at the smallest weight, monotonically negative after — the same absorbed-or-harmful
pattern FINDING 011 recorded for the cross-engine term against the old incumbent.

### The HEM checkpoint is a pocket classifier, not a pose scorer

`ATOMICA-Ligand/HEM_v1` was finetuned on **ligand-free protein binding sites** and asks
"does this pocket bind heme". It never sees the ligand; a pose can only move its output by
changing which residues the 8 Å cut selects. Measured anyway: within-ligand ρ **+0.023**,
selection **+0.0098**. Within-ligand sd (0.0224) is
*larger* than between-ligand sd (0.0174), i.e. most of its variation is pose-selection noise.
52.9% of ligands run in the right direction — a coin flip.
The "closest public checkpoint for a CYP project" is the wrong shape for this question.

## Result 2 — the pre-registered ablation, and it is the informative one

For each ligand, take the pool's **highest-LDDT-PLI pose** (mean 0.6975) as the near-native
reference, write it back through gemmi, and then write six rotated copies per axis × three
fixed pseudo-random axes: the LIG1 residue rigidly rotated about the **HEM iron** by
15/30/60/90/120/180°. Protein and heme untouched, so only ligand orientation differs. The
angle-0 reference goes through the identical gemmi round-trip, so no file-format confound.
1,653 structures, 0 processing failures.

`protheme`, AUC of the score separating rotated poses from their own unrotated reference:

| angle | 15° | 30° | 60° | 90° | 120° | 180° |
|---|---|---|---|---|---|---|
| `tr` AUC (rotated scores higher) | **0.951** | **0.969** | 0.931 | 0.813 | 0.751 | 0.739 |
| `ro` AUC | 0.185 | 0.076 | 0.048 | 0.042 | 0.039 | **0.040** |
| `ro` mean | 2.49 | 1.71 | 1.38 | 1.28 | 1.24 | 1.25 |

`ro` on the unrotated reference is **4.07 ± 1.62** and collapses monotonically to ~1.25 by
60°. Read as a separator (|AUC − 0.5| = 0.46), that is **AUC ≈ 0.96 at 30°** — a 30°
rotation of the ligand is detected almost perfectly. `prot` mode gives the same picture
(`ro` AUC 0.043 at 30°, i.e. 0.957 reversed; `tr` 0.858).

So, plainly:

> **ATOMICA separates a 30° rigid rotation of the ligand from the unrotated pose at
> AUC ≈ 0.96, and ranks the pool's real 0.00–0.96 LDDT-PLI spread at ρ ≈ 0.08.**

Two further details make this worse, not better:

* **`tr` changes sign between the two regimes.** On rigid decoys, rotated poses score
  *higher* (AUC 0.95). In the pool, higher `tr` goes with *better* poses (pooled ρ = +0.14;
  AUC separating pool poses with LDDT-PLI < 0.40 from > 0.70 is **0.378**, i.e. bad poses
  score *lower*). A feature whose direction flips between synthetic and real error is not a
  pose-quality feature.
* **The pool is not short of bad poses.** Its LDDT-PLI quantiles are
  0.00 / 0.244 / 0.446 / 0.607 / 0.722 / 0.825 / 0.959 (min, p5, p25, p50, p75, p95, max);
  338 poses are below 0.40. ATOMICA had genuine catastrophes to find and did not find them —
  on that bad-vs-good split the denoising heads run at AUC 0.35–0.42 (wrong side of 0.5),
  while the embedding medoid reaches 0.707, which is real but is again just consensus.

This is the outcome the dossier pre-registered: *"score native poses against deliberately
rotated ones; if the score separates those trivially, it is a crystallinity detector, not a
pose scorer."* It separates them trivially.

## Why, mechanistically

The errors in this pool are not rigid-body errors. FINDING 001-A: whole-protein CA fit is
~1.0 Å and the heme is placed correctly; FINDING 011/008: 84% of poses already sit inside
any reasonable Fe-coordination window. What varies between a 0.25 pose and an 0.85 pose is
**substituent placement around a correctly anchored core** — which ring flipped, which
halogen points at which residue. ATOMICA's global translation and rotation heads are, by
construction, blind to that: they predict a single rigid transform per graph. The one
channel that could see it, torsion denoising, needs `torch_cluster` and was disabled; it is
also a within-molecule term, not an interface term. So this is the fifth independent
confirmation of the repo's standing result — **the signal is in substituent placement, and
anchor-level or whole-body features cannot reach it.**

## Reproduce

```bash
ssh explorer
cd /scratch/shenoy.am/atomica_exp
./env/bin/python make_index.py                                   # 1,740 poses -> index_all.csv
sbatch run_prep.sh prot     index_all.csv all                    # 8 A interface graphs
sbatch run_prep.sh protheme index_all.csv all                    # ... with HEM in the receptor
sbatch run_score.sh proc_all_prot.pkl     scores_all_prot.npz    # tr, ro, z_graph, h_graph
sbatch run_score.sh proc_all_protheme.pkl scores_all_protheme.npz
sbatch run_hem.sh   proc_all_prot.pkl     scores_hem_all.npz     # ATOMICA-Ligand HEM_v1
srun  -p gyorilab -c 4 --mem=16G ./env/bin/python analyse.py \
      --truth truth.csv --xeng xeng_val87b.csv \
      --scores protheme=scores_all_protheme.npz --out results_protheme.json

# ablation
srun ./env/bin/python rotate.py --truth truth.csv --index index_all.csv \
     --outdir rotated --outindex index_rot.csv
sbatch run_prep.sh  prot index_rot.csv rot ; sbatch run_score.sh proc_rot_prot.pkl scores_rot_prot.npz
srun ./env/bin/python ablation.py --scores scores_rot_protheme.npz --index index_rot.csv \
     --tag protheme --out ablation_protheme.json
```

Scripts on Explorer: `prep.py`, `score.py`, `hem.py`, `rotate.py`, `analyse.py`,
`ablation.py`, `extra.py`. Nothing bulky was written to C: or D:; the 9.2 GB working set
lives on `/scratch`. **Do not reuse `/scratch/shenoy.am/cyp-finetune/env`** — this used a
separate venv, deliberately.

## Traps paid for here, for the next person

* **The login node silently SIGKILLs.** `conda create` and even `import scipy` were killed
  with no message. Everything — including two-second analysis snippets — must go through
  `srun`. This looked like a broken install three times before it was diagnosed.
* **`r.sample` on a pandas row is `DataFrame.sample`, the method**, not the column. It
  produced a filename containing a bound-method repr and a `FileNotFoundError` two frames
  away from the cause.
* **`remove_torsion_denoiser` lives on `InteractionModule`**, two levels down:
  `model.encoder.encoder.remove_torsion_denoiser()`.
* The ablation's angle-0 reference **must** go through the same writer as the rotated
  copies, or the "native vs decoy" AUC measures the file format.

## Standing

Logged as a **negative**, per the repo rule. ATOMICA is not added to the selector and
`cypstruct.xengine.select()` is unchanged. The dossier's §1.1 "✅ directly usable as a pose
scorer" should be read as falsified for this target: it is usable, it runs, and it does not
work. The wet-lab heme result in that paper stands; it is about *whether a protein binds
heme*, which is a different question from *where a ligand sits above one*.
