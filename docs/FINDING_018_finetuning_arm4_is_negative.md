# FINDING 018 — Fine-tuning Boltz-2 on arm4_mix is a small, consistent NEGATIVE

**Date:** 2026-09-20 · **Status:** measured, one arm, one dose · **Verdict so far:** negative

## What was run

Boltz-2 (`boltz2_conf.ckpt`, 506.7M params) fine-tuned on **arm4_mix** — the
novelty-matched arm, held-out NN Tanimoto 0.536 against the challenge's measured 0.587
(FINDING 017). 302 training structures after validation, 85 held out by whole Murcko
scaffold cluster.

350 steps, lr 3e-4 (0.3× the pretrain 1e-3), warmup 50, accumulate 4, 384-token crops
centred on the ligand chain, `diffusion_multiplicity` 1. Trained 3:30:09 on one H200.

Weights genuinely moved — this is not a null run:

    global_step 350
    structure_module   rel = 2.766e-02
    pairformer_module  rel = 1.941e-02
    msa_module         rel = 1.828e-02

Both checkpoints then predicted **the same 85 held-out pairs from the same 85 YAML
files**, 5 diffusion samples each, and were scored with `cypstruct.pose` — the same code
behind every pool number in this repo.

## The result

| metric | base | fine-tuned | Δ | ft wins | Wilcoxon p |
|---|---|---|---|---|---|
| LDDT-PLI, sample 0 | 0.3477 | 0.3352 | **−0.0125** | 15/85 | <1e-4 |
| LDDT-PLI, best of 5 | 0.3655 | 0.3420 | **−0.0235** | 15/85 | <1e-4 |
| BiSyRMSD, sample 0 | 5.85 Å | 6.06 Å | **+0.21 Å** | 20/85 | <1e-4 |

Median per-pair ΔLDDT-PLI is **−0.0007**. The effect is *tiny per pose* and *consistent in
sign*, which is what produces a significant p on a negligible magnitude.

## Why it is negative, and what it is not

The bimodal structure of the baseline is **completely unchanged**:

| | base | fine-tuned |
|---|---|---|
| catastrophes (LDDT-PLI < 0.1) | 39 / 85 | 38 / 85 |
| successes (LDDT-PLI > 0.5) | 23 / 85 | 23 / 85 |
| sub-2 Å poses | 44 / 85 | 42 / 85 |

That is the whole finding. Roughly half of these predictions are essentially correct and
half are catastrophic, and **fine-tuning moved zero poses between those two modes**. It
nudged nearly every pose very slightly in the wrong direction and left the failure
structure exactly where it was.

A 2.8% relative change in the structure trunk that reclassifies nothing is the signature
of a systematic perturbation, not of learning. Compare the PXR campaign, which spent
twelve days and five platforms on fine-tuning and returned **−0.0020**. That number now
looks less like bad luck and more like the same effect.

## What this does not yet establish

One arm, one dose. Before "fine-tuning does not work here" is safe to say:

- **Dose.** Checkpoints at 87 / 174 / 261 / 348 steps exist. If 87 is also negative and
  348 is more negative, the damage is monotonic and the lever is simply wrong. If there
  is an optimum, the recipe is wrong rather than the idea. *Running now.*
- **Gradient variance.** `diffusion_multiplicity` is forced to 1 by a released-code
  constraint (see FINETUNE_PLAN addendum) where pretraining used 32. Single-sample
  diffusion gradients on 302 examples may be too noisy to learn from at any dose.
- **The other three arms.** arm4_mix is the arm whose held-out number transfers; the
  others would test whether *composition* matters, but there is no reason to spend 8 h
  each on them while the mechanism on the best-matched arm is a null.

## Cost and comparison

~4.5 h GPU for training plus two 43-minute inference passes. Absolute LDDT-PLI is low for
both arms (no Cys-SG→heme-FE bond constraint, 85 different P450 targets rather than
CYP3A4, 5 samples) — but every one of those applies identically to both, which is the
point of a paired design.

Against the alternative: **cross-engine consensus selection gains +0.0381** on CYP3A4 and
+0.0357 over 81 held-out proteins, with no training, no GPU, and no fitted parameters
(FINDING 011). Fine-tuning has now cost more and returned less than nothing.
