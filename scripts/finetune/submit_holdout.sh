#!/bin/bash
#SBATCH --job-name=cypho
#SBATCH --partition=gyorilab
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/shenoy.am/cyp-finetune/runs/%x_%j.log
#
# Predict one arm's held-out pairs from ONE checkpoint. Run it twice - base and
# fine-tuned - over the same YAML directory, so the only difference between the two
# result sets is the weights.
#
#   sbatch submit_holdout.sh arm4_mix base
#   sbatch submit_holdout.sh arm4_mix /scratch/shenoy.am/cyp-finetune/runs/arm4_mix/last.ckpt
set -euo pipefail
ARM=${1:-arm4_mix}
CKPT=${2:-base}
cd /scratch/shenoy.am/cyp-finetune

if [ "$CKPT" = "base" ]; then
    CKPT=boltz_cache/boltz2_conf.ckpt
    TAG=base
else
    TAG=ft
fi

OUT=holdout_out/${ARM}_${TAG}
mkdir -p "$OUT"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
# --cache points at the local boltz_cache so nothing reaches for the network: GPU nodes
# here have no direct internet, and a silent download attempt just hangs the job.
srun ./env/bin/boltz predict "holdout_yaml/${ARM}" \
    --checkpoint "$CKPT" \
    --cache boltz_cache \
    --out_dir "$OUT" \
    --output_format mmcif \
    --diffusion_samples 5 \
    --num_workers 4 \
    --override
