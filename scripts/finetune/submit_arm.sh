#!/bin/bash
#SBATCH --job-name=cypft
#SBATCH --partition=gyorilab
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/shenoy.am/cyp-finetune/runs/%x_%j.log
#
# Submitted, not srun'd through an ssh session: a multi-hour run must not depend on a
# laptop staying connected. Poll with `sacct -j <id>` and the run's metrics.csv.
#
#   sbatch submit_arm.sh arm4_mix 350
set -euo pipefail
ARM=${1:-arm4_mix}
STEPS=${2:-350}
cd /scratch/shenoy.am/cyp-finetune
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
srun ./env/bin/python train_arm.py \
    --arm "$ARM" --steps "$STEPS" \
    --lr 3e-4 --warmup 50 --accum 4 \
    --workers 4 --max-tokens 384 --max-atoms 3072
