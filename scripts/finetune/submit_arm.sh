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
# Pre-flight: featurize every record and drop the ones that raise. Idempotent - it keeps
# manifest_*_full.json and re-filters from that each time. Costs ~2 minutes of the GPU
# allocation and removes the failure that killed the first launch three minutes in
# (KeyError 'C8' on one ligand out of 388, with four hours of GPU queued behind it).
srun ./env/bin/python validate_records.py --arm "$ARM" --workers 8

srun ./env/bin/python train_arm.py \
    --arm "$ARM" --steps "$STEPS" \
    --lr 3e-4 --warmup 50 --accum 4 \
    --workers 4 --max-tokens 384 --max-atoms 3072
