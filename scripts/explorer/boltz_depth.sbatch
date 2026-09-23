#!/bin/bash
#SBATCH --job-name=cypdep
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=07:00:00
#SBATCH --output=/scratch/shenoy.am/cyp-depth/logs/cypdep_%j.log
#
# Boltz-2 co-folding at conditioning MATCHED to the Modal pool.
#
#   sbatch --partition=gyorilab,gpu boltz_depth.sbatch <sub> <seed> <samples>
#
# --partition is passed on the command line, not baked in, because `gyorilab` is the
# priority partition but its four H200s are routinely held by multi-day interactive
# sessions. Listing `gyorilab,gpu` lets Slurm take whichever frees first. --time is
# 07:00:00 rather than days BECAUSE of that fallback: the general `gpu` partition caps at
# 8 hours and silently refuses anything longer, so a 24-hour request would pin the job to
# the one busy partition. Inference here takes ~45 minutes.
#
# Submitted with sbatch, never srun'd over ssh: three earlier attempts in this project
# were killed by client-side timeouts mid-run. A batch job survives the laptop.
set -euo pipefail
SUB=${1:-stratum}
SEED=${2:-101}
SAMPLES=${3:-20}

ROOT=/scratch/shenoy.am/cyp-depth
FT=/scratch/shenoy.am/cyp-finetune          # env + weights + CCD, from the finetune campaign
cd "$ROOT"

OUT="$ROOT/out/${SUB}_s${SEED}"
mkdir -p "$OUT"

# GPU nodes have no direct internet. Set the proxy so that IF anything reaches out it
# fails fast through a route that exists, rather than hanging on a dead socket for the
# whole walltime -- which is how the PXR campaign lost ~12 days. Nothing below should
# need it: weights, CCD and the alignment are all on /scratch.
export http_proxy=http://10.99.0.130:3128
export https_proxy=http://10.99.0.130:3128
export HTTP_PROXY=$http_proxy
export HTTPS_PROXY=$https_proxy
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

echo "=== node $(hostname) | $(date) | sub=$SUB seed=$SEED samples=$SAMPLES ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "=== inputs ==="
ls -la "$ROOT/msa/cyp3a4.a3m"
echo "yaml files: $(ls "$ROOT/yaml/$SUB"/*.yaml | wc -l)"
md5sum "$ROOT/msa/cyp3a4.a3m"

# --no_kernels is REQUIRED, not an optimisation. `boltz predict` sets
# use_kernels = not no_kernels and so defaults to the cuequivariance path; this env's
# cuequivariance-*-cu12 wheels were built against a torch that was later rolled back and
# the backend import dies on `libnvrtc.so.12: cannot open shared object file`. Two
# fine-tuning holdout jobs died in under 70 seconds on exactly that. The pure-PyTorch
# fallback is mathematically the same model, and it is the one confound this path carries
# against the Modal pool, which ran the kernel path.
#
# --cache points at the finetune campaign's local boltz_cache so no code path attempts a
# download. --override makes a resubmission idempotent.
srun "$FT/env/bin/boltz" predict "$ROOT/yaml/$SUB" \
    --checkpoint "$FT/boltz_cache/boltz2_conf.ckpt" \
    --cache "$FT/boltz_cache" \
    --out_dir "$OUT" \
    --output_format mmcif \
    --diffusion_samples "$SAMPLES" \
    --recycling_steps 3 \
    --sampling_steps 200 \
    --seed "$SEED" \
    --num_workers 4 \
    --no_kernels \
    --override

echo "=== produced ==="
find "$OUT" -name '*_model_*.cif' | wc -l
find "$OUT" -name '*_model_*.cif' | head -3
echo "=== done $(date) ==="
