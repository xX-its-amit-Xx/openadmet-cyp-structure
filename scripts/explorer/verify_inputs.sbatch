#!/bin/bash
#SBATCH --job-name=cypver
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/shenoy.am/cyp-depth/logs/verify_%j.log
#
# Verify the staged inputs INSIDE A JOB, with an explicit --mem.
#
# This exists because of a specific, measured trap: Explorer's LOGIN node runs under a
# memory cgroup that SIGKILLs long or memory-hungry reads, and the killed process can
# still leave a pipeline reporting SUCCESS. An md5sum of a 2.3 GB checkpoint on the login
# node is exactly that check -- it can print a hash of nothing and exit 0. Fetch on the
# login node (it has internet); verify in batch (it has memory).
#
# No GPU is requested: this is I/O and an import check, and asking for a GPU would put it
# behind the queue it is meant to de-risk.
set -uo pipefail
ROOT=/scratch/shenoy.am/cyp-depth
FT=/scratch/shenoy.am/cyp-finetune

echo "=== node $(hostname) | $(date) ==="
free -g | head -2

echo "=== MSA: must be the SAME alignment the Modal pool used ==="
ls -la "$ROOT/msa/cyp3a4.a3m"
md5sum "$ROOT/msa/cyp3a4.a3m"
echo "sequences: $(grep -c '^>' "$ROOT/msa/cyp3a4.a3m")"
echo "NUL bytes: $(tr -dc '\000' < "$ROOT/msa/cyp3a4.a3m" | wc -c)"

echo "=== checkpoint (2.3 GB -- the read that must not happen on a login node) ==="
ls -la "$FT/boltz_cache/boltz2_conf.ckpt"
md5sum "$FT/boltz_cache/boltz2_conf.ckpt"
echo "=== ccd ==="
ls -la "$FT/boltz_cache/ccd.pkl"
md5sum "$FT/boltz_cache/ccd.pkl"

echo "=== boltz version (must be 2.2.1, matching Modal's pin) ==="
"$FT/env/bin/python" - <<'PY'
import importlib.metadata as md
import torch
print("boltz  ", md.version("boltz"))
print("torch  ", torch.__version__)
try:
    import cuequivariance_ops_torch  # noqa: F401
    print("cueq   ", "IMPORTS -- kernels available")
except Exception as e:
    print("cueq   ", f"BROKEN ({type(e).__name__}) -- --no_kernels is required")
PY

echo "=== yaml sanity: does Boltz's own schema parser accept them? ==="
"$FT/env/bin/python" - <<'PY'
import glob, os, sys, yaml
# EVERY staged subdirectory, not a hardcoded pair. The list used to read
# ("smoke", "stratum"), which silently verified nothing at all for any later tag --
# a check that passes because it looked at no files is the repo's standing trap.
ROOT = "/scratch/shenoy.am/cyp-depth/yaml"
subs = sorted(d for d in os.listdir(ROOT) if os.path.isdir(os.path.join(ROOT, d)))
print("subdirectories found:", subs)
for sub in subs:
    fs = sorted(glob.glob(f"/scratch/shenoy.am/cyp-depth/yaml/{sub}/*.yaml"))
    if not fs:
        print(f"{sub}: 0 files -- NOTHING VERIFIED")
        continue
    ok = bad = 0
    for f in fs:
        try:
            d = yaml.safe_load(open(f))
            assert d["version"] == 1
            ids = [list(s)[0] for s in d["sequences"]]
            assert ids == ["protein", "ligand", "ligand"], ids
            b = d["constraints"][0]["bond"]
            assert b["atom1"][2] == "SG" and b["atom2"][2] == "FE", b
            import os
            assert os.path.exists(d["sequences"][0]["protein"]["msa"]), "msa path missing"
            ok += 1
        except Exception as e:
            bad += 1
            print("  BAD", f, type(e).__name__, e)
    print(f"{sub}: {ok} parse+assert OK, {bad} bad, {len(fs)} files")
PY
echo "=== done $(date) ==="
