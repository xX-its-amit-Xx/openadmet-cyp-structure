# RUNBOOK — Boltz-2 co-folding on Explorer (SLURM)

**Status: WORKING, exercised end to end on 2026-09-22.** This is the generation path to
use on drop day. Modal is over its spend cap and OpenProtein cannot serve `boltz2` with an
uploaded MSA (`FINDING_034`), so this is the only venue left that can produce poses at the
**same conditioning as the validation pool** — which `FINDING_034` measured to be the
difference between depth that pays and depth that costs 0.0065 of score.

Four GPUs on `gyorilab` with a 30-day walltime, no GPU cap, and a general `gpu` partition
behind it as overflow. One 14-ligand × 20-sample job is about one GPU-hour.

---

## What already exists — inventory before you build anything

The fine-tuning campaign (`scripts/finetune/`) left a complete **offline** Boltz-2 install
on `/scratch/shenoy.am/cyp-finetune`. Do not rebuild it.

| asset | path | verified |
|---|---|---|
| boltz **2.2.1** + torch 2.5.1 venv | `cyp-finetune/env/bin/boltz` | job `10528945` |
| Boltz-2 weights | `cyp-finetune/boltz_cache/boltz2_conf.ckpt` (2.3 GB) | md5 `2f0a1775bf8fc366a1a85e2019eca288` |
| CCD | `cyp-finetune/boltz_cache/ccd.pkl` | md5 `aa43b8ba95608ce179033e990e39b4f3` |
| component mols | `cyp-finetune/boltz_cache/mols/`, `mols.tar` | present |
| **185 P450 MSAs**, a3m | `cyp-finetune/msa/` | 185 files, 1.0 GB |
| the same 185 as Boltz `.npz` | `cyp-finetune/msa_npz/` | 185 files, 304 MB — **training only**; `boltz predict` reads a3m straight from the YAML |
| working sbatch templates | `cyp-finetune/submit_holdout.sh` | 5 COMPLETED inference jobs, ~43 min each |

**boltz 2.2.1 is exactly what Modal pinned** (`boltz[cuda]==2.2.1`), and
`boltz_cache/get.sh` pulls the checkpoint from
`huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt` — the same artefact
`boltz predict` downloads for itself. **The conditioning match is free.**

The CYP3A4 alignment is *not* in the 185 — those are the P450 superfamily targets. The one
the validation pool was folded against lives in the repo at `data/reference/cyp3a4.a3m`
(6,979 sequences, md5 `6de0ee1330ea9aa1aef1efb6ec6a2e3d`) and is staged by step 2.

---

## The procedure

Everything below is `scripts/explorer/boltz_depth.py`. Each stage is resumable and prints
counts rather than a status word.

### 0. Connectivity

```bash
ssh explorer 'hostname; squeue -u shenoy.am'
```

### 1. Build the inputs, locally

```bash
python scripts/explorer/boltz_depth.py plan            # the 14-ligand stratum
python scripts/explorer/boltz_depth.py plan --smoke    # one non-stratum ligand
python scripts/explorer/boltz_depth.py conditioning    # what must match, with hashes

# a blind test set: every id in the csv, under its own tag
python scripts/explorer/boltz_depth.py plan \
    --csv data/processed/test_ligands.csv --tag drop
```

`plan` does **not** re-implement the YAML. It lifts `build_yaml` out of
`scripts/cofold/modal_boltz.py` by AST extraction and executes it, so the input is written
by the same function that wrote the pool being extended.

**On the blind-set flag, and the correction that matters.** This runbook used to say
"point it at the new ligand CSV; nothing else changes" — **there was no flag to do that**
until 2026-09-23. `plan` read `validation_ligands.csv` and `binding_mode_labels_cyp3a4.csv`
and filtered to the 14-ligand predicted-Type-I stratum, so on a test set it would have
written the wrong inputs for the wrong ligands. `--csv` now takes **every** id in the file
(the Type-I filter is a property of the FINDING 035 pre-registration, not of the
generator), and `--tag` names the subdirectory that `stage`, `submit`, `collect` and
`score` all key off — on both sides. Without `--tag`, every stage behaves exactly as
before (`stratum`, or `smoke` with `--smoke`). `plan` records what it did in
`C:/cyp_struct/matched_depth/plans/<tag>.json`, which is where `collect` gets the expected
ligand list and `score` learns whether the set is blind. Exercised on the 87 validation
ligands presented as `id,smiles` only: 87/87 YAMLs, 0 skipped.

**`score` now refuses a blind tag** rather than returning zero rows. It needs a `pdb` per
ligand; without one every ligand landed in `skipped["no_crystal"]` and it wrote an empty
CSV while exiting 0 (playbook G2).

### 2. Stage to `/scratch` — from a **login** node

```bash
python scripts/explorer/boltz_depth.py stage --smoke
python scripts/explorer/boltz_depth.py stage
python scripts/explorer/boltz_depth.py push            # the sbatch templates
```

**This is the step that matters.** GPU nodes have no direct internet, so every input must
be on disk before the job starts. Login nodes *do* have internet — fetch there, compute
there never.

`stage` re-hashes the alignment on the far side and reports
`msa_md5_matches_local`. Do not proceed on a transfer that only claims to have worked.

### 3. Verify the big files — **inside a job**

```bash
python scripts/explorer/boltz_depth.py verify
```

Submits `verify_inputs.sbatch` (CPU, `--mem=16G`, `short` partition, ~20 s) and prints its
log. It md5s the 2.3 GB checkpoint, counts MSA sequences and NUL bytes, prints the boltz
version, reports whether the cuequivariance kernels import, and runs Boltz's own YAML
parser over every input with structural asserts (chain order, `SG`↔`FE` bond, MSA path
exists).

**Never md5 the checkpoint on the login node.** Its memory cgroup SIGKILLs long reads and
the killed process can still exit 0 — a hash of nothing, reported as success.

### 4. Submit

```bash
python scripts/explorer/boltz_depth.py submit --seed 101 --samples 20
```

`sbatch`, never `srun` over ssh — three earlier attempts in this project were killed by
client-side timeouts mid-run. The job id is appended to
`data/processed/matched_depth_jobs.json`.

### 5. Poll for ADVANCING progress

```bash
python scripts/explorer/boltz_depth.py poll
python scripts/explorer/boltz_depth.py tail --job <id>
```

`poll` prints the queue **and the count of mmCIFs on disk**. Liveness is not progress:
judge by the file count going up.

### 6. Collect and score

```bash
python scripts/explorer/boltz_depth.py collect
python scripts/explorer/boltz_depth.py score
```

`collect` pulls only `*_model_*.cif` and flattens them to `<ligand>__s<seed>_m<k>.cif`;
the per-pose confidence `.npz` blobs stay on `/scratch` (D: has ~13 GB and is not a place
to land a pool). `score` is `FINDING_034`'s scorer with a different input directory —
LDDT-PLI, BiSyRMSD, `xeng`, Fe geometry, renumbering offset and identity, md5.

---

## The sbatch templates

`scripts/explorer/boltz_depth.sbatch`, in full, with the three lines that are load-bearing:

```bash
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=07:00:00
...
export http_proxy=http://10.99.0.130:3128
export https_proxy=http://10.99.0.130:3128

srun "$FT/env/bin/boltz" predict "$ROOT/yaml/$SUB" \
    --checkpoint "$FT/boltz_cache/boltz2_conf.ckpt" \
    --cache "$FT/boltz_cache" \
    --out_dir "$OUT" --output_format mmcif \
    --diffusion_samples "$SAMPLES" --recycling_steps 3 --sampling_steps 200 \
    --seed "$SEED" --num_workers 4 --no_kernels --override
```

**`--partition` is passed on the command line, not baked in.** `gyorilab` is the priority
partition but its four H200s are routinely held by multi-day interactive `bash` sessions —
they were 4/4 allocated when this was written, to jobs 12 and 6 days old. Submitting to
`gyorilab,gpu` lets Slurm take whichever frees first, and the smoke test duly landed on
the general `gpu` partition (a V100 on `d1010`) in under a minute instead of waiting.

**`--time=07:00:00` is a consequence of that fallback.** The general `gpu` partition caps
at 8 hours; a 24-hour request would be silently unschedulable there and would pin the job
to the one busy partition. Inference is ~45 minutes, so 7 hours is ample.

**`--no_kernels` is required, not an optimisation.** `boltz predict` sets
`use_kernels = not no_kernels` and so defaults to the cuequivariance path; this env's
`cuequivariance-*-cu12` wheels were built against a torch that was later rolled back and
the backend import dies on `libnvrtc.so.12: cannot open shared object file`. Two
fine-tuning holdout jobs died in under 70 seconds on exactly that. It is also the one
genuine mismatch against the Modal pool, which ran the kernel path — see `FINDING_035`
confound K1.

`scripts/explorer/verify_inputs.sbatch` asks for **no GPU**: it is I/O and an import
check, and requesting a GPU would put the de-risking step behind the queue it exists to
de-risk.

---

## Every failure hit, and its fix

| # | symptom | cause | fix |
|---|---|---|---|
| 1 | `rsync: dup() in/out/err failed`, 0 bytes, on **every** invocation — through Python and typed straight into the shell alike | MSYS2 rsync ↔ OpenSSH stdio incompatibility on this box. Nothing to do with the cluster | **rsync is unusable from here.** Transport is `scp` for single files and `tar czf - \| ssh … tar xzf -` for directories. Both work |
| 2 | `echo "rsync rc=$?"` printed **`rc=0`** for a transfer that moved nothing | `$?` had captured `tail`'s status through the pipe, not rsync's | `${PIPESTATUS[0]}`, and **count the files on the far side** after every push. A failed copy reporting success is the trap in `docs/README.md` §2 |
| 3 | staging returned `ok: false` with an **empty stderr** and a null return code — three times | `subprocess` ran bare `bash`, which on this box resolves to `C:\Windows\System32\bash.exe`, the **WSL launcher**. It answers every command with "Windows Subsystem for Linux has no installed distributions" **on stdout** and exits 1 | pin the absolute Git Bash path (`_git_bash()` in `boltz_depth.py`) and refuse to fall back |
| 4 | (inherited, `submit_holdout.sh`) job dies in <70 s | cuequivariance backend import | `--no_kernels` |
| 5 | (inherited, PXR, ~12 days) jobs hang or die on download | GPU nodes have no direct internet; `--use_msa_server` cannot work | stage the MSA on a login node; `--cache` at a local dir; set the proxy so an accidental reach-out fails fast instead of hanging |
| 6 | `#!/bin/bash^M: bad interpreter` (anticipated) | CRLF from the Windows checkout | `push` runs `sed -i 's/\r$//'` on the far side and prints `head -1 \| cat -A` to prove it |
| 7 | `MSA file cyp3a4.a3m not found` although it is next to the YAML | Boltz resolves a relative `msa:` against the **process** cwd, not the YAML's directory | the YAML carries the **absolute** `/scratch/...` path |
| 8 | `collect` would have swept the smoke run's poses into the stratum set | the pull walked all of `out/`, and the flattening step keys on the ligand directory name, so it could not have noticed | scope the remote `find` to `./<sub>_s*`. Caught before it mattered; it is the same shape as `FINDING_034`'s "two input sets writing to one directory produce one result wearing two labels" |

---

## Do NOT

- **Do not `srun` over ssh** for anything that takes minutes. Killed three times already.
- **Do not read or hash a large file on the login node.** Fetch on login, verify in batch.
- **Do not trust a status.** `verify` exercises the parser; `poll` counts files; `stage`
  re-hashes on the far side. Every one of these replaced a check that could pass while
  broken.
- **Do not count jobs as poses.** `FINDING_034` found 196 of 196 complexes byte-identical
  across fourteen replicate jobs. Dedupe on geometry
  (`matched_depth_analysis.py distinct`, tol 0.05 Å) before quoting any depth.
  *Native* `boltz predict` does diversify the ligand with `--diffusion_samples` — that is
  an OpenProtein-wrapper limitation, not a Boltz one — but count it anyway, every time.
- **Do not bake `--partition=gyorilab` in.** It is frequently 4/4 allocated to long
  interactive sessions.
- **Do not land a pool on `D:`.** ~13 GB free. Bulk stays on `/scratch`; local
  intermediates go to `C:\cyp_struct` via `cypstruct.paths.SCRATCH`.

---

## Drop-day sequence

1. `ssh explorer 'sinfo -p gyorilab,gpu -o "%P %t %G"'` — pick the partition list.
2. `plan --csv data/processed/test_ligands.csv --tag drop`. The MSA is already staged and
   the target sequence does not change. Pass the same `--tag` to every stage below.
3. `stage --tag drop` → `verify` → `submit --tag drop`. Budget ~1 GPU-hour per 14 ligands
   × 20 samples; scale linearly and split into several jobs if the set is large, since
   each is independently resumable.
4. `poll` on the file count. `collect --tag drop` — it now also reports `missing_ligands`
   and `poses_per_ligand` against the plan, which is check 4 of the four signals.
   **Do NOT run `score`**: it needs a crystal per ligand and refuses a blind tag.
5. Label the test ligands prediction-side (`FINDING_033` item 1) before selecting.
6. Select with `cypstruct.xengine.select()` — unchanged (`FINDING_033` item 2).
