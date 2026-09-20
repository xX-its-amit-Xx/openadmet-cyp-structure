# Fine-tuning plan — and an honest case for doing it second

**Status:** infrastructure built, training not yet launched. 2026-09-10.

---

## The case against rushing into it

The PXR structure campaign attempted an OpenFold3 fine-tune and it cost enormous runway
for nothing measurable:

- Five platforms over roughly twelve days produced **fifty training steps**, all on Kaggle.
- The Explorer attempt never started training at all — three successive debug jobs died on
  weight download, ultimately because **GPU compute nodes have no direct internet**.
- The result that did get measured: **−0.0020** in a controlled 9-vs-8-model A/B, and
  0.5395 standalone against the campaign's 0.5640 best.
- The campaign's own recorded verdict: *"Not killed. Deferred because training data,
  compute, and validation needed more runway."*

So fine-tuning was never disproven — it was never really run. But with **15 days to the
interim leaderboard**, an unbounded-runway activity that has historically consumed
everything and returned nothing is the second bet, not the first.

**What goes first:** the steered pose pool and the physics scorer, both of which are cheap,
both of which are already measurable against 99 ligands of crystal ground truth.

---

## Why it is worth doing anyway, and why it looks better this time

Three things are genuinely different from the PXR attempt:

1. **We have far more data for one target.** RCSB holds **122 CYP3A4 entries**, of which
   99 give a usable (ligand, crystal pose) pair. The published low-N precedent fine-tuned
   OpenFold3 on **ten** complexes and improved interface lDDT and DockQ measurably. We have
   an order of magnitude more, for a single target, with the cofactor present in all of them.
2. **There is now a documented low-N recipe**, rather than pretraining defaults. Apheris's
   PDE10A case study publishes the deltas that matter, and they are all in the direction of
   "let a small dataset actually move the weights without destroying them":

   | hyperparameter | pretraining default | low-N value |
   |---|---|---|
   | learning rate | 1.8e-3 | **3e-4** |
   | warmup steps | 1000 | **50** |
   | EMA decay | 0.999 | **0.99** |
   | gradient steps | many | **~350** |
   | crop / token budget | large | **384** |

   Roughly 20 hours on one 80 GB GPU. The PXR campaign's three OF3 configs specified no LR
   at all; its Latch curriculum did, and that curriculum is the other reusable idea:
   **3e-4 → 1e-4 → 5e-5 → 2e-5** with `interface_weight` ramping **1.0 → 1.5 → 2.0 → 3.0**
   and crop **384 → 384 → 256 → 256**.

3. **LoRA is a real option now, with a published precedent for exactly our situation.**
   IntFold applied per-layer LoRA with a frozen base to a single target (CDK2): the base
   model captured **0 of 5** allosteric conformations, the fine-tune captured **4 of 5
   while keeping 35 of 35** on the common state. That is a measured no-catastrophic-
   forgetting result on a single-target adaptation. `wiwnopgm/boltz-finetune` (MIT) exposes
   `use_lora`, `lora_r: 8`, `lora_alpha: 16`, `freeze_all` with per-module overrides.

---

## What we would actually train, and on what

**Target of adaptation.** Not "CYP3A4 in general" — the specific failure mode. Co-folding
already places the heme and the backbone correctly; it gets **ligand orientation** wrong.
So the objective weight goes on the interface, and the held-out metric is BiSyRMSD and
LDDT-PLI on ligands the model never saw — not global lDDT, which is already fine and will
happily improve while the thing we care about does not.

**Dataset.** `scripts/finetune/build_dataset.py` emits, per usable CYP3A4 entry:
mmCIF with protein + HEM + ligand, the ligand SMILES, the precomputed target MSA (one
sequence, so one MSA shared across the whole set), and the measured binding class.

**Splitting — the part most likely to produce a fake win.** The 99 ligands are *not* 99
independent examples. They cluster hard by chemical series: 62 of the 83 coordinated
ligands present a pyridine donor, and the ritonavir-analog series alone contributes
several near-identical entries (1RD, 5AW, 6AW, 7AW …). A random split would put analogs on
both sides and report a large, meaningless gain.

So: **cluster by Murcko scaffold and Tanimoto, hold out whole clusters**, and report the
number of independent held-out clusters alongside the metric. If that number is small
enough that nothing clears significance, say so — that is the honest result.

**Forgetting check.** Every checkpoint is also evaluated on a set of unrelated
protein-ligand complexes. A model that wins on CYP3A4 by forgetting everything else is
worse than useless for a blind set whose composition we do not control.

---

## Order of operations

| step | what | cost | gate to proceed |
|---|---|---|---|
| 0 | **Steered Boltz-2, no training.** Heme bonded to Cys442, template from nearest holo, pocket restraint, forced Fe→donor contact at 2.4 Å. Steered vs unsteered arms. | ~5 GPU-h | does steering beat unsteered on the 99-ligand ground truth? |
| 1 | Add engines (AF3, Chai, Protenix, OpenFold3) as **tail rescuers**, per the PXR rule | ~15 GPU-h | does each rescue depth sweep show an interior optimum? |
| 2 | Physics scorer (tier 0 → tier 1 QM → tier 2 xTB) | ~10 CPU-h | does it beat the z-hybrid on held-out clusters? |
| 3 | **LoRA fine-tune of Boltz-2** on the 99 complexes, low-N recipe | ~20 GPU-h | does it beat step 0 on held-out clusters, without forgetting? |
| 4 | Full OpenFold3 low-N fine-tune, Apheris recipe | ~20 GPU-h | only if 3 shows signal |

Steps 0–2 are the submission. Steps 3–4 are the upside, and they are only worth their
runway if the earlier gates pass.

---

## Known risks

- **Protenix may be a poor base for this target.** Its default training set is named
  `weightedPDB_before2109_wopb_nometalc_0925` — the `nometalc` token suggests metal-cluster
  complexes were filtered out. For a heme enzyme that would matter. *Unverified*; the token
  is undocumented. Check before spending on Protenix fine-tuning.
- **Boltz-2 training code is not released.** `docs/training.md` upstream still reads
  "Coming soon updated training information for Boltz-2!", and the fine-tune issues are
  unanswered. The LoRA route depends on third-party forks built against Boltz-1.
- **AF3 weights are licence-restricted** (non-commercial, no redistribution). They can be
  fetched to Modal from Kaggle under the user's own credentials for inference; they cannot
  be redistributed, and fine-tuning them is outside the terms.
- **"Apo drift" is documented across Boltz-2, Chai-1 and Protenix** — models predicting the
  ligand-free state despite a ligand present. That is precisely the CYP3A4 pocket-expansion
  problem, so it is a plausible thing for fine-tuning to fix, and equally a plausible thing
  for it to fail at.

---

## CORRECTION, 2026-09-19 — Boltz-2 training code IS released

The "Known risks" section above says:

> **Boltz-2 training code is not released.** `docs/training.md` upstream still reads
> "Coming soon updated training information for Boltz-2!", and the fine-tune issues are
> unanswered. The LoRA route depends on third-party forks built against Boltz-1.

**That is no longer true.** `pip install boltz` (2.2.1, on PyPI) ships the full training
stack inside the package:

| component | location |
|---|---|
| `Boltz2(LightningModule)` | `boltz/model/models/boltz2.py:40` |
| `training_step` | `boltz2.py:793` |
| `validation_step` | `boltz2.py:1002` |
| `configure_optimizers` | `boltz2.py:1132` |
| `BoltzTrainingDataModule` | `boltz/data/module/trainingv2.py:467` |
| `TrainingDataset` / `ValidationDataset` | `trainingv2.py:168 / 317` |
| losses, optimisers | `boltz/model/loss/`, `boltz/model/optim/` |

The only thing absent is a `train` subcommand — `boltz/main.py` registers `predict` and
nothing else. So fine-tuning needs a script that wires `BoltzTrainingDataModule` to
`Boltz2` and calls `pl.Trainer.fit`, which is ordinary Lightning work, **not** a
third-party fork and **not** a fallback to OpenFold3.

### Why this is a materially better position than the plan assumed

We fine-tune **the same architecture that generated our 16,500-pose corpus**. A gain
transfers directly to the pool we already have and to the cross-engine reference built on
it, instead of requiring a new engine whose poses would have to be re-validated from
scratch. The plan's step 4 ("full OpenFold3 low-N fine-tune, only if step 3 shows signal")
can stay closed unless Boltz fails.

The LoRA precedent still applies and is still the first configuration to try: IntFold's
per-layer LoRA on a frozen base captured 4 of 5 held-out CDK2 conformations while keeping
35 of 35 on the common state — a measured no-catastrophic-forgetting result on exactly this
shape of single-family adaptation.

### Environment, as actually built

On `/scratch/shenoy.am/cyp-finetune`, all offline:

```
env/          micromamba-built Python 3.11.16, torch 2.5.1+cu121
msa/          185 a3m, 1009 MB          <- the PXR blocker, closed
rcsb/         506 crystal mmCIF, 516 MB
finetune_arms/  4 split CSVs
```

GPU verified through Slurm on `gyorilab`: `cuda avail True`, `NVIDIA H200 NVL`,
`bf16 True`. The PXR campaign never reached a running GPU job; this one has.

### The checkpoint is a full training checkpoint, not an inference export

Verified on a GPU node, 2026-09-19:

```
top-level keys: epoch, global_step, pytorch-lightning_version, state_dict,
                loops, callbacks, optimizer_states, lr_schedulers
tensors: 5102          parameters: 521.0M
hyper_parameters: atom_s, atom_z, token_s, token_z, num_bins, training_args,
                  validation_args, embedder_args, msa_args, pairformer_args,
                  score_model_args, diffusion_process_args
```

`optimizer_states`, `lr_schedulers` and `loops` survive in the released file, so training
state was not stripped. `hyper_parameters.training_args` means the original training
configuration is **recoverable rather than guessed** - which matters because the PXR
attempt's three OpenFold3 configs specified no learning rate at all and nobody could say
afterwards what had actually been run.

Combined with `Boltz2.training_step` and `configure_optimizers` shipping in the package,
every piece needed to resume or adapt training is present. Nothing here requires a fork.

### The actual pretraining recipe, recovered from the checkpoint

Not guessed from another paper's defaults — read out of `hyper_parameters.training_args`:

| parameter | pretrain value |
|---|---|
| `max_lr` | **1e-3** |
| `base_lr` | 0.0 |
| `lr_scheduler` | `af3` |
| `lr_warmup_no_steps` | 1000 |
| `lr_start_decay_after_n_steps` | 50000 |
| `weight_decay` | 0.003 (excluded on norm/bias) |
| adam β1 / β2 / eps | 0.9 / 0.95 / 1e-8 |
| `diffusion_loss_weight` | **4.0** |
| `confidence_loss_weight` | 0.3 |
| `distogram_loss_weight` | 0.03 |
| `affinity_loss_weight` | 0.003 |
| `recycling_steps` | 3 |
| `diffusion_multiplicity` | 32 |
| stopped at | epoch 37, global_step 23,750 |

Two things follow.

**The low-N LR is now relative to reality.** Apheris's 3e-4 is 0.3× this model's true
pretrain LR, and warmup drops 1000 → 50. The plan's table assumed a 1.8e-3 pretrain
default; the real value is 1e-3, so the intended reduction is milder than it looked.

**The loss we need is already the dominant term.** `diffusion_loss_weight` 4.0 governs
coordinates and outweighs confidence (0.3) and distogram (0.03) by more than an order of
magnitude. Our failure mode is ligand *orientation* — a coordinate error — so the gradient
signal is pointed at the right thing without reweighting. The plan's instinct to "put the
objective weight on the interface" is already satisfied by the stock configuration, and
changing it would be tuning against a recipe that produced a working model.

### Data path

Boltz ships its own mmCIF parser (`data/parse/mmcif.py`, `mmcif_with_constraints.py`) and
`Manifest` / `Record` as JSON-serialisable types. So the 506 crystal structures convert
through **Boltz's own parser**, not a reimplementation — which matters because a mismatch
between how training data is parsed and how inference parses it would degrade the
fine-tune silently rather than failing.

### Incident: a pip install silently replaced torch, and the dependency was never needed

Boltz inference failed with `ModuleNotFoundError: cuequivariance_torch`, raised from
`kernel_triangular_mult`. I installed the package. That was wrong twice over.

**It broke the environment.** `pip install cuequivariance-torch` pulled a CUDA 13 stack and
**upgraded torch 2.5.1+cu121 → 2.14.0+cu130**, after which boltz would not import at all.
The stack I had verified against an H200 was replaced by a resolver decision I never asked
for and did not check.

**And it was unnecessary.** `boltz/model/layers/triangular_mult.py` takes
`use_kernels: bool = False` on `forward`, with a pure-PyTorch path behind it. The kernel is
an optimisation. The correct fix was to leave kernels off, not to satisfy the import.

Recovery: `pip install --force-reinstall torch==2.5.1 --index-url .../cu121`, then
re-verified on a GPU node rather than trusting the import — `cuda avail True`, `NVIDIA
H200 NVL`, `use_kernels default: False`.

**Rule going forward:** this environment is pinned in `VERIFIED_VERSIONS.txt`:

```
boltz==2.2.1   torch==2.5.1+cu121   numpy==1.26.4   pytorch-lightning==2.5.0
gemmi==0.6.5   pandas==3.0.6        scipy==1.13.1   rdkit==2026.3.6
```

Any future install here passes those constraints, because a CUDA-adjacent package can
rewrite the framework underneath a working stack and report success while doing it. And a
missing *optional* accelerator is not a blocker — check for a fallback flag before
installing anything to satisfy an ImportError.

### Inference verified: the released checkpoint runs on the PyPI package

```
boltz predict t.yaml --cache ./boltz_cache --use_msa_server --no_kernels
  Number of failed examples: 0
  t_out/boltz_results_t/predictions/t/t_model_0.cif
  7 s on an H200
```

**`--no_kernels` is the whole fix.** `main.py:1321` sets `use_kernels=not no_kernels`, so
the CLI turns kernels ON by default while the module default is `False`. The
`ModuleNotFoundError: cuequivariance_torch` was one flag, never a missing dependency — and
chasing it is what upgraded torch and broke the environment.

**This also retires the version-gap worry.** Three ticks were spent on 128 missing
`pairformer_module.layers.N.attention.norm_s.*` tensors and a 506.8M-vs-521.0M parameter
gap, which looked like a genuine architecture mismatch. It was not: those were artifacts of
constructing `Boltz2(**filtered_kwargs)` by hand. Lightning's own loader reports only
`v2.5.0.post0` vs `v2.5.0`. The supported path works; the hand-rolled one was the problem.

Lesson for the trainer: **build the model through boltz's own loading path**, not by
filtering hyper_parameters into the constructor. The checkpoint and the package agree with
each other; they disagreed with me.

---

## Addendum 2026-09-20 — the released Boltz-2 cannot train as shipped

Wiring the data path turned up six defects, five of them silent. Recording them because
each one produces a run that *completes* and is worthless, which is the expensive kind.

**1. `trainingv2.py` is Boltz-1's trainer wearing a v2 name.** It imports `BoltzFeaturizer`
and the v1 `Structure`. `Boltz2Featurizer` is referenced by exactly one file in the entire
package — `inferencev2.py`. Handing v1 features to `Boltz2.training_step` does not raise.
The training path in `scripts/finetune/boltz2_data.py` is therefore built from the
*inference* v2 path — the one the released weights were exported against — with
`training=True` and a cropper added.

**2. `self.validate_structure` is never assigned.** `Boltz2.__init__` takes
`validate_structure` as an argument and uses it in four methods, but never stores it, so
`setup()` raises `AttributeError` the moment `stage != "predict"`. The released model
cannot enter a Lightning training loop without a patch. Same shape as (1): the training
path was not exercised on release.

**3. Symmetry correction raises on every batch, and the handler hides it.**
`training_step` routes ground truth through `minimum_lddt_symmetry_coords`, which reads
three ragged symmetry keys out of the batch. Without `compute_symmetries=True` every batch
throws, `training_step` returns `None`, and the run finishes with a full progress bar and
no gradient. Worse, the handler prints `batch['pdb_id']` — a key nothing in the released
pipeline puts there — so the real error is replaced by a `KeyError`.

**4. `collate` exempts its ragged keys by name.** The v2 featurizer emits at least one
ragged key that is not on the six-name list, so the worker dies with `'list' object has no
attribute 'shape'`. Replaced with a type check, which is the same rule and cannot go stale.

**5. Two incompatible RNG APIs in one item.** `BoltzCropper` calls `random.randint(n)` —
legacy `RandomState`. `Boltz2Featurizer` expects a `Generator`. One seed, two objects.

**6. Our own: `MOLDIR` pointed at `mols/mols`,** which survived the earlier flattening as
an *empty directory*, so `.exists()` stayed `True` and the parser got a mol dir with
nothing in it. 28 structures lost. An existence check on a directory is not a check that
it has contents; the check is now `any(MOLDIR.glob("*.pkl"))`.

### Data path, as built

    a3m (185)  --make_msa_npz-->  msa_npz/<target_key>.npz     185/185, 0 failures
    mmCIF      --make_records-->  records/<arm>/structures/<pdb>.npz + manifest*.json

**Splits are per-PDB, not per-pair.** `build_arms.py` splits ligand pairs, but one npz
holds every ligand in the entry — a PDB with one train pair and one test pair would put
the held-out ligand's coordinates into training. Any PDB touching test goes to test whole.
On arm4_mix this demoted 0 pairs, but the guard stays.

arm4_mix: 406 pairs → 404 PDBs → **388 parsed, 303 train / 85 test, 0 records without an
MSA**. The 16 failures are CCD components (`A1Axx`) newer than boltz's cached dictionary.

First real batch: 78 feature tensors, `coords [1,1,4096,3]`, `msa [1,1024,512]` — the MSA
depth confirms the a3m conversion actually reaches the model rather than defaulting to
single-sequence.

### Validation is deliberately OFF in Lightning

`Boltz2.validation_step` dispatches through `self.validator_mapper[batch["idx_dataset"]]`,
boltz's own harness, which reports its internal metrics. The gate this project
pre-registered is held-out **LDDT-PLI measured by `cypstruct.pose`**. A validator that
reports a different number than the gate is how a run looks healthy and fails the gate.
Held-out scoring runs as a separate inference pass from the saved checkpoint.

### Capacity, measured 2026-09-20 — the campaign is serial

`gyorilab` is one node, `d4079`, with four H200s. Three of them are held by another user's
interactive `bash` sessions, running 3 and 9 days respectively. That leaves **one GPU**.

A held-out control submitted alongside training went `PD / Resources` with an estimated
start three weeks out — past both the interim and final deadlines. Slurm's estimate is
pessimistic, but the constraint is real and is not something to fix by killing another
user's work.

So the campaign runs **serially, with dependencies**, not in parallel:

    sbatch submit_arm.sh arm4_mix 350                          # ~3.7 h
    sbatch --dependency=afterany:<train> submit_holdout.sh arm4_mix base
    sbatch --dependency=afterany:<base>  submit_holdout.sh arm4_mix <ckpt>

`afterany`, not `afterok`: the base control does not depend on training succeeding, and a
failed train should not silently cancel the control that would have explained it.

Budget at one GPU: ~3.7 h per arm plus two ~2 h inference passes ≈ **8 h per arm**, so all
four arms are roughly a day and a half of wall clock. That is affordable, but it means
arm4_mix has to be the one that runs first — it is the novelty-matched arm (NN 0.536
against the challenge's 0.587), so it is the only one whose held-out number transfers.

Second venue: `ssh discovery` fails host-key verification from this box, so it is not
available without a key the user would need to accept. Not pursued.
