"""Fine-tune Boltz-2 on one arm, with the gates that make the result mean something.

Runs ON EXPLORER. Wires `Boltz2` (which ships `training_step`, `validation_step` and
`configure_optimizers`) to a Lightning `Trainer`. No fork, no reimplementation.

**The recipe is set relative to the model's OWN pretraining**, read out of the checkpoint's
`hyper_parameters.training_args` rather than borrowed from a different paper:

| | pretrain | here | why |
|---|---|---|---|
| `max_lr` | 1e-3 | **3e-4** | 0.3x. Low-N adaptation should move weights, not rewrite them. |
| warmup steps | 1000 | **50** | 1000 warmup steps on a ~350-step run would never leave warmup. |
| `diffusion_loss_weight` | 4.0 | **4.0, unchanged** | It already dominates confidence (0.3) and distogram (0.03) by >10x, and it governs COORDINATES - which is exactly our failure mode (ligand orientation). Reweighting would tune away from a recipe that demonstrably produced a working model. |

**What this is trying to beat, and why it is not another helping of a failed idea.** Every
diversity experiment in this repo added oracle the selector could not reach: a union of
engines gave +0.0375 oracle / -0.0017 selected, a sampler sweep +0.0355 / -0.0038. Those
add poses to the TAIL. A fine-tune moves the whole distribution - it raises the MEAN, and a
better mean needs no selector at all. Different mechanism.

**Two gates, both pre-registered.**

1. *Held-out gain.* Scored on scaffold-cluster holdouts whose chemical novelty was matched
   to the real challenge split (median NN 0.587, FINDING 017). arm4_mix is the arm to
   trust at 0.536; arm3_cyp3a4_only sits at 0.758 and would validate against near-analogues
   of its own training set, so a win there means little.
2. *Forgetting.* A model that wins on CYP by forgetting everything else is worse than
   useless on a blind set whose composition we do not control. IntFold's LoRA precedent
   kept 35 of 35 unrelated complexes while gaining 4 of 5 held-out conformations; that is
   the shape to check for, not just the gain.

    ./env/bin/python train_arm.py --arm arm4_mix --steps 350
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path("/scratch/shenoy.am/cyp-finetune")


def main() -> int:
    import torch

    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="arm4_mix")
    ap.add_argument("--steps", type=int, default=350)
    ap.add_argument("--lr", type=float, default=3e-4)      # 0.3x the 1e-3 pretrain LR
    ap.add_argument("--warmup", type=int, default=50)      # vs 1000 in pretraining
    ap.add_argument("--accum", type=int, default=4)
    ap.add_argument("--clip", type=float, default=10.0)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--max-tokens", type=int, default=512)
    ap.add_argument("--max-atoms", type=int, default=4096)
    ap.add_argument("--diffusion-multiplicity", type=int, default=1,
                    help="pretrain used 32; see the note in main() for why 1 is forced")
    ap.add_argument("--dry-run", action="store_true",
                    help="build everything and report shapes without stepping the optimiser")
    a = ap.parse_args()

    from boltz.data.types import Manifest

    ck = ROOT / "boltz_cache" / "boltz2_conf.ckpt"
    arm_dir = ROOT / "records" / a.arm
    train = Manifest.load(arm_dir / "manifest_train.json").records
    test = Manifest.load(arm_dir / "manifest_test.json").records
    print(f"arm={a.arm}  train={len(train)}  test={len(test)}", flush=True)

    # Load the checkpoint's own hyper_parameters so the model is constructed exactly as it
    # was trained - guessing these is how a fine-tune silently trains a different model.
    sd = torch.load(ck, map_location="cpu", weights_only=False)
    hp = dict(sd["hyper_parameters"])
    print("pretrain max_lr:", hp.get("training_args", {}).get("max_lr"), flush=True)
    print("pretrain warmup:", hp.get("training_args", {}).get("lr_warmup_no_steps"),
          flush=True)
    print("diffusion_loss_weight:", hp.get("training_args", {}).get("diffusion_loss_weight"),
          flush=True)

    # Override only the low-N knobs, leave loss weights and architecture untouched.
    # This has to be passed to load_from_checkpoint as a kwarg: mutating the dict read
    # out of the file changes nothing, because Lightning re-reads hyper_parameters from
    # the checkpoint itself. An earlier version edited `hp` and trained at 1e-3.
    ta = hp.get("training_args")
    try:
        ta.max_lr = a.lr
        ta.lr_warmup_no_steps = a.warmup
    except (AttributeError, TypeError):
        ta = dict(ta)
        ta["max_lr"] = a.lr
        ta["lr_warmup_no_steps"] = a.warmup
    # diffusion_multiplicity 1, NOT the stored 32. forward() reshapes feats["coords"] to
    # (B*multiplicity, L, 3) for the diffusion loss, and the confidence branch then asserts
    # coords.shape[0] == 1 - so anything above 1 makes every batch raise "Validation is not
    # supported for batch sizes=N" and return None. Since confidence cannot be switched off
    # either (see below), 1 is the only self-consistent setting in the released code.
    # The cost is gradient variance, which accumulate_grad_batches buys back; the benefit
    # is that the memory freed lets the crop stay large enough to hold the whole pocket.
    if a.diffusion_multiplicity:
        try:
            ta.diffusion_multiplicity = a.diffusion_multiplicity
        except (AttributeError, TypeError):
            ta["diffusion_multiplicity"] = a.diffusion_multiplicity

    # Construct through boltz's OWN loading path, exactly as main.py:1314 does.
    #
    # My first attempt filtered the checkpoint's hyper_parameters into the constructor by
    # hand. That produced 128 missing pairformer tensors and a 506.8M-vs-521.0M parameter
    # gap which looked like a version mismatch, and was not: boltz overrides the stale
    # stored args with FRESH dataclass defaults from the installed version, rather than
    # replaying what the checkpoint recorded. Feeding the old values back is what broke it.
    #
    # strict=True is the point. It is a real assertion that every tensor loads - stronger
    # and simpler than the load-fraction heuristic I was using, which passed at 0.9755 on
    # a model that was quietly wrong.
    from dataclasses import asdict

    from boltz.main import (
        Boltz2DiffusionParams,
        BoltzSteeringParams,
        MSAModuleArgs,
        PairformerArgsV2,
    )
    from boltz.model.models.boltz2 import Boltz2

    diffusion_params = Boltz2DiffusionParams()
    # Activation checkpointing on both trunks. Inference leaves these False because it
    # runs under no_grad and stores nothing; training a 64-block pairformer does, and the
    # triangular-multiplication activations alone OOM a 140 GB H200 at 384 tokens. This
    # trades recompute for memory and changes no numbers.
    pairformer_args = PairformerArgsV2(activation_checkpointing=True)
    msa_args = MSAModuleArgs(subsample_msa=True, num_subsampled_msa=1024,
                             use_paired_feature=True, activation_checkpointing=True)
    steering_args = BoltzSteeringParams()
    steering_args.fk_steering = False
    steering_args.physical_guidance_update = False

    model = Boltz2.load_from_checkpoint(
        str(ck),
        strict=True,
        map_location="cpu",
        diffusion_process_args=asdict(diffusion_params),
        ema=False,
        use_kernels=False,          # pure-PyTorch path; the kernel needs cuequivariance,
                                    # and installing that upgraded torch and broke boltz
        pairformer_args=asdict(pairformer_args),
        msa_args=asdict(msa_args),
        steering_args=asdict(steering_args),
        training_args=ta,
        # THE CHECKPOINT STORES structure_prediction_training=False. boltz2_conf.ckpt is
        # the *confidence* training stage, and Boltz2.__init__ acts on that flag by
        # setting requires_grad=False on everything outside confidence_module /
        # affinity_module / out_token_feat_update. Inherit it and the run trains only the
        # confidence head - the one signal this repo has measured as useless for ranking
        # poses (rho -0.092, FINDING 011) - while training_step also skips the distogram
        # and diffusion losses entirely, so diffusion_loss_weight 4.0 never applies.
        # A 4-step run did exactly this and exited 0: param_norm_structure_module 0.0,
        # grad_norm identical to grad_norm_confidence_module.
        structure_prediction_training=True,
        checkpoint_diffusion_conditioning=True,
    )
    # `validate_structure` is a constructor argument of Boltz2 that is never assigned to
    # self, so `Boltz2.setup()` raises AttributeError on the first non-predict stage -
    # i.e. the released model cannot enter a training loop as shipped. (Same shape as
    # trainingv2.py being Boltz-1's: the training path was not exercised on release.)
    # Setting it False both fixes the attribute and is what we want: boltz's internal
    # validators report their own metrics, and our gate is LDDT-PLI from cypstruct.pose.
    model.validate_structure = False

    # Confidence stays ON, and that is forced by an upstream coupling rather than by
    # preference. forward() builds `diffusion_conditioning` under the gate
    # `(not self.training) or self.confidence_prediction`, but the structure-training
    # branch further down uses that variable unconditionally - so training the trunk with
    # confidence off raises UnboundLocalError on the first batch. Confidence contributes
    # at weight 0.3 against diffusion's 4.0, so it is a rounding error on the objective
    # we care about, and the head is not what we gate on (FINDING 011: rho -0.092).

    # The checkpoint's diffusion_loss_args carry `add_bond_loss`, which the installed
    # AtomDiffusion.compute_loss no longer accepts. training_step splats these as **kwargs
    # inside a try/except that prints and returns None, so EVERY batch was skipped while
    # the run exited 0 with a full progress bar. Version drift between a released
    # checkpoint's stored args and the installed package, hidden by a broad except.
    #
    # Dropping unknown keys is only safe when they are inactive, so this refuses rather
    # than guesses: add_bond_loss is False in the released checkpoint, and the remaining
    # five keys match the installed signature exactly.
    import inspect

    allowed = set(inspect.signature(model.structure_module.compute_loss).parameters)
    stored = dict(model.diffusion_loss_args)
    dropped = {k: v for k, v in stored.items() if k not in allowed}
    if any(bool(v) for v in dropped.values()):
        msg = ("diffusion_loss_args has ACTIVE keys the installed compute_loss does not "
               f"accept: {dropped}. Dropping them would change the objective.")
        raise SystemExit(msg)
    if dropped:
        print("dropped inactive diffusion_loss_args:", dropped, flush=True)
    model.diffusion_loss_args = {k: v for k, v in stored.items() if k in allowed}

    n_par = sum(p_.numel() for p_ in model.parameters())
    trainable = {}
    for name, p_ in model.named_parameters():
        if p_.requires_grad:
            trainable[name.split(".")[0]] = trainable.get(name.split(".")[0], 0) + p_.numel()
    print(f"model loaded strict=True: {n_par/1e6:.1f}M params", flush=True)
    print("effective max_lr:", model.training_args.get("max_lr"), flush=True)
    print("trainable by module (M):",
          {k: round(v / 1e6, 1) for k, v in sorted(trainable.items(), key=lambda kv: -kv[1])},
          flush=True)

    # A run that trains only the confidence head exits 0 and looks fine. Refuse to start
    # one by accident: the structure trunk must be receiving gradient.
    if not any(k in trainable for k in ("structure_module", "pairformer_module")):
        msg = ("structure trunk is frozen - this would fine-tune the confidence head "
               "only. Check structure_prediction_training.")
        raise SystemExit(msg)

    sys.path.insert(0, str(ROOT))
    from boltz2_data import Boltz2FinetuneDataModule

    dm = Boltz2FinetuneDataModule(
        arm_dir=arm_dir,
        msa_dir=ROOT / "msa_npz",
        mol_dir=ROOT / "boltz_cache" / "mols",
        batch_size=1,
        num_workers=a.workers,
        max_tokens=a.max_tokens,
        max_atoms=a.max_atoms,
    )

    if a.dry_run:
        dm.setup()
        batch = next(iter(dm.train_dataloader()))
        shapes = {k: tuple(v.shape) for k, v in list(batch.items())
                  if hasattr(v, "shape")}
        print(json.dumps({
            "arm": a.arm, "train": len(train), "test": len(test),
            "params_M": round(n_par / 1e6, 1),
            "lr": a.lr, "warmup": a.warmup, "steps": a.steps,
            "batch_keys": len(shapes),
            "coords": shapes.get("coords"),
            "token_index": shapes.get("token_index"),
            "msa": shapes.get("msa"),
            "status": "DRY RUN - weights load and one real batch featurizes; no optimiser step",
        }, indent=2))
        return 0

    import pytorch_lightning as pl
    from pytorch_lightning.callbacks import ModelCheckpoint
    from pytorch_lightning.loggers import CSVLogger

    out = ROOT / "runs" / a.arm
    out.mkdir(parents=True, exist_ok=True)

    # Lightning validation is OFF on purpose. Boltz2.validation_step dispatches through
    # `self.validator_mapper[batch["idx_dataset"]]`, boltz's own validation harness, which
    # reports its internal metrics - not LDDT-PLI against our crystal set. The gate this
    # project pre-registered is held-out LDDT-PLI measured by cypstruct.pose, computed by
    # running inference from the saved checkpoint. Wiring a validator that reports a
    # different number than the gate is how a run looks healthy and fails the gate.
    trainer = pl.Trainer(
        accelerator="gpu", devices=1, precision="bf16-mixed",
        max_steps=a.steps, accumulate_grad_batches=a.accum,
        limit_val_batches=0, num_sanity_val_steps=0,
        gradient_clip_val=a.clip,
        log_every_n_steps=1,
        logger=CSVLogger(str(out), name="csv"),
        callbacks=[ModelCheckpoint(dirpath=str(out), save_last=True,
                                   every_n_train_steps=max(a.steps // 4, 1),
                                   save_top_k=-1, filename="step{step}")],
        enable_progress_bar=True,
    )
    trainer.fit(model, datamodule=dm)
    print(json.dumps({
        "arm": a.arm, "steps": a.steps, "out": str(out),
        "status": "TRAINED",
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
