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
    ap.add_argument("--dry-run", action="store_true",
                    help="build everything and report shapes without stepping the optimiser")
    a = ap.parse_args()

    ck = ROOT / "boltz_cache" / "boltz2_conf.ckpt"
    manifest = ROOT / "records" / a.arm / "manifest.json"
    recs = json.loads(manifest.read_text())
    train = [r for r in recs if r["split"] == "train"]
    test = [r for r in recs if r["split"] == "test"]
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

    # override only the low-N knobs, leave loss weights and architecture untouched
    ta = dict(hp.get("training_args", {}))
    ta["max_lr"] = a.lr
    ta["lr_warmup_no_steps"] = a.warmup
    hp["training_args"] = ta

    import inspect

    from boltz.model.models.boltz2 import Boltz2
    from boltz.model.modules.diffusionv2 import AtomDiffusion

    # The RELEASED checkpoint was trained with a different boltz version than the one on
    # PyPI: its diffusion_process_args carry `mse_rotational_alignment` and
    # `step_scale_random`, which 2.2.1's AtomDiffusion.__init__ does not accept, and which
    # appear NOWHERE in 2.2.1's source. Dropping them changes nothing about how this
    # version behaves - the code never reads them - but passing them is a hard TypeError.
    #
    # Filter against the real signature rather than hardcoding a drop-list, so a future
    # version that adds parameters back keeps them automatically.
    def fit_kwargs(d: dict, fn) -> tuple[dict, list[str]]:
        allowed = set(inspect.signature(fn).parameters)
        keep = {k: v for k, v in d.items() if k in allowed}
        dropped = sorted(set(d) - set(keep))
        return keep, dropped

    dpa, dropped = fit_kwargs(dict(hp.get("diffusion_process_args", {})),
                              AtomDiffusion.__init__)
    if dropped:
        print(f"dropped {len(dropped)} diffusion args unknown to this boltz: {dropped}",
              flush=True)
    hp["diffusion_process_args"] = dpa

    kw, dropped_top = fit_kwargs(hp, Boltz2.__init__)
    if dropped_top:
        print(f"dropped {len(dropped_top)} top-level args: {dropped_top}", flush=True)
    model = Boltz2(**kw)

    # strict=False is required across the version gap, but a silent partial load would
    # train a randomly-initialised model and look fine. Report the counts.
    res = model.load_state_dict(sd["state_dict"], strict=False)
    missing, unexpected = list(res.missing_keys), list(res.unexpected_keys)
    n_par = sum(p.numel() for p in model.parameters())
    n_train = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"model built: {n_par/1e6:.1f}M params, {n_train/1e6:.1f}M trainable", flush=True)
    print(f"weight load: missing={len(missing)} unexpected={len(unexpected)}", flush=True)
    if missing[:3]:
        print("  first missing:", missing[:3], flush=True)
    if unexpected[:3]:
        print("  first unexpected:", unexpected[:3], flush=True)
    frac = 1 - len(missing) / max(len(list(model.state_dict())), 1)
    print(f"  fraction of model tensors filled from checkpoint: {frac:.4f}", flush=True)
    if frac < 0.95:
        print("REFUSING: more than 5% of the model did not load. Training this would be "
              "training a partly-random model that reports a plausible loss.", flush=True)
        return 2

    if a.dry_run:
        print(json.dumps({
            "arm": a.arm, "train": len(train), "test": len(test),
            "params_M": round(n_par / 1e6, 1),
            "lr": a.lr, "warmup": a.warmup, "steps": a.steps,
            "status": "DRY RUN - model constructs and weights load; no optimiser step taken",
        }, indent=2))
        return 0

    print("full training path not yet wired - run with --dry-run", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
