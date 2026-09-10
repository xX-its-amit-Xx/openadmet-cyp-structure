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
