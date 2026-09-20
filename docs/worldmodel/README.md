# The shared-embedding world model — plan, status, and kill criteria

**Started 2026-09-20.** Three asks, one object. Kept here so the chain of reasoning is
followable and so no agent's research is lost to a context window.

---

## The unifying idea

A JEPA (Joint Embedding Predictive Architecture) predicts *representations* of masked
parts rather than reconstructing inputs. Pointed at interaction data it becomes:

> given the embedding of entity A and a context, predict the embedding of the entity B
> that binds it.

That is a **generative prior over binders in embedding space**. It is the same object
ask 3 wants — "reason backwards from what CYP3A4 evolved to do, to what it should bind" —
approached from the data side instead of the causal side. Where the learned prior and the
causal derivation agree, we have a candidate. Where they disagree is the interesting part.

It is also a strict generalisation of proteochemometrics: PCM concatenates a drug
embedding with a target embedding and regresses affinity. A JEPA instead learns the
*map* between the two spaces, so it can be asked "what should bind here?" and not only
"do these two bind?".

---

## Honest constraints, recorded up front

| constraint | value | consequence |
|---|---|---|
| GPU | **1× H200**, contended (3 of 4 held by another user for 3–9 days) | no from-scratch pretraining; frozen foundation embeddings + light predictor head |
| Local disk | C: ~34 GB, D: ~20 GB | nothing bulk lands locally, ever |
| Cluster scratch | `/scratch` 298 TB | all corpora live here |
| Cold storage | OneDrive ~4.9 TB via `cypstruct.storage` | finished corpora, never the hot path |
| Interim deadline | **2026-09-24** (4 days) | this program contributes nothing to it; interim ships from the corrected pool |
| Final deadline | 2026-11-03 | realistic target for a proof of concept, not a trained world model |

**The part that cannot be crammed is the data.** Compute can be bought later; a corpus
assembled badly has to be reassembled. So reconnaissance and acquisition run first and
continuously, and modelling waits on them.

---

## The three asks

### Ask 1 — the world model
Two arms, deliberately opposed:
- **Implicit.** Give it everything — DNA, RNA, protein, ligand, every affinity readout —
  and see what structure it finds unprompted.
- **Explicit.** Train/fine-tune only on binding interfaces structurally similar to
  CYP3A4's, and see whether a narrow, well-chosen prior beats a broad one.

The comparison *is* the experiment. If broad beats narrow, biological priors are better
learned than imposed; if narrow wins, the opposite. Either answer is worth having.

### Ask 2 — causal decomposition of CYP3A4
Folds across apo/holo and isoforms; the phylogenetic tree and what each branch was
selected *for*; ancestral and current binding partners; PTMs as evidence of use-cases;
liver cell-type and transcriptional context; then reasoning forward from shape and
sequence to what *should* bind, including chemistry we have no assay for.

### Ask 3 — the dashboard
Real, interactive, continuously updated: the phylogenetic tree, the domains, the
structures, the embedding space as it forms. Merged into the Vercel dashboard.

---

## Kill criteria (pre-registered, per this repo's rules)

1. **Embedding transfer.** If frozen foundation embeddings cannot beat a Morgan-fingerprint
   + one-hot-target PCM baseline on a held-out affinity task, the shared space is not
   carrying biology and the program stops at the baseline.
2. **Partner prediction.** If the JEPA head cannot retrieve true partners above chance on
   held-out interactions (AUC > 0.6 on a difficulty-matched negative set), the predictive
   framing is wrong.
3. **Difficulty matching.** Per FINDING 022, any benchmark whose baseline is already
   >0.80 LDDT-PLI or >85% sub-2 Å is *uninformative* and must not be used to claim a gain.
4. **Cost.** If acquisition exceeds ~10 TB on `/scratch` before a single predictive result,
   stop and prune.

---

## Status

| wave | scope | state |
|---|---|---|
| 1 | model + data reconnaissance, 6 parallel agents | **launched 2026-09-20** |
| 2 | acquisition to `/scratch`, deduplication, splits | pending wave 1 |
| 3 | embedding extraction, JEPA head, baselines | pending wave 2 |
| 2b | CYP3A4 causal decomposition | **launched with wave 1** |

Findings land in `docs/worldmodel/*.md`, one file per domain, written by the agent that
did the work. This file is the index and the argument; those are the evidence.

---

## Revision 1 — 2026-09-20, after the protein/ligand reconnaissance

**The naive shared space is already falsified in the literature.** Separate unimodal
encoders, concatenated or attention-fused — the classical PCM architecture — loses to a
count fingerprint:

| evidence | result |
|---|---|
| 25 models × 25 datasets, scaffold split (arXiv 2508.06199) | **ECFP mean rank 7.52**; MoLFormer 9.50, Mol2Vec 10.36, Uni-Mol2 13.32, GraphMVP 14.64 |
| 62,820-model study (Nat Commun) | fingerprints competitive or better |
| 132-dataset peptide study | fingerprints beat ESM-2 and ProtT5 |
| Volkov 2022 | interaction descriptors add nothing over ligand-only |
| Graber 2025 (Nat MI) | ~49% of CASF-2016 leaks from PDBbind; **ligand-only matches structure-based DL** |
| TDC ADMET audit 2026 | only fingerprint + boosting entries reproduced |

So the program is redirected before it is built, not after.

### What survives, and why

`boltz predict --write_embeddings` is a shipped flag (`boltz/main.py:1038`). It writes
`embeddings_<id>.npz` containing `s [n, 384]` and `z [n, n, 128]`.

**`z` is a function of protein × ligand × pose.** That is categorically different from
gluing a protein vector to a ligand vector: the fingerprint critique applies to models
whose joint representation is assembled *after* the fact, and `z` is joint by
construction, computed by a network trained to place atoms in space. It is MIT, it runs
on the torch 2.5.1 stack we already validated, and it falls out of inference runs we
are already paying for.

`cuequivariance` is an *optional* extra of boltz, not a requirement — consistent with
what we found the hard way during fine-tuning.

### The pre-registered first experiment

> Can the Boltz-2 pair representation predict CYP3A4 affinity better than a fingerprint?

- **Signal:** pooled `z` from a Boltz-2 pose, plus `s`.
- **Control that must be beaten:** count-ECFP4 + LightGBM.
- **Two further controls, both mandatory:** a **ligand-only** model (catches the shortcut
  where target features carry everything) and a **nearest-neighbour-by-Tanimoto**
  predictor (catches memorisation).
- **Split:** target-cluster-held-out, never random. Report the ligand-free baseline
  alongside every number.
- **Kill:** if pooled `z` does not beat count-ECFP4 on a target-held-out split, the
  shared-space premise is dead and we say so.

This is days of work, not months, and it tests the load-bearing assumption first.

### Operational notes worth not rediscovering

- **MoLFormer-XL must be called with `deterministic_eval=True`** or its embeddings differ
  between calls — a silent corruption of exactly the kind this project keeps hitting.
- `z` is ~128 MB per pose. **Pool at write time**, never store raw.
- `HF_HOME` must point into `/scratch`; the model shortlist alone is ~30 GB and local
  disk has 20 GB.
- Dead ends confirmed: GearNet (TorchDrug pins Python <3.11), ESM-IF1 (torch-scatter),
  SMI-TED (2021 CUDA kernels), UMA/OMol25 (fairchem will not install on torch 2.5.1).
- Non-commercial licences to avoid unless the user says otherwise: ESM C 600M, ESM-3,
  Ankh/Ankh3, MACE-OFF, Chai-1 local weights. MIT substitutes exist for all of them.

---

## Revision 2 — 2026-09-20, after the multimodal/JEPA reconnaissance

Three independent reconnaissance lines now agree, which is worth more than any one of them.

### The convergent finding: joint-by-construction, or nothing

| source | finding |
|---|---|
| protein/ligand recon | learned molecular embeddings lose to count-ECFP4 across 25×25 scaffold-split benchmarks |
| multimodal recon | two leak-controlled ablations find protein LM embeddings contribute **~nothing**; a 21-token residue vocabulary *beats* ESM-2 on strict tiers |
| multimodal recon | AI-Bind: shuffling **all** SMILES *and* sequences moved AUROC 0.86 → 0.84 — the model was not learning the pairing at all |
| multimodal recon | DUD-E ligand-only vs receptor-ligand features correlate at **R² = 0.98** |

A shared space assembled by concatenating frozen unimodal encoders is the architecture
all of that evidence indicts. Two representations escape the critique because they are
**joint by construction** — computed by a network that had to place atoms in space:

1. **Boltz-2's pair representation** `z [n,n,128]`, a function of protein × ligand × pose,
   free from runs we already pay for.
2. **ATOMICA's interface graph** — atoms within 8 Å of the partner, both sides, one
   SE(3)-equivariant graph with explicit *intermolecular* edges.

Everything else in the program is downstream of testing these two.

### The split trap, quantified

| model | random split | honest split |
|---|---|---|
| Papyrus PCM | r = 0.79 | r = 0.42 (temporal) |
| DrugBAN | 0.960 | 0.575 (cluster) |

**And the method ranking flips.** A random-split comparison does not merely inflate every
number, it reorders which approach looks best. Any result in this program computed on a
random split is not a weak result, it is a meaningless one.

### What does not exist yet

Bio-JEPA is real and recent — Mol-JEPA (Boehringer; 14 modalities including Boltz-2
embeddings; benchmarks on PXR, our sibling target), ProtJEPA, JEPA-DNA, Cell-JEPA and
others, nearly all inside 12 months. **But every one masks *within* a single entity.
Nobody predicts a binding partner's embedding from a target's.** The object described at
the top of this file does not exist in the literature.

That is either a real gap or a silent graveyard. Treat it as unproven either way, and
note Mol-JEPA's own Wilcoxon test favours a *tabular baseline* 47% to 38%.

### ATOMICA: usable today, with a mandatory ablation

ATOMICA embeds an **already-bound interface**, so it cannot answer "what should bind
here" — there is no unbound encoder to query. It *can* score poses, and CYP3A4 is a heme
protein while ATOMICA-Ligand ships a heme checkpoint. Its headline result is wet-lab:
5 of 6 predicted proteins confirmed heme binding.

The catch, and the reason the ablation is not optional: **86% of its pretraining corpus is
CSD small-molecule crystal packing**, and it was trained by denoising rotations and
translations. It may be a crystallinity detector. So it is scored against native poses
*and* deliberately rotated ones; if it separates rotated decoys easily but cannot rank
real prediction errors, that is the finding.

**Launched 2026-09-20** against the existing 3,360-pose CYP3A4 pool. No new data.

### Scope decision

With one contended H200 between now and 2026-11-03 we can afford a leak-controlled
baseline **or** a JEPA prototype, not both. **The baseline goes first**, because without
it no JEPA number is interpretable — and because three separate lines of evidence say the
baseline may simply win.

### Licence default (override if wrong)

Proceeding **MIT/Apache-only** until told otherwise. Excluded on that basis: ESM C 600M,
ESM-3, Ankh/Ankh3, MACE-OFF, Chai-1 local weights, ProtJEPA, Top-DTI, DTIAM (all NC);
TransformerCPI2.0 (GPL-2.0); HyperAttentionDTI (no licence at all). MIT substitutes exist
for every one of these.

---

## Revision 3 — 2026-09-20, after the DNA reconnaissance

### Volume is not the constraint. Supervision density is.

The whole useful DNA corpus is **under 200 GB**, and the part that actually teaches
*interaction* is far smaller than that:

| modality | non-redundant interaction complexes |
|---|---|
| protein–ligand (HiQBind) | **>30,000** over >18,000 PDB entries |
| protein–DNA (PDB, 30% identity clusters) | **3,027** — from 10,733 raw entries |
| protein–protein with ΔΔG (SKEMPI) | 345 unique structures, everything else is mutations of those |

So the implicit "give it everything" arm has a problem the breadth conceals: **DNA and RNA
contribute enormous *sequence* corpora and almost no *supervised interaction* signal.**
8.8 Tbp of genome against ~10⁻⁹ supervision density. Adding them may improve
representation pretraining; it will not add pairs to learn binding from.

That does not kill the implicit arm — it reframes what it can possibly be testing. It is
a representation-transfer experiment, not a data-scaling one, and it should be reported
as such.

### The prior art we should read before writing any code

**JEPA-DNA (NVIDIA, arXiv 2602.17162)** ran a JEPA objective across five DNA backbones,
with **Apache-2.0 code and checkpoints**. That is the closest existing thing to what this
program proposes, in the open, on a permissive licence. Studying it is cheaper than
rediscovering it and it is the natural template for the masked-within-entity half of the
design. (It still does not do partner prediction — see Revision 2 — so the gap stands.)

### Model shortlist, MIT/Apache only

| model | licence | size | dim | note |
|---|---|---|---|---|
| **Evo 2 7B** (`arcinstitute/evo2_7b`) | Apache-2.0, ungated | 13.8 GB | 4096 | 1 Mb context; bf16 without Transformer Engine. The 40B is 82.2 GB and will not fit |
| **AlphaGenome-PyTorch** (`gtca/alphagenome_pytorch`) | open reimpl. | 450M | 3072 @128 bp | 40.8 GB peak on one H200, plain torch |
| **ModernGENA-large** | open | 377M | — | pure `transformers`; the one that works first try on torch 2.5.1 |

### A licence trap worth naming

**AlphaGenome's own weights are non-commercial *and restrict training on its outputs*.**
Not merely "don't sell the model" — a restriction that propagates into anything we train
on predictions it produced. Excluded under the MIT/Apache default, and worth remembering
as a class of term that does not show up in a headline licence label. The open PyTorch
reimplementation is the route if we want that architecture.

### Blocked, needs the user

AlphaGenome gate, NTv3 gate, AlphaFold3 weights application. BioLiP2 / DeepPBS / GraphBind
returned 403 or bad TLS from this network. ProNAB bulk may need an email.
