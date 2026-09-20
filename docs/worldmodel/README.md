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

---

## Revision 4 — 2026-09-20, after the RNA reconnaissance

### RNA–small-molecule binding data does not exist, and that is now measured

| quantity | count |
|---|---|
| aggregated RNA–small-molecule measurements worldwide (R-SIM) | ~2,500 |
| usable Kd | 1,480 |
| with a structure | **101** |
| non-redundant drug-like co-structures ever solved | **48** |

And nearly all of it is the same PDB entries and the same R-BIND/Inforna extract re-cut
under new names. So RNA joins DNA: enormous sequence corpora, essentially **zero**
ligand-binding supervision. Forty-eight complexes is not a modality, it is a footnote.

The "train on every modality's affinity data" framing is therefore settled: for *ligand
binding*, protein–ligand is not one modality among several, it is ~99% of the evidence
that exists anywhere.

### The architecture the recon actually found

Here is the thing worth keeping. **RhoFold+ exposes `c_s = 384` and `c_z = 128` — exactly
Boltz-2's trunk dimensions.** That is not luck; both follow the AlphaFold-style
single/pair convention. Which means:

> The shared space should be built by aligning **pair representations from structure
> predictors**, not by aligning sequence-model embeddings.

Every one of those pair representations is joint by construction — protein × ligand ×
pose for Boltz-2, RNA × RNA × geometry for RhoFold+ — which is precisely the property
that survives the fingerprint critique in Revision 2. And they already live in a common
dimensionality, so one projection head serves both rather than a bespoke adapter per
modality.

This is a concrete, cheap, falsifiable version of what the program set out to do, and it
came out of the evidence rather than being assumed at the start.

Supporting facts: **RiNALMo-giga is 1280-d, matching ESM-2-650M**, so the sequence side is
symmetric too if we want it. **ATOMICA ships finetuned HEM/HEC heme checkpoints** —
independently confirmed by two reconnaissance lines, and CYP3A4 is a heme protein.

### Prior art to read before designing anything

**RNAPro (NVIDIA × Das Lab, 2026)** is already a frozen RNA foundation model gated into a
Protenix co-folder — i.e. someone has built the "frozen foundation model feeds a structure
predictor" design and published what happened. Alongside JEPA-DNA (Revision 3), that is
two of the three legs of this program already standing in the open literature.

### Shortlist additions, MIT/Apache-compatible

| model | licence | size | dims | note |
|---|---|---|---|---|
| RiNALMo-giga | CC-BY-4.0 | 2.60 GB | 1280 | original repo pins `flash-attn==2.3.2`; use the MultiMolecule rewrap (plain PyTorch) |
| RhoFold+ | Apache-2.0 | 0.51 GB | s 384 / z 128 | same convention as Boltz-2 |
| RNAcentral r27 | CC0 | 10 GB | — | 58,558,809 sequences |
| RNA3DB (2026-01-05) | MIT | 2.15 GB | — | ships defensible non-redundant splits |
| Ribonanza | open | — | — | chemical mapping on ~2 M sequences; the one RNA layer richer than its protein analogue |

### A legal call that needs making

**MultiMolecule is AGPL-3.0.** That is not the usual non-commercial annoyance — AGPL is
viral across a network boundary, so serving a model through it can oblige us to release
our own source. Using it to *produce embeddings offline* is a different matter from
shipping it in anything. Flagged rather than decided.

### Blocked, needs the user

R-BIND 2.0 (URL dead; data only in paywalled ACS supporting information), Inforna 2.0
(signed licence), PDBbind registration, OpenFold3 HF gate, RIBOSPAN-10K (on request),
AlphaFold3 (non-redistributable). Given the 48-structure ceiling above, **none of the
RNA-specific ones are worth chasing.**


---

## Revision 5 — 2026-09-20, the CYP3A4 causal decomposition

Full argument in `CYP3A4_EVOLUTION.md` (~11,300 words, DOIs/PDB IDs/UniProt throughout,
inferences explicitly tagged). What it changes here:

### The convergence worth noting

Its structural conclusion — **"the iron anchor is saturated and uninformative;
discrimination lives under the F/G roof"** — is the same claim this repo reached
empirically and independently, twice: FINDING 008 (anchor-local features all fail;
recalibrating the coordination window moved selection by -0.0002) and FINDING 019
(unbonded predictions already reproduce crystal Fe geometry, so an explicit heme bond is
a null at p = 0.43).

One route is pose measurement on 85 held-out complexes; the other is comparative genomics
and mutagenesis literature. They agree. That is the strongest support any claim in this
project currently has, and it says: **stop working on the anchor.**

### Selection tuned regioselectivity, not breadth

The three sites under selection across the whole phylogeny are **codons 437, 478 and
479**, on the cavity floor — and mutagenesis shows they change **where the substrate gets
oxidised, not how tightly it binds**. Two lineage-specific bursts dominate (CYP3A7 at the
hominoid stem, omega = infinity on 15.5 nonsynonymous / 0 synonymous, simultaneous with
its restriction to fetal liver; human CYP3A4 after the chimp split, 6.0/0.0) against a
background where 89% of sites are purifying.

The standard "CYP3A expanded to detoxify plant secondary metabolites" story is therefore
**not supported by the selection signal** — nothing in it demonstrably tuned promiscuity.
CYP3A4 is also under hard physiological constraint: the single point mutation **I301T
causes a Mendelian disease** (vitamin D-dependent rickets type 3, doi:10.1111/febs.70277)
via novel 11-alpha-hydroxylated vitamin D metabolites. The "drug enzyme" framing is a
modern accident of which molecules we happen to hand it.

### The pocket, as an inference

Built to hold a **rigid, membrane-partitioned, mostly apolar 300-600 Da molecule loosely
enough to oxidise it at several positions**, polar end tethered at the rim
(Ser119 / Arg212 / Thr224 / backbone amides), entrance opening *sideways into the bilayer*
through a one-way F-F' gate.

Chemical space it should accept but that has never been assayed: sulfated lipid
conjugates, endocannabinoid and oxylipin chemistry, beyond-rule-of-5 macrocycles and
PROTACs, and **two-ligand pairs where a co-binder converts an uncoupled non-substrate into
a productive one** - a drug-drug interaction class with no standard assay.

### Two facts with direct operational consequences

1. **Ligands *increase* CYP3A4's global dynamics** (HDX-MS,
   doi:10.1016/j.jinorgbio.2023.112211) rather than rigidifying it on binding. Any scoring
   function assuming induced-fit rigidification has the sign wrong.
2. **All 122 CYP3A4 PDB entries are X-ray. There is no cryo-EM CYP3A4 structure.** That is
   exactly the gap the OpenADMET release fills, so the challenge references are not merely
   new data, they are the first of their kind for this target.

### Multi-occupancy: checked, and it is rare

See `FINDING_W001_multicopy_check.md`. 5 of 87 validation entries hold more than one copy
of the query ligand (max 3). Three of the four real cases rank in the bottom ten of 87 by
pool oracle - and the fourth, caffeine with three copies, ranks fourth. n = 4. Logged as a
hypothesis with a defined test, not as a result.

---

## Revision 6 — 2026-09-20, revised RNA reconnaissance

Two items that change the build, both verified in source rather than inferred:

**1. Boltz-2 already does RNA, and already writes embeddings.** `entity_type: rna` plus
`--write_embeddings` are both shipped (`main.py:1038`, `writer.py:250`,
`schema.py:1021`). So the cross-modal pair-representation test in Revision 4 needs **one
model we already run**, not a second predictor and an adapter. RhoFold+ becomes a
cross-check rather than a dependency, and the whole idea gets testable this week.

**2. A free blind test set, if we act now.** **CASP17 releases coordinates for 52 RNA
targets on 2026-11-29.** Freezing a PDB training cutoff *today* buys an uncontaminated
blind evaluation for nothing. Costs one line in a manifest now; cannot be recovered later.

The honest scale number also got smaller: RNA3DB's 15,441 chains collapse to **142
structurally independent components**. The whole RNA corpus fits in under 100 GB unless we
build MSAs, which adds 933 GB for NCBI nt — and at 142 independent components, it is hard
to argue the MSAs are worth 933 GB.


---

## Revision 7 — 2026-09-20, correcting Revision 3

### AlphaGenome is out entirely, not merely non-commercial

Revision 3 listed **AlphaGenome-PyTorch as recommended model #2** and described the
licence issue as a restriction on training. That was too soft, and the recommendation is
**withdrawn**.

AlphaGenome's Model Terms state its outputs "should not be used for the training of other
machine learning models." A JEPA head trained on AlphaGenome embeddings is exactly that.
So it is unusable as a feature source *and* as a distillation target — the two ways we
would plausibly have used it. It survives only as an **external comparator**: we may
compare against its published numbers, not build on its outputs.

Revised DNA shortlist:

| slot | model | licence | dims |
|---|---|---|---|
| 1 | **Evo 2 7B** (`arcinstitute/evo2_7b`) | Apache-2.0 | 4096 @1 Mb ctx |
| 2 | **borzoi-pytorch** (`johahi/borzoi-replicate-*`) — port claims exact parity | CC-BY-4.0 | ~1920 @32 bp |
| 2b | **Enformer** | CC-BY-4.0 | 3072 @128 bp |
| 3 | **ModernGENA-large** | open | — |

### The methodological finding, which matters more than the model list

**GUE and the Nucleotide Transformer benchmark share roughly half their datasets. BEND,
DART-Eval and LRB are largely ENCODE again.**

> Scoring well on three suites is one result reported three times.

This is the same disease as leaky splits, one level up: not leakage between train and
test, but leakage between *benchmarks*, so apparent independent confirmation is
correlated by construction. It generalises beyond genomics and belongs in this project's
standing rules — we already require difficulty-matched benchmarks (FINDING 022) and
honest splits (Revision 2); add **provenance-disjoint benchmarks** to that list.

Before claiming a result holds "across N benchmarks", check how many distinct datasets
those N actually contain.

### Sizes, now measured rather than estimated

| resource | size | verdict |
|---|---|---|
| ENCODE narrowPeaks: TF 71.0 GB · histone 44 GB · DNase 9 GB · ATAC 7 GB | ~131 GB | **peaks only** |
| ENCODE bigWigs | 43 TB | no |
| all of ENCODE | 1.71 PB | no |
| **SCREEN cCRE V4** — 2,348,854 elements | **129 MB** | best value-per-byte in the survey |
| HT-SELEX Jolma 2013 | 10.76 GB | yes if we go this way |
| HT-SELEX Yin methyl | 507.7 GB | no |
| Codebook GHT-SELEX (2026) — raw 235.7 GB / **peaks 73.3 MB** | peaks | best new TF-DNA resource since 2017; real genomic flanks |

### A judgment call for the user, not for me

**Roadmap Epigenomics (1.2 TB) sits on a single unmirrored server behind a dead domain.**
The agent's recommendation is "mirror now or lose it."

I am not doing that unilaterally. It is 1.2 TB of a *shared* 1.5 PB university filesystem,
spent on data this project has no use for — DNA contributes essentially no interaction
supervision (Revision 3), and epigenomics contributes none at all. Preserving it may
well be a public good, but it is a public good at someone else's expense and it is not
our call to make quietly. Flagged for the user to decide.

### New blockers

DART-Eval requires a Synapse DUA — worth noting that **its task 5 is the only
held-out-population variant test found anywhere in the survey**, which is exactly the
kind of honest evaluation this program keeps saying it wants. Enformer/Borzoi training
buckets are requester-pays (~$600 for 5 TB) and unnecessary.
