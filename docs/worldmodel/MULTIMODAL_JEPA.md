# Multimodal interaction models, JEPA, and honest proteochemometrics

**Research dossier, 2026-09-20.** Compiled for the shared-embedding world model
(`docs/worldmodel/README.md`, ask 1). Three parts: (1) what already exists for joint
protein-ligand / multi-biomolecule embedding, (2) what JEPA is and whether biology has
actually adopted it, (3) proteochemometrics and — the part that matters most — the
evaluation traps that make most of the DTI literature unbelievable.

**Reading posture.** Every number below is someone's self-report. Where I can name the
split, the negative set, or the ablation that would falsify a claim, I do. Where I
cannot, the claim is marked unverified rather than repeated as fact.

---

## Part 1 — Existing multimodal biological interaction models

### 1.1 ATOMICA (Zitnik lab) — the one the user named

| field | value |
|---|---|
| name | ATOMICA — *Learning Universal Representations of Intermolecular Interactions* |
| authors | Ada Fang, Michael Desgagné, Zaixi Zhang, Andrew Zhou, Joseph Loscalzo, Bradley Pentelute, Marinka Zitnik (HMS / MIT) |
| year | bioRxiv v1 2025-04-02, v2 later 2025. **Still a preprint as of 2026-09-20** — 17 months, no journal version indexed in Europe PMC |
| paper | https://www.biorxiv.org/content/10.1101/2025.04.02.646906 · PMC12026499 |
| code | https://github.com/mims-harvard/ATOMICA — **MIT license** ("The code in this package is licensed under the MIT License") |
| weights | https://huggingface.co/ada-f/ATOMICA — **CC-BY-4.0**, genuinely downloadable |
| size | **not disclosed.** The hyperparameter search reports 4–8 tensor-field layers, node dim 16–32, edge dim 16–32, k∈{4,8,16} neighbours. This is a *small* model — order 10⁶ params, not 10⁸. Treat the "foundation model" framing accordingly. |

**What it embeds.** Not a protein and not a ligand — an **interaction interface**. Atoms
within 8 Å of the partner molecule, from both sides, grouped into "chemical blocks"
(amino acids, nucleotides, small-molecule functional moieties). Five modalities: proteins,
small molecules, metal ions, lipids, nucleic acids. Eight interaction types.

**How it fuses.** It does not fuse two encoders — that is the architectural point. One
SE(3)-equivariant tensor-field network runs over a single graph containing *both*
partners, with intramolecular edges (k-NN within a molecule) and **intermolecular edges**
(k-NN across the interface). Hierarchy: atom-level TFN (λmax=1, filtered to λmax=2) →
multi-head cross-attention pooling into block nodes → block-level TFN → multi-head
self-attention graph pooling. Embeddings come out at atom, block and interface level.

**Pretraining objective.** Self-supervised, two terms:

- *Denoising* — perturb the complex by rigid rotation (σ∈{0.25,0.5,1}), translation
  (σ∈{0.5,1,1.5}) and torsion-angle noise on rotatable bonds (σ∈{0.25,0.5,1}); regress the
  score. This is diffusion-style denoising score matching, **not** a JEPA loss.
- *Masked block identity* — mask 10% of blocks, cross-entropy on block type.
- Combined `L = β_ω·l_ω + β_t·l_t + β_θ·l_θ + β_m·l_m`, with β_ω=β_t=1, β_m=0.1.
- 150 epochs, cosine annealing with warm restarts between 1e-4 and 1e-6, cycle 400k steps.

**Training data.** 2,037,972 complexes:

- 1,747,710 small-molecule pair complexes from the **CSD** (organic, non-polymeric,
  6–50 heavy atoms, no disorder, no metals) — i.e. **86% of the corpus is small-molecule
  crystal packing, not biology.** Hold onto this: the "universal" claim rests mostly on
  organic crystal contacts.
- 290,262 from PDB / Q-BioLiP: 124,541 protein–protein, 119,017 protein–small-molecule,
  74,514 protein–ion, 8,475 protein–peptide, 5,185 nucleic-acid–ligand, 3,511 protein–RNA,
  2,750 protein–DNA.
- Splits clustered at **30% protein sequence identity with MMseqs2**, plus chemical-similarity
  sampling for val/test (10,000 each). *This is a genuinely careful split* — better than
  most of the DTI literature in Part 3.

**What it was evaluated on, and how honest it looks.**

- *Zero-shot interface-residue identification (ATOMICAScore).* 2.7 of 10 correct amino-acid
  blocks retrieved at rank 10, vs ESM-2-3B at 2.4 and random at 2.0, over 5,691 complexes.
  **Believable but tiny.** Random is 2.0; ATOMICA is 0.7 above random and ESM-2 is 0.4
  above. That is a weak signal presented as a win over a 3B-parameter baseline.
- *Cross-modality transfer (masked-block AUPRC).* Multi-modality vs single-modality
  training: protein–DNA 0.24→0.71, protein–RNA 0.19→0.55, protein–peptide 0.32→0.67.
  **This is the most convincing result in the paper** and the one directly relevant to our
  broad-vs-narrow experiment: for *data-poor* modalities (2–7k complexes), training on
  everything else is worth 2–3× AUPRC. It is also the result most exposed to the obvious
  confound that the single-modality baselines are simply undertrained.
- *ATOMICANets disease-pathway analysis.* Connected-component sizes with permutation
  p-values. I would not build on this: it is network-biology enrichment, sensitive to
  degree distribution, and the effect sizes are small (largest component of 7 proteins,
  p = 0.037).
- *Dark proteome.* 2,646 predicted ligand-binding sites; **AlphaFold3 ipTM used as a proxy
  validator** (KS 0.11 ions, 0.54 small molecules). Using a co-folder's confidence to
  validate a binding prediction is circular, and this repo's FINDING 009/011 already
  showed that co-folder confidence ranks nothing. Discount this section.
- *Wet lab.* **5 of 6 expressed proteins confirmed heme binding**; four predicted covalent
  (scores 0.997, 0.690, 0.992, 0.874), and one non-canonical binder without a CXXCH motif.
  This is the real result. n=6, no negative controls reported, so it is a demonstration
  rather than a rate — but it is an actual experiment, which is more than almost any model
  in this dossier can say.

**Stated limitations (theirs).** Dependence on high-resolution structures; IDRs and
antibody CDRs poorly captured; structural datasets small relative to proximity-based
screens (Y2H, AP-MS, DEL, FRET, SPR).

**Can we use it? — the load-bearing constraint.** ATOMICA takes **PDB files of an
already-bound complex** (`atomica.data.process_pdbs`). It embeds a *pose*. It cannot take
a protein plus a SMILES and say whether or how they bind, because the 8 Å interface is
only definable once a pose exists. Consequences:

- ✅ **Directly usable as a pose scorer.** Embed the interface of each Boltz/Protenix
  CYP3A4 pose and use the embedding (or a light head on it) as a selection feature. This
  is a drop-in candidate against `cypstruct.xengine.select()`'s +0.0381 and, unlike Boltz
  `complex_ipde`, it was never trained to be a confidence — so it is not pre-doomed the way
  FINDING 011 showed confidence to be.
- ✅ ATOMICA-Ligand ships fine-tunes for 9 metal ions and 12 small molecules, **heme
  included**. That is the single most on-topic public checkpoint for a CYP project, and the
  wet-lab validation was specifically heme binding.
- ❌ **Not usable as the world model itself.** There is no unbound-protein encoder to query,
  so it cannot answer "what should bind here". To get a partner prior we would have to build
  the missing half.
- ⚠️ Because it is pose-conditioned, any scorer built on it inherits a leak risk: the
  embedding may encode "does this look like a deposited crystal structure" rather than "is
  this the right pose". Pre-register the ablation — score native poses against deliberately
  rotated/translated ones; if the score separates those trivially, it is a crystallinity
  detector, not a pose scorer. (Its pretraining objective is literally denoising rotations
  and translations, so this is not a hypothetical concern — it is what it was trained to do.)

### 1.2 Multi-modality models (three or more modalities)

| model | year | URL | what it does | weights / license |
|---|---|---|---|---|
| **ATOMICA** | 2025 | above | 5 modalities, one interface graph | HF, CC-BY-4.0 (code MIT) |
| **BioBridge** | ICLR 2024, arXiv 2310.03320 | Wang, Wang, Srinivasan, Ioannidis, Rangwala, Anubhai (AWS) | does **not** retrain encoders. Learns small transformations *between* frozen unimodal FMs (protein LM, molecule LM, text LM), supervised by a biomedical **knowledge graph**. Evaluated on cross-modal retrieval + biomedical QA + drug discovery | code public; confirm license before use |
| **NatureLM** | 2025-02, arXiv 2502.07527 | Microsoft Research | one sequence model over small molecules, proteins, DNA, RNA and materials; supports protein→molecule and protein→RNA *generation* | gated; verify |
| **Central Dogma Transformer** | 2026-01, arXiv 2601.01089 | — | directional cross-attention stitching pretrained DNA, RNA and protein models into a "Virtual Cell Embedding" | unverified |
| **ProtST** | ICML 2023 oral, arXiv 2301.12040 | Xu, Yuan, Miret, Tang | protein sequence ↔ biomedical **text**; alignment + multimodal masked prediction; zero-shot protein classification and retrieval | code public |
| **Mol-JEPA** | 2026-08 | Part 2 — 14 modalities, the most modality-rich molecular model found | **weights + code released, CC-BY-4.0** |
| **ProtJEPA** | 2026-08 | Part 2 — 10 protein modalities | CC-BY-**NC** |

**UniBioseq: I could not verify that it exists.** Zero hits on arXiv (`all:UniBioseq`) and
zero on Europe PMC. Either the name is wrong, it lives on a venue neither indexes, or it
was seen in a talk. Flagged in BLOCKED — do not plan around it until the user confirms.

**BioBridge is the architecturally interesting one for us.** It is the cheapest possible
version of ask 1: keep every frozen foundation encoder, learn only the maps between them,
supervise with a KG. Given "1× contended H200, no from-scratch pretraining", that is close
to the only tractable shape. Its weakness is that it needs a KG edge to supervise each
pair, so it inherits the KG's coverage bias — and for CYP3A4 substrate chemistry, KG
coverage is exactly what is missing. A KG will tell you CYP3A4 metabolises midazolam. It
will not tell you what CYP3A4 *should* bind.

### 1.3 The two-tower / shared-space protein–ligand family

The user's question — joint space vs concatenation — splits this field three ways:

1. **True shared space (contrastive co-embedding).** Two towers, one metric space, trained
   so that binding pairs are close. You can *retrieve* with these: embed a pocket, nearest-
   neighbour into a molecule library. ConPLex, DrugCLIP, S²Drug, AANet, ConGLUDe.
2. **Interaction-map fusion.** No shared space; instead an explicit pairwise tensor over
   substructures/residues (bilinear attention, cross-attention). DrugBAN, MolTrans,
   TransformerCPI, PerceiverCPI, HyperAttentionDTI. More expressive per pair, but you
   cannot query "what binds here" without scoring the whole library.
3. **Single joint graph.** ATOMICA. Most expressive, requires a pose.

Only class 1 gives the generative-prior-over-binders object that `README.md` describes.
Class 2 is a scorer. That is the key architectural takeaway of Part 1.

Recent class-1 entries worth reading:

- **DrugCLIP** (arXiv 2310.06367, Gao et al.) — reframes virtual screening as *dense
  retrieval*, aligning pocket and molecule encoders contrastively with no affinity labels.
  The reframing is right; the evaluation is the problem (see Part 3 on DUD-E decoy bias).
- **AANet** (arXiv 2506.05768, 2025) — **tri-modal** contrastive over ligand, pocket and
  cavity, explicitly built for *structural uncertainty*. Closest published spirit to our
  pose-ambiguity problem.
- **ConGLUDe** (arXiv 2601.09693, 2026) — aligns ligands to *multiple* candidate binding
  sites rather than assuming one pocket. Relevant to CYP3A4, which has a large promiscuous
  cavity and multiple sub-sites.
- **S²Drug** (arXiv 2511.07006, 2025) — fuses protein sequence and structure at residue
  level before contrastive alignment with ligands.
- SE(3)-equivariant contrastive pocket/ligand shared space (arXiv 2604.19562, 2026).

Per-model detail cards — weights, licenses, benchmarks and split honesty for ConPLex,
DrugBAN, PerceiverCPI, MolTrans, HyperAttentionDTI and TransformerCPI — are in §3.9.

---

## Part 2 — JEPA, and whether biology has actually adopted it

### 2.1 The core architecture

A JEPA predicts **representations of masked content from representations of visible
content**, in latent space. It never reconstructs the input. The claimed advantage is that
the encoder is free to discard unpredictable detail (pixel noise, side-chain jitter,
dropout in scRNA-seq) instead of being forced to model it — which is exactly the failure
mode of masked autoencoders and, arguably, of coordinate-space diffusion.

| model | year | arXiv | what |
|---|---|---|---|
| **I-JEPA** | 2023 | 2301.08243 | Assran, Duval, Misra, Bojanowski, Vincent, Rabbat, LeCun, Ballas. Predict representations of target blocks in an image from one context block. ViT-H/14 trained on ImageNet in <72 h on 16 A100s — the efficiency claim is the headline. |
| **MC-JEPA** | 2023 | 2307.12698 | Bardes, Ponce, LeCun. Motion (optical flow) + content in one encoder. |
| **A-JEPA** | 2023 | 2311.15830 | Fei et al. Audio; time-frequency-aware curriculum masking. |
| **Point-JEPA** | 2024 | 2404.16432 | Saito et al. Point clouds, no reconstruction, no auxiliary modality. |
| **V-JEPA 2** | 2025 | 2506.09985 | Meta. **>1 million hours of internet video.** 77.3 top-1 Something-Something v2; 39.7 R@5 Epic-Kitchens-100 anticipation; 84.0 PerceptionTest / 76.9 TempCompass at 8B params aligned to an LLM. V-JEPA 2-AC adds action conditioning from **<62 h** of Droid robot data and does zero-shot pick-and-place from goal images. Weights on HF/GitHub. |
| **LeJEPA** | 2025-11 | 2511.08544 | **Balestriero & LeCun.** The one to actually read. Proves the isotropic Gaussian is the optimal embedding distribution for downstream prediction risk, and enforces it with **SIGReg** (Sketched Isotropic Gaussian Regularization). Removes stop-gradient, EMA teacher and the rest of the heuristic stack; one trade-off hyperparameter, linear time and memory, ~50 lines. 79% ImageNet-1k frozen eval with ViT-H/14. |

**LeJEPA matters to us specifically** because the heuristic stack (EMA teacher, stop-grad,
careful masking schedules) is where low-data JEPAs collapse, and we will be low-data. A
method with one hyperparameter and a proof is the right starting point for a team with one
contended GPU. Note that Mol-JEPA below already adopted SIGReg and **dropped the EMA
teacher entirely** — so this is not theoretical.

**The world-model framing.** LeCun's argument is that a system that predicts in
representation space, conditioned on an action, can *plan* by searching actions whose
predicted latent state matches a goal. V-JEPA 2-AC is the existence proof. The analogue
our `README.md` wants — "given entity A and a context, predict the embedding of the entity
B that binds it" — is a **conditional** JEPA where the "action" is the choice of partner.
Nobody has built that for binding. That is the gap, and it is a real one.

### 2.2 Does bio-JEPA exist? — **Yes, and more than expected. It is ~12 months old.**

This is the honest answer to the user's question, and it is the opposite of what I expected
going in. As of September 2026 there is a small but real bio-JEPA literature, concentrated
in the last twelve months, with working code in several cases.

#### Molecules

**Mol-JEPA** — arXiv 2608.22642, Rottach, Schieferdecker, Rudman, Balestriero, Eickhoff
(**Boehringer Ingelheim** + Balestriero, i.e. the LeJEPA author). Submitted 2026-08-23.
Code and checkpoints: https://github.com/Boehringer-Ingelheim/mol-jepa · **CC-BY-4.0**.
**This is the most directly relevant paper in the entire dossier.**

- **14 modalities**, each from a pretrained encoder, masked at the *modality* level:
  UMA (atomistic FM, 128d), ChemGPT (2048d), CLOOME (cell painting, 512d), BioXMol
  (multi-assay phenotype, 1024d), **Boltz-2 embeddings (6912d)** and Boltz-2 predictions
  (4d), graph transformer (512d), MOE descriptors (227d), ECFP (2048d), xTB (7d),
  ∇²DFT (17d), ChEMBL (306d), PCBA (1328d), TDC ADMET (672d).
- **~50M params.** Transformer predictor. **No EMA teacher** — collapse prevented by
  SIGReg. Modalities masked Bernoulli(r), at least one kept.
- 4.69M unique compounds (ChEMBL 268k, TDC 403k, PCBA 1.12M, ∇²DFT 1.27M, Enamine REAL 2.01M).
- Results (MAE, lower better) across 23 datasets: ExpansionRx 9-endpoint ADMET **0.402** vs
  TabICLv2 0.453 / CheMeleon 0.471; ASAP-Polaris **0.524** vs 0.609 / 0.703; Biogen ADME
  **0.337** vs 0.369 / 0.384; **PXR activity 0.55** vs TabICLv2 0.87 / RF 0.81.
- Leave-one-modality-out: graph contributes most, then MOE and ECFP; the experimental
  assay modalities (ChEMBL/PCBA/TDC) contribute meaningfully. Performance rises
  monotonically with modality count; training time is linear in it.

**Skeptical read.** The paper's own Wilcoxon signed-rank analysis reports **38% of
comparisons favouring Mol-JEPA vs 47% favouring TabICLv2** — i.e. by their own significance
test, a tabular in-context-learning baseline wins more head-to-heads than the JEPA does,
even though the JEPA wins the headline means. That is unusually honest reporting and it is
the number to quote, not the MAE table. Their own limitations note that most ablations were
run on a much smaller dataset than the full pretraining corpus. Verdict: **a real, usable,
openly-licensed multimodal molecular JEPA whose advantage over a strong tabular baseline is
not established.** Note also that it uses Boltz-2 embeddings as a modality — which is
directly compatible with what this repo already generates — and that it benchmarks on PXR
activity, the sibling repo's exact task.

#### Proteins

**ProteinJEPA** — arXiv 2605.07554, Dan Ofer, Dafna Shahaf, Michal Linial (Hebrew
University), 2026-05-08. *"Latent prediction complements protein language models."*
Best variant is **masked-position MLM+JEPA**: predict latent targets only at masked
positions while keeping the MLM cross-entropy. On ESM2-35M: 10 wins / 3 losses / 3 ties;
on ESM2-150M: 11 / 2 / 3, across 16 downstream tasks (15 frozen linear probes + SCOPe-40
zero-shot fold retrieval; stability, β-lactamase fitness, variant effect, disorder, remote
homology, EC classification, fluorescence). Code/weights not stated in the abstract.

**Skeptical read.** The honest framing is in the title: *complements*, not replaces. Pure
JEPA did not beat MLM — the hybrid did, and only at 35–150M params, which is two orders
below ESM2-3B. A win count of 11/2/3 on frozen linear probes with no reported seed variance
is weak evidence; linear-probe deltas of this kind are routinely within noise. **I would
not claim JEPA beats MLM for proteins on this evidence.** I would claim it does not hurt,
which is itself mildly interesting.

**ProtJEPA** — bioRxiv 10.64898/2026.08.03.742606, Vaibhava Lakshmi Ravideshik, Jinha Kim,
**Manolis Kellis** (MIT), 2026-08. **CC-BY-NC** (no commercial use). A *sequence-only
student* encoder trained to predict joint embeddings spanning **ten modalities**: sequence,
structure, knowledge graph, protein interactions, literature, localization, tissue
expression, GO function, anatomy, disorder. Multi-teacher, modality-attentive fusion, with
**target whitening** to kill anisotropy and collapse. On held-out "dark" proteins: 58.07%
Hit@10 zero-shot GO retrieval, 69.99% enzyme-class accuracy, improved subcellular
localization; gains reported to transfer to **drug-target interaction** and disorder
prediction.

**Skeptical read.** This is architecturally the closest published thing to ask 1 — distil
many modalities into one queryable sequence-only space — and it is the right idea. But:
it is a preprint six weeks old, NC-licensed, the "dark protein" held-out set is not a
standard benchmark, and the DTI transfer claim is exactly the kind of claim Part 3 says to
disbelieve without a cold-target split. Read it for the architecture (multi-teacher +
target whitening), not for the numbers.

#### Genomics

**JEPA-DNA** — arXiv 2602.17162 (v3), Larey, Dahan, Bleiweiss, … Yoli Shavit (a large
Intel/Sheba/Mount Sinai collaboration), 2026-02. Model-agnostic **continual pretraining**:
it upgrades an *existing* genomic FM rather than training one. Applied to DNABERT-2, NTv3
and HyenaDNA it improved 6–8 of 9 supervised tasks regardless of architecture, tokenization
or original objective; SOTA across 17 genomic benchmarks. Code on GitHub. The stated
motivation — the "granularity trap", where token-level objectives over-fit local nucleotide
patterns and miss function — is the same argument we would make for interfaces over
coordinates.

**Verdict: this is the strongest bio-JEPA result in the dossier.** Not because the numbers
are biggest, but because the experiment is the right one: same backbone, same data, one
objective swapped, three architectures, consistent direction. That is a controlled
comparison, and almost nothing else here is.

**GenoJEPA** — 2026-04, built on **LeJEPA**. Reported to outperform NT-v2 across 55 tasks
at **52M vs 494M params** (~10× smaller). Not found on arXiv; indexed via the HF
"JEPA in Bioscience" collection. Unverified primary source — flagged in BLOCKED.

#### Single cell

**GeneJEPA** — bioRxiv 2025.10.14.682378, BiostateAI. Code
https://github.com/BiostateAI/GeneJEPA, weights `elonlit/GeneJEPA` on HF. Predicts
representations of masked gene sets from visible context; trained on the **Tahoe-100M**
atlas (100M cells); reported to surpass scGPT and Universal Cell Embeddings on drug-response
and genetic-perturbation prediction. ⚠️ **No license stated in the repo**, no parameter count,
and **no quantitative benchmarks in the README** — the numbers live only in the preprint.
Treat as unusable until the license is clarified.

**Cell-JEPA** — arXiv 2602.02093, ElSheikh, Wang, Wu, … Aly Khan, Han Liu, 2026-02.
**0.72 AvgBIO cell-type clustering vs scGPT 0.53** (+36% relative, zero-shot). But
**mixed on perturbation effect-size metrics** — the authors say so. The pattern is worth
noting: JEPA wins the *geometry* task (clustering) and does not clearly win the
*quantitative* task (effect size). If that pattern is general, it predicts what a bio-JEPA
would buy us: better organised embedding space, not better regression.

**BioM-JEPA** — arXiv 2608.05928, Yuhao Wang, Zelin Zang, Yuxuan Liu, Zhen Lei, Stan Z. Li
(Westlake), 2026-08-06. Predicts aggregate representations of **graph-connected gene
blocks** defined by protein-association and coexpression evidence, rather than individual
genes. Student/teacher. Lowest aggregate perturbation-response error among models compared
on CellBench; 5.75× fine-tuning throughput and 3.76× embedding throughput vs scFoundation.
Checkpoint "will be released" — i.e. not yet.

### 2.3 What bio-JEPA has and has not shown

**Has:**
- Latent prediction is a drop-in improvement over token reconstruction for genomics,
  demonstrated across three backbones with one variable changed (JEPA-DNA).
- Modality-level masking scales: more modalities monotonically helps (Mol-JEPA).
- Parameter efficiency is real — 52M beating 494M (GenoJEPA, unverified), 50M competitive
  with tabular SOTA (Mol-JEPA).
- Collapse is now a solved-ish engineering problem: SIGReg (Mol-JEPA) or target whitening
  (ProtJEPA) instead of the EMA-teacher stack.

**Has not:**
- **Nobody has built an interaction JEPA.** Every bio-JEPA above masks *within* one entity
  (genes of a cell, positions of a sequence, modalities of a molecule). None predicts the
  embedding of a *binding partner* from the embedding of a target. The object
  `docs/worldmodel/README.md` describes does not exist in the literature.
- No bio-JEPA has been shown to beat a strong classical baseline decisively. Mol-JEPA's own
  Wilcoxon test goes against it; ProteinJEPA only wins as an MLM hybrid; Cell-JEPA loses on
  half its metrics.
- No bio-JEPA has been evaluated on cold-target splits for anything interaction-related.

**So the answer to "does bio-JEPA exist" is: yes, about a dozen papers, mostly in the last
year, mostly honest, mostly showing modest gains — and the specific thing we want is the
one thing none of them did.** That is a genuine opening, and also a warning: if it were
easy and obviously good, this crowd would have done it.

---

## Part 3 — Proteochemometrics done properly, and the traps

This is the section that matters most, because `docs/worldmodel/README.md` kill criterion 1
is *"if frozen foundation embeddings cannot beat a Morgan-fingerprint + one-hot-target PCM
baseline on a held-out affinity task, the shared space is not carrying biology."* That
criterion is only meaningful if "held-out" is defined correctly. Under a random split it is
unfalsifiable — everything scores 0.9.

### 3.1 The state of the art, and the calibration number

The canonical group is **Gerard van Westen's (LACDR, Leiden)**. The infrastructure:

**Papyrus** — Béquignon, Bongers, Jespers, IJzerman, van der Water, van Westen,
*J Cheminform* 15:3 (2023), doi 10.1186/s13321-022-00672-x. **59,775,087 activity points,
1,270,570 unique 2D compounds, 6,926 proteins, 499 organisms.** Quality tiers:
**1,238,835 high / 335,661 medium / 58,200,591 low** — the high-quality fraction is ~2% of
the headline pile, which is the first thing to know about "60 million datapoints". Their PCM
protein descriptor is already a *learned* one: **UniRep 64+256+1900 concatenated → 6,660-d**,
not z-scales.

**The single most important number in this dossier:**

| split | QSAR | PCM | DNN |
|---|---|---|---|
| **random** | MCC 0.51, r 0.66, RMSE 0.73 | MCC 0.61, **r 0.79**, RMSE 0.75 | MCC 0.60, **r 0.81**, RMSE 0.75 |
| **temporal** (train ≤2013) | MCC 0.20, r 0.24, RMSE 1.19 | MCC 0.29, **r 0.42**, RMSE 1.17 | MCC 0.30, **r 0.36**, RMSE 1.44 |

The *same PCM model* goes from **r = 0.79 to r = 0.42** purely by changing the split.
RMSE 0.75 → 1.17 log units. And the **method ranking flips**: the DNN looks best under random
(r 0.81) and *loses* to PCM under temporal on RMSE (1.44 vs 1.17). Memorise this table before
designing any of our own evaluations.

**QSPRpred** — van den Maagdenberg, Šícho, Araripe et al., *J Cheminform* 16:128 (2024),
doi 10.1186/s13321-024-00908-y, https://github.com/CDDLeiden/QSPRpred. First Python QSPR
toolkit with first-class PCM support and the hard splits already implemented. **This is the
stack to use for our baseline** rather than rolling our own — it removes an entire class of
self-inflicted split bugs.

Other current Leiden work with a directly transferable lesson:

- **3DDPDs** (Gorostiola González et al., *J Cheminform* 2023, doi 10.1186/s13321-023-00745-5)
  — MD-derived *dynamic* protein descriptors beat static descriptors **on temporal splits but
  only match them on random splits**. The descriptor's value is *invisible* under the easy
  split. Direct evidence that descriptor ablations must be run under the hard split, and
  directly relevant to us, since CYP3A4's whole story is F/G-loop dynamics.
- **Assay-aware bioactivity models** (Schoenmaker et al., *JCIM* 2025, doi
  10.1021/acs.jcim.5c00603) — adding explicit assay/biological-context features moved
  average **R² 0.67 → 0.69**; a bag-of-words assay encoding gave **0.66, no improvement**.
  Small, correctly reported, honest.
- **Sparse kinase multitask** (Luukkonen et al., *JCIM* 2023, doi 10.1021/acs.jcim.3c00132)
  — multitask DL beats single-task and tree models on sparse data; **imputation did not
  help**; dissimilarity-clustered splits revealed poor generalisation vs random.
- **Lenselink et al.** *J Cheminform* 9:45 (2017) — the origin of both the "DNN_PCM wins"
  claim and the temporal-split convention. Read it as the thing everything above is
  reacting against.

### 3.2 Do PLM embeddings actually beat z-scales / ProtFP? — No published leak-controlled evidence that they do

This is where marketing and ablations disagree, and it is load-bearing for our kill
criterion 1.

**"All That Glitters Is Not Gold: Importance of Rigorous Evaluation of Proteochemometric
Models"** — Avdiunina, Jamal, Gusev, **Isayev**, *JCIM* 2025, doi 10.1021/acs.jcim.5c00395.
Kinase–ligand bioactivity, multiple protein/ligand descriptors including MSA-augmented
embeddings. Findings: **permutation testing consistently showed that the protein embeddings
contributed minimally**, and **splitting strategy and class imbalance dominated performance**
far more than the choice of protein representation.

**HonestAffinity** (Wei et al., 2026 preprint) — three variants of one architecture on 11,513
LP-PDBBind complexes: ESM-2 650M + pocket marker; ESM-2 only; and a plain **21-token residue
embedding** + pocket marker. Result is a **"split-conditioned reversal"**: ESM-2 + pocket wins
on validation and CASF-2016, but the **no-ESM 21-vocabulary variant wins Pearson R on every
strict LP-PDBBind no-leak tier (test_cl1–cl3)**. Both the pocket prior and the 1280-d ESM-2
input **help on leaky splits and hurt on strict ones.**

> **Defensible summary:** there is no published, leak-controlled ablation showing PLM
> embeddings beating z-scales/ProtFP for PCM. The two studies that ran the ablation properly
> — permutation testing, and leak-proof tiers — found the protein representation contributes
> little, or that the *simpler* representation wins once leakage is removed. Papers claiming
> "ESM embeddings improve DTA" almost universally report only random/warm splits and do not
> ablate against a classical descriptor.

**This is a direct hit on our program's premise.** Ask 1 assumes frozen foundation embeddings
carry biology that a fingerprint does not. Two careful papers say that under honest
evaluation they may not. Kill criterion 1 is therefore not a formality — it is the most likely
outcome, and we should expect to have to *earn* the embedding, not assume it.

### 3.3 Does deep learning beat gradient boosting for PCM? — Not under hard splits

- **DrugBAN** (*Nat Mach Intell* 2023), the authors' own words: *"RF achieves good performance
  and even consistently outperforms other deep learning baselines (DeepConv, GraphDTA and
  MolTrans) on the BindingDB dataset… deep learning methods are not always superior to shallow
  machine learning methods under the cross-domain setting."*
- **Schapin, Navarro, Bou, De Fabritiis** (Acellera), arXiv:2407.19073 (2024): *"simpler models
  can surpass more complex ones in specific tasks"*; 3D models only become competitive with
  large multi-target labelled datasets.
- **Janela & Bajorath**, *Sci Rep* 13 (2023) doi 10.1038/s41598-023-45086-3 and *Pharmaceuticals*
  16:530 (2023). Across **367 target-based activity classes**, potency predictions are
  dominated by compounds near the class median, whose potency is *"consistently predicted with
  high accuracy, without the need for learning."* ML and trivial control models score
  comparably and results are *"surprisingly resistant to modifications."* Their conclusion:
  **conventional benchmark settings are unsuitable for directly comparing potency methods.**
  **If a paper reports RMSE on a ChEMBL activity class without a median baseline, the number is
  uninterpretable.** This includes any number we would produce.

### 3.4 Trap 1 — the split. How far do scores fall from random to cold-pair?

**The cleanest single-paper ladder** — Bai, Miljković, Ge, Greene, John, Lu, *Hierarchical
Clustering Split for Low-Bias Evaluation of DTI Prediction*, IEEE BIBM 2021. One model
(MolTrans), one de-biased BindingDB (29,674 positives / 32,752 negatives, **all negatives
experimentally validated at IC50 > 10 µM** against positives < 100 nM):

| split | AUROC | AUPRC | accuracy | **sensitivity** |
|---|---|---|---|---|
| random | **0.946** | 0.935 | 0.874 | 0.838 |
| cold drug | 0.921 | 0.909 | 0.841 | 0.798 |
| scaffold | 0.893 | 0.874 | 0.804 | 0.736 |
| HDBSCAN density clustering | 0.821 | 0.778 | 0.724 | 0.581 |
| single-linkage clustering | **0.768** | 0.717 | 0.676 | **0.483** |

**−18.8% AUROC, −23.3% AUPRC** random → single-linkage. Sensitivity collapses 0.838 → 0.483
while specificity holds near 0.89 — **the model degenerates toward "predict negative"**, which
a headline AUROC hides completely. Always report sensitivity separately.

Two second-order findings that matter more than the headline:

- **The gap between methods shrinks to nothing.** MolTrans' advantage over DeepDTA falls from
  **4.4% → 1.1%**. And **DeepConv-DTI beats DeepDTA under random split and loses to it under
  HDBSCAN** — *the leaderboard order is an artefact of the split.*
- **91% of drugs in binary BindingDB have only one pair type** (all-positive or all-negative).
  That alone permits a protein-blind classifier to do well. They had to delete those drugs to
  build a low-bias set.

**DrugBAN's cold-pair numbers**, same architecture, both splits:

| dataset | random | clustering-based pair split |
|---|---|---|
| BindingDB | **0.960** | **0.575** |
| BioSNAP | **0.903** | **0.654** |

Domain adaptation recovers only to 0.604 / 0.684. On **Human** under a random split *every*
deep model scores **AUROC > 0.98**; under cold-pair they all drop sharply.

**Rule of thumb to carry into our own work:** expect **−0.15 to −0.25 AUROC** from random to
cluster/cold-pair on classification, and **roughly half the correlation** (r 0.79 → 0.42) on
regression. A paper reporting AUROC ≥ 0.95 on a random split has told you approximately
nothing about prospective performance.

### 3.5 Trap 2 — the leakage channels that make random splits lie

1. **Near-duplicate ligands.** *Data Leakage and Redundancy in the LIT-PCBA Benchmark*,
   arXiv:2507.21404: **2,491 unique inactives shared between train and validation**; **>350
   active train–validation analogue pairs at ECFP4 Tanimoto ≥ 0.6** (ALDH1 alone: 323, many
   >0.9); 2,945 duplicated inactives within training. The payoff: **a trivial ECFP4
   nearest-neighbour memorisation baseline reaches median EF1% = 4.15, matching
   state-of-the-art 3D encoders.** The authors state LIT-PCBA cannot be repaired without
   discarding most of it.
2. **Protein homology.** A cold-target split that picks proteins at random still leaves close
   homologues on both sides. **Schuh, Daniluk & Sieber**, *Chem Sci* 2026, doi
   10.1039/d6sc01799a, audited **>50 dataset configurations** and found *"pervasive cross-split
   contamination, unresolved label conflicts, and severe structural redundancies,"* worst *"in
   drug-target interaction datasets with extensive protein overlap"* — and showed
   **leading-model conclusions change when the audited test-set composition is altered.**
   (ATOMICA's MMseqs2 30%-identity clustering is the right remedy; most DTI papers do nothing.)
3. **Same ligand, different target.** The ligand is in train against target A and in test
   against target B. The model recognises the ligand, not the pair. This is exactly what the
   **label-reversal** experiment was invented to defeat.
4. **Structural train/test overlap.** Graber et al., *"Resolving data bias improves
   generalization in binding affinity prediction,"* *Nat Mach Intell* 2025, doi
   10.1038/s42256-025-01124-5: **~600 train–test similarities between PDBbind and CASF,
   affecting 49% of all CASF complexes.** Retrained on their de-leaked **CleanSplit**:
   Pafnucy Pearson **0.808 → 0.746**, RMSE **1.046 → 1.484**; GenScore 0.815 → 0.780,
   RMSE 1.273 → 1.362.

### 3.6 Trap 3 — the negative set, and the shortcut it creates

The mechanism: build negatives by random drug–protein pairing. Positives come from bioactivity
databases (drug-like, ChEMBL-shaped); negatives come from a different distribution or different
degree statistics. The model reaches near-perfect AUROC by answering **"is this a drug-like
molecule at all / is this a well-annotated hub node"** instead of "do these two bind".

- **Human and C. elegans** (Liu et al. 2015) negatives came from an **in-silico screening
  framework** — the negative labels are *another model's output*, never measured, drawn from a
  different chemical distribution. The tell: in a 2024 benchmark of **26 DTI models**
  (arXiv:2407.04055), **every single model scores ROC-AUC 0.905–0.980 on Human and 0.910–0.987
  on C. elegans.** Twenty-six architectures inside a 0.07 band, all above 0.90, is not
  measuring architecture — it is measuring the benchmark.
- **DUD-E decoys** are property-matched but topologically dissimilar *by construction*, so they
  are separable with no protein at all.
- **AVE bias** — Wallach & Heifets, *"Most Ligand-Based Classification Benchmarks Reward
  Memorization Rather than Generalization,"* *JCIM* 58:916 (2018), arXiv:1706.06619.
  Asymmetric Validation Embedding measures train–validation redundancy **accounting for
  clumping among inactives as well as actives**. Across **seven** widely used VS benchmarks,
  AVE bias **correlates strongly with measured model performance, irrespective of the predicted
  property, fingerprint, similarity measure, or previously applied unbiasing technique.** The
  reported score is predictable from dataset redundancy alone.
- **AI-Bind's annotation imbalance** is the numeric anatomy: of **10,416 ligands in
  DeepPurpose's BindingDB training data, only 793 have both binding and non-binding
  annotations**; of 1,391 proteins, only 667 do. Everything else is classifiable from identity.

### 3.7 Trap 4 — the shortcut: the protein-blind ablation, and why it is the only control that matters

Five independent papers, five different ablations, one verdict.

1. **Chen, Cruz, Ramsey, Dickson, Duca, Hornak, Koes, Kurtzman** — *"Hidden bias in the DUD-E
   dataset leads to misleading performance of deep learning in structure-based virtual
   screening,"* **PLOS ONE** 14(8):e0220113 (2019). The definitive protein-blind ablation.
   They retrained the CNN **replacing the entire protein with a single dummy atom**: across 102
   DUD-E targets, receptor-ligand and **ligand-only models give AUCs correlated at R² = 0.98,
   slope 0.99, both averaging AUC 0.98.** Feeding the *trained* receptor-ligand model ligands
   with vs without the protein present: **mean AUC 0.98 both ways, mean absolute difference
   0.0006.** Cross-target: **for 74 of 102 targets, actives were separated from decoys at
   AUC > 0.9 by models trained on a completely different target's ligands.** And AutoDock Vina
   performed comparably to Gnina, both beating Pafnucy.
2. **Chatterjee et al., AI-Bind**, *Nat Commun* 14:1989 (2023), doi 10.1038/s41467-023-37572-z:

   | model | transductive AUROC/AUPRC | semi-inductive (unseen targets) | inductive (unseen nodes) |
   |---|---|---|---|
   | DeepPurpose | 0.82 / 0.48 | 0.76 / 0.70 | **0.61 / 0.43** |
   | **configuration model** (degree only — no SMILES, no sequence) | **0.83 / 0.50** | **0.77 / 0.71** | 0.50 / 0.30 |
   | MolTrans | 0.952 / 0.887 | — | **0.572 / 0.432** |
   | AI-Bind VecNet | — | — | 0.75 / 0.718 |

   A null model seeing **only bipartite network degrees and no chemistry whatsoever** matches or
   beats DeepPurpose transductively and semi-inductively. Then the kill shot: **randomly
   reshuffling the SMILES *and* the amino-acid sequences in the training set changed AUROC from
   0.86 to 0.84 and AUPRC from 0.64 to 0.62.** The chemistry inputs were decorative.
3. **TransformerCPI** (*Bioinformatics* 2020) — the **label-reversal** experiment. CPI-GNN's
   protein-CNN weights were *"significantly concentrated in zero"*, and ligand-only training
   gave statistically indistinguishable performance (p > 0.05). On the **Kinase label-reversal
   set the reference models fell below AUC 0.5** — worse than chance — while TransformerCPI held
   ~0.75.
4. **Volkov, Turk, Drizard, Martin, Hoffmann, Gaston-Mathé, Rognan** — *"On the Frustration to
   Predict Binding Affinities from Protein–Ligand Structures with Deep Neural Networks,"*
   *J Med Chem* 65(11):7946 (2022), doi 10.1021/acs.jmedchem.2c00487. Modular message-passing
   GNNs over ligand, protein and interaction graphs: the **ligand-only model scored Pearson
   0.749 / RMSE 1.567, beating the protein–ligand interaction-graph model at 0.687 / 1.605.**
   Explicit non-covalent interaction description gave *no* advantage.
5. **Graber et al. 2025** — a model **trained on PDBbind with no protein information at all
   achieved RMSE 1.424, outperforming AutoDock Vina**; pure ligand memorisation reached
   RMSE 1.539 / Pearson 0.707. **Only after de-leaking did the protein-blind model degrade**
   (to 1.572) — i.e. the protein term only starts to matter once the leakage is gone.
6. **Sieg, Flachsenberg & Rarey** — *"In Need of Bias Control,"* *JCIM* 59(3):947 (2019), doi
   10.1021/acs.jcim.8b00712. Reimplemented published CNN scoring functions as **ligand-only
   variants**; concluded *"bias is learned implicitly and unnoticed from standard benchmarks"*
   and that not every benchmark dataset is suitable. Guidelines + code:
   https://github.com/rareylab/MLValidation

### 3.8 Best-practice references

- **Suay-García & Falcó**, *Brief Bioinform* 2026, doi 10.1093/bib/bbag370 — an explicit
  evaluation playbook and reporting checklist: minimum split standards (scaffold, temporal,
  target-wise/cluster), metrics including early-recognition for VS, **required baselines and
  required ablations**, plus task → protocol → failure-mode tables. Closest thing to a formal
  checklist.
- **Schuh, Daniluk & Sieber**, *Chem Sci* 2026 — audit evaluations, report chemically meaningful
  test-set composition, use leakage-resistant splits.
- **Sieg/Flachsenberg/Rarey 2019** — bias-controlled validation setup, with code.
- **Avdiunina/Isayev 2025** — permutation testing as a standard PCM diagnostic.
- **Bai et al. 2021 (BIBM)** — the three pitfalls named: inappropriate splitting, low-confidence
  negatives, drug-wise pair imbalance.


### 3.9 Model cards — the DTI/CPI roster the user named

Compiled from primary sources (papers, GitHub, GitHub API license fields). **No paper in
this list reports a parameter count**, which is itself a small tell about the field's
engineering culture. Estimates below are computed from published layer dimensions and are
marked as estimates.

#### ConPLex — the pure co-embedding pole

| field | value |
|---|---|
| name / year | ConPLex (Contrastive Protein Language-model Extrapolation), 2023, **PNAS** 120(24):e2220778120 |
| URL | https://www.pnas.org/doi/10.1073/pnas.2220778120 · https://github.com/samsledje/ConPLex |
| weights | **YES** — `conplex-dti download --to . --models ConPLex_v1_BindingDB`. One of only three in this roster with downloadable weights. |
| license | **MIT** |
| size | ~3.1 M trainable (est.): `Linear(2048→1024)` drug projector + `Linear(1024→1024)` target projector. ProtBert backbone (~420 M) **frozen**. |
| embeds | protein: frozen ProtBert mean-pooled → 1024-d (code also ships ESM, ProtT5-XL, Bepler–Berger, BindPredict21, FoldSeek featurizers). ligand: **Morgan fingerprint, 2048-bit**. No GNN, no chemical LM. |
| fusion | **Contrastive shared co-embedding.** Two MLP projectors → one shared 1024-d space; score = **cosine similarity** → sigmoid. No cross-attention, no interaction map. Training alternates BCE with a **triplet margin loss** on (protein, true drug, DUD-E decoy), margin annealed 0.25→0. |
| evaluated on | BIOSNAP, BindingDB, DAVIS; DUD-E specificity; prospective kinase panel |
| splits | BIOSNAP: random, **unseen-drug**, **unseen-target**. BindingDB and DAVIS: **random only**. |
| numbers | AUPR — BIOSNAP random 0.897 / unseen-drug 0.874 / **unseen-target 0.842**; BindingDB random 0.628; **DAVIS random 0.458** |

**Honest?** The most trustworthy here, with one caveat. It reports unseen-target and the
number degrades believably; it reports a genuinely bad DAVIS number (0.458) rather than
hiding it; and it ran **prospective wet-lab validation** — 19 kinase–drug predictions, 12
confirmed, one EPHB1 binder at K_D = 1.3 nM. No other model in this roster offers that
category of evidence. Caveat: the *hard* datasets got random splits only, and the
unseen-target split is reported only on the easiest dataset. The DUD-E specificity result
(median Cohen's d 0.730→4.716) is measured against the same decoy generator it was
*trained* on — partly a measurement of "learned DUD-E".

**Why it matters to us.** Cosine similarity between two independently-computed vectors is
O(1) per pair after embedding, so it scales to proteome × library. It also structurally
**cannot** model per-atom contacts. That is the trade.

#### DrugBAN — the most honest evaluation, and the most damning result

| field | value |
|---|---|
| name / year | DrugBAN, 2023, **Nature Machine Intelligence** 5:126–136 |
| URL | https://www.nature.com/articles/s42256-022-00605-1 · https://github.com/peizhenbai/DrugBAN |
| weights | **NO** (training code + data only) · license **MIT** |
| size | not stated; low single-digit millions (est.) |
| embeds | protein: amino-acid embedding → **1D-CNN**, *no PLM at all*. ligand: 2D graph → **GCN**. |
| fusion | **Bilinear attention network** — a bilinear attention map over (drug substructure × protein subsequence), then bilinear pooling → joint vector → FC. Genuinely pairwise-local. Plus **CDAN** conditional domain-adversarial adaptation. |
| splits | random 7:1:2; **cold-pair** on Human; **clustering-based cross-domain** split (single-linkage on ECFP4 for drugs, PSC for proteins, γ=0.5) |
| numbers | in-domain random: BindingDB AUROC **0.960**, BioSNAP **0.903**. **Cross-domain cluster split: 0.575 ± 0.025 (BindingDB), 0.654 ± 0.023 (BioSNAP)**; with CDAN 0.604 / 0.684. |

**Honest?** Yes, conspicuously. It cites the Human dataset's ligand bias and refuses to
rest on the 0.98 random-split number, then publishes the collapse **0.960 → 0.575 AUROC**
— barely above a coin flip. It also states in plain text that **Random Forest consistently
beats DeepConv-DTI, GraphDTA and MolTrans on cross-domain BindingDB**. An author-admitted
result that deep DTI does not generalise out of distribution. A *Nature Machine
Intelligence* Reusability Report (Xu et al. 2024) reproduced and extended it.

**Numbers to believe: the cross-domain ones.** The in-domain +0.008 AUROC over MolTrans is
not a meaningful advance.

#### PerceiverCPI

2023, *Bioinformatics* 39(1):btac731 · https://github.com/dmis-lab/PerceiverCPI · **MIT** ·
**no weights**. Protein: TAPE tokenizer → 1D-CNN with residual blocks and GLU (tokenizer
from TAPE, but **not** a pretrained PLM encoder). Ligand: **dual** — Morgan FP → MLP *and*
D-MPNN graph (chemprop) — fused by single-head cross-attention. **Fusion: nested
cross-attention** (compound = Query, protein = Key/Value; "nested" because cross-attention
is used once inside the compound encoder and once across modalities). Splits are above
average: novel-pair, **novel-compound**, **novel-protein**, and **novel-hard-pair**
(Tanimoto / sequence dissimilarity < 0.3). Evaluated on Davis, KIBA, Metz as **regression**
(MSE, CI) — *do not tabulate these against the classification AUROCs above*. Davis is 100%
dense with only 68 compounds, so "novel-compound" there means holding out a handful of
kinase inhibitors, and the metric is dominated by the mass of pKd≈5 non-binders.
Reproducibility poor.

#### MolTrans

2021, *Bioinformatics* 37(6):830–836 · https://github.com/kexinhuang12345/MolTrans ·
**BSD-3-Clause** · **no weights** (but the processed datasets, including unseen-drug and
unseen-protein variants, *are* shipped — unusually good). ~20 M params (est.), dominated by
embedding tables (23,532 drug + 16,693 protein substructure vocab × 384-d). Both protein
and ligand go through **FCS frequent-substructure mining** → tokens → 2-layer transformer.
**Fusion: outer-product interaction map + CNN** — pairwise dot products between every drug
and protein substructure form a 2-D tensor, then a 3-filter CNN over it.

**Honest?** The **"up to 25% improvement" claim is the least defensible number in this
dossier** — relative improvement, weak baselines (logistic regression, plain DNN, DeepDTI),
random split. DrugBAN later measured MolTrans at 0.895 AUROC vs its own 0.903, i.e. the
whole field sits in a 0.01 band. MolTrans does report unseen-drug/unseen-protein and there
it is only "competitive" — that is the honest signal. Also: the FCS vocabulary is mined over
the union of databases **including test compounds** — mild but real transductive leakage.

#### HyperAttentionDTI

2022, *Bioinformatics* 38(3):655–662 · https://github.com/zhaoqichang/HpyerAttentionDTI
(typo is in the real URL) · **NO LICENSE STATED** (GitHub API returns `license: null` —
default all-rights-reserved; **do not vendor or redistribute**) · **no weights**. Protein and
ligand both from **character** embeddings (64-d) → 3× 1D-CNN. No PLM, no graph, no
fingerprint. **Fusion: "HyperAttention"** — an attention *tensor* over (drug atom × amino
acid × **channel**), so weights are per-feature-channel rather than scalar per pair.
Benchmarks: DrugBank, Davis, KIBA; four settings E1 (both seen) through E4 (both novel),
5-fold CV, random partitioning. DrugBank E1 AUC 0.889.

**Honest?** The negative sampling sinks it. Negatives were drawn by **random sampling of
unlabeled drug–protein pairs** 1:1 against positives — the classic artifact, where the
classifier separates on marginal drug and protein frequency without modelling interaction
at all. Credit for defining E2/E3/E4; **do not take E1 as evidence of anything**, and do not
compare its DrugBank AUC to anyone else's, since that dataset + negative-sampling
combination is unique to this paper.

#### TransformerCPI — the paper that diagnosed the field

2020, *Bioinformatics* 36(16):4406–4414 · https://github.com/lifanchen-simm/transformerCPI ·
**Apache-2.0** · **weights YES** (committed in repo). Protein: overlapping 3-grams →
word2vec (100-d) pretrained on human UniProt → gated CNN. Ligand: RDKit 34-d atom features →
single-layer GCN. **Fusion: transformer decoder cross-attention with the roles inverted** —
the *ligand* attends over the protein, which is what makes atom-level attention
deconvolution interpretable. Human AUC 0.973, C. elegans 0.988, BindingDB 0.951.

**The important part is not the model, it is the diagnostic.** They introduced the
**label-reversal split**: construct test sets where ligands seen *only as negatives* in
training appear *only as positives* at test. A model that memorised "this ligand ⇒ this
label" scores at or below chance. And they report that **CPI-GNN trained on Human achieves
statistically indistinguishable AUC using ligand features alone versus ligand+protein
(p = 0.067)**. The Human dataset's 0.97–0.99 AUCs measure ligand classification, not
interaction. Origin of the artifact: the Human and C. elegans negatives come from Liu et al.
2015, where "highly credible negatives" were produced by **in-silico screening** — the
negative labels are another model's output, drawn from a different chemical distribution
than the positives, never measured.

#### TransformerCPI2.0 — the best-evaluated model in the roster

2023, ***Nature Communications*** 14:4217 ·
https://github.com/lifanchen-simm/transformerCPI2.0 · **GPL-2.0 (copyleft — cannot be linked
into proprietary code without triggering GPL obligations)** · **weights YES** (Google Drive),
but **inference only**, no training code, and inference speed deliberately limited.
~110–120 M params (est.). Protein: **TAPE-BERT** (12 layers, 768-d) — the big upgrade over
v1's word2vec. Ligand: RDKit atom features plus a novel **virtual atom** connected to every
atom, → 1-layer GCN. Fusion: encoder–decoder cross-attention with roles swapped vs v1.
Trained on ChEMBL_23 (3,348 proteins, 69,616 compounds, 117,513 pos / 134,611 neg).

Evaluation: DUD-E and DEKOIS 2.0 enrichment, an external set of 342,447 pairs over **1,192
unseen proteins**, a **time-split** ChEMBL27 set of 92,919 pairs over **637 new proteins**,
and label-reversal. Result vs docking: comparable to GOLD, slightly better than AutoDock
Vina, **worse than Glide SP**.

**Honest?** This is the protocol the rest of the field should be held to, and the conclusion
is appropriately modest: a sequence-only model reaches *mid-tier docking*, not better. One
pushback — DUD-E enrichment factors remain contaminated by analogue bias for *every* method
compared, docking included, so quote the **time-split ChEMBL27** numbers, not the EFs.

#### 2024–2026 successors

- **DTIAM** — *Nat Commun* 16 (2025), doi 10.1038/s41467-025-57828-0. Unified interaction +
  affinity + **mechanism of action (activation vs inhibition)**. Self-supervised pretraining
  of drug and target representations on label-free data, then task heads. Explicitly targets
  **cold-start** and claims its largest gains there. **CC-BY-NC-ND** (no derivatives). The MoA
  task is genuinely novel and much harder to game than binary DTI.
- **LEP-AD** — *J. Cheminform.* 18:101 (**2026**), https://github.com/reem12345/LEP-AD,
  **CC-BY 4.0**. Protein = **ESM-3**; drug = graph via GAT/TransformerConv; **fusion is
  configurable across concatenation / element-wise sum / multi-head cross-attention.**
  *This is the single most useful paper for the fusion question, because it ablates the
  fusion operator directly.* Split: "similar" vs "dissimilar" from Louvain clustering of
  ChemBERTa and ESM-3 embeddings. **Key finding: substantial performance drop in dissimilar
  splits across all methods, including theirs** — i.e. the split dominates the fusion choice.
- **Top-DTI** — *Bioinformatics* 41(Suppl_1):i133–i141 (2025, ISMB),
  https://github.com/bozdaglab/Top_DTI, **CC-BY-NC**. A true **PLM + molecular-LM fusion**:
  ProtT5 (1024-d) + MoLFormer (768-d, pretrained on 1.1 B molecules), plus **persistent
  homology** features (Betti curves, persistence landscapes; 1200-d) from molecular images
  and protein contact maps. **Fusion: a learned sigmoid-weighted element-wise gate** between
  topological and LM features. BioSNAP random AUROC 0.939; **Human cold AUROC 0.898 /
  AUPRC 0.837**. Reporting Human *cold* rather than Human *random* is the right call and
  makes this number comparatively credible.
- **DrugCLIP** — NeurIPS 2023, arXiv:2310.06367. Dense-retrieval reformulation; CLIP-style
  contrastive alignment of 3D **pocket** and molecule towers. Conceptually ConPLex in 3D.
  Zero-shot gains over docking come with the usual DUD-E / LIT-PCBA caveats.
- **ESP-DTI** — ESM-2 + CLIP-style contrastive alignment with progressive contrastive
  curriculum learning, AAAI. ⚠️ **Unverified** — the AAAI and Semantic Scholar records could
  not be loaded (429s); known only second-hand from a 2025 review. Re-check before citing.


---

## Which published results I would and would not believe

**Do not believe a DTI / DTA / PCM result if any of these hold:**

1. **Random / warm split only.** Expect −0.15 to −0.25 AUROC and roughly halved r under a
   cluster or cold-pair split.
2. **Evaluated on Human or C. elegans.** 26 models sit inside ROC-AUC 0.90–0.99 there; the
   benchmark cannot discriminate. Same for DUD-E enrichment with no ligand-only control, and
   LIT-PCBA with no memorisation baseline.
3. **Negatives generated by random pairing**, with no per-drug positive/negative balance
   reported. 91% of BindingDB drugs are single-class.
4. **No protein-blind (ligand-only) ablation.** Cheapest and most decisive control that exists;
   five independent papers show it usually matches the full model.
5. **CASF-2016 / PDBbind numbers with no leakage filter.** 49% of CASF complexes have a near
   neighbour in PDBbind.
6. **A PLM-descriptor claim with no classical-descriptor arm and no leak-proof split.**
7. **A regression result with no median / mean control.**

**Specific numbers I would not repeat:** any Human/C. elegans random-split AUC ≥ 0.95
(including TransformerCPI's own 0.973 / 0.988); MolTrans's "up to 25% improvement";
HyperAttentionDTI's DrugBank E1 figures; DUD-E enrichment factors for *any* method, ConPLex's
Cohen's d included; the entire 0.89–0.96 AUROC band on BindingDB/BioSNAP random splits (six
methods, 0.01 spread, Random Forest inside it — the ranking is noise); Davis CI/MSE gains
(68 compounds, 100% dense, dominated by the pKd≈5 mass).

**Numbers I would believe, and would treat as the real state of the art:**

- DrugBAN cross-domain cluster split: **AUROC 0.575 (BindingDB) / 0.654 (BioSNAP)**;
  0.604 / 0.684 with domain adaptation. This is what "novel chemotype × novel target family"
  actually costs.
- DrugBAN's own admission that **RF beats DeepConv-DTI, GraphDTA and MolTrans cross-domain**.
- ConPLex BIOSNAP unseen-target AUPR **0.842**, and DAVIS random AUPR **0.458**.
- TransformerCPI2.0 on a **time split with 637 unseen proteins**, reaching docking-comparable
  but **sub-Glide** performance.
- Papyrus PCM **r 0.79 → 0.42** random → temporal.
- Graber CleanSplit: Pafnucy **0.808 → 0.746 Pearson**, RMSE **1.046 → 1.484**.
- LEP-AD (2026): "substantial performance drop in dissimilar splits **across all methods**."

**On the fusion question the user asked:** the five distinct mechanisms in the literature are
contrastive shared co-embedding with a dot product (ConPLex, DrugCLIP, ESP-DTI); bilinear
attention + bilinear pooling (DrugBAN); outer-product interaction map + CNN (MolTrans);
transformer cross-attention (TransformerCPI/2.0, PerceiverCPI, LEP-AD); and learned gated
element-wise mixing (Top-DTI). **No published evidence establishes that any of these beats the
others once the split is honest.** LEP-AD is the only paper that ablates the fusion operator
directly, and its finding is that **the split dominates the fusion choice.** The practical
differentiator is cost: co-embedding is the only mechanism that amortises to one dot product
per pair, and therefore the only one that runs at proteome × library scale — and the only one
that supports the "what should bind here?" query `README.md` wants.

---

## What this means for the world model — the part to act on

1. **Kill criterion 1 is the likely outcome, not a formality.** Two leak-controlled studies
   (Isayev's permutation test; HonestAffinity's no-leak tiers) found PLM protein embeddings
   contribute ~nothing to PCM, and that the *simpler* representation can win once leakage is
   removed. Plan for having to earn the embedding.
2. **Adopt the Papyrus + QSPRpred stack for the baseline** rather than rolling our own splits.
   The split bugs are the expensive kind and they are already solved there.
3. **Pre-register the controls before the model.** Minimum honest protocol, and it should go
   into the repo's rules alongside "measure the pool oracle first":
   - cluster split (single-linkage on ECFP4 for ligands, MMseqs2 sequence-identity for
     proteins) **+** cold-pair **+** temporal where dates exist;
   - a **degree/frequency null model** (AI-Bind's configuration model);
   - a **ligand-only model** with the protein input ablated;
   - a **median baseline** for any regression;
   - a **permutation test on the protein descriptor**;
   - report **sensitivity separately** — AUROC hides the collapse to "predict negative";
   - report the **gap to the pool oracle**, per this repo's existing rule.
4. **This repo's own noise floor already encodes the lesson.** FINDING 007 put it at +0.0138
   and the memory note says anything under +0.020 is noise. The DTI literature's entire
   0.89–0.96 AUROC band is the same phenomenon at a different scale. Measure the null first.
5. **ATOMICA is the near-term actionable item**, as a pose scorer for the CYP3A4 pool — with
   the rotated-pose ablation pre-registered, because its pretraining objective *is* rotation
   denoising.
6. **The world-model gap is real.** No bio-JEPA predicts a binding partner's embedding from a
   target's. If we build it, we are first — and the honest reading of Part 3 is that the reason
   may be that nobody has a benchmark that could tell whether it worked.

---

## BLOCKED — needs the user

1. **UniBioseq — I cannot find it.** Zero hits on arXiv (`all:UniBioseq`) and zero on Europe
   PMC. Need the exact name, a link, or where it was seen. Do not plan around it until then.
2. **GenoJEPA has no verifiable primary source.** Reported (52M params beating NT-v2's 494M
   across 55 tasks, built on LeJEPA) only via a HuggingFace collection and a Substack. Not on
   arXiv. If the user has the paper, it changes the parameter-efficiency argument materially.
3. **Web search budget was exhausted mid-run** (200/200 calls, shared across this session's
   agents). The last third of all three parts was done by direct fetch and API queries.
   Coverage is good but not exhaustive; raising `CLAUDE_CODE_MAX_WEB_SEARCHES_PER_SESSION`
   would let a follow-up close the gaps below.
4. **Two paywalled papers I could only read at abstract level**, both load-bearing for Part 3:
   - **Avdiunina, Jamal, Gusev, Isayev**, *JCIM* 2025, doi 10.1021/acs.jcim.5c00395 — the
     single most on-point PCM-evaluation paper. Does the user have institutional access
     (Northeastern)?
   - **Schuh, Daniluk & Sieber**, *Chem Sci* 2026, doi 10.1039/d6sc01799a — the concrete
     contamination percentages could not be extracted.
   - Also **DrugBAN's Reusability Report** (Xu et al., *Nat Mach Intell* 2024,
     doi 10.1038/s42256-024-00822-w) — paywalled, no abstract in the Semantic Scholar record.
5. **ESP-DTI is unverified.** AAAI and Semantic Scholar records returned 429s; known only
   second-hand from a 2025 review. Re-check before citing.
6. **Licensing decisions the user must make before anything is vendored:**
   - **HyperAttentionDTI has no license at all** — default all-rights-reserved. Do not copy.
   - **TransformerCPI2.0 is GPL-2.0** — copyleft. Cannot be linked into anything proprietary.
   - **Top-DTI (CC-BY-NC)**, **DTIAM (CC-BY-NC-ND)**, **ProtJEPA (CC-BY-NC)** — non-commercial.
   - **GeneJEPA states no license** in its repo despite shipping weights.
   - Is this project's output commercial or academic? The answer changes which half of this
     dossier is usable.
7. **Decision needed on scope.** Part 3 says the honest evaluation of an interaction model is
   more expensive than the model. With the final deadline on 2026-11-03 and one contended
   H200, we can afford *either* a proper leak-controlled PCM baseline *or* a first JEPA
   prototype — probably not both. Which does the user want first? My recommendation is the
   baseline, because without it no JEPA number will be interpretable.
8. **ATOMICA pose-scoring is ready to start and needs no new data** (checkpoints are CC-BY-4.0
   on HuggingFace, code MIT, and the existing CYP3A4 pool is already on disk as PDBs). Say the
   word and it runs against `cypstruct.xengine.select()`'s +0.0381 incumbent.
