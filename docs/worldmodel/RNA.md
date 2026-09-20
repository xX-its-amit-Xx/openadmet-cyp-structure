# RNA — models and data for the shared embedding space

**Author:** RNA lead, wave 1. **Written 2026-09-20.** Sections 1, 2, 4, 5, 7 and 8 are this
lead's own work; section 3 and most of section 6 were gathered by parallel reconnaissance agents
and merged here. Every PDB count was re-derived live from the RCSB Search API on 2026-09-20 and,
where two sweeps used different drug-likeness thresholds, **both** numbers are reported.

Hard constraints this document is written against: **one contended H200 (80 GB)**,
**298 TB cluster scratch**, **no local disk**, **torch 2.5.1**, **no custom CUDA kernels**
(no flash-attn build, no triton kernels, no mamba-ssm, no transformer-engine, no cuEquivariance).
That last constraint eliminates more RNA models than any scientific consideration does, so it
is flagged in every row.

---

## 0. The one-paragraph answer

For a frozen embedding source, three models carry the field: **RiNALMo-giga** (650 M, the
largest honestly-open RNA LM), **AIDO.RNA-1.6B** (the only >1 B RNA LM with downloadable
weights), and **RNA-FM** (99 M, the most-used, MIT, and the backbone of the RhoFold structure
stack so its embedding space already has a 3D decoder attached). All three run in fp16 on one
H200 with room to spare — RNA LMs are *small*; the binding constraint is sequence length, not
parameters. For data, RNA has the same shape of problem as protein did in 2018: sequence is
abundant (**58.6 M RNAcentral sequences**), structure is scarce (**10,399 PDB entries with RNA,
2,375 RNA-only**), and **RNA–small-molecule affinity is essentially absent**. See §6 for the
honest tally; the short version is **101** RNA–ligand pairs with both a structure and a measured
affinity, **~2,500** RNA–ligand affinity measurements in aggregated form worldwide, and **48**
non-redundant drug-like RNA–ligand co-structures ever solved — against PDBbind's ~20,000 and
BindingDB's 3.2 M. The one place RNA is *data-rich* is chemical probing: **Ribonanza's ~2 M
experimentally probed sequences** have no protein analogue, and only one small model (11 M
parameters) has ever been trained on them.

---

## 1. RNA foundation models (sequence LMs) — the table

Hidden dim = the per-nucleotide embedding you would actually extract. "Size" = weights on disk
as distributed (mostly fp32; halve for fp16). All sizes verified from the HuggingFace API on
2026-09-20 unless marked *(est.)*.

| Model | Ver | Yr | URL | Weights open? | License | Params | Emb dim | Max ctx (nt) | Input | Size GB | 1× H200? | Exotic deps | Use in shared space |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **RNA-FM** | v1 | 2022/24 | [ml4bio/RNA-FM](https://github.com/ml4bio/RNA-FM) · [multimolecule/rnafm](https://huggingface.co/multimolecule/rnafm) | yes | **MIT** (orig) / AGPL-3.0 (MM rewrap) | 99.5 M | **640** | 1024 | nt string | **0.40** | trivially | none — pure PyTorch | Default RNA encoder; its space already feeds RhoFold/RhoDesign, so 3D decoders exist |
| **mRNA-FM** | v1 | 2024 | [multimolecule/mrnafm](https://huggingface.co/multimolecule/mrnafm) | yes | MIT/AGPL | 239 M | 1280 | 1024 codons | codons (len%3=0) | 0.96 *(est.)* | yes | none | Coding-sequence arm; pairs with CodonBERT for CDS tasks |
| **RiNALMo-micro** | 1.0 | 2024/25 | [lbcb-sci/RiNALMo](https://github.com/lbcb-sci/RiNALMo) | yes | Apache-2.0 code / **CC-BY-4.0 weights** | 33 M | 480 | 1024 | nt | 0.13 | yes | *orig repo pins flash-attn 2.3.2* | Cheap ablation control |
| **RiNALMo-mega** | 1.0 | 2024/25 | [multimolecule/rinalmo-mega](https://huggingface.co/multimolecule/rinalmo-mega) | yes | CC-BY-4.0 / AGPL | 148 M | 640 | 1024 | nt | 0.59 | yes | MM rewrap avoids flash-attn | Mid-size point on the scaling curve |
| **RiNALMo-giga** | 1.0 | Nat.Comm. 2025 | [multimolecule/rinalmo-giga](https://huggingface.co/multimolecule/rinalmo-giga) · [Zenodo 15043668](https://zenodo.org/records/15043668) | yes | CC-BY-4.0 / AGPL | **650.9 M** | **1280** | 1024 | nt | **2.60** | yes (~3 GB fp16) | **orig repo requires `flash-attn==2.3.2` + CUDA≥11.8 → use the MultiMolecule rewrap instead, which is plain HF attention** | **Primary RNA encoder.** Best or near-best on every published benchmark; 1280-d matches ESM-2-650M so a shared projection head is symmetric |
| **AIDO.RNA-650M** | 1.0 | 2024 | [genbio-ai/AIDO.RNA-650M](https://huggingface.co/genbio-ai/AIDO.RNA-650M) | yes | "other" (ModelScope-style, **read before use**) | 648 M | 1280 | 1024 | nt | ~5.2 *(est.)* | yes | `modelgenerator` package | Second opinion at matched size |
| **AIDO.RNA-1.6B** | 1.0 | 2024 | [genbio-ai/AIDO.RNA-1.6B](https://huggingface.co/genbio-ai/AIDO.RNA-1.6B) | yes, not gated | **"other"** | 1.6 B | **2048** | 1024 | nt | **14.36** | yes (~3.2 GB fp16 weights, fits easily) | `modelgenerator`/`bionemo`-ish stack | Largest downloadable RNA LM; use as the capacity ceiling probe |
| **ERNIE-RNA** | — | Nat.Comm. 2025 | [multimolecule/ernierna](https://huggingface.co/multimolecule/ernierna) | yes | AGPL (MM) | 86 M | 768 | 1024 | nt | 0.35 *(est.)* | yes | none | Its attention maps are base-pair-like — gives a *structure-aware* embedding at 86 M |
| **ERNIE-RNA-ss** | — | 2025 | [multimolecule/ernierna-ss](https://huggingface.co/multimolecule/ernierna-ss) | yes | AGPL | 86 M | 768 | 1024 | nt | 0.35 *(est.)* | yes | none | Free secondary-structure channel to concatenate |
| **RNAErnie** | — | Nat.Mach.Int. 2024 | [multimolecule/rnaernie](https://huggingface.co/multimolecule/rnaernie) · [CatIIIIIIII/RNAErnie](https://github.com/CatIIIIIIII/RNAErnie) | yes | AGPL (MM); orig is PaddlePaddle | 86 M | 768 | **512** | nt | 0.35 *(est.)* | yes | **original is PaddlePaddle, not torch** — use the MM port | Motif-aware pretraining; weakest context window of the majors |
| **RNA-MSM** | — | NAR 2024 | [yikunpku/RNA-MSM](https://github.com/yikunpku/RNA-MSM) · [multimolecule/rnamsm](https://huggingface.co/multimolecule/rnamsm) | yes | **MIT** (orig) / AGPL (MM) | **95.9 M** | 768 | 1024 | **MSA** | **0.38** | yes | ⚠️ **needs an MSA from RNAcmap3 (blastn vs `nt` + Infernal)** — GPU nodes have no internet; MSAs must be precomputed on a login node | Only evolutionary-coupling RNA model; best SS scores in two independent benchmarks. Treat the MSA pipeline as its own multi-day task |
| **RNABERT** | — | NARGAB 2022 | [multimolecule/rnabert](https://huggingface.co/multimolecule/rnabert) | yes | AGPL | 0.5 M | 120 | 440 | nt | <0.01 | yes | none | Baseline / sanity floor only. **Tiny.** |
| **UTR-LM (te_el, mrl)** | — | Nat.Mach.Int. 2024 | [multimolecule/utrlm-te_el](https://huggingface.co/multimolecule/utrlm-te_el) | yes | AGPL-3.0+ | **1.21 M** | **128** | 1022 | 5′UTR nt | <0.01 | yes | none | ⚠️ **A 1.2 M-parameter, 6-layer, 128-d model.** Do not call this a foundation model. Useful only as a 5′UTR-specific feature |
| **3UTRBERT** | 3/4/5/6-mer | 2024 | [multimolecule/utrbert-3mer](https://huggingface.co/multimolecule/utrbert-3mer) | yes | AGPL | ~86 M | 768 | 512 | 3′UTR k-mers | 0.35 *(est.)* | yes | none | 3′UTR arm; k-mer tokenisation makes alignment to nt-level spaces awkward |
| **SpliceBERT** | 510 / full | 2024 | [multimolecule/splicebert](https://huggingface.co/multimolecule/splicebert) | yes | AGPL | 19.4 M | 512 | 1024 | pre-mRNA nt | 0.08 | yes | none | Splice-context features; strongest zero-shot fitness model in the 2026 benchmark |
| **CodonBERT** | — | Genome Res. 2024 | [Sanofi-Public/CodonBERT](https://github.com/Sanofi-Public/CodonBERT) | **weights behind a Sanofi CDN link with a separate, non-standard model licence** | code: BSD-ish; **weights: bespoke** | ~87 M *(unverified — paper paywalled)* | 768 *(unverified)* | ~1024 codons | codons | ~0.35 *(est.)* | yes | poetry env | mRNA/CDS embeddings. ⚠️ licence and paper both gated — see BLOCKED |
| **CaLM** | — | Nat.Mach.Int. 2024 | [oxpig/CaLM](https://github.com/oxpig/CaLM) · [multimolecule/calm](https://huggingface.co/multimolecule/calm) | yes | AGPL (MM) | 86 M | 768 | 1024 codons | codons | 0.35 *(est.)* | yes | none | Codon LM trained on coding sequences — the natural bridge from RNA space to *protein* space |
| **ProtRNA** | — | Cell Syst. 2025 | [roxie-zhang/ProtRNA](https://github.com/roxie-zhang/ProtRNA) · [Zenodo 14795554](https://zenodo.org/records/14795554) | yes | **Apache-2.0** (verified) | 651 M | 1280 | **512** | nt | ~2.6 *(est.)* | yes | ⚠️ **`tensorflow==2.14.0`, not PyTorch** — needs a port or a separate env | ⭐ **Transfer-learned from ESM-2**: a protein LM re-tuned on RNA with 1/8 the trainable params and 1/6 the data of native RNA LMs. **The existing prior art for the protein↔RNA shared space — read it before designing ours** |
| **RibonanzaNet** | v1 | 2024 | [multimolecule/ribonanzanet](https://huggingface.co/multimolecule/ribonanzanet) · [Shujun-He/RibonanzaNet](https://github.com/Shujun-He/RibonanzaNet) | yes | AGPL (MM) | ~11 M | 256 | 457 (train) | nt | 0.05 | yes | none | **The only RNA model trained on millions of real experimental structure measurements.** Small, but its features are grounded in chemistry, not co-evolution |
| **RibonanzaNet-SS / -DEG / -DROP** | — | 2024 | [multimolecule/ribonanzanet-ss](https://huggingface.co/multimolecule/ribonanzanet-ss) | yes | AGPL | ~11 M | 256 | — | nt | 0.05 | yes | none | Ready-made SS / degradation heads |
| **Orthrus (large-6-track)** | — | Nat.Methods 2026 | [bowang-lab/Orthrus](https://github.com/bowang-lab/Orthrus) · [quietflamingo/orthrus-large-6-track](https://huggingface.co/quietflamingo/orthrus-large-6-track) | yes | **MIT** | **~10 M** (40.7 MB ckpt) | 256–512 | very long (SSM) | mature mRNA, 4- or 6-track one-hot | **0.04** | yes | ⚠️ **Mamba backbone — `mamba-ssm` ships custom CUDA kernels.** A slow pure-PyTorch fallback exists but is not the shipped path | Best *mature-mRNA property* embeddings per parameter; contrastive objective over isoforms + 400-species orthologs. ⚠️ **10 M params — "foundation model" is generous** |
| **RIBOSPAN-FM** | 1K-15 / 1K-40 public; **10K-15 on request** | 2026 | [SII-GAIR-NLP/RIBOSPAN-FM](https://huggingface.co/SII-GAIR-NLP/RIBOSPAN-FM) | **partly gated** | **non-commercial bespoke** | 1.61 B | **2048** | **10 240** (long-ctx ckpt) | nt | ~6.5 *(est., fp32 ≈ 6.4 GB)* | yes | — | Only RNA LM with a 10 k-nt context — the one that could embed a full mRNA. ⚠️ long-context weights **gated**, licence non-commercial |
| **NV-CodonFM (Encodon-1B)** | v1 | 2025 | [nvidia/NV-CodonFM-Encodon-1B-v1](https://huggingface.co/nvidia/NV-CodonFM-Encodon-1B-v1) | yes | **NVIDIA Open Model License** | 911 M | ~1536 *(unverified)* | **2046 codons** | codons | ~3.6 *(est.)* | yes — **NVIDIA explicitly states torch 2.5.1, Hopper** | BioNeMo-adjacent stack | Largest codon-level model; RefSeq CDS from 22 k species. Licence is permissive-ish but not OSI |
| **Evo 2 (7B / 40B)** | 2 | 2025/26 | [ArcInstitute/evo2](https://github.com/ArcInstitute/evo2) · [arcinstitute/evo2_7b](https://huggingface.co/arcinstitute/evo2_7b) | yes | **Apache-2.0** | 7 B / 40 B | 4096 | 1 M (7B) | DNA/RNA nt | **13.77** (7B) | 7B yes; **40B no** on one H200 | ⚠️ **StripedHyena-2; needs `vortex` + transformer-engine + FP8 kernels.** Violates the no-custom-kernels rule | Genome-scale context; would be the only model able to see an entire transcript's genomic neighbourhood. Park it |
| **LucaOne** | default | Nat.Mach.Int. 2025 | [LucaOne/LucaOneApp](https://github.com/LucaOne/LucaOneApp) | yes | check | ~1.8 B | 2560 | 1280 | **DNA+RNA+protein in one tokenizer** | ~7 *(est.)* | yes | custom repo, not HF-native | ⭐ **Already a shared nucleic-acid + protein space.** The single most relevant prior art for ask 1 — test it before building our own |
| **OmniGenome-186M** | — | AAAI 2025 | [yangheng/OmniGenome-186M](https://huggingface.co/yangheng/OmniGenome-186M) | yes | check | 186 M | 768 | 1024 | nt (+SS) | 0.75 *(est.)* | yes | none | Sequence↔structure aligned pretraining |
| **PlantRNA-FM** | — | Nat.Mach.Int. 2024 | [yangheng/PlantRNA-FM](https://huggingface.co/yangheng/PlantRNA-FM) | yes | check | 35 M | 480 | 512 | nt | 0.14 | yes | none | Plant-specific; out of scope but a clean domain-shift control |
| **BiRNA-BERT** | — | Comm.Bio 2025 | [buetnlpbio/BiRNA-BERT](https://github.com/buetnlpbio/BiRNA-BERT) | yes | check | 116 M | 768 | 1024 | nt **and** BPE | 0.47 *(est.)* | yes | none | Dual tokenisation (nt + BPE) — lets one model serve both fine and long-range views |
| **RNA-km** | — | 2024 | [gongtiansu/RNA-km](https://github.com/gongtiansu/RNA-km) | yes | check | 152 M | 1024 | 512 | nt | 0.61 *(est.)* | yes | none | — |
| **LAMAR (2k/4k)** | — | 2025 | search "LAMAR RNA language model" | yes | check | 86 M | 768 | **2048 / 4096** | nt | 0.35 *(est.)* | yes | none | Cheapest route to >1 k context |
| **RFamLlama (base/large)** | — | 2024 | [jinyuan22/RFamLlama-base](https://huggingface.co/jinyuan22/RFamLlama-base) | yes | check | 31 / 89 M | 512 / 768 | — | nt (causal) | 0.12 / 0.36 | yes | none | ⚠️ trained on **676 k Rfam seqs only** — family-memorising, not general |
| **Uni-RNA (L8/L12/L16)** | — | 2023 | [ComDec/unirna_tf](https://github.com/ComDec/unirna_tf) | ⚠️ **weights never released by DP Technology** | — | 21/71/168 M | 512/768/1280 | 1024 | nt | — | — | — | **Cite, do not plan around.** Frequently benchmarked from a private checkpoint |
| **ATOM-1** | — | 2023 | [preprint](https://doi.org/10.1101/2023.12.13.571579) | ❌ **closed, Atomic AI commercial** | — | undisclosed | — | — | — | — | — | — | Trained on proprietary chemical mapping. **Not obtainable.** |
| **DGRNA / HydraRNA / CodonMamba** | — | 2024–26 | various | yes | check | 108 M+ | 768 | 2048+ | nt | ~0.4 | yes | ⚠️ **Mamba / Hydra state-space — custom CUDA kernels** | Long-context alternatives, all blocked by the kernel constraint |
| **mRNABERT** | — | Nat.Comm. 2025 | [Taykhoom/mRNABERT-no-flashattention](https://huggingface.co/Taykhoom/mRNABERT-no-flashattention) | yes | check | — | — | — | mRNA nt | — | yes | the author ships an explicit **no-flash-attention** checkpoint — use it | Convenient given our constraint |
| **NucleicBERT** | — | Nat.Mach.Int. 2026 | [KIT-MBS/NucleicBERT](https://github.com/KIT-MBS/NucleicBERT) | yes | check | — | — | — | DNA+RNA | — | yes | — | Recent; unreplicated. Read before trusting |
| **METAGENE-1** | — | 2025 | [metagene-ai/METAGENE-1](https://huggingface.co/metagene-ai/METAGENE-1) | yes | Apache-2.0 | 7 B | 4096 | 512 | metagenomic nt | ~14 | yes | none | Wrong domain (metagenomics reads), listed for completeness |

### Notes on §1

1. **MultiMolecule is the single biggest practical win and also the biggest licensing trap.**
   [`multimolecule`](https://multimolecule.danling.org) re-implements ~25 RNA models on plain
   HuggingFace `transformers` with a uniform API, identical tokenizer conventions and **no
   flash-attn, no Paddle, no TensorFlow**. It is how we get RiNALMo without building flash-attn
   and RNAErnie without PaddlePaddle. **But every MultiMolecule checkpoint is redistributed
   under AGPL-3.0**, which is *more restrictive than the originals* (RNA-FM is MIT, RiNALMo
   weights are CC-BY-4.0, RNA-MSM is MIT). For internal research AGPL is fine; if anything is
   ever served, pull the original weights from Zenodo/GitHub and load them with the MM code
   under its licence terms, or re-derive. **Flag this to legal before any deployment.**
2. **RNA LMs are small.** The entire §1 table, every checkpoint at once, is <60 GB — less than
   one H200's memory budget for a single 7 B protein model. Compute is not the constraint;
   *sequence length* is. Almost every model caps at 1024 nt, which is fine for ncRNA and
   riboswitches and useless for mRNA. Only LAMAR (4 k), NV-CodonFM (2046 codons = 6138 nt) and
   RIBOSPAN (10 k) exceed it, and the best of those is gated.
3. **The field has not converged and the benchmarks disagree.** The 2026 zero-shot benchmark
   (21 models, [Brief. Bioinform.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12963973/)) found
   RNA-MSM best at secondary structure, RNA-FM best at family classification, SpliceBERT/MP-RNA
   best at fitness — **and median zero-shot fitness correlation under 0.27 across every model**.
   The [2025 SS benchmark](https://academic.oup.com/bib/article/26/2/bbaf137/8109668) found
   ERNIE-RNA/RiNALMo/RNA-MSM best. **No model wins everything, and every model is weak in
   absolute terms.** Plan on an ensemble of 2–3, not a champion.
4. **Skeptical flags.** (a) *Uni-RNA* is benchmarked everywhere and its weights were never
   released — treat every Uni-RNA number as unverifiable. (b) *UTR-LM* (1.2 M params) and
   *Orthrus* (10 M) are published in Nature-family journals as "foundation models"; they are
   small special-purpose networks. (c) *RFamLlama* trained on 676 k Rfam sequences will look
   excellent on Rfam-derived benchmarks for the wrong reason. (d) *ATOM-1* is cited as a
   chemical-mapping foundation model and is commercially closed. (e) Most 2026 entries have a
   single paper and no independent replication.

---

## 2. Which ones we would actually run

| Rank | Model | Why | Cost on one H200 |
|---|---|---|---|
| 1 | **RiNALMo-giga via MultiMolecule** | Best average across independent benchmarks; 1280-d matches ESM-2-650M exactly, so the RNA↔protein projection head is dimensionally symmetric; no flash-attn in the MM path | 2.6 GB weights; ~58 M RNAcentral seqs at 1024 nt ≈ 2–3 GPU-days at batch 64 fp16. Embed a **subset**, not all of RNAcentral |
| 2 | **RNA-FM** | MIT, 0.4 GB, and — uniquely — **a trained 3D decoder already reads this space** (RhoFold+, RhoDesign, RiboDiffusion). If the shared space must ever emit structure, start here | hours |
| 3 | **AIDO.RNA-1.6B** | Capacity ceiling probe; if 1.6 B does not beat 650 M, stop scaling RNA | 14.4 GB fp32 → 3.2 GB fp16 |
| 4 | **RibonanzaNet** | Grounded in 2 M *measured* chemical-mapping profiles rather than co-evolution. Orthogonal error mode to every LM in the table | minutes |
| 5 | **LucaOne** | Already a joint DNA+RNA+protein space. **Test it as the null hypothesis for ask 1** — if an existing multi-modal model already does what our JEPA is meant to do, that is the most important thing we could learn in week 1 | ~7 GB |
| — | **RNA-MSM** | Only if someone owns the MSA pipeline end-to-end. Per this repo's PXR lesson, MSA staging on an internet-less GPU partition is its own multi-day task and has already cost one campaign 12 days | + days of CPU |

**Kill criterion 1 applies directly:** if concatenating a frozen RiNALMo embedding to a Morgan
fingerprint does not beat one-hot-target PCM on a held-out RNA–ligand affinity task, the RNA arm
of the shared space is not carrying biology. Given §4, that test may not have enough data to run
at all, which is itself the finding.

---

## 3. RNA structure models that expose embeddings

The extractable representation is the point of this table, not the accuracy.

| Model | Ver / yr | URL | Weights open? License | Params | **Extractable emb dims** | Max ctx | Input | DL size GB | 1×H200 / torch 2.5.1 / no kernels | MSA? | Use in shared space |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **RhoFold+** | v1.0, Nat.Meth. 2024 | [ml4bio/RhoFold](https://github.com/ml4bio/RhoFold) · [cuhkaih/rhofold](https://huggingface.co/cuhkaih/rhofold) | yes, **Apache-2.0**, ungated | ~127 M (RNA-FM 99.5 M + Rhoformer) | **c_s=384 single, c_z=128 pair, c_m=256 MSA** (verified in `rhofold/config.py`) | LM branch 1024 | FASTA (+opt. a3m) | **0.51** | yes, pure torch | optional; full accuracy wants ~900 GB DBs | **Best single RNA structure target.** Its 384/128 dims are *identical* to Boltz-2's trunk — one adapter serves both |
| **DRfold2** | v2, PLOS Biol. 2026 | [leeyang/DRfold2](https://github.com/leeyang/DRfold2) | yes, **MIT** | RCLM 47.5 M + trunk | **RCLM seq 512 / pair 128** (18 blocks); trunk 64/64 | crop 256–384, no hard cap | FASTA only | ~1.3 | yes (tested 1.11/2.0.1/2.2.1) | **MSA-FREE** | The only MSA-free end-to-end 3D model. Default choice for an airgapped GPU node |
| **trRosettaRNA2** | v2.0.2, 2025/26 | [YangLab-SDU/trRosettaRNA2](https://github.com/YangLab-SDU/trRosettaRNA2) · [weights](https://huggingface.co/datasets/quailwwk/trRNA2) | weights downloadable; **no LICENSE file in repo** | not disclosed | **unverified** (bioRxiv 429'd; read the checkpoint) | — | a3m/FASTA/SS | — | yes, torch 2.0+, PyG | explicit single-seq mode; DB only 32 GB — cheapest MSA path | Generates *conformers*, rare and valuable for an ensemble-aware model. Best automated CASP16 server (Yang-Server) |
| **RNAPro** | v1, **2026-01-09** | [NVIDIA-Digital-Bio/RNAPro](https://github.com/NVIDIA-Digital-Bio/RNAPro) · [nvidia/RNAPro-Public-Best-500M](https://huggingface.co/nvidia/RNAPro-Public-Best-500M) | ungated; code Apache-2.0, **weights NVIDIA Open Model License** | **488,301,921** | RibonanzaNet2 frozen encoder -> learned gating -> Protenix Pairformer; per-block dims unpublished | **crops to 512**, skips >10 000 | CSV + MSA FASTA | ~2 *(est.)* | yes; Triton/cuEq/TE optional | user-supplied MSAs, no bundled server, stage-friendly | **Someone already built our architecture.** NVIDIA x Das Lab x Kaggle winners: a frozen RNA foundation model fused into a co-folder by learned gating. Read its gating module before designing ours. Caveat: "best" = **Kaggle private leaderboard**, not an independent benchmark |
| **Boltz-2** | 2025-06 | [jwohlwend/boltz](https://github.com/jwohlwend/boltz) | yes, **MIT code + MIT weights** | unpublished | **64-layer Pairformer, single 384, pair 128** | affinity head caps ligands at 128 atoms | YAML, `entity_type: rna` | auto | yes (cuEquivariance optional) | optional | **CORRECTION, verified in source 2026-09-20:** my structure agent reported "no `--write_embeddings` flag". It **is shipped** — declared at `src/boltz/main.py:1038`, written by `src/boltz/data/write/writer.py:250` as `embeddings_<id>.npz` holding `s [n,384]` and `z [n,n,128]`; and `src/boltz/data/parse/schema.py:1021` accepts `entity_type: rna`. **RNA + ligand + embeddings works today with no patch.** The deeper per-layer recipe (concat layers 16/32/48/64) from [hsjang0/boltz-as-FM](https://github.com/hsjang0/boltz-as-FM) (MIT) still needs a patch. **RNA transfer untested — that paper is protein/ligand only** |
| **Boltz-1 / 1x** | 2024/25 | same | MIT | — | single 384, pair 128 | — | YAML | auto | yes | `msa: empty` | Same trunk dims; 1x adds inference-time potentials |
| **Protenix / Protenix-v2** | v1.0.0 2026-02; **v2 464 M 2026-04-08** | [bytedance/Protenix](https://github.com/bytedance/Protenix) | yes, **Apache-2.0 code AND weights** | 368 M / **464 M** | AF3-style, dims not in README | — | JSON/FASTA | pip | yes (custom LayerNorm optional, torch fallback) | **RNA MSA support only from v1.0.0** | Apache-2.0 weights make it a first-class citizen. Note this repo's FINDING 009: Protenix is **deterministic** — samples do not diversify |
| **OpenFold3** | preview, 2026 | [aqlaboratory/openfold-3](https://github.com/aqlaboratory/openfold-3) · [openfold/openfold3](https://huggingface.co/openfold/openfold3) | Apache-2.0 **but the HF repo is click-through GATED** | not disclosed | AF3-style | low-mem mode | AF3 JSON | — | works without kernels; cuEq + DeepSpeed4Science are the documented path | jackhmmer/hhblits or ColabFold | Best-licensed AF3 reproduction with real RNA training. "Only model to match AF3 on monomeric RNA" is a **vendor claim on preview weights** |
| **AlphaFold3** | 3.0.x | [google-deepmind/alphafold3](https://github.com/google-deepmind/alphafold3) | **GATED.** Params by Google form, **non-commercial, no redistribution, no derived-weight sharing** | — | trunk 384/128 (architecture-standard, unconfirmed) | ~5000 tokens | JSON | params ~1 GB; DBs ~630 GB | **no — JAX + Triton kernels** | jackhmmer/nhmmer | **Do not build on this.** The licence forbids exactly the redistribution a shared model implies. Accuracy reference only |
| **Chai-1** | 0.6.1 | [chaidiscovery/chai-lab](https://github.com/chaidiscovery/chai-lab) | HF card now **apache-2.0** — **changed from the original non-commercial release**; stale pages still say otherwise. Verify the file you pull | unpublished | undocumented | — | FASTA/API, RNA supported | — | yes, Linux + bf16 | optional MMseqs2 | Patch the trunk as for Boltz |
| **Chai-2** | 2025 report | chaidiscovery.com | **NO WEIGHTS**, API/partnership only | — | — | — | — | — | — | — | Not usable. Do not plan around it |
| **NuFold** | Nat.Comm. 2025 | [kiharalab/NuFold](https://github.com/kiharalab/NuFold) | weights open, **GPL-3.0 — strong copyleft** | not disclosed | not disclosed | — | FASTA + a3m + ipknot SS | — | yes (deepspeed listed but usually optional) | **rMSA + ~2 TB DB — effectively disqualifying** | Nucleobase-centred frames = a genuinely different geometric bias; scored *best* on the NAR "easy" benchmark (TM 0.482). **GPL + 2 TB kills it for us** |
| **RoseTTAFold2NA** | v0.2, **Apr 2023 — stale** | [uw-ipd/RoseTTAFold2NA](https://github.com/uw-ipd/RoseTTAFold2NA) | 1.1 GB weights; code MIT, IPD weights conventionally non-commercial | not disclosed | 3-track 1D/2D/3D | — | FASTA + MSA | 1.1 + **~480 GB DBs** (nt alone 151 GB) | **no — SE3Transformer + DGL**; DGL wheels lag torch 2.5.1 | required | Stale, 74 open issues. Backbone inside RNAFlow. Completeness only |
| **RoseTTAFold-All-Atom** | 2024 | [baker-laboratory/RoseTTAFold-All-Atom](https://github.com/baker-laboratory/RoseTTAFold-All-Atom) | weights open; template DBs restricted | — | — | — | FASTA/SMILES | — | no — SE3Transformer/DGL | **README says it CANNOT build RNA MSAs** | RNA is second-class here. Deprioritise |
| **HelixFold3** | v3.2, 2025-07 | [PaddleHelix](https://github.com/PaddlePaddle/PaddleHelix/tree/dev/apps/protein_folding/helixfold3) | weights **non-commercial only** | — | — | ~1200 tokens on A100-40G | AF3-style | DBs 190 -> 530 GB | **no — PaddlePaddle 3.1.0, not PyTorch** | required | Framework-incompatible. **Exclude** |
| **RibonanzaNet** | v1, 2024 | [Shujun-He/RibonanzaNet](https://github.com/Shujun-He/RibonanzaNet) · [multimolecule/ribonanzanet](https://huggingface.co/multimolecule/ribonanzanet) | Kaggle account for originals; **original repo has NO LICENSE file**; MM rewrap AGPL-3.0+ | **11.37 M** | **hidden 256, pairwise 64**, 9 layers | **unbounded** (no learned pos cap) | AUGC+N | 0.05 | yes, trivially | **MSA-free** | Tiny, unbounded length, MSA-free, and gives a **native 64-d pair representation**. Chemical-mapping-pretrained — encodes a physically measured signal no other model has |
| **RibonanzaNet2** | 2025 | [Kaggle model](https://www.kaggle.com/models/shujun717/ribonanzanet2) | licence not stated on the model page | larger, unpublished | seq + pairwise | — | same | — | yes | none | Easiest route is **via RNAPro**, which wires it up and ships the checkpoint path |
| *"RibonanzaNet-3D"* | — | — | **no standalone release found** | — | — | — | — | — | — | — | DasLab GitHub has no such repo — only `k1_tools`. The 3D work landed as **RNAPro**. Treat the name as literature-only |
| **RNAformer** | ICLR-W 2024 | [automl/RNAformer](https://github.com/automl/RNAformer) | open, permissive | 8 M / 32 M | **2D latent pairing matrix** (axial attention) — the pair rep *is* the output | — | single sequence | small | yes | **none** | Secondary structure only, but the cleanest pure *pair-representation* RNA encoder if the shared space has a contact channel |
| **RiboSphere** | **ICML 2026**, [arXiv 2603.19636](https://arxiv.org/abs/2603.19636) | [Zhangz312/RiboSphere](https://github.com/Zhangz312/RiboSphere) | code released, CC-BY-4.0 | unpublished | **discrete FSQ codes** over an SE(3)-invariant geometric encoder; codebook size unpublished | — | 3D structure | small | likely yes (unverified) | none | **An RNA structure *tokenizer*** — the analogue of FoldSeek 3Di / ESM3 structure tokens. If the shared space wants one discrete vocabulary across modalities, this is the RNA one. Recon 1.25 A / TM 0.84 — **one group, unreplicated** |
| **ATOMICA** | 2025, bioRxiv v4 | [mims-harvard/ATOMICA](https://github.com/mims-harvard/ATOMICA) · [ada-f/ATOMICA](https://huggingface.co/ada-f/ATOMICA) | ungated; code **MIT**, weights **CC-BY-4.0** (both verified) | unpublished | **8 named reps at atom / block / interface / graph level**, selectable via `--guidance` | interface-local | PDB files -> Parquet | — | yes; explicitly needs **no** host CUDA toolkit, no torch-scatter/cluster | none | Zitnik lab. One embedding space over **protein–small molecule, protein–ion, small molecule–small molecule, protein–protein, protein–peptide, protein–RNA, protein–DNA and nucleic-acid–small-molecule** complexes, trained on 2,037,972 complexes. Ships ATOMICA-Interface, an **RNAGlib** variant, and **finetuned ligand models for HEM and HEC** plus Fe/Mg/Zn/etc. **Closest existing artifact to our target, and its heme head is directly on this repo's CYP3A4 problem.** Needs structures, so pair it with a folder |
| **RhoDesign / RiboDiffusion** | 2024 | [ml4bio/RhoDesign](https://github.com/ml4bio/RhoDesign) · [RiboDiffusion](https://github.com/ml4bio/RiboDiffusion) | open | — | **GVP node/edge features** | — | 3D backbone -> seq | small | yes | none | The *structure-side* RNA encoder to pair with sequence-side LMs; gives the decoder end of a shared space |
| **RNA-FrameFlow** | TMLR 2025 | [rish-16/rna-backbone-design](https://github.com/rish-16/rna-backbone-design) | open | small | SE(3) frames, 13 atoms/nt | 40–150 nt validated | de novo | small | yes | none | Synthetic structure generator to pretrain a geometric encoder. ">40% valid" is **self-consistency scTM, not experiment** |
| **RNAbpFlow** | **Nat.Meth. 2026** | [Bhattacharya-Lab/RNAbpFlow](https://github.com/Bhattacharya-Lab/RNAbpFlow) | open | small | SE(3) frames, **base-pair-conditioned** | ~70 nt; 10 samples ~25 s | seq + bp constraints | small | yes | none (needs an SS predictor) | Cleanest conditional **conformational-ensemble** generator; better vetted than most rows |
| **RNAFlow** | ICML 2024 | [divnori/rnaflow](https://github.com/divnori/rnaflow) | open | — | inherits RF2NA internals | — | protein target -> RNA | — | no — RF2NA/DGL stack | RF2NA DBs | Protein-conditioned RNA design; carries RF2NA's whole burden |
| **RNAJP / IsRNA / IsRNA2 / IsRNAcirc** | 2021–24 | Chen lab (Missouri) / bioRxiv | academic | n/a | **no learned embeddings** — coarse-grained MD | — | seq+SS | — | CPU | — | Ensemble / negative-sample generators only. **Not embedding models** |
| **DeepFoldRNA** | 2022–23 | [zhanggroup.org/DeepFoldRNA](https://zhanggroup.org/DeepFoldRNA/) | downloadable | — | restraint predictor, not a clean rep | — | MSA | — | yes | heavy | Superseded by DRfold2 from the same group |
| **E2Efold-3D** | 2022 | redirects to RhoFold | — | — | — | — | — | — | — | — | **Deprecated — it *is* RhoFold's ancestor** |
| **RNA-GPT** | NeurIPS-W 2024 | [arXiv 2411.08900](https://arxiv.org/abs/2411.08900) | code availability unverified | — | RNA encoder + **linear projection into an LLM token space** | — | RNA + text | — | ? | — | Architecturally *exactly* the adapter pattern we want (frozen encoder -> linear projector -> shared space), but it targets text QA. Design reference, weak component |
| **boltz-as-FM** (tooling) | [arXiv 2602.13249](https://arxiv.org/abs/2602.13249) | [hsjang0/boltz-as-FM](https://github.com/hsjang0/boltz-as-FM) | MIT | — | per-layer pair reps (16/32/48/64) + hybrid pooling | — | — | — | yes | — | Only needed if the shipped `--write_embeddings` final-layer `s`/`z` proves insufficient |

### Notes on section 3

**The dimensional coincidence to exploit.** RhoFold+ (`c_s=384, c_z=128`), Boltz-1/1x/2
(`single 384, pair 128`) and the AF3-family trunks all land on **384 single / 128 pair**. Three
independently trained RNA-capable models with identical rep shapes means **one shared adapter
(384→d, 128→d_pair) covers all of them with no per-model surgery.** RNA-FM (640), DRfold2
(512/128) and RibonanzaNet (256/64) need their own projections.

**The field-level negative result — read this before believing any accuracy number above.**
Two independent 2025–26 assessments say RNA 3D prediction does not generalise:

- *Limits of deep-learning-based RNA prediction methods*, [NAR 2026 54(16):gkag813](https://academic.oup.com/nar/article/54/16/gkag813/8769247),
  benchmarked AF3, Boltz-1, Chai-1, DeepFoldRNA, DRfold, HelixFold3, NuFold, RhoFold+, RF2NA,
  trRosettaRNA. **Easy targets TM 0.342–0.482; hard (training-dissimilar) targets TM 0.198–0.247
  for every method.** 12 of 79 structures passed all four metrics. Conclusion: *"methods
  recognize recurring motifs rather than generalize to novel RNA folds."* High TM with low DockQ
  was common — correct global fold, RNA docked to the wrong protein surface.
- [CASP16 nucleic-acid assessment](https://pmc.ncbi.nlm.nih.gov/articles/PMC12248019/): **no
  previously-unseen natural RNA exceeded TM 0.8 by anyone.** Of 36 monomer targets where any
  group beat TM 0.45, only **2** lacked a usable template. The AF3 *server* ranked poorly.

For a world model this is arguably good news — it says the embeddings encode *motif lookup*, and
motif lookup is what a shared retrieval space is for. But per this repo's kill criterion 3, do
not quote a single-model accuracy claim as evidence of generalisation.

**On the "Kong lab".** There is no Kong lab. RNA-FM and RhoFold+ are the **`ml4bio` group / CUHK
AIH Lab** at the Chinese University of Hong Kong; PI **Yu Li**, senior co-author **Irwin King**.
"Kong" is almost certainly a corruption of "King" or of "Hong Kong" in the affiliation string —
and **Liang Hong** *is* a real RNA-FM co-author. Not involved: Shuangjia Zheng (SJTU), Yi Xiong
(SJTU). **Yuedong Yang** (Sun Yat-sen) is real but belongs to **RNA-MSM**, not RNA-FM. The ml4bio
stack is one coherent four-tool ecosystem: **RNA-FM → RhoFold+ → RiboDiffusion → RhoDesign**.

**Stale forks outrank canonical repos in search.** Pin org names in automation. Canonical:
`jwohlwend/boltz`, `bytedance/Protenix`, `chaidiscovery/chai-lab`, `ml4bio/RhoFold`,
`aqlaboratory/openfold-3`. Search returned `BLUE-Flowing/boltz-2`, `fuad021/boltz2`,
`Maikuraky/Chai-1`, `Dharmogata/RhoFold` (a pre-"+" snapshot) and others *above* the real ones.

### MSA dependency — the airgap axis

This repo has already lost 12 days to `--use_msa_server` on an internet-less GPU partition. The
field splits cleanly:

| class | models | DB footprint |
|---|---|---|
| **MSA-free** | **DRfold2**, RibonanzaNet/2, RNA-FM, RiNALMo, ERNIE-RNA, RNAformer, RiboSphere, RhoDesign, RiboDiffusion, RNA-FrameFlow, RNAbpFlow | **0** |
| **degraded single-seq mode** | RhoFold+, trRosettaRNA2, Boltz-1/1x/2 (`msa: empty`), Chai-1, Protenix | 0 |
| **MSA required** | trRosettaRNA2 (full) | 32 GB — cheapest |
| | RNAPro | user-supplied, stage-friendly |
| | HelixFold3 | 190 → 530 GB |
| | RoseTTAFold2NA | ~480 GB (nt alone 151 GB) |
| | AlphaFold3 | ~630 GB |
| | RhoFold+ (full accuracy) | ~900 GB |
| | **NuFold** | **~2 TB — disqualifying** |

Counter-pressure from CASP16: MSA depth is what rescues template-free targets, so single-sequence
mode is not free. **Budget MSA staging as its own task, or restrict the pool to the MSA-free column.**

---

## 4. RNA structure and sequence data at scale

Counts pulled live on **2026-09-20** from the resource's own release notes or FTP listing, or
from the RCSB Search API, unless marked otherwise.

| Resource | Ver / date | URL | Open? | License | Size (exact counts) | Format | DL GB | Use in shared space |
|---|---|---|---|---|---|---|---|---|
| **RNAcentral** | **release 27, 2026-07-20** | [rnacentral.org](https://rnacentral.org/) · [FTP](https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/) | yes, no registration | **CC0** | **58,558,809 unique ncRNA sequences; 266,374,326 cross-refs; 59 expert databases** | FASTA, JSON, TSV, Postgres dumps | **10.0** (`rnacentral_active.fasta.gz`), +1.2 inactive, +11 species-specific | ⭐ The pretraining corpus every RNA LM in §1 was trained on. Stage once to `/scratch`; embed a deduplicated subset, not all 58 M |
| **Rfam** | **15.1, January 2026** | [rfam.org](https://rfam.org/) · [FTP](https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/) | yes | **CC0** | **4,227 families** with seed alignments, consensus SS and covariance models | Stockholm, CM, MySQL | seeds <0.5; full-region hits ~10s | Family labels are the only broad functional annotation RNA has. Use as the **held-out split key** — leave-one-family-out is RNA's analogue of leave-one-target-out |
| **PDB — RNA-containing** | live | [RCSB search API](https://search.rcsb.org) | yes | **CC0** | **10,399 entries** with ≥1 RNA entity; **2,375 RNA-only**; **8,024 protein–RNA**; 574 RNA-only entries released since 2023 | mmCIF | ~40 for the RNA subset | Ground truth for every structural claim. **Everything in §6 is a re-cut of this** |
| **RNA3DB** | **2026-01-05 full release** | [marcellszi/rna3db](https://github.com/marcellszi/rna3db) | yes, GitHub releases | **MIT** | Non-redundant clustering of all PDB RNA chains with **pre-built train/test splits**, plus Infernal cmscan results for every chain | mmCIF + JSON | **2.15** (`rna3db-mmcifs.tar.xz`), 0.088 cmscans, 0.007 jsons | ⭐ **The split we should use.** Someone has already defended a non-redundant RNA structural split by sequence *and* Infernal homology. Do not reinvent it — this repo's rules demand honest gating, and RNA3DB is the honest RNA split |
| **RNAsolo / RNAsolo2** | live | [rnasolo.cs.put.poznan.pl](https://rnasolo.cs.put.poznan.pl/) | yes | not stated on the front page | Cleaned PDB-derived RNA 3D structures, organised by **BGSU equivalence class**, with resolution-cutoff and redundancy-filtered subsets | PDB, mmCIF, extracted chains | ~GB | Convenient pre-cleaned download of the same PDB data. ⚠️ **Not independent of the PDB** — it is a cleaning layer. Front page returns almost no text to automated fetchers; check the version by hand |
| **BGSU RNA Structure Atlas / RNA3DHub** | live, weekly | [rna.bgsu.edu/rna3dhub](http://rna.bgsu.edu/rna3dhub/) | yes | academic | Representative sets and non-redundant lists by resolution; **RNA 3D Motif Atlas** (internal-loop and hairpin-loop motif groups) | lists + coords | small | The canonical **equivalence-class** definition. Use BGSU classes (or RNA3DB clusters) whenever reporting a "non-redundant" RNA number — otherwise the number is inflated |
| **bpRNA / bpRNA-1m** | 1.0 | [bprna.cgrb.oregonstate.edu](https://bprna.cgrb.oregonstate.edu/) | yes | academic | **102,318 sequences** = CRW 55,600 + Rfam 43,273 + PDB 669 + tmRNA 728 + SRP 959 + tRNA 623 + RNase P 466. Annotated: 2,075,928 stems, 708,144 hairpins, 538,670 internal loops, 517,672 bulges, 317,046 multiloops, **57,686 pseudoknots** | dbn/ST | <1 | Standard secondary-structure corpus. ⚠️ **Only 669 of the 102,318 come from experimental 3D structures** — the rest are covariance-model predictions from Rfam/CRW. A "100 k structure dataset" that is really ~700 experimental structures plus comparative-sequence inference |
| **bpRNA-1m(90)** | 1.0 | same | yes | academic | 90%-identity-filtered subset (~28 k) | dbn | <1 | The version to use; the unfiltered one leaks between train and test |
| **bpRNA-new** | 2021 | via SPOT-RNA2 / UFold papers | yes | academic | Rfam families **added after** bpRNA-1m — built specifically as a cross-family generalisation test | dbn | <1 | ⭐ The only SS test set that is family-disjoint from training. Use it, not ArchiveII |
| **ArchiveII / RNAStralign** | 2010s | via ML papers | yes | academic | ArchiveII ~3,975 sequences, 10 families; RNAStralign ~37,000, 8 families | ct | <1 | ⚠️ **Massively redundant and family-overlapping with every training set.** Near-perfect scores here are memorisation. Per kill criterion 3, an uninformative benchmark |
| **Ribonanza** | Das lab, bioRxiv 2024 / Kaggle 2023 | [paper](https://www.biorxiv.org/content/10.1101/2024.02.24.581671v2) · [DataPrepRibonanzaKaggle2023](https://github.com/DasLab/DataPrepRibonanzaKaggle2023) | yes | competition + open data | **Chemical mapping on ~2 million diverse RNA sequences** (Eterna + crowdsourced). Kaggle train: **214,831 sequences** with a 2A3 or DMS profile at SNR>1.0; **Ribonanza+ post-competition: 494,111** at SNR>1.0; +563 k noisy (SNR<1); RibonanzaNet trained on **2.1 M** | parquet/CSV reactivity profiles | ~10s of GB | ⭐⭐ **The single most valuable RNA dataset for a world model.** It is *experimental* per-nucleotide structure signal at a scale nothing else in RNA approaches — 2 M sequences vs 2,375 RNA-only PDB entries. 2A3 and DMS are two chemical probes reporting nucleotide flexibility/accessibility |
| **RMDB (RNA Mapping Database)** | live | [rmdb.stanford.edu](https://rmdb.stanford.edu/) | yes | open | The long-tail SHAPE/DMS/chemical-mapping archive that predates Ribonanza; hundreds of entries, thousands of constructs | RDAT | <1 | Older and much smaller than Ribonanza; use for probe/condition diversity, not scale |
| **OpenKnot / Eterna** | ongoing | Das lab / Eterna | yes | open | Pseudoknot-focused crowdsourced chemical mapping, layered on the Ribonanza pipeline | same | — | ⚠️ **Not independent of Ribonanza** — same platform, same lab, overlapping libraries |
| **CASP15 / CASP16 RNA targets** | 2022 / 2024 | [CASP16 NA assessment](https://pmc.ncbi.nlm.nih.gov/articles/PMC12248019/) · [prot.70072](https://doi.org/10.1002/prot.70072) | yes | open | CASP16: **42 nucleic-acid targets, 65 groups, 46 labs**; CASP15 had 12 RNA targets | PDB | tiny | The only genuinely blind RNA structural evaluation. **Tiny — 42 targets.** Any claim of "beating CASP" rests on a handful of structures |
| **RNA-Puzzles** | rounds I–V | [rnapuzzles.org](https://www.rnapuzzles.org/) | yes | open | Low tens of targets per round | PDB | tiny | Same caveat as CASP |
| **MSA resources: rMSA / RNAcmap3 / Infernal+nt** | — | [RNAcmap3 server](https://aigene.cloudbastion.cn/#/rna-msm) | tool open; DBs are the cost | — | rMSA needs **~2 TB**; RNAcmap3 needs blastn vs `nt` (151 GB) + Infernal vs Rfam | a3m/Stockholm | 0.15–2 TB | **Budget as its own task.** GPU nodes have no internet; MSAs must be built on a login node and staged |
| **GtRNAdb / SILVA / CRW** | live | various | yes | mostly academic | tRNA and rRNA specialist corpora, already largely inside RNAcentral | FASTA/alignment | 1–20 | ⚠️ Mostly **redundant with RNAcentral**, which aggregates 59 expert databases including these |

### Notes on §4

**The scale asymmetry, stated plainly.**

| layer | RNA | protein (for reference) |
|---|---|---|
| sequences | **58.6 M** (RNAcentral r27) | ~250 M (UniRef100) |
| families | **4,227** (Rfam 15.1) | ~20,000 (Pfam) |
| experimental 3D structures | **10,399 entries / 2,375 RNA-only** | ~230,000 |
| *measured* per-nucleotide structure signal | **~2 M sequences** (Ribonanza) | no direct analogue |

RNA's sequence layer is within an order of magnitude of protein's. Its **structure layer is two
orders of magnitude smaller.** Ribonanza is the one place where RNA has *more* experimental data
than protein does in kind — and it is exactly the layer no RNA LM except RibonanzaNet trains on.
**That gap is an opportunity, not a footnote.**

**Redundancy, explicitly.** PDB → RNAsolo → BGSU/RNA3DHub → RNA3DB → HARIBOSS → RNAmigos2 →
fpocketR are **one dataset presented seven ways.** Cite the PDB entry count once; cite everything
downstream as a *curation* of it. Similarly, bpRNA-1m's 102,318 "structures" are 669 experimental
ones plus ~101,600 comparative-sequence inferences, and GtRNAdb/SILVA/CRW are already inside
RNAcentral.

**Total acquisition footprint for everything above: well under 1 TB** (RNAcentral 10 GB, Rfam
seeds <1 GB, PDB RNA subset ~40 GB, RNA3DB 2.2 GB, Ribonanza tens of GB, bpRNA <1 GB) — **unless**
MSA databases are staged, which alone costs 0.15–2 TB. Kill criterion 4 (10 TB cap) is not
threatened by the RNA arm except through MSAs.

---

## 5. Protein–RNA interaction data

All counts below were pulled live on **2026-09-20** (RCSB search API, ENCODE portal API) or
from the resource's own release notes. Where a paper's number and the live site disagree, both
are given.

| Resource | Ver / yr | URL | Open? | License | Size (exact) | Format | DL GB | Use in shared space |
|---|---|---|---|---|---|---|---|---|
| **PDB protein–RNA complexes** | live | [RCSB search API](https://search.rcsb.org) | yes | CC0 | **8,024 entries** with ≥1 protein **and** ≥1 RNA entity | mmCIF | ~40 (full mmCIF subset) | The only source of 3D protein–RNA interfaces. Heavily ribosome-dominated — cluster before use |
| **PDB entries containing RNA (any)** | live | same | yes | CC0 | **10,399** | mmCIF | — | Denominator for every RNA structure claim below |
| **ENCODE eCLIP** | v4 portal, 2020–26 | [encodeproject.org](https://www.encodeproject.org/search/?type=Experiment&assay_title=eCLIP) | yes, no registration | **ENCODE terms (free use, cite)** | **290 experiments, 228 unique targets** (244 classed RBP); K562 164 / HepG2 124 / adrenal 2; **human only** | BED narrowPeak, bigWig, BAM | peaks ~2 GB; BAMs ~15 TB | Positive set for RBP↔transcript binding. **Binding sites, not affinities** |
| **ENCODE RNA Bind-n-Seq** | 2020 | [encodeproject.org](https://www.encodeproject.org/search/?type=Experiment&assay_title=RNA+Bind-n-Seq) | yes | ENCODE terms | **216 experiments, 116 unique targets** (103 RBPs), Burge lab | FASTQ + motif tables | small | *In vitro*, so it separates intrinsic sequence preference from cellular context — the cleanest RBP motif set |
| **POSTAR3** | v3, NAR 2022 (D287) | [paper](https://academic.oup.com/nar/article/50/D1/D287/6414586) · `http://postar.ncrnalab.org` | yes | not stated | **1,499 CLIP-seq datasets, 348 RBPs, 7 species, ~50 M binding sites**; + 300 Ribo-seq, 66 structure-seq, 6 icSHAPE, 83 degradome-seq | BED/TSV | tens of GB | Broadest CLIP compendium and the only multi-species one. ⚠️ **the site is HTTP-only and redirects to a bare IP (111.198.139.65) that refused HTTPS on 2026-09-20** — verify reachability before planning acquisition |
| **ENCORI / starBase** | v3, live | [rnasysu.com/encori](https://rnasysu.com/encori/) | yes, live, has an API | not stated | **2,725 CLIP-seq datasets, 23 species**; ~1.29 M RBP–mRNA + ~1.6 M RBP–ncRNA interactions | web API, TSV | tens of GB | Largest interaction count of any CLIP resource. ⚠️ **it is largely the same GEO/ENCODE raw data as POSTAR3, reprocessed** — do not treat the two as independent evidence |
| **ProNAB** | 2022, IIT Madras | [web.iitm.ac.in/bioinfo2/pronab](https://web.iitm.ac.in/bioinfo2/pronab/) | yes | academic | **20,219 binding-affinity entries** (protein–DNA **and** protein–RNA combined; the RNA share is a minority and the site does not break it out) | web + download | <0.1 | ⭐ The only sizeable protein–nucleic-acid **affinity** set (Kd, ΔG, ΔΔG). This is the protein-side analogue of what §6 shows does not exist for RNA–ligand |
| **CISBP-RNA / RBPDB / ATtRACT / oRNAment** | 2014–2019 | various | yes | academic | motif PWMs for a few hundred RBPs; **derived from RNAcompete/RBNS**, not new measurements | PWM/TSV | <1 | Motif priors only. **All four are the same upstream experiments rewrapped** |
| **RNAcompete** | 2013–2021 | CISBP-RNA | yes | academic | ~200 RBP motifs | PWM | <1 | Superseded by RBNS for depth; still the widest species coverage |

### Honest tally for protein–RNA

| quantity | number |
|---|---|
| 3D protein–RNA complexes (PDB, 2026-09-20) | **8,024 entries** — but heavily redundant (ribosomes, tRNA synthetases, spliceosomes) |
| High-throughput binding-site measurements | **~50 M sites** (POSTAR3) / **~2.9 M interactions** (ENCORI) — **no affinities** |
| RBPs with any binding data | **348** (POSTAR3) / **228** (ENCODE eCLIP) / **116** (RBNS) |
| Protein–nucleic-acid complexes with a **measured Kd/ΔG** | **20,219 entries in ProNAB, protein–DNA and protein–RNA combined** — the protein–RNA subset is low thousands at best |

**The structural point for a JEPA.** Binding-site data is abundant and affinity data is not,
in exactly the same ratio as everywhere else in biology. A JEPA that predicts *which* partner
binds (a retrieval task over 2.9 M positive pairs) is trainable today; a JEPA that predicts
*how tightly* is not, on RNA, in either direction.

**Redundancy warning.** POSTAR3, ENCORI, CLIPdb and doRiNA are largely the **same GEO
submissions reprocessed with different pipelines**. Counting them as four independent corpora
inflates the apparent data by ~4×. Deduplicate by GEO/SRA accession, not by database name.

---

## 6. RNA–small-molecule binding data — the priority section

> **The suspicion was right, and the reality is worse.** Two independent RCSB Search API sweeps
> were run today (mine and the reconnaissance agent's, at different drug-likeness thresholds) and
> they agree: **the entire RNA–small-molecule structural universe is 708 entries at the loose
> threshold and 540 at a stricter one, against 128,528 protein entries with a drug-like ligand.**
> The number of RNA–ligand pairs with **both** a structure and a measured affinity is **101**.
> The number of non-redundant, genuinely drug-like RNA–ligand co-structures ever solved is
> **48**, covering **38 unique RNAs**.

### 6.1 Live PDB counts — run independently twice on 2026-09-20

| Query | This lead's sweep (MW>250, ≥13 heavy atoms) | Agent sweep (MW 150–900, ≥10 heavy) |
|---|---|---|
| All experimental PDB entries | — | 259,693 |
| Entries with ≥1 RNA entity | **10,399** | 10,399 |
| Protein–RNA complexes | **8,024** | — |
| RNA-only entries (no protein) | **2,375** | 2,375 |
| RNA-only + any non-polymer entity (incl. ions/solvent) | 1,466 entries / **2,896 entities** | — |
| **RNA-only + drug-like ligand** | **540 entries / 595 ligand instances** | **708 entries** |
| Unique ligand chem comps, RNA-only | **497** | 402 |
| Unique ligand chem comps, any RNA-containing entry | **1,075** | 1,150 |
| …of the RNA-only drug-like ligand instances, **riboswitch** | **243** | 330 entries / 228 comps |
| …**aptamer** | **155** | 237 entries |
| RNA-only entries released since 2023-01-01 | **574** | — |
| Protein entries with a drug-like ligand (scale reference) | — | **128,528** |

**Read the riboswitch/aptamer rows carefully.** Of ~595 drug-like ligand instances bound to
protein-free RNA, **~two-thirds are riboswitches binding their natural cognate metabolite, or
engineered fluorogenic aptamers binding their designed dye.** These are not medicinal-chemistry
series. They are one-ligand-per-target co-crystals. The med-chem-style SAR series that make
PDBbind useful for protein–ligand ML essentially do not exist for RNA.

And of the *loose* 3,499 entries that contain RNA and a drug-like ligand, **1,617 (46%) are
ribosomes** — nearly half of all "RNA–ligand structural biology" is antibiotics bound to one
macromolecular machine.

### 6.2 The table

**SIZE legend:** **[A]** measured affinities · **[S]** structures only · **[H]** binary hits ·
**[X]** annotations only.

| Name | Ver | Yr | URL / DOI | Open? License | SIZE (precise) | Format | DL | Use |
|---|---|---|---|---|---|---|---|---|
| **HARIBOSS** | live | 2022→26 | [hariboss.pasteur.cloud](https://hariboss.pasteur.cloud/) · [btac483](https://doi.org/10.1093/bioinformatics/btac483) | free web; **licence not stated** | **2,077 complexes, 4,830 pockets, 537 unique ligands [S]** (live today). At publication: 716 redundant / **484 non-redundant**, 1,226 pockets, 267 ligands | web + per-entry files | <1 GB | ⭐ **The canonical RNA–ligand structure set. Start here.** Note: `www.` prefix does **not** resolve; use the bare host |
| **RNAmigos 2** | 2 | 2025 | [Nat.Comm.](https://doi.org/10.1038/s41467-025-57852-0) · [cgoliver/rnamigos2](https://github.com/cgoliver/rnamigos2) · [Zenodo 14803961](https://doi.org/10.5281/zenodo.14803961) | **open, MIT** | **1,740 RNA–ligand binding sites**, clustered by RMAlign ≥0.75 into **436 groups → 367 train / 69 test [S]**; + 1.3 M synthetic docking scores, 1.5 M unlabelled compounds | graphs + Zenodo | ~GB | ⭐ Best-engineered RNA virtual-screening benchmark with **pre-defended structural splits**. Prospectively validated on a 20,000-compound microarray (EF₁% = 2.93) |
| **ROBIN** | 1 | 2022/23 | [ky66/ROBIN](https://github.com/ky66/ROBIN) (MIT) · [Figshare 20401974](https://figshare.com/articles/dataset/Machine_Learning_Informs_RNA-Binding_Chemical_Space/20401974) (CC BY 4.0) · [Angew.](https://doi.org/10.1002/anie.202211358) | **fully open** | **24,572 compounds × 36 SMM screens (27 RNA + 9 DNA) = 1,627,072 assayed interactions [H]; 2,003 RNA-binding hits.** Hit rates 0.48% (hairpin) → 0.85% (pseudoknot). **ZERO Kd — binary only** | SDF + Mordred CSV | **1.28 GB** | ⭐ The **only genuinely large independent** open RNA-binding dataset. ⚠️ its negatives are BindingDB *protein* binders → documented physicochemical shortcut; pair with RNAdecoyDB or your AUROC measures molecular weight |
| **R-SIM** | 1 | 2023 | [web.iitm.ac.in/bioinfo2/R_SIM](https://web.iitm.ac.in/bioinfo2/R_SIM/index.html) · [JMB](https://doi.org/10.1016/j.jmb.2022.167914) | free web; no explicit licence; paper paywalled | **2,501 experimentally validated RNA–SM interactions with binding affinity [A]** — the largest RNA affinity set that exists. **Explicitly assembled from SMMRNA + R-BIND + literature** | web DB | MB | ⭐ **This is the RNA analogue of BindingDB — and it is ~1,300× smaller** |
| **RSAPred** train set | — | 2024 | [Brief.Bioinform. bbae002](https://doi.org/10.1093/bib/bbae002) · [server](https://web.iitm.ac.in/bioinfo2/RSAPred/) | CC BY-NC paper; server free | Kd-only subset of R-SIM: **aptamers 516 / miRNA 146 / repeats 97 / rRNA 294 / riboswitch 101 / viral 326 = 1,480 pairs [A]; 341 unique RNA targets; 897 unique small molecules** | CSV/SI | <10 MB | ⭐ The only usable RNA pKd **regression** benchmark. ⚠️ **11 rRNA targets carry 294 of the points** — stratify or you are fitting the ribosome |
| **PDBbind — RNA subset** | v2020 lineage | 2020– | [pdbbind.org.cn](http://www.pdbbind.org.cn) (**socket hang-up today**) / pdbbind-plus.org.cn | academic **registration** | **101 RNA–ligand pairs with affinity + structure [A+S]** (confirmed independently by the JCIM 2026 review and RNAmigos2) | PDBbind format | MB | The apples-to-apples number: **101 vs ~20,000 protein** |
| **fpocketR / Weeks drug-like set** | 1 | 2025 | [PNAS](https://doi.org/10.1073/pnas.2422346122) · [PMC12054788](https://pmc.ncbi.nlm.nih.gov/articles/PMC12054788/) | CC BY-NC-ND paper; tool open | HARIBOSS non-redundant (538) filtered to QED≥0.3, ≤3.5 Å, unique RNA, ≤200 nt → **48 drug-like RNA–ligand complexes; 38 unique RNAs (13 train / 25 test); 10 holo/apo pairs.** Authors: "**roughly 50** drug-like pockets have ever been visualised" | PDB/CSV | MB | ⭐ **The most honest number in the field. Quote this one.** Any RNA drug-likeness generalisation claim is being tested on ~25 structures |
| **FURNA** | 1 | 2024 | [PLOS Biol.](https://doi.org/10.1371/journal.pbio.3002476) · [server](https://seq2fun.dcmb.med.umich.edu/furna/) | **CC BY, bulk download provided** | All PDB RNA chains ≥10 nt (**>16,000 chains**) with ligand-site annotations, GO/EC/Rfam, cross-refs to ChEMBL/DrugBank/ZINC [S/X] | TSV + coords | ~GB | The RNA analogue of BioLiP2 (same group). **Most complete RNA annotation layer, free and bulk-downloadable** |
| **R-BIND** | 1.0 / 1.2 / **2.0** | 2017 / 2019 / **2022** | [2.0 PMC9343015](https://pmc.ncbi.nlm.nih.gov/articles/PMC9343015/) · [hargrovelab.org/rbind](https://hargrovelab.org/rbind) | ⚠️ **`rbind.chem.duke.edu` is DEAD (DNS fails).** Data now only via **paywalled ACS Supporting Information**. No licence stated | 1.0 = 104 ligands; 1.2 = 135; **2.0 = 188 (153 SM + 35 MW>700)**, RNA targets grouped into 5 structure classes. **The distinct-target count and the affinity-vs-binder-only split are NOT published** — they are inside the SI spreadsheet | XLSX | <1 MB | Best-curated "RNA-privileged chemistry" set. **R-BIND 3.0 does not exist** — the lab page says "update coming soon", undated; the lab moved Duke→Toronto |
| **DRTL** (Duke RNA-Targeted Library) | — | 2023 | [bioRxiv 2023.07.31.551350](https://doi.org/10.1101/2023.07.31.551350) — **still a preprint in 2026** | no public catalog | **>800 compounds**, a physical deck matched to R-BIND property space; 4 RNA targets screened [H] | physical | — | Buy/clone an RNA-biased screening deck |
| **Inforna** | **2.0** | 2016 | [PMC4912454](https://pmc.ncbi.nlm.nih.gov/articles/PMC4912454/) · [disney.scripps.ufl.edu/software](https://disney.scripps.ufl.edu/software/) | ⚠️ **GATED — a signed academic licensing agreement is required.** Reviews calling it "freely available" are wrong | **1,936 RNA-motif–small-molecule interactions; 244 unique small molecules; 1,331 motifs.** Scored by **StARTS Z-score from 2DCS, not per-entry Kd** | web tool, **no bulk download** | n/a | Second genuinely independent body of RNA-binding data — and it is licence-gated and not affinity-valued. **No version >2.0 exists.** "SMIRNAs" is lab terminology, not a database |
| **RNALigands** | 1 | 2022 | [PMC8906548](https://pmc.ncbi.nlm.nih.gov/articles/PMC8906548/) · [SaisaiSun/RNALigands](https://github.com/SaisaiSun/RNALigands) | CC BY-NC | **841 motif–ligand pairs = PDB 386 + R-BIND 67 + Inforna 388** (their own Fig. 1A) | web + code | small | ⭐ **The clearest single proof of the double-counting problem — it is a literal merge of three other databases** |
| **RPocket** | 1 | 2021 | [PMC8424408](https://pmc.ncbi.nlm.nih.gov/articles/PMC8424408/) | free | **240 RNA pockets from 94 non-redundant complexes [S]** | web | small | Pocket geometry descriptors |
| **RNASite TR60 / TE18** | — | 2021– | via RLBind/RNet/MultiModRLBP papers | open | **TR60 = 60 structures (train), TE18 = 18 (test)** — the field's de-facto binding-site benchmark, **78 structures total** | PDB | MB | ⚠️ A "benchmark" of 78 structures. Every RNA binding-site paper reports on it |
| **RNAdecoyDB + TRN** | 1 | 2026 | [JCIM 10.1021/acs.jcim.6c00360](https://doi.org/10.1021/acs.jcim.6c00360) — no PMCID yet | **availability unverified** | Property-matched **hard decoys** for RNA VS; demonstrates ROBIN's negatives are physicochemically mismatched | — | — | ⚠️ **Read before training on ROBIN** |
| **DRLiPS** | — | 2025 | [NAR gkaf239](https://doi.org/10.1093/nar/gkaf239) | open | druggable-pocket classifier with novel negative sampling; blind test on apo + modelled RNA | — | — | Druggability triage |
| **ChEMBL_37 RNA slice** | 37 | 2026-05 | [ebi.ac.uk/chembl](https://www.ebi.ac.uk/chembl/) | CC BY-SA | **129 nucleic-acid targets; 84 RNA-like; only 43 have ≥1 activity. 66,328 activities — but 64,732 (97.6%) are against miR-21 alone** (reporter/expression assays). **Excluding miR-21: ≈1,596 activities.** Top genuine entries: IRES 229, miR-30a 159, rRNA A-site 154, RRE 119, r(CUG) 80, PreQ1 riboswitch 57 | API/SQL | — | ⚠️ **Textbook "tiny dataset dressed up as large."** 66 k collapses to 1.6 k |
| **BindingDB** | live | 2026 | [bindingdb.org](https://www.bindingdb.org/) | free, CC | 3.2 M measurements, 1.4 M compounds, 11,500 targets — **no RNA target class at all** | — | GB | Scale reference only |
| **PLINDER** | live | 2024–26 | [plinder.sh](https://www.plinder.sh/) | open | >400 k protein–ligand systems, >50 k unique small molecules. **No RNA/nucleic-acid analogue exists** | — | ~TB | Scale reference |
| **SM2miR** | — | 2013 | [Bioinformatics bts698](https://doi.org/10.1093/bioinformatics/bts698) | free | 2,925 small-molecule–miRNA "interactions" — but these are **expression-modulation effects, not binding** | web/flat | small | ⚠️ **Do NOT treat as binding data.** It is routinely miscited as such |
| **Ribocentre-aptamer** | 1 | 2026 | [NAR gkaf1016](https://doi.org/10.1093/nar/gkaf1016) · [aptamer.ribocentre.org](https://aptamer.ribocentre.org/) | free, batch export; licence not stated | **191 aptamer types / 669 publications / 510 sequences / 123 small-molecule ligands / 344 protein ligands**, with structures, pockets and **affinity metrics [A/S/X]** | web + batch export | small | ⭐ Newest and best aptamer resource. **123 small-molecule ligands is the realistic ceiling for aptamer–SM pairs** |
| **Apta-Index** (Aptagen) | live | 2026 | [aptagen.com/apta-index](https://www.aptagen.com/apta-index/) | free to browse, **no bulk download**, commercial vendor | Hundreds of entries, **Kd per entry**; small-organic targets are a small minority | HTML only | n/a | Hand-mining only |
| UTexas Aptamer DB / AptamerBase / RNAapt3D / AptaDB | — | legacy | — | **liveness in 2026 unverified** | — | — | — | Treat as legacy; superseded by Ribocentre-aptamer |
| **HT-SELEX** | — | — | Jolma/Taipale lineage | open | ⚠️ HT-SELEX is **protein**-target. **Small-molecule SELEX is not high-throughput** — one target per campaign | — | — | Do not expect scale here |
| **Riboswitch affinities** | — | — | — | — | ⚠️ **NO aggregated riboswitch-affinity database exists.** The only aggregation is R-SIM's slice: **101 pairs / 34 riboswitch targets / 63 ligands**. ~55–60 riboswitch classes are experimentally validated (Breaker 2023); not all have a structure | — | — | **A real, novel gap. Curating this would be new work, not re-curation** |
| **CASP16 nucleic acids** | 16 | 2024 exp / 2026 papers | [Proteins 10.1002/prot.70072](https://doi.org/10.1002/prot.70072) | open | 42 NA targets, 65 groups, 46 labs. **No prediction of unseen natural RNA reached TM > 0.8.** **There was NO RNA–ligand scoring category** — the LG category was protein-ligand | — | — | Companion paper title: *"Ligand affinity prediction is totally unreliable"* |
| *"RNA-LigandBench" / "RLABench" / "PoseBusters-RNA"* | — | — | — | — | ⚠️ **NO SUCH RESOURCES EXIST.** These names return nothing. PoseBusters/PoseBench remain protein-ligand. **There is no RNA pose-validity suite and no named RNA-ligand benchmark suite** | — | — | Do not design anything that assumes one |
| Arrakis / Expansion / Skyhawk / Ribometrix / Remix | — | — | — | — | **No open compound–RNA affinity datasets found.** PEARL-seq / Chem-CLIP produce **target-ID maps, not quantitative affinity tables** | — | — | Nothing reusable |
| DEL vs RNA | — | 2026 review | [PMC13019106](https://pmc.ncbi.nlm.nih.gov/articles/PMC13019106/) | review only | DEL-vs-RNA is discussed as emerging; **no public dataset located.** DNA-barcode/nucleic-acid cross-reactivity is the structural obstacle | — | — | — |
| NMR fragment screening vs RNA | — | 2026 | [PMC13029659](https://pmc.ncbi.nlm.nih.gov/articles/PMC13029659/) | open review | Per-paper fragment hits; **no aggregated fragment-vs-RNA database** | — | — | — |
| **Boltz-2 affinity on RNA** | — | 2025/26 | — | MIT weights | ⚠️ **Boltz-2's affinity module was never evaluated on RNA by its authors.** The JCIM 2026 review ran it on 6 ROBIN targets (SAM_II, ZTP, TPP, PreQ1, NRAS, RRE2B) — the first published RNA assessment | — | — | Directly relevant to a co-folding team: **the affinity head we already use has no RNA validation** |

### 6.3 The honest tally

**(a) How many RNA–small-molecule 3D complexes exist**

| definition | count |
|---|---|
| Loose (any RNA entry + 150–900 Da ligand, incl. protein-bound ligands) | 3,499 — **of which 1,617 (46%) are ribosomes** |
| Curated with a real RNA contact (HARIBOSS, today) | **2,077** (4,830 pockets) |
| Strict: RNA-only entries with a drug-like ligand | **708** (loose threshold) / **540** (MW>250, ≥13 heavy) |
| Non-redundant ML-usable binding sites (RNAmigos2) | 1,740 sites → **436 structural clusters** |
| Non-redundant HARIBOSS (2025) | 538 |
| **Non-redundant, drug-like, ≤3.5 Å, non-ribosomal, unique RNA** | **48 complexes / 38 unique RNAs / ~50 pockets, ever** |

**(b) How many have a measured binding affinity**

| | count |
|---|---|
| RNA–ligand pairs with **both** structure and affinity (PDBbind RNA) | **101** |
| RNA–ligand affinity measurements in aggregated form (R-SIM, largest anywhere) | **2,501** |
| …usable Kd after filtering (RSAPred) | **1,480** |
| ChEMBL_37 activities vs RNA targets, excluding the miR-21 artefact | **≈1,596** across 83 targets (only 43 have any data) |
| Riboswitch–ligand pairs with an aggregated Kd | **101** (34 targets, 63 ligands) |
| Aptamer–small-molecule ligands with structure + affinity | **123** |
| ROBIN | **0 affinities** (2,003 binary hits) |

> **The single sentence: roughly 2,500 RNA–small-molecule affinity measurements exist worldwide
> in any aggregated form, of which about 101 come with a 3D structure, and only about 48
> non-redundant drug-like RNA–ligand co-structures have ever been solved.**

**(c) Unique targets and ligands**

- Unique protein-free RNA polymer entities with a drug-like ligand: **845**
- Unique ligand chem comps: **402–497** (RNA-only) / **537** (HARIBOSS) / **1,075–1,150** (any RNA entry)
- Unique RNA targets with a measured Kd: **341** — but **11 rRNA targets carry 294 of the measurements**
- Unique small molecules with a measured Kd vs RNA: **897**
- Unique molecules ever *reported* to bind RNA on any evidence: **~3,000–4,000** (union of ROBIN 2,003 + Inforna 244 + R-BIND 188 + PDB 537 + tail). Only a few hundred have a real Kd against a defined RNA; **only ~50 have a drug-like co-structure**

**(d) How much of it is the same data rewrapped — almost all of it**

Two lineages account for nearly the whole field:

- **Lineage 1, "the PDB re-cut N times."** HARIBOSS (2,077) → RNAmigos2 (1,740 / 436) → fpocketR
  (48) → RPocket (94 / 240) → RNASite TR60/TE18 (78) → RLBind / RNet / MultiModRLBP / CapBind /
  DRLiPS / DeepRSMA (all TR60/TE18 or HARIBOSS re-splits) → FURNA (all PDB RNA chains) →
  PDBbind-RNA (101). **Distinct underlying PDB entries behind all of it: ~2,000–3,100. There is
  no independent structural data anywhere in that list.**
- **Lineage 2, "R-BIND + Inforna + SMMRNA re-merged."** RNALigands states it outright:
  841 = 386 (PDB) + 67 (R-BIND) + 388 (Inforna). R-SIM collates SMMRNA + R-BIND + literature.
  RSAPred trains on R-SIM. DLRNA-BERTa uses RSAPred's subtypes. **R-BIND's 188 ligands therefore
  appear, whole or in part, inside RNALigands, R-SIM, RSAPred and DLRNA-BERTa.**

The **only genuinely independent large dataset is ROBIN** — and it has no affinities and a
documented negative-set shortcut. The second genuinely independent body is Inforna's 2DCS — and
it is licence-gated and Z-score-valued, not Kd-valued. **No cross-database deduplication has ever
been published.**

**(e) Scale gap against protein–ligand**

| axis | protein | RNA | gap |
|---|---|---|---|
| PDB entries with a drug-like ligand | **128,528** | 708 strict / 3,499 loose | **37–182×** |
| Curated complexes with affinity (PDBbind) | ~20,000 | **101** | **~200×** |
| Non-redundant drug-like co-structures | tens of thousands | **48** | **>400×** |
| Binding measurements (BindingDB vs R-SIM) | **3,200,000** | **2,501** | **~1,280×** |
| ChEMBL activities (total vs RNA ex-miR-21) | 24,527,044 | **≈1,596** | **~15,000×** |
| Benchmark infrastructure | PDBbind, CASF, PoseBusters, PLINDER, Binding MOAD | HARIBOSS + RNAmigos2 + ROBIN; **no PoseBusters-RNA, no PLINDER-RNA, no CASP RNA-ligand category** | qualitative |

> **RNA is to protein–ligand roughly what one well-studied protein family is to the whole PDB.
> The entire RNA affinity corpus is smaller than the bioactivity record of many single protein
> targets in ChEMBL.**

### 6.4 What this means for the shared embedding space

1. **Kill criterion 1 cannot be tested at scale on RNA.** "Beat a Morgan-fingerprint PCM baseline
   on a held-out RNA–ligand affinity task" has **1,480 data points**, of which 294 are 11 rRNA
   targets. Per FINDING 007's logic, the noise floor on n≈1,500 with heavy target imbalance will
   swallow most effects. Run it — but pre-register the null and report it as a bounded test.
2. **The usable framing is transfer, not RNA-native training.** RNA cannot supply the affinity
   data. What it *can* supply is (i) **1.6 M binary interaction labels** from ROBIN for a
   classification-flavoured JEPA, (ii) **2,077 structures** from HARIBOSS for geometry, and
   (iii) **~2.9 M protein–RNA interaction pairs** (§5) for partner-retrieval. Affinity regression
   should be learned on the protein side and *transferred* to RNA, which is precisely the
   argument for a shared space in the first place.
3. **The single most defensible RNA experiment available to us:** take the shared space trained
   on protein–ligand affinity, freeze it, and test zero-shot ranking on the **48-complex drug-like
   RNA set** and the **69 held-out RNAmigos2 clusters**. If it transfers at all, that is a real
   result on a real gap. If it does not, that is a clean negative and it costs days, not months.
4. **Boltz-2's affinity head has never been validated on RNA.** Since this repo already runs
   Boltz-2, running its affinity module across HARIBOSS/ROBIN targets is a cheap, novel,
   publishable-scale check that nobody else has done except one 6-target sample in a review.

---

## 7. The concrete RNA plan for the shared space

Ordered by evidence-per-day, and written against the README's Revision 1 finding that
*concatenated unimodal embeddings lose to a count fingerprint* — so nothing below proposes
gluing an RNA vector to a ligand vector and calling it a shared space.

| # | Step | Cost | What it decides |
|---|---|---|---|
| 1 | **Test LucaOne and ATOMICA as the null hypothesis.** LucaOne is already a joint DNA+RNA+protein space. ATOMICA is already one interface-embedding space over eight complex types including protein–RNA, protein–DNA and nucleic-acid–small-molecule, trained on 2,037,972 complexes, ungated, code MIT / weights CC-BY-4.0, **and it ships finetuned HEM and HEC ligand models** | ~2 days | Whether the artifact we are proposing to build already exists. **This is the single highest-value experiment in the whole RNA arm — and its heme head lands directly on this repo's CYP3A4 problem, not only on RNA.** If ATOMICA's interface embedding already separates binders at a heme site, the program's framing changes |
| 2 | **Read RNAPro's gating module.** NVIDIA × Das Lab fused a *frozen* RibonanzaNet2 into a Protenix Pairformer by learned gating — i.e. someone shipped the frozen-encoder→co-folder adapter we intend to design | hours | The adapter design, for free, from people who trained it |
| 3 | **Extract Boltz-2 `s` and `z` for RNA systems** with the shipped `--write_embeddings` flag. `z` is joint by construction — protein × ligand × pose — which is exactly what survives the fingerprint critique. **Verified in source on 2026-09-20:** the flag is declared at `src/boltz/main.py:1038`; `src/boltz/data/write/writer.py:250` writes `embeddings_<id>.npz` when `"s"` and `"z"` are in the prediction; and `src/boltz/data/parse/schema.py:1021` accepts `entity_type` in `{protein, dna, rna, ligand}`. **So RNA + ligand + embeddings works today with no patch** | days | Whether `z` behaves on RNA as it does on protein. **Nobody has published this** |
| 4 | **Zero-shot transfer test on the 48-complex drug-like RNA set + 69 RNAmigos2 held-out clusters.** Freeze the protein-trained space; rank RNA ligands | days | Whether the shared space transfers across the protein/RNA boundary at all. A clean positive *or* negative, cheap either way |
| 5 | **Embed RNAcentral subsets with RiNALMo-giga + RNA-FM + RibonanzaNet** and check whether the three spaces agree (CKA / Procrustes) | ~3 GPU-days | Whether RNA LMs have converged. Per §1, the benchmarks say they have not, so expect disagreement — and disagreement is a usable ensemble signal (cf. this repo's FINDING 011, where cross-engine agreement was the only selector that beat random) |
| 6 | **Ribonanza as the RNA pretraining signal.** 2 M experimentally probed sequences is the one RNA data layer that is *larger* than its protein analogue | weeks | Whether grounding in measured chemistry beats grounding in co-evolution. RibonanzaNet is 11 M params — there is obvious headroom |

**What we will not do:** build an RNA-native affinity regressor (§6: n≈1,480, imbalanced), stage
a 2 TB MSA database (§3), or depend on AlphaFold3, HelixFold3, Chai-2 or Inforna (licence or
availability, §8).

**Dimension bookkeeping for the shared projection.** 1280 (RiNALMo-giga, AIDO.RNA-650M,
mRNA-FM, ProtRNA) · 2048 (AIDO.RNA-1.6B, RIBOSPAN) · 768 (ERNIE-RNA, RNAErnie, RNA-MSM, CaLM,
3UTRBERT) · 640 (RNA-FM) · 512/128 (DRfold2) · **384 single / 128 pair (RhoFold+, Boltz-1/1x/2,
AF3 family — one adapter covers all three)** · 256/64 (RibonanzaNet). RiNALMo-giga's 1280 matches
ESM-2-650M exactly, so the RNA↔protein projection head can be made symmetric with no reshaping.

---

## 8. BLOCKED — needs the user

Ordered by how much it costs us.

| # | Item | What is blocked | What we need from Amit |
|---|---|---|---|
| 1 | **R-BIND 2.0 data** | The canonical curated RNA-binder set (188 ligands). **`rbind.chem.duke.edu` no longer resolves** (Hargrove lab moved Duke→Toronto). The data now exists only as **Supporting Information on a paywalled ACS article** ([ACS Chem. Biol. 2022, 10.1021/acschembio.2c00224](https://doi.org/10.1021/acschembio.2c00224)). The distinct-target count and the affinity-vs-binder-only split are *only* in that spreadsheet | **Institutional ACS access** to download the SI XLSX — or a direct request to the Hargrove lab |
| 2 | **Inforna 2.0** | 1,936 RNA-motif–small-molecule interactions, 244 unique compounds — the **second genuinely independent** body of RNA-binding data after ROBIN. **Requires a signed academic licensing agreement** (SCRPS-DisneyLabAdmin@mail.ufl.edu). No bulk download exists even after licensing | A decision on whether to sign an academic licence, and who signs it |
| 3 | **PDBbind / PDBbind+** | The 101 RNA–ligand pairs with structure *and* affinity — the direct apples-to-apples number. Site requires **academic registration**; `pdbbind.org.cn` returned a socket hang-up on 2026-09-20 | Register an academic account at [pdbbind-plus.org.cn](http://www.pdbbind-plus.org.cn), or confirm we should rely on RNAmigos2/HARIBOSS instead |
| 4 | **OpenFold3 weights** | Apache-2.0 licence but the **HuggingFace repo is click-through gated** — breaks unattended download to the cluster | One manual click on [huggingface.co/openfold/openfold3](https://huggingface.co/openfold/openfold3) with the team HF account, then `HF_TOKEN` on the cluster |
| 5 | **RIBOSPAN-10K** | The only RNA LM with a **10,240-nt context** — the one model that could embed a whole mRNA. The long-context checkpoints are **"access upon request"**, and the licence is **non-commercial** | An access request to SII-GAIR-NLP, plus a decision on whether a non-commercial licence is acceptable |
| 6 | **AlphaFold3 parameters** | Non-commercial, **by Google form, no redistribution, no derived-weight sharing** | Probably a firm *no* — the licence forbids exactly what a shared embedding model implies. Confirm we exclude it |
| 7 | **CodonBERT** | Paper is behind **Genome Research SSO**; weights sit on a **Sanofi CDN under a bespoke, non-standard model licence** separate from the code licence. Architecture specs in §1 are therefore unverified | Institutional access to *Genome Research*, and a read of the model licence before we use the weights |
| 8 | **MultiMolecule AGPL-3.0** | The most convenient dependency in this document — ~25 RNA models behind one HF API, no flash-attn, no Paddle — **redistributes every checkpoint under AGPL-3.0**, which is *more restrictive than the originals* (RNA-FM MIT, RiNALMo CC-BY-4.0, RNA-MSM MIT) | A legal call: AGPL is fine for internal research; if anything is ever served or shipped, we must pull original weights from Zenodo/GitHub instead. **Same hazard applies to NuFold (GPL-3.0)** |
| 9 | **POSTAR3 reachability** | 1,499 CLIP-seq datasets / 348 RBPs / 7 species / ~50 M binding sites — the broadest and **only multi-species** CLIP compendium. `http://postar.ncrnalab.org` is HTTP-only and redirects to a bare IP (`111.198.139.65`) that **refused HTTPS** on 2026-09-20 | Someone to open it in a browser and confirm it is alive; if it is dead, ENCORI is the fallback (but is largely the same reprocessed GEO data) |
| 10 | **Paywalled papers hit during this survey** | ACS (both R-BIND papers), Wiley (ROBIN/*Angew.*), RSC, PNAS direct, *J. Mol. Biol.* (R-SIM), *Proteins* (CASP16 assessments — PMC has abstract only), *Genome Research* (CodonBERT) | Institutional access, or accept the PMC/preprint mirrors already used above |
| 11 | **Web budget** | This session exhausted its **200-call WebSearch budget**. Several cells remain marked *unverified*: trRosettaRNA2 parameter count and rep dims, exact param counts for AF3/Boltz/Chai/OpenFold3, RNAPro and OpenFold3 download sizes, RiboSphere model and codebook size, RNAsolo's current release number, liveness of the legacy aptamer databases, and whether RNAdecoyDB's data is actually downloadable | Raise `CLAUDE_CODE_MAX_WEB_SEARCHES_PER_SESSION`, or accept the flagged gaps |

### Not blocked, but decide before acquiring

- **Kaggle account** is needed for the original RibonanzaNet / RibonanzaNet2 weights. The
  MultiMolecule mirror avoids this but carries AGPL.
- **`mamba-ssm`, `flash-attn`, `vortex`/transformer-engine, DGL, PaddlePaddle, JAX** are each
  disqualifying under the stated constraint. That single rule removes Orthrus (as shipped), DGRNA,
  HydraRNA, CodonMamba, Evo 2, RoseTTAFold2NA, RFAA, RNAFlow, HelixFold3 and AlphaFold3 — **ten
  models, on engineering grounds alone.** If any of them is scientifically necessary, the
  constraint is what needs revisiting, not the model list.
