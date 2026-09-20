# RNA — models and data for the shared embedding space

**Author:** RNA lead, wave 1. **Written 2026-09-20.** *Status: sections 1–2 and 5–6 complete;
sections 3–4 merged from parallel reconnaissance.*

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
2,375 RNA-only**), and **RNA–small-molecule affinity is essentially absent**. See §4 for the
honest tally; the short version is that the closest RNA analogue of PDBbind contains a few
hundred non-redundant complexes and low thousands of measured affinities, against PDBbind's
~20 k and BindingDB's ~2.8 M.

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
| **RNA-MSM** | — | NAR 2024 | [yikunpku/RNA-MSM](https://github.com/yikunpku/RNA-MSM) · [multimolecule/rnamsm](https://huggingface.co/multimolecule/rnamsm) | yes | **MIT** | 96 M | 768 | 1024 | **MSA** | 0.39 *(est.)* | yes | ⚠️ **needs an MSA from RNAcmap3 (blastn vs `nt` + Infernal)** — GPU nodes have no internet; MSAs must be precomputed on a login node | Only evolutionary-coupling RNA model; best SS scores in two independent benchmarks. Treat the MSA pipeline as its own multi-day task |
| **RNABERT** | — | NARGAB 2022 | [multimolecule/rnabert](https://huggingface.co/multimolecule/rnabert) | yes | AGPL | 0.5 M | 120 | 440 | nt | <0.01 | yes | none | Baseline / sanity floor only. **Tiny.** |
| **UTR-LM (te_el, mrl)** | — | Nat.Mach.Int. 2024 | [multimolecule/utrlm-te_el](https://huggingface.co/multimolecule/utrlm-te_el) | yes | AGPL-3.0+ | **1.21 M** | **128** | 1022 | 5′UTR nt | <0.01 | yes | none | ⚠️ **A 1.2 M-parameter, 6-layer, 128-d model.** Do not call this a foundation model. Useful only as a 5′UTR-specific feature |
| **3UTRBERT** | 3/4/5/6-mer | 2024 | [multimolecule/utrbert-3mer](https://huggingface.co/multimolecule/utrbert-3mer) | yes | AGPL | ~86 M | 768 | 512 | 3′UTR k-mers | 0.35 *(est.)* | yes | none | 3′UTR arm; k-mer tokenisation makes alignment to nt-level spaces awkward |
| **SpliceBERT** | 510 / full | 2024 | [multimolecule/splicebert](https://huggingface.co/multimolecule/splicebert) | yes | AGPL | 19.4 M | 512 | 1024 | pre-mRNA nt | 0.08 | yes | none | Splice-context features; strongest zero-shot fitness model in the 2026 benchmark |
| **CodonBERT** | — | Genome Res. 2024 | [Sanofi-Public/CodonBERT](https://github.com/Sanofi-Public/CodonBERT) | **weights behind a Sanofi CDN link with a separate, non-standard model licence** | code: BSD-ish; **weights: bespoke** | ~87 M *(unverified — paper paywalled)* | 768 *(unverified)* | ~1024 codons | codons | ~0.35 *(est.)* | yes | poetry env | mRNA/CDS embeddings. ⚠️ licence and paper both gated — see BLOCKED |
| **CaLM** | — | Nat.Mach.Int. 2024 | [oxpig/CaLM](https://github.com/oxpig/CaLM) · [multimolecule/calm](https://huggingface.co/multimolecule/calm) | yes | AGPL (MM) | 86 M | 768 | 1024 codons | codons | 0.35 *(est.)* | yes | none | Codon LM trained on coding sequences — the natural bridge from RNA space to *protein* space |
| **ProtRNA** | — | Cell Syst. 2025 | [roxie-zhang/ProtRNA](https://github.com/roxie-zhang/ProtRNA) | yes | check repo | 651 M | 1280 | **512** | nt | ~2.6 *(est.)* | yes | TensorFlow in the original release | **Transfer-learned from ESM-2.** Directly relevant: it is literally a protein LM re-tuned on RNA, i.e. an existing attempt at the shared space |
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
| **Boltz-2** | 2025-06 | [jwohlwend/boltz](https://github.com/jwohlwend/boltz) | yes, **MIT code + MIT weights** | unpublished | **64-layer Pairformer, single 384, pair 128** | affinity head caps ligands at 128 atoms | YAML, `entity_type: rna` | auto | yes (cuEquivariance optional) | optional **CORRECTION:** my structure agent reported "no `--write_embeddings` flag"; the sibling protein/ligand agent read the source and found it **is shipped** (`boltz/main.py:1038`), writing `embeddings_<id>.npz` with `s [n,384]`, `z [n,n,128]`. Use the flag, not a patch. The deeper per-layer recipe (concat layers 16/32/48/64) from [hsjang0/boltz-as-FM](https://github.com/hsjang0/boltz-as-FM) (MIT) still needs a patch. **RNA transfer untested — that paper is protein/ligand only** |
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
| **ATOMICA** | 2025, bioRxiv v4 | [mims-harvard/ATOMICA](https://github.com/mims-harvard/ATOMICA) · [ada-f/ATOMICA](https://huggingface.co/ada-f/ATOMICA) | yes, **MIT**, ungated | unpublished | **8 named reps at atom / block / interface / graph level**, selectable via `--guidance` | interface-local | PDB files -> Parquet | — | yes; explicitly needs **no** host CUDA toolkit, no torch-scatter/cluster | none | Zitnik lab. One embedding space over **protein-RNA, protein-DNA and nucleic-acid-small-molecule interfaces**, trained on 2,037,972 complexes — **and it already ships HEM/heme finetunes.** Closest existing artifact to our target on the interaction side. Needs structures, so pair it with a folder |
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
