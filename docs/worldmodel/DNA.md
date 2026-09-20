# DNA / genomics — frozen embedding sources and structure data at scale

**Researched 2026-09-20.** Domain report for the shared-embedding world model
(`docs/worldmodel/README.md`). Two jobs: (1) DNA foundation models we can freeze and
embed with, (2) DNA structure and protein–DNA interaction data we can actually get.

Target environment, and the thing every recommendation is filtered through:
**one contended H200 (80–140 GB), 298 TB of `/scratch`, no local disk, torch 2.5.1, no
`cuequivariance`.** Custom CUDA kernels pinned to a torch version are a real cost here,
not a footnote.

---

## 0. The short version

- **Evo 2 7B** is the default DNA encoder: Apache-2.0, ungated, 13.8 GB, 1 Mb context,
  single-nucleotide, 4096-d, and the only large model that runs bf16 on one GPU *without*
  Transformer Engine. The 40B is 82.2 GB and does **not** fit an 80 GB H200 in bf16.
- **AlphaGenome's weights are downloadable now — it is no longer API-only** — with a
  validated PyTorch port whose `encode()` returns (B, 1024, 3072) at 128 bp over 1 Mb in
  40.8 GB on one H200, no exotic kernels. **But its terms forbid using its outputs to train
  other ML models**, so it can only ever be an external comparator, never a component of
  the shared space. That is the most consequential licence fact in this report.
- **The function tower is therefore Borzoi (CC-BY-4.0, ~1920-d @32 bp) or Enformer
  (CC-BY-4.0, 3072-d @128 bp)**, both plain torch, both permissively licensed.
- **NTv3 (650M, 1 Mb, single-nucleotide, 1536-d)** is the best-matched encoder in the table
  and is **gated behind a non-commercial licence** — see *BLOCKED*.
- **JEPA-DNA (NVIDIA, 2026)** already ran the exact experiment this program proposes, on
  the DNA axis, across five backbones, and released code and checkpoints. Read it before
  writing any JEPA head. The checkpoints are demo-scale; the *recipe* is the value.
- **The structural half is tiny and that is the load-bearing finding.** 10,733 protein–DNA
  PDB entries → **3,027 non-redundant clusters at 30% identity**, total download **4.32 GB**.
  Everything else in section 3 — NAKB, DNAproDB, BioLiP2, twelve PDNA benchmarks — is the
  same PDB coordinates re-wrapped. Against 8.8 Tbp of sequence corpus that is a ratio near
  10⁻⁹.
- **A complete, useful DNA corpus is under 200 GB**, against a 298 TB scratch and the
  README's 10 TB cap — ENCODE's four narrowPeak sets are 131 GB, the SCREEN cCRE registry
  is 129 MB, and every motif bank in existence is under 7 MB. Data volume is not the
  constraint on this arm; supervision density is. **The only million-scale human
  sequence→function supervision that exists openly is SuRE-SNP (5.9 M SNPs, 119 GB).**
- **Benchmarks overlap more than they admit.** GUE and the NT suite share roughly half
  their datasets; BEND, DART-Eval and LRB are largely ENCODE again. Scoring well on three
  of them is one result reported three times.

---

## 1. DNA / genomic foundation models

### 1.1 Master table

Embedding dim = width of the residual stream you would pool. "1×H200" = fits and runs
inference on one 80 GB card.

| # | Model | Ver | Yr | URL | Weights open? | License | Params | Emb dim | Max ctx | Input | Size | 1×H200 | Exotic deps |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **Evo 2** | `evo2_7b` | 2026 | [HF](https://huggingface.co/arcinstitute/evo2_7b) · [GH](https://github.com/arcinstitute/evo2) · [Nature](https://www.nature.com/articles/s41586-026-10176-5) | yes, ungated | Apache-2.0 | 7B, 32 L, 32 H | **4096** (confirmed, repo `evo2-7b-1m.yml`) | **1 Mb** | raw nt, byte-level | **13.8 GB** | **yes**, bf16, no TE | `vortex`; flash-attn 2.8; torch 2.6/2.7 rec. |
| 2 | Evo 2 | `evo2_40b` | 2026 | [HF](https://huggingface.co/arcinstitute/evo2_40b) | yes | Apache-2.0 | 40B, 50 L, 64 H | **8192** (confirmed, `evo2-40b-1m.yml`) | 1 Mb | raw nt | **82.2 GB** (2 parts) | **no in bf16**; FP8 only (~41 GB) | **requires Transformer Engine 2.3 + FP8** |
| 3 | Evo 2 | `evo2_20b` | 2026 | [HF](https://huggingface.co/arcinstitute/evo2_20b) | yes | Apache-2.0 | 20B | ~6144 (inferred) | 1 Mb | raw nt | ~41 GB | yes w/ TE FP8 | TE required |
| 4 | Evo 2 | `evo2_7b_base` / `_262k` / `evo2_1b_base` | 2026 | [collection](https://huggingface.co/collections/arcinstitute/evo) | yes | Apache-2.0 | 7B / 7B / 1B | — | 8k / 262k / 8k | raw nt | 13.8 / 13.8 / ~2 GB | yes (1b_base wants TE) | as above |
| 5 | **Evo 1** | `evo-1-131k-base` | 2024 | [HF](https://huggingface.co/togethercomputer/evo-1-131k-base) · [GH](https://github.com/evo-design/evo) | yes | Apache-2.0 | 7B | 4096 | 131k | raw nt | ~14 GB | yes | StripedHyena-1 kernels |
| 6 | **NTv3** | `NTv3_650M_pre` / `_post` | 2026 | [HF](https://huggingface.co/InstaDeepAI/NTv3_650M_pre) · [paper](https://instadeep.com/wp-content/uploads/2025/12/NT_v3.pdf) | **GATED** | **custom, non-commercial** | 650M, 12L/24H | **1536** | **1 Mb**, single-nt | raw nt, len % 128 | ~2.6 GB | yes | none — `transformers>=4.55`, `trust_remote_code` |
| 7 | NTv3 | `NTv3_100M_*`, `NTv3_8M_*`, `_131kb`, `_8kb`, `NTv3_generative`, `NTv3_5downsample_*` | 2026 | [collection](https://huggingface.co/collections/InstaDeepAI/nucleotide-transformer-v3) | **GATED** | non-commercial | 7.7M / 100M / 600M | 1536 / ~768 | 8kb–1 Mb | raw nt | 30 MB–2.6 GB | yes | none |
| 8 | **NT v2** | `nucleotide-transformer-v2-500m-multi-species` | 2024 | [HF](https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-500m-multi-species) · [Nat Methods](https://www.nature.com/articles/s41592-024-02523-z) | yes, ungated | **CC-BY-NC-SA-4.0** | 500M, 29L | **1024** | 2048 tok = **12,288 bp** | 6-mer, vocab 4107 | ~2.0 GB fp32 | yes | none (ESM arch) |
| 9 | NT v2 | `-v2-50m/100m/250m-multi-species` | 2024 | [collection](https://huggingface.co/collections/InstaDeepAI/nucleotide-transformer) | yes | CC-BY-NC-SA-4.0 | 50–250M | 512–768 | 12,288 bp | 6-mer | 0.2–1.0 GB | yes | none |
| 10 | NT v1 | `-2.5b-multi-species`, `-2.5b-1000g`, `-500m-human-ref` | 2023 | [HF](https://huggingface.co/InstaDeepAI/nucleotide-transformer-2.5b-multi-species) | yes | CC-BY-NC-SA-4.0 | 2.5B | 2560 | 1000 tok = 6 kb | 6-mer | ~10 GB | yes | none |
| 11 | **HyenaDNA** | `hyenadna-large-1m-seqlen-hf` | 2023 | [HF](https://huggingface.co/LongSafari/hyenadna-large-1m-seqlen-hf) · [GH](https://github.com/HazyResearch/hyena-dna) | yes | **BSD-3** | 54.6M reported (paper says 6.6M — see note) | **256** (`d_model`) | **1,000,002** | single nt, vocab 12 | ~220 MB | trivially | **`-hf` variants are pure torch**; original repo wants flash-attn |
| 12 | HyenaDNA | `-tiny-1k`, `-small-32k`, `-medium-160k`, `-medium-450k` | 2023 | [collection](https://huggingface.co/LongSafari) | yes | BSD-3 | 0.4–28M | 128–256 | 1k–450k | single nt | <120 MB | trivially | none (`-hf`) |
| 13 | **Caduceus** | `caduceus-ph_seqlen-131k_d_model-256_n_layer-16` | 2024 | [HF](https://huggingface.co/kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16) · [arXiv](https://arxiv.org/abs/2403.03234) | yes | Apache-2.0 | **7.7M** | **256** | 131k | single nt, vocab 16 | **30.9 MB** | trivially | ⚠️ **`mamba-ssm` + `causal-conv1d` custom CUDA, torch-version-pinned wheels** |
| 14 | Caduceus | `caduceus-ps_…-131k_…` (RC-equivariant) | 2024 | [HF](https://huggingface.co/kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16) | yes | Apache-2.0 | ~1.9M (weight-tied) | 256 | 131k | single nt | ~8 MB | trivially | same mamba deps |
| 15 | **PlantCaduceus** | `PlantCaduceus_l32` (also l20/l24/l28) | 2025 | [HF](https://huggingface.co/kuleshov-group/PlantCaduceus_l32) · [PNAS](https://doi.org/10.1073/pnas.2421738122) | yes | Apache-2.0 | 225M | 1024 | 512 nt | single nt | ~0.9 GB | yes | same mamba deps |
| 16 | **DNABERT-2** | `DNABERT-2-117M` | 2024 | [HF](https://huggingface.co/zhihan1996/DNABERT-2-117M) · [GH](https://github.com/MAGICS-LAB/DNABERT_2) | yes | **code Apache-2.0; the HF weights card carries NO licence tag** ⚠️ | 117M, 12L | **768** | 512 BPE tok (ALiBi extrapolates) ≈ 2–3 kb | BPE, vocab 4096 | ~0.47 GB | trivially | ⚠️ **ships a Triton flash-attn that breaks on modern Triton**; use the drop-in [`quietflamingo/dnabert2-no-flashattention`](https://huggingface.co/quietflamingo/dnabert2-no-flashattention) (bit-identical output) |
| 17 | DNABERT-S | `DNABERT-S` | 2024 | [HF](https://huggingface.co/zhihan1996/DNABERT-S) | yes | Apache-2.0 | 117M | 768 | as DNABERT-2 | BPE | ~0.47 GB | trivially | same Triton caveat |
| 18 | DNABERT (v1) | `DNA_bert_6` | 2021 | [HF](https://huggingface.co/zhihan1996/DNA_bert_6) | yes | Apache-2.0 | 86M | 768 | 512 tok = 512 bp | 6-mer overlap | ~0.4 GB | trivially | none — **superseded, skip** |
| 19 | **GENA-LM** | `gena-lm-bert-large-t2t` | 2023/25 | [HF](https://huggingface.co/AIRI-Institute/gena-lm-bert-large-t2t) · [GH](https://github.com/AIRI-Institute/GENA_LM) | yes | open (per-repo; code Apache-2.0) | 336M, 24L | **1024** | 512 BPE tok ≈ **4.5 kb** | BPE, vocab 32k | ~1.4 GB | trivially | none |
| 20 | GENA-LM | `gena-lm-bigbird-base-t2t` | 2023 | [HF](https://huggingface.co/AIRI-Institute/gena-lm-bigbird-base-t2t) | yes | open | 110M | 768 | 4096 tok ≈ **36 kb** | BPE | ~0.45 GB | trivially | BigBird sparse attn (pure torch) |
| 21 | **ModernGENA** | `moderngena-large` / `-base` | 2026 | [HF](https://huggingface.co/AIRI-Institute/moderngena-large) · [OpenReview](https://openreview.net/forum?id=RB8LzQ6MIX) | **yes — weights confirmed on HF** | not stated on card ⚠️ | **377M** / 100M | ~1024 (ModernBERT-large) | long (ModernBERT rotary + local/global) | BPE 32k (GENA-LM tokenizer) | ~1.5 GB / ~0.4 GB | trivially | **none** — flash-attn optional, `torch.compile` supported ✅ |
| 21b | **GENAtator** (ModernGENA/Caduceus/GENA-RMT annotation heads) | 2026 | [HF org](https://huggingface.co/AIRI-Institute) | yes | not stated ⚠️ | 14M–0.4B | — | — | — | small | trivially | none |
| 22 | **JEPA-DNA** | `NV-JEPA-DNA-{DNABERT2,NTv3,HyenaDNA}` | 2026 | [GH](https://github.com/NVIDIA-Digital-Bio/JEPA-DNA) · [arXiv 2602.17162](https://arxiv.org/abs/2602.17162) · [HF](https://huggingface.co/collections/nvidia/jepa-dna) | yes, ungated | **code Apache-2.0; the NTv3 checkpoint card says NVIDIA non-commercial** | 100M (NTv3 var.) | inherits | **8,192 bp** | raw nt | small | yes | flash-attn *optional* (DNABERT-2 path) |
| 23 | **gLM2** | `tattabio/gLM2_650M` (+ `_150M`) | 2024/25 | [HF](https://huggingface.co/tattabio/gLM2_650M) · [ICLR'25 OMG](https://proceedings.iclr.cc/paper_files/paper/2025/file/dbb40375daa3c39f1c098b12608aed12-Paper-Conference.pdf) | yes | **Apache-2.0** | 650M | ~1280 | 4096 tok | **mixed: AA (upper) + nt (lower)** | ~2.6 GB | yes | `trust_remote_code`, bf16 |
| 24 | gLM (v1) | Hwang et al. | 2024 | [Nat Commun](https://www.nature.com/articles/s41467-024-46947-9) | yes | open | ~1B | 1280 | 30 genes | **genes-as-tokens (ESM2 embs)** | ~4 GB | yes | none |
| 25 | **GenomeOcean** | `GenomeOcean-4B` (+100M/500M) | 2025 | [HF](https://huggingface.co/pGenomeOcean/GenomeOcean-4B) · [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11838515/) | yes | **BSD** | 4B | ~2560 | 10,240 tok ≈ **50 kbp** | BPE, vocab 4096 | ~8 GB bf16 | yes | **flash-attn-2 required in the quickstart** |
| 26 | **METAGENE-1** | `metagene-ai/METAGENE-1` | 2025 | [HF](https://huggingface.co/metagene-ai/METAGENE-1) | yes | Apache-2.0 | 7B | 4096 | **512 tok only** | BPE | ~13 GB | yes | none (Llama-style) |
| 27 | **GROVER** | `PoetschLab/GROVER` | 2024 | [HF](https://huggingface.co/PoetschLab/GROVER) · [Nat Mach Intell](https://www.nature.com/articles/s42256-024-00872-0) | yes | open | ~90M (BERT-base) | 768 | 512 BPE tok | BPE | ~0.4 GB | trivially | none. ⚠️ re-tokenises unstably below 50 nt — pad 100 nt |
| 28 | **SegmentNT** | `InstaDeepAI/segment_nt` | 2024/25 | [HF](https://huggingface.co/InstaDeepAI/segment_nt) | yes, ungated | **CC-BY-NC-SA-4.0** | 562M (NT + 53M U-Net head) | 1024 | 30 kb (→50 kb) | 6-mer, vocab 4105 | ~2.3 GB | yes | `transformers` from source |
| 29 | ChatNT / Agro-NT / Isoformer / BulkRNABert / sCT / MOJO | — | 2024–25 | [GH](https://github.com/instadeepai/nucleotide-transformer) | yes | **CC-BY-NC-SA-4.0** | varies | varies | varies | varies | varies | yes | ChatNT is **JAX/Haiku** |

### 1.2 Notes, and where to be sceptical

- **Evo 2 is the one to start with, but mind the torch pin.** The repo recommends torch
  2.6/2.7 and flash-attn 2.8.0.post2. Our torch-2.5.1 environment will probably load the
  7B in plain bf16 (the README explicitly says 7B needs no Transformer Engine) but the
  `vortex` fast path may not build. Budget a day for this. Embeddings come out of an
  *intermediate* block (`blocks.28.mlp.l3`) — the repo says intermediate beats final,
  which matters because the obvious `last_hidden_state` pool is the wrong choice here.
- **`evo2_7b`'s `config.json` on HF contains three keys** — `{"name","architecture",
  "_name_or_path"}`, no hidden size, no layer widths. The real shapes are in the GitHub
  repo's YAML, not on the Hub: `evo2-7b-1m.yml` gives `hidden_size 4096, num_layers 32,
  num_attention_heads 32, vocab_size 512, max_seqlen 1048576`; `evo2-40b-1m.yml` gives
  `hidden_size 8192, num_layers 50, num_attention_heads 64`. The 20B width is still
  inferred. Anyone pulling Evo 2 config from the Hub alone will get nothing.
- **HyenaDNA's parameter count is inconsistent.** The HF card for
  `hyenadna-large-1m-seqlen-hf` reports 54.6M; the paper reports 6.6M for the 1M-context
  model. The gap is the materialised Hyena positional/filter tensors at `max_seq_len ≈
  10⁶`. Neither number is wrong; they measure different things. Don't quote either as
  "model size" without saying which.
- **Caduceus is tiny and good and will fight our environment.** 30.9 MB for 131k context
  is the best size/context ratio in the table, but it needs `mamba-ssm` and
  `causal-conv1d`, whose prebuilt wheels are pinned to specific torch/CUDA combos.
  On torch 2.5.1 this is a source build. There is a slow pure-PyTorch Mamba fallback;
  expect a large slowdown. This is the single biggest "exotic dep" risk in the table.
- **DNABERT-2's weights are, strictly, unlicensed.** The `MAGICS-LAB/DNABERT_2` GitHub
  repo is Apache-2.0, but the HuggingFace model card for `zhihan1996/DNABERT-2-117M`
  displays no licence tag at all. Same for GENA-LM's per-repo tags, which I could not
  confirm. For a research world model this is a non-issue; for anything downstream it is a
  question to ask before it matters.
- **DNABERT-2's shipped Triton kernel is broken on current Triton** — a well-documented,
  years-old open issue. The community fork removes flash-attention and is reported to
  produce identical output. Use the fork; do not debug the kernel.
- **The Nucleotide Transformer family is CC-BY-NC-SA-4.0 across v1/v2/SegmentNT/ChatNT,
  and NTv3 adds gating on top.** For a research world model this is usable; it forecloses
  anything commercial and the *share-alike* clause arguably reaches derived embeddings.
  Flagged below.
- **Redundancy to prune.** Evo 1 is strictly superseded by Evo 2 (same lab, same
  architecture family, Evo 1 is prokaryote-only). DNABERT v1 is superseded by DNABERT-2.
  NT v1 2.5B is superseded by NT v2 500M on the NT benchmark at 1/5 the size. Carrying all
  three costs 25 GB and buys nothing.
- **METAGENE-1 is 7B parameters with a 512-token context.** That ratio is the giveaway:
  it is a pathogen-detection classifier backbone, not a long-range genomic encoder. Skip
  it for our purposes despite the attractive size and licence.
- **Benchmarks disagree with the marketing.** [Nat Commun 2025,
  s41467-025-65823-8](https://www.nature.com/articles/s41467-025-65823-8) benchmarked
  DNABERT-2, NT v2, HyenaDNA, Caduceus-Ph and GROVER and found general-purpose DNA FMs
  *competitive on pathogenic-variant identification but worse than specialised models at
  gene-expression prediction and causal-QTL identification*. Treat any claim that a DNA FM
  "matches Enformer" as unreplicated until it is on a difficulty-matched split — the same
  trap FINDING 022 names in this repo.
- **ModernGENA is the quiet best-value entry.** A ModernBERT backbone retrained on
  **443 vertebrate genome assemblies / 353,574,093,776 bp**, sampled in a
  `[-16 kbp, +8 kbp]` window around transcription start sites — i.e. deliberately
  regulatory-biased pretraining, unlike every other model in this table. 377M params,
  pure `transformers`, flash-attention *optional*, `torch.compile` supported, no mamba and
  no Triton. On our torch-2.5.1 box it is the model most likely to just work on the first
  try. The ICLR-2026 workshop paper claims it "ranks among the top-performing models" on
  the NT benchmark; that is a workshop poster, not a replicated result.
- **Benchmark suites worth using, in rough order:** GUE (DNABERT-2's 28 tasks), the NT
  18-task suite, [BEND](https://github.com/frederikkemarin/BEND), Genomic Benchmarks,
  GenBench, [NABench 2025](https://arxiv.org/abs/2511.02888), DART-Eval, and GFMBench
  (what JEPA-DNA scores on). They overlap heavily — GUE and the NT suite share tasks.

### 1.3 The finding that matters most

**JEPA-DNA (NVIDIA Digital Bio, arXiv 2602.17162; v1 Feb 2026, final Aug 2026).** A
model-agnostic *continual-training* framework that adds a joint-embedding predictive
objective on top of existing genomic FMs — "supervises global sequence embeddings in a
latent space to predict functional representations of masked segments" — alongside the
usual generative loss. Backbones supported: DNABERT-2, NTv3 (8M and 100M), HyenaDNA,
Caduceus, DNABERT. Code Apache-2.0; checkpoints at `nvidia/NV-JEPA-DNA-DNABERT2`,
`nvidia/NV-JEPA-DNA-NTv3`, `nvidia/NV-JEPA-DNA-HyenaDNA`; claims SOTA on 17 genomic tasks.

This is the DNA-axis version of this repo's ask 1, already built, already open, with the
compute already spent. Two consequences: we should not re-derive it, and we should treat
its published numbers as the bar a home-grown DNA JEPA head has to clear. It also gives
us a free ablation — JEPA'd vs. base backbone from the *same* lab — which is exactly the
"broad vs. narrow prior" comparison the README pre-registers, for free.

**Caveats, and they matter.** (a) The released checkpoints are *small*:
`NV-JEPA-DNA-HyenaDNA` is **600k parameters** on a HyenaDNA-16k backbone with an 8,192 bp
context; `NV-JEPA-DNA-NTv3` is **100M** on NTv3-100M, also 8,192 bp. These are
demonstration artefacts, not a drop-in replacement for Evo 2. (b) Every checkpoint card
states an **NVIDIA non-commercial licence** even though the *code* is Apache-2.0 — one
card adds "subject to internal legal review", which is not a licence status anyone should
build on. (c) It wants torch 2.6.0. (d) 18 authors, an arXiv preprint with a six-month
revision window, no venue named in what I could fetch, SOTA claim not independently
replicated.

So: take the **recipe and the ablation design**, not the weights. The value is that
somebody has already measured what a JEPA objective adds on top of a genomic FM, across
five backbones, on 17 tasks — which is the exact experiment this repo's README
pre-registers, and it means our DNA arm starts from a published baseline instead of zero.

---

## 2. Regulatory-genomics models (supervised track predictors)

These are a **different species of model** from section 1 and, for anything to do with
what a genomic locus *does*, they are the better embedding source. DART-Eval (NeurIPS
2024, Kundaje lab) is blunt about it: DNA language models "do not offer compelling gains
over alternative baseline models for most tasks, while requiring significantly more
computational resources." Section 1 gives evolutionary constraint; section 2 gives
function.

### 2.1 Master table

| # | Model | Yr | URL | Weights open? | License | Gated | Params | Emb dim | Max ctx (bp) | Size | 1×H200 | Exotic deps |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | **AlphaGenome** (official, JAX) | 2025/26 | [HF `google/alphagenome-all-folds`](https://huggingface.co/google/alphagenome-all-folds) · [GH](https://github.com/google-deepmind/alphagenome_research) | **yes — weights are downloadable now, no longer API-only** | **AlphaGenome Model Terms — non-commercial only** | **yes, accept-terms** | **450M** | not exposed via the JAX API | **1 Mb (2²⁰)** | — | yes (≥H100; 34.6 GB @1 Mb) | **JAX + Haiku + JAXline** ❌ |
| 31 | **AlphaGenome-PyTorch** | 2026 | [GH](https://github.com/genomicsxai/alphagenome-pytorch) · [HF `gtca/alphagenome_pytorch`](https://huggingface.co/gtca/alphagenome_pytorch) | yes | code Apache-2.0; **weights inherit the NC terms** | port not HF-gated | 450M | **3072** @128 bp (1024 tokens) | 1 Mb | `model_all_folds.safetensors` + 4 folds | **yes — 40.8 GB peak @1 Mb, 260 ms/forward** | **none — plain torch** ✅ |
| 32 | **Enformer** | 2021 | [HF `EleutherAI/enformer-official-rough`](https://huggingface.co/EleutherAI/enformer-official-rough) · [GH](https://github.com/lucidrains/enformer-pytorch) | yes | **CC-BY-4.0** | no | ~250M | **3072** (896×3072 trunk) | **196,608** | ~1 GB | trivially | none |
| 33 | **Borzoi** (original) | 2023→2025 | [GH](https://github.com/calico/borzoi) | yes (GCS `.h5`) | code Apache-2.0 | no | ~0.2B | — | **524,288** (32 bp bins) | 4 folds | — | **TensorFlow 2.15 only** ❌ |
| 34 | **borzoi-pytorch** | 2024/25 | [GH](https://github.com/johahi/borzoi-pytorch) · `johahi/borzoi-replicate-{0..3}{,-mouse}` | yes | **CC-BY-4.0** | no | 0.2B | ~1920 trunk (`get_embs`) | 524 kb | F32 safetensors | yes | plain torch ✅ |
| 35 | **Flashzoi** | 2025 | `johahi/flashzoi-replicate-{0..3}` | yes | MIT / CC-BY | no | 0.2B | ~1920 | 524 kb | small | yes | **FlashAttention-2 + autocast required** ⚠️ |
| 36 | **Decima** | 2026 | [GH](https://github.com/Genentech/decima) · Zenodo 15092691 · `Genentech/decima-model` | yes | permissive (confirm) | no | ~0.2B (Borzoi trunk) | Borzoi trunk, avg-pooled | ~524 kb | 4 replicates | yes | `grelu` + torch |
| 37 | **scooby** | 2025 | [GH](https://github.com/gagneurlab/scooby) · `johahi/neurips-scooby` | yes | **MIT** | no | Borzoi + LoRA | needs a *precomputed cell embedding* as input | 524 kb | adapters only | yes (batch=1) | `snapatac2-scooby`, `peft` |
| 38 | **gReLU** (framework + zoo) | 2025 | [GH](https://github.com/Genentech/gReLU) · [HF zoo](https://huggingface.co/collections/Genentech/grelu-model-zoo) | yes | CC-BY | no | — | — | — | — | yes | torch-lightning ✅ |
| 39 | **ChromBPNet** | 2025 | [GH](https://github.com/kundajelab/chrombpnet) · `kundajelab/encode-chrombpnet-*` | yes | **MIT** | no | few M each | penultimate conv (not officially exposed) | **2114** | tiny each | trivially | **TensorFlow/Keras** ❌ (3rd-party torch port unvalidated) |
| 40 | **Sei** | 2022 | [GH](https://github.com/FunctionLab/sei-framework) · Zenodo 4906996 | yes | **academic / NC (Princeton OTL)** ❌ | no | — | 21,907 profiles → **40 sequence classes** | 4096 | — | trivially | PyTorch + **Selene** (stale) |
| 41 | **DeepSEA / Beluga / ExPecto** | 2015–18 | [GH](https://github.com/FunctionLab/ExPecto) | yes | **academic only, commercial prohibited** ❌ | no | small | 919 / 2002 | 1000 / 2000 | small | trivially | old torch + XGBoost |
| 42 | **Basenji / Basenji2** | 2018/20 | [GH](https://github.com/calico/basenji) | yes | Apache-2.0 | no | ~30M | 1536 | 131 kb | small | yes | **TensorFlow only** ❌ — superseded, skip |
| 43 | **Basset / DanQ** | 2016 | — | yes | permissive | no | ~4M | ~1000 | 600 / 1000 | tiny | yes | Torch7 / Keras-TF, **effectively dead** ❌ |
| 44 | **scBasset** | 2022 | [GH](https://github.com/calico/scBasset) | **no pretrained weights** | Apache-2.0 | no | small | **32** (bottleneck) | 1344 | — | yes | **TF/Keras** ❌ |
| 45 | **Orca** | 2022 | [GH](https://github.com/jzhoulab/orca) · Zenodo 6234936 | yes | **NC / academic** ❌ | no | — | — | **1 Mb / 32 Mb / 256 Mb** | 1.3 GB core + 34 GB data | yes | **custom Selene fork**, py3.9 ⚠️ |
| 46 | **Akita** | 2020 | [GH](https://github.com/calico/basenji) (`akita/`) | yes | Apache-2.0 | no | small | — | ~1 Mb | small | yes | **TensorFlow only** ❌ |
| 47 | **C.Origami** | 2023 | [GH](https://github.com/tanjimin/C.Origami) · Zenodo 7226561 | yes | not stated ⚠️ | no | — | — | **2 Mb** | — | yes | torch 1.12 pin ⚠️, pybigwig |
| 48 | **GET** | 2024 (Nature) | [GH](https://github.com/GET-Foundation/get_model) | yes (S3 + HF Space) | **CC-BY-NC-4.0** ❌ | no | — | — | — | — | yes | torch, bedtools, R/pcalg. **Input is ATAC peaks + motifs, not raw sequence** |
| 49 | **CREsted** | 2026 | [GH](https://github.com/aertslab/CREsted) | code yes, **no model zoo** | not stated ⚠️ | no | small | — | — | — | yes | **Keras 3** (torch backend OK; multi-GPU TF-only) |
| 50 | **Enformer-Celltyping** | 2024 | [GH](https://github.com/neurogenomics/EnformerCelltyping) | yes | **MIT** | no | ~Enformer | — | 196,608 | — | yes | **TensorFlow** ❌ |
| 51 | **regLM** | 2024 | [GH](https://github.com/Genentech/regLM) | weights unclear ⚠️ | MIT | no | small (HyenaDNA) | — | — | — | yes | **generative, not an encoder — skip** |
| 52 | **DeepPBS** | 2024 (Nat Methods) | [GH](https://github.com/timkartar/DeepPBS) · figshare 25678053 | yes | **BSD-3** ✅ | no | — | — | **structures (PDB/mmCIF), not sequence** | — | yes | torch-geometric 2.5 / torch 2.3; 3DNA/Curves preprocessing. **no `cuequivariance`** ✅ |
| 53 | **Geneformer** | V2, 2024 | [HF `ctheodoris/Geneformer`](https://huggingface.co/ctheodoris/Geneformer) | yes | **Apache-2.0** | no | 10M/104M/316M | 512–1024 | 2048/4096 **genes**, not bp | ≤1.3 GB | yes | HF transformers. **Not a DNA-sequence model** |
| 54 | **scGPT** | 2024 | [GH](https://github.com/bowang-lab/scGPT) | yes | MIT | no | ~50M | 512 | genes, not bp | small | yes | flash-attn optional. Same caveat |
| 55 | **DART-Eval** | 2024 | [GH](https://github.com/kundajelab/DART-Eval) | benchmark | — | no | — | — | — | — | — | **read before betting on any DNA LM** |

### 2.2 Notes

- **⚠️ AlphaGenome is comparator-only.** Its Model Terms say its outputs *"should not be
  used for the training of other machine learning models."* It is the best-matched encoder
  in this table and we cannot put it in the model. Use **Borzoi (CC-BY-4.0) or Enformer
  (CC-BY-4.0)** as the function tower instead; see §5.3.
- **The dimensions line up anyway.** Enformer 3072 @128 bp, Borzoi ~1920 @32 bp,
  AlphaGenome 3072 @128 bp. Enformer and Borzoi pool to a common 128-bp grid cleanly,
  which is what the shared space needs.
- **Port fidelity is not uniform, and this repo has already been burned by exactly this**
  (PXR lesson 1: validate the reference frame before trusting any cross-model number).
  AlphaGenome-PyTorch is the best validated — layer outputs match JAX to <1e-5 relative
  error, contact-map Pearson r = 0.9999. `borzoi-pytorch` claims exact parity with
  validation notebooks. **`enformer-pytorch` is the weak one**: its own README admits
  absolute error as high as 0.5 accumulating across layers, correlation preserved but not
  values. Enformer embeddings are correlation-valid, not value-valid — do not compare
  magnitudes across models.
- **Scepticism on AlphaGenome-PyTorch.** The port and its H200 benchmark both come from a
  single group publishing on a personal blog; ~161 stars, 9 open PRs, young. The 40.8 GB /
  260 ms numbers are theirs and are unreplicated. Time one forward pass before planning a
  corpus-scale extraction around them.
- **The TensorFlow/JAX friction list**, a real cost on our stack: original AlphaGenome
  (JAX/Haiku), original Borzoi (TF 2.15), **ChromBPNet (TF/Keras — and its 1,512-model
  ENCODE zoo across 408 biosamples is the thing we would most want)**, Basenji/Basenji2,
  Akita, scBasset, Enformer-Celltyping.
- **Fully-permissive fallback stack**, if the NC licences are unacceptable:
  Enformer (CC-BY-4.0) + borzoi-pytorch/Flashzoi (CC-BY-4.0/MIT) + ChromBPNet (MIT).
- **DeepPBS is the one entry bridging sections 2 and 3** — protein–DNA binding specificity
  predicted *from structure*, BSD-3, torch-geometric, no `cuequivariance`. It is the
  natural supervision signal linking a DNA embedding to a protein embedding.
- **Redundancy to prune here too.** Basenji2 → Borzoi (same lab, strict successor).
  Akita → AlphaGenome contact maps. DeepSEA/Beluga → Sei (same lab, Sei is the successor
  and subsumes the 919/2002 track sets). Carrying the ancestors buys nothing but
  TensorFlow dependencies.

### 2.3 One line each: how every model above would be used in the shared space

Numbers match the tables. "Encoder" = frozen, we pool its trunk and project into the
shared space; "channel" = an auxiliary vector concatenated to an encoder's output;
"supervision" = it produces targets, not embeddings.

| # | Model | Concrete use |
|---|---|---|
| 1 | Evo 2 7B | **Primary DNA encoder.** Mean-pool `blocks.28.mlp.l3` over a locus → 4096-d evolutionary-constraint vector; the DNA-side tower of the JEPA. |
| 2–4 | Evo 2 40B/20B/base | Ablation only — does the DNA tower improve with scale, at 3–6× the extraction cost? |
| 5 | Evo 1 | Skip. Strictly dominated by Evo 2. |
| 6–7 | NTv3 | Encoder if licensing clears: the only 1-Mb *single-nucleotide* 1536-d embedding; `_post` variants add 16k functional tracks as free supervision. |
| 8–10 | NT v1/v2 | Cheap 1024-d baseline encoder to sanity-check that Evo 2's extra 6.5B params buy anything. |
| 11–12 | HyenaDNA | Ultra-cheap long-context encoder for whole-gene-body sweeps where we cannot afford Evo 2 per window. |
| 13–14 | Caduceus | Reverse-complement-**equivariant** encoder — the only one whose embedding is invariant to strand, which matters because a TF site has no canonical strand. Worth the mamba build pain for that property alone. |
| 15 | PlantCaduceus | Skip for CYP3A4. Keep only if the world model is ever asked to cross into plant P450s. |
| 16–17 | DNABERT-2 / -S | 768-d baseline encoder; DNABERT-S is the contrastive variant, better for retrieval-style partner prediction. |
| 18 | DNABERT v1 | Skip. |
| 19–20 | GENA-LM | 36-kb-context BigBird encoder — covers a whole *CYP3A4* locus plus flanks in one pass. |
| 21 | ModernGENA | **Best first thing to actually run.** TSS-centred pretraining, pure torch, 377M — the cheapest model that is both likely to work and regulatory-relevant. |
| 21b | GENAtator | Supervision: per-nucleotide gene-structure labels to condition the DNA tower. |
| 22 | JEPA-DNA | **The prior art and the bar.** Use its checkpoints as a pre-JEPA'd DNA tower and as the control that says whether our own JEPA head adds anything. |
| 23 | gLM2 | **The bridge model.** One vocabulary over amino acids *and* nucleotides — the only pretrained model that already puts protein and DNA in a single embedding space. |
| 24 | gLM v1 | Conceptual reference for gene-neighbourhood context; not a per-base encoder. |
| 25 | GenomeOcean | 50-kbp-context encoder for operon/BGC-scale context; prokaryote-biased, so only for the microbial arm. |
| 26 | METAGENE-1 | Skip — 512-token context. |
| 27 | GROVER | Human-only 768-d baseline; useful only as a third opinion. |
| 28 | SegmentNT | Supervision: 14 per-nucleotide element classes (exon/intron/UTR/splice site) as auxiliary targets for the DNA tower. |
| 29 | ChatNT etc. | Not embedding sources. Isoformer/BulkRNABert are RNA-side; revisit for the RNA report. |
| 30–31 | AlphaGenome (+ torch port) | ⚠️ **Comparator only — its terms forbid training other models on its outputs.** Report against it; never embed with it. |
| 32 | Enformer | **Functional encoder (promoted).** 3072-d @128 bp, CC-BY-4.0 — the permissively-licensed tower now that AlphaGenome is comparator-only. |
| 33–35 | Borzoi / borzoi-pytorch / Flashzoi | **Primary functional encoder.** ~1920-d @32 bp — finest-grained available, only one carrying RNA-seq, CC-BY-4.0, and the port claims exact parity with the original. |
| 36 | Decima | Supervision: single-cell expression targets; gives the DNA tower a *liver-hepatocyte* readout, which is what CYP3A4 actually needs. |
| 37 | scooby | Channel: per-cell-type profile conditioned on a cell embedding — only if we go single-cell. |
| 38 | gReLU | **Infrastructure, not a model.** One torch API over Enformer/Borzoi/Decima — use it to avoid writing four loaders. |
| 39 | ChromBPNet | Supervision: base-resolution, bias-corrected accessibility over 408 biosamples. Blocked by TensorFlow; worth a conversion effort. |
| 40 | Sei | Channel: its 40 "sequence classes" are a ready-made 40-d interpretable regulatory embedding — cheap to concatenate. |
| 41 | DeepSEA/Beluga/ExPecto | Channel: 919-/2002-d chromatin vector, orthogonal and nearly free. Academic-only licence. |
| 42–44, 46, 50 | Basenji, Akita, Basset/DanQ, scBasset, Enformer-Celltyping | Skip — superseded and TensorFlow-bound. |
| 45, 47 | Orca, C.Origami | Channel: 3D-contact features at 1–256 Mb. Orthogonal modality, awkward deps, low priority. |
| 48 | GET | Not sequence-input. Skip unless we acquire matched ATAC. |
| 49 | CREsted | Train-your-own scATAC enhancer models — a wave-3 option, not a frozen source. |
| 51 | regLM | Generative; skip. |
| 52 | **DeepPBS** | **The section-2↔3 bridge.** Predicts protein–DNA specificity from structure → gives us (protein embedding, DNA embedding, affinity) triples with real coordinates behind them. |
| 53–54 | Geneformer, scGPT | Different space (genes, not bases). Only relevant if the world model gains an expression tower. |
| 55 | DART-Eval | Benchmark. The difficulty-matched split that keeps us honest per README kill-criterion 3. |

---

## 3. DNA structure data at scale

Every count below was **queried live against the RCSB Search API on 2026-09-20**, and
every size marked *(measured)* was obtained by sampling real files, not quoted from a
paper. Two of us ran the PDB queries independently and got identical numbers.

### 3.1 The ceiling, measured

| quantity | count | how |
|---|---|---|
| PDB experimental entries | **259,693** (holdings list: 260,089 IDs) | rcsb.org/stats |
| entries with **DNA** | **13,578** | `polymer_entity_count_DNA ≥ 1` |
| entries with **protein *and* DNA** | **10,733** | `count_DNA ≥ 1 AND count_protein ≥ 1` |
| ↳ X-ray / cryo-EM / NMR | **6,978 / 3,584 / 145** | method facet |
| ↳ polymer **entities** (chains) in those | **57,382** — 5.3× the entry count | `return_type: polymer_entity` |
| ↳ **non-redundant protein clusters @30% id** | **3,027** | seq-cluster facet over 34,740 protein entities |
| entries with DNA and **no** protein | **2,845** | `count_protein == 0` |
| protein + RNA | 8,024 | |
| all NA-containing entries | **22,364** | |

Reproduce with a POST to `https://search.rcsb.org/rcsbsearch/v2/query`:

```json
{"query":{"type":"group","logical_operator":"and","nodes":[
 {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.polymer_entity_count_DNA","operator":"greater_or_equal","value":1}},
 {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.polymer_entity_count_protein","operator":"greater_or_equal","value":1}}]},
 "return_type":"entry","request_options":{"return_counts":true}}
```

**The honest N is 3,027, not 10,733.** That is the leave-one-target-out budget for
protein–DNA — the same order of magnitude as the 185 targets in this repo's P450
universe, not the "tens of thousands" the raw entry count suggests. Against 8.8 trillion
bases of sequence corpus that is a ratio near 10⁻⁹: **the DNA arm is
structure-supervision-starved by four or five orders of magnitude.** This is the single
strongest argument for the frozen-encoder + light-head design the README already commits
to, and against anything that wants to learn DNA structure end to end.

### 3.2 Table

| # | Dataset | Ver | Yr | URL | Open? | License | Records | Size | Use in the shared space |
|---|---|---|---|---|---|---|---|---|---|
| D1 | **PDB protein–DNA subset** | live | 2026 | [search.rcsb.org](https://search.rcsb.org) → per-ID `files.rcsb.org/download/<id>.cif.gz` | yes | **CC0** | **10,733 entries / 3,027 clusters** | **4.32 GB gz *(measured*, 59-entry sample @393 KB); ~15–20 GB uncompressed** | The whole structural protein–DNA universe. Interface contact maps = the only true geometric supervision for a DNA↔protein link. |
| D2 | **PDB DNA-only** | live | 2026 | same | yes | CC0 | **2,845** (2,770 with no RNA either) | <1 GB | Sequence→shape targets: groove width, roll/twist/propeller from sequence alone. |
| D3 | Full PDB mmCIF archive | weekly | 2026 | `rsync://rsync.rcsb.org` :33444 · [wwPDB FTP sites](https://www.wwpdb.org/ftp/pdb-ftp-sites) | yes | CC0 | 260,089 | **~85 GB gz *(measured*, 249-entry sample)**; whole `ftp_data` tree ~787 GB | Don't mirror. Pull the 13,578 NA IDs. FTP is dead; HTTPS + rsync only. |
| D4 | PDB **NextGen Archive** | since 2023-02 | 2026 | `rsync://rsync-nextgen.rcsb.org` | yes | CC0 | same entries | **~103 GB** | mmCIF with UniProt/Pfam/SCOP/CATH mappings baked in. **Adds annotation, not coordinates** — the RCSB GraphQL API gives the same mappings for free. |
| D5 | **NAKB** (successor to the **NDB**) | 2024→ | 2026-09-16 | [nakb.org](https://nakb.org) · [NAR 52:D245](https://academic.oup.com/nar/article/52/D1/D245/7416380) | yes, no registration | open | 22,417 NA structures | **annotation dump = 6.9 MB JSON *(measured*)** | ⚠️ **The NDB was retired July 2023**; `ndbserver.rutgers.edu` 302s to nakb.org. Its 22,417 vs our 22,364 is a **0.24% difference** — it *is* the PDB. Worth exactly 7 MB: NA conformational class, entity typing, plus ideal A/B/Z fiber models. |
| D6 | **DNAproDB** | **2.0**, NAR 53:D396 | 2025, dump **2026-09-19** | [dnaprodb.usc.edu](https://dnaprodb.usc.edu) · [`dnaprodb2.tar.gz`](https://dnaprodb.usc.edu/data/dnaprodb2.tar.gz) · [pipeline](https://github.com/dnaprodb/dnaprodb-pipeline) | yes, no registration | **CC-BY-NC** (data); pipeline GPL-3.0 | **6,731 structures** (weekly auto-update) | **0.94 GB *(measured*)**; legacy v1.2 `.7z` **187 MB *(measured*)** | **The best derived layer over D1**: per-residue/per-nucleotide H-bonds incl. water-mediated, groove of contact, base-pair geometry, SASA. ⚠️ Its own page claims "several GB" for v1.2; it is 187 MB. **Run the pipeline yourself on predicted poses** so crystal and prediction get identical labels — that is what turns interface features into a *loss* rather than a metric. |
| D7 | **BioLiP2** | NAR 52:D404 | 2024, weekly | `zhanggroup.org/BioLiP/` · mirror `aideepmed.com/BioLiP/` | yes | BSD (code), data free | 385,160 protein chains / 781,684 chain–ligand interactions; **DNA 36,784** | ~GB, TSV per ligand type | ⚠️ **403s to non-browser clients.** The canonical residue-level binding-site labels — and the source almost every benchmark below quietly re-extracts. "36,784 DNA" is **chain–ligand pairs, not complexes.** |
| D8 | **DeepPBS dataset** | 2024 | Nat Methods 21:1674 | [figshare 25678053](https://doi.org/10.6084/m9.figshare.25678053) · [GH](https://github.com/timkartar/DeepPBS) | yes (figshare 403s to bots; browser OK) | BSD-3 | cross-family protein–DNA + measured specificity | not stated | **The only curated structure → binding-specificity (PWM) supervision.** Ideal auxiliary head: embed complex → predict PWM. |
| D9 | **ProNAB** | NAR 50:D1528 | 2021, **frozen Feb 2023** | [web.iitm.ac.in/bioinfo2/pronab](https://web.iitm.ac.in/bioinfo2/pronab/) | yes, `#download` | not stated ⚠️ | **20,219**: **14,606 protein–DNA** + 5,323 protein–RNA + 161 hybrid; 1,027 unique NA-binding proteins, lit. 1979–2021 | small (tabular Kd / ΔG) | **The affinity labels — the actual JEPA regression target.** Join to PDB by UniProt + sequence. Stale since Feb 2023. |
| D10 | **Deep DNAshape** | Nat Commun 15:1243 | 2024 | [GH](https://github.com/JinsenLi/deepDNAshape) · [server](https://deepdnashape.usc.edu/) | yes | **BSD-3** | 13 shape features + fluctuations, arbitrary flanks, MD-trained | repo 61 MB incl. weights | **Strictly better than DNAshapeR.** Differentiable per-bp geometry for *any* sequence — lets the DNA branch be trained at genome scale rather than on 13.5k PDB entries. |
| D11 | DNAshapeR / DNAshape | 2013/16 | — | [rohslab.usc.edu/DNAshape](https://rohslab.usc.edu/DNAshape/) · Bioconductor | yes | GPL-2 | pentamer MC lookup: MGW, HelT, ProT, Roll | <10 MB | Cheap input channel. Superseded by D10. |
| D12 | **BIGNASim** | NAR 44:D272 | 2016, **frozen** | [mmb.irbbarcelona.org/BigNASim](https://mmb.irbbarcelona.org/BigNASim/) | yes, browsable | — | MD trajectories: B-DNA (ABC/parmbsc0/bsc1), Z-, A-DNA, G-loops, quadruplexes, triplexes, + DNA–protein; 100 ns–10 µs | TB-scale, not published | The only large source of DNA conformational **ensembles** — a flexibility/entropy target the static PDB cannot give. Paper is 10 years old. |
| D13 | **DSSR-G4DB** | live, **2026-09-16** | 2026 | [g4.x3dna-dssr.org](https://g4.x3dna-dssr.org/DSSR-G4DB.html) | browsable; no bulk file | — | **608 G4-containing PDB entries**, DSSR-annotated | tiny | The cleanest non-duplex-DNA set. **Critical edge cases** so the DNA branch isn't a pure B-form-duplex memoriser. |
| D14 | **ONQUADRO** | NAR 50:D253 | 2022 | [onquadro.cs.put.poznan.pl](https://onquadro.cs.put.poznan.pl/) | yes | — | experimental quadruplexes, tetrads, geometric descriptors | small | Same role as D13, with per-tetrad geometry. |
| D15 | **non-B DB** / **non-B_gfa** | 2011 / pushed **2026-06-22** | — | [nonb-abcc.ncifcrf.gov](https://nonb-abcc.ncifcrf.gov/) · [GH](https://github.com/abcsFrederick/non-B_gfa) | yes | MIT (tool) | Z-DNA, G4, inverted/mirror repeats, cruciform, triplex, STRs, genome-wide | genome-scale GFF / small tool | ⚠️ **The database serves data built in June 2012.** Run the *tool* on a current genome instead. |
| D16 | PDNA-224 / -316 / -335 / -543 / -960, DNA-573/-129/-181, TransBind bundle | 2008–2024 | — | [ULDNA](https://github.com/yiheng-zhu/ULDNA) · [GraphSite](https://github.com/biomed-AI/GraphSite) · [Zenodo 10215073](https://zenodo.org/records/10215073) | yes | MIT / CC-BY-4.0 | 224–960 chains each | TransBind bundle **2.3 GB** (easiest one-stop) | Residue-level DNA-binding benchmarks. See the redundancy warning below. |
| D17 | **NPIDB** | 2013/16 | — | `npidb.belozersky.msu.ru` | — | — | — | — | ❌ **DEAD** — redirect loop → "temporarily unavailable" (checked 2026-09-20). |
| D18 | **PDIdb** | BMC Bioinf 11:262 | 2010 | `melolab.org/pdidb/` | — | — | 922 complexes ≤2.5 Å | — | ❌ **DEAD** — TCP connect refused (checked 2026-09-20). Survives only inside republishing papers. |
| D19 | **G4RNA** | 2015 | — | — | — | — | — | — | ❌ **DEAD** — both hostnames fail. Use G4Atlas or G4Bank. |
| D20 | *"hmmDB"* | — | — | — | — | — | — | — | ⚠️ **Could not verify that any such protein–DNA structural database exists.** Every search resolves to HMMER / HMM-profile *methods* papers. Flagged rather than invented. |

### 3.3 What the structure predictors trained on — and the split we should inherit

| Model | NA training data | PDB cutoff | Public split | License |
|---|---|---|---|---|
| **RoseTTAFold2NA** | RF2 protein data + all RNA/protein–RNA/protein–DNA in the PDB; 60:40 sampling | **2020-05** | described (224 complexes / 116 clusters) but **no split file** | code MIT; weights 1.1 GB. ⚠️ needs **~480 GB of sequence DBs** (UniRef30 46 GB + BFD 272 GB + nt 151 GB + RNAcentral 12 GB) |
| **AlphaFold3** | PDB + RNA distillation (Rfam/RNAcentral) + **protein–DNA distillation from JASPAR 2022 × SELEX (Jolma 2015, Yin 2017)** | **2021-09-30** | eval set = 8,856 complexes 2022-05-01→2023-01-12 | weights gated, non-commercial, application required |
| **Boltz-2** | PDB + **protein–DNA distillation built like AF3: JASPAR 2024 CORE × two SELEX sets, 10 ssDNA motifs per PFM + reverse complement, kept if PDE≤2.0 / iPDE≤1.0 / ipTM≥0.7, *no sequence clustering***. DNA–protein sampling weight **0.045** — the largest distillation modality after AFDB | train **2023-06-01**; validation to 2024-01-01, 398 structures, rule "retain all structures containing RNA or DNA entities" | procedure documented (App. A.1.5); files not shipped | **MIT, code + weights** |
| **Chai-1** | PDB + distillation | not extracted | no | **Apache-2.0, code *and* weights** — the most permissive of the four |
| Protenix / v2 | AF3 reproduction, weighted-PDB recipe | follows AF3 | pipeline documented | open |

**CASP is not a usable blind DNA set.** CASP15 (2022) was RNA-only. CASP16 (2024) had 42
NA targets: **35 RNA monomers, 1 DNA monomer, 11 RNA multimers, 18 NA–protein complexes**.
The assessors' own verdict: *"prediction accuracy for nucleic acid complexes was generally
poor unless 3D templates were available."* So **across CASP15+16 there is one blind DNA
monomer and 18 hybrid complexes** — we will be building our own held-out set from a
post-2024 PDB cutoff, which is fine, because 2024–2026 depositions are absent from every
published benchmark.

### 3.4 Skeptic's ledger

1. **There is essentially ONE protein–DNA structural dataset in the world: the PDB.**
   NAKB (0.24% apart, 7 MB of extra annotation), DNAproDB (a filtered subset + derived
   features), BioLiP2 (a re-extraction), NPIDB, PDIdb and all twelve PDNA/DNA-*n*
   benchmarks are re-wraps of it. **Budget one 4.3 GB download and one feature pipeline,
   not twelve datasets.**
2. **Per-chain vs per-complex inflates everything ~5×.** 10,733 entries → 57,382 polymer
   entities. BioLiP2's "36,784 DNA" is chain–ligand *pairs*. Never quote entity counts as
   structure counts. (This repo has a memory note for exactly this class of error.)
3. **The residue-level benchmarks are nested and stale.** PDNA-224 (2008) ⊂ PDNA-316
   (2011) ⊂ PDNA-543 (2014) by construction — they are the same chains binned by release
   date. Cross-testing between them **leaks**. DNA-573 cuts at 2016-01-06. None of them
   sees the **3,584 cryo-EM protein–DNA structures that are 33% of the field today** and
   are where nucleosomes, replisomes, CRISPR and transcription complexes live.
4. **Four resources are dead or zombie**: NPIDB (dead), PDIdb (dead), G4RNA (dead),
   non-B DB (serving 2012 data), ProNAB (frozen Feb 2023). Do not design around them.
5. **Beware AF3/Boltz-2 DNA pseudo-labels.** Both distil protein–DNA from **JASPAR ×
   SELEX**. Using their outputs as labels re-ingests JASPAR-derived *synthetic* DNA, not
   experimental structure — and Boltz-2 explicitly skipped the sequence clustering AF3
   applied, so that synthetic set is redundancy-heavy. A "structural" signal that traces
   back to a PWM is not independent evidence.
6. **Published size claims are unreliable.** DNAproDB says "several GB"; the file is
   187 MB. Every size above with *(measured)* was checked.

### 3.5 Practical minimum build

- `rsync`/HTTP the **10,733-ID protein–DNA subset — 4.3 GB, ~30 minutes.**
- Add the **NAKB annotation JSON (7 MB)** for NA conformational class. Skip NextGen
  (103 GB) — use the RCSB GraphQL API for UniProt/Pfam mappings.
- Run the **DNAproDB pipeline (GPL-3.0)** on crystal *and* predicted structures so
  interface labels are computed identically for both.
- Cluster at 30% → **3,027 folds**; hold out whole clusters, never entries.
- Auxiliary heads: DeepPBS (structure→PWM), ProNAB (structure→ΔG/Kd, 14,606 protein–DNA),
  Deep DNAshape (sequence→13 geometric features, distillable at genome scale).
- Edge cases so the DNA branch isn't a B-form memoriser: DSSR-G4DB's 608 G4 entries,
  ONQUADRO, BIGNASim ensembles.
- Held-out test: PDB entries released after each model's cutoff (Boltz-2 2023-06-01,
  AF3 2021-09-30, RF2NA 2020-05) plus the 18 CASP16 NA–protein targets.

**Total structural acquisition: well under 10 GB.**

---

## 4. Functional-genomics / binding data at scale

No coordinates here — these are *binding and activity* readouts, three to four orders of
magnitude larger than section 3. **Every size marked ✓ is a real byte sum or
`Content-Length` obtained live on 2026-09-20** against the primary API (ENA Portal,
ENCODE portal, JASPAR API, GEO/NCBI FTP, GTDB, GCS/S3, Zenodo), not quoted from a paper.

### 4.1 HT-SELEX and successors

| # | Dataset | Yr | Accession | Open? | License | Records ✓ | FASTQ size ✓ | Use |
|---|---|---|---|---|---|---|---|---|
| E1 | **Jolma 2013** human TF HT-SELEX | 2013 | ENA **PRJEB3289** | yes | ENA open | 2,726 runs, 490,439,803 reads; 830 profiles → **239 distinct specificities**, ~410 TFs | **10.76 GB** | The founding in-vitro TF↔DNA corpus, cheap enough to embed whole. Round-wise enrichment = a natural affinity-ranking objective. |
| E2 | Yang/Jolma deep re-sequencing of the **same** libraries | 2017 | ENA **PRJEB14744** | yes | open | 2,510 runs, **4.34 B reads** | **118.27 GB** | 9× the depth of E1 on identical libraries. Use this if you want quantitative counts; **do not use both**. |
| E3 | **Yin 2017 methyl-HT-SELEX** | 2017 | ENA **PRJEB9797** | yes | open | 7,443 runs, **8.79 B reads**; **542 TFs ±CpG methylation** | **507.74 GB** | The only large paired methylated/unmethylated corpus — lets the embedding carry a 5mC channel instead of treating DNA as 4-letter. |
| E4 | **Zhu 2018 NCAP-SELEX** (TF vs nucleosome) | 2018 | ENA **PRJEB22684** | yes | open | 3,962 runs, 4.65 B reads, 220 TFs | **456.54 GB** | The only in-vitro data encoding chromatin context. |
| E5 | **Codebook GHT-SELEX** (genomic HT-SELEX) | **2026** | ENA **PRJEB76622** · [codebook.ccbr.utoronto.ca](https://codebook.ccbr.utoronto.ca) · Nat Methods `10.1038/s41592-026-03177-9` | yes | ⚠️ **no licence posted** | 5,694 runs, 7.59 B reads; 1,534 GHT + 1,578 HT experiments; **139 of 331 uncharacterised TFs** succeeded | **235.71 GB** raw; **Triple-Optimized peaks 73.3 MB** | ★ **The most valuable new TF–DNA resource since 2017.** Uses fragmented *real genomic* DNA, so motifs come with native flanking context — directly comparable to ChIP-seq peaks in the same space. |
| E6 | Codebook ChIP-seq companion | 2026 | ENA **PRJEB78913** | yes | no licence | 1,107 runs, 52.5 B reads | **2.78 TB** raw; peaks only **94.1 MB** | In-vivo counterpart for the same dark TFs. |
| E7 | **Expanded codebook of human TF specificity** | 2026 | Nature `10.1038/s41586-026-10798-9` | yes | not stated | 332 putative + 61 control TFs; ~4,873 experiments; motifs for **177 (53%)** | PWM bundle **1.86 MB** | Eight years, five platforms, ~5k experiments → **1.9 MB of PWMs**. |
| — | *"ENCODE HT-SELEX"* | — | — | — | — | ❌ **Does not exist.** `assay_title=HT-SELEX` returns **0** ✓ | — | ⚠️ HT-SELEX motifs "in ENCODE" arrive via **Factorbook**, which *imports* Jolma 2013 / Yin 2017. Not an independent dataset. |

Total raw HT-SELEX-family FASTQ = **~4.17 TB**, of which 2.78 TB is Codebook ChIP-seq and
only ~1.39 TB is actual SELEX. **Nobody should re-derive motifs from raw reads** — the PWM
outputs of all of it are single-digit MB. The reads matter only for modelling *round-wise
enrichment* as a quantitative affinity signal, which is the one thing PWMs discard and the
one thing a shared embedding could exploit.

### 4.2 Motif banks — nested, not complementary, and tiny

| # | Resource | Ver | Yr | URL | Open? | License | Records ✓ | Size ✓ | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| E8 | **JASPAR** | **2026 (11th)** | 2026 | [jaspar.elixir.no](https://jaspar.elixir.no/) · NAR 54(D1):D184 · REST API, no key | yes | **CC-BY-4.0** | **CORE = 2,633** non-redundant (vertebrates 1,019, plants 927, insects 296…); +306 new; **1,259 BPNet DL models** (first DL collection) | **jaspar.zip 627 KB**, MEME 1.14 MB; genome-wide TFBS bed.tar.gz 298 MB | The only large motif bank with an unambiguous, commercially-safe licence. ⚠️ UNVALIDATED count: paper says 1,231, live API returns 1,031. ⚠️ **`JASPAR2026` 404s in Bioconductor release *and* devel**; pyJASPAR stops at 2024; **zero JASPAR datasets on HuggingFace**. |
| E9 | **HOCOMOCO** | **v14** | 2023/24 | [hocomoco14.autosome.org](https://hocomoco14.autosome.org/downloads_v14) | yes | **WTFPL** (treat as CC-BY) | **1,595 motifs / 1,107 human TFs**, 809 mouse; H14CORE-CLUSTERED 648 | **PWM tar.gz 868 KB** | Best-curated human/mouse bank. ⚠️ **v14 has FEWER motifs than v13** (1,595 vs 1,611) — a curation revision, not an expansion. ⚠️ **No v14 Zenodo deposit exists**; the cited DOI resolves to a record titled "v13". |
| E10 | **CIS-BP** | **Build 3.10** | **2026-04-26** | [cisbp.ccbr.utoronto.ca/bulk.php](https://cisbp.ccbr.utoronto.ca/bulk.php) | yes, no login | ⚠️ **no licence statement anywhere**; silently re-exports TRANSFAC | **13,030 motifs; 169,272 TFs with ≥1 motif — only 4,989 direct experimental**; 741 species | PWMs.zip **5.76 MB**; full archive **~6.4 GB** | ⚠️ **CIS-BP ⊇ JASPAR + HOCOMOCO + TRANSFAC + UniPROBE + Factorbook** (its own FAQ: >70 sources). "HOCOMOCO + JASPAR + CIS-BP" **is** CIS-BP. 164k of its 169k TF→motif edges are DBD-homology **inferences, not measurements**. |
| E11 | **UniPROBE** (PBM) | — | last new data **2021** | `thebrain.bwh.harvard.edu/uniprobe/` | partly | **academic click-through** | ~708 accessions, 34 publications | — | ⚠️ **Dying in public** — its homepage solicits a volunteer student. Only non-absorbed asset is continuous 8-mer E-scores. Mirror; take no runtime dependency. |
| E12 | **TRANSFAC** | — | rolling | genexplain.com/transfac | ❌ **no** | **proprietary, €3,200/yr academic** | — | — | Unusable in an open model — **and CIS-BP re-exports it**, a redistribution landmine. |
| E13 | SwissRegulon | matrices **2016-10-12** | — | swissregulon.unibas.ch | yes | none stated | ~190 human WMs | 340 KB | **Skip.** Ten years stale. |

**These are three estimators of one quantity.** JASPAR CORE vertebrates, HOCOMOCO H14CORE
and CIS-BP's 4,989 direct motifs derive from the *same* HT-SELEX / PBM / ChIP-seq
experiments, differing in inference algorithm and curation. Use their **disagreement as a
confidence weight**, never as three datasets. And note the scale: HOCOMOCO's entire
human+mouse PWM set is **868 KB**; JASPAR CORE 2026 is **627 KB**. **Motifs are a fixed
featurizer, not a data modality** — anyone budgeting them in GB has misread the resource.

### 4.3 ENCODE, and the peak/signal decision

Portal API, 2026-09-20 ✓: **28,642 Experiments indexed (27,043 released), 1,659,000 File
objects.** By phase: ENCODE4 **11,552** · ENCODE3 6,433 · ENCODE2 2,805 · Roadmap 2,763 ·
modERN 2,325. ⚠️ **There is no ENCODE5 RFA on the portal.**

Assays: **TF ChIP-seq 5,360 (5,002 released) across 2,168 distinct targets and 218
biosamples** · Histone ChIP-seq 3,992 · DNase-seq 3,582 · ATAC-seq 560 · intact Hi-C 366 ·
snATAC-seq 370. FunctionalCharacterizationExperiments: 869, of which **MPRA 125,
STARR-seq 42**.

| Slice (status=released) | Files ✓ | Size ✓ |
|---|---|---|
| **All TF ChIP-seq narrowPeak BEDs** | 40,662 | **71.0 GB** |
| All histone ChIP-seq narrowPeak | 26,846 | 44 GB |
| All DNase-seq narrowPeak | 6,552 | 9 GB |
| All ATAC-seq narrowPeak | 2,676 | 7 GB |
| All TF ChIP-seq **bigWigs** | 61,372 | **43.19 TB** |
| All TF ChIP-seq **FASTQ** | 14,941 | 21.55 TB |
| **ALL ENCODE** bam / fastq / bigWig / bed / bigBed | — | **≈1.71 PB** |

**Take the peaks (131 GB for all four assays), never the signal (>100 TB).** No account
needed; bulk via search → `files.txt` → `xargs curl`; API rate limit **10 GET/s**.
Licence, verbatim: *"External data users may freely download, analyze and publish results
based on any ENCODE data without restrictions."* **The cleanest licence in this report.**

**SCREEN / cCRE Registry V4** (Nature `10.1038/s41586-025-09909-9`): human GRCh38
**2,348,854 cCREs across 1,888 cell types — 129.1 MB** ✓ (PLS 47,532 · pELS 249,464 ·
dELS 1,469,205 · CTCF-bound 948,642); mouse 926,843 / 50.6 MB; T2T lift 123 MB. A **2.5×
expansion over V3** and **the single highest-value-per-byte file in this report: 129 MB
that indexes the regulatory genome.**

### 4.4 Aggregated ChIP-seq atlases — one dataset, four wrappers

| # | Resource | Real last update | Open? | License | Records ✓ | Size ✓ | Verdict |
|---|---|---|---|---|---|---|---|
| E14 | **ChIP-Atlas** | **actively maintained** — metadata 2026-09-09 | yes, open tree | **CC-BY-4.0** | 867,660 rows ≈ **454,500 unique SRX**; hg38 197,044 | metadata 353 MB; hg38 `allPeaks_light.05` **22.17 GB**, `.50` **3.39 GB**; bigWigs **~40 TB hg38** | ★ **Pick this one.** Widest coverage, permissive licence, only atlas with matched ATAC/DNase/bisulfite on one pipeline. ⚠️ hg19+hg38 are the same experiments twice; **only ~17% of hg38 rows are TF ChIP-seq**; peak bundles dated 2024-11-13. |
| E15 | **ReMap 2022** | ⚠️ **changelog's last entry 2021-09-14. There is no ReMap 2024/2025/2026** — those hosts do not resolve | yes | ⚠️ **CC-BY-NC-4.0** | human **8,103 datasets, 1,210 regulators, 182.4 M peaks** (68.2 M non-redundant, **3.4 M CRMs**) | hg38 all-peaks **4.52 GB**, non-redundant 1.46 GB, CRMs 200 MB | Best-QC'd set, and its **3.4 M CRMs** are the right granularity for a regulatory-element token vocabulary. **Five years frozen and NC-licensed.** ENCODE is **25.3%** of ReMap-human ✓. |
| E16 | **GTRD v21.12** | ⚠️ **abandoned** — `downloads/current/` serves files dated **2020-09-20** | yes | ⚠️ none posted | 60,285 experiments (human 27,503) | human MACS2 8.3 GB, meta-clusters 6.6 GB | Only non-duplicated asset is the **4-caller meta-cluster consensus** (a peak-confidence prior). Otherwise superseded. |
| E17 | **Cistrome DB v3.0** | 2024; **front-end degrading 2026** — maintenance banner, Django `DEBUG=True` in production, all `/api/*` routes 404 | ⚠️ not programmatically | none stated | ~45k human + ~44k mouse samples | not obtainable anonymously | ❌ **Do not build on it.** Harvest its per-sample QC table if you can, nothing else. |

**Quantified redundancy:** ChIP-Atlas hg38 TF (33,368) vs GTRD human (27,503) vs ReMap
human (8,103) — all ingesting the same GEO/SRA submissions, differing only in peak caller,
QC filter and genome build. Union them and you weight the same GSE two or three times:
a silent duplication weight, not extra signal. De-duplication by metadata text is
impossible — **only 1,750 of 197,044 ChIP-Atlas hg38 rows even contain the string
"ENCODE."** Pick one.

### 4.5 Hi-C / 3D genome

| # | Dataset | Yr | URL | Open? | License | Records ✓ | Size ✓ | Use |
|---|---|---|---|---|---|---|---|---|
| E18 | **Akita training tensors** | 2020 | `gs://basenji_hic/1m/data/tfrecords/` | **yes — plain public GCS, NOT requester-pays** | Apache-2.0 | **5 targets**: HFF, H1hESC, GM12878, IMR90, HCT116; 1 Mb windows, 2048 bp bins | **10.71 GB**, 32 tfrecords | ★ **The single most reusable artefact in this section.** Five cell types already binned, clipped and sharded. A drop-in Hi-C head for 10 GB and $0. |
| E19 | **4D Nucleome** | live | `data.4dnucleome.org`; mirror **`s3://4dn-open-data-public`** | ⚠️ portal `@@download` **403s anonymously**; **the S3 mirror is open** | AWS Open Data, no restrictions; portal's own `/about/data-use-policy` **404s** | 3,392 ExperimentSetReplicates (1,988 Hi-C; **only 26 Micro-C**), 66,440 files | **365.7 TB total; 324.2 TB open — but 91% is reads.** Consumable part: **.hic 3.31 TB + .mcool 2.71 TB ≈ 6.2 TB** | ⚠️ Live banner: *"under review for potential modification in compliance with Administration directives."* **Mirror what you need.** |
| E20 | **Rao 2014 GM12878** | 2014 | GEO **GSE63525** | yes | public domain | 200 GSMs, 4.9 B contacts, 1 kb res | **1.66 TB total — but 1.0 TB is one `RAW.tar`**; the `.hic` you want is **51 GB** (MAPQ0) / 37 GB (MAPQ30) | Still the densest single human map. |
| E21 | **Micro-C H1/HFFc6** | 2020 | 4DN `4DNES21D8SP8`, `4DNESWST3UBH` | yes via S3 | as 4DN | 10 / 14 replicates | processed **81.6 GB / 130.0 GB** | The Akita/Orca targets. |
| E22 | **Orca resources** | 2022 | [Zenodo 6234936](https://zenodo.org/records/6234936) | yes | **CC-BY-4.0** | H1 + HFF Micro-C, 1 kb→1 Mb, to 256 Mb input | core **1.33 GB**, mcools 35.71 GB | Chromosome-scale contact supervision, pre-binned. |
| E23 | **C.Origami data** | 2023 | [Zenodo 7226561](https://zenodo.org/records/7226561) | yes | CC-BY | IMR-90 + CTCF ChIP + ATAC, 2 Mb windows | **1.81 GB** | Cheapest worked example of *conditioning contacts on cell-type tracks* — exactly the shared-embedding pattern. |
| E24 | **HiCFoundation** | **2026** | [GH](https://github.com/Noble-Lab/HiCFoundation) · HF `wang3702/hicfoundation_models` | weights yes; **corpus not released** | **Apache-2.0** | "hundreds of Hi-C assays, 118 M patches", 316 species, 7 checkpoints | weights only | ★ **A free pretrained Hi-C encoder you can graft in.** Its epigenomic head already does Hi-C→ATAC/ChIP — a partial shared-embedding model that already exists. |
| E25 | **Human body single-cell 3D + methylation atlas** | **2026-07-31** | GEO **GSE326618**, Science `adx0673` | **yes, fully open** | public domain | **86,689 single nuclei, 16 tissues, 206 cell subtypes** (snm3C-seq) | **`_RAW.tar` = 1.3 TB** | ★ **The biggest genuinely new 3D resource of 2026**, and the only one with paired methylation at single-cell scale. Not a re-cut of Rao/Krietenstein. |
| E26 | **Evo2HiC** | 2025-11 preprint | `10.1101/2025.11.18.689171` | preprint | — | multimodal DNA-FM embeddings ⊕ Hi-C | — | ⚠️ **Read before designing anything** — someone is already building the sequence+Hi-C half of this thesis, and already hit the "genome-wide sequence embeddings are too expensive" wall. |

**"Hi-C as structure" — the honest answer: no.** A bulk contact map is a
population-averaged frequency matrix whose signal is ~80% genomic-distance decay P(s);
inverting it to coordinates is non-unique and averages over an ensemble containing
**mutually exclusive conformations** — structurally the same failure mode as CYP3A4's F/G
loop, which this repo already knows how to lose to. Real 3D supervision comes from
single-cell Hi-C / snm3C (E25) and chromatin-tracing imaging. **Use bulk `.hic`/`.mcool`
as a 2D auxiliary head (as Akita and AlphaGenome do); never merge it into a geometry
loss.** Note also: **AlphaGenome's contact head is 2048 bp — the same resolution Akita
used in 2020. Sequence→contact resolution has not improved in six years.**

### 4.6 Reporter assays, expression, variant effects

| # | Dataset | Yr | URL / accession | Open? | License | Records ✓ | Size ✓ | Use |
|---|---|---|---|---|---|---|---|---|
| E27 | **lentiMPRA, Agarwal 2025** | 2025 | ENCODE **ENCSR022GQD / ENCSR382BVV / ENCSR244FWB** + Zenodo 13908857 | yes | **ENCODE — unrestricted** | **>680,000 cCRE sequences**, 41.7% active, 3 cell types | raw 17.8/34.8/16.6 GB; **element-quantification TSVs only 0.148/0.195/0.039 GB** | ★ **Use this as the MPRA backbone.** Permissive, and cCRE-anchored so it joins directly to SCREEN and the ENCODE tracks. |
| E28 | **lentiMPRA, Gosai 2024** (Malinois) | 2024 | [Zenodo 10698014](https://zenodo.org/records/10698014) | yes | ⚠️ **CC-BY-NC-4.0** | **776,474 sequences** in K562/HepG2/SK-N-SH | Zenodo record 72.74 GB — **but 62 GB is immunofluorescence images**; the table is **0.28 GB** | Largest single MPRA set, **NC-licensed**. ⚠️ boda2's README warns **the bioRxiv supplementary Table S2 was WRONG** — refetch from the GCS URL. |
| E29 | **SuRE-SNP** | 2019 | GEO **GSE128325** | yes | public domain | **~5.9 M SNPs assayed**, 4 haplotype libraries | **118.9 GB** | ★ **The only million-scale human allelic reporter supervision that exists openly.** A natural contrastive / variant-effect objective. |
| E30 | **de Boer 2020 / Vaishnav 2022 yeast random promoters** | 2020 | GEO **GSE104878** | yes | public domain | **~31.4 M sequence→expression pairs** in the pTpA file alone | **~3.1 GB** | ★ **Highest record count per GB here by two orders of magnitude.** Yeast, so no human transfer — but the ideal scaling-law testbed for a sequence→expression head. |
| E31 | **Kircher 2019 saturation-mutagenesis MPRA** | 2019 | GEO **GSE126550** | **yes, ungated** | public domain | 20 disease-associated elements at single-bp resolution | RAW 12 GB; tables MB-scale | ★ **CAGI5 without the DUA.** Dense, quantitative, single-bp — ideal regression target on embedding deltas. |
| E32 | Sharpr-MPRA | 2016 | GEO GSE71279 | yes | public domain | ~487k constructs | **82 MB** | ⚠️ **Eval only. Far too small to train on.** |
| E33 | **SuRE** | 2017 | GEO GSE78709 | yes | public domain | >10⁸ native fragments | counts 1.1 GB | Native-fragment promoter activity; complements oligo MPRA. |
| E34 | **FANTOM5 CAGE** | frozen 2017/2021 | `fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/` | yes | **CC-BY-4.0** | Enformer's CAGE head = **exactly 638 `CNhs*` libraries** of its 5,313 ✓ | peaks 5.2 MB; **TPM matrix 799 MB** | Transcription-initiation modality; small and permissive. |
| E35 | **Enformer training data** | 2021 | `gs://basenji_barnyard/data` | ⚠️ **REQUESTER-PAYS** (HTTP 400 anonymously) | Apache-2.0 code | human head **5,313 tracks**: CHIP 3,991 · DNASE 674 · CAGE 638 · ATAC 10 ✓ | "multiple TB", ~**$600** egress for 5 TB | Reference track inventory. Don't pay for it — use the weights. |
| E36 | **Borzoi training data** | 2023–25 | `gs://borzoi-paper/data` | ⚠️ requester-pays | Apache-2.0 code | human head **7,611 tracks**: CHIP 3,886 · **RNA 1,543** · CAGE 1,276 · DNASE 674 · ATAC 232 ✓ | multiple TB | Weights and QTL benchmarks are free; only the tfrecords cost money. |
| E37 | **GTEx v10 / v11** | v11 **2026-01-15** | `gs://adult-gtex/bulk-gex/v11` | **yes — open GCS, no auth** | open | **19,788 samples, 68 tissues, 946 donors** | v10 **42.55 GB** / v11 66.38 GB; gene TPM 2.26 GB | ⚠️ **v11 is a reprocessing, not new data — v10 and v11 `SampleAttributesDS.txt` are byte-identical (same MD5).** Do not re-download. Individual-level data is dbGaP-gated (2–4 months). |
| E38 | **Roadmap Epigenomics** | Release 9, files **2013–2016** | `egg2.wustl.edu/roadmap/data/byFileType/` | yes | ⚠️ **no machine-readable licence anywhere** | **127 consolidated epigenomes**; 1,032 pval bigWigs, 32 marks | **pval signal 605 GB** ✓, ~1.2 TB with fold-change | ⚠️ **`roadmapepigenomics.org` is dead (HTTP 000).** 1.2 TB of irreplaceable signal on one unmirrored university web server. ENCODE hosts 2,763 Roadmap experiments but not all 127 epigenomes. **Mirror it now if you want it.** |
| E39 | **ClinVar** | **2026-09-13** | `ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/` | yes | US-gov public domain | **4,471,671 records; 858,435 non-coding (19.2%)**, of which **55,597 P/LP**; intron 1,005,249 · 5′UTR 157,292 · splice-donor 60,909 | **193.7 MB** | The only non-coding label set big enough to matter, public-domain, weekly-recomputable, and **you control the split**. |
| E40 | **gnomAD v4.1.1** | **2026-03-30** | `gs://gcp-public-data--gnomad/release/4.1/` | yes, no requester-pays | MIT + terms; commercial OK | **807,162 individuals** (730,947 exomes + 76,215 genomes) | **1.64 TB** sites VCFs | Negative/common-variant control and AF covariate — essential to avoid a frequency-confounded probe. ⚠️ **gnomAD v5 does not exist.** |
| E41 | **BEND** | 2024 | [GH](https://github.com/frederikkemarin/BEND) | yes, no login | **code BSD-3; data CC-BY-4.0** | 7 tasks: chromatin_accessibility **2,062,129** · cpg_methylation 959,039 · histone 625,229 · variant_effects_disease 295,495 · gene_finding 5,977 · **enhancer_annotation 285** | **0.23 GB** (~1 GB working set) | ★ **Best-designed suite for frozen-embedding probing** — ships a webdataset embedding-precompute pipeline. **Drop `enhancer_annotation` (n=285).** |
| E42 | **Genomics Long-Range Benchmark (LRB)** | 2024 | HF `InstaDeepAI/genomics-long-range-benchmark` | yes, ungated | ⚠️ **CC-BY-NC-SA-4.0 — share-alike contaminates** | 9 tasks incl. variant_effect_causal_eqtl 88,717 · pathogenic_clinvar 38,634 · **pathogenic_omim 2,321,473** · chromatin 2,203,689 · enhancer 1,914,575 | streams; multi-GB | ★ **Best fit for variant-effect probing.** Three ready-made tasks with arbitrary `sequence_length` and pre-built ref/alt pairs. The eQTL task is Enformer's GTEx fine-mapped set — **and it gets you Enformer's eQTL slice free**, avoiding the requester-pays bucket. |
| E43 | **GUE** (DNABERT-2) | 2024 | mirror HF `leannmlindsey/GUE` | ⚠️ official data is a **bare Google Drive link — no version, no checksum, no licence** | code Apache-2.0; **data none** | mirror: 37 configs, **1,045,150 rows** — ⚠️ the mirror **adds `phage_fragments` + `fungi_species_20` which are NOT official GUE** | **1.06 GB** | Comparability with the DNABERT-2 literature only. |
| E44 | **NT benchmark (revised)** | 2024 | HF `..._downstream_tasks_revised` | yes | ⚠️ none declared | **532,064 rows**, chromosome-held-out | 0.199 GB rows / 0.59 GB repo | The original 18-task version is **deprecated by its own authors**; use this. |
| E45 | **Genomic Benchmarks** | 2023 | [GH](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks) | yes | **Apache-2.0** | 9 datasets, **890,705 sequences** | **0.147 GB** | Cheap linear-probe smoke test. **Not a result.** |
| E46 | **DART-Eval** | 2024 | [GH](https://github.com/kundajelab/DART-Eval) · Synapse **syn59522070** | ❌ **Synapse account + DUA** | ⚠️ **repo has no LICENSE (404)** | 5 tasks; task 5 = **African caQTLs + Yoruban dsQTLs** | not disclosable anonymously | ★ Task 5 is the **only held-out-population variant-effect test** in this list — the least ENCODE-recycled signal available. **Registration-gated: see BLOCKED.** |
| E47 | **DNALongBench** | HF, **2026-05-11** | HF `andyjzhao/dnalongbench` | yes | ⚠️ **none declared, no README at all** | 4 tasks: contact map · eQTL · enhancer–target gene · transcription initiation | multi-GB | Only 2026-dated suite found. **Treat as pre-release.** |
| E48 | MFASS / Vex-seq | 2018/19 | GH `KosuriLab/MFASS` · GEO GSE113163 | yes | ⚠️ MFASS has no LICENSE | 32,669 SNVs / ~2,000 variants | 41 MB / **1.6 MB** | MFASS is small but genuinely independent. **Vex-seq is too small to conclude anything.** |

### 4.7 Pretraining corpora

| # | Corpus | Yr | URL | Open? | License | Records ✓ | Tokens (bp) | Size ✓ |
|---|---|---|---|---|---|---|---|---|
| E49 | **OpenGenome2** (Evo 2) | 2025 | HF `arcinstitute/opengenome2` | yes, **ungated** | **Apache-2.0** | 28,177 GTDB prokaryotes + 15,044 NCBI eukaryotes + 32,240 organelles + metagenomes | **8.8 T claimed / 8.60 T itemized** | **5.52 TB — but `fasta/` 2.735 TB and `json/` 2.788 TB are the SAME sequences twice.** Pull one subtree. ⚠️ 76% of tokens are animal+plant from only 15,044 assemblies. |
| E50 | **GTDB R11-RS232** | **2026-04-15** | `data.gtdb.ecogenomic.org/releases/latest/` | yes | ⚠️ no explicit licence (INSDC-derived) | **901,341 genomes in 199,923 species clusters** | reps ≈0.6–0.7 T | **`gtdb_genomes_reps.tar.gz` = 192.2 GB** | ★ **The right prokaryote dedup. Use the 199,923 representatives, never raw RefSeq bacteria** — exactly what Evo 1/2 did. |
| E51 | **T2T-CHM13 v2.0 + GRCh38 no-alt** | 2022 | NCBI / UCSC `hs1`, `hg38` | yes | public domain | 24 / 195 seqs | **3,117,275,501 / 3,099,734,149 bp** | **933 MB / 873 MB** | Anchor in hg38 (every functional track is there); project T2T in for centromeric/segdup sequence. |
| E52 | **UCSC multiz100way (hg38)** | current | `hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/` | yes | free for any use | 100 species | — | **74.7 GB** | ★ **Start here for conservation, not multiz470way (1.22 TB) or the 447-way Cactus HAL (1.26 TB).** 16× smaller for most of the usable signal. |
| E53 | **HPRC Release 2** | 2025 | `s3://human-pangenomics` (`--no-sign-request`) | yes, no egress fee | public domain (466/466 GenBank-accessioned) | **466 assemblies / 234 samples** | ~1.4 T nominal | **421 GB as gz — but 3.295 GB in AGC. A 127× reduction.** | ⚠️ **Proof that pangenomes add almost no tokens**: human haplotypes are 99.9% identical, so 1.4 T bp carries maybe 5–10 Gbp of new information. |
| E54 | **1000 Genomes 30x** | 2022 | `ftp.1000genomes.ebi.ac.uk` | yes, **no DAC, no registration** | fully open | 3,202 samples, 125 M variants | — | **VCF chr1-22,X = 29.8 GB**; raw reads 230.2 TB | **Use the VCF.** The CRAMs are near-worthless as pretraining tokens. ⚠️ NT's `-1000g` corpus is **99.9% redundant** and its models underperform the 850-species ones: 6,400× the tokens for 1.00004× the information. |
| E55 | **Earth BioGenome / DToL** | live | goat.genomehubs.org | yes | INSDC public domain | **69,472 eukaryotic assemblies** (17,300 chromosome-level) | ≈50–100 T | via ENA/NCBI | ★ **The most under-used token source in genomics** — Evo 2 used only 15,044 eukaryotic genomes. |
| E56 | RefSeq 237 / GenBank 273 | 2026-08 | NCBI FTP | yes | public domain | 645 M / 6.7 B records | 6.70 T / **60.07 T** | 13.6 TB / ~11.6 TB+ | ⚠️ 20× species-redundant among bacteria; the worst signal-per-token ratio available. |
| E57 | MGnify / IMG-M / Tara | live | ebi.ac.uk/metagenomics · img.jgi.doe.gov | MGnify yes; **IMG/M requires JGI SSO** ❌ | varies | 56,782 species-rep MAGs / 240k datasets | 0.17 T / **~30 T** | tens of GB | IMG/M is **the biggest closed pool in genomics** — why Evo 2 used MGD/IMG-VR extracts instead. |

### 4.8 Skeptic's ledger

- **Four ChIP-seq "atlases" are one pile of GEO/SRA data.** ENCODE is 25.3% of ReMap-human
  and ~10% of ChIP-Atlas-human. Pick one. §4.4.
- **Three resources are dead or dying; mirror today if wanted.** ReMap (frozen 2021, and
  the "2022" branding is the most misleading recency signal in this report), GTRD (2020-era
  files), Roadmap (`roadmapepigenomics.org` returns nothing; 1.2 TB on one unmirrored
  server), Cistrome (`DEBUG=True` in production).
- **"Size" is misleading in six specific places**, all measured: 4DN's 365 TB is 91% reads
  (consumable: 6.2 TB); OpenGenome2's 5.52 TB is `fasta/` + `json/` of the *same*
  sequences; Gosai's 72.74 GB Zenodo record is 62 GB of microscopy images (the table is
  0.28 GB); HPRC R2 is 421 GB as gz and 3.295 GB in AGC; **the NT 850-genome HF repo is
  289 KB — a loader script, not data**; GSE63525's 1.66 TB is 1.0 TB of one `RAW.tar`.
- **Benchmarks reshuffle the same data.** Genomic Benchmarks' 8 core datasets all derive
  from Ensembl 97/100 and nothing else. **GUE's `emp_*` configs are the same 2006 yeast
  histone data as NT's; GUE's promoter splits are the same DeePromoter splits as NT's —
  GUE and the NT benchmark overlap on roughly half their datasets.** BEND
  histone/accessibility, DART-Eval tasks 1/3/4 and LRB `chromatin_features_*` are ENCODE
  again. **Scoring well on GUE *and* NT *and* Genomic Benchmarks is one result reported
  three times.**
- **Too small to matter — the FINDING-007 problem verbatim.** BEND `enhancer_annotation`
  n=285 (an 11 KB BED: any delta between two models is noise); Vex-seq ~2,000 variants;
  Genomic Benchmarks `dummy_mouse_enhancers` n=1,210 (its own authors named it "dummy" and
  papers still report it); Sharpr-MPRA 82 MB; 4DN has **26 Micro-C experiment-sets** against
  1,988 Hi-C, so "Micro-C at scale" does not exist. **Measure the permutation null before
  reporting any delta on these.** And watch saturation: HyenaDNA hits 96.6% on
  `demo_human_or_worm`, which is solvable by GC content, where a plain CNN gets 93.3% —
  that is a ceiling, not a benchmark, and README kill-criterion 3 forbids claiming on it.
- **Firm negatives.** gnomAD **v5 does not exist** (v4.1.1 is current). **ReMap 2025 does
  not exist.** **There is no ENCODE HT-SELEX dataset.** **There is no ENCODE5 RFA on the
  portal.** **JASPAR2026 is not in Bioconductor**, release or devel, despite the paper.
  **CAGI6 has no regulatory saturation-mutagenesis challenge.** GTEx v11 sample metadata
  is byte-identical to v10.
- **No licence at all** — "downloadable" ≠ "you may redistribute derived embeddings":
  CIS-BP, GTRD, Cistrome, SwissRegulon, Codebook, the DART-Eval repo, the official GUE
  data, the NT task datasets, DNALongBench, MFASS, Roadmap.
- **The scale reality.** The only million-scale human *sequence→function* supervision that
  exists openly is **SuRE-SNP (~5.9 M SNPs, 119 GB)**; everything else is sub-million. If
  you want scaling-law behaviour on sequence→expression today, **the data is in yeast**
  (de Boer GSE104878, 31.4 M pairs, 3.1 GB).

---

## 5. Recommendation

### 5.1 Models — take three

1. **Evo 2 7B** (`arcinstitute/evo2_7b`, Apache-2.0, ungated, 13.8 GB, 1 Mb ctx, 4096-d).
   The evolutionary-constraint tower, and the only genuinely unencumbered large DNA model.
   Budget a day for the environment: the repo wants torch 2.6/2.7 + flash-attn 2.8 +
   `vortex`, and we have 2.5.1. The 7B is the only size that runs bf16 without Transformer
   Engine.
2. **borzoi-pytorch** (`johahi/borzoi-replicate-{0..3}`, **CC-BY-4.0**, ~0.2B, ~1920-d
   @32 bp, 524 kb ctx, plain torch) — *or* **Enformer** (`EleutherAI/enformer-official-rough`,
   CC-BY-4.0, 3072-d @128 bp, 196 kb). The function tower. **This replaces AlphaGenome in
   the recommendation — see 5.3.** Borzoi is the finer-grained of the two and carries
   RNA-seq; its port claims exact parity with the original, where `enformer-pytorch` does
   not.
3. **ModernGENA-large** (`AIRI-Institute/moderngena-large`, 377M, pure `transformers`).
   TSS-centred pretraining on 443 vertebrate genomes, no mamba, no Triton, flash-attention
   optional. **Run this first** to get the extraction pipeline working before spending
   contended H200 hours on Evo 2.

### 5.2 Datasets — three, and they are small

1. **PDB protein–DNA subset — 10,733 entries, 4.32 GB, CC0** (→ 3,027 non-redundant
   clusters), plus **DNAproDB 2.0, 0.94 GB** for precomputed interface geometry.
2. **ENCODE narrowPeak BEDs — 71.0 GB TF ChIP-seq + 44 GB histone + 9 GB DNase + 7 GB
   ATAC = ~131 GB**, plus **SCREEN cCRE Registry V4 at 129 MB** (2,348,854 elements — the
   highest value-per-byte file in this report). Licence: *"without restrictions."*
3. **HT-SELEX: Jolma 2013 (`PRJEB3289`, 10.76 GB)** for the cheap version, or **Codebook
   GHT-SELEX (`PRJEB76622`, 235.71 GB raw / 73.3 MB peaks)** for the 2026 one that uses
   real genomic DNA with native flanks. The only quantitative (TF, sequence) → enrichment
   signal, which is exactly the shape a partner-prediction JEPA head consumes.

**Total under 200 GB** against a 298 TB scratch and the README's 10 TB cap. Volume is not
the constraint on this arm; supervision density is.

### 5.3 The correction that matters: AlphaGenome cannot be *in* the model

AlphaGenome is the best-matched functional encoder in this report — 450M, 3072-d @128 bp
over 1 Mb, 40.8 GB on one H200, plain torch via the port. **And its terms forbid exactly
what we would use it for:** the Model Terms state its outputs *"should not be used for the
training of other machine learning models."* A JEPA head trained on AlphaGenome embeddings
is precisely that. Not as features, not as distillation targets, not as pseudo-labels.

So AlphaGenome is an **external comparator only** — a number we report our model against,
never a component of it. That is a demotion from where this file had it an hour ago, and
it is the single most consequential licence fact in the report.

### 5.4 Alignment with README revision 1

The protein/ligand reconnaissance found that concatenating unimodal encoders loses to a
count fingerprint, and that what survives is a *jointly computed* representation
(Boltz-2's `z`). The same discipline applies here, and the DNA arm has an exactly
analogous trivial control: **a k-mer count vector**. Before any DNA embedding is believed
it must beat k-mer counts + LightGBM on a leave-one-target-cluster-out split.

This is not a formality. §1.2 records that published benchmarking finds general-purpose
DNA LMs lose to specialised models on gene expression and QTLs; §4.8 records that GUE, the
NT suite and Genomic Benchmarks overlap so heavily that a good score on all three is **one
result reported three times**; and several of the standard benchmarks (BEND
`enhancer_annotation` n=285, Vex-seq n≈2,000, `dummy_mouse_enhancers` n=1,210) are small
enough that any delta is inside the permutation null — the FINDING-007 problem verbatim.
**Measure the null first, on the specific benchmark, before reporting any gain.**

### 5.5 Three things that already exist and change the plan

- **JEPA-DNA** (NVIDIA, arXiv 2602.17162) already ran this experiment on DNA across five
  backbones. Start from its published baseline, not from zero.
- **HiCFoundation** (Nat Methods 2026, **Apache-2.0 weights on HuggingFace**) is a
  pretrained Hi-C encoder whose epigenomic head already maps Hi-C→ATAC/ChIP — a partial
  shared-embedding model, free.
- **Evo2HiC** (bioRxiv 2025.11.18.689171) is someone building the sequence+Hi-C half of
  this thesis, and they already hit the "genome-wide sequence embeddings are too
  expensive" wall. Read it before committing architecture.

---

## BLOCKED — needs the user

Ordered by cost if it stays blocked. All are obtainable; none is a technical problem.

### B1. AlphaGenome — accept the gate, **and accept that it can only be a comparator**
`google/alphagenome-all-folds` is **gated**: log in, complete the fields, Accept. Two
clauses matter, both verified on the model card:

> derivatives inherit the identical non-commercial terms, **and training a new model on
> AlphaGenome's outputs or predictions is itself restricted** — the terms state outputs
> *"should not be used for the training of other machine learning models."*

That forecloses AlphaGenome as a component of the shared embedding space (§5.3). Accepting
the gate is still worth doing — it lets us *evaluate against* AlphaGenome, which is the
strongest published baseline on 24 of 26 variant-effect evals. **Decide explicitly whether
we accept comparator-only status or drop it entirely.**

*Honest note on a laundering hazard:* the community port `gtca/alphagenome_pytorch`
carries the same weights and is **not** HF-gated. Downloading from the mirror does not
change what the terms permit.

### B2. NTv3 — accept the gate (one click)
Every `InstaDeepAI/NTv3_*` repo is **gated** behind a **custom non-commercial licence**
("Licensed Models are only available under this License for Non-Commercial Purposes"), and
the HF licence field reads `other`/NOASSERTION. **Both the weights and the NTv3 benchmark
show "Access restricted."** One acceptance covers the collection. It is the only open 1-Mb
single-nucleotide 1536-d encoder, and the `_post` checkpoints carry ~16,000 functional
tracks across 24 species.

### B3. DART-Eval — Synapse account + Data Use Agreement
`syn59522070` requires **registration and a DUA**; anonymous REST returns metadata only,
and the GitHub repo has **no LICENSE file**. Its **task 5 (African caQTLs + Yoruban
dsQTLs) is the only held-out-population variant-effect test in this entire report** — the
least ENCODE-recycled signal available, and therefore the most honest gate we could use
per README kill-criterion 3. Worth the registration.

### B4. AlphaFold3 weights — a written application, not a click
Gated, non-commercial, **requires an application to Google DeepMind with institutional
details**; turnaround days to weeks. **Boltz-2 (MIT) and Chai-1 (Apache-2.0, code *and*
weights) do the same job with no gate**, so this is optional. Start the application only
if there is a specific reason AF3 must be in the comparison.

### B5. Four downloads that block on bot-detection or a bad cert, not on permission
None needs an account; they need a browser or a real user-agent.
- **BioLiP2** — `zhanggroup.org/BioLiP/` returns **HTTP 403** to non-browser clients.
  Mirror: `aideepmed.com/BioLiP/`.
- **DeepPBS dataset** — figshare `10.6084/m9.figshare.25678053` **403s to bots**; fine from
  a browser. Our best structure→PWM supervision.
- **GraphBind benchmark sets** — `csbio.sjtu.edu.cn` serves a **bad TLS certificate**;
  `curl -k`, or use the GraphSite GitHub mirror (MIT).
- **4D Nucleome portal** — `@@download` **403s anonymously**; the **`s3://4dn-open-data-public`
  mirror is open** and serves identical bytes. Use S3.

### B6. Two requester-pays buckets — this one costs actual money
`gs://basenji_barnyard/data` (Enformer training tensors) and `gs://borzoi-paper/data`
(Borzoi tensors) are **requester-pays** — confirmed HTTP 400 anonymously. A 5 TB pull is
**~$600 of GCP egress on a billing project we would have to supply.** We do **not** need
them: the weights are free, and the eQTL benchmark slice is available free through the
Genomics Long-Range Benchmark. **`gs://basenji_hic` (Akita, 10.71 GB) is NOT
requester-pays.** Only spend here if we decide to retrain rather than freeze.

### B7. ProNAB bulk file — possibly one email
`web.iitm.ac.in/bioinfo2/pronab/` exposes a `#download` anchor, but one read of the page
says *"If you would like to access the entire dataset, please contact us."* **14,606
protein–DNA affinity measurements** is the actual regression target for a binding JEPA. If
the anchor does not yield the full table, one email to the IIT-Madras group is worth
sending. Last updated **February 2023**.

### B8. A decaying resource worth mirroring *now*, before it needs a decision
**Roadmap Epigenomics.** `roadmapepigenomics.org` is dead (HTTP 000). **1.2 TB of
irreplaceable signal — 127 consolidated epigenomes, 1,032 p-value bigWigs across 32 marks
— sits on one unmirrored university web server (`egg2.wustl.edu`) with no
machine-readable licence.** ENCODE hosts 2,763 Roadmap experiments but does not reproduce
all 127 epigenomes. Same story, smaller stakes: **GTRD** (serving 2020 files),
**Cistrome** (`DEBUG=True` in production), **ReMap** (frozen since 2021), **UniPROBE**
(homepage soliciting a volunteer student). If any of these is wanted, the decision is
"mirror this month", not "download later".

### B9. Two closed pools, for completeness
- **JGI IMG/M** requires JGI SSO — **~30 × 10¹² bp, larger than all of RefSeq, and the
  biggest closed pool in genomics.** It is why Evo 2 used MGD/IMG-VR extracts instead.
- **GTEx individual-level data** (phs000424.v11.p2) needs a dbGaP DAR plus an institutional
  signing official: **budget 2–4 months.** The summary-level data is open and is almost
  certainly enough.

### B10. The licensing policy call — a decision only the user can make
If the program must stay commercially clean, these drop out: **AlphaGenome and anything
trained on its outputs**, NTv3, every Nucleotide Transformer (CC-BY-NC-SA — *share-alike
arguably reaches derived embeddings*), SegmentNT, ChatNT, Sei, DeepSEA/Beluga/ExPecto,
Orca, GET, the NVIDIA JEPA-DNA checkpoints, **ReMap (CC-BY-NC-4.0)**, **DNAproDB data
(CC-BY-NC)**, **Gosai 2024 MPRA (CC-BY-NC)**, and **the Genomics Long-Range Benchmark
(CC-BY-NC-SA — share-alike contaminates)**.

Separately, a long list of resources have **no licence at all**, which is not the same as
permissive — "downloadable" ≠ "you may redistribute derived embeddings": CIS-BP, GTRD,
Cistrome, SwissRegulon, Codebook, the DART-Eval repo, the official GUE data, the NT task
datasets, DNALongBench, MFASS, Roadmap.

**What survives clean:** Evo 2 + OpenGenome2 (Apache-2.0), Enformer and borzoi-pytorch
(CC-BY-4.0), Flashzoi (MIT), Caduceus / HyenaDNA / PlantCaduceus, ModernGENA, gLM2,
GenomeOcean, ChromBPNet (MIT), DeepPBS and Deep DNAshape (BSD-3), HiCFoundation and Akita
(Apache-2.0), Orca resources (CC-BY-4.0), the PDB (CC0), **ENCODE ("without
restrictions")**, JASPAR (CC-BY-4.0), HOCOMOCO (WTFPL), ChIP-Atlas (CC-BY-4.0), FANTOM5
(CC-BY-4.0), BEND data (CC-BY-4.0), Genomic Benchmarks (Apache-2.0), gnomAD (MIT), and
ClinVar/GEO/GenBank/RefSeq/1000G/HPRC (public domain). **That is a complete stack.**
Knowing which one we are building is worth more than any individual model.

### B11. Two environment limits that cost this survey coverage
- **`CLAUDE_CODE_MAX_WEB_SEARCHES_PER_SESSION` hit 200/200** partway through, for all three
  agents. Everything after that was direct HTTP against primary APIs — which is *why* most
  numbers here are measured rather than quoted, but it did cost breadth. Not covered:
  **Delphi**, DeepBind / BindSpace / ProBound zoos (`proboundmodels.org` does not resolve),
  Kipoi model-zoo liveness, and dedicated protein–DNA *co-embedding* models.
- **The concurrent-subagent cap** was hit while fanning out. Raising both closes those gaps
  in one pass.

---

## Provenance

Compiled 2026-09-20 by the DNA/genomics lead plus two parallel reconnaissance agents
(structure data; regulatory models). Conventions used throughout:

- **✓** = queried live today against a primary API (RCSB Search, ENCODE portal,
  HuggingFace config/file listings), not recalled from training.
- ***(measured)*** = file size obtained by sampling real files, not quoted from a paper.
- **⚠️** = unverified, stale, inflated, or licence-ambiguous.
- **❌** = dead, or incompatible with the torch-2.5.1 / no-`cuequivariance` stack.

Numbers without a mark are quoted from the cited source and were not independently
checked. The PDB counts in §3.1 were derived twice, independently, and agreed exactly.

**Sections 3 and 4 are the best-evidenced parts of this file.** Nearly every count and
size in them is a live API response or a real byte sum taken on 2026-09-20 — ENA Portal,
ENCODE portal, RCSB Search, JASPAR REST, GEO/NCBI FTP, GTDB, public GCS/S3, Zenodo, and
HuggingFace tree-size. The PDB counts in §3.1 were derived twice, independently, and
agreed exactly. Section 1's parameter counts and licences were read from live
`config.json` files and model cards; section 2's port-fidelity and H200-memory claims are
the weakest, since they come from the porting authors themselves and are unreplicated.

**Known gaps.** The session's WebSearch budget was exhausted at 200/200 partway through,
for all three agents; the remainder was direct HTTP. Not covered for want of search
breadth: **Delphi**, DeepBind / BindSpace / ProBound model zoos (`proboundmodels.org` does
not resolve), Kipoi model-zoo liveness, and dedicated protein–DNA *co-embedding* models.
