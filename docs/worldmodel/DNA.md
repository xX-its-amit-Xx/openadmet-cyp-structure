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
  single-nucleotide, and the only large model that runs bf16 on one GPU *without*
  Transformer Engine. Its 40B sibling does not fit an 80 GB H200 in bf16.
- **NTv3 (650M, 1 Mb, single-nucleotide)** is technically the best-matched encoder and is
  **gated behind a non-commercial licence** — see *BLOCKED*.
- **JEPA-DNA (NVIDIA, 2026)** already did the exact thing this program proposes, for DNA,
  and released weights. Read it before writing any JEPA head.
- Structure data is the scarce half. Protein–DNA co-crystals are ~10⁴ complexes, not 10⁶;
  everything else in section 2 is either the same PDB rewrapped or is *binding* data
  (motifs, ChIP, SELEX) that carries no coordinates.

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

- **The dimensions line up.** AlphaGenome 3072 @128 bp, Enformer 3072 @128 bp, Borzoi
  ~1920 @32 bp. AlphaGenome and Enformer are *directly concatenable at matched 128-bp
  resolution*, which is a real gift for a shared space.
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
| 30–31 | AlphaGenome (+ torch port) | **Primary *functional* DNA encoder.** `encode()` → (B, 1024, 3072) at 128-bp resolution across 1 Mb; concatenate with Enformer at matched resolution. |
| 32 | Enformer | Second functional encoder, 3072-d @128 bp; the permissively-licensed stand-in if AlphaGenome's terms fail. |
| 33–35 | Borzoi / borzoi-pytorch / Flashzoi | ~1920-d @32-bp encoder — the finest-grained functional embedding available, and the only one carrying RNA-seq coverage. |
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


*(pending)*

---

## 4. Functional-genomics / binding data at scale

*(pending)*

---

## BLOCKED — needs the user

*(pending)*
