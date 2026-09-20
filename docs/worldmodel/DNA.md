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
- **AlphaGenome's weights are downloadable now — it is no longer API-only** — and there is
  a validated PyTorch port with an explicit `encode()` giving (B, 1024, 3072) at 128 bp
  over 1 Mb, running in 40.8 GB on one H200 with **no exotic kernels**. It is
  non-commercial, and the terms restrict training on its outputs.
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
  README's 10 TB cap. Data volume is not the constraint on this arm; supervision density is.

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

No coordinates here — these are *binding and activity* readouts. They are three to four
orders of magnitude larger than section 3 and are where the DNA arm's data volume
actually comes from.

**Counts marked ✓ were queried live against the ENCODE portal API on 2026-09-20.**

### 4.1 Table

| # | Dataset | Ver | Yr | URL | Open? | License | Records | Size | Use in the shared space |
|---|---|---|---|---|---|---|---|---|---|
| E1 | **ENCODE** | ENCODE4 | 2026 | [encodeproject.org](https://www.encodeproject.org) | yes, no registration | **fully open** | **27,043 released experiments** ✓ — of which **TF ChIP-seq 5,002** ✓, **histone ChIP-seq 3,832** ✓, **DNase-seq 3,490** ✓, **ATAC-seq 559** ✓ | narrowPeak BEDs for all TF ChIP: **~20–50 GB**; bigWigs: **~10–30 TB**; FASTQ: **petabyte-scale** | The primary DNA-side supervision. **Pull peaks, not signal.** Peak BEDs give (sequence window, bound/unbound, which TF) triples — the exact form a partner-prediction JEPA head consumes. |
| E2 | **ENCODE SCREEN / cCRE Registry** | **V4** | 2026 | [screen.wenglab.org](https://screen.wenglab.org/) | yes | open | **2,348,854 human cCREs / 1,888 cell types** ✓ (1,718,669 enhancers, 47,532 promoters); **926,843 mouse / 366 cell types** | **129.1 MB human, 50.6 MB mouse** ✓ | A ready-made regulatory-element vocabulary — the candidate set for "what binds here". Tiny and high-value. |
| E3 | **JASPAR** | **2026 (11th release)** | 2026 | [jaspar.elixir.no](https://jaspar.elixir.no/downloads/) | yes | **CC-BY-4.0** | ~2–3k CORE profiles (non-redundant + redundant variants) | **<100 MB** | The cheapest possible DNA-side "label": a PWM per TF. Use as a *baseline* the learned embedding must beat, not as a feature. |
| E4 | **HOCOMOCO** | **v14** | 2023/24 | [hocomoco14.autosome.org](https://hocomoco14.autosome.org/) | yes | **WTFPL** ("treat as CC-BY") | **1,595 models**; 1,107 human + 809 mouse TFs | <100 MB | Same role as E3, better coverage. ⚠️ built from **ChIP-Seq (via GTRD) + HT-SELEX + GHT-SELEX + SMiLE-Seq + PBM** — i.e. it is a *meta-analysis of E1/E5/E7*, not independent data. |
| E5 | **HT-SELEX (Jolma 2013)** | — | 2013 | [ENA `PRJEB3289`](https://www.ebi.ac.uk/ena/browser/view/PRJEB3289) — "DNA-binding specificities of human transcription factors" ✓ | yes | ENA open | ~500 human TFs, hundreds of millions of reads | raw FASTQ **~10² GB** (est.) | The only high-throughput assay giving a *quantitative binding curve per TF over random sequence space*. Ideal JEPA training signal: (TF, sequence) → enrichment. |
| E6 | **HT-SELEX / methyl-HT-SELEX (Yin 2017), GHT-SELEX (Codebook, 2024–25)** | — | 2017–25 | GEO/ENA (per-paper accessions) | yes | open | ~540 TFs (Yin); ~hundreds more (Codebook) | ~10² GB | Adds methylation-conditioned specificity — a dimension no motif DB carries. |
| E7 | **CIS-BP** | **v3.10, 2026-04-26** ✓ | 2026 | [cisbp.ccbr.utoronto.ca](https://cisbp.ccbr.utoronto.ca/) | yes | not stated ⚠️ | **13,030 motifs; 169,272 TFs — but only 4,989 from direct experiment** ✓; 741 species; 321 DBD families | <1 GB | ⚠️ **The 169k figure is 97% inferred by DBD similarity, not measured.** Use the 4,989 direct set; treat the rest as a prior, never as a label. |
| E8 | **UniPROBE (PBM)** | — | 2015+ | [uniprobe.org](http://the_brain.bwh.harvard.edu/uniprobe/) | yes | academic | ~700 TFs | <10 GB | Orthogonal assay (universal PBM). Small, and substantially absorbed into E4/E7. |
| E9 | **ReMap** | **2022 (4th)** | 2022 | [remap.univ-amu.fr](https://remap.univ-amu.fr/) | yes | ⚠️ **CC-BY-NC-4.0** | **1,210 regulators, 8,103 ChIP-seq datasets, 182M peaks** ✓ | ~10–30 GB BED | Broader than ENCODE alone. ⚠️ **Assembled from GEO + ENCODE + ENA** — ENCODE is a *subset*, so ReMap ∪ ENCODE double-counts. Pick one as primary. |
| E10 | **ChIP-Atlas** / GTRD / Cistrome DB | ongoing | 2024–26 | [chip-atlas.org](https://chip-atlas.org/) · gtrd.biouml.org · cistrome.org | yes | **ChIP-Atlas CC-BY-4.0** ✓ (permissive — unlike ReMap) | ChIP-Atlas: **433,000 ChIP-seq/ATAC-seq/Bisulfite-seq experiments** ✓ | 10–100 GB peaks | ⚠️ **These three plus E9 are four reprocessings of the same public GEO/SRA deposits.** Their union is not 4× the data. HOCOMOCO's ChIP input comes from GTRD, closing the loop. Choose **one** aggregator. |
| E11 | **4D Nucleome (Hi-C, Micro-C)** | ongoing | 2026 | [data.4dnucleome.org](https://data.4dnucleome.org/) | yes, free | open (coordinated-publication courtesy) | Hi-C, Micro-C, chromatin tracing, microscopy | 10–100 TB if taken whole | 3D *genome* contact, not molecular structure. Relevant only if we add the Orca/Akita contact modality. **Low priority** — it does not touch protein–DNA interfaces. |
| E12 | **MPRA / reporter assays** (Sharpr-MPRA, lentiMPRA, Agarwal 2025, Gosai 2024) | — | 2016–25 | GEO per-paper | yes | open | 10⁵–10⁶ sequences each | 1–50 GB | Direct sequence→activity labels at scale. The best *quantitative* regression target available for a DNA tower. |
| E13 | **FANTOM5 CAGE / Roadmap / GTEx** | — | 2014–24 | fantom.gsc.riken.jp · gtexportal.org | yes | open (GTEx needs dbGaP for individual-level) | 10³ samples | 10–100 GB (summary) | These are Enformer's and Borzoi's *training targets*. Redundant if we use those models frozen. Needed only if we retrain. |
| E14 | **OpenGenome2** | — | 2025 | [HF](https://huggingface.co/datasets/arcinstitute/opengenome2) | yes, ungated | **Apache-2.0** | **8.8 trillion bp**; GTDB v220, IMG/PR, IMG/VR, MGD, NCBI, Ensembl, RNAcentral, Rfam, EPD | **5.52 TB** ✓ | Evo 2's pretraining corpus. We do **not** need it to use Evo 2 frozen — and at 5.52 TB it eats half the README's 10 TB acquisition cap. Pull only if we fine-tune. |
| E15 | **NT downstream-task benchmark (revised)** | — | 2024 | [HF](https://huggingface.co/datasets/InstaDeepAI/nucleotide_transformer_downstream_tasks_revised) | yes | not stated ⚠️ | **532,064 rows, 18 tasks**, chromosome-held-out ✓ | **591 MB** ✓ | The standard evaluation. Chromosome-held-out splits make it a legitimate difficulty-matched test per README kill-criterion 3. |
| E16 | **Genomic Benchmarks** | — | 2023 | [HF `katarinagresova/Genomic_Benchmarks_*`](https://huggingface.co/datasets/katarinagresova/Genomic_Benchmarks_human_enhancers_cohn) | yes | not stated ⚠️ | e.g. human_enhancers_cohn **27,791 rows** ✓ | **6.65 MB** per task ✓ | ⚠️ **Too small to distinguish models.** At 7 MB and ~28k examples, differences between strong DNA FMs here are inside the noise — the same trap as this repo's +0.0138 selector noise floor. Report it, don't decide on it. |
| E17 | **GUE**, **BEND**, **DART-Eval**, **GenBench**, **NABench** | — | 2023–25 | [DART-Eval](https://github.com/kundajelab/DART-Eval) · [NABench](https://arxiv.org/abs/2511.02888) | yes | varies | 20–60 tasks each | <10 GB each | Overlapping benchmark suites — GUE and the NT suite share tasks. **DART-Eval is the one with a real negative result**; use it as the honest gate. |
| E18 | **ClinVar (non-coding) / gnomAD** | ongoing | 2026 | ncbi.nlm.nih.gov/clinvar · gnomad.broadinstitute.org | yes | open | 10⁶ / 10⁸ variants | 1–500 GB | Zero-shot variant-effect evaluation for a frozen DNA encoder — the one task where DNA LMs demonstrably hold up. |

### 4.2 Notes

- **The big redundancy, stated plainly.** ENCODE (E1) → GEO/SRA → reprocessed four times
  as ReMap, ChIP-Atlas, GTRD and Cistrome (E9, E10) → meta-analysed into HOCOMOCO (E4)
  and CIS-BP (E7). Taking all of them is not more data; it is the same experiments at
  five levels of processing, with correlated errors. **Take ENCODE peaks as primary and
  exactly one aggregator for coverage beyond it — and make that aggregator ChIP-Atlas
  (CC-BY-4.0, 433k experiments) rather than ReMap (CC-BY-NC-4.0, 8,103 datasets):
  more data, no licence problem.** That alone cuts the acquisition budget
  by most of an order of magnitude.
- **Licences to watch on the data side.** ReMap is **CC-BY-NC-4.0**, which is the only
  non-commercial restriction in this whole section — ENCODE, JASPAR, HOCOMOCO and the PDB
  are all clean. If ReMap's NC terms matter, ChIP-Atlas or GTRD covers the same ground.
- **CIS-BP's headline number is inflated ~34×.** 169,272 TFs with a motif, of which
  **4,989** come from a direct experiment. Everything else is inferred from DNA-binding-
  domain similarity. Quoting 169k as "TFs with known specificity" would be the exact
  too-clean-number failure this repo has a memory note about.
- **Size reality check against the README's 10 TB cap.** ENCODE peak BEDs (~tens of GB),
  JASPAR/HOCOMOCO/CIS-BP (<2 GB total), HT-SELEX raw (~hundreds of GB), the PDB NA subset
  (single-digit GB), the benchmarks (<10 GB). **A complete, useful DNA corpus is well
  under 1 TB.** The only things that would blow the cap are OpenGenome2 (5.52 TB), ENCODE
  bigWigs (10–30 TB) and 4DN (10–100 TB) — and we need none of the three to run frozen
  encoders. Acquire none of them in wave 2.

---

## 5. Recommendation

**Models — take three, not fifteen.**

1. **Evo 2 7B** (`arcinstitute/evo2_7b`, Apache-2.0, 13.8 GB, 1 Mb ctx, 4096-d). The
   evolutionary-constraint tower. Budget a day for the environment: the repo wants torch
   2.6/2.7 + flash-attn 2.8 + `vortex`, and we have 2.5.1. The 7B is the only size that
   runs bf16 without Transformer Engine.
2. **AlphaGenome via the PyTorch port** (`gtca/alphagenome_pytorch`, 450M, 3072-d @128 bp,
   1 Mb, **40.8 GB peak on one H200**). The function tower. Plain torch, no exotic kernels,
   `encode()` hands back (B, 1024, 3072) directly. Non-commercial — see BLOCKED.
3. **ModernGENA-large** (`AIRI-Institute/moderngena-large`, 377M). The cheap one that will
   actually run on our stack on the first try, pretrained TSS-centred on 443 vertebrate
   genomes. Use it to get the pipeline working before spending H200 hours on Evo 2.

Permissive fallback if non-commercial is ruled out: **Evo 2 + Enformer (CC-BY-4.0) +
borzoi-pytorch (CC-BY-4.0)**. Still a good stack.

**Datasets — three, and they are small.**

1. **PDB protein–DNA subset — 10,733 entries, 4.32 GB, CC0.** 3,027 non-redundant clusters
   after 30% clustering. Plus DNAproDB 2.0 (0.94 GB) for precomputed interface geometry.
2. **ENCODE TF ChIP-seq peak BEDs — 5,002 released experiments, tens of GB, fully open.**
   Peaks only. Not bigWigs, not FASTQ.
3. **HT-SELEX (`PRJEB3289` + Yin 2017 + Codebook GHT-SELEX) — ~10² GB.** The only
   quantitative (TF, sequence) → enrichment signal, which is the exact shape a
   partner-prediction JEPA head consumes.

**Total under 200 GB against a 298 TB scratch and a 10 TB cap.** Data volume is not the
constraint on the DNA arm; supervision *density* is.

**Alignment with README revision 1.** The protein/ligand reconnaissance found that
concatenating unimodal encoders loses to a count fingerprint, and that what survives is a
*jointly computed* representation (Boltz-2's `z`). The same discipline applies here, and
the DNA arm has an exactly analogous trivial control: **a k-mer count vector**. Before any
DNA embedding is believed, it must beat k-mer counts + LightGBM on a
leave-one-target-cluster-out split. Per section 1.2, published benchmarking already says
general-purpose DNA LMs lose to specialised models on gene expression and QTLs — so this
control is not a formality.

---

## BLOCKED — needs the user

Ordered by cost if it stays blocked. All are obtainable; none is a technical problem.

### B1. AlphaGenome weights — accept the terms (one click), *and* make a licensing decision
`google/alphagenome-all-folds` is **gated**: log in, complete the fields, Accept. Licence
is the **AlphaGenome Model Terms — non-commercial only**, and the scope is unusually broad:

> fine-tuned derivatives inherit the identical non-commercial terms, **and training a new
> model on AlphaGenome's outputs or predictions is itself restricted.**

That second clause reaches this program directly: a JEPA head trained on AlphaGenome
embeddings is a restricted derivative. For an academic project that is fine — but it must
be a decision, not a discovery.

*Honest note on a laundering hazard:* the community port `gtca/alphagenome_pytorch`
carries the same weights and is **not** HF-gated. Downloading from the mirror does not
change what the terms permit. Use the official gate.

### B2. NTv3 — accept the gate (one click)
Every `InstaDeepAI/NTv3_*` repo is **gated** behind a **custom non-commercial licence**
("Licensed Models are only available under this License for Non-Commercial Purposes").
One acceptance covers the collection. It is the only open 1-Mb single-nucleotide 1536-d
encoder, and the `_post` checkpoints carry ~16,000 functional tracks across 24 species.
Worth having even if we end up not using it.

### B3. AlphaFold3 weights — a written application, not a click
Gated, non-commercial, **requires an application to Google DeepMind with institutional
details**; turnaround is days to weeks. Only needed if we want AF3 as a protein–DNA
structure source. **Boltz-2 (MIT) and Chai-1 (Apache-2.0, code *and* weights) do the same
job with no gate**, so this is optional — start the application only if there is a
specific reason AF3 must be in the comparison.

### B4. Three downloads that block on bot-detection, not on permission
None of these needs an account; they need a browser or a `curl` with a real UA.
- **BioLiP2** — `zhanggroup.org/BioLiP/` returns **HTTP 403** to non-browser clients.
  Mirror: `aideepmed.com/BioLiP/`.
- **DeepPBS dataset** — figshare DOI `10.6084/m9.figshare.25678053` **403s to bots**;
  downloads fine from a browser. This is our best structure→PWM supervision.
- **GraphBind benchmark sets** — `csbio.sjtu.edu.cn` serves a **bad TLS certificate**;
  `curl -k` works, or use the GraphSite GitHub mirror (MIT).

### B5. ProNAB bulk file — possibly an email
`web.iitm.ac.in/bioinfo2/pronab/` exposes a `#download` anchor, but one read of the page
says *"If you would like to access the entire dataset, please contact us."* **14,606
protein–DNA affinity measurements** is the actual regression target for a binding JEPA, so
if the anchor does not yield the full table, one email to the IIT-Madras group is worth
sending. Note it has not been updated since **February 2023**.

### B6. The licensing policy call — a decision only the user can make
If the program must stay commercially clean, these all drop out: AlphaGenome **and
anything trained on its outputs**, NTv3, every Nucleotide Transformer (v1/v2/SegmentNT/
ChatNT are CC-BY-NC-SA-4.0, and *share-alike* arguably reaches derived embeddings), Sei,
DeepSEA/Beluga/ExPecto, Orca, GET, the NVIDIA JEPA-DNA checkpoints, **ReMap
(CC-BY-NC-4.0)** and **DNAproDB data (CC-BY-NC)**.

What survives clean: **Evo 2, Enformer, borzoi-pytorch/Flashzoi, Caduceus, HyenaDNA,
PlantCaduceus, ModernGENA, gLM2, GenomeOcean, ChromBPNet, DeepPBS, Deep DNAshape, the
PDB (CC0), ENCODE, JASPAR (CC-BY-4.0), HOCOMOCO (WTFPL).** That is a complete stack.
Knowing which one we are building is worth more than any individual model.

### B7. Two environment limits that cost this survey coverage
- **`CLAUDE_CODE_MAX_WEB_SEARCHES_PER_SESSION` was exhausted at 200/200** partway through.
  Everything after that point was done by direct HTTP fetch against canonical URLs, which
  is why most numbers here are *measured* rather than quoted — but it did cost breadth.
  Not covered for lack of search: **Delphi**, DeepBind/BindSpace, Kipoi model-zoo
  liveness, and protein–DNA co-embedding models.
- **The concurrent-subagent cap** was hit while fanning out. Raising both would close
  those gaps in one pass.

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

**Known gaps.** The session's WebSearch budget was exhausted at 200/200 partway through;
the remainder was direct HTTP. Not covered for want of search breadth: **Delphi**,
DeepBind / BindSpace / ProBound model zoos (`proboundmodels.org` does not resolve),
Kipoi model-zoo liveness, and dedicated protein–DNA *co-embedding* models. §4 (functional
genomics) is the least independently cross-checked section — its ENCODE, JASPAR,
HOCOMOCO, CIS-BP, ReMap, ChIP-Atlas and SCREEN figures were fetched directly and carry ✓,
but the HT-SELEX and MPRA size estimates are order-of-magnitude, marked (est.), and should
be measured before any acquisition is budgeted against them.
