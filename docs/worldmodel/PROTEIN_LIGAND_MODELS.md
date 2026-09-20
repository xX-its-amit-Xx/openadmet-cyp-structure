# Protein + small-molecule representation models — frozen-embedding dossier

**Compiled 2026-09-20** for the shared-embedding world model (see `README.md` in this
directory). Scope: models whose **frozen** embeddings we can extract on **one H200**,
under **torch 2.5.1**, without installing anything that upgrades torch or compiles CUDA.

**Hard constraints this document is filtered against**

| constraint | consequence |
|---|---|
| torch **2.5.1**, cannot change | anything pinning `torch<2.3`, `torch>=2.7`, or shipping prebuilt kernels against another ABI is **out** |
| 1× H200 (80 GB), contended | ≤7B params in bf16 is fine; no pretraining |
| C: ~34 GB / D: ~20 GB free | weights land on `/scratch` (Explorer) or OneDrive, **never** locally |
| no custom-CUDA builds | `cuequivariance_ops_cu12`, `Uni-Core`, `torch-scatter`/`torch-sparse` wheels, DGL → treat as blockers until proven wheel-available |

Verification method: HuggingFace `config.json` + `/api/models?blobs=true` for dims, sizes,
gating and license; PyPI JSON for dependency pins and last-upload dates; GitHub API for
license and last push. Numbers below are measured, not recalled, except where marked ~.

**The shortlist, if you read nothing else**

| rank | protein | ligand |
|---|---|---|
| 0 (control) | one-hot target / 3Di-bigram over pocket | **count-ECFP4 (2048) + RDKit descriptors + LightGBM** |
| 1 | **Boltz-2 trunk** — `--write_embeddings`, MIT, already running here, s=384 / z=128 (§2) | **Boltz-2 pair block** — the same forward pass; the only pose-conditioned ligand vector |
| 2 | **ESM-2 650M** — MIT, 7.8 GB, 1280-d, zero friction | **MoLFormer-XL** — Apache-2.0, 0.19 GB, 768-d, pure torch |
| 3 | **SaProt 650M** — MIT, 1280-d, loads as `EsmModel`, needs only the `foldseek` binary | **Uni-Mol v1/v2 via `unimol_tools` 0.1.6** — MIT, 512/768-d, **no Uni-Core any more** |

Everything else in §1 and §3 is either a decorrelation experiment, a licence compromise, or a
dependency trap. **§4 is the part that should change behaviour**: the honest literature says most
of these lose to a fingerprint.

---

## 1. Protein models — summary table

Size = total HF repo bytes (usually fp32 + duplicate formats; bf16 download is ~half).
"Emb dim" = per-residue hidden width you get from the frozen trunk.

| Model | Ver / Year | URL | Weights / license | Params | Emb dim | Input | Size GB | H200 / deps | How we'd use it |
|---|---|---|---|---|---|---|---|---|---|
| **ESM-2** | t33_650M_UR50D, 2022 | [hf](https://huggingface.co/facebook/esm2_t33_650M_UR50D) | open, **MIT**, ungated | 650 M | **1280** (33 layers, vocab 33) | AA sequence | 7.82 | trivial; `transformers` `EsmModel`, pure torch | **Default protein axis.** Mean-pool per chain; per-residue for pocket residues. |
| ESM-2 (small) | t12_35M / t30_150M / t6_8M | [hf](https://huggingface.co/facebook/esm2_t12_35M_UR50D) | MIT | 35 M | 480 | AA seq | 0.41 | trivial | Ablation control — how much of the gain is just "a PLM" vs "a big PLM". |
| ESM-2 (3B) | t36_3B_UR50D, 2022 | [hf](https://huggingface.co/facebook/esm2_t36_3B_UR50D) | MIT | 3 B | 2560 | AA seq | 22.72 | fits bf16; 3× slower | Only if 650M plateaus. 15B exists, not worth the disk. |
| **ESM C 600M** | esmc-600m-2024-12 | [hf](https://huggingface.co/EvolutionaryScale/esmc-600m-2024-12) | ungated on HF, **Cambrian NON-COMMERCIAL** | 600 M | 1152 | AA seq | 2.30 | `esm` pypi pkg, pure torch | Better-than-ESM2 per FLOP. **License blocks any commercial path** — see BLOCKED. |
| ESM C 300M | esmc-300m-2024-12 | [hf](https://huggingface.co/EvolutionaryScale/esmc-300m-2024-12) | **Cambrian Open** (permissive, incl. commercial) | 300 M | 960 | AA seq | 1.33 | trivial | The one EvolutionaryScale model with a clean license. Good default if 600M's licence bites. |
| ESM C 6B | 2024-12 | API only (Forge) | **no weights** | 6 B | — | — | — | — | Not usable frozen-offline. Ignore. |
| ESM-3 (open) | esm3-sm-open-v1, 2024 | [hf](https://huggingface.co/EvolutionaryScale/esm3-sm-open-v1) | ungated, **Cambrian NC** | 1.4 B | 1536 | seq + struct tokens + function tracks | 5.50 | `esm` pkg; needs 3Di-ish struct tokenizer for the structure track | Multi-track: can condition on structure. NC licence; treat as research-only. |
| **ProtT5-XL-U50 (enc)** | half_uniref50-enc, 2021 | [hf](https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc) | open, no explicit license tag on repo (paper: academic-friendly) | 1.2 B enc-only | **1024** (24 layers) | AA seq | 2.42 | trivial, `T5EncoderModel` fp16 | The classic non-ESM sequence axis — decorrelated from ESM by training corpus/objective. |
| **ProstT5** | 2023 | [hf](https://huggingface.co/Rostlab/ProstT5) | open, **MIT** | 1.2 B | 1024 | AA seq **↔ 3Di** (bilingual) | 11.28 | trivial | Translates seq→3Di and back; gives a *structure-flavoured* embedding **without a structure**. |
| Ankh-large | 2023 | [hf](https://huggingface.co/ElnaggarLab/ankh-large) | open, **CC-BY-NC-SA-4.0** | 1.15 B | **1536** (48 layers) | AA seq | 7.52 | trivial | Strong small-data embedder; NC licence. |
| **Ankh3-large** | 2025 | [hf](https://huggingface.co/ElnaggarLab/ankh3-large) | open, **CC-BY-NC-SA-4.0** | 1.15 B | 1536, vocab 256 | AA seq, `[NLU]` prefix for embeddings | 7.52 | trivial; **must prepend `[NLU]`** or you get the generation head's space | Newest Ankh; multi-task (denoise + completion). NC. |
| Ankh3-XL | 2025 | [hf](https://huggingface.co/ElnaggarLab/ankh3-xl) | CC-BY-NC-SA-4.0 | ~3 B | 2560 | AA seq | 22.92 | fits | Only if Ankh3-large wins the bake-off. |
| **AMPLIFY 350M** | 2024 | [hf](https://huggingface.co/chandar-lab/AMPLIFY_350M) · [gh](https://github.com/chandar-lab/AMPLIFY) | open, **MIT**, ungated | 350 M | **960** (32 layers) | AA seq | 1.42 | `trust_remote_code=True`; pure torch; repo pushed 2026-05 | Best quality-per-GB. MIT. Nvidia also mirrors it (BioNeMo). Cheap third axis. |
| ProGen2-base | 2022 | [hf](https://huggingface.co/hugohrban/progen2-base) | open, **BSD-3** | 764 M | 1536 | AA seq (decoder-only) | 3.06 | community HF port; original repo is JAX-ish/stale | Autoregressive → gives a **likelihood** per sequence, not just an embedding. Weak as an embedder. |
| ProGen2-xlarge | 2022 | [hf](https://huggingface.co/hugohrban/progen2-xlarge) | BSD-3 | 6.4 B | 4096 | AA seq | 12.89 | fits bf16 | Low priority: decoder-only PLMs underperform encoders on embedding transfer. |
| **SaProt 650M** | SaProt_650M_AF2, 2024 | [hf](https://huggingface.co/westlake-repl/SaProt_650M_AF2) · [gh](https://github.com/westlake-repl/SaProt) | open, **MIT**, ungated | 650 M | **1280**, **vocab 446** (= 20 AA × 20 3Di + specials) | AA **+ Foldseek 3Di** token per residue | 5.21 | ESM-2 architecture → loads with `EsmModel`. **Needs `foldseek` binary** (GPL-3, conda/static binary, no CUDA) | **The structure-aware axis we can actually run.** Feed our Boltz/crystal CYP structures → 3Di → SaProt. Same code path as ESM-2. |
| FoldSeek 3Di tokens | v10, 2026 | [gh](https://github.com/steineggerlab/foldseek) | **GPL-3.0** binary | n/a | 20-letter alphabet | backbone PDB/mmCIF | <0.1 | CPU only, fast; GPL is viral for linked code but we shell out | Free structural alphabet: a bag-of-3Di or 3Di-kmer vector is a *shockingly strong, nearly free* structure baseline. Also the input to SaProt/ProstT5. |
| **ProteinMPNN (encoder)** | v_48_020, 2022 | [gh](https://github.com/dauparas/ProteinMPNN) | open, **MIT** | ~1.7 M | **128** hidden, 3 enc layers | backbone coords (N,CA,C,O) | **0.0067** | trivially runs anywhere, pure torch, no deps | Per-residue *structural-environment* embedding at 6.7 MB. Ideal cheap pocket descriptor; run on the pocket residues of every pose. |
| ESM-IF1 | esm_if1_gvp4_t16_142M_UR50, 2022 | [gh](https://github.com/facebookresearch/esm) | open, **MIT** | 142 M | 512 | backbone coords | ~0.5 | ⚠ needs `torch_geometric` **+ `torch-scatter`/`torch-sparse`** and `biotite`; scatter/sparse need matching prebuilt wheels for torch 2.5.1+cu12x | Inverse-folding embedding. Only if ProteinMPNN's 128-d proves too thin. `fair-esm` PyPI last shipped **2022-11**. |
| GearNet / ESM-GearNet | ICLR'23, weights on [zenodo](https://zenodo.org/record/7593637) | [gh](https://github.com/DeepGraphLearning/GearNet) | MIT, weights open | ~22 M | 3072 (concat of 6×512) | residue graph from PDB | ~0.1 | ❌ **TorchDrug 0.2.1, last released 2023-07, `requires_python >=3.7,<3.11`**, needs torch-scatter/cluster | **Do not attempt.** Dependency stack is dead; cost of resurrection exceeds the value of a 22 M-param encoder. |
| **Boltz-2** | 2.2.1, 2025 | [gh](https://github.com/jwohlwend/boltz) · [hf](https://huggingface.co/boltz-community/boltz-2) | open, **MIT**, ungated | ~1 B trunk | **s = 384**, **z = 128** | seq + ligand SMILES/CCD + MSA | 6.20 (`boltz2_conf.ckpt` + `boltz2_aff.ckpt` + `mols.tar`) | `torch>=2.2` → **2.5.1 OK**. `cuequivariance` is an **optional extra** (`pip install boltz[cuda]`); `--no_kernels` runs without it | **Already running here.** See §2 — a single flag dumps the trunk embeddings. This is our only *joint* protein+ligand representation. |
| Boltz-1 / 1x | 2024 | same repo | MIT | ~0.9 B | s=384, z=128 | same | ~2 | same | Only as a decorrelation check against Boltz-2. |
| Chai-1 | 0.6.1, 2024 | [gh](https://github.com/chaidiscovery/chai-lab) · [hf](https://huggingface.co/chaidiscovery/chai-1) | HF repo tagged **Apache-2.0**, but chai-lab README restricts local weights to **non-commercial**; commercial use needs their server/licence — **contradictory, treat as NC** | ~0.5 B | pair/single trunk (not exposed by CLI) | seq + ligand | 1.18 | **`torch<2.7,>=2.3.1` → 2.5.1 OK**; `numpy~=1.21` pin fights modern stacks | Second co-folding opinion (repo already uses it). No embedding export; needs a forward-hook patch. |
| Protenix / Protenix-v2 | 2025–26 | [gh](https://github.com/bytedance/Protenix) | open, **Apache-2.0** | ~0.4 B | trunk s/z (not exported) | seq + ligand | ~2 | Apache-2.0 is the cleanest licence of any AF3 reimpl | Third co-folding opinion; repo already uses it via OpenProtein. Embeddings need a hook. |
| ESMFold | v1, 2022 | [hf](https://huggingface.co/facebook/esmfold_v1) | open, **MIT** | 690 M folding trunk on ESM-2-3B | 2560 (LM) + 1024 (structure module s) | AA seq, **no MSA** | 8.44 | `transformers` `EsmForProteinFolding`, pure torch | MSA-free folding trunk. Fast per-sequence structural embedding for the *whole* P450 universe (493 pairs / 185 targets) where MSAs are unaffordable. |
| OpenFold | 2.x, 2025 | [gh](https://github.com/aqlaboratory/openfold) | **Apache-2.0**, weights open (AF2 params are CC-BY-4.0) | 93 M | Evoformer s = 384, z = 128 | seq + MSA | ~5 + MSA DBs (TBs) | ⚠ custom **DeepSpeed / triton attention kernels**; build pain is real | AF2 Evoformer embeddings with a *usable* licence. Only if AF2-space is specifically wanted. |
| AlphaFold2 (DeepMind) | 2021 | [gh](https://github.com/google-deepmind/alphafold) | code Apache-2.0, **params CC-BY-4.0** | 93 M | s=384, z=128 | seq + MSA | ~4 + DBs | JAX, not torch — separate env | Reference only; OpenFold gives the same numbers in torch. |
| **AlphaFold3** | 2024 | [gh](https://github.com/google-deepmind/alphafold3) | code open, **weights by individual request, non-commercial, no redistribution** | ~1 B | s=384, z=128 | seq + ligand | ~1 (params) + DBs | JAX | **BLOCKED** — weights require a signed request to DeepMind. Boltz-2 is the drop-in. |

### Protein notes, one to three lines each

- **ESM-2 650M is the control, not the champion.** Everything else must beat it on a held-out
  task before it earns disk. It is MIT, ungated, 7.8 GB, and loads in four lines of
  `transformers`. Use `output_hidden_states=True` — the *last* layer is over-specialised to
  the MLM head; layers ~30–33 (or a mean over the top third) usually transfers better.
- **ESM C vs ESM-2 is a licence decision, not a quality decision.** ESM C 600M beats ESM-2 650M
  on EvolutionaryScale's own evals, but carries the **Cambrian non-commercial** licence. ESM C
  **300M** is under the permissive Cambrian Open licence and is the version to standardise on if
  anything downstream could ever be commercial. ESM C 6B is API-only — no weights, ever.
- **ESM-3's value is that it eats structure,** via a structure-token track. That makes it the
  only sequence-family model that can be conditioned on our predicted CYP3A4 poses. Non-commercial.
- **ProtT5 and Ankh are the decorrelation play.** They are T5 encoders trained on different
  corpora with different objectives, so their errors are less correlated with ESM's than another
  ESM checkpoint would be. Per FINDING-J in the root CLAUDE.md, *decorrelation is the thing that
  has repeatedly failed to materialise* in this project — measure it (CKA / per-target Spearman
  between residuals) before assuming an ensemble helps.
- **Ankh3 requires the `[NLU]` prefix token** to get the representation space rather than the
  sequence-completion space. Easy to get wrong silently.
- **AMPLIFY 350M is the efficiency pick**: MIT, 1.42 GB, 960-d, actively maintained (last push
  2026-05), and its paper's whole claim is matching much larger PLMs. Needs `trust_remote_code`.
- **SaProt is the one structure-aware PLM with no dependency tax.** It is literally ESM-2's
  architecture with a 446-token vocabulary (20 AA × 20 3Di + specials), so `EsmModel` loads it.
  The only new dependency is the **`foldseek` binary** (GPL-3, CPU, static binary available).
  Concretely: our poses → `foldseek structureto3didescriptor` → interleave with the sequence →
  SaProt. This is the highest-value protein item after Boltz-2 because it lets the *predicted
  pose* change the protein embedding.
- **3Di on its own is a free baseline.** Before SaProt, try a 400-dim 3Di-bigram count vector
  over pocket residues. It costs milliseconds and is a legitimate control for "does the fancy
  structure PLM add anything over the structural alphabet it consumes?"
- **ProteinMPNN's encoder is 6.7 MB and gives 128-d per-residue structural context.** It is
  almost free and is an honest "structure only, no evolution" axis. Pair it with the pocket
  residue list we already compute in `cypstruct.targets`.
- **ESM-IF1 and GearNet are dependency traps.** ESM-IF1 needs `torch-scatter`/`torch-sparse`
  (prebuilt wheels must match torch 2.5.1 + CUDA exactly, else a 40-minute source build).
  GearNet needs **TorchDrug 0.2.1**, which caps Python at <3.11 and was last released in
  **July 2023**. GearNet is a *no* — 22 M params is not worth resurrecting a dead stack.
- **AF3 is blocked and it does not matter.** Boltz-2 is MIT, runs here already, and exports
  embeddings with a flag. AF3 needs a signed non-commercial request, JAX, and cannot be shared.

---

## 2. Boltz-2 trunk embeddings — the concrete recipe

**Answer: yes, and it is a one-flag operation. We do not need to patch anything.**

`boltz` 2.2.1 (`src/boltz/main.py:1038`) exposes:

```
--write_embeddings    " to dump the s and z embeddings into a npz file. Default is False."
```

and `src/boltz/data/write/writer.py:249-257` does exactly:

```python
# Save embeddings
if self.write_embeddings and "s" in prediction and "z" in prediction:
    s = prediction["s"].cpu().numpy()
    z = prediction["z"].cpu().numpy()
    path = struct_dir / f"embeddings_{record.id}.npz"
    np.savez_compressed(path, s=s, z=z)
```

| item | value | source |
|---|---|---|
| command | `boltz predict <yaml_dir> --out_dir <out> --model boltz2 --write_embeddings` | `main.py` |
| file written | `<out>/predictions/<id>/embeddings_<id>.npz`, keys `s`, `z` | `writer.py` |
| **single** repr `s` | `[n_tokens, 384]` (`token_s: 384`) | `scripts/train/configs/structure.yaml:81` |
| **pair** repr `z` | `[n_tokens, n_tokens, 128]` (`token_z: 128`) | `structure.yaml:82` |
| atom-level dims | `atom_s: 128`, `atom_z: 16` | `structure.yaml:79-80` |
| trunk depth | Pairformer, `num_blocks: 48` in the shipped structure config (the *boltz-as-FM* paper describes a 64-layer trunk — **check the checkpoint you actually load**) | `structure.yaml:108` |
| torch | `torch>=2.2` → **2.5.1 fine** | `pyproject.toml` |
| `cuequivariance` | **optional extra only** (`boltz[cuda]`); `--no_kernels` disables. **Do not install it** — the `cuequivariance_ops_cu12` wheels are the thing that would fight torch 2.5.1 | `pyproject.toml` `[project.optional-dependencies] cuda` |
| weights | `boltz-community/boltz-2` on HF, **MIT**, ungated, 6.20 GB total (`boltz2_conf.ckpt`, `boltz2_aff.ckpt`, `mols.tar`) | HF API |
| extra artefact | `pre_affinity_<id>.npz` is also written — the affinity-module input representation | `writer.py:179` |

**Caveats that will bite:**

1. `z` is **O(n²·128)**. CYP3A4 is ~480 residues; with heme + ligand that is ~500 tokens →
   500·500·128 float32 = **128 MB per pose, before compression**. At 20 poses × 87 ligands
   that is ~220 GB of `z`. **Pool on the fly, never archive raw `z`.** `s` is 500×384 = 0.77 MB
   and is cheap to keep.
2. The published recipe that works (*Boltz is a Strong Baseline for Atom-level Representation
   Learning*, [arXiv:2602.13249](https://arxiv.org/html/2602.13249v1), code
   [hsjang0/boltz-as-FM](https://github.com/hsjang0/boltz-as-FM)) uses the **pair** repr, not the
   single one: take layers **{16, 32, 48, 64}** of the Pairformer, concatenate → 512-d per pair,
   then "hybrid pooling" = mean+std statistic pooling over (diagonal entries ∥ bonded atom-pair
   entries ∥ all entries) → a **3072-d** fixed-length vector. The shipped `--write_embeddings`
   gives you only the **final** layer, so reproducing their multi-layer version needs a hook on
   `pairformer_module`. Start with the final layer; escalate only if it underperforms.
3. That paper's headline: Boltz-2 representations were best on **9 of 22 TDC ADMET tasks** among
   non-ensemble methods, and showed **low CKNNA alignment with existing small-molecule models** —
   i.e. it is *complementary* to Uni-Mol/MolFormer rather than redundant. It did **not** report a
   head-to-head against ESM or ECFP4, and the repo README says "code is currently being
   refactorized" (17 stars, no licence file). Treat as a promising recipe, not a proven one.
4. For the **ligand** side specifically, the sub-block of `z` indexed by ligand tokens (and the
   protein×ligand cross block) is a genuinely *interaction-aware* ligand embedding — something no
   ligand-only model can give. This is the single most differentiated asset we have, because we
   already pay for the forward pass.
5. `--write_embeddings` runs the full structure prediction. If we only want representations, the
   marginal cost over what the pose pipeline already does is ~zero; standalone it is ~1–3 min/pose
   on an H200 with MSA precomputed.

**Recommendation:** add `--write_embeddings` to `scripts/cofold/modal_boltz.py` now, keep `s`,
pool `z` into the three blocks (diag, ligand×ligand, protein×ligand) at write time, and discard
raw `z`. Every pose we generate from here on then also produces a joint embedding for free.

---

## 3. Small-molecule / ligand models — summary table

Sizes and dims below are read from HF `config.json`, the HF blob API, PyPI `requires_dist`, or
the repo's own model zoo table. `~` means the source did not state it.

| Model | Ver / Year | URL | Weights / license | Params | Emb dim | Input | Size GB | H200 / deps | How we'd use it |
|---|---|---|---|---|---|---|---|---|---|
| 🎯 **ECFP4 / Morgan (count)** | RDKit 2026.03 | [rdkit](https://www.rdkit.org/) `rdFingerprintGenerator.GetMorganGenerator(radius=2)` | **BSD-3**, `pip install rdkit` | 0 | **2048** bits (1024/4096 also standard); use **counts**, not bits | SMILES → 2D graph | 0.0 | zero deps, no torch, 10⁴–10⁵ mol/s/core | **The control.** Nothing ships unless it beats this on target-held-out folds. |
| 🎯 **RDKit 2D descriptors** | 2026.03 | `rdkit.Chem.Descriptors.descList` | BSD-3 | 0 | **~210–220** (count `len(descList)` in-env; it drifts by version) | SMILES | 0.0 | zero | Second control, and the *direct ablation* for CheMeleon (which is descriptor distillation). |
| ⭐ **MoLFormer-XL-both-10pct** | 2022 paper / HF port | [hf](https://huggingface.co/ibm-research/MoLFormer-XL-both-10pct) · [gh](https://github.com/IBM/molformer) | open, **Apache-2.0**, ungated, 160k dl/mo | **46.8 M** | **768** (12 layers, 12 heads, vocab 2362, max_pos 202) | SMILES | **0.19** | 🟢 lowest risk: HF port is pure `torch.nn`, no `fast_transformers`, no CUDA ext | `AutoModel(trust_remote_code=True, deterministic_eval=True)` → `pooler_output` [B,768]. First learned ligand axis to try. |
| ⭐ **Uni-Mol v1** (via `unimol_tools`) | 2022/23; tools **0.1.6, 2026-06** | [gh](https://github.com/deepmodeling/Uni-Mol) · weights [hf](https://huggingface.co/dptech/Uni-Mol-Models) | open, **MIT**, ungated | ~47.6 M | **512** CLS + per-atom reprs | SMILES → RDKit ETKDG **3D conformer** | 0.19 (`mol_pre_no_h_220816.pt` = 190 MB) | 🟢 `unimol_tools` 0.1.6 `requires_dist` **verified**: `torch>=2.4.0` floor, **no Uni-Core, no CUDA ext, no apex** | `UniMolRepr(model_name='unimolv1').get_repr(smiles, return_atomic_reprs=True)` → `cls_repr` + atom reprs aligned to `mol.GetAtoms()`. **Only clean 3D-aware ligand route.** |
| **Uni-Mol2** | 2024, NeurIPS | [hf](https://huggingface.co/dptech/Uni-Mol2) | open, **MIT**, ungated | 84 M / 164 M / 310 M / 570 M / **1.1 B** | **768** (84M, 164M) · **1024** (310M) · **1536** (570M, 1.1B); pair dim 512 | SMILES → 3D conformer, two-track atom+pair transformer | 0.34 / 0.66 / 1.24 / 2.28 / **4.52** (zoo total 9.04) | 🟢 **via `unimol_tools`** (`model_name='unimolv2', model_size='84m'…'1.1B'`); 🔴 via the `unimol2/` repo folder, which hard-requires **Uni-Core** | Same `get_repr()` API. Start at 84m/164m — and note §4: Uni-Mol2 ranked **below ECFP** in the largest scaffold-split benchmark. |
| Uni-Mol+ | 2023-24 | [gh](https://github.com/deepmodeling/Uni-Mol/tree/main/unimol_plus) | MIT, release ckpts | 27.7/52.4/77 M | not documented | 2D graph → refined 3D | 0.1–0.3 | 🔴 hard-requires Uni-Core | **Don't.** Ships no embedding extractor; only PCQM4Mv2 HOMO–LUMO inference. |
| **ChemBERTa-2 (77M-MTR)** | 2022 | [hf](https://huggingface.co/DeepChem/ChemBERTa-77M-MTR) | open, **no license tag on HF** (legally unspecified), ungated | tiny: **3 layers**, hidden 384 | **384** | SMILES | **0.014** | 🟢 none; arch class is `RobertaForRegression` (non-stock) but `AutoModel` loads the encoder | MTR variant (multi-task regression on 200 RDKit descriptors) usually embeds better than MLM. ⚠ **"77M" is the *corpus* size, not parameters.** |
| ChemBERTa-2 (77M-MLM) | 2022 | [hf](https://huggingface.co/DeepChem/ChemBERTa-77M-MLM) | no license tag | 3 layers / 384 | 384 | SMILES | 0.014 | 🟢 none | MLM sibling of the above. Cheap ablation pair. |
| ChemBERTa v1 | 2020 | [hf](https://huggingface.co/seyonec/ChemBERTa-zinc-base-v1) | **no license tag**, ungated, 139k dl/mo | ~44 M (RoBERTa 6L/12H) | **768** | SMILES | 0.36 | 🟢 stock `transformers` | Historical baseline. CLS or mean-pool. |
| ⭐ **CheMeleon** | 2025 | [gh](https://github.com/JacksonBurns/chemeleon) · Zenodo [10.5281/zenodo.15426600](https://doi.org/10.5281/zenodo.15426600) | open, **MIT**, ungated | small D-MPNN | not documented | SMILES / RDKit Mol | <0.1 | 🟢 via `chemprop>=2.2.0`; chemprop 2.3.1 pins only `torch>=2.1`, `lightning>=2.0` — **no PyG, no torch-scatter, nothing compiled** | `CheMeleonFingerprint()` is a real public embedding class. Pretrained to predict **Mordred descriptors** — i.e. descriptor distillation, not contrastive SSL. Best pain-to-value ratio of the 2025 graph models. ⚠ name collision: `hspark1212/chemeleon` is **crystal** diffusion, unrelated. |
| **ChemFM-1B / 3B** | 2024-25 | [hf](https://huggingface.co/ChemFM/ChemFM-3B) · [gh](https://github.com/TheLuoFengLab/ChemFM) | ⚠ **license conflict**: HF tag says `mit`, GitHub README says **CC-BY-NC-4.0**. Resolve before use. Ungated | **970 M** / **3.00 B** | **2048** / **3072** | SMILES, `max_position_embeddings` 512, vocab 320 | ~2 / ~6 | 🟢 stock `LlamaForCausalLM` | Decoder-only → `output_hidden_states` + mean/last-token pool. Also ~20 task-specific `ChemFM/admet_*` heads. |
| GP-MoLFormer-Uniq | 2025 | [hf](https://huggingface.co/ibm-research/GP-MoLFormer-Uniq) | open, Apache-2.0 | 46.8 M | 768 | SMILES | 0.19 | 🟢 same profile as MoLFormer-XL | Generative/causal. Strictly worse as a *representation* than bidirectional MoLFormer-XL. Off-label. |
| ⚠ **SMI-TED (smi-ted-Light)** | 2024 | [hf](https://huggingface.co/ibm-research/materials.smi-ted) · [gh](https://github.com/IBM/materials) | open, Apache-2.0, ungated | **289 M** (+ an 8×289M MoE) | **768** (12L/12H, max_len 202) | SMILES | 1.16 per file; repo **4.75** (same weights in 4 formats) | 🔴 **highest risk**: `load.py` imports `fast_transformers.*` → **`pytorch-fast-transformers` 0.4.0, last released 2021-04**, sdist only, no wheels, compiles CUDA, predates sm_90 | **Skip under our constraint** unless someone patches it to pure-torch linear attention. |
| MiniMol | 2024 | [gh](https://github.com/graphcore-research/minimol) | open, **MIT**, weights ship with the pip package | **10 M** | **512** | SMILES | <0.1 | 🔴 `requires_dist` **verified**: `graphium==2.4.7` (exact pin) → torch-geometric + torch-sparse + torch-cluster + torch-scatter **and `torchmetrics<0.11`**; README installs PyG wheels for torch 2.3.0+cu124 | Nicest API of the graph set, nastiest install. Isolate in its own venv or skip. It is the model Boltz-2 was ensembled with in arXiv:2602.13249. |
| MolCLR | 2022 | [gh](https://github.com/yuyangw/MolCLR) | MIT; **checkpoints are in-repo** (`ckpt/pretrained_gin`, `pretrained_gcn`) | ~1.8 M (GIN, 5 layers) | 300 node/graph (512 projection head) | 2D graph | <0.05 | 🔴 pins torch 1.7.1+cu110, PyG 1.6.3, torch-sparse 0.6.9, torch-scatter 2.0.6, rdkit 2020.09 | Last push **2023-11**, 10 commits, **no embedding script ships**. Port is feasible (standard GIN) but is real work for 1.8 M params. Deprioritise. |
| GraphMVP | 2022 | [gh](https://github.com/chao1224/GraphMVP) | MIT; ckpts on a **Google Drive folder**, liveness unverified | ~1.8 M (GIN, 5 layers) | ~300 | 2D graph (3D used at pretrain only) | <0.05 | 🔴 README pins torch 1.9.1, PyG 1.7.2, torch-scatter 2.0.9, torch-sparse 0.6.12, **Python 3.7** | Last push **2022-09**. Ranked **14.64 mean rank vs ECFP's 7.52** in §4. Skip. |
| Mol2vec | 2018 | [gh](https://github.com/samoturk/mol2vec) | **BSD-3**; repo **ARCHIVED 2023-01, read-only** | n/a (word2vec) | **300** | SMILES → Morgan substructure "sentences" | <0.02 | 🟡 CPU-only; risk is **gensim** — no download URL for `model_300dim.pkl` in the repo, and gensim-3.x pickles commonly fail on gensim 4.x | Historical baseline only. Ranked 10.36 vs ECFP 7.52 in §4. |
| MolE | 2024, Nat. Commun. | [gh](https://github.com/recursionpharma/mole_public) | **CC-BY-NC-4.0** | not documented | not documented | 2D graph (DeBERTa disentangled attention) | ? | 🟡 2 commits total | ⚠ README's "Download pretrained models" is **still a TODO** — `mole_predict.encode()` exists but the weights it needs are **not published**. Treat as unavailable. (It is the model that won the Polaris potency track — see §4.) |
| SmilesT5 | 2025 | [gh](https://github.com/hothousetx/smiles_t5) · [hf](https://huggingface.co/hothousetx/smiles_t5) | ⚠ **GPL-3.0** (viral) | not documented | not documented | SMILES (CSV column) | ? | 🟢 likely low (T5 via `transformers`), torch pin unverified | Ships `extract_embeddings.py` → mean-pooled `.pt` files. Small repo, niche. GPL matters if we ship code. |
| **MACE-MP-0 / MPA-0** | 2023-24 | [docs](https://mace-docs.readthedocs.io/en/latest/guide/foundation_models.html) · [gh](https://github.com/ACEsuit/mace) | **MIT** (MP-0/MPA-0 only; OMAT-0, MATPES-0, MH-1/0, MDP are **ASL**) | small/med/large | invariant node features (128-ch/layer medium) | **3D conformer** (ASE `Atoms`) | 0.02–0.06 | 🟢 `mace-torch` 0.3.16: `torch>=1.12` floor, **`e3nn==0.4.4`** exact but pure-python; `cuequivariance` is an opt-in `[cueq]` extra | `mace_mp(model="medium")` → **`get_descriptors(atoms)`**. **89 elements incl. Fe** — the only MACE that can see the heme. Materials-domain training is a real distribution mismatch. |
| MACE-OFF23 / OFF24 | 2023/24 | [gh](https://github.com/ACEsuit/mace-off) | ⚠ **ASL — academic only, explicitly no commercial use**; weights are plain files in the repo | small/med/large | ~256 (concat over 2 layers) | 3D conformer | **0.007 / 0.018 / 0.053** | 🟢 same as above | Organic-chemistry-trained, so better-matched than MP-0 — but **H,C,N,O,F,P,S,Cl,Br,I only: no Fe**, so it cannot describe the heme. Academic licence. |
| SchNet / PaiNN (SchNetPack) | 2017 / 2021; pkg **2.2.0, 2025-12** | [gh](https://github.com/atomistic-machine-learning/schnetpack) | code open; ⚠ **no general-purpose pretrained checkpoints ship** | ~1 M | 128 | 3D conformer | n/a | 🟢 `torch>=2.5.0` floor, `requires_python>=3.12`, **no torch-scatter** | **Not a foundation model.** Include only as a trainable 3D baseline; the "pretrained SchNet embedding" people cite is usually a QM9 checkpoint someone made themselves. |
| EquiformerV2 | 2023 | [gh](https://github.com/atomicarchitects/equiformer_v2) | code **MIT**; **checkpoints only via the FAIR-Chem zoo** | ≤153 M | irreps-based | 3D periodic (OC20/OC22) | ~0.6 | 🔴 fairchem torch pin (below); standalone needs e3nn + PyG + torch-scatter | Trained on catalyst surfaces — **domain mismatch** for drug-like ligands. Skip. |
| DimeNet++ / GemNet-OC / eSEN | 2020-25 | FAIR-Chem zoo | open, various | 2–200 M | — | 3D periodic | 0.1–1 | 🔴 fairchem torch pin | Catalysis domain. Skip. |
| **UMA / OMol25** (Meta FAIR) | 2025 | [hf](https://huggingface.co/facebook/UMA) · [hf](https://huggingface.co/facebook/OMol25) | FAIR Chemistry License v1 (commercial OK) but **GATED — `gated: manual`**, needs name/DOB/org + terms; unavailable in CN/RU/BY | uma-s / uma-m (MoE) | — | 3D, all elements | several | 🔴🔴 **hard blocker.** `fairchem-core` `requires_dist` across its whole history: 1.10.0→`torch~=2.4.0`; 2.0.0b0–2.5.0→`torch~=2.6.0`; 2.8.0–2.21→`torch~=2.8.0`; 2.22.0→`torch~=2.13.0`. **No release is installable on torch 2.5.1** | Scientifically the most interesting 3D entry (OMol25 is genuinely biomolecular, unlike OC20) and therefore the most painful loss. Run it in its own container on Explorer/Modal, or not at all. |
| ⚠ **UME** (Universal Molecular Encoder, Lobster) | 2025-26 | [gh](https://github.com/prescient-design/lobster) | code Apache-2.0, active (2026-09) — but **weights are NOT public** | 12 M / 90 M / 480 M / 870 M | NeoBERT `hidden_size` per checkpoint | **SMILES + amino acids + nucleotides + 3D coords in ONE shared encoder** | ? | 🟡 unverified | **Exactly the model this program wants.** But `_ume_models.py` points at `s3://prescient-lobster/...` with a `TODO: currently, these will work for internal users`, and the `prescient-design` HF org publishes **zero models**. Their published checkpoints (`asalam91/lobster_*`) are protein-only. → **BLOCKED**, worth an email. |
| ether0 | 2025-06 | [hf](https://huggingface.co/futurehouse/ether0) | open, **Apache-2.0**, ungated | **23.6 B** | 5120 (Mistral-Small-24B base) | SMILES-in-text | **~47** bf16, 21 shards | 🟡 fits 80 GB but tight | **A reasoning/generative model, not an embedder.** No embedding API. Listed for completeness; don't build on it. |
| ChemDFM v1.5-8B / v2.0-14B | 2024-25 | [hf](https://huggingface.co/OpenDFM/ChemDFM-v2.0-14B) | open weights | 8 B / 14 B | LLM hidden | chemistry text + SMILES | ~16 / ~28 | 🟢 stock transformers | Chemistry **dialogue LLM**. No embedding API. Completeness only. |
| MolPILE | 2025 | [hf dataset](https://huggingface.co/datasets/scikit-fingerprints/MolPILE) | mixed per source (CC0 / CC-BY / **CC-BY-NC** for Mcule, ChemSpace) | — | — | SMILES | **5.16**, 238 M molecules | none (data) | ⚠ **It is a DATASET, not a model** — no pretrained models released with it. Only relevant if we pretrain. |

### Ligand notes, the ones that will actually bite

- **Uni-Core is the trap, and it is now avoidable.** Uni-Core's last release is **0.0.3 (June 2023)**;
  prebuilt wheels stop at `cu118torch2.0.0-cp{38,39,310}`, and a source build compiles nine fused
  CUDA extensions. **`unimol_tools` dropped the Uni-Core dependency in v0.1.0 (June 2024)** and
  reaches the same MIT weights for both Uni-Mol v1 and v2. Use it; never touch Uni-Core. Its one
  real side effect is `numpy<2.3,>=2.0` — a forced numpy-2 upgrade that breaks anything compiled
  against numpy 1.x. **Own venv.** Also set `HF_ENDPOINT` explicitly; it silently falls back to
  `hf-mirror.com`.
- **MoLFormer has two silent footguns.** (1) `main` targets transformers v5
  (`transformers.masking_utils.create_bidirectional_mask`); on transformers 4.x pass
  `revision="compat-v4"`. (2) **Pass `deterministic_eval=True`** — the linear-attention random
  features otherwise make embeddings non-reproducible *between calls*, which would silently poison
  any consensus or agreement metric, exactly the class of bug FINDING 021 already cost this repo.
  Also: trained on isomeric-info-stripped SMILES ≤202 tokens → **blind to stereochemistry**.
- **"ChemBERTa-77M" is a corpus size, not a parameter count.** hidden 384, **3 layers**, vocab 600.
  Neither the DeepChem nor the seyonec repos carry a licence tag on HF — legally unspecified,
  which may matter at publication.
- **The 2022-era GNNs (GraphMVP, MolCLR) are the weakest entries here** — ~10-commit repos pinned
  to torch 1.7–1.9 and PyG 1.6–1.7, neither shipping an embedding-extraction script, and both
  ranking below ECFP in §4. A day of porting for 1.8 M parameters. No.
- **The FAIR-Chem pin is absolute** and takes out UMA, OMol25, eSEN and the whole OCP zoo in one
  stroke: `~=` is a compatible-release clause and **every fairchem-core release ever published
  excludes torch 2.5.1**.
- **MACE is the one 3D family that installs cleanly** — and the licence/element split matters:
  MACE-**OFF** is organic-trained but **ASL (academic-only) and has no Fe**; MACE-**MP-0** is MIT
  with 89 elements including Fe but is materials-trained. For a heme protein, MP-0 is the only
  option that can see the iron at all.
- **Three "models" that are not what the name suggests:** MolPILE is a dataset;
  `hspark1212/chemeleon` is crystal-structure diffusion; ether0 and ChemDFM are chat LLMs with no
  embedding API.
- **Two models we would most want are unobtainable:** **UME** (the only published
  shared SMILES+protein+3D encoder — private S3) and **MolE** (the Polaris potency winner —
  "download pretrained models: TODO"). See §5.

### Recommended order for the ligand axis

1. **ECFP4 (count) + RDKit descriptors + LightGBM** — the control. Get this number *first*; §4 says
   several rows above will not beat it.
2. **Boltz-2 `--write_embeddings`** — free, MIT, no new deps, and the only pose- and
   protein-conditioned ligand representation available to us.
3. **MoLFormer-XL** — 190 MB, Apache-2.0, pure torch, real embedding API, `deterministic_eval=True`.
4. **Uni-Mol v1 (then v2-84m) via `unimol_tools`** — the only clean 3D-conformer embedding.
5. **CheMeleon** via chemprop — low-risk, and the natural ablation against the RDKit-descriptor
   control it was distilled from.
6. **MACE-MP-0 `get_descriptors`** — only if a physics-flavoured 3D embedding that includes Fe is
   wanted.
7. **Do not attempt:** SMI-TED (2021 CUDA kernels), native Uni-Core paths, UMA/OMol25/EquiformerV2
   via fairchem (torch pin), MolE and UME (weights unpublished), GraphMVP/MolCLR (dead stacks).

## 4. What actually beats a fingerprint

**Short version: almost nothing, reliably.** The strongest, most-replicated finding in the
2019–2026 literature is that **count-based ECFP/Morgan (or RDKit descriptors) + a gradient-boosted
tree** is a baseline that pretrained molecular embeddings have still not decisively beaten on
property prediction — and on *protein–ligand binding affinity* the deeper problem is that
structure-based deep models are largely **not using the structure at all**.

### 4.1 Property prediction — head-to-head

| Study | Year | Task / dataset | Comparison | Result | URL |
|---|---|---|---|---|---|
| Praski, Adamczyk, Czech — *Benchmarking Pretrained Molecular Embedding Models* | 2025 (v2 2026) | 25 datasets (7 MoleculeNet + 18 TDC), **scaffold splits** | 25 pretrained embedders vs ECFP, hierarchical Bayesian test | ECFP mean rank **7.52**, AUROC **79.89%**. Beaten only by CLAMP (5.40/82.55), R-MAT (6.08), MolBERT (6.92). **MoLFormer 9.50/79.80, Mol2Vec 10.36, Uni-Mol2 13.32, GraphMVP 14.64 — all worse than ECFP.** Only CLAMP (itself fingerprint-based) is statistically significant | [arXiv:2508.06199](https://arxiv.org/abs/2508.06199) |
| Sun, Dai, Yu — *Does GNN Pretraining Help Molecular Representation?* | 2022, NeurIPS | MoleculeNet + split/feature ablations | SSL-pretrained GNN vs none vs hand-crafted features | **No statistically significant gain** from self-supervised pretraining; gains shrink with richer features and balanced splits | [arXiv:2207.06010](https://arxiv.org/abs/2207.06010) |
| Deng et al. — *Key elements underlying molecular property prediction* | 2023, Nat. Commun. | MoleculeNet + opioid/activity sets, **62,820 models trained** | Fixed FP/descriptors vs SMILES LM vs graph | "Representation learning models exhibit **limited performance**… in most datasets"; dataset **size** is what lets them win | [10.1038/s41467-023-41948-6](https://doi.org/10.1038/s41467-023-41948-6) |
| van Tilborg et al. — MoleculeACE / activity cliffs | 2022 (corr. 2023, 2024), JCIM | 30 ChEMBL targets, 48.7k mols, cliff-aware splits | 24 methods: ECFP/descriptors + SVM/RF/GBM vs MPNN/GAT/SMILES-transformer | "**ML approaches based on molecular descriptors outperformed more complex deep learning methods**"; all degrade on cliffs. The 2023 correction did not change the conclusion | [PMID 36456532](https://pubmed.ncbi.nlm.nih.gov/36456532/) |
| Adamczyk et al. — *Fingerprints are strong models for peptide function* | 2026, Bioinformatics 42(5) | **132 datasets**, 6 suites | Count ECFP/TT/RDKit + LightGBM vs GNNs **and protein LMs (ESM2, ProtBERT, ProtT5)** | Fingerprints SOTA or tied: Peptides-func AUPRC **74.60** vs 73.11; AMPBenchmark AUROC **97.37** vs 96.79 — with ~22k params and **19 s** vs 60 GPU-hours | [10.1093/bioinformatics/btag179](https://doi.org/10.1093/bioinformatics/btag179) |
| Notwell & Wood (MapLight) | 2023 → TDC leaderboard | 22 TDC ADMET benchmarks | CatBoost on ECFP+Avalon+ErG+200 RDKit props vs published DL | "RF or SVM paired with **ECFP consistently outperformed recently developed methods**"; still top-ranked on TDC | [arXiv:2310.00174](https://arxiv.org/abs/2310.00174) |
| Koleiev et al. — *Critical assessment of TDC ADMET leaderboards* | 2026, bioRxiv | Top-3 models × 22 endpoints | Reproducibility + leakage audit | Only **3 audited methods reproduced** (CaliciBoost, MapLight, MapLight+GNN — all fingerprint/descriptor + boosting). Direct and indirect **leakage** found; deliberate contamination inflated metrics and distorted ranks | [10.64898/2026.02.26.708193](https://doi.org/10.64898/2026.02.26.708193) |
| ⚑ *counter-example:* Burns et al. — **CheMeleon** | 2025 | 58 Polaris + MoleculeACE assays | D-MPNN pretrained to predict **Mordred descriptors** on 1M PubChem mols, vs RF/fastprop/Chemprop | Polaris win rate **75%** vs RF 68%, Chemprop 32%; MoleculeACE **97%** vs RF 50%. Note *what* won: classical descriptors distilled into a GNN — not language-model pretraining | [arXiv:2506.15792](https://arxiv.org/abs/2506.15792) |
| ⚑ *partial:* ASAP/Polaris Antiviral Blind Challenge (65+ teams) | 2025, JCIM 65(24) | Blinded **time-split** Mpro pIC50 + ADME | Classical FP/descriptor ML vs modern DL / large pretrained | **DL significantly beat classical on ADME**; on **potency, classical stayed competitive and SVM was among the top**. Biggest lever = **extra task-specific data**, not architecture | [10.1021/acs.jcim.5c01982](https://doi.org/10.1021/acs.jcim.5c01982) |
| Inductive Bio post-mortem, same competition | 2025 | Polaris ADMET + potency | Fingerprint baseline vs MolMCL, MolE, MolGPS, finetuned LLMs | Fingerprint baseline ~20th–24th on ADMET but only **1.7% worse than the best pretrained model (MolE)** on potency. "Massive non-ADMET pretraining: mixed results" | [inductive.bio](https://www.inductive.bio/blog/lessons-from-the-polaris-admet-competition) |
| Yang et al. (Chemprop) — the pro-learned datapoint | 2019, JCIM | 19 public + 16 proprietary sets | D-MPNN vs fixed descriptors | D-MPNN "consistently **matches or outperforms**" — and even here none reach experimental reproducibility | [10.1021/acs.jcim.9b00237](https://doi.org/10.1021/acs.jcim.9b00237) |

### 4.2 Binding affinity specifically — the structure is not doing the work

This is the part that matters most for us, because the whole program assumes a *joint*
protein–ligand space is worth building.

| Study | Year | Dataset | Comparison | Result | URL |
|---|---|---|---|---|---|
| Volkov et al. — *On the Frustration to Predict Binding Affinities…* | 2022, J. Med. Chem. | PDBbind, modular MP-GNNs over ligand / protein / complex | Explicit noncovalent-interaction descriptors vs ligand-only vs protein-only | "Explicit description of protein–ligand noncovalent interactions does **not** provide any advantage with respect to ligand or protein descriptors." Nearest-neighbour baselines already do well ⇒ "**memorization largely dominates true learning**." Ligand-graph-only ≈ Pearson **0.749 / RMSE 1.567**, beating the interaction-graph model | [10.1021/acs.jmedchem.2c00487](https://doi.org/10.1021/acs.jmedchem.2c00487) |
| Boyles, Deane, Morris — *Learning from the ligand* | 2020, Bioinformatics | PDBbind 2007/2013/2016 core | RF-Score / NNScore2 / Vina ± 185 RDKit **ligand-only** descriptors | Pearson 0.790→**0.836** (2007), 0.746→**0.780** (2013), 0.814→**0.821** (2016). Ligand-only alone hits **R ≈ 0.71–0.74** — i.e. it predicts "the **mean affinity of a ligand for its binding partners**" | [10.1093/bioinformatics/btz665](https://doi.org/10.1093/bioinformatics/btz665) |
| Graber et al. — *Resolving data bias…* (PDBbind **CleanSplit**) | 2025, Nature Mach. Intell. | PDBbind → CASF-2016 | Structure-based DL (Pafnucy, GenScore, GEMS) vs **ligand-only GEMS** | **~49% of CASF test complexes have structurally similar counterparts in PDBbind train.** Ligand-only GEMS reaches CASF-2016 **RMSE 1.424**, comparable to structure-based DL (Pafnucy R 0.746/1.484; GenScore 0.780/1.362) ⇒ predictions are "clearly **not** based on an understanding of protein–ligand interactions." On CleanSplit Pafnucy collapses; GEMS holds at R 0.803/RMSE 1.308 | [10.1038/s42256-025-01124-5](https://doi.org/10.1038/s42256-025-01124-5) |
| Chen et al. — *Hidden bias in DUD-E* | 2019, PLOS ONE | DUD-E SBVS | 3D-CNN vs AutoDock Vina | CNN enrichment traced to **analogue and decoy bias**; PDBbind-trained CNNs **not superior to Vina** | [10.1371/journal.pone.0220113](https://doi.org/10.1371/journal.pone.0220113) |
| Sieg, Flachsenberg, Rarey — *In Need of Bias Control* | 2019, JCIM | 3 SBVS benchmark sets | Re-evaluation of published ML SBVS | "Bias is learned **implicitly and unnoticed** from standard benchmarks" | [10.1021/acs.jcim.8b00712](https://doi.org/10.1021/acs.jcim.8b00712) |
| Li et al. — **LP-PDBbind** (Leak Proof) | 2023 | PDBbind re-split by seq + ligand similarity | Standard vs leak-controlled | Train/test "cross-contaminated with proteins and ligands with high similarity" | [arXiv:2308.09639](https://arxiv.org/abs/2308.09639) |
| Zhang et al. — *Rethinking generalization of DTA prediction* | 2025, **ICLR oral** | 5 DTA methods × 4 datasets | Random vs similarity-aware split | Random test sets "dominated by samples with high similarity to training"; performance **severely degraded** at low similarity | [arXiv:2504.09481](https://arxiv.org/abs/2504.09481) |
| Kopko, Graber et al. — *Generalization Beyond Benchmarks* | 2025 | Unseen-target splits | Learnable scoring functions, target-held-out | Benchmarks "do **not** reflect the true challenge of generalizing to novel targets"; a little test-target data beats architecture | [arXiv:2512.05386](https://arxiv.org/abs/2512.05386) |
| *Robustness of affinity models to docked and predicted structures* | 2026, JCIM | CASF-2016, 5 pipelines × {crystal, GNINA-docked holo/apo, **AF3 co-folded**} | Crystal-input benchmark vs deployment-realistic input | Real **benchmark→deployment gap**: accuracy drops off crystal inputs. Apo receptors gave **minimal advantage over AF3-predicted**; **co-folding matched or beat rigid-receptor docking**; multipose averaging did **not** recover crystal-level performance | [10.1021/acs.jcim.6c00592](https://doi.org/10.1021/acs.jcim.6c00592) |

### 4.3 Caveats that decide most of these comparisons

- **Split type is the whole ballgame.** Random splits flatter learned models; scaffold, time and
  similarity-aware splits close or reverse the gap. For affinity, ~49% of CASF-2016 leaks from
  PDBbind.
- **Dataset size.** Representation learning needs scale; on the 1k–10k-row assays that dominate
  drug discovery, fixed representations win.
- **Leakage masquerading as interaction learning.** Volkov's nearest-neighbour control and
  Graber's ligand-only ablation are the cheapest and most damning tests available, and are still
  rarely reported. Hermes (arXiv 2602.13503, 2026), a DEL-trained transformer claiming
  held-out-target generalisation, reports no fingerprint baseline in its abstract at all.
- **Baseline quality is usually the confound.** Where a pretrained model "wins", the fingerprint
  baseline is often *binary* ECFP with default hyperparameters. **Count**-based ECFP + tuned
  LightGBM/CatBoost + RDKit descriptors is a materially stronger baseline, and it is the one that
  keeps winning when someone bothers to run it.

### 4.4 Bottom line, and what it means for this program

Learned molecular embeddings have not earned their reputation. Across the largest head-to-heads —
25 models × 25 scaffold-split datasets, 62,820 trained models, 132 peptide datasets, 30
activity-cliff targets — pretrained SMILES/graph/3D embeddings are at best indistinguishable from
count-ECFP + a boosted tree, and **MoLFormer, Mol2Vec, GraphMVP and Uni-Mol2 rank below it.** The
two credible exceptions are instructive rather than encouraging: **CheMeleon** wins by distilling
**classical Mordred descriptors** into a GNN, and the Polaris ADME result was driven mainly by
**extra task-specific data**.

For affinity it is worse than "no gain": ligand-only features reproduce nearly all of the reported
structure-based accuracy, nearest-neighbour memorisation explains the leaderboards, and explicit
interaction descriptors add nothing.

**This is directly a kill criterion for the world model.** `README.md` kill criterion 1 already
says frozen embeddings must beat a Morgan-FP + one-hot-target PCM baseline. The literature says
that bar is *higher than it sounds*, and that most published wins over it are leakage. So:

1. Baseline = **count**-ECFP4 (2048 bits) + RDKit descriptors + LightGBM, tuned — not binary ECFP
   with defaults.
2. Always run the **two controls from the affinity literature**: ligand-only, and
   nearest-neighbour-in-training. Report them next to every number, the way this repo already
   reports the pool oracle.
3. Split by **target** (we have 185 from PF00067) and by scaffold, never randomly.
4. Report **within-target** correlation, per the repo's standing rule that between-target variance
   is not selection signal.
5. The one item in this dossier with a *mechanistic* reason to beat a fingerprint is the
   **Boltz-2 pair representation**, because it is the only representation that is a function of
   the protein *and* the ligand *and* the pose. If even that does not beat count-ECFP + LightGBM
   under target-held-out splits, the shared-embedding thesis is falsified cheaply and early —
   which is the desired outcome of running it.

---

## 5. BLOCKED — needs the user

Nothing in the shortlist (§1 items 1–8, §3 items 1–6) is blocked. These are the items where an
external gate, a licence decision, or a dependency quarantine needs Amit.

| # | Item | What is blocked | What we need from you |
|---|---|---|---|
| 1 | **Meta UMA / OMol25** (`facebook/UMA`, `facebook/OMol25`) | HF repos are **`gated: manual`** — name, date of birth, organisation and licence acceptance required before download. Separately, **no `fairchem-core` release in its entire history is installable on torch 2.5.1** (`torch~=2.4` / `~=2.6` / `~=2.8` / `~=2.13`). | (a) Accept the FAIR Chemistry License on HF under your account if you want OMol25 at all; (b) decide whether it is worth an **isolated container on Explorer or Modal** with its own torch. This is the single most scientifically interesting 3D model we cannot run locally. |
| 2 | **UME / Lobster** (`prescient-design/lobster`) | The checkpoint registry points at **`s3://prescient-lobster/ume/checkpoints.json`**, a private bucket, with the code comment *"TODO: currently, these will work for internal users"*. The `prescient-design` HF org publishes **zero models**. | This is *exactly* the shared SMILES + amino-acid + nucleotide + 3D-coordinate encoder our world model wants. Worth an **email to Prescient Design / Genentech** asking for weights or a collaboration. Do not plan around it until they answer. |
| 3 | **MolE** (`recursionpharma/mole_public`) | README's "Download pretrained models" section is **still a TODO**; `mole_predict.encode()` has no weights to load. MolE won the Polaris potency track (§4), so this is not an idle ask. | Ask Recursion, or accept that MolE is unavailable. Licence is **CC-BY-NC-4.0** regardless. |
| 4 | **AlphaFold3 weights** | Weights are **not downloadable**: they require an individual, signed, non-commercial request to Google DeepMind, with no redistribution. | Only worth pursuing if we specifically want AF3-space embeddings. **Recommendation: don't.** Boltz-2 is MIT, already running, and exports embeddings with a flag. |
| 5 | **Non-commercial licences — a decision, not a gate** | Freely downloadable but **NC**: ESM C 600M and ESM-3 (Cambrian NC), Ankh / Ankh3 (CC-BY-NC-SA-4.0), MACE-OFF23/24 (ASL, academic only), MolE (CC-BY-NC), Chai-1 local weights (README restricts to non-commercial despite an Apache-2.0 tag on HF). ChemFM's licence is **self-contradictory** (HF `mit` vs GitHub CC-BY-NC). | Tell us whether this work must stay commercially clean. If yes, the substitutions are: **ESM C 300M** (Cambrian *Open*) instead of 600M, **ESM-2 / AMPLIFY** (MIT) instead of Ankh, **MACE-MP-0** (MIT) instead of MACE-OFF. All are already in the tables. |
| 6 | **`torch-scatter` / `torch-sparse` wheels** | ESM-IF1, MiniMol, MolCLR and GraphMVP all need them. Prebuilt wheels must match torch 2.5.1 + the exact CUDA build, otherwise it is a long source compile. | Confirm whether we may stand up a **second venv** (or a container) for this family. If not, all four are simply out — and none of them is on the shortlist anyway. |
| 7 | **Scratch space on Explorer** | The shortlist alone is ~30 GB of weights (ESM-2 650M 7.8 + SaProt 5.2 + ProtT5 2.4 + AMPLIFY 1.4 + Boltz-2 6.2 + MoLFormer 0.2 + Uni-Mol 0.2 + ESM C 300M 1.3 + spares). Local C:/D: have ~34/~20 GB and the repo rule is that **bulk artefacts never land locally**. | Confirm a `/scratch` path and `HF_HOME` on Explorer so every download goes there. Nothing in this dossier should be fetched to this box. |

### Not blocked, but flag before anyone relies on it

- **`boltz-as-FM`** (the multi-layer pair-pooling recipe in §2) has **17 stars, no licence file, and
  a README that says "the code is currently being refactorized."** The *flag* (`--write_embeddings`)
  is upstream and solid; their *3072-d hybrid pooling* is not yet a load-bearing dependency.
- **Boltz-2 `z` is O(n²)** — ~128 MB per CYP3A4 pose uncompressed. Pool at write time or this fills
  `/scratch` faster than the P450 corpus does.
- **MoLFormer without `deterministic_eval=True`** returns different embeddings on repeated calls.
  This repo has already lost a finding to a silent determinism/indexing bug (FINDING 021); assume
  it will happen again unless someone asserts `allclose` across two forward passes as a unit test.
