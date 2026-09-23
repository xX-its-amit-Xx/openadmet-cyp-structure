# PRE-REGISTRATION — depth on the predicted-Type-I stratum at MATCHED conditioning

**Written:** 2026-09-22, before a single stratum pose was generated.
**Committed before submission.** Not to be edited afterwards; the result goes in
`docs/FINDING_035_matched_depth.md`.

---

## Why this exists

`FINDING_034` measured the depth law and could not buy it:

> **The depth actually purchasable is not depth.** Losing the MSA — forced, because
> `boltz2` with an uploaded MSA fails server-side on OpenProtein — produced poses whose
> ceiling matches and whose average is 0.057 lower. Mixed into the pool they add **+0.0118
> to +0.0171 of oracle and subtract 0.0065 to 0.0082 of score.**

and closed by naming the experiment that would settle it:

> It does not test depth at matched conditioning, which is the experiment that would
> actually settle it and which needs Modal or a staged-MSA Explorer run.

Modal is still over cap. The staged-MSA Explorer run is now possible, because the
fine-tuning campaign left behind a working offline Boltz-2 install on `/scratch`. So this
is `FINDING_034`'s own stated follow-up, run at the venue it named.

> **The question.** When the new poses are generated at the SAME conditioning as the pool
> they join, does depth move SELECTED LDDT-PLI — or does it join FINDINGs 013, 016, 027,
> 031 and 034 as a sixth consecutive oracle-up / selection-flat result?

---

## The ligands — fixed, and already fixed once

The **predicted-Type-I stratum**: median `fe_donor_dist` over a ligand's own 20 existing
pool poses **> 2.6 Å**. Prediction-side, so nothing from any crystal enters the selection.
`COORD_MAX = 2.6` is the repo constant from `build_reference_set.py`, unchanged.

**n = 14**, reproduced today from `binding_mode_labels_cyp3a4.csv` and identical to
`FINDING_034`'s list:

```
08J  08Y  A1CIW  CFF  D0R  ERY  MF8  MWS  MWV  MWY  PG0  PG4  TCI  YNV
```

Filter counts, which will be re-printed with the result: 87 validation ligands → **14**
pass → **73** excluded as predicted Type II.

The smoke test that exercises the path uses **`1RD`** — the first predicted-Type-II ligand
alphabetically, chosen by that rule *a priori* and **deliberately not in the stratum**, so
no stratum pose exists before this document is committed.

---

## The conditioning I will match, and how I verified it matches

The existing pool is the **`unsteered`** arm of `poses_scored_val87b.csv`: 87 ligands ×
20 poses, Boltz-2 on Modal. Verified today: 1,740 rows, 87 ligands, exactly 20 samples
each.

**The YAML is not re-implemented. It is imported.** `build_yaml` is lifted out of
`scripts/cofold/modal_boltz.py` by AST extraction and executed, so the new arm's input is
written by the same function that wrote the pool's. Source sha1 **`aefda513…`**, recorded
here so a later edit to that function is detectable.

| | existing pool (Modal) | new arm (Explorer) | matches? | how checked |
|---|---|---|---|---|
| engine / version | Boltz-2, `boltz[cuda]==2.2.1` | Boltz-2, **2.2.1** | **yes** | `importlib.metadata.version` in job `10528945` |
| checkpoint | `boltz2_conf.ckpt`, auto-downloaded by boltz | `boltz_cache/boltz2_conf.ckpt` | **yes** | `get.sh` pulls `huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt` — the same artefact boltz fetches. md5 `2f0a1775bf8fc366a1a85e2019eca288`, hashed **inside a Slurm job** |
| MSA | `/cache/msa/cyp3a4.a3m`, 6,979 seqs | the **same file**, staged | **yes** | md5 `6de0ee1330ea9aa1aef1efb6ec6a2e3d` local **and** on `/scratch`; 6,979 `>` lines; 0 NUL bytes both sides |
| heme bond | `bond: [A,442,SG] ↔ [H,1,FE]` | identical | **yes** | asserted on all 15 YAMLs by Boltz's own parser in job `10528945` |
| chain ids | A protein / H HEM / L ligand | identical | **yes** | same assert |
| steering | none (unsteered arm) | `steer=False` | **yes** | no `pocket`, no `contact`, no `templates` block emitted |
| `diffusion_samples` | 20 | **20** | **yes** | |
| `recycling_steps` | 3 | **3** | **yes** | |
| `sampling_steps` | 200 | **200** | **yes** | |
| ligand SMILES | `cypstruct.chem.standardize` | same call | **yes** | |

**Depth is bought with `diffusion_samples`, not with replicate jobs, and that is
deliberate.** `FINDING_009`/`034` established that `diffusion_samples` does not sample the
ligand **on OpenProtein** — that is a property of their wrapper. Native `boltz predict`
does: `FINDING_034`'s own control found **20 unique `xeng` values on every one of 87
ligands** in this very pool, median within-ligand LDDT-PLI sd 0.062, zero ligands with
sd 0. The pool being extended is itself proof that the mechanism works here. It will still
be **counted, not assumed** (control N4 below).

Per `FINDING_034` item 3b, the whole purchase is **one submission wave**: a single sbatch
job, seed 101.

---

## The confounds — named before any pose is scored

| id | confound | why it is unavoidable | expected direction |
|---|---|---|---|
| **K1** | **`--no_kernels`.** Modal ran the default cuequivariance kernel path; this env's `cuequivariance-*-cu12` wheels were built against a torch that was later rolled back and the backend dies on `libnvrtc.so.12: cannot open shared object file` (reproduced in job `10528945`). The pure-PyTorch fallback is the same model and the same weights, evaluated through different kernels. | rebuilding the env risks the whole path for a numerically equivalent code path | small, unsigned. **This is the largest named risk to the "matched" claim.** |
| **K2** | **GPU.** Modal: A100-40GB / L40S / A100-80GB. Explorer: H200 or V100, whichever the partition frees. | no control over either | small, unsigned |
| **K3** | **torch build** cu121 vs cu124 | env is prebuilt | negligible |
| **K4** | **seed differs** (101 vs the pool's). | **intentional** — an identical seed would regenerate the pool, which is the opposite of buying depth | none; this is the mechanism |
| **K5** | **batching.** Modal ran one YAML per job; this runs a directory of 14. The RNG stream position per ligand therefore differs. | batching is what makes it one wave | none for validity — independent samples are wanted, not reproducible ones |

K1–K3 are all *numerical-path* differences at identical weights, MSA, constraints and
sampler settings. **If the new arm's pool mean falls materially below the existing pool's
(as OpenProtein's single-sequence arm did, by 0.057), the matched-conditioning claim has
failed and the result must be read as `FINDING_034` again, not as a matched test.** That
check is A0 below and is reported **before** the primary endpoint.

---

## The depth target

**One doubling: 20 → 40.** 14 ligands × 20 new samples = **280 new poses**, one sbatch
job, ~1 GPU-hour. This is the interval `FINDING_034`'s §2 curve actually extends
(+0.0127 selected per doubling whole-curve, **+0.0333** in the 10→20 octave) and the one
its §5 projection priced at +0.006 to +0.017 on a 50/50 test set.

---

## What will be measured, oracle before selection, in this order

Everything mirrors `FINDING_034` so the numbers are comparable.

**A0 — the new arm's quality, alone, before it is mixed into anything.** Pool mean
(= random) and oracle at matched depth, new arm vs existing pool, on the same 14 ligands.
This is the confound check, not the endpoint.

**A1 — controls, with counts, including the zeros.**
- **N2** the shipped column reproduces: selected **0.6164**, oracle **0.6975**, random
  **0.5769**, gain **+0.0395**, within-ligand ρ **−0.2582**, correct sign **75.86%**.
- **N2b** `argmin(xeng)` equals `cypstruct.xengine.select()` to 1e-9 on all 87.
- **N1** numbering (`FINDING_021`): renumbering offset and residue-name identity on every
  new pose; any pose under the 0.95 identity bar is counted and reported.
- **N4** **distinct** poses, not jobs: md5 over every collected mmCIF **and** distinct
  `xeng` values per ligand. Every depth number is quoted on the deduplicated pool.
- **N5** filters, every one with its count, including skips that are zero.
- **C7** any exact-zero LDDT-PLI rows, adjudicated as ejection vs numbering the way
  `FINDING_033` did — BiSyRMSD and Fe distance, not assumption.

**A2 — the noise floor INSIDE this stratum.** A random feature, **4,000 draws at n = 14**.
`FINDING_034` measured **+0.0431 / +0.0631** (p95/p99) here against the pooled +0.0137.
It will be recomputed, not quoted.

**A3 — the primary endpoint.** Augmented pool, all 20 existing poses kept plus `k`
distinct new poses, 256 draws, matched depth. Reported at k = 0 (baseline), and at the
matched `k` set by the minimum distinct count over the 14 ligands, plus "all available".

| statistic | form |
|---|---|
| Δ oracle, Δ selected, Δ random | 20 → 20+k, paired by ligand |
| 95% CI | 10,000-resample percentile bootstrap over the 14 ligands |
| **Wilcoxon signed-rank** | paired, on per-ligand selected score |
| **tie count** | ligands where the selector keeps an old pose and nothing changes |
| within-ligand ρ, correct-sign fraction | on the augmented pool |
| poses taken from the new arm | count out of 14 |

**A4 — the rate.** Selected LDDT-PLI per doubling across the augmented curve, against
`FINDING_034`'s Type I **+0.0127** (whole-curve) and **+0.0333** (top octave). This is the
prediction being tested: matched depth should convert at ~29%.

---

## Acceptance rule — fixed now

Let **Δ_sel** = selected LDDT-PLI of the 40-pose pool minus the 20-pose pool, on the 14
stratum ligands, paired.

| verdict | condition |
|---|---|
| **SHIPS** | Δ_sel ≥ **+0.020** *and* its 95% CI excludes 0 *and* A0 shows the new arm's pool mean within **0.02** of the existing pool's (i.e. the conditioning really matched) |
| **REFUTED** | Δ_sel ≤ 0 *and* the CI upper bound < +0.020 — a result that excludes a useful effect, not merely one that fails to find one |
| **MEASURED** | anything between: notably Δ_oracle ≥ +0.020 with Δ_sel short of the SHIPS bar, which is the oracle-up / selection-flat shape |
| **NOT REACHED** | generation did not complete, or A0 fails by more than 0.02 so the arm is not matched and the comparison is `FINDING_034` over again |

**+0.020 is not chosen today.** It is `FINDING_007`'s standing bar, restated in
`docs/README.md` as "anything under +0.020 is noise", and it is the same bar
`FINDING_034` pre-registered.

A secondary, weaker clause, also fixed now: **the rate clause fires** if selected
LDDT-PLI per doubling over 20→40 is ≥ **+0.0125** (`FINDING_004`'s benchmark, the same
clause `FINDING_034`'s SHIPS (3) used). It can fire while the primary fails, and it will
be reported either way.

---

## What is NOT being done

- No threshold, rung, ligand or acceptance bar will move after a score is seen.
- The selector is **not** being changed — `FINDING_033` item 2.
- No crystal information enters stratum selection, pose generation, or the `xeng` feature.
- The smoke test is on `1RD`, outside the stratum, and its poses are **not** scored against
  any crystal before this file is committed.
- Modal is not touched. OpenProtein is not touched.

---

## Artefacts this will produce

| what | where |
|---|---|
| runner (plan/stage/verify/submit/poll/collect/score) | `scripts/explorer/boltz_depth.py` |
| sbatch templates | `scripts/explorer/boltz_depth.sbatch`, `verify_inputs.sbatch` |
| analysis | `scripts/structure/matched_depth_analysis.py` |
| the new poses, scored | `data/processed/matched_depth_poses.csv` |
| controls | `data/processed/matched_depth_controls.json` |
| primary endpoint | `data/processed/matched_depth_augment.json` |
| job ledger | `data/processed/matched_depth_jobs.json` |
| procedure | `docs/RUNBOOK_explorer_boltz.md` |
| result | `docs/FINDING_035_matched_depth.md` |
