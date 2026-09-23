# PRE-REGISTRATION — does the SELECTED depth curve turn over? A second matched-depth stratum

**Written 2026-09-23, committed before a single Boltz job was submitted.** Nothing below —
no ligand, no rung, no threshold, no acceptance clause — may be edited after the first
`submit`. Deviations, if any are forced, are recorded in `FINDING_040` as deviations, not
folded back into this file.

---

## 1. The question, and why it is worth buying

`FINDING_035` doubled a 14-ligand predicted-Type-I stratum from 20 poses to 40 at matched
conditioning: Δ selected **+0.0116**, Wilcoxon **p = 0.715**, **10 of 14 ligands tied**.
`FINDING_039` §3b recomputed that same curve in **closed form**, with zero sampling error,
and found something nobody had remarked on:

> the exact **selected** score peaks at pool **33** (0.5424) and **falls** to 0.5386 by
> pool 40, while the oracle stays monotone — and the entire effect is **one ligand**
> (08J, +0.1721 between pool 30 and 40, against CFF's −0.1195), with **ten of fourteen
> ligands unable to move at all**.

039 named its own successor explicitly:

> What would settle it: the same exact closed-form curve on a **second matched-depth
> stratum of ≥ 30 ligands**.

This is that experiment. It is not for the record. On release day we choose a generation
depth. If depth keeps paying past 20 poses we buy more; if the selected curve genuinely
turns over, buying past ~30 is wasted wall clock and may actively cost score.

---

## 2. The stratum, and the rule that picks it — fixed here, applied after

**Rule (stated before any pose exists):**

> **Take every validation ligand in the prediction-side Type II majority — every ligand of
> the 87 in `data/processed/binding_mode_labels_cyp3a4.csv` with
> `pred_fe_donor_median <= 2.6 Å`. This is, by construction, the exact complement of
> `FINDING_035`'s 14-ligand predicted-Type-I stratum: overlap is zero by definition, not
> by inspection.**

Three reasons this is the right rule and not a convenient one:

1. **It cannot overlap 035.** The two strata partition the 87 on a single prediction-side
   threshold. A second stratum that shared ligands with the first would not be a second
   stratum.
2. **It is where the test set will sit.** `FINDING_033` measured the validation set at
   83% Type II and identified the Type I/II mix as the one unmeasured exposure. 035 bought
   depth on the 17% and has no coverage of the 83% at all.
3. **It is prediction-side**, so it is computable on a blind set at inference time (the
   repo's no-leaky-features rule), and it is the same labelling `FINDING_034`/`035` used.

**Applied: n = 73.** (87 validation ligands → 73 pass `pred_fe_donor_median <= 2.6`, 14
excluded as the 035 Type-I stratum.) That is 2.4× the ≥ 30 that 039 asked for.

Reference depth, checked before committing: **all 73 have ≥ 6 reference poses**
(6:18, 7:43, 8:5, 9:2, 10:1, 11:1, 12:3), every one above the 4 that `FINDING_011`
requires for `xeng`. **No ligand is dropped for reference depth**, and no ligand-level
filter beyond the rule above is applied.

Ligand list: `data/processed/depth_turnover_stratum.csv` (id, smiles, pdb — the same three
columns `validation_ligands.csv` carries, subset by the rule).

**Batching, fixed in advance.** 73 ligands × 20 samples is ~8 GPU-hours, over the 7-hour
walltime the `gpu` fallback partition allows, so generation is split into **four jobs**.
The split is **round-robin over the alphabetically sorted ids** (`i % 4`), tags
`dt2a` (19), `dt2b` (18), `dt2c` (18), `dt2d` (18) — *not* contiguous blocks. This is a
direct response to `FINDING_038`, where an **ordered** csv under `--batch` packing put all
three Type I picks into the single job that failed, manufacturing a chemistry signal out of
one job's death. Round-robin makes any job-level failure chemistry-neutral.
Manifest: `data/processed/depth_turnover_batches.json`.

---

## 3. The conditioning to match, and how the match is verified

Generated on Explorer with the **existing** path — `scripts/explorer/boltz_depth.py`,
stages `plan / conditioning / stage / push / verify / submit / poll / collect / score` —
documented in `docs/RUNBOOK_explorer_boltz.md`. Nothing is rebuilt.

| what | value, and how it is checked |
|---|---|
| `build_yaml` | **not re-implemented** — AST-lifted out of `scripts/cofold/modal_boltz.py`, source sha1 **`aefda51308eadcbd68872e7b8934df01374ac1b0`**, the same sha1 `FINDING_035` recorded |
| MSA | `data/reference/cyp3a4.a3m`, **6,979 sequences**, md5 **`6de0ee1330ea9aa1aef1efb6ec6a2e3d`**, 0 NUL bytes — re-hashed on the far side **inside a Slurm job** |
| checkpoint | `boltz2_conf.ckpt`, md5 **`2f0a1775bf8fc366a1a85e2019eca288`**, hashed **in a job** with explicit `--mem`, never on the login node |
| boltz / torch | 2.2.1 / 2.5.1 — the version Modal pinned |
| heme bond | `[A,442,SG] ↔ [H,1,FE]`, asserted by Boltz's own YAML parser on every input |
| chains / steering | A,H,L / unsteered — no pocket, contact or template block |
| samples / recycling / sampling steps | **20 / 3 / 200** |
| seed | **202** (Modal pool 1, `FINDING_035` 101) |

**A0 — the conditioning-match check, and it is reported BEFORE any curve.**

> A0 passes iff | mean LDDT-PLI of the 1,460 new poses − mean LDDT-PLI of the existing
> 20-pose Modal pool on the same 73 ligands | **≤ 0.020**.

This is `FINDING_035`'s own bar, on its own statistic (035 measured −0.0019;
`FINDING_034`'s unmatched OpenProtein arm measured −0.057). **If A0 fails, the depth curve
is not reported as a matched-depth curve and the finding is NOT REACHED** — an unmatched
arm makes the whole curve meaningless and that is decided here, not afterwards.

Secondary, reported but not gating: fraction coordinating and median Fe–donor distance in
both arms (Type II character should be *preserved*, i.e. mostly coordinating here, which is
the mirror of 035's Type I stratum).

---

## 4. Confounds — every conditioning difference we are forced to accept

Named here, before the result, because a confound named afterwards is an excuse.

| # | difference | why it is forced | why it is judged non-material |
|---|---|---|---|
| **K1** | `--no_kernels` (pure PyTorch) vs the Modal pool's cuequivariance path | this env's `cuequivariance-*-cu12` wheels die on `libnvrtc.so.12` | mathematically the same model; `FINDING_035` measured A0 = −0.0019 with K1 in force, and the smoke ligand `1RD` scored 0.288/0.302 against its Modal pool's mean 0.283 |
| **K2** | GPU is a V100-SXM2 or H200, Modal ran A100/L40S | whatever Slurm gives | inference only; no training-time nondeterminism claim |
| **K3** | torch 2.5.1+**cu121** vs Modal's 2.5.1+**cu124** | the staged venv | same torch minor |
| **K4** | seed **202** vs the pool's 1 | intentional — identical seeds would risk redrawing the same poses | dedupe at 0.05 Å in the heme frame is run regardless |
| **K5** | 18–19 YAMLs per job vs Modal's 1 YAML per job | walltime | `FINDING_035` already carried this (14/job) at A0 = −0.0019 |
| **K6** | **four** jobs, not one — a job-level effect (node, GPU model) could align with a batch | walltime again | round-robin assignment (§2) decorrelates batch from chemistry; per-batch pool means are reported |
| **K7** | `1RD` already has **two** poses from `FINDING_035`'s smoke run | it is the first Type II id alphabetically and the rule takes it | those two poses live in `matched_depth_smoke_poses.csv` and are **not** used here; 1RD gets 20 fresh poses at seed 202 |
| **K8** | the new arm is scored by `boltz_depth.score`, the pool by the shipped board | different entry points to the same scorer | controls **N2** and **N2b** reproduce the shipped board to four decimals and assert `argmin(xeng) == cypstruct.xengine.select()` to 1e-9 |

---

## 5. Controls — each with counts, run before the endpoint

| id | control | bar |
|---|---|---|
| **N2** | the shipped column reproduces on all 87 | selected 0.6164, oracle 0.6975, random 0.5769, gain +0.0395, ρ −0.2582, correct sign 75.86% — to 5e-4 |
| **N2b** | `argmin(xeng)` **is** `cypstruct.xengine.select()` | max abs diff < 1e-9 on all 87 |
| **N1** | residue numbering, `FINDING_021` | per-pose renumber offset and identity printed; **0** poses below 0.95 identity; `mapped` true on all |
| **N4** | distinct poses, not jobs | md5 uniqueness **and** 0.05 Å heme-frame dedupe (`FINDING_009`/`034`); per-ligand distinct counts printed |
| **N5** | every filter, firing, with a count | 87 → 73 by the rule; scoring skips (`no_crystal`, `no_ligand`, `load_failed`, `mapping_failed`, `no_heme_frame`, `low_renumber_identity`, `no_files`) all printed even when zero |
| **C7** | exact-zero LDDT-PLI rows | counted and attributed (ejection vs numbering) rather than dropped |
| **REF** | reference depth per ligand | min ≥ 4 |

**The noise floor.** `FINDING_039` established floors are **population-specific** and that a
4,000-draw p95 carries ±0.0009, so anything within 0.003 of a floor is undecided. This
stratum is a **new population** (prediction-side Type II, n = 73) and its floor is
**recomputed here at ≥ 200,000 draws**, at depth 20 and at depth 40, with the standard
error of the p95 stated. No floor is quoted from another finding.

---

## 6. The curves, and the rungs

Both curves are computed **in closed form**, not by subsampling. `argmin` over a uniformly
random subset has a combinatorial distribution, so the exact expectation is available and
has **zero sampling error** — this is `FINDING_039` §3b's method, applied to a stratum with
5.2× the ligands.

And the reason it matters, from `FINDING_035`'s own methods correction: **a subsample
curve's last rung is the pool's deterministic maximum, with zero variance, while every
earlier rung is an expected maximum over random subsets — so the final slope is inflated
and over-prices new samples.** The union curve below has the same property at its last
rung for the *oracle*; the closed form removes the Monte-Carlo error but not that
structural fact, and it is stated in the finding rather than glossed.

**Primary — the UNION subsample curve.** Draw `d` poses uniformly from each ligand's
40-pose union (20 Modal + 20 Explorer), exactly. Report **oracle**, **selected** and
**random** separately at every rung:

> **d = 5, 10, 15, 20, 25, 28, 30, 33, 35, 38, 40**, plus the full sweep `d = 1…40` to
> locate the argmax.

**Secondary — the AUGMENT curve.** All 20 existing poses kept, `k = 0…20` new poses added
(pool 20 → 40). This is the *purchase* curve — what a depth buy actually does to a pool
that already exists — and it is the one `FINDING_035` reported as its primary endpoint.

**Per-ligand decomposition**, so a one-ligand artefact cannot hide again: the per-ligand
selected score at **d = 20, 30, 33, d\*** (the argmax) and **40**, and the per-ligand delta
peak−40, written to csv in full and tabulated for every ligand whose |delta| > 0.

**And the bound on what any depth purchase can buy:** the count and fraction of ligands
with **at least one new pose beating their incumbent on `xeng`** — the only ligands depth
can move at all. Reported plainly, before the curve.

---

## 7. Acceptance — when a turnover is called real

Let `d*` = argmax over `d ∈ [1,40]` of the exact union **selected** curve, and
`Δ_peak = E[sel(d*)] − E[sel(40)]`.

| clause | requirement |
|---|---|
| **T1 — existence** | `d* < 40` (the peak is strictly interior) |
| **T2 — magnitude** | `Δ_peak ≥ 0.0050` |
| **T3 — not one ligand** | the sign of `Δ_peak` survives **leave-one-ligand-out**: deleting any single one of the 73 ligands leaves `Δ_peak > 0` |
| **T4 — inference** | a ligand-level bootstrap (10,000 resamples) 95% CI on the per-ligand `sel(d*) − sel(40)` **excludes 0**, *and* Wilcoxon signed-rank over those per-ligand deltas gives **p < 0.05** |

| verdict | fires when |
|---|---|
| **TURNOVER REAL** | T1 ∧ T2 ∧ T3 ∧ T4 |
| **TURNOVER IS ONE LIGAND** | T1 ∧ ¬T3 — i.e. `FINDING_039`'s 08J result reproducing at larger n, which is a curiosity, not a law |
| **TURNOVER UNDECIDED** | T1 ∧ T3 ∧ (¬T2 ∨ ¬T4) — an interior peak too small or too noisy to act on |
| **NO TURNOVER (monotone)** | ¬T1 — the selected curve does not turn over on this stratum |
| **NOT REACHED** | A0 fails, or fewer than 30 ligands yield 20 distinct new poses |

**The depth recommendation is derived from the verdict, not chosen:**
under **TURNOVER REAL** → generate to `d*`; under **NO TURNOVER** or
**TURNOVER IS ONE LIGAND** or **TURNOVER UNDECIDED** → generate to the depth at which the
selected curve's marginal gain per doubling falls below `FINDING_004`'s **+0.0125**
benchmark, reported from the same curve.

**Ligand-yield fallback, fixed here.** The primary curve is computed on ligands yielding
**exactly 20 distinct** new poses. Any ligand yielding fewer is reported with its count and
excluded from the primary curve (a ragged depth axis is not a depth axis); a sensitivity
curve over all ligands with ≥ 16 distinct is reported beside it. If the primary set falls
below 30 ligands, the verdict is **NOT REACHED**.

---

## 8. Traps this run must not walk into

- **Explorer's login node SIGKILLs long reads and the killed process can exit 0.** Every
  big-file check runs **inside a Slurm job** with explicit `--mem` (`verify`).
- **`sbatch`, never `srun` over ssh.** Three earlier attempts in this project were killed
  by client-side timeouts mid-run.
- **Poll for ADVANCING progress** — the count of mmCIFs on disk climbing — never for
  liveness.
- **A Slurm job can report `FAILED` with all its work written** (`find | head` takes
  SIGPIPE under `pipefail`, `FINDING_035` job `10529025`). Judge by four independent
  signals: step exit code, boltz's own failure count, the file count, and per-ligand depth.
- **Count distinct poses, never jobs** (`FINDING_009`, `FINDING_034`).
- **Nothing lands on `D:`** — pools stay on Explorer `/scratch` and in `C:\cyp_struct`.
- Nothing that was not created by this run is deleted, on either machine.

---

## 9. Deliverables

| what | where |
|---|---|
| this pre-registration | `docs/PREREG_depth_turnover.md`, committed before `submit` |
| the finding | `docs/FINDING_040_depth_turnover.md` |
| ligand list and batches | `data/processed/depth_turnover_stratum.csv`, `depth_turnover_batches.json` |
| scored new poses | `data/processed/depth_turnover_poses.csv` |
| controls incl. A0, N1, N4, N5, C7, the floor | `data/processed/depth_turnover_controls.json` |
| exact curves, both | `data/processed/depth_turnover_curves.json` |
| per-ligand decomposition | `data/processed/depth_turnover_per_ligand.csv` |
| release-day depth recommendation | `docs/DROP_DAY_PLAYBOOK.md`, generation step |
