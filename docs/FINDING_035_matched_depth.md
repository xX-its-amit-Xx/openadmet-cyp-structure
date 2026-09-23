# FINDING 035 — depth at MATCHED conditioning: the purchase is real, the conversion is not

**Date:** 2026-09-23 · **Status:** measured. **3 Slurm jobs, 282 new complexes, 0 failed
examples, 1 h 34 m of one V100, $0** · **Verdict: MEASURED** — and for the first time in
six attempts, **selection did not go backwards**.

Pre-registered in `docs/PREREG_matched_depth.md`, committed at `4a0a99f` before the
stratum job was submitted. No threshold, rung, ligand or acceptance clause was moved
afterwards.

Generation path and its traps: `docs/RUNBOOK_explorer_boltz.md`.

---

## What was asked

`FINDING_034` measured the depth law on the predicted-Type-I stratum and could not buy it,
because the only venue open sold MSA-less poses whose average was 0.057 LDDT-PLI worse. It
closed by naming its own successor:

> It does not test depth at matched conditioning, which is the experiment that would
> actually settle it and which needs Modal or a staged-MSA Explorer run.

> **The question.** When the new poses come from the SAME generator at the SAME
> conditioning as the pool they join, does depth move SELECTED LDDT-PLI — or does it join
> FINDINGs 013, 016, 027, 031 and 034 as a sixth oracle-up / selection-flat result?

---

## The headline, in three numbers

| | `FINDING_034` (OpenProtein, single-sequence) | **this (Explorer, matched)** |
|---|---|---|
| new arm's pool mean vs the pool it joins | **−0.057** | **−0.0019** |
| Δ oracle, 20 → 40 | +0.0118 | **+0.0097** [+0.0012, +0.0195] |
| Δ **selected** | **−0.0065** | **+0.0116** [−0.0687, +0.1034] |

**The conditioning matched** — that is the thing `FINDING_034` could not get, and it cost
0.0019 instead of 0.057. **The sign flipped.** And **nothing is resolved**: ten of fourteen
ligands are unchanged, Wilcoxon p = 0.715, and two ligands moving in opposite directions
carry the entire mean.

---

## Goal A first — the generation path works, and that is the durable part

Modal is over cap; OpenProtein cannot serve `boltz2` with an uploaded MSA. Explorer now
can, end to end, and almost none of it needed building: the fine-tuning campaign had
already left a complete **offline** Boltz-2 install on `/scratch`.

| job | what | node | result |
|---|---|---|---|
| `10528945` | verify inputs **inside a job** | c0587 (CPU, 16 G) | COMPLETED 0:20 |
| `10528965` | smoke, `1RD` (non-stratum), 2 samples | d1010 V100 | COMPLETED **4:54** |
| `10529025` | **stratum, 14 ligands × 20 samples, seed 101** | d1010 V100 | **1:34:19**, 280/280 mmCIFs, "Number of failed examples: 0" |

**Job `10529025` reports `FAILED` and its work is complete.** The `srun` step exited
`0:0`; the batch step exited 13 on the trailing `find … | head -3` under `set -o
pipefail` — `head` closes the pipe, `find` takes SIGPIPE, *after* every pose is written.
The repo's standing trap is a **passing check over broken work** (`docs/README.md` §2);
this is its mirror, and it would have thrown away a 94-minute job. Verified by four
independent signals before being believed: step exit code, boltz's own failure count,
280 files, and 20 per ligand on all 14. Fixed in the template with `|| true`.

Full procedure, sbatch templates, staging steps, the proxy, and all eight failures hit
(rsync unusable from this box; `$?` after a pipe reporting success for a no-op transfer;
bare `bash` resolving to the WSL launcher and complaining on **stdout**; CRLF shebangs;
`--no_kernels`; no internet on GPU nodes; collect scoping) are in
`docs/RUNBOOK_explorer_boltz.md`.

---

## The conditioning match, verified rather than asserted

`build_yaml` is **not re-implemented** — it is lifted out of `scripts/cofold/modal_boltz.py`
by AST extraction and executed, so the new arm's input is written by the same function
that wrote the pool it joins (source sha1 `aefda513…`).

| | Modal pool | Explorer arm | matched | evidence |
|---|---|---|---|---|
| boltz | `boltz[cuda]==2.2.1` | **2.2.1** | ✔ | job `10528945` |
| checkpoint | auto-downloaded | `boltz2_conf.ckpt`, same HF path | ✔ | md5 `2f0a1775bf8fc366a1a85e2019eca288`, hashed **in a job** |
| **MSA** | `cyp3a4.a3m`, 6,979 seqs | **the same file** | ✔ | md5 `6de0ee1330ea9aa1aef1efb6ec6a2e3d` on **both** sides, 0 NUL bytes |
| heme bond | `[A,442,SG] ↔ [H,1,FE]` | identical | ✔ | asserted by Boltz's own parser, 14/14 |
| chains / steering | A,H,L / unsteered | identical | ✔ | no pocket, contact or template block |
| samples / recycling / sampling steps | 20 / 3 / 200 | 20 / 3 / 200 | ✔ | |
| **kernels** | cuequivariance | **`--no_kernels`** | ✘ **K1** | cueq dies on `libnvrtc.so.12` |
| GPU | A100/L40S | V100-SXM2-32GB | ✘ K2 | |
| torch | 2.5.1+cu124 | 2.5.1+cu121 | ✘ K3 | |
| seed | 1 | 101 | intentional K4 | |
| batching | 1 YAML/job | 14 YAML/job | K5 | |

**K1–K3 were pre-registered as the largest risk to the "matched" claim, with a numeric
test attached (A0): if the new arm's mean fell more than 0.02, the comparison would be
`FINDING_034` again.**

| A0 | |
|---|---|
| new arm pool mean | **0.4577** |
| existing pool mean, same 14 ligands | **0.4597** |
| **Δ** | **−0.0019** — inside the 0.02 bar by an order of magnitude |
| for contrast, `FINDING_034`'s arm | −0.057 |
| fraction coordinating / median Fe–donor | 1.1% / 5.26 Å vs existing 1.8% / 5.30 Å — Type I character preserved |

**The kernel path is not a material confound.** It was also visible on the smoke ligand
before any stratum pose existed: `1RD` scored 0.2876 / 0.3022 against its own Modal pool's
mean 0.2826 (IQR 0.249–0.304), Fe–donor 2.17/2.14 against 2.19.

---

## Controls, with counts, including the zeros

| id | control | result |
|---|---|---|
| **N2** | shipped column reproduces | selected **0.6164**, oracle **0.6975**, random **0.5769**, gain **+0.0395**, ρ **−0.2582**, correct sign **75.86%** — `FINDING_011`/`033` to four decimals |
| **N2b** | `argmin(xeng)` **is** `cypstruct.xengine.select()` | max abs diff **0.0**, equal to 1e-9 on all 87. Asserted, not assumed |
| **N1** | numbering (`FINDING_021`) on 280 new poses | offset **0 on 280 of 280**, min residue-name identity **0.9979**, **0** poses below the 0.95 bar, `mapped` True on **280 of 280** |
| **C7** | 2 exact-zero LDDT-PLI rows | PG0 ×1, PG4 ×1 — BiSyRMSD **23.2 / 25.9 Å**, ligand **20.9–22.4 Å from the iron**, offset 0, mapped True. **Ejections, not numbering**, and the same two PEG fragments that eject in the existing pool |
| **N4** | distinct poses, not jobs | 280 rows → **280 distinct md5** → **280 distinct at 0.05 Å in the heme frame**. **20.0 per ligand** |
| **N5** | filters | 87 validation ligands → **14** pass `median fe_donor_dist > 2.6` → **73** excluded. Scoring skips: no_crystal 0, no_ligand 0, load_failed 0, **mapping_failed 0**, no_heme_frame 0, low_renumber_identity 0, no_files 0 |
| | existing pool diversity in the stratum | **20 unique `xeng` on every one of 14**, median within-ligand sd 0.0804, zero ligands with sd 0 |
| | reference depth | 6–11 per ligand, min **6**, all ≥ the 4 `FINDING_011` requires |

### A venue fact worth as much as the endpoint

**One job bought 20.0 distinct poses per ligand. `FINDING_034`'s 30 OpenProtein replicate
jobs bought 10.79**, and its replicates 16–29 bought **zero**.

`FINDING_009`'s "`diffusion_samples` does not sample the ligand" is a property of
**OpenProtein's wrapper**, not of Boltz. Native `boltz predict --diffusion_samples 20`
diversifies the ligand fully — visible even on the 2-sample smoke run, whose two poses
differ by **4.07 Å** in the heme frame. Depth on Explorer is bought with
`diffusion_samples`, in one job, at 1/30th the job count.

---

## The primary endpoint

All 20 existing poses kept, `k` distinct new poses added, 256 draws. `k_matched = 20`
(every ligand yielded 20 distinct), so this is a **true, full doubling: 20 → 40.**

| pool | **oracle** | selected | random |
|---|---|---|---|
| **20** (existing only) | **0.6499** | **0.5270** | 0.4597 |
| 22 | 0.6509 | 0.5316 | 0.4595 |
| 24 | 0.6524 | 0.5354 | 0.4594 |
| 28 | 0.6542 | 0.5393 | 0.4590 |
| 30 | 0.6553 | **0.5411** | 0.4589 |
| **40** (+20 new) | **0.6595** | **0.5386** | 0.4587 |

| | value | 95% CI | |
|---|---|---|---|
| **Δ selected (20 → 40)** | **+0.0116** | **[−0.0687, +0.1034]** | **2 better / 2 worse / 10 unchanged** |
| **Δ oracle (20 → 40)** | **+0.0097** | **[+0.0012, +0.0195]** | excludes zero |
| Δ random | −0.0010 | | |
| **Wilcoxon signed-rank** | | | **p = 0.715** |
| ligands where the selector ever took a new pose | **4 of 14** | | |
| within-ligand ρ | **−0.2883** (was −0.2357) | | correct sign **78.57%** (was 71.43%) |
| augmented-pool gain over random | **+0.0799** | | vs its own null p95 **+0.0456** |
| noise floor **inside n=14**, 4,000 draws | **p95 +0.0453 · p99 +0.0648** | | (`FINDING_034`: +0.0431 / +0.0631) |

### The pre-registered verdict

| clause | required | measured | fires? |
|---|---|---|---|
| SHIPS (1) | Δ_sel ≥ +0.020 | **+0.0116** | no |
| SHIPS (2) | 95% CI excludes 0 | [−0.0687, +0.1034] | no |
| SHIPS (3) | A0 within 0.02 | **−0.0019** | **yes** |
| **REFUTED** | Δ_sel ≤ 0 **and** CI upper < +0.020 | **+0.0116**, upper +0.1034 | **no** |
| NOT REACHED | generation incomplete or A0 fails | 280/280, A0 passed | no |
| **MEASURED** | anything between | | **YES** |
| *rate clause* | selected ≥ +0.0125 per doubling | **+0.0116** | **no**, narrowly |

**MEASURED.**

---

## Why it did not convert — and it is not the reason FINDING 034 gave

### 1. Two ligands are the entire result, and both turn on a razor-thin margin

| ligand | selected before | after | Δ | the xeng margin that decided it |
|---|---|---|---|---|
| **CFF** | 0.4721 | **0.9547** | **+0.4826** | new winner `xeng` **0.4151** vs old **0.4720** |
| **08J** | 0.4806 | **0.1364** | **−0.3442** | new winner `xeng` **1.5159** vs old **1.5604** |
| YNV | 0.6093 | 0.6448 | +0.0354 | |
| PG0 | 0.4824 | 0.4706 | −0.0118 | |
| *the other ten* | | | **exactly 0.0000** | |

| leverage | Δ selected |
|---|---|
| all 14 | **+0.0116** |
| drop 08J | +0.0389 |
| drop CFF | **−0.0247** |
| **drop both** | **+0.0020** |

A **0.045 Å** Chamfer margin cost 0.34 LDDT-PLI on 08J; a **0.057 Å** margin gained 0.48 on
CFF. That is the mechanism: **a deeper pool does not give the selector a better pose to
find, it gives it more near-ties — and a near-tie is a coin flip with an enormous
payoff.** 08J's new winner is *not* a coordinating pose (Fe–donor 3.14 Å,
`is_coordinated` False) and its BiSyRMSD is 5.10 Å, while the same arm's best pose for
08J reaches 0.4924. The selector had a good new pose available and did not take it.

### 2. The oracle law over-predicts the real purchase — and the reason is methodological

`FINDING_034` measured, by subsampling the 20-pose pool, a Type I oracle rate of
**+0.0428 per doubling** and **+0.0501 in the 10→20 octave**. The real 20→40 purchase
delivered **+0.0097**.

That gap is largely an artefact of where a subsample curve ends, and it is worth stating
plainly because this repo will subsample again. **The last rung of a subsample curve is
the pool's own maximum, with zero variance**, while every earlier rung is an *expected*
maximum over random subsets. So the final octave's slope is inflated, and extrapolating it
past the pool's size overstates what new samples will buy.

Subsampling the **40-pose union** instead — where no rung is anybody's deterministic
maximum — gives a clean curve:

| depth | oracle | selected | random |
|---|---|---|---|
| 1 | 0.4597 | 0.4597 | 0.4587 |
| 2 | 0.5034 | 0.4800 | 0.4587 |
| 5 | 0.5485 | 0.4906 | 0.4587 |
| 10 | 0.5879 | 0.5064 | 0.4587 |
| 20 | 0.6264 | 0.5291 | 0.4587 |
| 30 | 0.6478 | 0.5393 | 0.4587 |
| **40** | **0.6595** | **0.5386** | 0.4587 |

| octave | oracle | **selected** |
|---|---|---|
| 10 → 20 | +0.0385 | **+0.0227** |
| **20 → 40** | **+0.0331** | **+0.0095** |

**The oracle law survives — +0.033 per doubling at this depth, barely below the previous
octave's +0.0385. The SELECTION rate does not: it falls by more than half, +0.0227 →
+0.0095, over the same interval where the ceiling keeps rising almost linearly.** Depth
keeps putting better poses in the pool; the selector keeps not finding them, and finds
proportionally fewer of them the deeper the pool gets.

### 3. The two arms are equal in quality and disagree about which pose is best

| 20-pose arm, alone | oracle | mean | selected |
|---|---|---|---|
| existing (Modal) | **0.6499** | 0.4597 | 0.5270 |
| **new (Explorer)** | 0.6206 | **0.4577** | **0.5501** |

The new arm's *mean* matches to 0.002 and its *maximum* is 0.029 lower — ordinary
extreme-value noise at n=20, and it beats the existing arm's oracle on **5 of 14**
ligands. Notably the new arm **selects better alone** (0.5501 vs 0.5270) and the union
lands between the two, which is the same non-additivity `FINDING_013` found for a second
engine.

---

## Verdict — **MEASURED**

1. **Matched conditioning is achievable and was achieved.** The new arm's pool mean is
   **0.0019** from the pool it joins, against `FINDING_034`'s 0.057. `--no_kernels`,
   a different GPU and a different CUDA build are **not** material. This retires the
   confound that made `FINDING_034` unable to answer its own question.
2. **Depth at matched conditioning raises the ceiling — +0.0097, CI [+0.0012, +0.0195],
   excluding zero — and the oracle law holds at +0.033 per doubling** on an unbiased
   curve.
3. **It does not convert. Δ selected +0.0116, CI [−0.0687, +0.1034], Wilcoxon p = 0.715,
   ten of fourteen ligands unchanged, and +0.0020 once the two swing ligands are
   dropped.** The rate clause misses its +0.0125 bar at +0.0116.
4. **But the sign flipped, and that is new.** Five consecutive pool expansions —
   013, 016, 027, 031, 034 — moved selection *negative*. This one did not. The
   oracle-up / selection-flat pattern is now **oracle-up / selection-flat-but-not-worse**
   once the added poses are genuinely the same kind. The distinction `FINDING_034` drew
   between DEPTH and EXPANSION is supported: expansion costs score, matched depth does
   not.
5. **The conversion rate itself decays with depth.** +0.0227 per doubling at 10→20,
   **+0.0095** at 20→40, while the oracle rate barely moves. Extrapolating a selected-score
   rate measured at shallow depth to a deeper purchase will over-promise, and
   `FINDING_034`'s §5 projection (+0.006 to +0.017 at a 50/50 test set) should now be read
   at the **bottom** of its range.
6. **A subsample curve's last rung is a deterministic maximum.** Do not read its final
   slope as the price of new samples. This is a methods correction to `FINDING_034` §2 and
   to `FINDING_004`.

**What this does NOT license.** n = 14, and the CI is 0.17 wide — this cannot exclude a
real +0.05 effect any more than it can establish one. It tests one doubling, 20→40, on one
stratum; 40→80 is unmeasured. It does not say depth is worthless — the ceiling
demonstrably rises. And the two swing ligands mean the point estimate should be quoted
with `drop-both = +0.0020` attached, always.

---

## What to do at submission time

`FINDING_034` item 3 said: *buy Type I depth at matched conditioning or not at all.* That
is now testable and the answer is narrower than hoped.

1. **The Explorer path is the generation venue. Use it** — it is the only one that
   produces poses at the pool's conditioning, and it is tested
   (`docs/RUNBOOK_explorer_boltz.md`). Budget ~6.6 min per ligand at 20 samples on a
   V100, much less on an H200.
2. **Buy depth with `--diffusion_samples`, in one job, not with replicates.** 20.0
   distinct poses per ligand from one job, against 10.79 from thirty OpenProtein jobs.
3. **Do not expect depth to raise the submitted score.** At matched conditioning a
   doubling is worth **+0.0116 ± 0.09** on selection, and **+0.0095 per doubling** on the
   unbiased curve. Spend it for the ceiling and for robustness, not for the leaderboard.
4. **Do not change the selector on this** — `FINDING_033` item 2 stands. The one
   encouraging signal is that the deeper pool *ranks* better (ρ −0.2357 → −0.2883, correct
   sign 71.4% → 78.6%, gain over random +0.0673 → +0.0799 against a floor of +0.0456).
5. **The near-tie is the target, not the pose.** Both swing ligands turned on a <0.06 Å
   Chamfer margin producing a >0.3 LDDT-PLI swing. A tie-break for the top-2 `xeng` poses
   is a better-evidenced next experiment than more depth — and it is measurable on this
   very pool, with no new inference.

---

## Artefacts

| what | where |
|---|---|
| pre-registration | `docs/PREREG_matched_depth.md` (`4a0a99f`) |
| generation path + traps | `docs/RUNBOOK_explorer_boltz.md` |
| runner | `scripts/explorer/boltz_depth.py` |
| sbatch templates | `scripts/explorer/boltz_depth.sbatch`, `verify_inputs.sbatch` |
| analysis | `scripts/structure/matched_depth_analysis.py` |
| the 280 new poses, scored | `data/processed/matched_depth_poses.csv` |
| controls incl. A0, N1, C7, N4 | `data/processed/matched_depth_controls.json` |
| distinct-pose dedupe | `data/processed/matched_depth_distinct.json` |
| primary endpoint + curve + nulls | `data/processed/matched_depth_augment.json` |
| per-ligand deltas | `data/processed/matched_depth_per_ligand.csv` |
| smoke poses (`1RD`, non-stratum) | `data/processed/matched_depth_smoke_poses.csv` |
| job ledger | `data/processed/matched_depth_jobs.json` |
| the 282 mmCIFs | `C:\cyp_struct\matched_depth\` (local) and `/scratch/shenoy.am/cyp-depth/out/` (Explorer) — never `D:` |

---

# Correction, 2026-09-23 — the two floors here are two populations, and the pool-30 peak is real but is one ligand

Appended by `FINDING_039`. The body is unchanged.

## 1. The floors

This finding prints two, and they are correctly two different quantities:

| printed | population | authoritative (2 × 10⁶ draws, SE ±0.00004) | this finding's error |
|---|---|---|---|
| `noise_floor_in_stratum` **+0.0453 / +0.0648** | prediction-side Type I, **depth 20** | **p95 +0.04490 · p99 +0.06398** | +0.44 sd — correct |
| `noise_floor_augmented` **+0.0456 / +0.0669** | prediction-side Type I, **depth 40** | **p95 +0.04400 · p99 +0.06377** | +1.64 sd — high |

A 4,000-draw p95 at n = 14 has **±0.0009** of Monte-Carlo error. That, and nothing else,
separates +0.0453 here from `FINDING_034`'s +0.0431 on the *identical* population, and
+0.0456 here from `FINDING_036`'s +0.0435 on the identical population. The depth-40 floor to
use is **+0.0440**.

**The verdict strengthens.** The augmented gain **+0.0799** against the authoritative +0.0440
gives an exact one-sided **p = 0.0020** (2 × 10⁶ draws).

## 2. The two 20-pose baselines (playbook C8) — **two estimands, both right**

- **0.5270** — `argmin(xeng)` over **the 20 Modal poses that exist**. Deterministic, no draw.
  Recomputed exactly: **0.5270**. This is the primary endpoint's baseline and it is the
  correct one, because a purchase is measured against the pool that existed before it.
- **0.5291** — **E**[ `argmin(xeng)` over 20 poses drawn from the **40-pose union** ], a pool
  half of which is the new Explorer arm. Exact closed-form value **0.5293**.

Different populations, so this was never an inconsistency — but the two should not appear in
one sentence about "the" conversion.

## 3. The selected column peaks at pool 30 and falls at 40 — **real, and it is 08J**

Unremarked in the body. Both curves here are 256-draw Monte Carlo and both are exactly
computable in closed form. Computed exactly, with **zero sampling error**:

| pool | 20 | 24 | 28 | 30 | **33** | 36 | 38 | **40** |
|---|---|---|---|---|---|---|---|---|
| **exact** selected | 0.5270 | 0.5349 | 0.5400 | **0.5415** | **0.5424** | 0.5418 | 0.5406 | **0.5386** |
| published (256 draws) | 0.5270 | 0.5354 | 0.5393 | 0.5411 | | | | 0.5386 |

The union curve peaks the same way: exact **0.5410 at depth 34**, 0.5386 at 40. Every
published rung is within ±0.0007 of its exact value against a ±0.0012 256-draw SE. **So the
non-monotonicity is not noise.** The oracle stays monotone throughout.

**But ten of the fourteen ligands do not move at all** — they have no new pose that beats
their incumbent on `xeng`, so depth cannot touch them. Exact, pool 30 minus pool 40:
**08J +0.1721**, CFF −0.1195, YNV −0.0177, PG0 +0.0059, the other ten exactly 0.0000. 08J has
one new pose that wins on `xeng` and scores **0.1364** where the incumbent scored 0.4806; at
pool 30 it is drawn half the time, at pool 40 always.

**What it licenses:** at n = 14 a single adversarial pose can invert the selected curve while
the oracle rises. That is this finding's own thesis at its smallest instance. It does **not**
establish a turning point in depth, and it is **not** a reason to stop buying depth — it is
one more reason the depth question cannot be settled on this stratum.
