# DROP-DAY PLAYBOOK — OpenADMET CYP3A4 structure track

**Execute this document on release day.** It is the consolidation of thirty-seven findings
into one ordered path. Every claim carries its finding number. Where something is not
measured, it says **not measured** — do not fill the gap under deadline pressure.

**Upstream state at time of writing (2026-09-23, `data/processed/upstream_state.json`):**
`STRUCTURE_TRACK_LIVE = False` · `STRUCTURE_DATASET_SIZE = 20` · dataset last modified
`2026-08-27T23:16:39Z` · 7 files, none of them a structure set. `20` replaced the earlier
`184` PXR placeholder but is still marked TODO upstream — **re-read it on the day.**

**Read before executing:** §9. It used to list eight gaps that made the printed path
unrunnable on a blind set; seven are now closed and the whole sequence below §1 has been
**rehearsed end to end on the 87 validation ligands with their crystals hidden**
(`python scripts/ops/blind_rehearsal.py run`, 2026-09-23). Unblinded, it reproduced the
board to four decimals — selected **0.6164**, oracle **0.6975**, random **0.5769**, gain
**+0.0395**, ρ **−0.2582**, correct sign **75.86%** — and picked the **same pose on 87 of
87** ligands as the crystal-dependent path. **All eight gaps are now closed:** G6, the
reference-pose volume, was measured on 2026-09-23 (FINDING 038) — 9 of 9 ligands clear
depth 4 *and* depth 6, and the remaining risk is a 6.25% per-job failure rate, not
chemistry. §9 now carries the recommended sweep and the schedule.

---

## 1. The executable sequence

Run the steps in order. Each has a *done* condition and a *check*. Do not advance on a
status word — exercise the thing (`README` §Traps 2).

### Step 0 — detect the release

```bash
python scripts/ops/readiness.py
```

* **Done:** exit code `2` and `*** UPSTREAM CHANGED ***` printed.
* **Check:** `STRUCTURE_TRACK_LIVE` flipped to `True`, `STRUCTURE_DATASET_SIZE` changed,
  or a new dataset file appears. Any one is the signal.
* **Then, by hand:** re-read the Space config. `STRUCTURE_DATASET_SIZE` and the example
  compound-id format (`x00011-1`) were both PXR placeholders. The file count is what the
  portal gates on first, with no diagnosis on failure, so get it right here.
* Download the test CSV to `data/processed/test_ligands.csv` with columns
  `id,smiles` (the validator accepts `structure` or `id`).

### Step 1 — preflight-parse the new ligands (CPU, no GPU allocated)

```bash
python scripts/cofold/preflight_parse.py \
    --csv data/processed/test_ligands.csv --tag drop --samples 20
```

* **Done:** every ligand parses; the script names any that do not.
* **Check:** organometallics. Boltz cannot parse Ir/Ru dative-bond SMILES — 12 were
  excluded from the validation set (`data/processed/excluded_organometallic.json`). If any
  appear in the test set they must be **reported as uncovered, not quietly dropped**, and
  Protenix will parse them where Boltz will not (RUNBOOK, venue facts).
* **Why this is a gate, not a formality:** it has already caught two batch-killers,
  including a `contact` constraint with wrong list arity that fails as *rc=0, zero
  structures, 25 seconds* — indistinguishable from success at a glance (RUNBOOK Step 0).

### Step 2 — generate the pool (Explorer; see §2 for the venue decision)

```bash
ssh explorer 'sinfo -p gyorilab,gpu -o "%P %t %G"'      # pick the partition list
python scripts/explorer/boltz_depth.py plan --csv data/processed/test_ligands.csv --tag drop
python scripts/explorer/boltz_depth.py conditioning     # prints the hashes that must match
python scripts/explorer/boltz_depth.py stage  --tag drop    # from a LOGIN node
python scripts/explorer/boltz_depth.py push             # sbatch templates, CRLF stripped
python scripts/explorer/boltz_depth.py verify           # inside a job, not on the login node
python scripts/explorer/boltz_depth.py submit --tag drop --seed 101 --samples 20 \
    --partition gyorilab,gpu
python scripts/explorer/boltz_depth.py poll
```

* `--csv` takes **every** id in the file; the 14-ligand Type-I filter is the FINDING 035
  pre-registration and applies only when `--csv` is absent. `--tag` flows through
  `plan`/`stage`/`submit`/`collect`/`score` on both sides and defaults to `stratum`, so an
  existing campaign is untouched.

* **Done:** `poll` reports `n_ligands × 20` mmCIFs on disk and boltz's own log says
  `Number of failed examples: 0`.
* **Check, in this order, and believe nothing until all four agree (FINDING 035):**
  1. the `srun` step's exit code (not the job state — see §8 T2),
  2. boltz's own failure count,
  3. the file count,
  4. 20 poses per ligand on every ligand.
* **Budget:** ~6.6 min/ligand at 20 samples on a V100; 14 ligands × 20 = **1 h 34 m**,
  0 failures, **$0**, 280/280 written, **280 distinct** md5 and 280 distinct at 0.05 Å in
  the heme frame — 20.0 per ligand (FINDING 035). Scale linearly; split into several
  independently resumable jobs if the set is large.
* `stage` re-hashes the alignment on the far side and reports `msa_md5_matches_local`.
  Do not proceed on a transfer that only claims to have worked (RUNBOOK_explorer_boltz §2).

### Step 3 — collect

```bash
python scripts/explorer/boltz_depth.py collect --tag drop
```

* **Done:** a flat directory of `<ligand>__s<seed>_m<k>.cif` under `POSE_DIR/drop_flat`.
  Confidence `.npz` blobs stay on `/scratch` — D: has ~12 GB and is not a place to land a
  pool.
* **Check:** the returned JSON now carries `planned`, `missing_ligands`,
  `poses_per_ligand` and `uniform_depth` against the plan manifest — that is check 4 of
  Step 2's four signals, and a file count alone never was. Also that the remote `find` was
  scoped to `./<tag>_s*` so a smoke run's poses were not swept in
  (RUNBOOK_explorer_boltz failure 8).
* **DO NOT run `boltz_depth.py score`.** It requires a crystal per ligand. It used to skip
  every blind ligand as `no_crystal` and return zero rows with exit 0; it now **refuses**
  a blind tag by name (G2, closed).

### Step 4 — build the reference set and the selection feature

The shipped selector scores a pose by its mean Chamfer distance, in the heme frame, to
poses of the SAME ligand from **independent** engines (FINDING 011). Those reference poses
do not exist for new ligands — they must be generated. See §9 G6; this is the longest pole.

```bash
# one submission WAVE per sampler setting - waves are the unit of independence, and
# --replicates is byte-identical on both engines (FINDING 015/016/034).
# The ligand order is SHUFFLED before --batch packing (seed 38, recorded in jobs.json)
# and --samples now defaults to 1. Both are defaults; you do not do anything.
python scripts/cofold/openprotein_cofold.py submit --engine protenix_v2 \
    --csv data/processed/test_ligands.csv --tag drop \
    --sweep 3x200,10x200,3x50 --batch 5
python scripts/cofold/openprotein_cofold.py submit --engine protenix \
    --csv data/processed/test_ligands.csv --tag drop \
    --sweep 3x200,10x200,3x50 --batch 5
python scripts/cofold/openprotein_cofold.py collect --engine protenix_v2 --tag drop
python scripts/cofold/openprotein_cofold.py collect --engine protenix    --tag drop

# IF `collect` REPORTS `"failed": n > 0`: run the IDENTICAL submit command again.
# That is the whole recovery. It retries the failed work in place, at the same wave
# index and the same setting, and touches nothing else.

# freeze the reference set. This REFUSES below depth 4 and names the ligands.
python scripts/cofold/openprotein_cofold.py refset --tag drop \
    --engines protenix_v2 protenix --csv data/processed/test_ligands.csv \
    --out data/processed/reference_set_drop.npz

python scripts/structure/build_xeng_feature.py --tag drop \
    --pool <flat pose dir> --pool-flat --ref-npz data/processed/reference_set_drop.npz
```

* **Done:** `data/processed/xeng_drop.csv` with `{ligand, sample, xeng}`.
* `refset` applies three filters and prints the count for each, because a filter that
  never fires beats any passing check (T8): **one pose per (engine, wave)** — within a job
  every diffusion sample shares one ligand conformation (FINDING 009); **md5 of the file
  bytes** — deterministic engines return byte-identical replicates (FINDING 015);
  **coordinates in the heme frame at 0.05 Å** — two files can differ in a header and be
  the same pose. Exercised offline on a constructed pool: 24 files → 6 dropped as
  same-wave, 9 as same-md5, depth 3, **refused**; and 18 files → 3 same-md5, 3
  same-coordinates-only, depth 4, written.
* **RECOVERING A FAILED JOB — re-run the identical command.** FINDING 038 measured a
  **1-in-16** job failure (opaque `internal server error` after 24.7 min) and recorded
  that `submit` could not recover from it: its resume set was `{(rep, ligand) for EVERY
  batch}` with no check of `b["done"]`, so the lost ligand-waves stayed claimed and the
  identical command printed "0 ligands". **Fixed 2026-09-23**; FINDING 038's §7 and the
  "fresh wave at a new setting" advice in its §10 item 4 are superseded by this bullet.
  The semantics now are:
  * a batch whose `done` is `"failed"` **releases** its `(wave, ligand)` claims;
  * a batch that is pending, collected, or already superseded **keeps** them;
  * anything already **on disk** for that wave keeps them regardless, so a half-landed
    job resubmits only the half that is missing.
  So the recovery is: **`collect` first** (a failure is only known once polled), then
  **re-run the same `submit --sweep ...` command**. It retries **in place at the same
  wave index and setting** — which is what keeps `refset`'s depth accounting
  `engines × settings`. Do **not** allocate a fresh wave to recover: reusing a setting at
  a new index makes `refset` count one opinion as two, and a genuinely new setting
  changes what the sweep means. Re-run the **same** `--sweep` string, too — wave indices
  are `enumerate(sweep)` positions, so a shortened sweep shifts them and trips the
  setting-collision guard (which now says so).
  * The retry batch records `retry_of` / `retried_ligands`, `submit` returns
    `claims_released` / `ligand_waves_retried` / `failed_batches_superseded`, and the
    superseded failure moves out of `collect`'s `failed` count into `superseded` — so
    `failed` means *outstanding*, not *historical*.
  * Exercised on FINDING 038's own pool with the failed job reconstructed: 68 poses on
    disk, depth **7** on `A1APA/PG0/MWS/ERY` and 8 on the rest. Old code: **0 jobs**.
    New code: **1 job, 4 ligand-waves retried**, and after it lands, depth **8 on 9 of 9**
    with `below_min_depth: []`. A further re-run submits 0. Pending batches are untouched.
* **The csv is shuffled for you before `--batch` packing.** `--batch` chunks in order, so
  an ordered csv aligns job boundaries with strata: in FINDING 038 all three Type I picks
  sat in chunk 2, chunk 2 was the job that failed, and the result was a Type I "depth
  deficit" at p = 0.0369 that was pure packing artifact. `submit` now shuffles with a
  seed derived from `(tag/engine, wave, --shuffle-seed)`, default **38**, recorded in
  `jobs.json` — so the two engines and the successive waves pack *differently* and a
  failure cannot take the same stratum out twice. On that 9-ligand csv: file order puts
  all 3 Type I in one chunk in **4 of 4** waves; seeded, **1 of 4** (9.2% of chunks over
  200 seeds). The same seed reproduces the same packing exactly. `--no-shuffle` restores
  file order **for debugging only**.
* **`--samples` now defaults to 1.** `refset` keeps exactly one file per `(engine, wave)`
  and always takes the `__r<w>s0.cif` file, so the rest are never read: 340 downloaded →
  272 discarded unread (**80%**) in FINDING 038, because `diffusion_samples` does not
  sample the ligand (FINDING 009). Nothing else consumes them — the only other artefact is
  `confidence_<tag>_<engine>.csv`, which no script in this repo reads. At `--samples 20`
  over 100 ligands that is ~5 GB fetched to keep ~250 MB on a D: with 12 GB free. Pass
  `--samples 20` explicitly if you want the confidence block.
* `--pool-flat` is what Explorer's `collect` produces. Without it `build_xeng_feature`
  globbed for per-ligand *directories* and read zero poses, reporting it as a `--pattern`
  problem (G4, closed).
* **Check — this is the one that decides whether the selector is usable:**
  `build_xeng_feature` prints `reference depth` and **REFUSES below 4 independent poses
  per ligand**. At one reference pose the feature measures **−0.0055**: a thin reference
  set does not give a weaker selector, it gives a **harmful** one (FINDING 011).
* **Both OpenProtein engines are deterministic** (FINDING 015): a 12→24 replicate doubling
  moved the oracle on **0 of 489 pairs**. `--replicates` buys nothing at any count. Depth
  comes from **submission waves** (FINDING 034 — distinct poses track waves, and replicates
  16–29 returned 196 of 196 byte-identical) or from the **sampler**
  (`num_recycles`/`num_steps`, FINDING 016 — 4 distinct placements from 4 settings, 1.2%
  catastrophic; never `num_recycles=1`, which degrades by −0.0596).
* **If depth cannot be reached:** pass `--skip-thin`. Those ligands fall back to the
  FINDING 003 selector (+0.0265) rather than losing the feature for every ligand. Record
  which ligands fell back and report it with the entry.
* **More engines is NOT better.** Best measured reference set is the two Protenix
  checkpoints (+0.0380); adding esmfold2 at matched depth drops it to +0.0178 (FINDING 011).
* **Never add a reference engine to the POOL.** A union pool adds +0.0375 of oracle that
  selection cannot reach (FINDING 013), reproduced by a second mechanism in FINDING 016
  (+0.0355 oracle, −0.0038 selected).

### Step 5 — type the test ligands (do this before selecting; §3)

```bash
python scripts/ops/pool_diagnostics.py --xeng data/processed/xeng_drop.csv

python scripts/structure/binding_mode_robustness.py predict \
    --pool <flat pose dir> --pool-flat \
    --ligands data/processed/test_ligands.csv \
    --out data/processed/binding_mode_pred_drop.csv
```

The rule, from FINDING 033: **median `fe_donor_dist` over that ligand's own poses, cut at
2.6 Å**. Needs no crystal. Re-run on the validation pool it reproduces
`pred_fe_donor_median` to 9e-16 and the crystal agreement to the ligand — **84 of 87 =
96.6%**, composition 73 Type II / 13 Type I / 1 peripheral, i.e. the 14-ligand
non-coordinating stratum (G7, closed).

* **Done:** a per-ligand `pred_mode` and the composition (n Type I / n Type II).
* **Check:** the label is computable from the prediction alone on **84 of 87 = 96.6%**;
  median fraction of poses coordinating is **1.00** for Type II and **0.00** for Type I;
  class separation AUC = 1.00 with the closest Type I sitting 0.244 Å beyond the furthest
  Type II (FINDING 033). A ligand landing inside [2.59, 2.82] Å is unlabelled by this rule
  in either direction — flag it, do not guess.

### Step 6 — select

```python
from cypstruct import xengine as X
import pandas as pd
picks = X.select(pd.read_csv("data/processed/xeng_drop.csv"))   # -z(xeng), one row per ligand
```

* **Done:** one `sample` per ligand.
* **Check:** `argmin(xeng)` is exactly `cypstruct.xengine.select()` — verified at max abs
  difference **0.0** on all 87 in FINDINGs 035, 036 and 037 independently.
* **Do not tune it.** Unweighted and alone on purpose: adding the sibling-consensus term
  measured **worse** at full depth (+0.0183 against +0.0336) and fitting weights collapses
  it into noise (+0.0149) (FINDING 011). The board it reproduces: selected **0.6164**,
  oracle **0.6975**, random **0.5769**, gain **+0.0395**, ρ **−0.2582**, correct sign
  **75.86%**.

### Step 7 — build the submission

```bash
python scripts/submit/build_submission.py build --blind \
    --tag drop --arm unsteered --pool-dir <flat pose dir> --pool-flat \
    --ligands data/processed/test_ligands.csv \
    --heme keep --out submissions/02_drop.zip
```

* **`--blind` is the drop-day path and it is EXERCISED** (G5, closed). It calls
  `choose_poses_blind()`, which selects with `cypstruct.xengine.select()` over
  `xeng_<tag>.csv` alone — no `poses_scored_*`, no `consensus_features_*`, nothing from a
  crystal — asserts `argmin(xeng) == select()` before using either, and writes the picks
  and the per-ligand rule into `submission_report_<tag>.json`. Run without `--blind`,
  `build` behaves exactly as before.
* **Ligands that fell back are NAMED, not absorbed.** A ligand missing from
  `xeng_<tag>.csv` (one `--skip-thin` dropped) takes the FINDING 003 rule from
  `orientation_features_*` + `consensus_features_*` — both prediction-side — and the
  report lists it under `fell_back`; a ligand with no pose at all is listed under
  `uncovered`. Verified by thinning 10 ligands' reference sets to depth 1: all 10 were
  named, all 10 got a **different** pose from the one the full feature picks, and the
  other 77 were **bit-identical** to the full run.
* **Format — the OFFICIAL spec**, read from the Space's own `submission.py`, not inferred:
  a flat `.zip` of exactly `STRUCTURE_DATASET_SIZE` `.pdb` files, named
  `<compound_id>.pdb`, each a **full protein–ligand complex with the ligand residue named
  literally `LIG`**.
* **Keep the heme.** "Full protein–ligand complex" resolves it; `--heme keep` is the
  default and is now the documented answer, not a judgement call.
* Boltz emits `LIG1`, Chai emits `LIG3`. A zip of raw co-folder output is **rejected**.
  The converter renames and then asserts it.
* **Done:** `wrote N PDBs`, `errors: 0`, `files without exactly one LIG residue: 0`.

### Step 8 — validate before uploading

```bash
python scripts/submit/build_submission.py validate \
    --zip submissions/02_drop.zip \
    --ligands data/processed/test_ligands.csv \
    --expect-n <STRUCTURE_DATASET_SIZE>
```

* **Done:** `VALID: all checks passed`.
* **Check:** file count (the portal refuses on count **before any other check**, with no
  diagnosis), one `LIG` residue per file, and the ligand heavy-atom count against the
  requested SMILES. Heavy atoms are counted by atomic number, not via `RemoveHs` — some
  challenge SMILES carry explicit stereo-defining hydrogens (metformin/MF8) and `RemoveHs`
  keeps them, which would reject a good pose and send you to debug the co-folder.

### Step 9 — report alongside the entry

* **Pool oracle before selection, always.** A selection number without its ceiling is
  uninterpretable (RUNBOOK standing rule). On a blind set the oracle is unavailable — say
  so explicitly rather than omitting it.
* Report: binding-mode composition (§3), pre-announced expected absolute score (§3),
  reference depth per ligand and which ligands fell back to the FINDING 003 selector,
  any uncovered organometallics, and the pool regime from `pool_diagnostics.py` (§7).

---

## 2. Generation: where from, in what order of preference

### 1st — **Explorer (SLURM).** The tested path.

`docs/RUNBOOK_explorer_boltz.md` is the full procedure. Exercised end to end 2026-09-22;
FINDING 035 measured it: Boltz-2 **2.2.1** offline, checkpoint md5
`2f0a1775bf8fc366a1a85e2019eca288`, staged **6,979-sequence** MSA md5
`6de0ee1330ea9aa1aef1efb6ec6a2e3d`, heme bond `[A,442,SG] ↔ [H,1,FE]` asserted by Boltz's
own parser 14/14 — **280 poses at 20/ligand, 0 failures, ~6.6 min/ligand, $0, all 280
distinct.** The new arm's pool mean was **−0.0019** from the pool it joins: the conditioning
match is real, and free.

**MSA staging is not optional.** GPU nodes have no direct internet. Every input must be on
disk before the job starts; login nodes do have internet — **fetch there, compute there
never.** `--use_msa_server` is structurally incompatible with this and cost the PXR campaign
~12 days.

The sbatch template (`scripts/explorer/boltz_depth.sbatch`), with the load-bearing lines:

```bash
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=07:00:00
set -euo pipefail
export http_proxy=http://10.99.0.130:3128
export https_proxy=http://10.99.0.130:3128

srun "$FT/env/bin/boltz" predict "$ROOT/yaml/$SUB" \
    --checkpoint "$FT/boltz_cache/boltz2_conf.ckpt" \
    --cache "$FT/boltz_cache" \
    --out_dir "$OUT" --output_format mmcif \
    --diffusion_samples "$SAMPLES" --recycling_steps 3 --sampling_steps 200 \
    --seed "$SEED" --num_workers 4 --no_kernels --override
```

* `--partition` is passed on the command line, **never baked in**. `gyorilab`'s four H200s
  are routinely 4/4 allocated to multi-day interactive sessions. `gyorilab,gpu` lets Slurm
  take whichever frees first; the smoke test landed on a general-partition V100 (`d1010`)
  in under a minute.
* `--time=07:00:00` is a consequence of that fallback: the general `gpu` partition caps at
  8 h and silently refuses anything longer, which would pin the job to the busy partition.
* `--no_kernels` is **required**, not an optimisation. `boltz predict` defaults to the
  cuequivariance path and this env's `cuequivariance-*-cu12` wheels die on
  `libnvrtc.so.12: cannot open shared object file`. Two fine-tuning jobs died in <70 s on it.
  It is the one genuine mismatch against the Modal pool (confound K1) and FINDING 035
  measured it as immaterial: new-arm pool mean 0.4577 vs existing 0.4597.

**The five documented failure modes and their fixes** (RUNBOOK_explorer_boltz, plus 035):

| # | symptom | cause | fix |
|---|---|---|---|
| 1 | **Job reports `FAILED` with all 280 poses written.** `srun` step exited `0:0`; batch step exited **13** | the trailing `find … \| head -3` under `set -o pipefail`: `head` closes the pipe, `find` takes SIGPIPE — *after* every pose is written. Job `10529025` did exactly this and it would have thrown away a **94-minute** run (FINDING 035) | `\|\| true` on the summary `find`. **Already fixed in the template.** Before believing a FAILED state, check all four signals in Step 2 |
| 2 | `rsync: dup() in/out/err failed`, 0 bytes, on every invocation | MSYS2 rsync ↔ OpenSSH stdio incompatibility on this box. Nothing to do with the cluster | **rsync is unusable from here.** `scp` for single files, `tar czf - \| ssh … tar xzf -` for directories |
| 3 | staging returns `ok: false`, **empty stderr**, null return code | `subprocess` ran bare `bash`, which resolves to `C:\Windows\System32\bash.exe`, the **WSL launcher**; it answers on **stdout** and exits 1 | pin the absolute Git Bash path (`_git_bash()` in `boltz_depth.py`) and refuse to fall back |
| 4 | `#!/bin/bash^M: bad interpreter` | CRLF from the Windows checkout | `push` runs `sed -i 's/\r$//'` far side and prints `head -1 \| cat -A` to prove it |
| 5 | `MSA file cyp3a4.a3m not found` although it sits next to the YAML | Boltz resolves a relative `msa:` against the **process cwd**, not the YAML's directory | the YAML carries the **absolute** `/scratch/...` path |

Also: `echo "rc=$?"` after a pipe captures the *last* command's status, not the transfer's —
use `${PIPESTATUS[0]}` and **count files on the far side** after every push.

### 2nd — **nothing.** There is no tested second venue.

* **Modal is over its spend cap.** $37.58 billed, launches blocked. It built the validation
  pool (`--tag val87 --samples 20 --seeds 1`) and is kept only for provenance.
* **OpenProtein cannot serve Boltz-2 at matched conditioning.** `boltz2` **rejects an
  uploaded MSA** — verified on a stored probe: `boltz_2__msa → JobStatus.FAILURE, "internal
  server error"` while single-sequence `boltz2 → JobStatus.SUCCESS` (FINDING 034). The only
  depth it can sell is single-sequence, whose **ceiling matches and whose average is 0.057
  lower**; mixed into the pool it adds **+0.0118 oracle and −0.0065 selected** — the fifth
  time an expansion raised the ceiling and lowered the score (013, 016, 027, 031, 034).
* OpenProtein remains the venue for the **reference** engines (Step 4), never for the pool.
* Single-sequence is disqualifying for unfamiliar P450s generally: on the 4 pairs with both,
  MSA gives pool mean 0.5825 / oracle 0.8760 with **zero** catastrophic poses; single-sequence
  gives 0.0505 / 0.0505 with **75%** catastrophic (RUNBOOK).

### How much depth to buy

**Very little, and only at matched conditioning.** FINDING 035 ran the true 20→40 doubling
on Explorer: Δ oracle **+0.0097 [+0.0012, +0.0195]**; Δ selected **+0.0116 [−0.0687,
+0.1034]**, **Wilcoxon p = 0.715**, **10 of 14 unchanged**, and **+0.0020** once the two
swing ligands are dropped. The selected rate decays with depth (**+0.0227** at 10→20,
**+0.0095** at 20→40) while the oracle rate barely does (+0.0385 → +0.0331). Twenty samples
per ligand is the right default; a second seed is cheap insurance, not a lever.

---

## 3. What to do about the test ligands' composition

**Type the test set immediately — before selecting, before submitting.**

**The rule (FINDING 033, item 1):** median `fe_donor_dist` over that ligand's own pool
poses; **≤ 2.6 Å = Type II (iron-coordinated), > 2.6 Å = Type I (active-site)**. Needs
nothing from a crystal. Accurate on **84 of 87 = 96.6%** of the validation set.

**Report the composition with the entry, and pre-announce the expected absolute score.**

| | CYP3A4 Type II (n=73) | CYP3A4 Type I (n=14) | cost |
|---|---|---|---|
| random (pool mean) | 0.5985 | **0.4644** | −0.134 |
| **oracle** | **0.7068** | **0.6488** | **−0.058** |
| selected (shipped) | 0.6376 | **0.5055** | −0.132 |

**A Type I-rich set scores nearer 0.51 than 0.64, and over half that shortfall is
unreachable by any selector** — the oracle itself is 0.058 lower. That is a *generation*
exposure, not a selection failure. The fraction of available headroom captured is
statistically indistinguishable between strata (CYP3A4 36% vs 22%, difference CI
[−0.17, +0.44]; family 65% vs 58%, CI [−0.08, +0.23]).

**Do not change the selector on the strength of the composition.** The only adequately
powered test — 342 family pairs — puts the headroom-matched Type II − Type I difference at
**−0.0030, CI [−0.0132, +0.0075]**, with both strata clearing their own nulls at p = 0.0000.
CYP3A4 alone is underpowered: n = 14, CI [−0.063, +0.050], MDD 0.066. The heme-frame
explanation fails (partial ρ = −0.033). Changing a selector on an n=14 stratum is exactly
the post-hoc pick FINDING 011 retracted `-zm - zx` for.

**Do not buy Type I depth from a degraded configuration.** FINDING 033's own item 3 said to;
**FINDING 034 replaced it** — buy at matched conditioning or not at all — and FINDING 035
then measured the matched purchase as a near-null (§2).

**Watch QDY's failure mode.** On 1 of 14 Type I validation ligands the co-folder forced
coordination on a non-coordinator: all 20 poses at 2.23 Å against a crystal 2.82 Å. A ligand
whose poses all coordinate but whose scaffold has no accessible sp²/sp³ nitrogen donor
(`cypstruct.chem.coordinating_atoms`, 98.8% recall) is a flagged case. **Not measured as a
selector** — flag it in the report, do not act on it.

**On chemistry novelty (FINDING 017):** the released activity test set is interpolation —
median nearest-neighbour Tanimoto **0.587**, only **0.4%** of compounds below 0.4, 89% in
0.45–0.70. In that band the selector gains **+0.0427**, above its own overall +0.0383, and
it never collapses (+0.0227 on the most novel decile, above the +0.0138 floor). This removes
a worry; it is **not** a licence to expect more, because the gain tracks pool quality, not
chemical novelty (FINDING 012).

---

## 4. THE DEAD LIST

**Every lever below has been refuted. Do not re-run any of them under deadline pressure.**
One line each: what it was, why it died, its finding.

### Generation-side

| lever | verdict | mechanism, in one sentence | # |
|---|---|---|---|
| **Fine-tuning Boltz-2** | dead | monotonically negative from the smallest dose — true Δ **−0.0135 / −0.0369 / −0.0301** at 87/174/350 steps against a true base of 0.8110, so it is the lever that is wrong, not the recipe; composition cannot rescue it | 018 |
| **The Cys442→FE heme-bond constraint** | null | the **unbonded** predictions already reproduce the crystal's Fe coordination distribution, so there was nothing to fix: **Δ −0.0069 at p = 0.43** | 019 |
| **Holo templates** | dead | a template carries **no ligand and no side chains** (backbone frames + CB only), and blind holo over the F/G span 210–216 is **4.25 Å median / 7.20 Å p90** against the model's own **1.03 / 4.98** — 4.1× worse at exactly the residues that block | 024, 028 Add. 1 |
| **Apo templates** | dead | all six deposited apo entries are **one conformer** (apo–apo median 0.330 Å) and fail both pre-registered legs — **1.084 Å > 1.03 Å** median, **5.154 Å > 4.98 Å** p90; the model is closer on **46 of 58**, Wilcoxon p = 4.9e-08 | 028 Add. 2 |
| **Cytochrome b5 as a co-folded partner** | dead | it could only act by stabilising the protein, and the protein is already right to **0.73 Å** with ρ(pocket CA, ligand RMSD) = **+0.03**; the interface is 9.7 Å away on the **far side of the heme** and no CYP3A4–b5 complex is deposited | 024 |
| **Animal orthologs** | dead | **zero** ligand-bound structures exist for rat/dog/macaque CYP3A — the whole subfamily in a 493-pair harvest is human. There is nothing to run | 024 |
| **Two-ligand co-folding** (second copy of the query, or a decoy) | refuted | a second copy is **filler, not a constraint**: it lands in the active site on **12 of 15** and the scored copy gets **−0.066 LDDT-PLI** (CI [−0.142, +0.005]), rotation **+4.5°**; an unrelated second ligand is worse (−0.093, ejecting the query on 10 of 60) | 031 |
| **Pool expansion by rigid orientation grid** (1.78 M rotations) | dead | the model's own pocket **excludes the true pose** — the crystal ligand clashes below 2.2 Å on **71%** of poses where **0 of 87** crystals do in their own protein; oracle **+0.0108**, selection **−0.0145**, 31 of 41 failures unreachable by any rotation | 027 |
| **Pool expansion by degraded-conditioning depth** (OpenProtein single-seq) | refuted | the only depth the venue sells costs **0.057** of pool mean; mixed in it is **+0.0118 oracle, −0.0065 selected** | 034 |
| **Replicate jobs as a depth purchase** | dead | both OpenProtein engines are **deterministic** — a 12→24 doubling moved the oracle on **0 of 489 pairs**; replicates 16–29 returned **196 of 196** byte-identical | 015, 034 |

### Selection-side

| lever | verdict | mechanism, in one sentence | # |
|---|---|---|---|
| **Boltz's own confidence** | worse than random | within-ligand ρ between `complex_ipde` and true LDDT-PLI is **−0.092 / −0.033**, and selecting the most confident pose scores **−0.0417** against random | 001, 023 |
| **Protenix's confidence** | worse than random | seven confidence fields rank poses within a ligand at `frac(ρ>0) = 0.50` for four of five, every p > 0.6 — two independently trained architectures both at chance | 009 |
| **ATOMICA as a pose scorer** | loses | best feature **+0.0148** against a +0.0138 floor; the pre-registered ablation says why — it separates a synthetically rotated ligand at **AUC 0.95–0.97** and ranks real errors at chance. A **rigid-perturbation detector, not a pose-quality scorer** | worldmodel/EXPERIMENT_atomica_selector |
| **Boltz-2 trunk embeddings (`z`)** | null | the trunk runs **once per input, before diffusion**, so `z` is identical across every pose of the same ligand — pose-to-pose spread **1e-8**; it lands exactly on the ECFP4 ligand-only control (−0.0047, ρ 0.011) | 023 |
| **The CYP3A4 physics scorer** (34 terms above the iron) | null | a fixed vocabulary of pocket features averages over the one thing that varies; **not one term clears +0.0138**, and the fitted ensemble's +0.0252 is the **exact observed maximum** of its own 200-draw matched-dimension null, and loses to the incumbent by 0.0130 | 025 |
| **Induced-fit demand + all steric/clash/strain terms** | refuted | the repack probe fires on synthetic 90° rotations (p = 3e-142) and selects nothing (best single **+0.0012** against a +0.0140 floor); the clash terms are **anti-selective** (−0.0632) because the data prefer the snugger pose. **There is nothing left to try on the steric axis** | 029 |
| **Fragment placement transfer** from the P450 superfamily | refuted | **there is nothing to transfer** — donor fragment centroids scatter **4.48 Å** where the error to fix is 2.5 Å, the query's own crystal ranks at the **56th percentile** of its 20 predictions, and it is **−0.0464 paired with only 5 of 84 tied** | 030 |
| **Ligand torsion / internal-conformer priors** (superfamily, ETKDG, MMFF strain) | refuted | the prior is **sharp (25.0°)** and points the **wrong way** — CYP3A4 binds torsionally unusual conformers and the crystal ranks at the **31st percentile**, p = 2.6e-04; and the **torsion ORACLE is worth only +0.0213** and loses to the incumbent | 032 |
| **Top-k tie-breaking** (5 candidates × 3 k's) | refuted in substance | the near-tie is the **normal case** (median margin 0.0194 Å, 70 of 87 under 0.06 Å) and the stake does **not** concentrate in it (ρ = −0.14, p = 0.18) — a **variance mechanism, not a prize**. A **perfect** top-2 tie-break is **+0.0129** against a **+0.0134** floor; 15 of 16 configurations negative, best **+0.0011, p = 0.68, 71 of 87 tied** | 036 |
| **Minority-cluster rescue** ("the consensus is wrong, take the minority") | refuted in substance | the mechanism is **confirmed** on all four predictions and is still not a prize: the perfect rescue is **+0.0427**, its ceiling is a monotone function of the cut converging on the pose oracle at singletons, and the orientation-specific excess over a size-matched shuffle is **+0.004 / +0.015 median** against the +0.0134 floor. **Choosing the MAJORITY mode is worth +0.0546 (p = 1.3e-06)** — more than a perfect rescue of the minority — with a **−0.151** downside | 037 |
| **Anchor-local / population-level features** (~30 of them: occupancy and orientation priors, contact fingerprints, azimuth consensus, crystallographic contact prior, tier-1 donor prior, classical interaction energy) | all null | the iron anchor is **saturated** — 84% of poses already sit inside any reasonable coordination window — so the signal is in substituent placement, and nothing that scores a pose from its own geometry has ever beaten random here | 002, 006, 008, 010, 014 |
| **Splitting the shipped Chamfer into translation + orientation** | closed | the Chamfer was **already** ranking by orientation (ρ(mix,ori) = +0.879, same pose on 49 of 87); orientation alone is the strongest single unfitted term (+0.0424) and still **ties** paired (+0.0041, p = 0.61) | 026 |
| **A second engine as a pool member** | dead | architectural diversity has not bought decorrelation, 0 for 2 — Chai ρ = +0.45, Protenix ρ = +0.474–0.60 against Boltz's per-ligand oracle | 005, 013 |

**The pattern, five to six deep (CLAIM F):** oracle rises, selection falls — 013 (+0.0375
unreachable), 016 (−0.0038), 027 (+0.0108 / −0.0145), 031 (+0.0237 / −0.0027), 034
(+0.0118 / −0.0065). **A filter that improves the ceiling can make the selector worse.
Report both, always, and in that order.**

**What is left, from FINDING 037's close:** *"Every route that reads the pool's own
geometry — its mean (011), its higher moments (036), its modal structure (037) — is now
measured and spent. Anything that moves it must bring information the pool does not
contain."*

---

## 5. The two diagnostic rungs, and their scope

If a new idea appears on drop day, run the cheap diagnostic **before** building anything.
Both are free, both are leaky by construction, and neither ever enters anything scored.

### Rung 1 — the answer-recognition gate (R1), from FINDING 030

**How:** score the query's own **crystal ligand** as a 21st pose on the candidate feature,
and report where it ranks among that ligand's 20 predictions. Percentile = the fraction of
the 20 predicted poses whose feature value is **worse** than the crystal's.
**Pass rule, fixed a priori: mean percentile ≥ 0.65 and binomial p < 0.05 in the correct
direction.**

**⚠ SCOPE — this is the correction that matters (FINDING 036 §3):**

* **Apply R1 to a PRIOR** — a claim that some external record says where the ligand goes.
  That is what 030 and 032 both were; **both verdicts stand unchanged.**
* **Do NOT apply R1 to a within-ligand COMPARATOR** — a claim that these predictions can be
  ranked against each other. **The shipped selector fails the gate**: crystal at the
  **34th percentile**, median 0.20, 28 of 87, **binomial p = 0.0012**, while scoring
  **+0.0395** on the very pool it is gated on. The co-folders all share CYP3A4's 30° error
  (024; Chai ρ = +0.45, Protenix ρ = +0.47–0.60), so the consensus is **displaced from the
  truth while remaining informative about the ordering** — a biased but informative estimator.
* **Always score the incumbent on the same gate first, and read the branch that gives.**

### Rung 2 — the term oracle, from FINDING 032. **The safer instrument.**

**How:** hand the answer to the feature's **own scoring function** — score each pose by its
distance to the *true* value of whatever the term measures — and select with it. That is the
ceiling of every possible version of that feature, under every prior, weighting and
calibration.

**What it bounds:** 032's torsion oracle came in at **+0.0213**, at the noise floor's p99
and below the incumbent — which closes the internal conformer as a selection term *no matter
what prior is used*. 037's cluster oracle (+0.0427) and 036's perfect top-2 tie-break
(+0.0129) are the same instrument.

**The gate says "this prior fails"; the term oracle says "no prior can succeed."** It is
what closed 036 before a single candidate was scored, and it is unaffected by the R1 scope
correction. **Use it first.**

---

## 6. The statistical bar — stated once, not to be relitigated

A candidate ships only if it clears **both** bars, in order.

**Bar 1 — beat a null recomputed on THAT pool at THAT n.** Not an inherited constant.
Measured floors (random-feature selection, 4,000 draws):

| pool | n | p95 | p99 |
|---|---|---|---|
| CYP3A4 full | 87 | **+0.0134** | +0.0193 |
| CYP3A4 Type-I stratum | 14 | **+0.0435** | +0.0623 |

The n=87 floor has now been re-derived from scratch four times and lands on **+0.0134 to
+0.0141** every time (007, 025, 026, 029, 036). **Anything under +0.020 at n=87 is noise.**
*Caveat: the n=14 floor is quoted as +0.0435 (036), +0.0453 (035) and +0.0456 (035) in
different places — see §10 C1. Use the largest when judging a candidate.*

**Bar 2 — a PAIRED test against the incumbent**, reporting Δ, bootstrap **CI**, **Wilcoxon
p**, and the **tie count**. The tie count is not optional: 036's best candidate was
+0.0011 at p = 0.68 with **71 of 87 tied**, and 029's ensemble was −0.0167 with 5 of 87 tied.

**Nine-plus candidates have cleared Bar 1 and died at Bar 2.** 025's ensemble (+0.0252,
loses by 0.0130), 026's orientation term (+0.0424, ties at p = 0.61), 029's 44-column
ensemble (+0.0228, −0.0167 paired), 030's fragment consensus, 032's torsion family, 034's
depth purchase, 036's five tie-breakers, 037's minority rescue, 023's fitted z-confidence.

**Bar 3, where applicable — complementarity.** The candidate's per-ligand gain must correlate
with where the incumbent errs. 025 died at r = −0.047; 029 died at r = −0.078. Anything
absorbed by the incumbent is a **negative** — log it.

**Fix every sign a priori.** Leave-one-out sign selection manufactures confident negatives
from nulls: one term came back at ρ = −0.316 with 0 of 57 ligands positive and **p = 9e-09**
and its raw mean ρ was **+0.0005 on 51% of ligands** (025).

**Do not read ρ as selection value.** Rank correlation and top-1 selection have moved in
opposite directions **three** times (011, 012, 026/027, 036). Judge on the metric that will
be used.

**Report the pool oracle before any selection number, always.**

---

## 7. The one thing that is still live

**Our validation pool is unimodal. OpenADMET states several ligands fit multiple mutually
exclusive conformations. Those two facts have never met.**

What FINDING 037 established: the clustering is *not* uninformative — the orientation-specific
excess over a size-matched shuffle is **positive at 11 of 14 cuts** and reaches **+0.0191**,
which is at the +0.0134 floor, not below it. And, quoted:

> *"A pool with genuinely bimodal orientations — the `FINDING_012` catastrophe regime, which
> the blind cryoEM release may be, and which OpenADMET's 'multiple mutually exclusive
> conformations' explicitly describes — would have more clusters and a larger stake. **There
> is no drop-day proxy for it and none is claimed.**"*

**Be explicit about that last sentence. No drop-day proxy is claimed here either.**

### What to measure on the real pool

```bash
python scripts/ops/pool_diagnostics.py --xeng data/processed/xeng_drop.csv
python scripts/structure/minority_report.py cluster      # ⚠ hardwired to val87b — see §9
```

* `pool_diagnostics.py` **already reports the regime**: median within-ligand relative spread
  of `xeng`, interpolated between exactly **two** measured anchors — spread **0.51 → gain
  +0.0381** (CYP3A4 Boltz, 0.1% catastrophic) and spread **0.83 → gain +0.1448** (P450
  protenix, 12% catastrophic). It prints its own calibration next to the answer because two
  points define a line and nothing more. It does **not** extrapolate far beyond them, and
  its first version was wrong precisely because its thresholds were guessed before either
  anchor was measured.
* **Expect the +0.04 regime for a well-behaved release.** CYP3A4 wide-spread ligands gained
  +0.0704, narrow-spread ones +0.0063 (null) — the gain is a **catastrophe rate**, so it
  falls as pools improve (FINDING 012: 19.3% → 7.74% catastrophic giving +0.3006 → +0.0697).
  A falling number there is the mechanism working, not a regression.
* **Cluster-count comparison, if you want the regime answer:** cluster each ligand's poses
  in its own heme frame, **average linkage**, pose-to-pose **Chamfer** at **1.00 Å** (the
  headline cut) and **rotation at 15°**. Our unimodal pool gives **3.55** and **4.13** mean
  clusters respectively. A materially higher count on the blind pool is the signal that this
  is a different regime. **This comparison is a description, not a licensed action:** no rule
  keyed on it has been measured, and `minority_report.py` is hardwired to the validation set.
* **Two prediction-side quantities have been tested as regime indicators and FAILED.** The
  top-2 margin distribution does not predict the stake (ρ = −0.14, p = 0.18) and pool
  tightness points the *wrong* way (TIGHT pools are wrong 20.7% of the time, WIDE 29.3%).
  Do not revive either.
* **If the pool does look bimodal: still take the majority mode.** At Chamfer 1.00 Å,
  choosing the majority cluster is worth **+0.0546 (p = 1.3e-06)** against random cluster,
  more than a *perfect* rescue of the minority (+0.0427), with a **−0.151** downside if you
  choose wrong. `argmin(xeng)` already does this — it is above chance at **14 of 14** cuts.
  **Do not build a rescue rule on the day.**

---

## 8. Known traps, each of which has cost real time

**T1 — The residue-numbering offset that silently scores exactly 0.0 (FINDING 021).**
`lddt_pli` pairs protein atoms by residue number. Boltz numbers its output `1..N` from the
input sequence; a crystal uses auth numbering (CYP3A4 starts near 29; 3NA0 is offset by 43).
Where they disagree, every contact is compared against the wrong residue and the score
collapses to **exactly 0.0**. It failed silently on **28 of 85 pairs** and invented a 46%
catastrophe rate; the true base was **0.8110**, not 0.3477, and **1/84** catastrophes, not
39/84. FINDING 020 was retracted in full.
**Guard:** scan integer offsets and take the one maximising **residue-NAME agreement**, then
renumber. `best_offset()` / `shift()` in `scripts/finetune/score_holdout.py` (`--min-identity`,
default 0.80), `renumber_to_reference()` in `scripts/cofold/two_ligand_cofold.py`, or
`cypstruct.qmscore.pocket.resolve_offset()` against the canonical UniProt sequence for a
**blind** target. **Emit `resnum_offset` and `seq_identity_at_offset` on every pair and refuse
anything under 0.80.** Match on residue NAME, **never on overlap count** — that variant
produced **16.8 Å** CA RMSDs and was hit again inside one session after 021 was written.
**Any score of exactly 0.0000 is this bug until proven otherwise.**

**T2 — Explorer's login-node cgroup SIGKILLs long reads while reporting success.** A killed
`md5sum` of the 2.3 GB checkpoint can still exit 0 — a hash of nothing, reported as success.
**Never md5 or read a large file on the login node.** Fetch on login, verify in a batch job:
`boltz_depth.py verify` (CPU, `--mem=16G`, `short`, ~20 s) hashes the checkpoint, counts MSA
sequences and NUL bytes, prints the boltz version, checks the kernel import, and runs Boltz's
own YAML parser with structural asserts.

**T3 — `bash` resolves to the WSL launcher, with empty stderr.** `C:\Windows\System32\bash.exe`
answers every command with "Windows Subsystem for Linux has no installed distributions"
**on stdout** and exits 1, producing `ok: false`, empty stderr and a null return code — three
times before it was diagnosed. Pin the absolute Git Bash path and refuse to fall back.

**T4 — rsync is unusable from this box.** `rsync: dup() in/out/err failed`, 0 bytes, on every
invocation, through Python and typed into the shell alike — an MSYS2/OpenSSH stdio
incompatibility, nothing to do with the cluster. Use `scp` for single files and
`tar czf - | ssh … tar xzf -` for directories. And `$?` after a pipe is the *last* command's
status: use `${PIPESTATUS[0]}` and count files on the far side.

**T5 — a subsample curve's last rung is a deterministic maximum (FINDING 035).** The final
rung is the pool's own maximum with **zero variance**; every earlier rung is an *expected*
maximum over random subsets. So the final octave's slope is inflated and extrapolating it
overstates what new samples buy. This is exactly why FINDING 034's predicted **+0.0428** per
doubling met a **+0.0097** reality. It is a methods correction to 034 §2 and to 004.

**T6 — a job status, or a passing check, over broken work.** Three engines were written off
on a misread status or an unread `failure_message`; a readiness check reported **11/11** while
the submission was **unbuildable** because `pool_dir` was a dead parameter. **Exercise the
thing, do not stat it.** And its mirror: a job that reports FAILED with all its work complete
(§2, failure 1).

**T7 — counting artifacts instead of independent opinions.** 20 identical models read as 20
poses; 39% of "replicates" were duplicates; Protenix-v1 and v2 both fully deterministic; a
5,892-pose doubling turned out to be one file repeated. **Count distinct poses, never jobs —
and apply it to the POOL, not only the reference**, which is the half that went unchecked for
a whole campaign. Dedupe on geometry (`matched_depth_analysis.py distinct`, tol 0.05 Å).

**T8 — too-clean numbers.** An sd of exactly 0, or a filter that never fires, beats any
passing check. A seed-resampling sd of 0.0000 was **vacuous, not clean** — the `random_state`
only drove a binning subsample that does not trigger below 10,000 rows.

**T9 — disk.** C: ~19 GB, D: ~12 GB. Do not land a pool on D:; bulk stays on `/scratch`,
local intermediates go to `C:/cyp_struct` (`cypstruct.paths.SCRATCH`). **Never write bulk
data through the `O:` drive letter** — its VFS cache is unbounded and lives on C:. Use
`cypstruct.storage.push()`. Note `build_submission.py` hardcodes its temp dir to
`D:/cyp_scratch` in two places (§9 G8).

**T10 — background jobs die between turns.** Run long local work blocking or make it
resumable. Remote work (Slurm, OpenProtein) survives; local does not. And poll for
**advancing** progress, never for liveness — a 168-job batch was killed 2.5 minutes after
launch because it had produced no output directories yet, which is what a healthy run looks
like at that point.

**T11 — the released structures are a symmetric TRIMER.** Pair each ligand to the heme **in
its own chain**, or every distance is nonsense. Re-run reference-frame validation against the
new cryoEM structures **first**, before any cross-model metric (RUNBOOK Step 4).

---

## 9. What does not exist yet — read before Step 2

A code audit on **2026-09-23** found eight gaps that between them made the path above
unrunnable on a blind test set: `readiness.py` printed a sequence that routed through
Modal (over cap) and assumed crystals existed, and the runbook claimed a `plan --csv` flag
that did not exist. **All eight are now closed and exercised.** G6 — the one that needed a
network and a queue — was measured the same day (FINDING 038) and is below.

The acceptance test is in the repo and is re-runnable:

```bash
python scripts/ops/blind_rehearsal.py run      # the full sequence, crystals hidden
python scripts/ops/blind_rehearsal.py gates    # the REFUSALS, against broken input
```

`run` takes the 87 validation ligands, reduces them to `id,smiles`, re-presents the
existing Boltz pool in the flat layout Explorer's `collect` writes (hardlinks, not a copy —
T9), and executes plan → refuse-to-score → type → xeng → select → build → validate through
the same entry points drop day will use. Then it unblinds. Result, 2026-09-23:

| | expected | blind rehearsal |
|---|---|---|
| selected | 0.6164 | **0.6164** |
| oracle | 0.6975 | **0.6975** |
| random | 0.5769 | **0.5769** |
| gain | +0.0395 | **+0.0395** |
| ρ (within-ligand) | −0.2582 | **−0.2582** |
| correct sign | 75.86% | **75.86%** |

**All six to four decimals, and the same pose on 87 of 87 ligands as the crystal-dependent
path.** `argmin(xeng) == cypstruct.xengine.select()` is asserted inside
`choose_poses_blind` before either is used. Zero poses scored exactly 0.0000, which is the
FINDING 021 numbering control (T1). `validate --expect-n 87` passed; 87 PDBs, 0 errors, 0
files without exactly one `LIG` residue.

*One definition to keep straight, because getting it wrong makes a working pipeline look
broken:* **"correct sign 75.86%" is the fraction of ligands whose WITHIN-LIGAND ρ has the
expected negative sign (66 of 87)** — not the fraction of ligands where the pick beat the
pool mean, which is 67.8%. Both are in `blind_rehearsal.json`.

---

### G6 — reference poses for NEW ligands. **MEASURED. No longer the longest pole.**

**FINDING 038, 2026-09-23** replaced the one-ligand probe with the real thing: 9 validation
ligands, stratified in advance by heavy-atom tercile, rotatable-bond extremes and
prediction-side binding mode (3 of the 13 Type I), run through the shipped `submit --sweep`
(4 settings) → `collect` → `refset` on **both** Protenix checkpoints. **16 jobs, 25.7 min
wall clock, $0.**

**The rate: 9 of 9 = 100% clear depth 4, and 9 of 9 = 100% clear depth 6.** Median depth
**8**, minimum **7**, maximum **8**. `refset` wrote the `.npz` with `below_min_depth: []`.

**The mechanism, which is what makes this schedulable:**

> **depth = (engines) × (sampler settings), minus one per failed job containing the ligand.**

There is **no attrition from duplication at all**. 340 files → 272 dropped as same-wave
(FINDING 009, 4 of every 5 files discarded unread) → **`dropped_same_md5` = 0** and
**`dropped_same_coords` = 0**. Both zeros are real, not dead code: five positive controls
were run and all five fired (exact byte copy, identical vector, +1e-4 Å header-only change,
+5.0 Å genuinely distinct kept as 2, 272 same-wave). Because every wave is a *different*
sampler setting, FINDING 015's per-configuration determinism never gets a chance to bite.
Minimum pairwise Chamfer anywhere in the set is **0.288 Å** against a 0.05 Å tolerance —
~6× clear, not marginal passes.

**Depth does NOT depend on chemistry.** ρ(depth, n_heavy) = −0.087 (p = 0.82),
ρ(depth, n_rot) = +0.088 (p = 0.82), and pose spread is null on both plus mode (p = 0.90).
The apparent Type I deficit (7.00 vs 7.83, MWU p = 0.037) is **entirely** the one failed
job — see the trap below. Type I ligands are **not** harder to reference, so this does
*not* compound with FINDING 033's worse Type I pool.

**RECOMMENDED SWEEP — 2 engines × 3 settings, `--sweep 3x200,10x200,3x50`, `--batch 5`,
`--samples 1`.** Measured by re-deduping every sub-configuration of the poses on disk:

| sweep | jobs / 9 ligands | min depth | frac ≥4 | frac ≥6 |
|---|---|---|---|---|
| 2 eng × 1 setting | 4 | 1–2 | 0.00 | 0.00 |
| 2 eng × 2 settings | 8 | **4** | **1.00** | 0.00 |
| **2 eng × 3 settings** | **12** | **6** | **1.00** | **1.00** |
| 2 eng × 4 settings | 16 | 7 | 1.00 | 1.00 |
| 1 eng (protenix) × 4 | 8 | 4 | 1.00 | 0.00 |
| 1 eng (protenix_v2) × 4 | 8 | 3 | 0.56 | 0.00 |

Two settings is the **arithmetic** minimum — depth exactly 4, the refusal threshold with
**zero margin**. Three is the **viable** minimum: depth 6 survives one job failure
(6 → 5 ≥ 4) and it drops `3×400`, the slowest setting on both engines and the only one
that failed. The fourth setting is slack, not baseline: +33% jobs, +8 min, depth 6 → 8.
**A single engine cannot do the job.** (Every `0.56` and every low `min depth` in the grid
is the same four ligands from the one failed job, not a weak setting.)

**THE 6.25% JOB FAILURE — `submit` NOW RECOVERS FROM IT (fixed 2026-09-23).** 1 of 16 jobs
died with an opaque `internal server error` after 24.7 min at `progress_counter=75`. At the
time `submit` would not resubmit it: its resume set was `{(rep, ligand) for every batch}`
with **no check of `b["done"]`**, so all 4 lost ligand-waves stayed claimed and re-running
the identical command bought nothing. FINDING 038 §7 left that unfixed deliberately (it was
a measurement) and its §10 item 4 advises a fresh wave at a new setting. **Both are now
superseded:** failed batches release their claims and **re-running the identical command
IS the recovery**, retried in place at the same wave and setting. Full semantics in Step 4
above. Do not allocate a fresh wave — that would make `refset` count one opinion as two.
The 6.25% rate is still the reason the recommendation carries margin: recovery costs
another ~18 min of wall clock, and depth 6 means you do not have to spend it.

**THE CSV IS SHUFFLED FOR YOU BEFORE PACKING (fixed 2026-09-23).** `--batch` packed the CSV
in **file order**, so an ordered CSV aligned job boundaries with strata. Here the stratified
CSV put all three Type I picks in chunk 2 — and chunk 2 was the job that failed, taking out
**3 of 3** Type I ligands at once and manufacturing the p = 0.0369 above. `submit` now
shuffles per `(tag/engine, wave)` from `--shuffle-seed` (default **38**, recorded in
`jobs.json`); on this csv that drops "one chunk holds all 3 Type I" from **4 of 4 waves** to
**1 of 4**. `--no-shuffle` restores file order for debugging.

**Cost model — `jobs = 6 × ceil(N/5)` at the recommended sweep.** Concurrency measured at
**≥16 jobs with no queueing penalty**; above that unmeasured, so the schedule is a band
(lower = slowest single wave, 17.4 min; upper = serialised in blocks of 16 at ~18 min).

| test set | jobs (3 settings) | jobs (4 settings) | wall clock, optimistic → conservative |
|---|---|---|---|
| 20 ligands | 24 | 32 | ~20 min → ~36 min |
| 50 ligands | 60 | 80 | ~20 min → ~1 h 12 m |
| 100 ligands | 120 | 160 | ~20 min → **~2 h 25 m** |

Add ~5 min for `collect`, **under 1 min** for `refset` (2.9 s for 340 files). All fit the
2,000/month cap with **610 already spent**. **Even the conservative column is far shorter
than Step 2's pool generation (~6.6 min/ligand × 100 = 11 h), so reference depth is no
longer the schedule risk.**

**Use `--samples 1`, not the default 20.** 80% of downloaded files were discarded unread —
`diffusion_samples` does not sample the ligand (FINDING 009) and `refset` keeps one file
per (engine, wave). At `--samples 20` on 100 ligands that is ~5 GB downloaded to keep
~250 MB, onto a D: drive with 12 GB free.

Still true, and unchanged:

* `submit --sweep` records the setting beside the wave index and **refuses to reuse a wave
  index at a different setting** (two settings counted as one replicate is the T7 error).
  `--replicates` is kept but dead: byte-identical on both engines at any count (015/034).
* `refset` freezes the waves into an `.npz` and **exits 2** naming any ligand below depth 4.
* **If depth cannot be reached:** `refset --allow-thin`, then `build_xeng_feature
  --skip-thin`, then `build_submission --blind`, which names every fallback. Measured cost
  shape, from thinning 10 of 87 to depth 1: those 10 all got a different pose, the other
  77 were bit-identical, and the population difference between the two rules is 0.0395 −
  0.0265 ≈ **0.013 LDDT-PLI per fallback ligand** (at n=10 the realised difference was
  inside noise — do not read a single small stratum as the effect).
* **More engines is NOT better** (+0.0380 for two Protenix checkpoints, +0.0178 adding
  esmfold2), and **never add a reference engine to the POOL** (FINDING 013/016).
* **Not claimed by 038:** that 8 reference poses select better than 4. It measures supply,
  not value. n = 9, all CYP3A4, all known to fold; the 6.25% failure rate is **one
  observation** (95% CI ≈ 0.2–30%).

---

### Closed, with what closed it

**G1 — `plan` had no `--csv`. CLOSED.** `plan --csv <file> --tag <name>` takes every id in
the file; `--tag` flows through `plan`/`stage`/`submit`/`collect`/`score` on both sides and
defaults to the historical `stratum`/`smoke`, so an existing campaign is untouched. `plan`
writes `plans/<tag>.json`; `collect` reads it and now reports `missing_ligands`,
`poses_per_ligand` and `uniform_depth`. Exercised: 87/87 YAMLs from an `id,smiles` CSV, 0
skipped; the default `plan` still returns the same 14-ligand stratum.

**G2 — `score` returned zero rows on a blind set. CLOSED.** It refuses a tag whose plan
CSV has no `pdb` column, naming the CSV, instead of counting every ligand into
`skipped["no_crystal"]` and exiting 0. Exercised (exit 1, message quoted in
`blind_rehearsal.json`).

**G3 — three scripts in the printed path imported `modal` at module scope. CLOSED.**
`collect_and_score.py`, `orientation_features.py` and `test_consensus_selector.py` now
import it lazily, inside the branch that uses it, and all three import cleanly on a box
with no `modal` installed (verified). The two prediction-side ones take `--pool-dir`,
`--ligands` and `--no-eval` and were run against the local pool with no Modal and no
crystals. `collect_and_score.py` is importable but **has no blind mode and cannot get
one** — it scores against deposited structures by definition; that is G2, not a bug.

**G4 — `build_xeng_feature --pool` needed per-ligand DIRECTORIES. CLOSED.** `--pool-flat`
reads Explorer's flat `<lig>__*.cif` directly. The old workaround was a full copy of the
pool on a disk with ~12 GB free. It also prints `pool units` and how many have no
reference pose, so a silent zero-row read is no longer possible.

**G5 — `choose_poses()` could not run blind. CLOSED and EXERCISED.** `build --blind` calls
`choose_poses_blind()`: `xengine.select()` over `xeng_<tag>.csv` alone, the
`argmin == select` assertion, an explicit FINDING 003 fallback for thin ligands from two
prediction-side feature files, and `fell_back`/`uncovered` lists in
`submission_report_<tag>.json`. The ~10-line bypass sketched in the old §9 is superseded —
it had no fallback and no report.

**G7 — no prediction-side binding-mode labeller. CLOSED.**
`binding_mode_robustness.py predict --pool <dir> [--pool-flat]`. On the validation pool it
reproduces the stored `pred_fe_donor_median` to 9e-16 and the 84/87 = 96.6% crystal
agreement, and gives 73 Type II / 13 Type I / 1 peripheral. It **flags** any ligand whose
median lands in the empty band **[2.59, 2.82] Å** between the two observed classes rather
than guessing a side.

**G8 — two hardcoded `D:/cyp_scratch` temp dirs in `build_submission.py`. CLOSED.** Both
now use `cypstruct.paths.SCRATCH` (`C:/cyp_struct`) behind `guard_scratch()`, which raises
*before* the first write.

### Still stale, and cheap

`readiness.py`'s printed drop-day path (steps 3–5) still names `detached.py launch`,
`collect_and_score`, `orientation_features` and `test_consensus_selector` — the Modal,
crystal-dependent sequence. The correct sequence is §1 of this document. The readiness
*checks* are sound; only the five printed lines are wrong.

## 10. Contradictions between findings — resolve before citing

These are live inconsistencies found while assembling this document. They do not change any
verdict, but do not quote both sides of one in the same table.

> **Status, 2026-09-23 — `FINDING_039` closed C1, C2, C3, C5 and C8**, each rewritten below
> with its resolution. **C4, C6, C7 and C9 remain open**; all four are supersessions or
> reconciliations and none is decision-relevant on drop day. **One published claim changed
> sign** (FINDING 036 §2e/§5, see C1). Every affected finding carries a dated correction note.

**C1 — RESOLVED. The "n=14 noise floor" is three populations plus Monte-Carlo error.**
The five quoted values (+0.0431, +0.0433, +0.0435, +0.0453, +0.0456) are **one** estimator —
*the p95 of the gain from picking one pose per ligand uniformly at random, over that pool's
mean* — evaluated on **three different sets**, each estimated at 4,000 draws, which carry
**±0.0009** of slop at n=14. Authoritative, 2,000,000 draws, SE ±0.00004:

| population | **p95** | p99 |
|---|---|---|
| **crystal-side** Type I, pool A depth 20 — FINDING 033's | **+0.0431** | +0.0616 |
| **prediction-side** Type I, pool A depth 20 — 034, 035 `in_stratum` | **+0.0449** | +0.0640 |
| **prediction-side** Type I, **pool B depth 40** — 035 `augmented`, 036 | **+0.0440** | +0.0638 |

Crystal-side and prediction-side differ by one ligand (`D0R` in, `QDY` out) — FINDING 034's
own headline — and the random baseline moves 0.4644 → 0.4597. **Use the floor for the
population you are on, not "the largest".** Pooled at n=87 it is +0.0137. **Recompute on the
actual pool (§6), with ≥200,000 draws if a decision turns on it**: inside ±0.003 of a floor a
4,000-draw comparison is undecided.

**One sign flipped.** FINDING 036 §2e and §5 say pool B's perfect top-2 ceiling **+0.0447**
"lands **on**" the floor. Against the authoritative **+0.0440** it is **above**, exact
one-sided **p = 0.0475**. 036's *verdict* is unchanged — its stage-1 stop is the n=87
measurement (+0.0129 vs +0.0134) and its best realisable pool-B candidate is +0.0392, still
below the floor. 033's, 034's and 035's floor verdicts all stand (exact p = 0.0581 / 0.0073 /
0.0020). `FINDING_039` §1.

**C2 — RESOLVED. The measurement stands, the mechanism is refuted, the result is now
unexplained.** `fe_centroid_dist` scoring **+0.0408** on the arm4 holdout **reproduces to four
decimals**, and was always computed on post-fix scores (0 LDDT mismatches in 425 poses; 61 of
85 pairs at a non-zero residue offset, none below the 0.80 identity bar). The post-fix
catastrophe rate is **1 of 85** pairs / 7 of 425 poses — the 39 does not exist under any
definition. The attribution to FINDING 012's catastrophe scaling is not merely unverified, it
is **dead**: a *perfect* catastrophe-avoider is worth **+0.0002**; deleting all 7 catastrophic
poses leaves the gain at **+0.0409**; 83 clean pairs supply +0.0396 of the +0.0408. The
fallback account — pool uncertainty — fails and runs backwards (CYP3A4's pool has +0.1174
headroom and 34.2% sub-0.5 poses against the holdout's +0.0652 and 14.6%, and the feature is
worth +0.0013 there). **Still do not build on it** — now because nobody can say why it works.
`FINDING_039` §2.

**C3 — RESOLVED by pointer.** FINDING 028's body section *"LICENSED, and it reverses a call in
FINDING 024"* is dead: its own Addendum 1 withdraws it and Addendum 2 removes the last
qualifier. A dated warning box now sits at the **top of 028's body**, so a reader cannot reach
the superseded claim without passing it. Templates of every kind stay on the DEAD LIST.

**C4 — STILL OPEN (supersession, not conflict). FINDING 018 closes all fine-tuning arms; FINDING 022 conditionally reopens one.**
018: "arm1_promiscuous, arm2_family and arm3_cyp3a4_only stay built and unused … composition
cannot rescue a lever that damages the model at its smallest dose", reaffirmed by 021.
022: "**Re-examine arm3_cyp3a4_only** … worth one run *if* a difficulty-matched holdout can be
built for it", on a mechanism 022 itself calls "plausible, not demonstrated". **Arms 1 and 2
are closed; arm 3 is conditionally reopened. It is not a drop-day action.**

**C5 — RESOLVED. Two statistics on the same six pairs, not a disagreement.** 022's
per-target column is the mean of **`lddt_sample0`** — the engine's own rank-0 pose — and 024's
parenthesised scores are the mean over **all five samples**. Both reproduce exactly on five
targets and both aggregate rows: aromatase **0.8629 / 0.6313**, Q2IU02 0.9306 / 0.9241, CYP3A4
0.5550 / 0.5526, CYP51A1 0.9093 / 0.9076, CYP3A5 0.4758 / 0.5059. Aromatase is the only target
where they diverge materially (**0.232**, against ≤0.007 elsewhere): Boltz-2 puts a very good
aromatase pose at rank 0 and much worse ones behind it. 024's `entries` column is a **third**
denominator — crystals in the 406-entry cavity set (aromatase 7), not scored pairs (6).
**Prefer the all-5 mean for cross-target comparison**, because rank 0 is the engine's own
confidence ordering and FINDING 001 is that it does not rank poses. `FINDING_039` §4.

**C6 — STILL OPEN (reconciled in 037's own text). FINDING 037's "a bimodal pool would have more clusters and a larger stake" collides
with its own degeneracy result.** §5a establishes that the cluster oracle is **monotone in
cluster count by construction**, converging on the full pose oracle at singletons — which is
why "a single headline number for this ceiling does not exist". The honest quantity is the
**excess over the size-matched shuffle**, which is *not* shown to grow with cluster count: it
peaks mid-grid (+0.0158 Chamfer at 0.75 Å, +0.0191 rotation at 15°) and goes negative at the
loose end. §7 above is worded to that correction.

**C7 — STILL OPEN (supersession, not conflict). FINDING 035 item 5 recommends the experiment FINDING 036 then kills.** 035: "a
tie-break for the top-2 `xeng` poses is a better-evidenced next experiment than more depth."
036 measured it and said no. Supersession, not conflict — but 035's submission-time item 5 is
**stale as written; do not action it.**

**C8 — RESOLVED. Two estimands, and the peak is real but is one ligand.** **0.5270** is
`argmin(xeng)` over **the 20 Modal poses that exist** — deterministic, no draw, recomputed
exactly. **0.5291** is **E**[`argmin(xeng)` over 20 poses drawn from the **40-pose union**],
a pool half of which is the new Explorer arm — exact value 0.5293. Different populations, so
both are right; **the primary endpoint's baseline is 0.5270**, because a purchase is measured
against the pool that existed before it.

The non-monotonicity is **real, not Monte-Carlo noise**: computed in closed form, with zero
sampling error, the augment curve peaks at **pool 33 (0.5424)** and falls to 0.5386 at 40, and
the union curve peaks at **depth 34 (0.5410)**; every published rung is within ±0.0007 of its
exact value against a ±0.0012 256-draw SE. The oracle stays monotone. **But it is one
ligand** — ten of fourteen have no new pose that beats their incumbent on `xeng` and cannot
move at all; per-ligand the drop is 08J **+0.1721** net of CFF −0.1195. It does **not**
establish a turning point in depth and is **not** a reason to stop buying it: at n=14 one
adversarial pose can invert the selected curve while the oracle rises. `FINDING_039` §3.

**C9 — STILL OPEN (reconcilable, unreconciled in the text). FINDING 017's analogue ladder vs FINDING 022's "superfamily benchmarking is nearly
worthless."** Reconcilable but unreconciled in the text: 017's splits are **chemical** (ligand
novelty), 022's objection is to **protein** analogues chosen by family similarity. 024 R2
supplies the replacement criterion 022 lacked (cavity ≥ 600 Å³, target sub-2 Å rate near
CYP3A4's 33%). Treat 017's ladder as valid for **ligand-side** matching only.

**C10 — a non-contradiction worth recording.** The shipped board reproduces to four decimals
in 035, 036 and 037 independently — selected 0.6164, oracle 0.6975, random 0.5769, gain
+0.0395, ρ −0.2582, correct sign 75.86% — and `argmin(xeng)` equals
`cypstruct.xengine.select()` at max abs difference 0.0 on all 87 in all three.

---

## 11. One-page summary

1. `readiness.py` → exit 2. Re-read the Space config; **the file count gates first**.
2. `preflight_parse.py` on the test CSV. Organometallics are reported, not dropped.
3. Generate on **Explorer only**: `plan --csv <test csv> --tag drop`, then every stage
   with the same `--tag`. Stage the MSA on a login node; `verify` in a job; `sbatch`,
   never `srun`. **A FAILED state may be a complete run — check four signals.**
4. `collect --tag drop`, and check `missing_ligands` / `uniform_depth`. **Do not run
   `score`** — it refuses a blind tag, by design.
5. Buy reference poses for the new ligands:
   `openprotein_cofold.py submit --sweep 3x200,10x200,3x50 --batch 5` on **both** Protenix
   checkpoints, `collect`, then `refset`. The csv is shuffled before packing and
   `--samples` defaults to 1 — both automatic. Measured 9 of 9 at depth 6 (FINDING 038);
   `jobs = 6 × ceil(N/5)`, ≤2 h 25 m even at 100 ligands. `refset` and `build_xeng_feature`
   both **refuse below depth 4**; `--allow-thin`/`--skip-thin` rather than force, and
   **name the fallbacks**. **If `collect` reports `failed > 0`, re-run the IDENTICAL submit
   command** — failed batches release their claims and are retried in place at the same
   wave. Never allocate a fresh wave to recover.
6. Type the set: `binding_mode_robustness.py predict --pool <flat dir> --pool-flat`
   (median `fe_donor_dist`, 2.6 Å). Report the composition and pre-announce **~0.51 if
   Type I-rich, ~0.64 if not**. A median inside [2.59, 2.82] Å is flagged, not guessed.
7. Select with `cypstruct.xengine.select()`. **Do not tune it. Do not change it on an n=14
   stratum.**
8. Build the zip with `build_submission.py build --blind --pool-flat`: one flat
   `<id>.pdb` per compound, residue `LIG`, **keep the heme**. `validate --expect-n`.
9. Report the oracle before the selection number, or say it is unavailable.
10. **If a new idea appears: run the term oracle (§5 rung 2) first. If it cannot clear
    +0.0134 with perfect knowledge, it cannot clear it without.**
