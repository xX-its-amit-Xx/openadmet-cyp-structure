# Runbook — what to launch next, and what each step has to prove

Ordered so that each step's gate decides whether the next one is worth its runway.
`scripts/ops/watchdog.py` runs every two hours; the ops tick reads this file to decide
what to launch when nothing is in flight.

**Deadline: interim leaderboard 2026-09-24.** Budget backwards from that.

---

## Step 0 — the pool that scores itself  ✅ infrastructure ready

We have **87 CYP3A4 ligands with deposited crystal poses**
(`data/processed/validation_ligands.csv`, 72 iron-coordinated + 15 active-site; 12
organometallic ligands were excluded, see below). Every
pose we generate for these can be scored immediately with `cypstruct.pose.lddt_pli` and
`bisy_rmsd`. No waiting for the challenge release to start measuring.

**Always dry-parse on CPU first.** This is a gate, not a formality:

```bash
python scripts/cofold/preflight_parse.py \
    --csv data/processed/validation_ligands.csv --tag val87 --samples 20
python scripts/cofold/modal_boltz.py submit \
    --csv data/processed/validation_ligands.csv --tag val87 --samples 20 --seeds 1
python scripts/structure/score_pool.py --pool <collected dir>
```

The gate has already caught two batch-killers: a `contact` constraint with the wrong list
arity (which fails as **rc=0, zero structures, 25 seconds**, indistinguishable from success
at a glance), and **12 organometallic ligands** (Ir/Ru dative-bond SMILES) that Boltz
cannot parse. Those 12 are excluded and recorded in
`data/processed/excluded_organometallic.json`; they must be reported as uncovered, not
quietly dropped. The set is now 87 ligands, 72 coordinated and 15 active-site.

Watch the budget: 570 jobs at 10 samples estimates ~40 GPU-h and **the preflight will
refuse it** as more than half the monthly cap. More samples per job is cheaper than more
jobs, because Boltz runs the trunk once then diffuses: measured ~122 s for a 3-sample job
on an A100, and the running 168-job batch at 20 samples estimates 18.5 GPU-h.

**Gate:** does the steered arm beat the unsteered arm on LDDT-PLI, on ligands where the
donor prediction was correct? Report the **pool oracle** (best achievable in the pool)
alongside the selected score — a selection number without its oracle is uninterpretable.

**What would make me abandon steering:** if steered ≈ unsteered, the forced contact is
either being ignored or is landing on the wrong atom. Check the realised Fe–donor distance
distribution in the steered arm before concluding anything about selection.

---

## Step 1 — more engines, as tail rescuers only

```bash
python scripts/cofold/modal_af3.py weights          # one-time Kaggle pull into a Modal volume
python scripts/cofold/modal_af3.py submit --csv data/processed/validation_ligands.csv --tag val87
```

Then Chai-1 and Protenix on the same pattern. **Do not add a new engine as a z-hybrid pool
member.** In PXR, Protenix beat the pool on 6 of 8 holo cases and still regressed the board
0.5551 → 0.5241 as a member, while an 8-ligand tail swap gained +0.017.

```python
from cypstruct.select import z_hybrid, pool_confidence, sweep_tail_rescue
base, wins = z_hybrid(boltz_poses)
sweep_tail_rescue(base, af3_best, pool_confidence(all_poses), truth)   # sweep N, don't assume 8
```

**Gate:** does the rescue-depth sweep show an interior optimum? A flat curve means the new
engine adds nothing and should be dropped, not tuned.

---

## Step 2 — the physics scorer

Tier 0 (`cypstruct.qmscore.geometry`) is implemented and calibrated against the reference
set. Tiers 1–2 are specified in `docs/QM_SCORER_DESIGN.md`.

Order matters: tier 1 is **per ligand**, not per pose, so 99 ligands is 99 QM jobs even
though the pool is tens of thousands of poses. That is what makes this affordable.

**Gate:** beats the z-hybrid on **leave-one-scaffold-cluster-out** folds, and its gain
correlates with where the z-hybrid errs. Absorbed = negative = logged as one.

**Expected negative:** the strain term. PXR built an MMFF strain validator, tested it
against holo ground truth and killed it. If ours gates null too, that is a replication.

---

## Step 3 — fine-tuning, only if steps 0–2 have landed

See `docs/FINETUNE_PLAN.md`. Short version: LoRA on Boltz-2 first, Apheris low-N recipe
(LR 3e-4, warmup 50, EMA 0.99, ~350 steps), cluster-held-out splits, forgetting check on
unrelated targets. The PXR campaign spent twelve days here and got −0.0020; do not repeat
that without the earlier gates passing.

---

## Step 4 — the actual submission

The released structures are **split across the challenge train and test sets**, so part of
the release is scoreable ground truth and part is the blind target. When the data lands:

1. Re-run reference-frame validation against the new cryoEM structures **first**. They are
   a symmetric trimer — pair each ligand to the heme in its own chain, or every distance is
   nonsense.
2. Check whether the released ligands are type II or type I. Our steering and our scorer
   both lean on coordination; if the release is mostly substrates, say so and reweight.
3. Fit the selector on the released train portion; predict the test portion.

---

## Standing operational rules

- **Never write bulk data through `O:\`.** Its VFS cache is unbounded and lives on the
  near-full C:. Use `cypstruct.storage.push()` / `Batch()`.
- **Launch with `run_in_background`, not `nohup ... &`.** A backgrounded shell job dies
  when the tool call's shell exits — this already killed one Modal launch silently.
- **Do not pipe an unbuffered log through `grep` without `--line-buffered`.** It hides all
  output until the buffer fills, which looks exactly like a hung job.
- **Poll for ADVANCING progress**, not liveness - but check the run is OLD ENOUGH for
  "no progress" to mean anything first. Call `budget.may_judge_stalled(kind, started)`.
  A co-folding job takes about five minutes; a 168-job batch was killed 2.5 minutes
  after launch because it had produced no output directories yet, which is exactly
  what a healthy run looks like at that point. A liveness check that fires faster than
  one unit of work can complete does not detect stalls, it causes them.
- Every launch goes through `budget.preflight_hours` and is recorded in the ledger
  **before** it starts, so the watchdog can clean up even if the launcher dies.

---

## The whole ops tick, in four commands

```bash
python scripts/ops/watchdog.py            # runaways, real dollar spend, disk
python scripts/ops/readiness.py           # upstream state + drop-day path, EXERCISED
python scripts/ops/topup_p450.py --submit # learns skip lists, fills pool+reference depth
python scripts/structure/score_p450_pool.py score &&   python scripts/structure/score_p450_pool.py validate   # the FINDING 012 test
```

All four are idempotent and safe to run on an empty queue. Only the last is slow (it
rescores new poses), so run it when `topup` reports new pairs, not every tick.

**What "nominal" looks like right now:** watchdog 0 runaways with Modal blocked (expected -
everything is on OpenProtein), readiness 14/14, topup reporting 0 thin pairs, and validate
reporting a gain that *shrinks* as the pools improve. That last one is the counter-intuitive
part: a falling number there is the mechanism working, not a regression. It tracks the
catastrophe rate, which falls as pools deepen (FINDING 012 - eight measurements,
19.3% -> 7.74% catastrophic giving +0.3006 -> +0.0697).

---

## Every ops tick, run this first

```bash
python scripts/ops/watchdog.py        # runaways, real dollar spend, disk
python scripts/ops/readiness.py       # upstream state + drop-day readiness
```

`readiness.py` replaces polling the Space by hand. It reports `STRUCTURE_TRACK_LIVE`,
`STRUCTURE_DATASET_SIZE` and the dataset file list, diffs them against the last run, and
exits 2 when anything upstream moves. It also checks that every piece of the drop-day path
exists and that the selector still imports — a pipeline that has never been run is not a
pipeline.

**As of 2026-09-12: `STRUCTURE_TRACK_LIVE = False`, data unreleased, 11/11 readiness
checks pass.** `STRUCTURE_DATASET_SIZE = 184` is the PXR count with a `TODO` beside it, so
treat both it and the example id format as placeholders until the track goes live.

### Drop-day path, in order

1. **Re-read the Space config.** Dataset size and identifier format were placeholders.
2. `scripts/cofold/preflight_parse.py` on the new ligands — CPU, catches schema errors
   before any GPU is allocated. It has already caught two batch-killers.
3. `scripts/cofold/detached.py launch --engine boltz --arms unsteered` — detached so it
   survives the client; unsteered because FINDING 001 retired the steered arm.
4. `collect_and_score.py` → `orientation_features.py` → `test_consensus_selector.py`.
4b. **`build_xeng_feature.py --tag <tag> --pool <dir> --refs protenix_v2 protenix`**
   — writes `xeng_<tag>.csv`, which `build_submission.py` then prefers automatically.
   It REFUSES on a thin reference set rather than degrading, because at one reference
   pose the feature measures -0.0055: it would make selection worse, not weaker.
   Needs >= 4 independent reference poses per ligand, i.e. replicate jobs, deduplicated.
5. `build_submission.py build`, then `validate --expect-n <size>`.

Budget note: Boltz is about $0.156 per job at 20 samples. The primary Modal workspace is
over its spend limit; `xx-its-amit-xx` has ~$7.61 left, and the ~$140 hackathon workspace
is not configured on this box.

---

## Current state, 2026-09-14 — what the next tick should do

**The selector changed.** Cross-engine agreement (FINDING 011/012) replaces the FINDING 003
rule: **+0.0381** on CYP3A4 against its +0.0265, and it **generalises** - +0.1098 across 27
held-out P450 PROTEINS (44 construct sequences), positive on 13 of 16. No fitted parameters. A validated submission
already exists at `submissions/01_xeng_val87b.zip` (87 PDBs, mean LDDT-PLI 0.6164 against
PXR's winning 0.5640).

**The mechanism matters more than the number.** It is a *catastrophe detector*: the gain
tracks how bad the pool's worst poses are, not which protein it is. CYP3A4 wide-spread
ligands +0.0704, narrow-spread +0.0063 (null). Expect the **+0.04 regime** for a
well-behaved release. Which pool we submit from now matters as much as the selector.

### Venue facts, all measured

- **Six** OpenProtein engines fold protein+HEM+ligand: protenix, protenix_v2, esmfold2,
  rosettafold_3, boltz2, boltz_1x. Four need `--single-sequence`; an UPLOADED msa makes
  them fail server-side. Only alphafold2 genuinely cannot (it discards ligands).
- **`--samples` never diversifies the ligand here, for any engine.** Only `--replicates`.
  And replicates are not automatically distinct - dedupe (39% were duplicates).
- **protenix-v1 is DETERMINISTIC**: sd 0.0000 across replicates on 6/6 standard ligands
  and 14/14 organometallics. One replicate is all it will ever give, and that single pose
  is worth +0.007 in the reference set. Never buy more than one. rosettafold_3 in
  single-sequence mode is deterministic for about half of ligands; esmfold2 is
  properly diverse.
- **protenix_v2 is ALSO DETERMINISTIC** (FINDING 015, 2026-09-15). This list used to say
  it was "properly diverse". A 12 -> 24 replicate doubling moved the oracle on **0 of 489
  pairs**, delta +0.000000: replicates 12-23 are one file written twelve times, and the
  older half holds 3-4 distinct files out of 12. Nominal depth 12 is a REAL depth of ~4.
  Replicate count is not a lever for this engine; `num_recycles` / `num_steps` are the
  untested ones. Dedupe the POOL before believing any depth number, not just the
  reference - checking only the reference is how this went unnoticed for a whole campaign.
- **More engines is NOT better.** Best reference set is the two Protenix checkpoints
  (+0.0380); adding esmfold2 at matched depth drops it to +0.0178.
- Single-sequence mode unblocks all 185 P450 targets at once against ~2 days of serial
  MSA queue. Quality cost not yet measured - 30 probe jobs are running for that.
- Protenix parses the Ir/Ru organometallics Boltz cannot; esmfold2 does not (skip list at
  `p450_universe/esmfold2_unsupported.json`). rosettafold_3 parses them but places them
  badly in single-sequence mode (n=2, LDDT-PLI 0.002/0.034 - watch, do not yet conclude).

### Do next, in order

1. ~~Measure what single-sequence costs~~ **MEASURED - single-sequence is NOT usable for
   the P450 set.** On the 4 pairs with both, MSA gives pool mean 0.5825 / oracle 0.8760
   with **zero** catastrophic poses; single-sequence gives 0.0505 / 0.0505 with **75%**
   catastrophic. MSA wins 4/4, mean oracle delta **-0.8255**. Wait for the MSA queue;
   do not take the shortcut.

   This is exactly the trap that was flagged in advance: a degraded pool has more
   catastrophes, which INFLATES the cross-engine gain while producing far worse
   submissions. Judged on absolute oracle, as intended, it is disqualifying.

   Note single-sequence DOES work for one well-studied target - the CYP3A4 probe placed
   the ligand at Fe 2.23 A. It is diverse and unfamiliar P450s where removing the MSA
   destroys the prediction. Do not generalise from the CYP3A4 result.
2. **Re-run the P450 generalisation** as targets accumulate (`score_p450_pool.py score`
   then the cross-engine block). It has held at 17 targets; 185 is the goal.
3. **Keep P450 depth building** - `python scripts/ops/topup_p450.py --submit`. One
   command, idempotent, and it reports the gaps before queueing. The gap that matters is
   not visible from job counts: a pair with fewer than 3 POOL poses is dropped from the
   generalisation test entirely, so thin pairs cost whole PROTEINS - 48 of them were
   gating 10 proteins, and filling them moved the count 22 -> 27 -> 30. Run with no flags
   to see the gaps without submitting. It skips rosettafold_3 and protenix-v1 on purpose:
   both are deterministic enough that replicates produce byte-identical files.
3b. **Disk — DONE, 2026-09-14.** All four op1 pools are now on OneDrive
   (`onedrive:rclone-offload/cyp-structure/pool/op1`, 9,234 objects / 2.92 GiB), which
   took D: from 6.2 GB back to 9.6 GB. The selector is unaffected and that was verified
   against the genuinely archived state, not a simulated one: it reproduces the feature
   across all 1,740 poses to 2.8e-08 from the frozen reference. Pull any pool back with
   `cypstruct.storage.pull("pool/op1/<engine>", <dest>)` if a reference set ever needs
   re-deriving with different dedup settings.

   Original note kept for the method: the reference set is frozen at
   `data/processed/reference_set_cyp3a4.npz`
   (241 KB, verified to reproduce xeng to 2.8e-08). So the 2.2 GB `op1/protenix_v2` and
   589 MB `op1/protenix` pools ARE archivable - `build_xeng_feature` falls back to the
   frozen copy automatically. Archive order when D: approaches ~3 GB: `rosettafold_3` and
   `esmfold2` op1 first (both measured NOT to help the CYP3A4 reference set), then the
   Protenix pools. Push with `cypstruct.storage.push`, never through the `O:` drive letter.
4. **Rebuild the submission** whenever the pool changes:
   `build_xeng_feature.py` then `build_submission.py build --pool-dir <dir>`.

### Do NOT

- Add single-feature selection candidates without a mechanism. ~30 tested; everything
  anchor-local or population-level has failed.
- Trust `rho` as a proxy for selection value. It moved OPPOSITE to the gain twice.
- Trust a job status or a passing check without exercising the thing it names. Three
  engines and one unbuildable submission path were lost to exactly that this session.

