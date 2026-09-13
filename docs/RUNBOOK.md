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
5. `build_submission.py build`, then `validate --expect-n <size>`.

Budget note: Boltz is about $0.156 per job at 20 samples. The primary Modal workspace is
over its spend limit; `xx-its-amit-xx` has ~$7.61 left, and the ~$140 hackathon workspace
is not configured on this box.

---

## Current state, 2026-09-13 — what the next tick should do

Venue is **OpenProtein**, not Modal (over cap, reserved for fine-tuning). Only `protenix`
and `protenix_v2` run protein+HEM+ligand there. **Read FINDING 009 before launching**:
`diffusion_samples` does not sample the ligand, `--replicates` does.

### In flight

| campaign | state | command |
|---|---|---|
| CYP3A4 replicates, 6 of 12 | ~419 folds queued | `openprotein_cofold.py collect --engine protenix_v2 --tag op1` |
| 11 organometallic + 3 new CYP3A4 ligands | queued | same, `--tag recover` |
| P450 MSAs | 11 of 185 done, ~1 per 15 min | `p450_campaign.py msa-status` |
| P450 folds | 12 jobs, 120 pairs unblocked | `p450_campaign.py collect` |

### Do these in order, when the gate opens

1. ~~Resume replicate depth to 12 once the fold queue is under ~100 pending~~
   **CORRECTED and already done, 2026-09-13.** The premise was wrong. MSA throughput is
   serial and independent of fold load: 7 -> 9 in ~1 h while 490 folds were queued, then
   9 -> 12 in ~2 h with the queue draining - about one per 20-25 min either way.
   Throttling the folds bought nothing and cost six replicates of depth. Do not throttle
   folds for the MSA queue's sake again; the two do not compete. Command, if depth is
   needed beyond 12:
   `openprotein_cofold.py submit --engine protenix_v2 --samples 5 --batch 4
   --replicates 12 --tag op1` (resumes on `(rep, ligand)`; re-running is safe).
   It was stopped at 6 because ~490 queued folds were starving the MSA searches.
   Justified by FINDING 004: the oracle is still climbing and the selector tracks it at
   +0.0125 per doubling, which is larger than any feature gain measured so far.

1b. **BEFORE using esmfold2 or rosettafold_3 in the FINDING 011 reference set, verify
   their replicates actually differ.** rosettafold_3 runs in single-sequence mode, which
   removes the MSA as a source of stochasticity, so its replicates could be deterministic -
   the FINDING 009 trap in a new engine. One command:

       python - <<'PY'
       # per-atom sd across replicates must be >> 0; FINDING 009 was sd = 0.0000
       PY

   Expect ~2-4 A per-atom sd, as Chai and Protenix-between-replicates show. If it is
   0.0000, replicates of that engine are worthless as independent opinions and only ONE
   pose per ligand can enter the reference set.

2. **Re-run the FINDING 011 test** once most ligands have >= 4 replicates:
   `python scripts/structure/cross_engine_agreement.py`. This is a pre-registered test
   with a stated prediction, not a fishing expedition. Cross-engine agreement gave
   rho = -0.202 at p = 0.00053 using a Protenix side of ONE pose per ligand. If the
   signal scales with independent poses, the combination may stop being subtractive; if
   it does not, retire the feature.

3. **Fold the rest of the P450 set** as MSAs land:
   `p450_campaign.py submit --samples 3 --batch 4 --replicates 3` (raise `--limit` as the
   queue allows), then `score_p450_pool.py score / features / validate` for the
   leave-one-TARGET-out result.

### Do NOT

- Add more single-feature selection candidates without a mechanism. Roughly 30 have been
  tested; everything anchor-local or population-level has failed (FINDINGS 002, 006, 008
  addendum, 010) and only within-ligand comparisons have ever worked. Six were tried on
  2026-09-13 alone, and best-of-six random scores about +0.012, which is most of why the
  best of them did not count.
- ~~Submit large fold batches while MSA searches are pending~~ — measured false, see
  item 1. Folds and MSA searches do not compete; the MSA queue is serial at ~1 per
  20-25 min regardless of what else is running.
