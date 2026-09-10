# Runbook — what to launch next, and what each step has to prove

Ordered so that each step's gate decides whether the next one is worth its runway.
`scripts/ops/watchdog.py` runs every two hours; the ops tick reads this file to decide
what to launch when nothing is in flight.

**Deadline: interim leaderboard 2026-09-24.** Budget backwards from that.

---

## Step 0 — the pool that scores itself  ✅ infrastructure ready

We have **99 CYP3A4 ligands with deposited crystal poses**
(`data/processed/validation_ligands.csv`, 83 iron-coordinated + 16 active-site). Every
pose we generate for these can be scored immediately with `cypstruct.pose.lddt_pli` and
`bisy_rmsd`. No waiting for the challenge release to start measuring.

**Always dry-parse on CPU first.** This is a gate, not a formality:

```bash
python scripts/cofold/preflight_parse.py --csv data/processed/validation_ligands.csv     --tag val87 --samples 20            # exits non-zero if anything fails to parse
python scripts/cofold/modal_boltz.py submit --csv data/processed/validation_ligands.csv     --tag val87 --samples 20 --seeds 1
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
python scripts/cofold/modal_af3.py submit --csv data/processed/validation_ligands.csv --tag val99
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
- **Poll for ADVANCING progress**, not liveness. A run that is alive and not advancing is a
  runaway; kill it and record why.
- Every launch goes through `budget.preflight_hours` and is recorded in the ledger
  **before** it starts, so the watchdog can clean up even if the launcher dies.
