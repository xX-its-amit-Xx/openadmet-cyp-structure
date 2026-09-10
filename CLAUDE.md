# OpenADMET CYP3A4 — Structure Prediction Track

Entry for the structure-prediction track added to the OpenADMET CYP blind challenge,
announced 2026-09-08. First scored at the **interim leaderboard, 2026-09-24**; final
2026-11-03. Sibling repo `../OpenADMET-cyp-challenge` holds the activity track.

---

## What the challenge actually is

OpenADMET solved **CYP3A4 by cryoEM** (Orta, with the Fraser and Lander labs) and is
releasing ligand-bound structures for chemical matter already inside the blind challenge.
Predictions go through the same HuggingFace portal. Preprint:
`https://www.biorxiv.org/content/10.64898/2026.09.01.748687`.

From their announcement, three facts that shape everything here:

1. The construct is a **symmetric trimer** of the same truncated protein used for
   crystallography, solved unliganded and ligand-bound.
2. For several ligands the density fits **multiple mutually exclusive conformations**,
   not one pose.
3. The pocket remodelling is concentrated in the **F/G loop** — the region left
   unmodelled in many X-ray structures — and they name this as why **co-folding does
   poorly here**.

Structures are split across the train and test sets. Metric follows the PXR structure
track: **LDDT-PLI** primary, BiSyRMSD and lddt_lp secondary.

---

## The thesis this repo is built on

Three measurements, made here, that compose into one plan.

**A. Co-folding fails on CYP3A4 by getting the ligand's ORIENTATION wrong, not its
location.** OpenADMET's own co-folding analysis reports the heme is placed correctly and
the backbone is accurate (no Cα RMSD above 0.75 Å), while only **57.8% of Boltz-2 poses**
and **31.5% of OpenFold3 poses** reach BiSyRMSD < 2 Å on CYP3A4. Errors include
near-180° flips.

**B. Most CYP3A4 ligands coordinate the iron, and that fixes orientation.** Parsing 116
deposited CYP3A4 entries (`scripts/structure/build_reference_set.py`): of 101 unique
ligands, **83 coordinate the heme iron directly**, through nitrogen in 103 of 104 chain
observations. The measured geometry is tight:

| observable | p5 | p50 | p95 |
|---|---|---|---|
| Fe–donor distance (Å) | 1.94 | 2.20 | 2.38 |
| S(Cys442)–Fe–donor angle (°) | 159.5 | 171.2 | 177.6 |
| Fe out-of-plane (Å) | −0.42 | −0.05 | +0.06 |
| Fe–SG(Cys442) (Å) | 2.12 | 2.37 | 2.51 |

Only **0.7%** of ligand atoms across the whole set sit on the proximal heme face, making
the distal-side test a near-perfect validity filter.

**C. We can predict which atom coordinates.** `cypstruct.chem.coordinating_atoms` names a
competent donor for **82 of the 83** coordinated ligands (98.8% recall,
`scripts/structure/validate_donor_prediction.py`). Its weakness is honest and known: it
over-calls coordination for **5 of 15** type I substrates that carry an unused aromatic
nitrogen (33% false-positive rate).

**Therefore:** for the dominant binding mode, naming the donor atom converts an
orientation failure into a constraint we can hand the model — and into a check the
scorer can run. Every co-folding run submits a **steered and an unsteered arm** so this
is measured, not assumed.

And a prior pilot in the sibling repo cofolded 24 CYP ligands × 5 samples without
steering and put **zero of 120 poses** inside the Fe-coordination window (closest 2.80 Å).
That is the baseline the steered arm has to beat.

---

## Storage — read before writing anything

**Measured 2026-09-09: C: ~1 GB free, D: ~1.8 GB free, O: 4.9 TB free.** The sibling
repo's CLAUDE.md says "C: ~16 GB" and the global instructions say "~85 GB". Both are
stale by an order of magnitude.

**`O:\` is not a bulk data path.** It is an rclone mount with `--vfs-cache-mode full`,
`--cache-dir C:\Temp\rclone-cache`, and **no `--vfs-cache-max-size`**. Every byte read or
written through it is cached on C: indefinitely. Merely *listing* OneDrive during
orientation grew that cache to 8 GB and drove C: to zero, breaking unrelated tools.

- Bulk artifacts → `cypstruct.storage.push()` / `storage.Batch()`, which use
  `rclone copy` against the `onedrive:` **remote** and never touch the VFS cache.
- Reference mmCIFs stay local (`data/reference/`, a few hundred MB, re-read constantly).
- Heavy *hot* compute → Modal or Explorer. **This box has no GPU.**
- `cypstruct.paths.guard_scratch()` raises before a write rather than filling the disk.

**Worth fixing at the source:** add `--vfs-cache-max-size 8G --vfs-cache-max-age 24h` to
`C:\Temp\tools\mount-onedrive.vbs`. One line, needs a remount, and permanently removes
this hazard.

---

## Compute venues

- **Modal** — primary. Monthly credits. All remote functions declare `timeout`,
  `retries=1`, `max_containers`. `cypstruct.budget` refuses launches over the monthly cap
  **before** they start; `scripts/ops/watchdog.py` kills anything that runs away.
- **Explorer (SLURM)** — the `gyorilab` partition has 4× H200, 30-day walltime, no GPU cap.
  ⚠️ **GPU nodes have no direct internet** (proxy `http://10.99.0.130:3128`). This is
  structurally incompatible with `--use_msa_server`; MSAs must be precomputed on a login
  node and staged to `/scratch`. The PXR campaign lost ~12 days to this and never started
  training. Budget MSA staging as its own task.
- **Kaggle** — AF3 weights live here; also free GPU with a weekly quota.
- **Boltz hosted API** — authed, but accepts only templates, **not `constraints`**, which
  is why co-folding is self-hosted on Modal.
- **OpenProtein.ai**, **Colab Pro+**, **molab** — burst venues.

---

## Lessons inherited from the PXR structure campaign (2nd/50, 0.5640 LDDT-PLI)

1. **Validate the reference frame before trusting any cross-model metric.** Some engine
   exports landed ~20 Å off the crystal and inflated every agreement number about twofold
   before anyone noticed. `cypstruct.pose.align_by_residue` matches by residue number,
   never by array position, and `lddt_pli` is superposition-free by construction.
2. **Select with a cross-model z-hybrid, not with pLDDT.** Per engine take the best sample
   by its native confidence, z-score within that engine, then take the best engine per
   ligand. Sign matters: PDE and PAE are *errors*, so a raw `max()` picks the worst pose.
3. **A new engine joins as a TAIL RESCUER, not a pool member.** Protenix beat the pool on
   6 of 8 holo cases but *regressed* the board 0.5551 → 0.5241 as a full member; swapping
   it in for only the 8 least-confident ligands gained **+0.017** and won the campaign.
4. **Fine-tuning cost that campaign enormous runway and returned −0.0020.** Five platforms,
   twelve days, fifty training steps. It was deferred for lack of runway, not disproven —
   but with 15 days to the interim deadline it is not the first bet. See
   `docs/FINETUNE_PLAN.md` for the better-evidenced low-N recipe if we do spend on it.
5. **Pool oracle first.** Selection gains are meaningless without the best achievable score
   in the pool. Report it alongside every selection result.

---

## Layout

```
src/cypstruct/
  paths.py      paths + guard_scratch()
  storage.py    OneDrive push/pull that bypasses the C: VFS cache
  budget.py     spend ledger, cost model, preflight that REFUSES
  targets.py    sequences, heme, Cys442, pocket prior, F/G span, RCSB queries
  chem.py       standardise + coordinating_atoms() (98.8% recall)
  pose.py       structure loading, alignment, lddt_pli / bisy_rmsd / lddt_lp
  select.py     z-hybrid + tail rescue + conformer clustering
  qmscore/geometry.py   tier-0 physics terms, calibrated on the reference set
scripts/
  structure/build_reference_set.py        harvest + measure all PDB CYP3A4
  structure/validate_donor_prediction.py  honest check on the steering mechanism
  cofold/modal_boltz.py                   steered/unsteered Boltz-2 on Modal
  ops/watchdog.py                         kill runaways, watch spend and disk
docs/QM_SCORER_DESIGN.md                  the CYP-specific physics scorer
```

---

## Rules

- **Honest gating.** A term or feature ships only if it beats the incumbent z-hybrid on
  **leave-one-ligand-cluster-out** folds, not random CV, and its gain correlates with
  where the incumbent errs. Anything absorbed by the incumbent is a negative — log it.
- **No leaky features.** Anything derived from the reference structure, *including the
  pocket definition*, must be computable from the prediction at inference time.
- **Resumable.** Skip-done markers everywhere; poll for ADVANCING progress, never assume
  liveness means progress.
- **Measure the pool oracle before reporting any selection number.**
