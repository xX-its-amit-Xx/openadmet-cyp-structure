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

**Superseded in part by `docs/FINDING_001_selection_is_the_whole_problem.md` (n=87,
3,360 poses). Read that first — two claims below were falsified by it, and are kept here
with corrections attached because the errors are instructive.**

**A. Co-folding fails on CYP3A4 by getting the ligand's ORIENTATION wrong, not its
location.** OpenADMET's own analysis reports the heme is placed correctly and the
backbone is accurate, while only 57.8% of Boltz-2 poses reach BiSyRMSD < 2 Å. Our own
run agrees on the mechanism: whole-protein CA fit is ~1.0 Å and the errors are in the
ligand. **Still stands.**

**B. Most CYP3A4 ligands coordinate the iron.** Parsing 116 deposited entries: of 101
unique ligands, **83 coordinate the heme iron directly**, through nitrogen in 103 of 104
chain observations. Fe–donor 1.94 / **2.20** / 2.38 Å (p5/p50/p95), S(Cys442)–Fe–donor
159.5 / **171.2** / 177.6°, and only **0.7%** of ligand atoms on the proximal face.
**Still stands** — and it is the calibration the scorer uses.

**C. ⚠️ FALSIFIED: "co-folders never reach coordination geometry."** Inherited from a
sibling-repo pilot that put zero of 120 poses inside the window, closest 2.80 Å. With an
explicit `bond` from Cys442 SG to the heme FE and a 6,979-sequence MSA, Boltz-2 reaches a
**median Fe–donor distance of 2.23 Å with 84% of poses inside 1.90–2.45 Å** — essentially
the crystallographic distribution. The pilot's failure was its setup, not the model.

**D. ⚠️ FALSIFIED: "naming the donor atom converts an orientation failure into a
constraint."** Measured across 87 ligands, forcing the Fe→donor contact is a **null,
slightly negative**: oracle −0.0023, selected −0.0278, sub-2 Å ligand rate −4.7 points.
It was redundant, because the heme bond alone already produces the right geometry.
**Drop the steered arm.**

**E. THE ACTUAL FINDING: Boltz's confidence does not rank poses.** Within-ligand Spearman
between `complex_ipde` and true LDDT-PLI is **−0.092 (steered) and −0.033 (unsteered)**,
positive for about half the ligands, i.e. chance. **Selecting the highest-confidence pose
is worse than selecting at random** (0.5428 vs 0.5585; 0.5706 vs 0.5769). The cross-model
z-hybrid that won PXR ranks within an engine by exactly this signal, so it has no
foundation here.

**F. The prize, quantified.** Pool oracle **0.6975** against selection **0.5706** —
**0.127 LDDT-PLI per ligand** unclaimed in a pool already paid for. The PXR winning entry
scored 0.564, below our pool's current selection. **Generation is not the bottleneck.**

**G. The coordination term survives as a SCORER.** Poses that coordinate the iron score
**+0.14 LDDT-PLI** over those that do not (0.5992 vs 0.4588). Real CYP-specific signal —
used as a selection feature, not as a generation constraint.

**H. Confidence fails in BOTH engines, so this is the problem, not a Boltz quirk.**
Protenix-v2's seven confidence fields rank poses within a ligand at `frac(rho>0) = 0.50`
for four of five, every p > 0.6 (FINDING 009). Two independently trained architectures,
both at chance ranking their own samples.

**I. The iron anchor is saturated; the signal is in substituent placement.** 84% of Boltz
poses already sit inside any reasonable coordination window, so recalibrating that window
on the whole P450 superfamily moved selection by -0.0002 and the angle term is inert to
four decimals (FINDING 008 addendum). Fourth independent confirmation.

**J. Architectural diversity has not bought decorrelation - 0 for 2.** Chai rho = +0.45,
Protenix rho = +0.474 against Boltz's per-ligand oracle. The hard ligands are hard for
everyone, which points the remaining upside at scoring rather than at another engine.

**Therefore:** the physics scorer in `docs/QM_SCORER_DESIGN.md` is now the whole project.
It does not have to beat a strong incumbent; it has to beat random, which the incumbent
fails to do - and which a second engine's incumbent also fails to do.

## The data situation changed on 2026-09-13

`n = 87` was the ceiling on every negative result: FINDING 007 puts the noise floor at
+0.0138, so a real +0.015 effect is unprovable there. The whole P450 superfamily is now
harvested (PF00067): **493 foldable pairs, 367 ligands, 185 distinct targets**, with the
active-site atoms cached in `data/processed/p450_universe/p450_atoms.parquet` so no
geometric question needs another 585 downloads. 185 targets turns
leave-one-ligand-cluster-out into **leave-one-TARGET-out**.

**Venue note.** Modal is over its cap and reserved for fine-tuning; pools now come from
**OpenProtein**, where `protenix`, `protenix_v2` and **`esmfold2`** run protein+HEM+ligand
(boltz-1/-1x/-2 and rosettafold-3 fail; alphafold2 discards ligand chains).
Read FINDING 009 before launching anything there: `diffusion_samples` does NOT sample the
ligand - it varies only the protein, and pose diversity requires `--replicates`.

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
