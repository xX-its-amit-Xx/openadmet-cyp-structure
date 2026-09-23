# PRE-REGISTRATION — does reference depth beyond 4 buy SELECTION on CYP3A4?

Written and committed **before any job was submitted**. Nothing below may be edited
afterwards; deviations go in the finding, labelled as deviations.

---

## 1. The conflict this settles

The shipped selector (`cypstruct.xengine.select`, FINDING 011) scores a Boltz pose by its
mean Chamfer distance, in the heme frame, to independent-engine reference poses of the
**same** ligand. `build_xeng_feature` and `openprotein_cofold refset` both **refuse** below
**4** reference poses, because at one pose the feature measured **−0.0055**.

Two findings now give opposite release-day advice about going deeper than 4:

* **FINDING 016** ("reference depth saturates at 4"): within-pair on **487 P450 pairs**,
  4 sampler settings gave **+0.0315** and 8 gave **+0.0345**. Delta **+0.0030**, under a
  +0.0138 floor. Verdict: *"Buy exactly 4 settings. Eight costs twice the jobs for a
  difference indistinguishable from zero."*
* **FINDING 040 §7**: *"Spend the wall clock on reference depth instead ... Depth in the
  pool is worth +0.0027; depth in the reference is what makes the +0.0395 gain exist at
  all."*

**FINDING 038** measured that depth 8 is cheaply *achievable* on CYP3A4 (9 of 9 ligands,
16 jobs, 25.7 min) and states explicitly that it **does not claim** 8 reference poses
select better than 4 — it measured **supply, not value**.

CYP3A4 is a known family outlier (FINDING 022: family 0.86 / 89% sub-2 Å, CYP3A4 alone
0.555 / 33%), so FINDING 016's P450-universe answer is not automatically CYP3A4's answer.
This experiment measures the **value** of reference depth on **CYP3A4**.

---

## 2. The ligand rule — fixed here, applied without reading any score

Source: `data/processed/refdepth_chemistry_all87.csv` (FINDING 038). Every column used is
computable without a crystal structure: `n_heavy` / `n_rot` from SMILES, `pred_mode` from
the median `fe_donor_dist` over the ligand's own **predicted** pool poses cut at 2.6 Å
(FINDING 033), `tercile` from `n_heavy`.

1. Take **all 13 `type_I`** ligands. Type I is the scarce stratum and the one FINDING 033
   names as the real exposure; it is taken whole rather than sampled.
2. From the `type_II` ligands take **8 per heavy-atom tercile** (`small` / `medium` /
   `large`), drawn with `numpy.random.default_rng(41).choice(sorted(ids), 8,
   replace=False)`.
3. Drop the single `peripheral` ligand.

**n = 37** (13 Type I + 24 Type II; heavy atoms 8–51, rotatable bonds 0–17).
Written to `data/processed/refvalue_ligands.csv` **before** submission. No LDDT-PLI,
oracle, gain, `xeng` or selection column is read anywhere in the rule.

---

## 3. What is generated, and why ten waves and not eight

**Reference:** `submit --sweep 3x200,10x200,3x50,5x200,3x150` on **both** Protenix
checkpoints (`protenix_v2`, `protenix`) under tag `refvalue`, `--batch 5`, `--samples 1`.
That is **5 settings × 2 engines = 10 waves**, and by FINDING 038's mechanism
(`depth = engines × settings`, minus one per failed job) the target depth is **10**.

Three reasons the sweep is 5 settings rather than 4:

* **The top rung must not be a deterministic maximum.** FINDING 035 §"methods
  corrections": *a subsample curve's last rung is a deterministic maximum, so its final
  slope over-prices new samples.* If available depth were exactly 8, there would be
  exactly one subset of size 8 and the depth-8 point would be a max where depths 4–7 are
  averages — a bias in favour of the very conclusion being tested. At available depth 10
  every rung from 4 to 8 is a genuine average over `C(10, d) ≥ 45` subsets.
* **Failure margin.** FINDING 038 measured a 6.25% per-job failure rate. At depth 10 a
  ligand can lose two waves and still reach 8.
* `3x400` is dropped (slowest setting on both engines, the only one that failed in 038);
  `5x200` and `3x150` are taken from FINDING 016's second four. `num_recycles < 2` is
  refused by `parse_sweep` and is not attempted.

**Job count:** `10 × ceil(37/5) = 80` jobs. `data/processed/openprotein/jobs.json` holds
**610** against the 2,000/month cap; this takes it to **690 (34.5%)**.
`cypstruct.budget.preflight('openprotein_cofold', 80, venue='openprotein',
cap_key='openprotein_jobs')` must return `ok=True`. **If it refuses, the run stops and the
refusal is reported.**

**Pool: held completely fixed.** The existing 20-pose **unsteered Boltz-2** pool
(`D:/cyp_scratch/val87b_unsteered`, scored in `data/processed/poses_scored_val87b.csv`).
No pose is added, removed or re-scored. The **only** thing that varies across the primary
endpoint is which reference poses the feature sees.

---

## 4. The endpoint, the subsampling scheme and the draw count

For a ligand with `D` available reference poses and depth `d ∈ {4, 5, 6, 7, 8}`:

* the 20 × D Chamfer matrix `C[p, r]` is computed once;
* `xeng(p | S) = mean_{r ∈ S} C[p, r]` for a subset `S`, `|S| = d`;
* the pick is `argmin_p xeng(p | S)` (first index on a tie, which is what
  `select()`'s `idxmax` does);
* the ligand's score at depth `d` is the **exact mean over all `C(D, d)` subsets** of the
  picked pose's LDDT-PLI. `C(10, d) ≤ 252`, so **every subset is enumerated; there is no
  Monte-Carlo error in the depth curve itself.**

**Spread across draws.** A "draw" is one independent subset per ligand. The pooled mean's
exact variance is `Σ_l Var_l / n²` with `Var_l` the exact per-ligand across-subset
variance; the full distribution (p5 / p50 / p95, min, max) is sampled at **200,000 draws**
per depth with `default_rng(20260923)`. Both are reported. The spread *is* part of the
answer to "how much does which-4-you-pick matter".

**Reported at every depth, beside every pooled mean:** pool **oracle** (mean per-ligand
max LDDT-PLI over the 20 poses) and **random** (mean per-ligand pool mean) first; then
selected, gain over random, mean within-ligand Spearman ρ(`xeng`, LDDT-PLI), and the
fraction of ligands with the correct sign (ρ < 0).

**Noise floor, recomputed for THIS population**, using the FINDING 039/040 construction
(`noise_floor` in `depth_turnover_analysis.py`): a random feature's argmin is a uniform
draw, gain = `mean(drawn) − mean(pool means)`, at **2,000,000 draws** (≥ the 200,000
required), with the SE of the p95 reported. No floor is quoted from another finding.

---

## 5. The acceptance rule — fixed a priori

Primary quantity: **Δ = selected(depth 8) − selected(depth 4)**, exact, same 37 ligands,
same 740 poses.

| clause | requires |
|---|---|
| **A0 — gate** | the new sweep reference at full depth clears this population's floor; if it does not, the instrument is broken and no depth claim is made |
| **T1 — magnitude** | `Δ ≥ +0.0125` — FINDING 004's per-doubling benchmark, and 4 → 8 is exactly one doubling |
| **T2 — inference** | paired bootstrap CI over ligands (10,000 resamples) excludes 0 |
| **T3 — not one ligand** | the sign of Δ survives leave-one-ligand-out on all 37 |
| **T4 — monotone** | selected(d) is non-decreasing over d = 4…8 |

**Verdict map:**

* **T1 ∧ T2 ∧ T3 → "reference depth beyond 4 PAYS on CYP3A4."** FINDING 040's advice
  survives; FINDING 016 is family-specific.
* **0 < Δ < +0.0125 → "detectable, not worth the wall clock."**
* **|Δ| < 0.0030 → FLAT. FINDING 016 replicates on CYP3A4**, the depth-4 refusal rule is
  the right number and not a lower bound, and the release-day spare hour does **not** go
  to reference depth. This is the pre-declared conclusion for a flat curve, and it is
  reported as the result rather than as a failure to find one.
* **Δ < −0.0030 → deeper is worse**, which would be the sixth "more diversity, less
  selection" result in this repo (FINDINGs 013 / 016 / 027 / 034 / 035).

**Leverage, reported whichever way it lands:** drop-1 sensitivity on Δ (min and max over
the 37 leave-one-outs), and **how many ligands actually change pick** between depth 4 and
depth 8 (measured as the modal pick over subsets, and as the fraction of (d=4, d=8) subset
pairs that disagree).

---

## 6. Controls, all of which must pass before any number is believed

| # | control | must show |
|---|---|---|
| **C1** | numbering control (FINDING 021) | no ligand has LDDT-PLI identically 0.0 across its pool; the pool's per-ligand means match `poses_scored_val87b.csv` exactly |
| **C2** | reproduce the shipped `xeng` column | recomputing `xeng` for the 740 poses from `reference_set_cyp3a4.npz` matches `data/processed/xeng_val87b.csv` to < 1e-6 |
| **C3** | `argmin(xeng)` **is** `xengine.select()` | the per-ligand pick from `argmin` is identical to `select()`'s on all 37 ligands, on the shipped column |
| **C4** | `refset` filter counts | all eight counts reported **including the zeros**, with the FINDING 038 positive controls (md5 on a byte copy, coords on an identical vector, coords on +1e-4 Å, coords keeps +5.0 Å, same-wave on duplicate samples) |
| **C5** | distinctness margin | min / median pairwise Chamfer inside each ligand's reference set against the 0.05 Å dedupe tolerance |
| **C6** | waves, not jobs | depth is counted per `(engine, wave)`; the job count and the distinct-pose count are reported as two different numbers and never conflated |

**The answer-recognition gate (R1) is NOT run.** FINDING 036 established it applies to
*priors*, not to within-ligand comparators, and the shipped selector already fails it by
construction (crystal at the 34th percentile) while working.

---

## 7. What is *not* claimed

This measures reference depth 4 → 8 for **one** reference recipe (the two Protenix
checkpoints under a sampler sweep), on **one** pool (20 unsteered Boltz-2 samples), on
**37** CYP3A4 validation ligands with known crystals. It says nothing about depth beyond
10, about a third engine, or about a pool from a different generator. Depths below 4 are
not tested: the refusal there is FINDING 011's and is not under review.
