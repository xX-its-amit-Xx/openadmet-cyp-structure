# PRE-REGISTRATION — pool DEPTH on the predicted-Type-I stratum

**Written:** 2026-09-22, before a single new job was submitted and before any new pose was
scored. **Committed before launch.** Nothing in this document may be edited afterwards;
corrections go in `docs/FINDING_034_type_i_depth.md` and are labelled as deviations.

---

## 1. Why

`FINDING_033` measured that the shipped selector `-z(xeng)` is **not** biased by binding
mode where the data can say so (headroom-matched Type II − Type I difference **−0.0030**,
CI [−0.0132, +0.0075], 342 P450 pairs), and that the real exposure is the **pool**: on
CYP3A4 the predicted-Type-I stratum has a random baseline **0.134** below Type II and an
**oracle 0.058 below**, so a Type I-rich blind test set scores about **0.51** instead of
**0.64**, with more than half of that gap unreachable by any selector. Its recommendation
(§"What to do at submission time", item 3) was to spend the remaining budget on **pool
depth for the predicted-Type-I ligands**.

That recommendation is an assumption until it is measured, and this repo has a hard-won
distinction that makes the assumption non-trivial:

* **Pool EXPANSION** — adding a *different kind* of pose — has raised the oracle and left
  selection flat or worse **four times**: a second engine (FINDING 013, +0.0375 of
  unreachable oracle), a sampler sweep (FINDING 016, −0.0038), 1.78 M rigid ligand
  rotations (FINDING 027, oracle +0.0108 / selection **−0.0145**), and a second copy of
  the query ligand (FINDING 031, oracle +0.0237 / selection **−0.0027**).
* **Pool DEPTH** — more samples of the *same kind* — is the one lever that has paid, at
  **+0.0125 of selected score per doubling** (FINDING 004), with the oracle climbing at
  +0.024 per doubling and no sign of saturation at 20.

**This experiment is DEPTH, not EXPANSION, and the distinction is itself under test.**
If the oracle climbs and selection does not, that is expansion behaviour wearing depth's
clothes, and that is the finding.

> **The question.** Does adding poses to the predicted-Type-I stratum convert into
> SELECTED LDDT-PLI, and at what rate per doubling, against the +0.0125 benchmark
> established for the pool at large?

---

## 2. The ligand list, and exactly how it was chosen

**Rule, fixed a priori, prediction-side only.** A ligand is **predicted Type I** if the
**median `fe_donor_dist` over its own 20 existing pool poses is > 2.6 Å**. `COORD_MAX =
2.6 Å` is the repo's existing constant from `scripts/structure/build_reference_set.py`,
unchanged, and is the rule `FINDING_033` validated at **96.6% accuracy** against the
crystal label (84 of 87). Nothing from any crystal enters the selection of this list.

This yields **n = 14 ligands**, computed from
`data/processed/binding_mode_labels_cyp3a4.csv` before this document was written:

| ligand | pdb | crystal mode | predicted median Fe–donor (Å) | frac coordinating | pool mean LDDT-PLI |
|---|---|---|---|---|---|
| MWS | 6OOB | type_I | 3.848 | 0.00 | 0.668 |
| CFF | 8SO1 | type_I | 3.952 | 0.00 | 0.481 |
| MF8 | 5G5J | type_I | 3.965 | 0.00 | 0.649 |
| PG0 | 9BV6 | type_I | 4.499 | 0.00 | 0.413 |
| PG4 | 9YK4 | type_I | 4.870 | 0.00 | 0.501 |
| ERY | 2J0D | type_I | 5.012 | 0.00 | 0.611 |
| 08Y | 3UA1 | type_I | 5.250 | 0.00 | 0.704 |
| **D0R** | 3TJS | **type_II** | 5.430 | 0.25 | 0.325 |
| 08J | 5TE8 | type_I | 5.820 | 0.00 | 0.282 |
| TCI | 9PLJ | type_I | 5.884 | 0.00 | 0.365 |
| MWY | 6OOA | type_I | 6.077 | 0.00 | 0.372 |
| A1CIW | 9PLK | type_I | 6.144 | 0.00 | 0.378 |
| YNV | 7LXL | type_I | 6.491 | 0.00 | 0.289 |
| MWV | 6OO9 | type_I (predicted **peripheral**, 7.70 Å) | 7.695 | 0.00 | 0.398 |

**Where the two labels disagree — reported, not hidden.** Over all 87 ligands the
prediction-side and crystal labels disagree on exactly three, and two of them touch this
list:

| ligand | pdb | crystal Fe | predicted median | what it means |
|---|---|---|---|---|
| **D0R** | 3TJS | 2.107 Å (**II**) | 5.430 Å | model fails to coordinate a genuine Type II → **included** here, because the prediction-side rule is what we would have at submission time |
| **QDY** | 6UNJ | 2.824 Å (**I**) | 2.231 Å | model **forces** coordination on a non-coordinator, all 20 poses → **excluded** here, for the same reason |
| MWV | 6OO9 | 3.763 Å (I) | 7.695 Å | right mode, ejected past the active site → included (median > 2.6 Å) |

**Deliberately NOT excluded:** PG0 and PG4 are PEG fragments and `FINDING_031` excluded
that class from a different experiment. They are kept here because `FINDING_033`'s Type I
stratum contains them and this experiment must deepen *that* stratum. A sensitivity
analysis dropping them (n = 12) is pre-registered as a secondary, not a substitute.

---

## 3. Venue, engine, and what "depth" is forced to mean

**Venue: OpenProtein only.** Unmetered for this account; `cypstruct.budget.preflight` is
called with `cap_key="openprotein_jobs"` (cap 2,000/month) and **if it refuses, this
experiment stops and reports the refusal**. Modal is over its spend cap and will not be
touched. Explorer is not reachable inside the interim deadline (MSA staging, CLAUDE.md).

**Engine: `boltz2`, in single-sequence mode.** The existing pool is Boltz-2, so this is
the same architecture — the closest thing to "more of the same kind" that the permitted
venue can supply. Measured facts fixing the configuration:

* `diffusion_samples` does **not** diversify the ligand on OpenProtein, for any engine
  (FINDING 009). So `diffusion_samples=1`, and depth comes from **replicate jobs**.
* An **uploaded MSA makes `boltz2` fail server-side.** Verified on the stored probe jobs
  before writing this: `boltz_2__msa` → `JobStatus.FAILURE, "internal server error"`,
  while the single-sequence `boltz2` probe → `JobStatus.SUCCESS`. This is not a choice.
* `num_recycles=3`, matching every other OpenProtein run in this repo.

### Confounds, listed before any number exists

| # | confound | size, and what is known about it |
|---|---|---|
| **C1** | **No MSA.** The existing pool used a 6,979-sequence CYP3A4 alignment; the new arm cannot use one. | The largest confound, and unavoidable at this venue. It is expected to make the new poses *worse*, which biases the depth test **against** the hypothesis, not toward it. |
| **C2** | **No explicit Cys442-SG→heme-FE bond.** The existing pool has one. | `FINDING_019` measured this as a **null**: Δ = −0.0069 LDDT-PLI, p = 0.43 over 84 paired P450 pairs, and unbonded predictions already reproduce the crystal Fe-donor distribution. Weak confound. |
| **C3** | **Different serving stack / checkpoint.** OpenProtein-hosted Boltz-2 vs self-hosted on Modal. | Unmeasurable from here. Reported. |
| **C4** | **Sampling mechanism differs.** The existing 20 poses are 20 diffusion samples of one Modal job (verified distinct: 20 unique `xeng` values on every one of 87 ligands, median within-ligand LDDT-PLI sd 0.062). The new poses are one sample per replicate job. | Both produce genuinely distinct ligand poses; the *mechanism* of diversity differs. |
| **C5** | **Pool heterogeneity.** Mixing two configurations changes the within-ligand distribution that `-z(xeng)` standardises over. | This is precisely the expansion-vs-depth question, and is why the pure-subsample curve (§4.1) is reported as the primary internally-consistent measurement. |
| **C6** | **n = 14.** `FINDING_007`'s noise floor is n-dependent; at n = 14 a *random* feature's selection gain reaches **+0.0433** at the 95th percentile and +0.0622 at the 99th, against the pooled **+0.0144**. | Every threshold below is set against the floor recomputed **inside this stratum**, never against the pooled one. |
| **C7** | The four exact-zero LDDT-PLI rows in the existing pool (PG0 ×1, PG4 ×3) are genuine large-RMSD ejections, confirmed by `FINDING_033` C2 (BiSyRMSD 23–27 Å, `mapped` True, other poses of the same ligand scoring 0.44/0.58). They are **kept**. | Not a numbering failure. |

**If the measured distinct-pose yield of this configuration makes depth unbuyable, that
is a reportable result and the experiment terminates at §4.1 + the yield number.**
`FINDING_031` measured protenix_v2 as *bimodal* in its configuration — exactly two
distinct poses from four replicate jobs — and `FINDING_015` measured it as fully
deterministic in another. **Determinism is a property of the configuration, not of the
engine**, so it is measured here before depth is bought.

### The replicate-count rule, fixed now

1. **Yield probe first.** 2 ligands × 6 replicate jobs, `diffusion_samples=1`. Count
   **distinct poses** by `cypstruct.xengine._dedupe` (per-atom max deviation ≥ 0.05 Å in
   the ligand's own heme frame). Let `y` = mean distinct poses per replicate.
2. **Buy** `R = ceil(20 / max(y, 0.1))` replicates per ligand, **capped at R = 30** and at
   **500 total complexes**. If `y < 0.25` (i.e. fewer than ~5 distinct poses reachable in
   20 replicates) the configuration is declared **depth-incapable** and no bulk run is
   made.
3. `R` is chosen from the **yield**, which is a property of the configuration, never from
   any LDDT-PLI, oracle or selection number.

---

## 4. The measurements

### 4.1 The subsample depth curve — primary, free, internally consistent

Subsample the existing 20-pose pool. **Rungs fixed now: 1, 2, 3, 5, 8, 10, 14, 20.**
(1, 2, 3, 5, 8, 10, 20 are `FINDING_004`'s ladder; 14 is added to fill the top octave.)

For each rung `k` and each of **≥ 64 draws** (we will use 256), draw `k` of the 20 poses
per ligand without replacement, then report, per stratum:

* **oracle** = mean over ligands of max LDDT-PLI in the drawn subset — **reported first, always**;
* **selected** = mean over ligands of the LDDT-PLI of the pose chosen by
  `cypstruct.xengine.select()` on the drawn subset, with **random tie-breaking**;
* **random** = mean over ligands of the mean LDDT-PLI of the drawn subset (exact
  expectation), which must be flat in `k` — a failure of that flatness invalidates the
  subsampling and the run is void.

Strata: **predicted Type I (n = 14)**, **predicted Type II (n = 73)**, **all (n = 87)**.
Per-ligand tables are written out alongside the pooled numbers.

**Rate per doubling** = OLS slope of the stratum's mean vs `log2(k)` over all eight rungs,
fitted separately for oracle and for selection, with a bootstrap CI over ligands
(10,000 resamples). The benchmark is `FINDING_004`'s **+0.0125 per doubling** for
selection and **+0.024** for the oracle.

### 4.2 The augmented depth curve — the new poses

With `N` new distinct Boltz-2 poses per Type I ligand, repeat §4.1 on the **union** pool
at rungs **20 + {0, ⌈N/4⌉, ⌈N/2⌉, N}**, always keeping all 20 existing poses and drawing
the new ones, ≥ 64 draws. Also report the **new-poses-only** pool at matched depth, so
that the quality of the new arm is visible separately from its contribution.

### 4.3 Projection onto a Type I-rich test set

For a test set that is fraction `f` Type I, expected selected score is
`f · S_I + (1 − f) · S_II`. Report at **f = 0.5** (and f = 0.0, 0.25, 1.0 for context),
before and after the depth spend, with a 95% bootstrap CI propagated from the per-ligand
gains of both strata. `S_II` is unchanged by this experiment by construction.

### 4.4 Controls that must pass before any score is believed

* **N1 — numbering (FINDING 021).** Every new pose is renumbered onto the crystal by
  residue-NAME agreement; the offset and the identity fraction are recorded per pose. Any
  pose with identity < 0.95 is reported and excluded with a count.
* **N2 — the shipped column reproduces.** `select()` over the existing unsteered pool must
  return selected **0.6164**, oracle **0.6975**, random **0.5769**, gain **+0.0395**,
  within-ligand ρ **−0.2582**, 75.9% correct sign. *(Verified before writing this
  document: all six match to four decimals.)*
* **N3 — the within-stratum noise floor**, recomputed here: a random feature substituted
  for `xeng`, **4,000 draws inside the 14-ligand stratum**, at every depth rung.
* **N4 — distinct poses, not jobs.** Dedupe every new pose at tol 0.05 Å; report
  distinct/collected. Every depth number is quoted on the **deduplicated** pool.
* **N5 — every filter reports its count**, including the ones that fire zero times.

---

## 5. The acceptance rule

Primary endpoint **Δ_sel** = selected LDDT-PLI on the 14 predicted-Type-I ligands at the
deepest achieved rung minus selected at depth 20, both from §4.2, paired per ligand.
Secondary **Δ_oracle**, same form.

| verdict | fires when |
|---|---|
| **SHIPS** | **Δ_sel ≥ +0.020**, **and** its 95% paired bootstrap CI over the 14 ligands excludes 0, **and** the §4.1 Type I *selection* rate per doubling is ≥ **+0.0125** |
| **MEASURED — expansion in depth's clothing** | **Δ_oracle ≥ +0.020** while Δ_sel < +0.020 or its CI includes 0 |
| **REFUTED** | **Δ_sel ≤ 0 and the CI upper bound < +0.020** |
| **MEASURED — depth unbuyable here** | the §3 yield rule declares the configuration depth-incapable (`y < 0.25`), or fewer than 5 distinct new poses per ligand are obtained |

If none of the four fires, the verdict is **MEASURED** and the curve is the deliverable.

**Explicitly pre-committed:** the rungs in §4.1 and §4.2, the 2.6 Å threshold, the ligand
list in §2, the 256/4,000/10,000 draw counts, and the +0.020 bar are all fixed by this
document. **No rung, threshold or subset may be chosen after seeing a scored result.**

---

## 6. Deliverables

`docs/FINDING_034_type_i_depth.md`; scripts under `scripts/cofold/` and
`scripts/structure/`; results under `data/processed/type_i_depth_*`. mmCIF to `C:\Temp`,
never to `D:`. Job count actually used is reported in the finding.
