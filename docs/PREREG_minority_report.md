# PRE-REGISTRATION — the minority report

**Written:** 2026-09-23, before any endpoint in this line exists. Committed before the
analysis script produces a single number. Not edited afterwards.

**Successor to:** `FINDING_036`, which closed tie-breaking by its own term oracle
(a **perfect** top-2 tie-break is worth **+0.0129** against a **+0.0134** floor; top-5
reaches 39.2% of the 0.0811 oracle gap) and named the open question:

> **60% of the gap sits outside the top 5.** The question is not *"which of these two"*
> but *"why is the right pose ranked 8th?"*

---

## 1. The hypothesis, and where it comes from

`FINDING_036` §3 established that the shipped selector fails the answer-recognition gate —
the query's own crystal sits at the **34th percentile** of its own 20 predictions,
binomial **p = 0.0012** — while still selecting at **+0.0395**. Its stated mechanism is
`FINDING_024`: the co-folders **share** CYP3A4's 30° orientation error, so the
cross-engine consensus is *displaced from the truth* but still informative about the
*ordering* of the predictions.

**The hypothesis under test.** If every engine makes the same mistake, then the best pose
in the pool is a **minority report** — an outlier that the consensus actively down-ranks —
and `xeng` is majority-voting for a shared error. Concretely: predicted poses fall into
orientation clusters in the heme frame; `xeng` ranks by cluster membership (large cluster →
low `xeng`); and the oracle pose sits disproportionately in a **small** cluster.

**What would make this worth building on.** A perfect *"the consensus is wrong here, take
the minority cluster instead"* rule would have to be worth more than a random feature is,
on this pool, at this n. `FINDING_036` is the standing warning that a **perfect** rule can
be worth **less than the noise floor** — which is why step 4 is an oracle and not a
selector, and why nothing is built before it clears.

---

## 2. Data — everything is on disk, zero new inference

| what | source | n |
|---|---|---|
| pool A | `D:/cyp_scratch/val87b_unsteered`, heme-frame cache `C:/Temp/cyp_tiebreak/frames_A.npz` | 87 ligands × 20 unsteered Boltz-2 poses = 1,740 |
| truth | `data/processed/poses_scored_val87b.csv`, `arm == "unsteered"` | 1,740 rows |
| shipped feature | `data/processed/xeng_val87b.csv` | 1,740 rows |
| pose–reference Chamfer matrix | `C:/Temp/cyp_tiebreak/dmat_A.npz` (rebuilt from `reference_set_cyp3a4.npz` if absent) | per ligand |
| prediction-side covariates | `data/processed/binding_mode_labels_cyp3a4.csv` | 87 |

Pool B (the 14-ligand, depth-40 `FINDING_035` stratum) is **out of scope as an endpoint**.
`FINDING_036` established that n = 14 cannot grade a rule (floor +0.0435; its best number
was one ligand). It may be quoted as a descriptive appendix and **nothing is concluded from
it**.

CPU only. Intermediates go to `C:/Temp/cyp_minority`, never to `D:`.

---

## 3. Controls that must pass before any endpoint is read

| id | control | pass condition |
|---|---|---|
| **C-XENG** | recompute `xeng` from `dmat_A` and join to the shipped column | max abs diff < 1e-6, 1,740/1,740 joined, 0 left-only, 0 right-only |
| **C-SEL** | `argmin(xeng)` **is** `cypstruct.xengine.select()` | max abs diff < 1e-9 on all 87. **Re-verified, not inherited** |
| **C-N2** | the shipped board reproduces | selected 0.6164, oracle 0.6975, random 0.5769, gain +0.0395, ρ −0.2582, correct sign 75.86% — to 5e-4 |
| **C-NUM** | `FINDING_021` numbering | `mapped` True on 1,740/1,740; exactly **4** exact-zero LDDT-PLI rows, all at BiSyRMSD > 10 Å → ejections, not numbering |
| **C-REF** | reference depth ≥ 4 per `FINDING_011` | min ≥ 4 |
| **C-ATOM** | every pose of a ligand has the same atom count (required for the rotation metric's identity correspondence) | 87/87 ligands homogeneous |
| **C-FILT** | every filter fires, with counts, **including the zeros** | reported as counts; 0 ligands dropped for no reference, 0 for no truth |

A control that reads zero rows is a failure, not a pass (`too-clean-numbers-are-the-tell`).

---

## 4. The clustering definition, fixed here

Poses are clustered **per ligand**, in that pose's **own heme frame** (`in_heme_frame`,
superposition-free).

**Primary distance — pose-to-pose Chamfer, `dC`.** The same permutation-invariant metric
the shipped selector uses. No atom correspondence is required, so `cyp-pose-atom-mapping-trap`
cannot apply. Chosen as primary because the hypothesis is about *what `xeng` sees*, and
`xeng` sees exactly this quantity against the reference set.

**Secondary distance — rotation angle after centroid alignment, `dR`.** `FINDING_024`'s
rotation measure: centre both ligands on their own centroids, Kabsch with **identity atom
correspondence** (legal within one ligand's pool — all 20 poses come from one generator
with one atom order, verified by C-ATOM), take
`angle = arccos((tr(R) − 1)/2)` in degrees. This isolates orientation, which is the 30°
error `FINDING_024` measured, and is reported as a robustness check on every endpoint.

**Linkage:** average, on the precomputed condensed distance matrix
(`scipy.cluster.hierarchy.linkage(..., method="average")`), flat clusters by distance
threshold (`fcluster(..., criterion="distance")`).

**Cut grid, fixed a priori and reported in full — no cut is tuned:**

* Chamfer: **0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 2.50, 3.00 Å**
* Rotation: **15, 30, 45, 60, 90, 120 degrees**

Every endpoint in steps 2 and 4 is reported **at every cut**. If the answer depends on the
cut, that is stated plainly as the result. Where a single number is needed, the headline is
the **median over the Chamfer grid** and the **maximum** is quoted separately as an upper
bound.

---

## 5. The measurements, in order

### Step 1 — where the oracle pose sits in the ranking

Sort each ligand's 20 poses ascending by `xeng` (lower is better; `argsort(kind="stable")`,
matching `FINDING_036`). Define:

* `rank_oracle` = 1-indexed position of the pose with the **maximum** LDDT-PLI. If several
  poses tie at the maximum within 1e-9, take the **smallest** rank — **conservative
  against this hypothesis**, since it makes the incumbent look better.
* `stake_all` = per-ligand (oracle − incumbent's pick). This is the per-ligand share of the
  0.0811 gap.

Report the **full distribution**: the count at each rank 1..20, and the fractions at
rank 1, ≤ 3, ≤ 5, > 10 (bottom half) and = 20 (dead last). Report the mean and median
`rank_oracle`, and the mean `rank_oracle` **weighted by `stake_all`**.

**Stake strata, fixed a priori.** HIGH = `stake_all` > **0.05** (the threshold
`FINDING_036` used for "ligands with stake > 0.05"); LOW = ≤ 0.05. Terciles of `stake_all`
are reported beside them. A ligand whose poses are all equivalent must not dilute one that
matters, so **every step-1 and step-2 endpoint is reported pooled AND within the HIGH
stratum.**

Null for the rank distribution: under no information, `rank_oracle` is uniform on 1..20
(mean 10.5). Report a one-sample test against uniform.

### Step 2 — is the oracle pose a minority report?

At every cut, for every ligand, record for each pose its cluster label, its cluster size,
and its `xeng` rank. Endpoints:

* **E2a** — fraction of ligands whose **oracle pose** lies in the **largest** cluster,
  against the matched chance rate `mean(size_largest / n_poses)`.
* **E2b** — the same for the **incumbent's pick** (`argmin xeng`).
* **E2c** — within-ligand Spearman ρ(`cluster_size`, `xeng_rank`), mean over ligands, the
  fraction of ligands with the predicted sign, and a binomial test on that fraction.
* **E2d** — mean cluster size of the oracle pose minus mean cluster size of the incumbent's
  pick, paired, with a bootstrap CI and a Wilcoxon p.
* **E2e** — E2a–E2d recomputed **within the HIGH-stake stratum** (conditioning on how much
  is at stake, as required).

**Pre-registered predictions, signs fixed now:**

| id | prediction | direction |
|---|---|---|
| **P1** | the oracle pose is in the largest cluster **less** often than chance | observed < matched chance |
| **P2** | `xeng` ranks by cluster membership | ρ(cluster_size, xeng_rank) **< 0**, correct sign on > 50% of ligands |
| **P3** | the oracle pose sits in a **smaller** cluster than the incumbent's pick | E2d **< 0** |
| **P4** | the effect **survives** conditioning on stake | P1–P3 hold in the HIGH stratum |

A failure of P2 falsifies the mechanism outright: if `xeng` does not rank by cluster size,
there is no majority vote to be a minority against.

### Step 3 — when is the consensus confidently wrong?

Prediction-side tightness only, nothing that needs the answer:

* `spread` = mean pairwise pose-to-pose Chamfer within the ligand's pool.
* `n_clusters` at the median Chamfer cut (1.25 Å).
* `margin_12` = `xeng(2) − xeng(1)`.

Fixed a priori: **TIGHT** = bottom tercile of `spread`; **WRONG** = incumbent-selected
LDDT-PLI < **0.50** (below this pool's random baseline of 0.5769, rounded down to a round
number chosen before looking). Report the full 2×2 (TIGHT/WIDE × WRONG/RIGHT) with counts,
and the rate of TIGHT-and-WRONG.

Then: do any **prediction-side** covariates separate TIGHT-and-WRONG from TIGHT-and-RIGHT?
Fixed list, no additions after the fact: `spread`, `n_clusters`, `margin_12`,
`mean xeng` (the absolute agreement level), `n_heavy` (atom count from the frames),
`pred_fe_donor_median` and `pred_frac_coordinated` from `binding_mode_labels_cyp3a4.csv`.
Mann-Whitney per covariate, **Holm-corrected across the seven**, reported whether or not
anything survives.

This is descriptive. **No selector is built from it in this sitting**, and any covariate
that separates is a hypothesis for a later pre-registration, not a result.

### Step 4 — the term oracle, which is the number that decides everything

The rule under test is *"the consensus is wrong here — take the minority cluster instead."*
Its perfect form gets to choose the **cluster** with full knowledge of the truth, and then
must take the pose that a **prediction-side** rule names inside it. It does **not** get to
choose the pose; that would be the pool oracle (0.0811) and would test nothing.

* **O-A (the term oracle)** — per ligand, at each cut:
  `max over clusters c of  lddt[ argmin_{i in c} xeng_i ]`, minus the incumbent's pick.
  The representative of each cluster is fixed by `xeng`; only the cluster choice is oracular.
* **O-B** — the same restricted to a **binary** choice: keep the incumbent, or switch to the
  single best **non-incumbent** cluster. (Identical to O-A by construction; reported so the
  framing is explicit.)
* **O-C** — the full pose oracle, **0.0811**, quoted for scale.

**Two bars, both recomputed on this pool, neither quoted:**

1. **F87** — the `FINDING_007` random-feature floor: argmin of a random within-ligand
   feature over all 20 poses, **4,000 draws**, p95 and p99. Expected ≈ +0.0134;
   it must **reproduce** or the run is void.
2. **N-A, the matched null** — O-A is a **best-of-`n_clusters`** oracle and is therefore
   positive by construction. The matched null **shuffles the cluster labels within each
   ligand, preserving the cluster-size multiset**, and recomputes O-A identically,
   **4,000 draws**, p95. This is the bar that asks whether clustering by *orientation*
   beats clustering at *random into the same-shaped groups*.

### Acceptance rule — decided now

| outcome | condition | verdict | consequence |
|---|---|---|---|
| **CLOSED** | **max** over the Chamfer cut grid of O-A **< F87 p95** | **REFUTED** | the ceiling is below the floor. No selector is built, now or later, on this pool. The line closes and that is a complete answer |
| **MEASURED** | O-A clears F87 at some cuts, but the **median-over-grid** O-A < 2 × F87, **or** O-A fails N-A p95 at more than half the cuts | **MEASURED** | report the ceiling and stop. No selector |
| **LICENSED** | median-over-grid O-A ≥ 2 × F87 **and** O-A > N-A p95 at ≥ half the cuts | **MEASURED, licensed** | a follow-up pre-registration may build a prediction-side minority-rescue rule. **Still nothing is built in this sitting** |

**No selector is built in this sitting under any outcome.** Step 4 is a ceiling, and
`FINDING_036`'s lesson is that a ceiling below the floor closes a line more cheaply and
more honestly than any candidate sweep.

---

## 6. What I will conclude if the ceiling is below the floor

That the minority-report account is **unrecoverable on this pool even if true**. It would
remain possible that the oracle pose *is* a minority report — steps 1–3 can still show
that, and it is worth knowing — while the corresponding rescue rule is worth less than
noise, exactly as a perfect top-2 tie-break was in `FINDING_036`. I will say so plainly,
record the descriptive result, and **close the line**. "The mechanism is real and the prize
is not" is the outcome this pre-registration expects most, and it is a complete answer.

---

## 7. Things this experiment will NOT do

* **No answer-recognition (R1) gate on any comparator.** `FINDING_036` §3 established that
  R1 tests whether a **prior** recognises the truth and is invalid for a **within-ligand
  comparator**; the shipped selector fails it at p = 0.0012 while working. Applying it here
  would condemn a working feature. The crystals are used in step 1–4 only as **truth for
  the oracle**, never as an object to be ranked.
* **No re-test of a refuted family.** Clash and strain (`FINDING_029`, anti-selective at
  −0.0632), the physics ensemble (025), Boltz confidence (001-E), ATOMICA and trunk
  embeddings (023), fragment placement (030), ligand torsion (032), and top-k tie-breaks
  (036) are all closed and none is touched.
* **No fitted parameters and no tuned cut.** The cut grid is fixed above and reported in
  full.
* **No conclusion from n = 14.**
* **No new inference, no downloads, no writes to `D:`.**

---

## 8. Outputs

| what | where |
|---|---|
| script | `scripts/structure/minority_report.py` |
| controls | `data/processed/minority_controls.json` |
| step 1, rank of the oracle pose | `data/processed/minority_rank.json` |
| steps 2 + 4, clustering and the term oracle, every cut | `data/processed/minority_clusters.json` |
| step 3, confidently wrong | `data/processed/minority_confident.json` |
| per-ligand table | `data/processed/minority_per_ligand.csv` |
| write-up | `docs/FINDING_037_minority_report.md` |

Seed **20260923**, `numpy.random.default_rng`, fixed for every null.
