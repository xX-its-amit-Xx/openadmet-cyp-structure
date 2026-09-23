# FINDING 031 — a second copy of the ligand goes into the active site, not the peripheral groove, and it makes the first copy worse

**Date:** 2026-09-22 · **Status:** measured, 48 OpenProtein jobs, ~22 minutes wall clock,
**$0** · **Verdict:** **REFUTED**, by the first pre-registered clause — the effect on the
scored copy is negative and its 95% CI excludes anything useful.

Pre-registered in `docs/PREREG_two_ligand_cofold.md`, committed at `f6ec85c` before a
single two-ligand pose was scored. No threshold below was moved afterwards.

---

## What was asked

> co-fold CYP3A4 together with the things it normally interacts with, to see whether that
> reveals cryptic pockets or different dynamics.

`docs/worldmodel/CYP3A4_BIOLOGY_MAP.md` §6.9/§9 ranks eight candidate partners and its
verdict is unambiguous: the highest-value one is **a second copy of the query ligand
itself**. It is the only "partner" with CYP3A4 structural precedent — 6 of 122 entries
carry more than one copy of one ligand per chain (2V0M, 4K9T, 4K9U, 8SO1, 8SO2, plus
8SG5/7SV2/8GK3 in the paralogs), all clash-free, minimum inter-copy heavy-atom distance
**2.90 Å** — the only one whose interface touches the **F/G roof** where `FINDING_028`
measured 83% of the pose exclusion, and it costs one extra ligand entity in the request.

The protein partners were not run and no budget went to them. `FINDING_024` Part 3 already
measured why: CYB5A and POR bind the **proximal face**, 9.7 Å from the iron and 10.7 Å from
the nearest modelled F/G residue, **across the porphyrin**, on a target where
ρ(pocket CA, ligand RMSD) = **+0.03**. There is no measured route from that face to this
error.

`FINDING_024` R3 proposed exactly this experiment, ranked it **last of three**, and asked
for it "as a probe with a pre-registered null, not as a pool change". That is what this is.

## What was done

`scripts/cofold/two_ligand_cofold.py` (submit / collect / score) and
`scripts/structure/two_ligand_analysis.py`.

| | |
|---|---|
| venue | **OpenProtein**. Modal is over its cap and was not touched. `cypstruct.budget.preflight` passed at 23 of a 2,000-job monthly cap |
| engine | `protenix_v2`, uploaded 6,979-sequence CYP3A4 MSA, `diffusion_samples=1` (FINDING 009: samples do not diversify the ligand), `num_recycles=3` |
| ligands | **15**, by the pre-registered rule: exclude the curated cryoprotectant/polymer component IDs (**2 fired** — PG0, PG4), then take the 15 smallest by RDKit heavy-atom count. 87 → 85 eligible → 15 |
| arms | **single** (protein+HEM+1×query) · **double** (protein+HEM+**2×query**) · **decoy** (protein+HEM+1×query+1×**ibuprofen**) |
| depth | 4 replicate jobs per (ligand, arm) = **180 complexes, 48 jobs, 0 failures** |
| scored copy | **the copy nearest the heme iron** — prediction-only, as FINDING 024 R3 specified. The higher-LDDT-PLI copy is reported only as a ceiling |

**Feasibility was established first, with one job.** The server accepts a four-chain
complex, and for caffeine it placed the two copies at **3.51 Å and 12.7 Å** from the iron,
6.43 Å apart — the 8SO1 arrangement, unprompted.

### Cost, stated before anything is proposed

**49 OpenProtein jobs** (48 + 1 feasibility probe), 181 complexes, **~22 minutes** of wall
clock from first submission to last collection, **$0 billed** — OpenProtein is unmetered
for this account, and the guard that matters is the job cap, of which this pilot used
**2.5%**. 61 MB of mmCIF on disk. **Zero Modal spend.** Scaling this to all 85 eligible
ligands would cost roughly **270 jobs and two hours**, which is why the result below is
worth having before spending it.

---

## Controls, with counts, including the zeros

| id | control | result |
|---|---|---|
| **C1** | FINDING 021 numbering | renumbering offset **0 on all 180 poses**, residue-name identity **1.000 minimum**. 4 exact-0.0 LDDT-PLI values, **all 4** with BiSyRMSD 17.3 Å and the ligand 20.0 Å from the iron — genuine ejections, not the numbering bug. All 4 are one ligand (MF8) in the **decoy** arm, one per replicate |
| **C2** | distinct poses, not jobs | 4 nominal replicates are **1.96 distinct poses** per (ligand, arm), min 1, max 2. 180 rows → **88 distinct**. Every paired test is repeated on the deduplicated pool |
| **C3** | both copies really present | **120 of 120** two-entity rows carry two entities; **60 of 60** double rows carry two copies of the query. Inter-copy minimum heavy-atom distance: median **3.11 Å**, p5 **2.47 Å**, and **4 of 60 below 2.5 Å** — below the crystallographic floor, where **no** deposited entry has an inter-copy pair under 2.5 Å |
| **C5** | symmetry mapping | `best_ligand_mapping` succeeded on **180 of 180** |
| **C6** | the two copy-choice rules | the nearest-iron copy is **also** the better-matching copy on **60 of 60** double poses, maximum gap **0.000000** |

**C6 is the reason two rows of this document are numerically identical, and it was checked
rather than reported.** The pre-registration defined a prediction-only rule and a 2-way
oracle over the copies; they never disagree, so the oracle buys exactly nothing and the
"ceiling" is the result. That is a fact about the prediction — the far copy is never the
one that matches the crystal — not a copy-paste.

**C2 is FINDING 015 in a new costume and it is the number to carry forward.** protenix_v2
was recorded as fully deterministic across replicates. In *this* configuration it is
**bimodal**: it returns one of exactly two answers, so four jobs buy two opinions. Not
deterministic, not diverse. Every depth figure here is quoted after deduplication.

---

## Result 1 — the primary endpoint: the first copy gets **worse**

Paired over the 15 ligands, per-ligand mean within each arm, on the copy nearest the iron.
The sign columns are `n_pos / n_neg`, not "better/worse": for RMSD and angle a positive
difference is worse.

**Arm B (second copy of the same ligand) minus arm A (single):**

| metric | arm A | arm B | **Δ** | bootstrap 95% CI | Wilcoxon p | +/− |
|---|---|---|---|---|---|---|
| **LDDT-PLI** | 0.6182 | 0.5518 | **−0.0664** | **[−0.1424, +0.0047]** | 0.135 | 6 / 9 |
| BiSyRMSD (Å) | 2.819 | 3.061 | +0.242 | [−0.402, +0.877] | 0.303 | 10 / 5 |
| **rotation (°)** | 64.8 | 69.3 | **+4.45** | [−22.4, +31.8] | 0.330 | 11 / 4 |
| translation (Å) | 1.065 | 1.204 | +0.139 | [−0.184, +0.430] | 0.277 | 11 / 4 |
| conformer (Å) | 1.442 | 1.527 | +0.084 | [−0.177, +0.322] | 0.303 | 11 / 4 |
| pocket CA (Å) | 0.744 | 0.853 | +0.109 | [−0.003, +0.249] | 0.762 | 6 / 9 |

**Arm C (second copy of an unrelated ligand) minus arm A:**

| metric | Δ | bootstrap 95% CI | Wilcoxon p | +/− |
|---|---|---|---|---|
| **LDDT-PLI** | **−0.0929** | [−0.1995, −0.0018] | 0.107 | 4 / 11 |
| BiSyRMSD (Å) | **+1.851** | [+0.332, +3.881] | **0.030** | 12 / 3 |
| rotation (°) | +1.54 | [−22.6, +27.5] | 0.890 | 9 / 6 |
| translation (Å) | **+1.989** | [+0.457, +4.260] | **0.012** | 11 / 4 |

On the **deduplicated** pool the picture is the same and the decoy arm sharpens:
Δ_B(LDDT-PLI) = **−0.0633** [−0.1357, +0.0003], p = 0.151; Δ_C = **−0.1205**
[−0.2300, −0.0138], **p = 0.041**.

### The pre-registered verdict

| clause | required | measured | fires? |
|---|---|---|---|
| SHIPS (1) | Δ_B ≥ +0.020 | **−0.0664** | no |
| SHIPS (2) | Wilcoxon p < 0.05 | 0.135 | no |
| SHIPS (3) | Δ_B − Δ_C ≥ +0.010 | +0.0265 | **yes** |
| SHIPS (4) | Δ_B(rotation) ≤ −5.0° | **+4.45°** | no |
| REFUTED (a) | Δ_B ≤ 0 **and** CI upper < +0.020 | −0.0664, CI upper **+0.0047** | **yes** |
| REFUTED (b) | Δ_C ≥ Δ_B | −0.0929 < −0.0664 | no |

**REFUTED**, by clause (a). The point estimate is negative and the confidence interval
excludes the smallest effect that would have been worth having.

Clause (b) does **not** fire, and that is the one nuance worth keeping: **the same ligand's
second copy is less damaging than an unrelated second entity** (−0.066 against −0.093, and
the decoy is the arm that reaches significance). The biology is not *absent* — a second
copy of the same molecule is better tolerated than a foreign one, exactly as the crystal
record would predict. It is simply not *helpful*.

**The within-ligand spread says how much of this to believe.** Mean within-cell sd is
**0.054 LDDT-PLI** and **15.5°** of rotation over ~2 distinct poses. Against that, a
−0.066 paired mean over 15 ligands is a real direction and a weak significance, which is
what the CI and the Wilcoxon between them say.

### The coordination statistic, which is the cleanest number here

Scored copy more than 5 Å from the iron — i.e. the query ligand pushed off the heme:

| arm | poses out of coordination | LDDT-PLI < 0.05 |
|---|---|---|
| **single** | **0 / 60** | 0 |
| **double** | **2 / 60** | 0 |
| **decoy** | **10 / 60** | **6** |

Ibuprofen **evicts the query ligand outright** on a sixth of the poses — MF8 goes to
LDDT-PLI 0.000 with BiSyRMSD 17.3 Å, TCI to 0.007 at 18.1 Å. A second copy of the *same*
molecule almost never does this. The control is doing its job: "any second entity perturbs
the prediction" is true, and it perturbs it **harder** than the hypothesis does.

---

## Result 2 — where the second copy goes: the **active site**, not the peripheral groove

Classified by the second copy's closest heavy atom against the three-site geometry measured
in biology map §7.1–7.3 (active site 0–11 Å, channel 11–15 Å, peripheral groove > 15 Å with
≥ 3 contacts on the measured groove set).

| class | arm B (2×query) | arm C (query + ibuprofen) |
|---|---|---|
| **active site** | **12 / 15** | 12 / 15 |
| channel | 2 / 15 | 1 / 15 |
| **peripheral groove** | **1 / 15** | 2 / 15 |

Median second-copy iron distance **8.81 Å** (arm B, per pose) / **9.75 Å** (per ligand).
The pre-registered bar for calling a cryptic peripheral site OBSERVED was **≥ 5 of 15**.
**Not observed.**

**But the 9.75 Å is not a null — it is a reproduction.** 2V0M's second ketoconazole sits at
**9.3–9.8 Å**, in the upper cavity at the channel mouth, and the biology map records
explicitly that it is "**not** the 1W0F peripheral site". The model puts the second copy
where the only multi-copy CYP3A4 X-ray structure of a large ligand puts it, to within
half an ångström of the median, without being told. What it does **not** do is populate the
17–21 Å progesterone/caffeine groove.

**And the F/G surface *is* engaged — from the inside.** 45 of 60 second copies make at
least one contact on the measured groove residue set, dominated by **Phe213 (31/60)** and
**Phe220 (28/60)**:

| residue | 212 | 213 | 214 | 217 | 219 | 220 | 238 | 239 | 240 |
|---|---|---|---|---|---|---|---|---|---|
| arm B contacts (of 60) | 3 | **31** | 0 | 5 | 6 | **28** | 0 | 3 | 8 |
| arm C contacts (of 60) | 14 | 18 | 3 | 7 | 8 | 4 | 8 | 4 | 4 |

This is the biology map's own point made by the model: the peripheral groove, the oligomer
interface and the F/G loop are **one surface**, so a ligand at 9 Å from the iron can stack
on Phe213/Phe220 while still being an active-site ligand by distance. The engagement is
real. It does not help.

### The pre-registered follow-up: does an active-site second copy help or hurt?

| where the second copy landed | n | Δ LDDT-PLI | Δ rotation |
|---|---|---|---|
| **active site** | 12 | **−0.0794** | **+4.57°** |
| channel | 2 | −0.0158 | −0.56° |
| peripheral groove | 1 | −0.0116 | +13.03° |

**Hurts, and the damage is concentrated exactly where the competition is.** The three
ligands whose second copy stayed out of the active site lose almost nothing; the twelve
whose second copy competes for the cavity lose 0.08 and turn 4.6° further wrong. This is
the direct answer to the cryptic-pocket question, and it is the one the mechanism predicted
if a second copy were being used as filler rather than as a constraint.

---

## Result 3 — the F/G loop moves, and not enough to call it

Backbone `N, CA, C, O` over residues **210–216**, 28 atoms on 14 of 15 ligands.

| quantity | median | mean | max |
|---|---|---|---|
| arm A vs **arm B**, superposed on the shared pocket CA | **0.401 Å** | 0.847 | 2.326 |
| arm A vs **arm C** | **0.267 Å** | 0.472 | 1.539 |
| arm A vs arm B, fit **excluding** 210–216 (sensitivity) | 0.416 Å | 0.899 | 2.485 |
| arm A vs arm C, fit excluding 210–216 | 0.274 Å | 0.507 | 1.704 |
| arm A → crystal | 0.723 Å | 1.702 | 6.358 |
| arm B → crystal | 0.825 Å | 1.836 | 6.082 |
| arm C → crystal | 0.599 Å | 1.557 | 6.259 |

Pre-registered rule: MOVED requires median ≥ **0.50 Å** and greater than the decoy arm's.
The second condition holds (**0.401 > 0.267**); the first does not. **NOT MOVED**, by
0.099 Å.

Two honest readings, both belong here. The direction is right and the specificity is right
— a second copy of the *same* ligand moves this span **1.5×** as far as a foreign one, and
three ligands (08J 2.33 Å, A1ASQ 2.31 Å, A1ASV 2.15 Å) move it a great deal. And it does
not clear a bar that was set from FINDING 028's own measurement of how rigid this pocket
is, and **the movement goes the wrong way against the crystal**: 0.723 → 0.825 Å. The one
intervention that could plausibly have loosened residues 211–216 nudges them, and nudges
them away from the answer.

---

## Result 4 — oracle first, then selection

On the **deduplicated** pool, because an oracle over duplicates is an oracle over one pose
written twice (FINDING 015).

| pool | distinct poses / ligand | **oracle** | random | **selected** | gain vs random | own null p99 |
|---|---|---|---|---|---|---|
| arm A only | 1.93 | **0.6608** | 0.6296 | **0.6357** | +0.0061 | +0.0242 |
| **arm A + arm B** | 3.93 | **0.6845** (+0.0237) | 0.5948 | **0.6330** (−0.0027) | +0.0382 | +0.0478 |

**The oracle rises +0.0237 and selection falls 0.0027.** That is CLAIM F of `CLAUDE.md`, and
FINDINGS 013 / 016 / 027, for the fourth time — a second engine, a sampler sweep, 1.78 M
rigid rotations, and now a second ligand copy have each added ceiling that selection cannot
reach. Note also that **neither** pool's selection gain clears its own random-draw null at
the 99th percentile, at n = 15 with ~2–4 poses per ligand: this pilot is not powered to say
anything about the selector, and does not.

**A pool that merely gets bigger is not a result.** It is not reported as one.

---

## Verdict

**REFUTED.** The lever does not work, the pre-registered clause that fires is the one that
excludes a useful effect rather than the one that merely fails to find it, and the
mechanism behind the failure is visible in the subgroup table.

- **SHIPS** — fails on 3 of 4 clauses, including the sign of the primary endpoint.
- **REFUTED** — fires on clause (a): Δ_B = −0.0664 with a CI upper bound of +0.0047.
- **cryptic peripheral site** — **not observed**: 1 of 15 against a bar of 5.
- **F/G loop moved** — **no**: 0.401 Å against a bar of 0.50 Å, though specifically
  greater than the decoy arm's 0.267 Å.

## What this changes

1. **The biology map's top-ranked co-folding partner is spent, and cheaply.** Eight
   candidates were ranked; the one with structural precedent, F/G contact and the lowest
   cost has now been run and is negative. Ranks 3, 4, 7 and 8 were already DO NOT SPEND on
   the proximal-face geometry. **What remains from §6.9 is rank 2, membrane lipid** — the
   same shape of experiment (extra ligand entities, no second chain), aimed at restraining
   211–216 from outside rather than filling the cavity from inside. This result does not
   test it and does not forbid it, but it lowers the prior: the one thing measured here is
   that *adding mass to the distal side makes the pose worse*, and lipid adds mass to the
   distal side.
2. **Filling the cavity is not the same as constraining it.** FINDING 024 found an
   881 Å³ cavity and a 30° rotation error and asked whether a second copy would take up the
   slack. It takes up the slack — the second copy lands in the active site on 12 of 15 — and
   the first copy turns **4.6° further wrong**. A competitor for the same volume is not a
   constraint on orientation; it is one more body with no term determining its own
   orientation either.
3. **The model already knows the 2V0M arrangement.** Median second-copy iron distance
   **8.81–9.75 Å** against 2V0M's measured 9.3–9.8 Å, and Phe213/Phe220 stacking on 45 of
   60 poses. The engine reproduces multi-occupancy geometry it was never asked for. That
   is a point in favour of the *structural* prior in the biology map and against the
   *predictive* value of acting on it.
4. **"Any second entity perturbs" is true, and it perturbs harder than the hypothesis.**
   Ibuprofen ejects the query ligand from the heme on 10 of 60 poses and to LDDT-PLI < 0.05
   on 6. Without the decoy arm, arm B's −0.066 would have looked like a specific biological
   effect. It is the *smaller* of two perturbations. **This is the control earning its
   place**, and any future multi-entity co-folding experiment here needs one.
5. **protenix_v2 on OpenProtein is bimodal, not deterministic.** FINDING 015 recorded it as
   deterministic on the P450 campaign; in this configuration four replicate jobs return
   **exactly two** distinct answers per cell (1.96 mean, max 2). Neither "replicates buy
   nothing" nor "replicates buy depth" is right. **Count distinct poses — and re-count them
   per configuration, because the answer changed.**
6. **The two copy-choice rules never disagree.** On 60 of 60 double poses the copy nearest
   the iron is also the copy that matches the crystal. Any future two-copy work can use the
   prediction-only rule without paying an oracle, and does not need to argue about it.

## What it forecloses

- **A second copy of the query ligand as a pool member or a generation-side intervention on
  CYP3A4.** Closed. It is negative on the score, positive on the rotation error, and its
  oracle gain is unreachable by the selector.
- **"A ligand in the cavity will push the F/G loop"** as an approach to residues 211–216.
  Measured at 0.401 Å of backbone motion, in the wrong direction against the crystal.
- It does **not** close: membrane lipid (§6.9 rank 2, untested), nor the same experiment on
  larger ligands — this pilot is the 15 *smallest* by construction, and
  `FINDING_024` measured ρ(n_heavy, score) = **+0.243** inside CYP3A4, so nothing here
  generalises to the other 70 without a second run. The cost of that run is ~270 jobs and
  two hours, and on this evidence it is not worth it.

## Files

| | |
|---|---|
| pre-registration | `docs/PREREG_two_ligand_cofold.md` (`f6ec85c`) |
| runner | `scripts/cofold/two_ligand_cofold.py` |
| analysis | `scripts/structure/two_ligand_analysis.py` |
| per-pose scores (180) | `data/processed/two_ligand_poses_protenix_v2.csv` |
| per-ligand table | `data/processed/two_ligand_per_ligand_protenix_v2.csv` |
| within-ligand spread | `data/processed/two_ligand_within_ligand_protenix_v2.csv` |
| F/G span | `data/processed/two_ligand_fg_protenix_v2.csv` |
| everything above, as JSON | `data/processed/two_ligand_analysis_protenix_v2.json` |
| the 180 mmCIFs (61 MB) | `data/processed/two_ligand/protenix_v2/` |
