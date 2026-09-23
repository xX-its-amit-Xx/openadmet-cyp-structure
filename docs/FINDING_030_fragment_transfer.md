# FINDING 030 — the superfamily does not place a shared fragment tightly enough to rank a pose

**Date:** 2026-09-22 · **Status:** measured, zero new inference, zero new downloads, CPU-only ·
**Verdict:** **REFUTED**, on the pre-registered rule, and refuted with a mechanism rather
than with a p-value.

Pre-registered in `docs/PREREG_fragment_transfer.md`, committed at `fa6257f` before a
single number existed. Nothing below was tuned afterwards. Everything post-hoc is labelled
post-hoc.

---

## What was asked

Dock fragments and functional groups, dock to paralogs and homologs, and use those poses to
make medicinal-chemistry-style corrections to the pose on the real target.

FINDING 024 refuted every intervention that conditions the **protein** on CYP3A4 —
templates, cytochrome b5, orthologs — because the pocket is already right to 0.73 Å and
ρ(protein error, ligand error) = **+0.03**. This idea is different in kind: a **ligand-side**
prior. Nothing measured so far ruled it out, so it was measured.

The cheapest strong version, and the one that has to work first:

> A functional group that a P450 binds sits somewhere specific relative to the iron and the
> porphyrin plane. That placement is an empirical crystallographic fact over 185 distinct
> targets, the heme frame is the one landmark that exists in all of them, and it is
> available for any query ligand that shares a substructure with something already
> crystallised. **Score a pose by how well its copy of a shared fragment matches the
> crystallographic consensus.**

Scoring prior only. No pose was edited, rotated or rebuilt — FINDING 027/028 say the
co-folded pocket is rigid (0.077 Å of ligand-to-ligand side-chain motion against the
crystals' 0.723 Å) and already excludes the true ligand on 71% of poses, so a transferred
fragment dropped into it would clash. Establish signal as a score first.

## What was done

`scripts/structure/fragment_transfer.py`, four CPU stages, ~25 min on 10 cores.

| | |
|---|---|
| queries | 87 CYP3A4 ligands, `data/processed/validation_ligands.csv` |
| pool | `val87b_unsteered` — 87 × 20 = **1,740** Boltz-2 poses |
| donors | `p450_universe/p450_atoms.parquet` — **1,002** crystal ligand-chain observations; **936 usable** (354 CCD codes, 475 entries, **183 target keys**) |
| frame | `cypstruct.xengine.heme_frame`, unchanged — Fe at the origin, z along the porphyrin normal flipped away from the thiolate, x from the propionate oxygens |
| chemistry | RDKit mol built **in coordinate order** (3D connectivity perception, then `AssignBondOrdersFromTemplate` with the CCD SMILES), MCS at `CompareElements` / `CompareOrderExact` / rings-only, **≥ 6 heavy atoms** |
| feature | `d_i` = symmetry-minimised RMSD, **no superposition**, between the shared fragment of the pose in the pose's own heme frame and of donor *i* in donor *i*'s own heme frame. `frag_cons` = mean over legal donors (primary), `frag_best` = min, `frag_med` = median. **Lower is better**, fixed a priori |
| selector | within-ligand z, argmax, ties broken at random over 64 draws. **No fitted parameters** |

Donor mol construction dropped **66 of 1,002** observations: 42 template-assignment
failures, 21 atom-count mismatches, 3 with no CCD SMILES. One query (**MF8**) failed the
same way and is excluded throughout.

---

## Coverage first, as pre-registered

**Coverage is not the problem.** This was the expected failure mode and it did not happen.

| | ≥ 1 legal donor | ≥ 3 legal donors | median donors | p90 | max | zero |
|---|---|---|---|---|---|---|
| **primary (no CYP3A donors)** | **86 / 87** | **84 / 87** | **522** | 544 | 613 | 1 |
| with-3A (optimistic bound) | 86 / 87 | 85 / 87 | 644 | 670 | 742 | 1 |
| loose MCS arm | 86 / 87 | 86 / 87 | 728 | 729 | 730 | 1 |

Median shared fragment **7 heavy atoms** (mean 7.1) in the primary arm, 12 in the loose arm.
Donors are drawn from **170 distinct P450 targets** in the primary arm. Every query that has
a molecule has a crystallographic precedent for part of it, several hundred times over.

### Every leakage filter fires — C2

Counts are summed over the 87 × 936 = **80,496** query-donor cells.

| filter | rule | cells removed |
|---|---|---|
| L1 | the query's own PDB entry | **108** |
| L2 | the same CCD ligand code, any entry | **61** |
| L3 | count-ECFP4 Tanimoto ≥ 0.90 (fixed a priori) | **153** |
| L5 | `closest_fe` > 10 Å — not an active-site copy | **4,887** |
| L0 | target not classifiable to a UniProt | **258** |
| **L4** | **the whole CYP3A subfamily** (P08684/P20815/P24462/Q9HB55) | **12,265** |

L4 is the expensive one and it is the one that makes this leave-one-**TARGET**-out rather
than leave-one-ligand-out. The with-3A variant lifts only L4 and is reported below **as the
optimistic bound it is**, never as a headline.

---

## The headline

Primary arm, 84 queries with ≥ 3 legal CYP3A-free donors. Random baseline **0.5789**,
pool oracle **0.6947**.

| selector | selected | **gain vs random** | fraction of the null ≥ it |
|---|---|---|---|
| `frag_cons` (PRIMARY) | 0.5723 | **−0.0066** | 0.790 |
| `frag_best` | 0.5488 | **−0.0301** | 1.000 |
| `frag_med` | 0.5783 | **−0.0006** | 0.536 |
| `frag_combo` = `−z(xeng) − z(frag_cons)` | 0.6066 | +0.0277 | 0.0005 |
| **incumbent** `cypstruct.xengine.select()` | **0.6188** | **+0.0399** | 0.000 |

**Random-feature null recomputed on this pool and this subset** (2,000 draws, as
pre-registered, not FINDING 007's numbers reused): p95 **+0.0140**, p99 **+0.0193**,
max +0.0334.

**Bar 1 fails.** The primary feature is *below random*, and the best of the three is
indistinguishable from it. The only row that clears the null is `frag_combo`, and
`frag_combo` is **worse than the incumbent alone** (0.6066 against 0.6188) — it is the
incumbent diluted, and reading its +0.0277 as a result would be reading the incumbent's
+0.0399 through a fog. It is included precisely because a naive table would report it as a
win.

Within-ligand, the feature is at chance: mean ρ = **−0.019**, median −0.030,
**correct sign on 54.8%** of ligands, binomial p = **0.45**.

### The with-3A optimistic bound, and the loose arm

| arm | n | `frag_cons` gain | null p99 | exceedance |
|---|---|---|---|---|
| primary, no 3A | 84 | −0.0066 | +0.0193 | 0.790 |
| primary, **with 3A** | 85 | **+0.0048** | +0.0195 | 0.306 |
| loose MCS, no 3A | 86 | −0.0221 | +0.0178 | 0.996 |
| loose MCS, with 3A | 86 | −0.0218 | +0.0178 | 0.995 |

Letting the query's own subfamily donate — which is the most leakage-flavoured thing this
experiment can legally do — buys **+0.011** and still does not reach the 95th percentile of
the null.

**One cell of twenty is positive**: loose / with-3A / `frag_best` at **+0.0220**
(exceedance 0.002). Twenty cells were computed (5 selectors × 2 leakage arms × 2 chemistry
settings), so one at that level is expected; more to the point it dies at bar 2 anyway
(**Δ = −0.0187 against the incumbent, Wilcoxon p = 0.60, 38 better / 37 worse**). It is the
"find the nearest CYP3A4 crystal ligand" statistic wearing a fragment costume, and it is
reported here so it is on the record rather than in a drawer.

### Bar 2 — the paired test that decides it

Primary arm, same 84 ligands, identical poses, against the shipped
`cypstruct.xengine.select()`.

| selector | mean Δ | bootstrap 95% CI (10,000) | Wilcoxon | ligands tied | better / worse |
|---|---|---|---|---|---|
| `frag_cons` | **−0.0464** | **[−0.0779, −0.0158]** | **0.0019** | **5** | 26 / 53 |
| `frag_best` | −0.0699 | [−0.1022, −0.0396] | 0.00042 | 0 | 30 / 53 |
| `frag_med` | −0.0405 | [−0.0692, −0.0122] | 0.0025 | 6 | 27 / 51 |
| `frag_combo` | −0.0122 | [−0.0330, +0.0096] | 0.189 | 30 | 23 / 31 |

This is not the usual death by ties. **Only 5 of 84 ligands select the same pose**, 79 pairs
are informative, and the feature loses on 53 of them. It is a genuine, significant,
*negative* result, not an underpowered one — the opposite failure mode from FINDING 026 and
027, where the tie column was the whole story.

**Whole board** (86 ligands, uncovered queries falling back to the incumbent's pick):
0.5714 against the incumbent's 0.6168, Δ = **−0.0454**; with-3A 0.5818, Δ = −0.0350.
Oracle 0.6959.

**Bar 3, complementarity:** r = +0.602 (p = 1.4e-09) between per-ligand Δ and the
incumbent's headroom — and it is **not** evidence of anything, because Δ and headroom share
the term −(incumbent selection) by construction. It is reported for comparability with
FINDING 025's −0.047 and should be read as moot for a feature whose Δ is negative
everywhere.

**Leave-one-ligand-out is vacuous** and is said to be vacuous rather than performed
decoratively: no cross-ligand quantity is fitted anywhere in this work. Ties are broken at
random over 64 draws throughout.

### The controls

| control | required | measured |
|---|---|---|
| **C4** frame identity | this frame reproduces the shipped `xeng` | **1,740 / 1,740 rows, max \|Δ\| = 2.45e-07** |
| **C5** crystal re-run | the diagnostic run reproduces the primary `d` matrix | **max \|Δ\| = 0.0 exactly, 82 ligands** |
| **C3** non-constant feature | > 95% of queries | **100%**, smallest within-ligand sd 0.024 Å |
| **C2** every filter fires | non-zero counts | **all six**, table above |
| **C1** scrambled donors | matched must beat scrambled | **fires: −0.0066 vs −0.0241** |

**C1 fires, and it is the one interesting thing in the negative.** Replacing each matched
donor with a random legal donor and a random atom subset of the same size — same donor
count, same fragment size, same number of alternatives minimised over, 8 redraws, seed
20260922 — scores **−0.0241** (sd over redraws 0.0037). The matched feature beats it by
**+0.0175, about 4.7 sd**. So the substructure match is *not* decorative: knowing which
atoms correspond is worth something real. It is simply worth far less than the 0.0193 the
null demands, and both arms sit below random pose selection.

---

## Why it fails — four post-hoc measurements

Labelled post-hoc; none of these is a candidate feature and none may become one without a
fresh pre-registration. `data/processed/frag_transfer_diagnostics.json`.

### D4. There is no consensus to match

| quantity | value |
|---|---|
| RMS spread of **donor fragment centroids** about their own mean, in the heme frame | **4.48 Å** (median 4.40) |
| RMS spread of one ligand's **own 20 predicted centroids** | **0.58 Å** |
| mean transfer distance `d` | **5.01 Å** |

The prior scatters **eight times more widely than the thing it is asked to discriminate**.
CYP3A4's error is a 30° rotation of a ligand with a 4.78 Å radius of gyration — about
**2.5 Å** of ligand RMSD. A distribution whose own width is 4.5 Å cannot resolve 2.5 Å.
This is the physical content of the result: *P450s do not place a shared 7-atom fragment in
a common position relative to the heme.* The superfamily binds the same group in many
orientations, and averaging over 500 crystals returns the cavity, not a pose.

### D5. The true answer is not closer to the consensus than the wrong guesses are

The decisive one. The query's **own crystal ligand** — the answer — was scored on exactly
the same feature, by appending it as a 21st pose (atom order transferred by substructure
match, never by array position; 83 of 87 available, 82 evaluable).

| | |
|---|---|
| mean percentile of the crystal among its own 20 predicted poses | **0.558** |
| median percentile | 0.65 |
| crystal better than the median predicted pose | **45 of 82** |
| binomial p | **0.44** |
| mean Δ (crystal − mean predicted) | −0.072 Å |

**Chance.** The feature cannot tell the right pose from a wrong one *even when handed the
right pose*. No selector, weighting, calibration or donor-choosing scheme can rescue a prior
that ranks the truth at the 56th percentile — and this single number is what closes the
whole family of variants, not just the three that were run.

### D3. What the feature actually ranks by

Within ligand, `frag_cons` correlates with the pose's **centroid radius in the heme frame**
at ρ = **+0.363**, positive on **82.9%** of ligands (and with the radius of gyration at
+0.327, positive on 78%). It is, to first order, a *hug-the-iron-and-be-compact* term.
FINDING 024 measured the CYP3A4 centroid as already correct to 1.07 Å with isotropic
scatter and no direction — so the feature's dominant axis is precisely the coordinate that
carries no information here. It is FINDING 026's lesson inverted: there, the shipped
Chamfer turned out to be ranking by orientation when it looked like it ranked by position;
here, a term that was meant to capture orientation turns out to rank by position.

### D1. "Pick a better donor" is not the lever

Donor oracle — for each query, the single best of its ~470 legal donors, chosen *knowing the
answer*: **0.6942**. Best-of-N random rankers at the same N: mean **0.6999**, p95 **0.7000**.
**The donor oracle is below its own matched null.** Any reported "best donor" gain here
would be a counting artifact, which is trap 1 of `docs/README.md` in a new costume.

### D2. Nor is "use only the similar donors"

Keeping only the *k* most chemically similar legal donors (post-hoc, Tanimoto-ranked):

| k | 1 | 3 | 10 | 30 | 100 | 300 | all |
|---|---|---|---|---|---|---|---|
| gain | −0.0362 | −0.0298 | −0.0192 | −0.0266 | −0.0389 | −0.0254 | −0.0064 |

Every restriction is **worse** than using all of them. The failure is not dilution by
irrelevant donors; the closest analogues are no better placed than the far ones.

---

## Verdict

**REFUTED**, on the pre-registered rule: bar 1 fails (`frag_cons` −0.0066 against a
recomputed p99 of +0.0193), and bar 2 fails in the strong direction (−0.0464, Wilcoxon
p = 0.0019, only 5 of 84 ligands tied). The scrambled-donor control C1 *did* fire, so the
chemistry is doing something — it is simply an order of magnitude too little.

This is a **ligand-side** negative, which is new. Everything killed in FINDING 024 was
killed for conditioning a protein that was already correct. This one was killed for a
different reason: **the crystallographic prior it depends on does not exist at the required
resolution.**

## What it licenses, and what it closes

**Closed — fragment-transfer pose editing, for now.** The user's endpoint was med-chem-style
correction of the pose using transferred fragment geometry. FINDING 027/028 said it must come
second because the co-folded pocket is rigid and would clash. This says something stronger and
prior to that: **there is nothing to transfer.** An editor that moved a fragment toward the
superfamily consensus would be aiming at a 4.5 Å-wide cloud to fix a 2.5 Å error, and D5 says
it would move poses away from the answer as often as toward it. Do not build it on this prior.

**What would have to be true to reopen it.** The donor-centroid spread must fall below about
1 Å — a factor of 4.5. Nothing in the *ligand* can do that, because the spread was measured
over an exact substructure match. It would require conditioning the placement on something
else: the surrounding pocket. That is the protein, and FINDING 024's ρ = +0.03 says the
protein is already right and uninformative here. **The two vetoes meet**, which is why this
is written as closed rather than as unfinished.

**What is NOT refuted.** (a) No docking was run: this tests transfer of *deposited*
fragment geometry, not a fragment docking calculation. (b) It tests *absolute* placement in
the heme frame, not *internal* geometry — torsion or conformer priors transferred from
crystals remain untested, and FINDING 024 left 26% of CYP3A4's error (1.60 Å) in the internal
conformer, which is the largest untouched term in the decomposition. **That is the version of
this idea worth the next sitting**, and it is a different measurement, not a re-tuning of
this one. (c) The with-3A bound is not a usable predictor but it does say the subfamily's own
crystals are only marginally more informative than the superfamily's, which is consistent
with FINDING 024's 1.00 Å crystal-to-crystal pocket spread.

**One methodological carry-forward.** D5 — scoring the *truth* on a candidate feature — cost
one extra stage and closed a whole family of variants that bar 1 and bar 2 alone would only
have wounded. Any future prior-shaped feature should be handed the answer and asked to
recognise it **before** a selector is built around it. Add it to the RUNBOOK's checks.

---

## Artefacts

| what | where |
|---|---|
| script, four stages + diagnostics | `scripts/structure/fragment_transfer.py` |
| pre-registration | `docs/PREREG_fragment_transfer.md` (committed `fa6257f`) |
| primary + with-3A results | `data/processed/frag_transfer_primary.json` |
| loose sensitivity arm | `data/processed/frag_transfer_loose.json` |
| post-hoc diagnostics D1–D5 | `data/processed/frag_transfer_diagnostics.json` |
| per-ligand table (82 rows) | `data/processed/frag_transfer_per_ligand.csv` |

Intermediates (framed poses, donor table, per-query match JSONs) are in the session
scratchpad on `C:`, never on the full `D:` drive. Nothing was downloaded and no inference
was run.
