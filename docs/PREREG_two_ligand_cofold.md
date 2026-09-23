# PREREG — co-folding CYP3A4 with a second copy of the query ligand

**Written:** 2026-09-22, before any two-ligand job was submitted and before any
two-ligand pose was scored. Committed before launch. **Not to be edited afterwards.**

---

## The question, in the user's words

> co-fold CYP3A4 together with the things it normally interacts with, to see whether that
> reveals cryptic pockets or different dynamics.

`docs/worldmodel/CYP3A4_BIOLOGY_MAP.md` §6.9 and §9 rank eight candidate partners and the
verdict is not ambiguous: the highest-value partner is **a second copy of the query ligand
itself**. It is the only "partner" with CYP3A4 structural precedent (6 of 122 entries
carry more than one copy of one ligand per chain — 2V0M, 4K9T, 4K9U, 8SO1, 8SO2, plus
8SG5/7SV2/8GK3 in the paralogs; all clash-free, minimum inter-copy heavy-atom distance
2.90 Å), the only one whose interface touches the **F/G roof** where `FINDING_028` located
83% of the pose exclusion, and it costs one extra ligand entity in the request — no second
chain, no new protein degrees of freedom.

The protein partners are **DO NOT SPEND**, for a geometric reason already measured:
CYB5A and POR bind the **proximal face**, 9.7 Å from the iron and 10.7 Å from the nearest
modelled F/G residue, across the porphyrin (`FINDING_024` Part 3). No budget goes there.

`FINDING_024` R3 proposed exactly this experiment, ranked it **last**, and called for it to
be run "as a probe with a pre-registered null, not as a pool change". This is that
pre-registration.

## What is already known that constrains the design

| fact | source | consequence for this design |
|---|---|---|
| `diffusion_samples` does **not** diversify the ligand on OpenProtein, for any engine | FINDING 009 | submit `samples=1`; buy depth only with replicate jobs |
| protenix_v2 is **deterministic** across replicate jobs | FINDING 015 | replicates are a *control*, not a source of depth. Determinism is an advantage here: a paired single-vs-double comparison has no sampling noise to average away |
| an uploaded MSA makes four OpenProtein engines fail server-side | FINDING 009 / CLAUDE.md | any engine other than protenix_v2 / protenix / esmfold2 runs with `Protein.NullMSA` |
| predicted and crystal residue numbering differ; the symptom is an **exact 0.0** LDDT-PLI | FINDING 021 | every pose is renumbered by residue-name agreement before scoring, and the offset and identity are reported per pose |
| CYP3A4's error is 74% rigid-body: **30.1° rotation**, 0.94 Å translation, 1.60 Å internal conformer, in a pocket built to 0.73 Å at CA | FINDING 024 | the endpoint that matters is the **rotation**, not the raw score |
| a pool that merely gets bigger raises the oracle and leaves selection flat or worse — three times | FINDINGS 013, 016, 027 | oracle and selection are both reported, **oracle first**, and a bigger pool is explicitly not a result |
| Modal is over its spend cap | CLAUDE.md | venue is **OpenProtein only**. `cypstruct.budget.preflight` is called before launch and its refusal is final |

## Pilot scope, fixed in advance

**15 ligands.** Not 87.

**Ligand selection rule — prediction-only, stated before any score was consulted.** From
`data/processed/validation_ligands.csv` (87 CYP3A4 ligand/entry pairs):

1. Exclude PDB chemical-component IDs in the curated cryoprotectant / polymer set
   `{PG0, PG4, PEG, P6G, 1PE, PE4, EDO, GOL, MPD, DMS, SO4, PO4, ACT, TRS, IMD}` — the
   same exclusion the biology map §7.1 applied over the whole PDB. A second copy of a
   polyethylene-glycol fragment is not the hypothesis under test.
2. Of what remains, take the **15 with the fewest heavy atoms**, counted from SMILES with
   RDKit, ties broken by ligand ID ascending.

Rationale, and it is the repo's own: `FINDING_024` R3 proposed "the 8 smallest-ligand
CYP3A4 pairs"; the smallest ligands are the ones that leave room for a second copy, and
the clearest multi-occupancy precedent in the PDB (8SO1/8SO2, three and six caffeines) is
the smallest drug-like CYP3A4 ligand there is. This is 024's rule, widened from 8 to 15.

**The set, fixed here:**

`MF8` (9), `CFF` (14), `MYT` (17), `A1ASS` (22), `A1ASV` (22), `08J` (23), `PK9` (23),
`TCI` (23), `A1ASU` (25), `CL6` (25), `A1ASP` (26), `A1ASR` (26), `D7Y` (27), `MWS` (27),
`A1ASQ` (28). Heavy-atom counts in parentheses.

**Known limitation, declared now:** the 15 smallest are, family-wide, *not* a random
sample — `FINDING_024` measured ρ(n_heavy, score) = **+0.243** inside CYP3A4, so small
CYP3A4 ligands score worse. That biases the *level*, not the *paired difference*, because
both arms see the same 15 ligands. Generalisation to all 87 is therefore not claimed by
this pilot and will not be claimed in the write-up.

## Arms

All three arms are the identical request except for the ligand entities.

| arm | chains | what it tests |
|---|---|---|
| **A — single** | protein + `HEM` + 1 × query | the paired baseline. Same engine, same settings, same day |
| **B — double-same** | protein + `HEM` + **2 × query** | the hypothesis |
| **C — double-decoy** | protein + `HEM` + 1 × query + 1 × **ibuprofen** | the null: does *any* second entity perturb the prediction? |

**Decoy, fixed in advance: ibuprofen**, `CC(C)Cc1ccc(cc1)C(C)C(=O)O`, 15 heavy atoms, a
CYP2C9 substrate rather than a CYP3A4 one, drug-like and lipophilic so it is a fair
competitor for the cavity rather than a straw man. **Verified absent from the 87-ligand
set** by canonical-SMILES comparison before this document was written (`False`). It is the
same decoy for every ligand, so the control is a constant.

**Engine:** `protenix_v2` on OpenProtein, with the uploaded 6,979-sequence CYP3A4 MSA,
`samples=1`, `num_recycles=3`, `replicates=2` — the same settings as the existing `op1`
pool, so arm A is comparable to work already on disk.

**Second engine, if and only if it accepts a four-chain complex:** `boltz2` on OpenProtein
in **single-sequence mode** (`Protein.NullMSA`), same three arms. Its purpose is to keep
the result from being a protenix quirk. If it fails server-side, that is reported as a
platform fact and the pilot stands on one engine.

**Budget.** `cypstruct.budget.preflight(..., venue="openprotein", cap_key="openprotein_jobs")`
is called before submission. If it refuses, the launch does not happen and the refusal is
reported. The pilot is ~23 jobs per engine per replicate-pass at 4 complexes per job.
Actual job count and wall clock are reported in the finding.

## Which copy is scored — decided now, not after looking

A two-copy prediction has two candidate poses for "the ligand". Two rules are
pre-registered and both will be reported:

- **PRIMARY — `nearest_fe`:** the copy whose closest heavy atom is nearest the heme iron.
  This is computable from the prediction alone, uses no crystal information, and is what
  `FINDING_024` R3 specified ("score only the copy nearest the iron"). **All headline
  numbers and all acceptance tests use this rule.**
- **SECONDARY — `best_match`:** the copy with the higher LDDT-PLI. This is a 2-way oracle
  and is reported only as a ceiling, never as a result.

Arm A has one copy, so both rules coincide there and the comparison is paired and fair.

## What is measured

### 1. Primary — does the FIRST copy's pose change, and change for the better?

Paired over the 15 ligands, arm B minus arm A, and separately arm C minus arm A:

- **LDDT-PLI** (`cypstruct.pose.lddt_pli`, symmetry-mapped, the shipped metric)
- **BiSyRMSD** (`cypstruct.pose.bisy_rmsd`)
- the **FINDING 024 decomposition**: pocket-lining CA RMSD, ligand **rotation error in
  degrees**, centroid **translation** in Å, and internal **conformer** RMSD, obtained by
  superposing on the binding site by residue number and then optimally superposing the
  ligand on itself under the scorer's own symmetry mapping.

**The rotation is the endpoint that matters.** 024 measured it at 30.1° against 6.1° for
the rest of the P450 family, and 027 showed the raw score can move for reasons that have
nothing to do with it.

### 2. Where does the second copy go?

Classified by the closest heavy atom of the *other* copy (the one not scored), using the
three-site geometry measured in biology map §7.1–7.3:

| class | rule |
|---|---|
| **active site** | Fe → closest atom **< 11 Å** |
| **channel** | 11–15 Å |
| **peripheral groove** | > 15 Å **and** ≥ 3 contacts within 4.5 Å of the measured peripheral set `{R212, F213, D214, D217, F219, F220, I238, C239, V240}` |
| **elsewhere / surface** | anything else |

Distances to the peripheral set are computed in the **prediction's own** numbering after
the FINDING 021 renumbering, so the residue identities are checked, not assumed.

If the second copy lands in the **active site**, the pre-registered follow-up is:
does the scored copy's rotation error go **down** or **up** on those ligands specifically?
That is the cryptic-pocket / competitive-displacement question and it is answerable
directly.

### 3. Does the F/G loop move?

Backbone atoms `N, CA, C, O` over residues **210–216** — the span `FINDING_028` identified
as the one where a repack fails on CB or main chain, and the span the model already
predicts to 1.03 Å median.

- `d_AB` = backbone RMSD of 210–216 between arm A and arm B, after superposing the two
  predictions on their shared pocket-lining CA set by residue number
- `d_AC` = the same for arm A vs arm C
- and each arm's own 210–216 backbone RMSD **to the crystal**

### 4. Pool effects, reported in the mandated order

**Oracle first, then selection.** Arm A alone (2 replicates) vs arm A + arm B pooled:
per-ligand maximum LDDT-PLI (oracle) and `cypstruct.xengine.select()` top-1 (selection),
against a random-draw baseline over the same pool. A pool that merely gets bigger is **not
a result** and will not be reported as one.

## Acceptance rules — fixed before any two-ligand pose was scored

Let Δ_B = mean paired (arm B − arm A) and Δ_C = mean paired (arm C − arm A), both on the
`nearest_fe` copy, over the 15 ligands.

**SHIPS** requires **all four**:

1. Δ_B (LDDT-PLI) ≥ **+0.020**
2. Wilcoxon signed-rank p < **0.05** on the 15 paired differences
3. Δ_B − Δ_C ≥ **+0.010** — the effect must be specific to *this ligand's* second copy
4. Δ_B (rotation error) ≤ **−5.0°**

**REFUTED** if **either**:

- Δ_B (LDDT-PLI) ≤ 0 **and** the 95% bootstrap CI upper bound on Δ_B is < +0.020 — i.e.
  a useful effect is excluded, not merely unobserved; **or**
- Δ_C ≥ Δ_B — any second entity does as much, so the biology adds nothing.

**MEASURED** in every other case, and the substantive reading is stated alongside the
formal verdict, as in FINDING 027.

**The F/G loop is called MOVED** only if median `d_AB` ≥ **0.50 Å** and median
`d_AB` > median `d_AC`. The 0.50 Å bar is set from FINDING 028's own numbers: the model's
pocket moves **0.077 Å** from one ligand to the next where the crystals move **0.723 Å**,
so anything under 0.5 Å is inside the model's measured rigidity and is not a dynamics
claim.

**A cryptic peripheral site is called OBSERVED** only if **≥ 5 of 15** second copies
classify as *peripheral groove* by the rule above. Fewer than that is compatible with
noise at n = 15 and will be reported as such.

Thresholds are round numbers chosen from prior measurements, not from any result in this
experiment. No cutoff in this document may be moved after a scored result is seen.

## Controls that must pass and will be reported with counts

| id | control | requirement |
|---|---|---|
| **C1** | FINDING 021 numbering | renumbering offset and residue-name identity reported for **every** pose; identity ≥ 0.95; no exact 0.0 LDDT-PLI accepted unless accompanied by BiSyRMSD > 10 Å |
| **C2** | determinism / distinct poses | md5 and per-atom coordinate sd across the 2 replicates of every (arm, ligand). **Distinct poses counted, never jobs** (README trap 1) |
| **C3** | both copies are really there | every arm-B and arm-C structure must contain **two** non-heme ligand entities with the expected heavy-atom counts; minimum inter-copy heavy-atom distance reported and compared against the crystallographic 2.90 Å floor |
| **C4** | every filter fires with a count | ligand exclusions, submission failures, collection failures, unscoreable poses — each reported as a number, including zeros |
| **C5** | symmetry mapping | `best_ligand_mapping` success rate reported; failures are excluded and counted, not silently index-mapped |

## What this cannot show

- n = 15 and the ligands are the smallest 15, so nothing here generalises to the 87-ligand
  set or to the challenge's cryoEM ligands without a second run.
- There is **no CYP3A4 crystal with two *different* drug-like ligands** (0 of 122, biology
  map §7.1), so arm C's geometry has no ground truth. It is a perturbation control, not a
  model of anything real.
- The engine is not the engine our pool came from. Every number here is arm B against
  **arm A of the same engine on the same day**, never against the Boltz-2 pool.

---

*Committed before submission. Scripts: `scripts/cofold/two_ligand_cofold.py`. Results:
`data/processed/two_ligand_*`. Write-up: `docs/FINDING_031_two_ligand_cofold.md`.*
