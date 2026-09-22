# PRE-REGISTRATION — fragment pose transfer from the P450 superfamily

**Written:** 2026-09-22 · **Committed before any score existed.** Nothing in this file may
be edited after the first number is computed. The write-up is
`docs/FINDING_030_fragment_transfer.md`.

---

## The idea

Dock-free, inference-free, ligand-side prior. FINDING 024 refuted every intervention that
conditions the **protein** on CYP3A4 (templates, cytochrome b5, orthologs) because the
pocket is already correct to 0.73 Å and ρ(protein error, ligand error) = +0.03. This is a
different kind of lever: it conditions the **ligand**.

> A functional group that a P450 binds sits somewhere specific relative to the iron and
> the porphyrin plane. That placement is an empirical, crystallographic fact, measurable
> across 185 distinct P450 targets in the heme frame, and it is available for any query
> ligand that shares a substructure with something already crystallised.

Tested here **as a scoring prior only**. Pose editing is explicitly out of scope — FINDING
027/028 say the co-folded pocket is rigid (side-chain ligand-to-ligand motion 0.077 Å
against the crystals' 0.723 Å) and already excludes the true ligand on 71% of poses, so a
transferred fragment dropped into that pocket would clash. Establish signal first.

---

## Data, fixed

| role | source | size |
|---|---|---|
| query ligands | `data/processed/validation_ligands.csv` | 87 CYP3A4 |
| query poses | `D:/cyp_scratch/val87b_unsteered` | 87 × 20 = 1,740 Boltz-2 |
| truth | `data/processed/poses_scored_val87b.csv`, `arm == "unsteered"` | LDDT-PLI |
| incumbent | `data/processed/xeng_val87b.csv` → `cypstruct.xengine.select()` | shipped |
| donors | `data/processed/p450_universe/p450_atoms.parquet` | 1,002 ligand-chain observations, 369 CCD codes, 494 entries |
| donor chemistry | `p450_cofold_set.csv` (`id` → `smiles`, `uniprot`) | 367 codes, 87 UniProts, 185 target keys |
| donor active-site geometry | `p450_geometry.parquet` (`closest_fe`) | 1,002 rows |

No new downloads, no new inference, CPU only.

## The common frame

`cypstruct.xengine.in_heme_frame` / `heme_frame`, unchanged: origin at Fe, z along the
porphyrin normal flipped away from the proximal thiolate, x from the propionate oxygens.
`p450_atoms.parquet` already stores `fe`, `normal`, `sg`, `xaxis` per row, so donors need
no recomputation; query poses are framed from their own predicted heme. **The heme is the
only landmark that exists in all 185 targets, and it is what makes cross-target transfer
meaningful.**

## Chemistry, fixed a priori

1. Build an RDKit mol for each ligand **in coordinate order**
   (`rdDetermineBonds.DetermineConnectivity` on the 3D heavy atoms, then
   `AssignBondOrdersFromTemplate` with the CCD SMILES). Coordinate order is preserved, so
   an MCS atom index is an xyz row index. A donor or a query whose bond-order assignment
   fails is **dropped**, and the count is reported.
2. Shared fragment = maximum common substructure, `rdFMCS.FindMCS`, **primary settings**:
   `atomCompare=CompareElements`, `bondCompare=CompareOrderExact`,
   `ringMatchesRingOnly=True`, `completeRingsOnly=True`, `timeout=10`.
   **Minimum 6 heavy atoms.**
3. **Sensitivity arm (`loose`), declared now, reported separately:** `bondCompare=CompareAny`,
   `ringMatchesRingOnly=False`, `completeRingsOnly=False`, minimum 5 heavy atoms.
4. Symmetry is handled by enumerating up to 64 substructure matches per side and taking
   the minimum over pairings (cap 4,096). No truth enters this minimisation.

## Features, signs fixed a priori

For pose *p* of query *q* and legal donor *i* with shared fragment *F*:

`d_i(p)` = RMSD, **without superposition**, between *F*'s atoms of *p* in *p*'s own heme
frame and *F*'s atoms of donor *i* in donor *i*'s own heme frame, minimised over symmetry
pairings.

| feature | definition | sign |
|---|---|---|
| **`frag_cons`** (PRIMARY) | mean over legal donors of `d_i` | **LOWER is better** |
| `frag_best` | min over legal donors of `d_i` | **LOWER is better** |
| `frag_med` | median over legal donors of `d_i` | **LOWER is better** |
| `frag_combo` (secondary) | `-z(xeng) - z(frag_cons)`, unweighted | HIGHER is better |

No fitted parameters anywhere. Terms z-scored **within** the ligand, because selection only
ever compares poses of the same molecule (CLAUDE.md, FINDING 011).

## Leakage rules — enforced and counted

Every filter below is applied per query, and the script reports how many donor rows each
one removes. **A filter that never fires is a bug, not a pass** (memory: too-clean numbers
are the tell).

| # | rule | rationale |
|---|---|---|
| L1 | drop donor rows whose `pdb` equals the query's own PDB entry | the answer itself |
| L2 | drop donor rows whose CCD `lig` code equals the query's | the answer under another entry |
| L3 | drop donors with **count-ECFP4 Tanimoto ≥ 0.90** to the query (Morgan radius 2, counts) — fixed a priori, not tuned | a near-identical analogue is the answer |
| L4 | **PRIMARY: drop every CYP3A donor** — UniProt `P08684` (3A4), `P20815` (3A5), `P24462` (3A7), `Q9HB55` (3A43), plus any entry whose target key maps to them. This is leave-one-**TARGET**-out, which 185 targets is what makes possible | the query's own subfamily |
| L5 | drop donor rows with `closest_fe > 10 Å` — a surface- or tunnel-bound copy is not an active-site pose and carries no placement prior | not the phenomenon |

The **with-3A variant** (L4 lifted, L1–L3 and L5 kept) is computed and reported separately
and **only** as the optimistic bound it is. It is never the headline.

## Coverage is reported FIRST

Before any selector number: how many of the 87 queries have **≥ 1** legal donor, how many
have **≥ 3**, and the full donor-count distribution. The primary analysis runs on queries
with **≥ 3 legal donors** (a mean over fewer is not a consensus; `xengine` demands ≥ 4
independent poses for the same reason). **If most queries have no legal donor, that is the
finding and it is stated plainly before any selector number.**

Two evaluation scopes, both reported:

* **covered subset** (primary) — queries with ≥ 3 legal donors;
* **whole board** — all 87, uncovered queries falling back to the incumbent's pick.

## Controls that must fire

| # | control | pass condition |
|---|---|---|
| **C1** | **scrambled donors** — for each query, replace each matched donor with a donor drawn uniformly at random from the legal pool and a random *k*-atom subset of it (same *k*, same donor count), averaged over 8 redraws, seed `20260922` | matched must beat scrambled. **If it does not, the chemistry is doing nothing and that is reported early and as the result.** |
| **C2** | every filter L1–L5 fires with a non-zero count | reported as a table |
| **C3** | `frag_cons` has non-zero within-ligand sd on > 95% of queries | a constant feature selects at random and would read as a null for the wrong reason |
| **C4** | the pool truth used here reproduces `poses_scored_val87b.csv` exactly on the 1,740 rows, and the incumbent reproduces the shipped `xeng` selection | the thing under test is the shipped metric and the shipped selector |

## Bars, in order — no skipping to the last

1. **Versus random.** Random-pose baseline = mean over ligands of mean LDDT-PLI. Gain must
   exceed the **99th percentile** of a **2,000-draw random-feature null computed on THIS
   pool and THIS ligand subset**. The +0.0138/+0.0196 figures of FINDING 007 are for other
   features and are recomputed here, not reused.
2. **Versus the incumbent, PAIRED, on identical poses.** `cypstruct.xengine.select()` on the
   same ligands. Report mean Δ, **10,000-draw bootstrap 95% CI**, **Wilcoxon signed-rank on
   the untied pairs**, and **the number of ligands tied**. Seven candidates have cleared
   bar 1 and died here. A large raw gain with 80 of 87 ligands tied is not a result.
3. **Complementarity.** Pearson r between the per-ligand Δ (frag − incumbent) and the
   incumbent's headroom (per-ligand oracle − incumbent selection). FINDING 025 died at
   r = −0.047.
4. **Pool oracle** reported alongside every selection number. **Within-ligand Spearman ρ**
   and **correct-sign fraction** reported alongside every pooled number. Ties broken at
   random, averaged over **64 draws**.
5. **Leave-one-ligand-out** is *vacuous* for an unfitted feature and will be said to be
   vacuous rather than performed decoratively (FINDING 027 precedent). If any variant does
   fit a parameter, it is scored leave-one-ligand-out and labelled as fitted.

## Verdict rules, fixed now

* **SHIPS** — bar 1 passes, bar 2 passes (bootstrap CI excludes 0 **and** Wilcoxon p < 0.05
  **and** fewer than 80 of the evaluated ligands tied), and C1 fires (matched > scrambled).
* **MEASURED** — bar 1 passes, bar 2 fails.
* **REFUTED** — bar 1 fails, or C1 fails (scrambled ≥ matched), or coverage is so thin that
  no bar is testable.

## Out of scope, deliberately

No pose is edited, rebuilt or rotated. No new inference. No new network fetch. If the
transferred geometry carries signal **as a score**, the write-up states exactly what pose
editing would then require; if it does not, pose editing is refuted with it.
