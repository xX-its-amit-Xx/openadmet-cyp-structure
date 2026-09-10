# A CYP-specific physics scorer for pose selection

**Status:** design, 2026-09-09. Tier 0 implemented; tiers 1–3 specified here.
**Objective it optimises:** LDDT-PLI against the released CYP3A4 experimental structures.
**Claim being tested:** for CYP, pose *selection* is the bottleneck, not pose *generation* —
and the selection signal the co-folders are missing is chemically specific enough to write down.

---

## 1. Why a generic docking score cannot do this job

Every widely used scoring function — Vina, Glide, GNINA, and the confidence heads inside
Boltz/AF3/Chai — treats the heme iron as a sphere with a charge and a Lennard-Jones radius.
That representation has no way to express the one interaction that decides CYP binding mode
for a large fraction of inhibitors: a **dative bond from a ligand sp2 nitrogen lone pair into
the Fe(III) d(z²) orbital**, at 2.0–2.3 Å, aligned with the proximal thiolate.

This is not a theoretical worry. In the sibling activity repo, an untemplated Boltz-2 pilot
cofolded 24 CYP ligands, 5 samples each, and **not one pose of 120 landed inside the true
Fe-coordination window** — the closest approach across the whole run was 2.80 Å. The failure
was geometric and systematic, in exactly the coordinate a generic score cannot see.

OpenADMET's own cryoEM announcement makes the complementary point from the experimental
side: the density for several ligands is consistent with **multiple mutually exclusive
conformations**, and the remodelling that reshapes the pocket is concentrated in the **F/G
loop** — the region left unmodelled in many X-ray structures. They name this as the reason
co-folding does poorly here.

Put together: the generators sample a plastic pocket badly *and* score the metal coordination
not at all. The second half is the one we can fix cheaply, so it is where the effort goes.

---

## 2. What is actually CYP-specific, term by term

Nothing below is a generic docking term rebadged. Each one encodes a fact about
cytochrome P450 chemistry that a general-purpose scorer has no way to represent.

### A. The sixth coordination site is the only one available

Fe(III) in P450 is already five-coordinate: four pyrrole nitrogens of protoporphyrin IX plus
the **proximal thiolate of Cys442** (CYP3A4 numbering; verified against UniProt P08684 in
`targets.py`). Only the distal, sixth site can accept a ligand.

Consequences we can score:

| observable | physical meaning | discriminating window |
|---|---|---|
| `d(Fe···N_lig)` | dative bond length | 1.9–2.4 Å for type II |
| `angle(S_Cys442–Fe–N_lig)` | *trans* to the thiolate | near 180° |
| `angle(lone-pair vector, N→Fe)` | the lone pair must point at the iron | < ~30° |
| `side(heme plane)` | ligand must be distal | hard filter — proximal is protein interior |
| `Fe out-of-plane displacement` | 6-coordinate low-spin pulls Fe into the plane | ~0–0.3 Å when coordinated |

The **S–Fe–N angle is the sharp one and essentially nobody checks it.** A pose can put a
nitrogen at a perfect 2.1 Å and still be wrong if it approaches off-axis; the trans
relationship to the thiolate is what makes the interaction real.

### B. Which nitrogen can coordinate is an electronic question

Type II binding is not "has an aromatic nitrogen". Imidazole coordinates strongly; an
ortho-substituted pyridine often does not; oxazole is weak. The discriminators are
**σ-donor strength** and **steric accessibility of the lone pair**, and both are computable
once per *molecule* rather than once per pose:

- **Proton affinity at each candidate N** (GFN2-xTB, ΔE of protonation) as a σ-donor proxy.
  Cheap, well-behaved, and monotone with coordinating ability across azole/azine families.
- **Condensed Fukui f⁻ / local softness** at each candidate N — nucleophilicity toward the metal.
- **Percent buried volume (%V_bur)** in a 3.5 Å sphere around the N, measured on the ligand
  alone — captures the ortho-substituent effect that kills coordination sterically.

These produce a **per-atom coordination prior** that the geometric terms are then read
against: a pose that coordinates through the molecule's *best* donor scores higher than one
coordinating through a hindered or weak nitrogen, even at identical distance.

### C. For substrates, the reactive atom should face the iron

Type I ligands do not coordinate; they sit over the iron at ~4.0–4.5 Å, positioned so the
**site of metabolism** is presented to the ferryl oxygen. The corresponding QM observable is
the **C–H bond dissociation energy** for H-abstraction (or an aromatic-carbon Fukui index for
oxidation), computed per atom, once per molecule.

The term is deliberately **rank-based, not absolute**: reward poses where the rank of an atom
by Fe-distance agrees with its rank by ease of abstraction. Rank form makes it robust to the
systematic error in any cheap BDE estimate, and it side-steps needing an absolute energy scale.

### D. Strain, clash, and stacking

- **Internal strain** — E(pose conformer) − E(relaxed), same method (GFN2-xTB, or MMFF for
  the fast pass). Above roughly 10 kcal/mol the pose is not a real bound state.

  ⚠️ **Strain has a negative precedent and should not be assumed to work.** The PXR campaign
  built an MMFF strain validator (#302), tested it against holo ground truth, and **killed
  it** — the artifacts record `"ACCEPT": false`, and a torsion-strain follow-on was rejected
  transitively. The campaign's own note is that *DFT would have computed the same energy
  ordering*, i.e. the failure was not the cheapness of the method but the premise that
  bound-conformer strain discriminates good poses from bad ones at all.

  Strain therefore enters this scorer as a **candidate term that must earn its place**,
  ranked below the coordination terms, and it is expected to gate null. If it does, that
  replicates PXR and gets logged as a replication rather than quietly dropped.
- **Hard-sphere clash** with protein and heme — a filter, not a soft term.
- **Porphyrin π-stacking** — CYP3A4 stacks many ligands flat on the heme. Scored as ring-centroid
  distance to the heme plane (3.3–4.0 Å) together with interplanar angle (< 20°).
- **Dispersion-dominated burial** — CYP3A4's cavity is large and hydrophobic; for substrates
  the binding energy is mostly dispersion. Captured by the tier-2 interaction energy, which
  must therefore use a **dispersion-corrected** method (GFN2-xTB includes D4-type dispersion).

### E. The polar anchors that CYP3A4 actually uses

Ser119, Arg106, Arg212, Asp214 recur across the CYP3A4 holo literature as the H-bond partners.
Scored as satisfied/unsatisfied donor-acceptor geometry, restricted to those positions rather
than counted generically.

---

## 3. Cost architecture — the part that makes this feasible

Scoring 10⁴–10⁵ poses with quantum chemistry is not affordable. It is also not necessary,
because **the expensive quantities are properties of the molecule, not of the pose.**

| tier | what runs | cost scales with | when |
|---|---|---|---|
| **0** | pure-geometry terms (§2A, §2D clash/stack, §2E) | poses | every pose, always |
| **1** | ligand-intrinsic QM: per-atom donor strength, BDE/Fukui, charges, %V_bur | **ligands** | once per molecule |
| **2** | GFN2-xTB interaction energy + strain on a pocket cutout | top-K poses | K≈5 per ligand |
| **3** | DFT on a truncated Fe-porphine + SH⁻ + ligand model | handful | only where coordination decides |

Tier 1 is the leverage. A pose pool of 200 ligands × 6 engines × 20 samples is 24,000 poses
but only **200 QM jobs**, because the per-atom donor and reactivity profiles are computed on
the free molecule and then *read* against each pose's geometry.

Tier 2 runs on Modal (GPU not required; xTB is CPU-bound and embarrassingly parallel).
Tier 3 is reserved and may never be needed — it is listed so the decision to skip it is explicit.

---

## 4. How the terms get combined — fitted, not hand-weighted

Hand-tuned weights are how a physics score quietly becomes a fitted model with unstated
parameters. Instead:

**Training data.** RCSB holds **122 entries for CYP3A4** (P08684, queried live). Each holo
entry with a real ligand gives a labelled example: run the same pool-generation pipeline on
its ligand, compute each pose's true LDDT-PLI against the deposited structure, and that is
the label. This yields thousands of (pose, score) pairs from public data alone, before the
challenge structures are released.

**Model.** A pairwise ranker (LambdaRank objective, gradient-boosted trees) over
`[tier-0 geometry] ⊕ [tier-1 ligand QM read against the pose] ⊕ [tier-2 energies] ⊕ [engine confidence]`,
trained to rank poses within a ligand. Ranking, not regression: we only ever need the argmax.

**Validation — the discipline inherited from the PXR campaign.**
- **Leave-one-ligand-cluster-out**, never random CV. Random splits were ~0.1 optimistic there.
- The scorer must beat the **cross-model z-hybrid** incumbent on held-out clusters, not beat nothing.
- Its gain must **correlate with where the z-hybrid errs**. A score that improves the cases
  already handled is absorbed and adds nothing at inference.
- **Beware leaky features.** Anything derived from the reference structure — including the
  pocket definition — must be computed from the *predicted* structure at inference time.
  This is the failure mode that inflates dev metrics and vanishes on the blind set.

**Deploy rule.** The physics scorer ships only if it wins on held-out clusters *and* the win
survives seed resampling. If it does not, that is a publishable negative and gets logged as one.

---

## 5. Handling the plasticity OpenADMET warned about

Their post says several ligands show **mutually exclusive conformations** rather than one pose,
and that the F/G loop does the remodelling. Two design consequences, both already reflected in
`targets.py`:

1. **The F/G span (residues 202–260) is excluded from the rigid alignment frame.** It moves
   legitimately between holo forms, so including it in a superposition contaminates every
   comparison with real biology mistaken for error.
2. **Multi-conformer scoring.** Where the pool contains two well-separated clusters that both
   score well, the right answer may be that both are real. LDDT-PLI is scored against a single
   reference, so we still submit one pose — but the *cluster structure* is a first-class output,
   and picking the cluster with more members and better physics beats picking a lone high-confidence outlier.

---

## 6. Honest statement of what could sink this

- **The coordination term only helps type II ligands.** If the released structures are mostly
  substrates, term A is inert and the scorer rests on strain/stacking/dispersion, which are
  much less CYP-specific. Mitigation: type the ligand set early and report the split.
- **Cheap QM may not resolve donor strength finely enough.** GFN2-xTB proton affinities are good
  for ordering across families, shakier within one. Tier 3 exists for this.
- **122 PDB entries is not 122 independent examples** — they cluster heavily by ligand series.
  Cluster-aware splitting will shrink the effective sample size a lot, and the honest number of
  independent test clusters may be small enough that no result clears significance.
- **If the pool never contains a near-native pose, no scorer can rescue it.** Pool oracle
  (best achievable LDDT-PLI in the pool) is measured first and reported alongside every
  selection result. Selection gains are meaningless without it.
