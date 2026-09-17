# Handoff: custom scoring functions for CYP poses

**Paste everything below the line into your agent as its opening brief.** It is written to
be self-contained: an agent that reads it should be able to start work without asking
anyone for context, and — more importantly — without re-deriving results that already cost
us several thousand GPU jobs to learn.

---

## Who you are and what you are building

You are building a **custom scoring function for protein–ligand poses in cytochrome P450s**,
for the OpenADMET CYP blind challenge. You are working with the BU AI & ML team
(`buaiml/openadmet-cyp`), on what their task manifest calls **T6 — structure-derived
features**. Task T5 (generating the poses) is done; you are scoring them.

**Your input data is already public and needs no access request:**

- **Poses + labels:** https://huggingface.co/datasets/xX-its-amit-Xx/cyp-cofold-poses
- **Methods, code, and the full findings record:** https://github.com/xX-its-amit-Xx/openadmet-cyp-structure
- **Live campaign dashboard:** https://openadmet-cyp-dashboard.vercel.app
- **Pose inspector (3D, interactive):** https://openadmet-cyp-dashboard.vercel.app/structures

Reload the HuggingFace dataset whenever you need fresh data; it is the shared surface.
Upload your own poses there too, in the same layout, so the pipeline stays one bucket.

## What is in the dataset

| path | what it is | why you care |
|---|---|---|
| `p450_poses_scored.csv` | **~15,000 P450 poses scored against crystal truth** (LDDT-PLI, BiSyRMSD), across 87 proteins | This is your training/validation signal. It is the most valuable file in the set. |
| `t5_bakeoff/poses/{CYP3A4,CYP2D6}/*.cif` | 400 co-folded complexes, protein + HEM + ligand | Challenge compounds, both isoforms |
| `t5_bakeoff/t5_manifest.csv` | per-pose completion, physical validity, Fe geometry, model confidence | Ready-made features + a validity label |
| `reference_set_p450_sweep.npz` | frozen cross-engine reference, 492 ligands | Needed to reproduce the shipped selector |
| `pool_scorecard.json` | every headline number with its date and n | Check any claim below against this |

## The single most important thing to understand

**Selection is the bottleneck, not generation.** On our CYP3A4 set the *oracle* (best pose
available in the pool) scores **0.6975** while *selection* picks **0.5706**. That gap —
0.127 LDDT-PLI per ligand — is sitting in a pool already paid for. Your scoring function's
entire job is to close it.

The corollary is that **a better generator does not help you.** Two independent attempts to
add pose diversity (a second engine; a sampler sweep) each added ~+0.037 of oracle and
selection captured **−0.0017 and −0.0038** of it. More poses did not become findable poses.

## The rules of evidence here — read before you measure anything

1. **The noise floor is +0.0138.** That is the 95th percentile of a *random* feature on this
   data. Anything below **+0.020** is indistinguishable from noise. We confirmed this two
   ways: the direct null, and a formal max-of-N multiple-comparisons correction at N=30
   trials. They agree to 0.0008.
2. **Check within-ligand variance FIRST.** Selection compares poses *of the same molecule*.
   A feature with a huge between-compound effect and no within-ligand variance cannot rank
   anything. `is_coordinated` had a within-ligand CV of 0.023 and died on exactly this,
   after looking excellent between compounds (+0.14 LDDT-PLI). **One line of code, run it
   before anything else.**
3. **Rank correlation and top-1 selection have moved in opposite directions twice.** Judge
   on the metric that ships, which is the score of the pose you actually pick.
4. **Hold out proteins, not ligands.** We have 87; leave-one-target-out is the honest split.
5. **Write the kill criterion before the experiment.** Every idea in `docs/IDEAS.md` has one
   attached, and three ideas have been closed by them.

## What is already dead — do not spend time here

| approach | result | why |
|---|---|---|
| **Model confidence as a ranker** | within-ligand ρ = −0.033 / −0.092; **worse than random** | Protenix's seven confidence fields sit at frac(ρ>0)=0.50, all p>0.6. Two architectures, both at chance. If you weight poses by confidence you weight by noise. |
| Iron-coordination *as a filter* | −0.0002 | 84% of poses already sit inside any reasonable coordination window on crystal ligands — saturated, so the filter is inert. |
| QM scorer tier 1 (per-atom donor prior) | no within-ligand variance | 34 of 40 ligands coordinate through the *same* atom in every pose. |
| QM scorer tier 2 (AutoDock Vina terms) | +0.0205 vs null p99 +0.0248 | Plenty of variance, still cannot discriminate. |
| ~30 assorted features (occupancy priors, contact fingerprints, azimuth consensus) | all under the noise floor | See FINDING 002/003/006/010. |
| Portfolio-variance weighting (mean/std etc.) | +0.0008 over incumbent | `std` alone scores −0.0682: pose-to-pose spread is not a quality signal. |

## What works, and the bar you have to beat

**Cross-engine agreement.** Score each pose by its mean Chamfer distance, computed in the
heme frame, to poses of the *same ligand* from an independent source. Permutation-invariant
(no atom correspondence needed) and superposition-free (no reference frame to get wrong).
**No fitted parameters** — fitting weights collapses it into the noise.

- **+0.0381** on CYP3A4 vs the incumbent's +0.0265
- **+0.0357** across 81 held-out P450 proteins, p=0.0, positive on 26 of 30
- Survives multiple-comparisons deflation at N=100 trials

**Its mechanism is a catastrophe detector**, not a fine-grained ranker: the gain tracks how
bad the *worst* poses in a pool are. As catastrophe rates fell 19.3% → 4.13%, the gain fell
+0.3006 → +0.0357. **This is your opening.** A scorer that discriminates among *good* poses
would be complementary rather than competing, because the shipped one demonstrably does not.

Requires **≥4 genuinely independent reference poses**. At one it measures −0.0055 — actively
harmful, not merely weak. Depth saturates at 4; going to 8 buys +0.0030.

## Where the remaining opportunity actually is

1. **Substituent placement.** The iron anchor is saturated; the discriminating information
   is in where the rest of the molecule sits. Nothing anchor-local has ever worked.
2. **Orientation, not location.** Our failures are not misplacements. Open the pose
   inspector: two caffeine poses, centroids 0.38 Å and 0.53 Å from truth — a rounding error
   apart — scoring 0.908 and 0.320. The molecule is in the right place facing the wrong way.
   A scorer targeting orientation specifically is aimed at the real failure mode.
3. **Coordination has variance on YOUR compounds.** This reverses our own finding and is the
   most actionable thing here. On curated crystal ligands, 84% coordinate — saturated,
   useless. On the challenge compounds it is **41.0% (CYP3A4) and 28.5% (CYP2D6)**, because
   most are not Type II binders. **A feature that is dead on our data may be alive on
   yours.** Verify with rule 2 above before trusting it.
4. **Mode-specific scoring.** The catastrophe detector implies a bimodal population. Classify
   the failure mode first (wrong face / wrong rotamer / right pose), then apply a scorer per
   mode. Untested.

## CYP-specific physics you need to know

- **The iron has five ligands already:** four pyrrole nitrogens plus the Cys442 thiolate
  (Cys443 in 2D6, Cys435 in 2C9, Cys458 in 1A2). Only the **sixth, distal** site is free.
  **Any ligand density on the proximal face is physically impossible** regardless of what
  the model reports — we use this as a hard validity filter and it catches real failures.
- **Type II inhibitors coordinate the iron:** an sp2 nitrogen (imidazole, triazole, pyridine)
  donates into Fe(III) d_z². Crystal Fe–donor distances are **1.94 / 2.20 / 2.38 Å** at
  p5/p50/p95 over 116 deposited CYP3A4 entries; 83 of 101 unique ligands coordinate directly.
  It is a dative bond with a hard directional requirement, not a soft contact.
- **The pocket is plastic in one place: the F/G loop.** OpenADMET's own cryoEM work names
  F/G-loop remodelling as *why* co-folding does poorly here. It is also the region left
  unmodelled in many X-ray structures, so exclude it from rigid alignment.

## Four traps that corrupt results silently rather than failing loudly

These cost us real time. Each produces a plausible-looking wrong number.

1. **Ligand atom order differs between prediction and crystal.** A Boltz mmCIF may list
   C,C,C,C,N,C where the crystal lists O,C,O,C,C,S. Index-for-index comparison gives 11.4 Å
   where the truth is far smaller, silently. Use graph-isomorphism matching with element
   labels (symmetry-aware) — `cypstruct.pose.best_ligand_mapping` does this.
2. **`(pdb, chain, ligand)` is not a unique key.** Crystals deposit multiple copies — 8SO1
   has three caffeine molecules, only one in the pocket. We manufactured 125 duplicate rows
   this way. Include `seqid`, and pick the copy nearest the iron.
3. **Residue numbering differs.** Predictions number from 1; crystals use author numbering
   (3TK3 starts at 28). Align by residue *number*, never array position, or whole targets
   score exactly 0.000 — which looks like a catastrophic model failure rather than a join bug.
4. **Replicates are not independent opinions.** Both OpenProtein engines currently return a
   **byte-identical** pose per input, so `--replicates` buys nothing at any count, and
   `diffusion_samples` never varied the ligand. Historical pools are ~85% duplicates.
   **Deduplicate before counting reference depth.** Vary `num_recycles` / `num_steps`
   instead — that still produces genuine diversity (4 distinct placements from 4 settings).
   Keep `num_recycles ≥ 2`; at 1 it degrades by −0.0596.

## Navigating the dashboard

- **https://openadmet-cyp-dashboard.vercel.app** — campaign board. Approach cards for both
  tracks with status and headline metric; blockers; the metric tiles. Driven by
  `site/campaign.json`, so it is data rather than prose.
- **/structures** — the pose inspector. Rotate, toggle crystal vs best vs worst prediction,
  show pocket residues and the Cys442 thiolate, split side-by-side, and drag the morph
  slider to watch a prediction slide off the truth. Three detail levels (orientation →
  practitioner → full record) gate the jargon. **Start here** — it makes the
  orientation-not-location point in about ten seconds.

## Suggested first week

1. Pull `p450_poses_scored.csv`. Reproduce the random baseline and the +0.0138 null yourself.
   Do not take our numbers on faith; the repo records two claims we later retracted.
2. Reproduce the cross-engine selector at +0.0357 from the frozen reference. That is your
   incumbent and the thing to beat.
3. Compute within-ligand CV for every feature you are considering, *before* testing any.
   Discard the flat ones. This is 20 minutes and will save you weeks.
4. Target **orientation** and **substituent placement**, on **leave-one-protein-out** splits.
5. Report gains with the trial count attached — quote the N=30 deflated bar, not the
   single-trial p-value.

## How to reach us

Open an issue on `xX-its-amit-Xx/openadmet-cyp-structure`, or comment on
`buaiml/openadmet-cyp#4` where the T5 handoff and its corrections live. We are affiliated
with the BU team rather than on its roster, and we run the structure track ourselves — so
treat our results as an input to check, not an authority.
