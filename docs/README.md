# Findings index

Forty findings, most of them negative. Read them in this order if you are arriving cold;
the numbering is chronological, not logical.

> # ▶ ON RELEASE DAY, READ `DROP_DAY_PLAYBOOK.md` FIRST — AND ONLY THAT
>
> **`docs/DROP_DAY_PLAYBOOK.md`** is the single document someone executes when the
> structure-track test set lands. It carries the executable sequence with real commands
> and a *done* condition per step, the generation venue decision (Explorer only, with the
> MSA-staging requirement, the sbatch template and all five documented failure modes),
> what to do about the test set's binding-mode composition, **the DEAD LIST** of every
> refuted lever with its mechanism and finding number, the two diagnostic rungs and their
> scope, the statistical bar, the one thing still live, the traps that have each cost real
> time, and **§9 — five pieces of the shipped path that are stale or blind-incompatible**.
>
> Everything below is the *why*. The playbook is the *what*. Under deadline pressure the
> expensive mistake is re-deriving or re-trying something already refuted — that is what
> the DEAD LIST exists to prevent.

## Start here

| # | one line | status |
|---|---|---|
| **001** | Selection, not generation, is the bottleneck — pool oracle 0.6975 vs selection 0.5706, and Boltz's own confidence ranks poses *worse* than random | stands |
| **011** | **Cross-engine agreement selects: +0.0381**, beating the incumbent's +0.0265. No fitted parameters. This is what ships | stands |
| **012** | It **generalises** — positive on 13 of 16 held-out P450 proteins — and it is a **catastrophe detector**, so its payoff is predictable from pool quality | stands |
| **024** | **CYP3A4 fails by ROTATION (30°) in a pocket the model builds right (0.73 Å)** — protein accuracy does not predict ligand accuracy here (ρ=+0.03) where it does everywhere else (ρ=+0.51). Kills templates, b5 and orthologs | stands |
| **027** | **The model's pocket excludes the true pose** — the crystal ligand clashes below 2.2 Å inside the co-folded protein on **71%** of poses, while **0 of 87** crystals violate that cutoff in their own protein. Only 2% of orientations fit; 31 of 41 failures are unreachable by any rigid rotation | stands |
| **028** | **Three residues do 83% of the exclusion — Phe215, Arg212, Phe304** — and the model's pocket is **rigid, not wrongly adaptive**: 0.077 Å of ligand-to-ligand side-chain motion against the crystals' 0.723 Å. Where a repack works it is **one residue turning 30°**; where it fails it is **CB or main chain at 211–216**, in the F/G loop | stands |
| **RUNBOOK** | what to run, in order, with each step's trap attached | live |
| **PLAYBOOK** | **`DROP_DAY_PLAYBOOK.md` — the release-day path, the DEAD LIST, and the gaps in the shipped pipeline** | **live** |

## The shape of the whole problem

Everything **anchor-local or population-level has failed**; only **within-ligand** comparisons
work. That is not a slogan, it is ~30 measured features:

| # | what was tried | result |
|---|---|---|
| 002 | first sweep of selectors | nothing beat random |
| 003 | pocket contacts + consensus medoid | **+0.0220 held-out** — the first that worked |
| 006 | contact-fingerprint consensus, azimuth consensus, crystallographic contact prior | three dead ends |
| 007 | **the noise floor**: a random feature scores +0.0138 at the 95th pct | anything under +0.020 is noise |
| 010 | occupancy / orientation priors from 747 crystal poses | null (+0.0018, +0.0008) |
| 014 | QM scorer tier-1 donor prior; then classical interaction energy | both fail the gate |
| 025 | **the whole CYP3A4 physics scorer above the iron** — Ser119/Arg106/Arg212/Asp214/Thr224 anchors, the Phe roof, F/G engagement, MMFF strain (34 terms) | every term null; the fitted ensemble (+0.0252) lands on its own 200-draw null's **maximum**, loses to the incumbent and fails the complementarity gate |
| 026 | splitting the shipped Chamfer into translation + orientation (R1 of 024) | orientation alone is the strongest single unfitted term (+0.0424) and still **ties the incumbent paired** (+0.0041, p = 0.61); translation alone is *not* a passenger (+0.0309) |
| 027 | expanding the pool with 1.78 M rigid ligand rotations inside the model's own protein | oracle **+0.0108**, selection **−0.0145**; within-ligand ρ *improved* −0.258 → −0.390 while top-1 got worse. The incumbent rejected **99.83%** of the injected poses |
| 028 | attributing 027's exclusion residue by residue | wrongness and blocking are **different residues** (ρ = +0.34, p = 0.06): Leu210 is off 100° and blocks nothing, Phe215 is off 8.7° and blocks most. Side-chain error correlates with pose error **within** a ligand (ρ = −0.31, 87% correct sign) and not **between** (−0.08) — 024's null again |
| 029 | **induced-fit DEMAND** — repack the pocket around each pose and make the repack cost the feature, the first scoring experiment run against a receptor that is allowed to move | the probe works (90° rotation: 0→4 residues, 2.52→0.79 Å, p=3e-142) and selects nothing: best single term **+0.0012** against a +0.0140 floor, the 44-column ensemble **+0.0228** at its own null's 99th pct, **−0.0167 paired** vs the incumbent with 5 of 87 tied, complementarity r = −0.078 |
| 030 | **fragment pose transfer from the P450 superfamily** — score a pose by how well its shared substructure matches where 183 other targets' crystals put that fragment in the heme frame | coverage is fine (522 legal donors per query, CYP3A excluded) and the prior is **empty**: donor fragment centroids scatter **4.48 Å** where the error to fix is 2.5 Å, and the query's own **crystal** ranks at the **56th percentile** of its own 20 predicted poses on the feature. −0.0066 vs random, **−0.0464 paired with only 5 of 84 tied**. Closes pose editing on this prior |
| 031 | **co-folding a SECOND COPY of the query ligand** — the biology map's top-ranked partner, the only one with CYP3A4 precedent and the only one touching the F/G roof | the second copy lands in the **active site on 12 of 15** (median 8.8 Å from Fe, reproducing 2V0M's 9.3–9.8 Å unprompted) and **not** the peripheral groove (1 of 15). The first copy gets **worse**: **−0.066 LDDT-PLI**, CI [−0.142, **+0.005**], rotation **+4.5°**, and the damage is concentrated on the 12 whose second copy competes for the cavity (−0.079). An **unrelated** second ligand is worse still (−0.093, and it ejects the query from the heme on 10 of 60 poses). F/G 210–216 moves 0.40 Å against a 0.50 Å bar. **REFUTED** |
| 032 | **the ligand's INTERNAL CONFORMER** — 024's last untouched term (26% of the error, 1.60 Å), transferred as torsions from the superfamily and, separately, as prior-free ETKDG plausibility | the prior is **sharp** this time (donors agree to **25.0°** where the error to fix is **59.3°**) and it points the **wrong way**: the crystal ranks at the **31st percentile** of its own 20 poses, p = 2.6e-04, and the same happens with a general small-molecule prior. **The torsion ORACLE — the true torsions handed over — is worth only +0.0213 and loses to the incumbent.** Gate failed, no selector built |
| 033 | **is the shipped selector biased by BINDING MODE?** — the validation set is 83% Type II, so a Type I-rich blind test set would be an unmeasured exposure | **ROBUST**, on the set with the power: family-wide the raw Type II − Type I gap is +0.0318 and **vanishes to −0.0030, CI [−0.0132, +0.0075], once matched on pool headroom**; both strata clear their own nulls at p=0.0000 over 342 pairs. CYP3A4 alone is **underpowered** — n=14, CI [−0.063, +0.050], MDD 0.066. The heme-frame explanation **fails** (partial ρ −0.033). The real exposure is the POOL: Type I random 0.464 vs 0.599 and **oracle 0.649 vs 0.707**. One label was wrong: **73/14, not 72/15** |
| 034 | **is FINDING 033's own recommendation right — buy DEPTH on the predicted-Type-I ligands?** The rate is real and the purchase is not | **the depth law holds and is steepest here**: Type I converts depth at **+0.0127 selected per doubling** (+0.0333 in the top octave) against Type II's +0.0069, because its oracle climbs **1.8× faster** (+0.0428 vs +0.0242) — while the *fraction* of depth converted is **29.6% vs 28.5%**, identical, and **22.8% vs 22.7%** across two different generators. But the only depth OpenProtein can sell is MSA-less `boltz2` (an uploaded MSA fails server-side), whose ceiling matches and whose average is 0.057 lower; mixed in it adds **+0.0118 oracle and −0.0065 score**, the fifth such case. **REFUTED** as a purchase, **MEASURED** as a law. Also: 30 replicates bought 10.8 distinct poses and replicates **16–29 bought zero** |

| 035 | **the experiment 034 asked for: depth at MATCHED conditioning.** Boltz-2 on Explorer, staged 6,979-sequence MSA, same checkpoint, same heme bond — 14 ligands × 20 new samples, a true 20→40 doubling | **the match is real** (new arm's pool mean **−0.0019** from the pool it joins, against OpenProtein's −0.057) and **the conversion is not**. Δ oracle **+0.0097** [+0.0012, +0.0195]; Δ selected **+0.0116** [−0.0687, +0.1034], Wilcoxon **p = 0.715**, **10 of 14 unchanged**, **+0.0020** once the two swing ligands are dropped. **But the sign flipped** — the first expansion in six that did not make selection *worse*. Two methods corrections: the **selected** rate decays with depth (+0.0227 at 10→20, **+0.0095** at 20→40) while the oracle rate barely does (+0.0385 → +0.0331); and **a subsample curve's last rung is a deterministic maximum**, so its final slope over-prices new samples — which is why 034's +0.0428 predicted rate met a +0.0097 reality. Also: **one job bought 20.0 distinct poses** where 034's thirty OpenProtein jobs bought 10.79 |

| 036 | **the experiment 035 asked for: break the near-tie.** Characterise the regime first, then five tie-breakers native to the shipped selector — medoid of the top-k, 2nd-nearest reference, Borda rank, within-pose sd across references, plurality vote | **the prize does not exist.** A **PERFECT** top-2 tie-break is worth **+0.0129** against a floor of **+0.0134** recomputed here; top-3 +0.0210, top-5 +0.0318. The near-tie is the *normal* case (median margin **0.0194 Å**, 70 of 87 under 0.06 Å) and the stake does **not** concentrate in it (ρ = −0.14, p = 0.18) — so it is a **variance mechanism**, not a prize. `argmin(xeng)` beats a coin flip over its own top 2 by **+0.0161** and beats **4,000 of 4,000** random tie-breaks. Fifteen of sixteen candidate configurations are negative; the best is **+0.0011, p = 0.68, 71 of 87 tied**. **60% of the 0.081 oracle gap sits outside the top 5** |

| 037 | **the experiment 036 asked for: why is the right pose ranked 8th?** Locate the oracle pose in the `xeng` ranking, cluster each pool by orientation in the heme frame, and price a perfect *"the consensus is wrong, take the minority cluster"* rule | **the median oracle pose IS ranked 8th**, and the ranking **inverts where the money is**: LOW stake median rank **3.5** and 33% at rank 1, HIGH stake median **11**, **2 of 45** in the top 3, mean rank **12.04 against a lottery's 10.5** (ρ(rank, stake) = +0.601, +0.327 de-confounded, p = 0.005). The **mechanism is confirmed on all four pre-registered predictions** — `xeng` ranks by cluster membership (ρ = −0.25 to −0.34, 69–85% correct, every cut of both metrics), the oracle pose is in a majority cluster **below** chance while the incumbent's pick is **above** it (14/14 cuts), and conditioning on the stake **doubles** it (**46.7% minority vs 23.1% chance**). And it is **not a prize**: the perfect rescue is **+0.0427**, its ceiling is a monotone function of the cut that converges on the pose oracle at singletons, the orientation-specific excess over a size-matched shuffle is **+0.004 / +0.015 median** against the **+0.0134** floor — and **choosing the MAJORITY mode is worth +0.0546 (p = 1.3e-06), more than a perfect rescue of the minority**, with a −0.151 downside. `FINDING_036` §2c one level up |
| 038 | **the last drop-day gap: reference-pose VOLUME.** 9 validation ligands stratified in advance by heavy-atom tercile, rotatable-bond extremes and prediction-side binding mode (3 of the 13 Type I), run through the shipped `submit --sweep` × 4 settings × 2 Protenix checkpoints → `collect` → `refset` | **it is no longer the longest pole.** **9 of 9 clear depth 4 AND depth 6**; median depth **8**, min 7. 16 jobs, **25.7 min**, $0. The mechanism is arithmetic: **depth = engines × settings**, minus one per failed job — **zero** duplicate attrition (`dropped_same_md5` = 0, `dropped_same_coords` = 0, both proved real by five positive controls; min pairwise Chamfer **0.288 Å** against a 0.05 Å tolerance). **Depth does not depend on chemistry** (ρ = −0.087 / +0.088, p = 0.82); the apparent Type I deficit is **one failed job**, because an ordered csv under `--batch` packing put all 3 Type I picks in it. Minimum viable sweep: **2 engines × 3 settings** (depth 6, survives one failure) — the 4th is slack. `jobs = 6 × ceil(N/5)`; 100 ligands in **≤ 2 h 25 m**. Live risk: a **6.25%** job failure rate that `submit` **cannot** resubmit (failed batches stay marked claimed) |

| 039 | **record repair, not science: the two decision-relevant contradictions in `DROP_DAY_PLAYBOOK.md` §10** — five quoted values of "the n=14 noise floor", and FINDING 025 explaining a live result with FINDING 021's retracted catastrophe count | the five floors are **three populations** (crystal-side d20 **+0.0431**, prediction-side d20 **+0.0449**, prediction-side d40 **+0.0440**, at 2×10⁶ draws) plus the **±0.0009** Monte-Carlo error of a 4,000-draw p95 — 034's low value replays bit-for-bit from its own seed. **One published claim flips:** 036's pool-B top-2 ceiling **+0.0447 clears** its floor at **p = 0.0475**, it does not "land on" it; 036's verdict, and 033/034/035's, stand. FINDING 025's **+0.0408 reproduces exactly** and its **mechanism is refuted** — a perfect catastrophe-avoider is worth **+0.0002**, deleting all 7 catastrophic poses leaves **+0.0409**, and pool uncertainty runs backwards; the result is now **unexplained**. Also settled: 035's two baselines are two estimands and its pool-30 peak is **exact, not noise, and is 08J alone**; aromatase's 0.863 and 0.631 are **sample-0 vs all-5 means** of the same 6 pairs |

| 040 | **the experiment 039 asked for: does the SELECTED depth curve turn over?** A second matched-depth stratum — the **73** prediction-side Type II ligands, the exact complement of 035's 14 — generated on Explorer at the same MSA, checkpoint and heme bond, and the curve computed in **closed form** | **NO TURNOVER. The exact union selected curve is STRICTLY increasing at all 39 steps** (min increment +0.00036, argmax depth **40**), and so is the combined 87-ligand curve: **039's pool-33 peak was 08J and does not generalise.** The conditioning match is the repo's best — **A0 = +0.0005** over 1,460 poses, against 035's −0.0019 and 034's −0.057 — and **depth still does not convert**: Δ oracle **+0.0172 [+0.0102, +0.0251]**, Δ selected **+0.0027 [−0.0130, +0.0182]**, **p = 0.870**, **38 of 73 unchanged**, **−0.0012** dropping the one biggest mover. **No octave clears +0.0125** (+0.0070 / +0.0039 / +0.0050 / +0.0077), the oracle-to-selection gap **widens** 0.0770 → 0.0876, and **37 of 73 ligands cannot move at any depth**. 035's "the selected rate decays with depth" is **Type-I-specific** — Type II rises. Floors for this population, 2×10⁶ draws: **+0.01383** (d20) / **+0.01367** (d40) |

**Read 007 before proposing any new feature.** It is the arithmetic that makes "promising"
a meaningless word here.

## One method correction, from 036, that changes how 030 and 032 are applied

**The answer-recognition gate does not apply to within-ligand comparators.** The shipped
selector — +0.0395 here, +0.0357 across 81 held-out P450 proteins — puts the query's own
**crystal at the 34th percentile** of its own 20 predictions, binomial p = **0.0012**:
`FINDING_032`'s exact kill signature, on the one feature in this repo that works. The
co-folders share CYP3A4's 30° error (024), so the consensus is displaced from the truth
while remaining informative about the *ordering* of the predictions.

Apply R1 to a **prior** — a claim that some external record says where the ligand goes,
which is what 030 and 032 both were; **both verdicts stand unchanged**. Do **not** apply it
to a **comparator** — a claim that these predictions can be ranked against each other.
`FINDING_032`'s **term oracle** rung is unaffected and is the stronger instrument: it is
what closed 036 before a single candidate was scored.

## Engine and venue facts, all measured the hard way

| # | fact |
|---|---|
| 004 | sampling still pays — the oracle keeps climbing, and the selector tracks it at +0.0125 per doubling |
| 005 | a second engine does not decorrelate: Chai ρ=+0.45, Protenix ρ=+0.60 against Boltz |
| 008 | the P450 superfamily replicates the coordination thesis; a p5–p95 window is **not** an acceptance test and was discarding 10% of true coordination |
| 009 | **`diffusion_samples` does not sample the ligand on OpenProtein** — for any engine. Only replicate jobs do, and replicates are not automatically distinct either. **035 bounds this to OpenProtein's wrapper**: native `boltz predict --diffusion_samples 20` gives 20 of 20 distinct poses in ONE job |
| 035 | **Explorer is a working Boltz-2 venue** — offline, MSA staged on a login node, `--no_kernels`, ~6.6 min/ligand at 20 samples on a V100. `docs/RUNBOOK_explorer_boltz.md`. A job can report **FAILED with all its work complete** (`find \| head` under `pipefail`) — the mirror of trap 2 |
| 013 | a union pool adds **+0.0375 of oracle that selection cannot reach** — keep a second engine as *reference*, never as a pool member |
| 016 | **the sampler sweep is a renewable REFERENCE** — `num_recycles`/`num_steps` diversify at 1.2% catastrophic; worthless as a pool expansion (−0.0038, FINDING 013 again). As a reference it is **equal** to esmfold2 in quality (a +0.0078 edge at n=80 reversed to −0.0028 at n=428) but references **59 pairs esmfold2 cannot**, and can be regenerated for any target |
| 034 | **replicate jobs are not a purchase order for poses** — distinct poses track submission WAVES. 30 replicates -> 10.79 distinct; replicates 16-29 returned 196 of 196 complexes byte-identical to replicate 15; a later wave with the same payload returned a *new* pose shared exactly by all four of its jobs |
| 015 | **protenix_v2 is deterministic too** — a 12→24 doubling moved the oracle on 0 of 489 pairs. Nominal depth 12 is real depth ~4, and the POOL was never deduped, only the reference |

## The three traps that cost the most

1. **Counting artifacts instead of independent opinions.** 20 identical models read as 20 poses; 39% of "replicates" were duplicates; RF3 half-deterministic; Protenix-v1 fully deterministic — and now protenix_v2 too, where a 5,892-pose doubling turned out to be one file repeated. Count *distinct poses*, never jobs — and apply it to the POOL, not only the reference, which is the half that went unchecked. (009, 011, 015)
2. **Trusting a status or a passing check.** Three engines were written off on a misread status or an unread `failure_message`; a readiness check reported 11/11 while the submission was unbuildable. Exercise the thing, do not stat it. (see RUNBOOK "Do NOT")
3. **Reading ρ as selection value.** Rank correlation and top-1 selection moved in *opposite* directions twice. Judge on the metric that will actually be used. (011, 012)

## Superseded or retracted, kept on purpose

- `-zm - zx` looked best at n=63 (+0.0305) and fell to +0.0183 at full depth — **retracted**.
- "Checkpoint diversity beats replicate count" rested on one data point — **retracted**; adding esmfold2 at matched depth halves the gain.
- 009's scope widened twice: Protenix-only → all OpenProtein engines; two working engines → six.
- **"The reference must be better than the pool"** — proposed to explain 012's outlier
  failures, then **contradicted by CYP3A4**, whose reference averages a third of its pool's
  quality and works anyway. Two successor hypotheses (reference self-consistency, reference
  oracle) also failed. P20815 remains an unexplained failure; the thread is stopped until
  20+ proteins have known outcomes.
- The "shallow pools" and "reference depth" caveats on 012 were **tested and retired** rather than repeated.

Each is left next to the reasoning that produced it, because the reasoning is what repeats.
