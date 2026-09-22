# Findings index

Nineteen findings, most of them negative. Read them in this order if you are arriving cold;
the numbering is chronological, not logical.

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

**Read 007 before proposing any new feature.** It is the arithmetic that makes "promising"
a meaningless word here.

## Engine and venue facts, all measured the hard way

| # | fact |
|---|---|
| 004 | sampling still pays — the oracle keeps climbing, and the selector tracks it at +0.0125 per doubling |
| 005 | a second engine does not decorrelate: Chai ρ=+0.45, Protenix ρ=+0.60 against Boltz |
| 008 | the P450 superfamily replicates the coordination thesis; a p5–p95 window is **not** an acceptance test and was discarding 10% of true coordination |
| 009 | **`diffusion_samples` does not sample the ligand on OpenProtein** — for any engine. Only replicate jobs do, and replicates are not automatically distinct either |
| 013 | a union pool adds **+0.0375 of oracle that selection cannot reach** — keep a second engine as *reference*, never as a pool member |
| 016 | **the sampler sweep is a renewable REFERENCE** — `num_recycles`/`num_steps` diversify at 1.2% catastrophic; worthless as a pool expansion (−0.0038, FINDING 013 again). As a reference it is **equal** to esmfold2 in quality (a +0.0078 edge at n=80 reversed to −0.0028 at n=428) but references **59 pairs esmfold2 cannot**, and can be regenerated for any target |
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
