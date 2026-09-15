# Findings index

Sixteen findings, most of them negative. Read them in this order if you are arriving cold;
the numbering is chronological, not logical.

## Start here

| # | one line | status |
|---|---|---|
| **001** | Selection, not generation, is the bottleneck — pool oracle 0.6975 vs selection 0.5706, and Boltz's own confidence ranks poses *worse* than random | stands |
| **011** | **Cross-engine agreement selects: +0.0381**, beating the incumbent's +0.0265. No fitted parameters. This is what ships | stands |
| **012** | It **generalises** — positive on 13 of 16 held-out P450 proteins — and it is a **catastrophe detector**, so its payoff is predictable from pool quality | stands |
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
| 016 | **the sampler sweep is a renewable REFERENCE** — `num_recycles`/`num_steps` give 5 distinct poses from 5 settings at 1.2% catastrophic; worthless as a pool expansion (−0.0038, FINDING 013 again) but references 100 pairs where esmfold2 manages 80, and can be generated for any target on demand |
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
