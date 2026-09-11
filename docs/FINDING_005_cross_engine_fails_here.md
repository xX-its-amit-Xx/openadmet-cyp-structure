# Finding 005 — a second engine adds nothing, because its errors are correlated

**Date:** 2026-09-11 · **Runs:** `val87b` (Boltz-2, 20 samples) + `chai87` (Chai-1, 10 samples)
**n = 87 ligands, 2,610 poses** · cost of the Chai pool: **$22.39**

FINDING 002 named cross-engine agreement as the most promising untested signal, and
FINDING 003 left it open. It is a null, and four separate tests agree.

---

## 1. Chai-1 is a legitimate engine, not a straw man

| at matched 10 samples | oracle | random pose | coordinated | median Fe–donor |
|---|---|---|---|---|
| Boltz-2 | 0.6739 | 0.5779 | 84.1% | 2.23 Å |
| Chai-1 | 0.6625 | 0.5323 | 80.8% | 2.29 Å |

87 of 87 jobs succeeded, 10 structures each. Chai reaches the crystallographic Fe–donor
distance without any bond constraint, as Boltz does. It is a little behind, not broken.

## 2. Cross-engine agreement does not select

| feature | Δ vs random | p |
|---|---|---|
| distance to nearest other-engine pose | +0.0007 | 0.89 |
| mean distance to the other engine | −0.0054 | 0.82 |
| count of other-engine poses within 2 Å | +0.0066 | 0.70 |

And adding it to the FINDING 003 selector **makes it worse**: +0.0279 → +0.0241 (w=0.25)
→ +0.0191 (w=0.5). It is not neutral, it is dilution.

## 3. Why: the two engines fail on the same ligands

Per-ligand oracle, Boltz vs Chai: **Spearman ρ = +0.451, Pearson r = +0.513** (p ≈ 1e-5).
Of the 20 ligands Boltz handles worst, **11 are also in Chai's worst 20** — chance would
be 4.6.

Agreement between two models with correlated errors measures **shared bias**, not
correctness. That is the whole result.

## 4. The union raises the ceiling and lowers the score

| pool | oracle | selected by FINDING 003 |
|---|---|---|
| Boltz only | 0.6975 | **0.6048** |
| Chai only | 0.6625 | 0.5430 |
| **union** | **0.7423** | 0.5875 |

The union's ceiling is **+0.0448** above Boltz, nearly twice what doubling Boltz samples
buys, and Chai's best pose beats Boltz's on **32 of 87 ligands**. Yet selecting over the
union scores **0.017 lower** than Boltz alone, or 0.045 lower without per-engine z-scoring.

**Pool diversity only pays if the selector is good enough to exploit it.** Chai's poses are
worse on average, so adding them mostly gives a weak selector more chances to err.

## 5. Tail rescue does not transfer

The PXR campaign's winning rule was to add a second engine only for the least-confident
ligands: Protenix regressed the board 0.5551 → 0.5241 as a pool member but gained +0.017
as an 8-ligand tail swap. Swapping Chai into the least-confident tail here:

| N swapped | 0 | 4 | 8 | 12 | 20 | 35 | 50 | 87 |
|---|---|---|---|---|---|---|---|---|
| score | 0.6048 | 0.6074 | 0.6022 | 0.5947 | 0.5937 | 0.5814 | 0.5634 | 0.5430 |

No interior optimum. N=4 is +0.0026 at p=0.27, i.e. noise, and everything past it declines
monotonically.

**This unifies the inherited rule rather than contradicting it.** Tail rescue works when the
second engine's errors are *decorrelated* from the first's — that is the mechanism, and PXR
simply never had to state it because it happened to hold there. Here ρ = +0.45, so where
Boltz is unsure Chai is unsure about the same thing, and there is nothing to rescue with.

---

## What to do

- **Stop spending on Chai.** $22.39 bought a clean negative, which is worth having, but
  there is no second increment to buy.
- **Spend on Boltz samples instead.** FINDING 004: +0.0125 per doubling, and Boltz is about
  5× cheaper per sample ($0.156/job vs $0.411).
- **A third engine is only worth trying if its errors would be decorrelated.** Test that
  cheaply — per-ligand oracle correlation on a 20-ligand pilot — before buying a full pool.
- **The selector remains the binding constraint.** It cannot exploit a ceiling that is
  already 0.7423.
