# W001 — Multi-ligand occupancy in the CYP3A4 validation set: real, rare, and not yet a finding

**Date:** 2026-09-20 · **Status:** measured, n too small to conclude · **Verdict:** hypothesis with a defined test, not a result

## Why we looked

The evolutionary decomposition turned up that CYP3A4 is **functionally multi-occupancy**:
a single testosterone in the site makes essentially no product — it uncouples — and
maximum rate requires **two**, with coupling improving again at three
(doi:10.1074/jbc.m609589200). That is now visible structurally: 3 and 6 caffeines in
8SO1/8SO2, and 4 DHEA-S in CYP3A7's 8GK3.

Our pipeline predicts **one** ligand per structure. If the crystals routinely hold
several, prediction and scoring are both answering a different question from the one the
crystal asks. So we counted, rather than assuming either way.

## What the crystals contain

Copies of the query ligand within a single chain, across all 87 CYP3A4 validation entries
(0 missing):

| copies in one chain | entries |
|---|---|
| 1 | 82 |
| 2 | 4 |
| 3 | 1 |

**5.7% multi-copy, maximum 3.** The five: `2V0M_KLN` (ketoconazole x2 — the classic
cooperativity structure), `8SO1_CFF` (caffeine x3 — exactly as the literature describes),
`4K9T_1RD`, `4K9U_5AW`, and `9YK4_PG4` (PG4 is a PEG fragment, i.e. a cryoprotectant, not
a real ligand case).

So multi-occupancy is real and it is rare. It is not a systematic scoring problem.

## Do they score worse? Suggestive, and contradicted

Per-ligand pool oracle, rank out of 87 (1 = best):

| ligand | copies | oracle | rank |
|---|---|---|---|
| 5AW | 2 | 0.2164 | **87** |
| 1RD | 2 | 0.4013 | **82** |
| KLN | 2 | 0.4684 | **79** |
| CFF | 3 | **0.9081** | **4** |
| PG4 | 2 | 0.8185 | 23 |

Three of the four real multi-copy ligands sit in the bottom ten of eighty-seven. That is
the pattern the mechanism predicts — we place one molecule where the crystal has two, so
even the best pose in the pool is penalised.

**And then caffeine, with the most copies of any of them, ranks fourth.** At n = 4 that is
an anecdote, not an effect, and the single strongest counterexample is the entry with the
highest occupancy. Reporting the first three without the fourth would be selection.

## The test that would settle it

Co-fold the five with the crystallographic number of ligand copies rather than one, score
against the same references, and compare paired. If 5AW, 1RD and KLN move substantially
while CFF does not, the mechanism is occupancy. If nothing moves, they are simply hard
ligands and the rank pattern was coincidence.

Cheap — five ligands, one inference run. Worth doing, not worth believing first.

## What this does NOT license

It does not license predicting multiple copies generally. 82 of 87 entries hold exactly
one, and adding copies there would inject a second ligand into a site that has none,
which is a reliable way to make good poses worse.
