# FINDING 016 — the sampler sweep is a renewable REFERENCE, not a pool expansion

**Status: stands.** Measured 2026-09-15, 100 P450 pairs / 29 proteins / 500 poses.

Context: FINDING 015 killed `--replicates` — both OpenProtein engines now return a
byte-identical pose per input, so pose depth was frozen at what was already on disk, and
the esmfold2 *reference* was stuck at a median of 3 distinct poses against a threshold of
4. This finding is the way out, and it is not the way that was expected.

## The sweep produces real, non-degrading diversity

Varying `num_recycles` and `num_steps` (3×200 default, 1×200, 10×200, 3×50, 3×400):
**5 distinct ligand placements from 5 settings on 100 of 100 pairs**, zero unretrievable.

Scored against crystal ground truth, this is diversity rather than damage — the check
that disqualified single-sequence mode (75% catastrophic, oracle 0.0505):

| setting | mean LDDT-PLI | vs default |
|---|---|---|
| 10×200 | 0.7353 | +0.0079 |
| 3×50 | 0.7307 | +0.0033 |
| 3×200 (default) | 0.7274 | — |
| 3×400 | 0.7251 | −0.0023 |
| **1×200** | 0.6677 | **−0.0596** |

**Catastrophic: 1.2%.** `num_recycles=1` is genuinely degrading and is excluded
everywhere below; the rest match the default. No single setting is meaningfully better —
the value is in the *union*, which is what a diversity lever should look like.

## As a POOL expansion it is worthless — FINDING 013, reproduced

The union of settings raises the oracle **+0.0355** over the default and the selector
captures **none** of it:

```
default 3x200 only : 0.7739
random over settings: 0.7766
SELECTED (xeng)     : 0.7701
oracle              : 0.8094
gain over random    : -0.0065   (null p99 +0.0100, p=0.9003)
gain over DEFAULT   : -0.0038   <- what shipping it would actually buy
```

Set against FINDING 013:

| experiment | oracle added | selection captured |
|---|---|---|
| union of engines (013) | +0.0375 | −0.0017 |
| sampler sweep (this) | +0.0355 | −0.0038 |

Two unrelated mechanisms — a different architecture, and different sampler settings on
the *same* architecture — add almost exactly the same oracle and both yield nothing.
**New diversity enlarges what is achievable without making it findable.** Treat that as
the rule now, not the coincidence it looked like once.

## As a REFERENCE it beats the incumbent and renews it

Matched: same 80 pairs, same pool, only the reference source varies.

| reference | random | selected | gain |
|---|---|---|---|
| esmfold2 (incumbent) | 0.7405 | 0.7764 | +0.0359 |
| **sweep only** | 0.7405 | **0.7842** | **+0.0436** |
| esmfold2 + sweep | 0.7405 | 0.7788 | +0.0383 |

Three things, in decreasing order of confidence:

1. **Coverage — the solid result.** The sweep can reference **100 pairs where esmfold2
   manages 80**, because esmfold2's determinism leaves 20 pairs below the 2-distinct
   minimum. 25% more pairs get a working selector, and this needs no statistical claim.
2. **Renewability — the strategic result.** esmfold2 is frozen and cannot be deepened at
   any price. A sweep can be generated for *any* pair on demand, including test ligands
   nobody has folded yet. The reference source stops being a fixed asset.
3. **⚠️ RETRACTED: "+0.0078 over the incumbent".** It did not replicate. See below.

Combining the two references looked **worse than either alone** here (+0.0383 vs
+0.0436). That did not replicate either — at full scale all three arms are equal.

## What to do with it

- **Do not** add sweep poses to the pool. It costs 4× the jobs to lose 0.0038.
- **Do** use the sweep as the reference source, especially for pairs where the reference
  engine is too shallow to score at all.
- Drop `1×200` from the sweep. Four settings, not five.
- On drop day this is the reference recipe: four sampler settings on the target, which
  needs no second engine and no MSA for the reference at all.

## Cost

500 jobs to close the pool question and open the reference one. The pool arm was decided
by one number — gain over *default*, not over random — because random over a sweep is not
the alternative anyone would ship.


---

## Correction — the quality edge was noise. Coverage is the whole result.

The matched comparison above used 80 pairs. Repeated on **428** — the full set where both
reference sources are usable, same pool, only the reference varying:

| reference | pairs | proteins | selected | gain | within-pair ρ |
|---|---|---|---|---|---|
| esmfold2 | 428 | 81 | 0.7479 | **+0.0270** | −0.170 |
| sweep | 428 | 81 | 0.7451 | **+0.0242** | −0.202 |
| esmfold2 + sweep | 428 | 81 | 0.7474 | +0.0264 | −0.193 |

**sweep − esmfold2 = −0.0028.** At n=80 it was +0.0078. The effect reversed sign, and
both values sit inside this pool's null p99 of ~+0.0100. That is FINDING 007's arithmetic
doing exactly what it was written to do: **under +0.020 is noise**, and a number quoted as
"suggestive" at n=80 was simply noise with a sign.

So the honest statement is that **the sweep and esmfold2 are equivalent as reference
sources.** Not better. Equivalent.

Note also that the sweep has the *better* rank correlation (−0.202 vs −0.170) while
losing on top-1 selection. Rank correlation and selection moving in opposite directions is
the third trap in `docs/README.md`, met for a third time. Judge on the metric that ships.

### What survives, and it is still worth having

- **Coverage.** The sweep references **488 pairs where esmfold2 manages 432**, including
  **59 pairs esmfold2 cannot reference at all**. At full-set scale that is 487 pairs / 87
  proteins in the generalisation test against 430 / 81 — six whole proteins that had no
  working selector now have one. This needs no statistical claim and does not depend on
  the retracted edge.
- **Renewability.** esmfold2 is deterministic and frozen; a sweep can be generated for any
  target, including test ligands nobody has folded. That is what makes it a drop-day
  reference recipe rather than a fixed asset.

The pool-expansion result is untouched: as a pool member the sweep is still −0.0038, and
FINDING 013 still holds twice over.
