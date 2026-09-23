# FINDING 039 — the five noise floors are three floors, and FINDING 025's mechanism is dead while its measurement stands

**Date:** 2026-09-23 · **Status:** record repair, zero new inference, CPU-only, one sitting ·
**Verdict:** **one published claim flips sign** (FINDING 036 §2e / §5), **one published
mechanism is withdrawn** (FINDING 025's catastrophe attribution) while the measurement it
explained reproduces to four decimals, and **two of the three cheaper conflicts turn out to
be different quantities rather than disagreements.**

This finding produces no new science. It exists because `DROP_DAY_PLAYBOOK.md` §10 logged
ten contradictions in the published record and logging is not resolving. Two of them were
decision-relevant. Both are now settled, from data already on disk.

---

## Controls, first

Nothing below is trusted until the shipped board reproduces and FINDING 021's numbering fix
is shown to be in force on each set that is touched.

| control | result |
|---|---|
| **C-XENG** — `argmin(xeng)` on all 87, pool A | selected **0.6164**, oracle **0.6975**, random **0.5769**, gain **+0.0395** — the shipped board, to four decimals |
| **C-NUM (pool A)** — FINDING 020's signature was 28 of 85 pairs scoring *exactly* 0.0 | **4** of 1,740 pool rows score exactly 0.0; `mapped` is true on 1,740/1,740 |
| **C-LDDT (arm4 holdout)** — the LDDT column in the feature file vs the post-fix score file | **0 mismatches in 425 poses.** The +0.0408 was always measured on post-fix scores |
| **C-NUM (arm4 holdout)** — FINDING 021's offset control | **61 of 85** pairs need a non-zero residue offset; minimum residue-name identity at the chosen offset **0.809**; **0** pairs below the 0.80 bar |
| **C-BOARD** — the three n=14 boards | pool B d20 `0.6499/0.5270/0.4597`, pool B d40 `0.6595/0.5386/0.4587`, crystal-side d20 `0.6488/0.5055/0.4644` — each matches its finding |

Every filter and its count is printed by the three scripts named under *Reproduce*.

---

## 1. C1 — the n = 14 noise floor. **Five numbers, three quantities, one estimator.**

### 1a. What the statistic is

Every one of the five quotations computes the same thing, under two different names
("a random feature's `argmin`", "random pose selection" — they are identical, because the
`argmin` of a random feature over a ligand's poses *is* a uniform draw from those poses):

> **gain** = mean over ligands of ( LDDT-PLI of **one pose drawn uniformly at random** from
> that ligand's pool ) − mean over ligands of ( that ligand's pool mean ).
> The **floor** is the p95 (and p99) of that gain over independent draws.

So the estimator was never in dispute. What differs is the **population** it is evaluated on
and the **Monte-Carlo stream**. There are exactly three populations in play:

| | population | n | depth | random baseline |
|---|---|---|---|---|
| **D1** | **crystal-side** Type I (`crystal_mode == type_I`), pool A | 14 | 20 | **0.4644** |
| **D2** | **prediction-side** Type I (`pred_fe_donor_median > 2.6 Å`), pool A | 14 | 20 | **0.4597** |
| **D3** | prediction-side Type I, **augmented** pool B | 14 | **40** | **0.4587** |

D1 and D2 differ by **one ligand** — `D0R` in, `QDY` out — which is FINDING 034's own
headline about the label swap, and it moves the random baseline by 0.0047. D2 and D3 differ
by the 20 Explorer poses FINDING 035 bought. They are **not** the same quantity and never were.

### 1b. The authoritative floors

2,000,000 draws each, chosen so the Monte-Carlo error is **twenty times smaller** than the
spread being explained.

| | **p95** | SE of that p95 | **p99** | sd of the gain |
|---|---|---|---|---|
| **D1** crystal-side Type I, depth 20 | **+0.04306** | ±0.00004 | +0.06160 | 0.0273 |
| **D2** prediction-side Type I, depth 20 | **+0.04490** | ±0.00004 | +0.06398 | 0.0281 |
| **D3** prediction-side Type I, depth 40 | **+0.04400** | ±0.00004 | +0.06377 | 0.0274 |

### 1c. Where each published number came from, and how far off it was

A **4,000-draw** p95 at this n has a standard error of **±0.0009**, so two independent
4,000-draw runs of the *same* quantity differ by ±0.0013 at one sd of their difference and
±0.0026 at two. That is the entire residual spread.

| quoted | in | population | authoritative | deviation, in sd of a 4,000-draw estimate |
|---|---|---|---|---|
| **+0.0433** | FINDING 033 | **D1** | +0.04306 | **+0.27 sd** — correct |
| **+0.0431** | FINDING 034 (stratum) | **D2** | +0.04490 | **−1.96 sd** — low |
| **+0.0432** | FINDING 034 (depth curve, rung 20) | **D2** | +0.04490 | **−1.85 sd** — low |
| **+0.0453** | FINDING 035 (`noise_floor_in_stratum`) | **D2** | +0.04490 | **+0.44 sd** — correct |
| **+0.0456** | FINDING 035 (`noise_floor_augmented`) | **D3** | +0.04400 | **+1.64 sd** — high |
| **+0.0435** | FINDING 036 (pool B depth 40) | **D3** | +0.04400 | **−0.51 sd** — correct |

**It is Monte Carlo, not a code difference, and that is shown rather than assumed.**
FINDING 034's `stratum()` null was replayed from its own seed, `default_rng(SEED + 11)` with
`SEED = 20260922`, and returns **p95 +0.0431 / p99 +0.0631 bit-for-bit** — the published
pair. The estimator is right; the stream was a low draw. The same replay confirms the
depth-curve rung at +0.0432.

**FINDING 036's "reproduces FINDING 035 (+0.0453)" was wrong twice over.** 036's own
+0.0435 is on **D3** and 035's +0.0453 is on **D2** — different pools — and 036 then compared
its number to 035's *other* number, the +0.0453 rather than the +0.0456 that shares its
population. The right statement is that 036's +0.0435 and 035's +0.0456 are two 4,000-draw
estimates of D3 = +0.0440, 0.5 sd low and 1.6 sd high respectively.

### 1d. **What changes sign**

| claim | as published | against the authoritative floor | verdict |
|---|---|---|---|
| **036 §2e / §5**: pool B's perfect top-2 ceiling **+0.0447** "lands **on** the floor … nothing is provable" | vs +0.0435 → already above it, by +0.0012, which 036 did not notice | vs **D3 +0.04400** → **above by +0.00070**; exact one-sided p against its own null, 2 × 10⁶ draws, **p = 0.0475** | **FLIPS.** The ceiling **clears** p95, by a hair |
| 033: Type I gain **+0.0411** does not clear its floor, p = 0.0595 | vs +0.0433 | vs **D1 +0.04306**, exact **p = 0.0581** | **stands** |
| 034: Type I gain **+0.0673** clears its floor, p = 0.0063 | vs +0.0431 | vs **D2 +0.04490**, exact **p = 0.0073** | **stands** |
| 035: augmented gain **+0.0799** clears its floor | vs +0.0456 | vs **D3 +0.04400**, exact **p = 0.0020** | **stands**, more comfortably |
| 036 §4g: best pool-B candidate **+0.0392** is below the floor | vs +0.0435 | vs **D3 +0.04400** | **stands** |
| 036 stage-1 stop: pool A ceiling **+0.0129** vs pooled floor **+0.0134** | n = 87, not in dispute | untouched | **stands** |

**The verdict of FINDING 036 does not change.** Its stage-1 stop is the n = 87 measurement
(+0.0129 against +0.0134), which nothing here touches; its n = 14 stratum was explicitly
"quoted, not concluded from"; and the realisable candidates on pool B (best **+0.0392**)
remain below the floor. What changes is a **statement about a ceiling**: the perfect,
oracular top-2 tie-break on the 14-ligand depth-40 stratum is **marginally above** its floor
at p = 0.0475, not "on" it. At n = 14, with a perfect-knowledge ceiling and no candidate
that reaches it, that is not a prize — but the record should say what it is.

### 1e. The operating rule this replaces

`DROP_DAY_PLAYBOOK.md` §10 C1 said *"use the largest when judging a candidate"* — i.e.
+0.0456. Superseded. **Use the floor for the population you are on:** +0.0449 for a 14-ligand
depth-20 stratum, +0.0440 for a 14-ligand depth-40 stratum, +0.0137 pooled at n = 87. And the
playbook's second half was always the real rule and still is: **recompute on the actual pool**
— at 4,000 draws, which is what every script here uses, ±0.0009 of slop is built in, so
anything inside ±0.003 of a floor is a coin flip about a coin flip. **Use ≥ 200,000 draws when
a decision turns on it**; on this data 2 × 10⁶ draws cost eight seconds.

---

## 2. C2 — FINDING 025's arm4 holdout. **The measurement stands. The mechanism does not.**

FINDING 025 reported `fe_centroid_dist` ("closer to the iron is better") scoring **+0.0408**
on the 85-pair arm4 P450 holdout against +0.0013 on the CYP3A4 pool, and explained the
difference by saying the holdout is a pool "where **39 of 85 pairs are catastrophes**" —
FINDING 012's catastrophe-detector scaling reproducing on a new feature. FINDING 021
re-measured that catastrophe count at 1 of 84 after the residue-numbering fix.

### 2a. The measurement reproduces exactly

The per-pose feature file was recovered from Explorer (`zexp/depth_terms_holdout.csv`,
35 KB — no inference, no re-scoring) and `depth_replication.evaluate` re-run verbatim:

| feature, "lower is better" | recomputed | published |
|---|---|---|
| `fg_depth` | −0.0051, ρ 0.023, p 0.6705 | −0.0051, ρ 0.023, p 0.6705 |
| **`fe_centroid_dist`** | **+0.0408**, ρ 0.138, p 0.0005 | **+0.0408**, ρ 0.138, p 0.000 |
| `fe_min_dist` | +0.0291, ρ 0.133, p 0.0075 | +0.0291, ρ 0.133, p 0.006 |

C-LDDT confirms the LDDT column in that file is the **post-fix** score, pose for pose, on all
425. **The +0.0408 was never contaminated by the numbering bug.**

### 2b. The catastrophe census, post-fix

| definition | count |
|---|---|
| poses scoring < 0.1 | **7 of 425 = 1.6%** |
| pairs with **any** pose < 0.1 | **2 of 85** |
| pairs with **all** poses < 0.1 (FINDING 021's definition) | **1 of 85** |
| pairs with mean < 0.1 | **1 of 85** |
| **asserted in FINDING 025** | **39 of 85** |

There is no reading of the post-fix holdout on which 39 pairs are catastrophes. The
population the explanation rests on does not exist, exactly as the playbook suspected.

### 2c. The mechanism is not merely unverified — it is **refuted**

Three independent ways, all decisive:

| test | result |
|---|---|
| **Ceiling of a perfect catastrophe-avoider** — a rule that never picks a pose < 0.1 and otherwise picks at random | **+0.0002.** That is **0.5%** of +0.0408. Avoiding every pose below 0.2 is worth +0.0076; below 0.3, +0.0351 |
| **Delete the catastrophes and re-run** — drop all 7 sub-0.1 poses, rescore the feature on the 418 that remain | **+0.0409** (ρ 0.136, p 0.0010). The gain does not move by one part in four hundred |
| **Decompose by pair** | the 2 pairs containing a catastrophic pose contribute **+0.0011** of the +0.0408; the 83 clean pairs contribute **+0.0396** |

### 2d. The obvious replacement fails too

If not catastrophes, then FINDING 012's broader claim — the gain scales with **pool
uncertainty**. It does not, and it points the wrong way:

| | headroom (oracle − random) | poses below 0.5 | `fe_centroid_dist` gain |
|---|---|---|---|
| arm4 holdout | **+0.0652** | 14.6% (62/425) | **+0.0408** |
| CYP3A4 pool | **+0.1174** | **34.2%** (595/1740) | **+0.0013** |

The CYP3A4 pool has **nearly twice the headroom and more than twice the bad poses**, and the
same feature is worth nothing there. Whatever separates the two pools, it is not pool quality
in the direction FINDING 012 describes.

One descriptive difference is on the record for whoever picks this up: within a pair, the
holdout's `fe_centroid_dist` varies by **0.51 Å** (sd 0.21) while CYP3A4's varies by
**1.77 Å** (sd 0.46). The holdout's poses are nearly all near-correct (76% score above 0.7)
and small distance differences track quality; CYP3A4's scatter widely in distance and it
tracks nothing. **That is a description, not a mechanism, and it is not offered as one.**

### 2e. Verdict

Of the three options the audit named — both stand, measurement only, neither —
**the evidence supports the second: the measurement stands and the mechanism does not.**
`fe_centroid_dist` really does select on the arm4 holdout, at +0.0408 with p = 0.0005 against
its own null, and **why it does is now unexplained.** It is not catastrophe detection and it
is not pool uncertainty. Nothing in this repo builds on it, and nothing should until someone
explains it.

---

## 3. C8 — FINDING 035's two 20-pose baselines, and the pool-30 peak

### 3a. Two estimands, not two estimates of one

| number | what it is |
|---|---|
| **0.5270** | `argmin(xeng)` over **the 20 Modal poses that exist**. Deterministic; no draw, no variance. Recomputed here: **0.5270** |
| **0.5291** | **E**[ `argmin(xeng)` over 20 poses **drawn from the 40-pose union** ]. Half the poses are Explorer poses. Exact closed-form value: **0.5293** |

Different populations, so both are right and neither is "the" baseline. **The primary
endpoint's baseline is 0.5270** — it is the pool that existed before the purchase, which is
what a purchase is measured against. 0.5291 belongs only to the union curve, where it is the
correct rung. 035 labelled both correctly in context and never says which is "the"
conversion; the playbook's C8 is a reading hazard, not an error.

### 3b. The peak at pool 30 is **real**, and it is **one ligand**

Both of 035's curves are 256-draw Monte Carlo, and both are exactly computable in closed
form, because `argmin` over a random subset has a combinatorial distribution. Computed
exactly — **zero sampling error**:

| pool | 20 | 24 | 28 | 30 | **33** | 36 | 38 | **40** |
|---|---|---|---|---|---|---|---|---|
| exact selected (augment curve) | 0.5270 | 0.5349 | 0.5400 | 0.5415 | **0.5424** | 0.5418 | 0.5406 | **0.5386** |
| 035 published (256 draws) | 0.5270 | 0.5354 | 0.5393 | 0.5411 | | | | 0.5386 |

The union curve does the same thing: exact peak **0.5410 at depth 34**, falling to 0.5386 at
40. **The non-monotonicity is not Monte-Carlo noise** — it is an exact property of this pool,
and the published numbers are all within ±0.0007 of their exact values (the 256-draw SE is
±0.0012). The oracle, as expected, is monotone throughout.

**But it is one ligand.** Per-ligand, exact, pool 30 minus pool 40:

| ligand | pool 30 | pool 40 | Δ |
|---|---|---|---|
| **08J** | 0.3085 | **0.1364** | **+0.1721** |
| CFF | 0.8352 | 0.9547 | −0.1195 |
| YNV | 0.6271 | 0.6448 | −0.0177 |
| PG0 | 0.4765 | 0.4706 | +0.0059 |
| the other **10** | — | — | **0.0000** |
| **mean** | 0.5415 | 0.5386 | **+0.0029** |

Ten of fourteen ligands have **no** new pose with a better `xeng` than their best existing
pose, so depth cannot move them at all. The decline is 08J net of CFF — the same two swing
ligands 036 §4g traced by hand. 08J has exactly one new pose that beats its incumbent on
`xeng` and scores **0.1364** where the incumbent scored 0.4806; at pool 30 there is a 50%
chance that pose is not drawn, at pool 40 it is always there and always wins.

**What this licenses.** It is a real, exactly-computed inversion, and it is effectively
n = 1. It does **not** establish "selected score turns down past ~33 poses" as a law. It
establishes that **at n = 14 a single adversarial pose can invert the selected curve while
the oracle rises monotonically** — which is 035's own thesis (the oracle law survives, the
selection rate does not) reduced to its smallest possible instance. It is one more reason
that the depth question cannot be settled on this stratum, and it is not a reason to stop
buying depth.

---

## 4. C5 — aromatase. **Two statistics on the same six pairs.**

FINDING 022 prints P11511 at **0.863** over 6 pairs; FINDING 024 prints **0.631**. Computed
from `scores_arm4_mix_base.json` and `arm4_mix.csv`:

| target | n pairs | **mean of `lddt_sample0`** | FINDING 022 | **mean of all 5 samples** | FINDING 024 |
|---|---|---|---|---|---|
| P08684 CYP3A4 | 15 | **0.5550** | 0.555 | 0.5526 | — |
| **P11511 aromatase** | **6** | **0.8629** | **0.863** | **0.6313** | **0.631** |
| P10614 CYP51A1 | 5 | **0.9093** | 0.909 | 0.9076 | — |
| Q2IU02 | 27 | **0.9306** | 0.931 | **0.9241** | 0.924 |
| P20815 CYP3A5 | 2 | 0.4758 | — | **0.5059** | 0.506 |
| every other P450 | 70 | **0.8648** | 0.8648 | 0.8466 | — |

**Settled, and it is not a disagreement.** FINDING 022 reports the **model's own rank-0
sample**; FINDING 024 reports the **mean over all five samples**. Every quoted figure in both
findings reproduces to three decimals under its own rule, on five targets plus both aggregate
rows. The **"7 entries"** on 024's aromatase row is a third denominator again — the count of
P11511 **crystals in the 406-entry cavity set**, not pairs scored — so that one line carries
a cavity from 7 crystals beside a score from 6 pairs.

One substantive thing falls out. Aromatase is where the two statistics diverge most by a wide
margin: **0.232**, against ≤ 0.007 for CYP3A4, CYP51A1 and Q2IU02. Boltz-2 puts a very good
aromatase pose at rank 0 and much worse ones behind it. **Prefer the all-5 mean** when
comparing targets — rank-0 is the engine's own confidence ordering, and FINDING 001 is that
this ordering does not rank poses. 022's table is not wrong, but it reads the engine's
self-assessment as if it were quality.

---

## 5. C3 — FINDING 028's licensed holo template

No measurement needed; the addenda already did it. A dated pointer has been added at the top
of 028's body so the superseded "LICENSED, and it reverses a call in FINDING 024" section
cannot be acted on by a reader who stops before line 514.

---

## What is still open

- **Why `fe_centroid_dist` selects on the arm4 holdout at all.** Not catastrophes, not pool
  uncertainty. Unexplained, and flagged as such in 025.
- **Whether the pool-33 peak generalises.** It cannot be answered on this stratum. What would
  settle it: the same exact closed-form curve on a second matched-depth stratum of ≥ 30
  ligands, which needs a depth purchase that 035 already priced as not worth making. Until
  then it is a curiosity about 08J.
- **Playbook §10 C4, C6, C7, C9** — not touched here. C6 and C9 are reconciliations the text
  already contains; C4 and C7 are supersessions, not conflicts.
- **The 4,000-draw convention itself.** Every null in FINDINGs 033–037 carries ±0.0009 of
  Monte-Carlo error at n = 14 and ±0.0003 at n = 87. No published verdict turns on it now,
  but the next one might.

---

## Reproduce

```
python scripts/ops/corrections_floor.py       # C1  -> corrections_c1_floor.json
python scripts/ops/corrections_arm4.py        # C2  -> corrections_c2_arm4.json
python scripts/ops/corrections_depth.py       # C8  -> corrections_c3_depth.json
python scripts/ops/corrections_aromatase.py   # C5  -> corrections_c5_aromatase.json
```

| output | file |
|---|---|
| three floors, their SEs, every quoted value placed | `data/processed/corrections_c1_floor.json` |
| exact p-values for the four gain-vs-floor verdicts | `data/processed/corrections_c1_pvalues.json` |
| arm4 reproduction, catastrophe census, avoidance ceiling, decomposition | `data/processed/corrections_c2_arm4.json` |
| exact augment and union depth curves, both baselines | `data/processed/corrections_c3_depth.json` |
| per-target sample-0 vs all-5 means | `data/processed/corrections_c5_aromatase.json` |

`corrections_arm4.py` needs `depth_terms_holdout.csv` and `depth_terms_cyp3a4.csv`, pulled
once from `explorer:/scratch/shenoy.am/zexp/` (276 KB total) into `C:\Temp\cyp_corrections`.
No inference, no GPU, nothing written to `D:` but the JSON outputs and these documents.

---

## The method note worth keeping

**A floor is a number *about a population*, not a constant.** Five quotations of "the n = 14
floor" turned out to be three populations differing by one ligand and by twenty poses, plus
±0.0009 of Monte-Carlo error that nobody had sized. The repo already knew the floor is
n-dependent — every script recomputes it "inside the stratum" for exactly that reason — and
still cross-quoted across populations, because the quantity had a short name and the name did
not carry the population. **Name the population in the number, or the number will be quoted
where it does not apply.**

And the second one, which is FINDING 021's lesson in a new coordinate: **a retraction has to
be chased into every document that cited it.** FINDING 021 retracted the 39 catastrophes in
full and FINDING 025 went on explaining a live result with them four findings later, because
nothing walks the citations. The measurement survived; the explanation was load-bearing and
false for three days.
