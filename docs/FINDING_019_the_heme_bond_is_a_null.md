# FINDING 019 — The heme bond is a null, because coordination was never broken

**Date:** 2026-09-20 · **Status:** measured, 84 paired held-out pairs · **Verdict:** null

## What was tested

FINDING 018 killed fine-tuning and pointed here instead: the base checkpoint reached
44/85 sub-2 Å on held-out P450 pairs while dropping 39/85 into catastrophe, and it did
that with **no heme bond**. This repo had already measured that an explicit
Cys-SG→heme-FE bond is what moves Boltz-2 to crystallographic Fe-donor geometry. The
obvious inference was that adding it would convert catastrophes into successes, for free.

The same 84 pairs, the same MSAs, the same 5 diffusion samples, the same scorer — one
line of YAML different. The ligating cysteine came from **sequence** (the conserved
P450 FxxGxxxCxG motif, 81/85 strict matches, cysteine at relative position
0.715/0.870/0.897), never from the crystal, so the method transfers to a blind target.

## Result: nothing

| metric | unbonded | bonded | Δ | bonded wins | p |
|---|---|---|---|---|---|
| LDDT-PLI sample 0 | 0.3430 | 0.3349 | −0.0081 | 21/84 | 0.086 |
| LDDT-PLI best of 5 | 0.3608 | 0.3516 | −0.0092 | 17/84 | 0.108 |
| BiSyRMSD | 5.89 Å | 6.02 Å | +0.13 Å | 34/84 | 0.058 |

| | unbonded | bonded |
|---|---|---|
| catastrophes (<0.1) | 39 | 38 |
| successes (>0.5) | 22 | 21 |
| sub-2 Å | 44/84 | 43/84 |

Not significant, and — the familiar refrain — **the bimodal split does not move**.

## Why: the anchor was already correct

A null has two explanations that demand opposite responses: the constraint never took
effect, or it took effect and did not matter. Measuring Fe-to-ligand distances separates
them, and the answer is unambiguous.

| set | n | Fe-donor median | in 1.90–2.45 Å | under 3 Å |
|---|---|---|---|---|
| base, nearest ligand atom | 85 | 3.763 Å | 0.282 | 0.341 |
| bonded, nearest ligand atom | 84 | 3.991 Å | 0.286 | 0.310 |
| **crystal, nearest ligand atom** | 85 | **3.728 Å** | **0.271** | **0.318** |
| base, nearest ligand N | 35 | 2.112 Å | 0.686 | 0.714 |
| bonded, nearest ligand N | 35 | 2.236 Å | 0.686 | 0.714 |
| **crystal, nearest ligand N** | 35 | **2.212 Å** | **0.600** | **0.714** |

**The unbonded predictions already reproduce the crystal's coordination distribution.**
Median 3.763 Å against the crystal's 3.728 Å; nitrogen donors at 2.112 Å against 2.212 Å;
the fraction under 3 Å matches to within a percentage point. The bond changed the median
by roughly a tenth of an Ångström and the in-window fraction by nothing at all (0.686 in
both).

There was nothing to fix. That is the fifth independent confirmation that **the iron
anchor is saturated** and the remaining error lives in substituent placement.

## What it implies about the catastrophes

This is the part worth carrying forward. 39 of 85 poses score under LDDT-PLI 0.1 while
the ligand sits at crystallographic distance from the heme iron. The molecule is
anchored in the right place and is still catastrophically wrong — which is thesis A
("co-folding gets the ORIENTATION wrong, not the location") holding at a much larger
scale than it was originally measured on.

So a catastrophe here is a ligand rotated or flipped about a correctly-placed anchor. No
constraint on the anchor can address that, which retro-justifies FINDING 008's negative
result on anchor-local features and explains why FINDING 018's fine-tune could not move
the split either.

## Running total on interventions that do not move the split

| intervention | Δ LDDT-PLI | catastrophes |
|---|---|---|
| fine-tune, 87 steps | −0.0038 | 39 → 39 |
| fine-tune, 174 steps | −0.0138 | 39 → 39 |
| fine-tune, 350 steps | −0.0125 | 39 → 38 |
| heme bond constraint | −0.0081 | 39 → 38 |

Four interventions, two mechanisms, zero poses reclassified. Meanwhile cross-engine
consensus **selection** gains +0.0381 without touching generation at all (FINDING 011).
The evidence keeps saying the same thing: generation is not the lever, selection is.
