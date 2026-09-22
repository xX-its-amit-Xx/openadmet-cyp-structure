# FINDING 026 — splitting the shipped selector into translation and rotation: the term separates, the selector does not

**Date:** 2026-09-22 · **Status:** measured, zero new inference · **Verdict:** no change to
`cypstruct.xengine.select()`. The orientation-only term is the best single unfitted term
ever measured here (+0.0424) and it **does not beat the incumbent in a paired test**
(+0.0041, Wilcoxon p = 0.61, Holm p = 1.00).

This is R1 of FINDING 024, run exactly as pre-registered.

## What was tested

The shipped selector scores a pose by `xeng_score` — the mean symmetric Chamfer distance,
in the pose's own heme frame, to the frozen independent-engine reference poses of the same
ligand. Chamfer mixes **where** the ligand sits with **how it is turned**. FINDING 024
measured the centroid as already correct to 1.07 Å with isotropic scatter and no direction,
and the orientation as wrong by 30.1° (family: 6.1°). The hypothesis was therefore that the
incumbent spends its dynamic range on a coordinate that carries almost no signal.

Four parameter-free terms, all in the **same** heme frame (`xengine.in_heme_frame`,
unmodified), each averaged over that ligand's references:

| term | definition |
|---|---|
| `mix` | `chamfer(pose, ref)` — **the incumbent, the control to beat** |
| `cen` | `‖centroid(pose) − centroid(ref)‖` — translation only |
| `ori` | `chamfer(pose − c_p, ref − c_r)` — orientation (and internal conformer), translation removed |
| `pax` | mean angle between the three principal axes of the two inertia tensors (cheap orientation proxy) |
| `pax1` | the same for the longest axis alone |

**Nothing is fitted.** No weight, no threshold, no held-out fold, in any row below.

**Pool:** the 1,740 poses of `/scratch/shenoy.am/zexp/out/run5_s4{2,3,4,5}`, 87 ligands,
truth from `truth_run5_s4*.csv`, references from
`cyp-finetune/reference_set_cyp3a4.npz` (87/87 ligands, depth min 6 / median 7 / max 12,
all ≥ 4). Scripts `zexp/dec.py` and `zexp/dec2.py`, both under `sbatch`; results in
`zexp/decomp_terms.csv` and `zexp/decomp_results.json`.

**Controls before any result.** The recomputed `mix` matches the shipped `xeng_4seed.csv`
to **max |Δ| = 7.1e-15** over all 1,740 rows, so the control genuinely *is* the incumbent
and not a lookalike. Pose and reference heavy-atom counts agree on **1,740 of 1,740** rows
(median 37), so the centroid term is comparing the same atom set.

**Oracle 0.6931 · random 0.5757 · incumbent 0.6140 (+0.0383).**

## The null, re-derived rather than borrowed

A matched random feature through the identical harness — same ligands, same per-ligand
argmax, same 64-draw random tie-breaking — 4,000 draws:

| | |
|---|---|
| mean gain | +0.0002 |
| sd | 0.0085 |
| **p95** | **+0.0141** |
| **p99** | **+0.0197** |

This reproduces FINDING 007's +0.0138 / +0.0196 almost exactly, which is a check on the
harness rather than a new number. It is used below only as the gain-versus-random bar.

## Part 1 — the diagnostic, before any selection

Within-ligand spread of each term, and how the incumbent's within-ligand *ranking* relates
to the two halves it decomposes into:

| | within-ligand sd | within-ligand mean |
|---|---|---|
| `mix` | 0.340 Å | 1.214 Å |
| **`cen`** | **0.469 Å** | 1.085 Å |
| **`ori`** | **0.162 Å** | 1.131 Å |
| `pax` | 7.12° | 25.1° |

| | within-ligand ρ with `mix` |
|---|---|
| `cen` | +0.667 |
| **`ori`** | **+0.879** |
| `pax` | +0.566 |

**Half the premise is confirmed and half is not.** In absolute Å the centroid term does
dominate the incumbent's spread — 0.469 Å against 0.162 Å, a factor of 2.9 — which is what
FINDING 024 predicted. But in **rank**, which is all selection uses, the incumbent already
tracks the orientation term far more closely (ρ = +0.879) than the translation term
(ρ = +0.667). The Chamfer was never ranking poses by translation, whatever its units did.

## Part 2 — per-term selection

Ties broken at random over 64 draws throughout. `rho` is the within-ligand Spearman
between the **distance** and LDDT-PLI, so negative is correct; `f(ρ<0)` is the fraction of
ligands on which the term points the right way; `f > rnd` is the fraction of ligands where
the selected pose beats that ligand's own mean.

| selector | selected | gain | ρ | f(ρ<0) | f > rnd | p vs null |
|---|---|---|---|---|---|---|
| `mix` **(incumbent)** | 0.6140 | **+0.0383** | −0.231 | 0.713 | 0.621 | 0.0000 |
| `cen` (translation) | 0.6067 | **+0.0309** | −0.195 | 0.713 | 0.575 | 0.0000 |
| **`ori` (orientation)** | **0.6181** | **+0.0424** | **−0.241** | **0.759** | 0.655 | 0.0000 |
| `pax` (3-axis angle) | 0.6107 | +0.0350 | −0.209 | 0.724 | 0.655 | 0.0000 |
| `pax1` (long axis) | 0.6010 | +0.0253 | −0.203 | 0.690 | 0.609 | 0.0015 |
| `cen+ori` rank-avg | 0.6086 | +0.0329 | — | — | 0.609 | 0.0000 |
| `cen+ori` z-sum | 0.6117 | +0.0360 | — | — | 0.598 | 0.0000 |
| **`mix+ori` rank-avg** | **0.6233** | **+0.0476** | — | — | 0.655 | 0.0000 |
| `ori+pax` rank-avg | 0.6057 | +0.0300 | — | — | 0.632 | 0.0000 |
| `cen+ori+pax` rank-avg | 0.6005 | +0.0247 | — | — | 0.609 | 0.0020 |

Two readings, and the second is the one that matters.

**The orientation term is the strongest single parameter-free term this project has
measured** — +0.0424 against random, ρ = −0.241, correct on 75.9% of ligands, three times
the p99 null. It beats the incumbent on the headline number.

**And the translation term is not the passenger the hypothesis said it was.** `cen` alone
scores +0.0309, i.e. **81% of the incumbent's whole gain**, from nothing but the distance
between two centroids. A coordinate that carries almost no signal does not do that.

The reason is a distinction the hypothesis elided. FINDING 024 measured the pose centroid
against the **crystal**, and found it right to 1.07 Å with isotropic scatter. This measures
the pose centroid against the **consensus of six independent engines**, where the
within-ligand sd is 0.47 Å and disagreement is informative: a pose whose centre six other
engines do not reproduce is in a different sub-basin, whether or not the *population*
centroid error is unbiased. **Agreement with a consensus is not the same quantity as error
against truth, and being unbiased against truth does not make a coordinate uninformative
for consensus.** That is the transferable lesson here.

## Part 3 — the paired test, which is what decides it

Per-ligand paired differences over the 87 ligands, 20,000 bootstrap draws, plus Wilcoxon,
plus Holm correction over the nine tests actually run:

| selector | diff vs `mix` | bootstrap 95% CI | Wilcoxon | Holm | better/worse/tied |
|---|---|---|---|---|---|
| **`mix+ori` rank-avg** | **+0.0094** | **[−0.0004, +0.0214]** | 0.069 | **0.62** | 19 / 10 / 58 |
| **`ori` alone** | **+0.0041** | **[−0.0106, +0.0190]** | 0.606 | **1.00** | 21 / 17 / 49 |
| `pax` | −0.0033 | [−0.0247, +0.0174] | 0.792 | 1.00 | 37 / 34 / 16 |
| `cen+ori` z-sum | −0.0023 | [−0.0203, +0.0141] | 0.763 | 1.00 | 22 / 16 / 49 |
| `cen+ori` rank-avg | −0.0054 | [−0.0248, +0.0124] | 0.757 | 1.00 | 22 / 25 / 40 |
| `cen` | −0.0073 | [−0.0273, +0.0117] | 0.733 | 1.00 | 29 / 32 / 26 |
| `ori+pax` rank-avg | −0.0083 | [−0.0284, +0.0102] | 0.549 | 1.00 | 29 / 33 / 25 |
| `pax1` | −0.0130 | [−0.0332, +0.0056] | 0.212 | 1.00 | 34 / 42 / 11 |
| `cen+ori+pax` rank-avg | −0.0135 | [−0.0352, +0.0057] | 0.715 | 1.00 | 26 / 23 / 38 |

**Nothing clears the bar.** Every confidence interval crosses zero. The best raw p is
0.069 and dies at 0.62 under Holm over the nine comparisons that were made.

The decisive statistic is the tie column: **`mix` and `ori` select the same pose on 49 of
87 ligands (56.3%)**. The orientation term is, on more than half the pool, the incumbent
under another name — which is exactly what ρ(`mix`,`ori`) = +0.879 predicted. On the 38
ligands where they differ, `ori` wins 21 and loses 17. That is a coin flip, and this is
the third time in this project that a +0.004 to +0.009 headline difference has evaporated
under a paired test (FINDING 023 addendum 3 is the previous one).

### The consistency that is real, and is not significance

Each of the four seeds is an independent 5-pose pool of the same 87 ligands:

| seed | random | oracle | `mix` | `ori` | diff | Wilcoxon |
|---|---|---|---|---|---|---|
| run5_s42 | 0.5747 | 0.6620 | 0.6096 | 0.6099 | +0.0003 | 0.69 |
| run5_s43 | 0.5761 | 0.6476 | 0.5937 | 0.5997 | +0.0060 | 0.17 |
| run5_s44 | 0.5767 | 0.6399 | 0.5948 | 0.6051 | +0.0103 | 0.11 |
| run5_s45 | 0.5754 | 0.6517 | 0.6019 | 0.6037 | +0.0018 | 0.36 |

`ori` is ahead in **4 of 4**, and `mix+ori` likewise in 4 of 4 (+0.0001, +0.0044, +0.0057,
+0.0023). Four sub-pools of one pool are not four experiments — they share the ligands, the
truth and the references, so this is a consistency check, not a replication, and no
individual seed reaches p < 0.10. It is enough to say the sign is stable and not enough to
say the effect exists.

Splitting by pool spread, where FINDING 012 says a cross-engine gain should concentrate:

| | n | `mix` | `ori` | diff | Wilcoxon |
|---|---|---|---|---|---|
| low-spread ligands | 44 | 0.6156 | 0.6144 | −0.0011 | 0.70 |
| high-spread ligands | 43 | 0.6124 | 0.6218 | **+0.0095** | 0.61 |

The direction is right — whatever `ori` adds, it adds where the pool disagrees with
itself — but at n = 43 and p = 0.61 this is a hypothesis for the next pool, not a result
on this one.

## Is this FINDING 006's azimuth null again?

**No, and the distinction matters.** FINDING 006 measured azimuth-about-the-Fe-axis
consensus at **−0.0078, p = 0.57** — a term with no signal, correctly retired. The
orientation term measured here has a great deal of signal: **+0.0424, three times the p99
null, ρ = −0.241, correct on 75.9% of ligands.** It is the strongest single unfitted term
in the project.

What fails is not the term, it is the *decomposition*. The incumbent already contains
almost all of it. So the honest statement is narrower than "orientation consensus is
closed": **orientation consensus works, and separating it out of the Chamfer buys nothing
measurable, because the Chamfer was already ranking by it.** FINDING 024's premise — that
the incumbent wastes its dynamic range on translation — is refuted by ρ(`mix`,`ori`)
= +0.879 and by `cen` alone scoring +0.0309.

## What this changes

1. **`cypstruct.xengine.select()` is unchanged.** Shipping `ori` would trade a
   long-replicated selector (+0.0395 / +0.0357 / +0.0383 on three pools, never retuned)
   for a paired difference of +0.0041 at p = 0.61.
2. **R1 of FINDING 024 is answered and closed.** The one proposal whose mechanism followed
   directly from the error decomposition returns a null against the incumbent. That is the
   **fifth** consecutive line — fine-tuning, the heme bond, ATOMICA, embeddings, now the
   geometric decomposition — where the honest answer is that the incumbent still wins.
3. **`pax` is worth keeping as a diagnostic, not a selector.** An inertia-tensor angle with
   no distance computation at all reaches +0.0350, 91% of the incumbent, at a fraction of
   the cost. That is useful for screening large pools cheaply; it is not an improvement.
4. **One pre-registered follow-up, and only one.** `mix + ori` rank-average is positive in
   4 of 4 seeds with a CI whose lower edge is −0.0004. If a *genuinely* new pool is ever
   scored — new ligands or new references, not new seeds of these poses — test that one
   combination, once, and accept only a paired Wilcoxon p < 0.05. Testing it again on this
   pool would be testing it on the data that generated the hypothesis.
5. **For the world model:** consensus agreement and error-against-truth are different
   quantities, and FINDING 024's isotropic-centroid result does not license the inference
   that centroid agreement is uninformative. It scores +0.0309 on its own.

## Method notes worth keeping

- **Ties were broken at random over 64 draws.** With continuous distances no ties occur,
  but the harness was built with it anyway because a constant score otherwise selects
  `_model_0` and reports Boltz confidence under another name (FINDING 023).
- **The null was re-derived, not borrowed**, and landed on +0.0141 / +0.0197 against the
  inherited +0.0138 / +0.0196. That agreement is the evidence the harness is sound.
- **The gain-versus-random null is the wrong bar here and was not used as one.** Both `ori`
  (+0.0424) and `mix+ori` (+0.0476) sail past +0.0197 and neither beats the incumbent. Any
  write-up quoting only the first number would have been a confident, wrong result — the
  precise failure FINDING 023 addendum 3 was written to prevent.
- **Nine comparisons were made and nine are reported**, with Holm applied to all nine. The
  best raw p of 0.069 is what one expects from nine draws on a null.
