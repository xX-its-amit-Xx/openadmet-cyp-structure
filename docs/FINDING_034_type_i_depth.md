# FINDING 034 — depth on the predicted-Type-I stratum: the law holds, the purchase does not

**Date:** 2026-09-22 · **Status:** measured, **124 OpenProtein jobs**, 424 complexes,
0 failures, ~36 minutes wall clock, **$0** · **Verdict:** **REFUTED** for the intervention
as it could actually be bought, **MEASURED** for the depth law underneath it — which is
intact, and is 1.8× steeper on Type I than on Type II.

Pre-registered in `docs/PREREG_type_i_depth.md`, committed at `686487a` before a single
job was submitted. No threshold, rung or ligand below was moved afterwards.

---

## What was asked

`FINDING_033` closed with a recommendation rather than a result:

> **3. Spend the remaining budget on POOL DEPTH for the Type I ligands, not on scoring.**
> Type I oracle is 0.058 below Type II, and FINDING 004 measured the oracle still climbing
> with sampling at +0.0125 of selected score per doubling. Doubling the pool on the
> predicted-Type-I subset is a small number of jobs and attacks the part of the gap that
> selection cannot.

That deserves testing rather than assuming, because this repo holds a hard-won
distinction. **Pool EXPANSION** — adding a *different kind* of pose — has raised the
oracle and left selection flat or worse four times running: a second engine (013, +0.0375
of unreachable oracle), a sampler sweep (016, −0.0038), 1.78 M rigid rotations (027,
oracle +0.0108 / selection **−0.0145**), a second copy of the query ligand (031, oracle
+0.0237 / selection **−0.0027**). **Pool DEPTH** — more samples of the same kind — is the
one lever that has paid, at +0.0125 per doubling (004).

> **The question.** Does depth on the predicted-Type-I ligands convert into SELECTED
> LDDT-PLI, and at what rate per doubling, against the +0.0125 benchmark?

**This is depth, not expansion — and that claim is itself under test.** If the oracle
climbs and selection does not, this is expansion behaviour wearing depth's clothes.

---

## What was done

| | |
|---|---|
| stratum | **prediction-side**: median `fe_donor_dist` over a ligand's own 20 pool poses **> 2.6 Å**. The rule `FINDING_033` validated at 96.6%. **n = 14**, and nothing from a crystal enters it |
| part 1 | **subsample** the existing 20-pose Boltz-2 pool at rungs 1/2/3/5/8/10/14/20, 256 draws, per stratum. Free, exact, internally consistent |
| part 2 | **generate** new poses for those 14 ligands and re-measure at 20 + k |
| venue | **OpenProtein**, `boltz2`, single-sequence, `diffusion_samples=1`, `num_recycles=3`. `cypstruct.budget.preflight` passed on every submission (`openprotein_jobs`, cap 2,000) |
| cost | **124 jobs** — 6 yield probe + 58 (replicates 0–15) + 56 (replicates 16–29) + 4 cache probe — **424 complexes, 0 failures, $0, 144 MB on C:, never on D:** |
| Modal | **not touched.** It is over its spend cap |

### Why the new arm is not quite the same kind of pose, and why it could not be

The existing pool is Boltz-2 on Modal with a 6,979-sequence MSA and an explicit
Cys442-SG→heme-FE bond. **`boltz2` on OpenProtein rejects an uploaded MSA** — verified on
the stored probe job before the pre-registration was written: `boltz_2__msa` →
`JobStatus.FAILURE, "internal server error"`, while single-sequence `boltz2` →
`JobStatus.SUCCESS`. Modal is over cap and Explorer cannot be MSA-staged inside the
interim deadline. So the single available depth purchase carries a confound, it was named
**C1** in the pre-registration before any pose was scored, and it is the mechanism behind
the result below. `FINDING_019` measured the missing heme bond (**C2**) as a null,
Δ −0.0069 at p = 0.43.

---

## Controls, with counts, including the zeros

| id | control | result |
|---|---|---|
| **N2** | the shipped column reproduces | selected **0.6164**, oracle **0.6975**, random **0.5769**, gain **+0.0395**, within-ligand ρ **−0.2582**, correct sign **75.86%**. Every one matches `FINDING_011`/`033` to four decimals |
| **N2b** | the `argmin(xeng)` shortcut equals `cypstruct.xengine.select()` | 0.616391 vs 0.616391, **equal to 1e-9** on all 87. z-scoring is monotone within a ligand, so this is an identity, and it is asserted rather than assumed |
| **N1** | numbering (FINDING 021) on the 420 new poses | renumbering offset **0 on 420 of 420**, minimum residue-name identity **0.9979**, **0** poses below the 0.95 bar |
| **C7** | the 4 exact-zero rows in the existing pool | PG4 ×3, PG0 ×1, BiSyRMSD **23.1 / 25.9 / 26.9 / 26.1 Å**, `mapped` True on all four, the same ligands' other poses at 0.441 / 0.577. Genuine ejections, as `FINDING_033` found |
| **C7′** | 24 exact-zero rows in the **new** arm | PG4 and ERY only, BiSyRMSD **22.5–28.0 Å**, ligand **21.7–25.2 Å from the iron**, offset 0, `mapped` True. Ejections again, not numbering |
| **N4** | distinct poses, not jobs | 420 collected rows → **151 distinct**, 10.79 per ligand. Every depth number below is quoted on the deduplicated pool |
| **N5** | filters | 87 validation ligands → 14 pass `median fe_donor_dist > 2.6` → 73 excluded as predicted Type II. Skips during scoring: no_crystal 0, no_ligand 0, load_failed 0, **mapping_failed 0**, no_heme_frame 0, low_renumber_identity 0 |
| | reference depth per stratum ligand | 6–11 independent poses (min 6), all ≥ the 4 that `FINDING_011` requires |
| | pool diversity of the **existing** pool | **20 unique `xeng` values on every one of 87 ligands**, median within-ligand LDDT-PLI sd 0.062, zero ligands with sd 0. The 20 samples are 20 poses |

---

## 1. The stratum, on the label we would actually have

The pre-registered rule is prediction-side, because that is what exists at submission
time. Against the crystal label it disagrees on **three of 87**, two of which touch this
list:

| ligand | pdb | crystal Fe | predicted median | in the stratum? |
|---|---|---|---|---|
| **D0R** | 3TJS | **2.107 Å (Type II)** | 5.430 Å | **yes** — the model fails to coordinate a genuine coordinator |
| **QDY** | 6UNJ | 2.824 Å (Type I) | **2.231 Å** | **no** — the model forces coordination on a non-coordinator, all 20 poses |
| MWV | 6OO9 | 3.763 Å (Type I) | 7.695 Å | yes — right mode, ejected past the active site |

So the prediction-side split is **14 / 73**, not `FINDING_033`'s crystal-side 14 / 73 over
a different membership: D0R swaps in and QDY swaps out. That one swap matters:

| at depth 20 | Type I (pred, n=14) | Type II (pred, n=73) |
|---|---|---|
| random (exact) | **0.4597** | 0.5994 |
| **oracle** | **0.6499** | 0.7066 |
| selected | **0.5270** | 0.6335 |
| gain | **+0.0673** | +0.0342 |
| its own null p95 / p99, 4,000 draws at this n | **+0.0431 / +0.0631** | +0.0140 / +0.0192 |
| **empirical p** | **0.0063** | 0.0000 |
| drop best 1 / best 2 | +0.0462 / +0.0234 | +0.0316 / +0.0291 |
| ligands with positive gain | 57.1% (8/14) | 69.9% |
| within-ligand ρ, correct sign | −0.2357, 71.4% | −0.2625, 76.7% |

**On the label we would actually have, the selector clears its own n=14 noise floor
(p = 0.0063)** where `FINDING_033`'s crystal-side stratum could not (+0.0411, p = 0.0595).
That is not a new selector and not a new pool — it is one ligand moving. It should be read
with the leverage row next to it: drop the two best ligands and +0.0673 becomes +0.0234,
below the floor. The honest statement is that the prediction-side stratum is *consistent
with* the selector working on Type I and still cannot prove it at n = 14. The family-wide
test in `FINDING_033` (169 Type I pairs, matched difference −0.0030 ± 0.010) remains the
evidence that it does.

The absolute gap `FINDING_033` identified is confirmed on this labelling and is slightly
larger: **random −0.140, oracle −0.057, selected −0.107**.

---

## 2. The depth curve — oracle first, then selection

Subsampling the existing pool. 256 draws per rung, random tie-breaking, null recomputed
**inside each stratum at each rung** because the FINDING 007 floor is n-dependent.

### Type I (predicted), n = 14

| depth | **oracle** | selected | random (exact) | gain | its own null p95 |
|---|---|---|---|---|---|
| 1 | 0.4600 | 0.4600 | 0.4597 | +0.0003 | +0.0473 |
| 2 | 0.5058 | 0.4775 | 0.4597 | +0.0178 | +0.0438 |
| 3 | 0.5280 | 0.4806 | 0.4597 | +0.0209 | +0.0461 |
| 5 | 0.5587 | 0.4852 | 0.4597 | +0.0256 | +0.0447 |
| 8 | 0.5847 | 0.4897 | 0.4597 | +0.0301 | +0.0445 |
| 10 | 0.5998 | 0.4937 | 0.4597 | +0.0340 | +0.0432 |
| 14 | 0.6223 | 0.5054 | 0.4597 | +0.0457 | +0.0459 |
| **20** | **0.6499** | **0.5270** | 0.4597 | **+0.0673** | +0.0432 |

### Type II (predicted), n = 73 · and all 87

| depth | oracle II | selected II | oracle all | selected all |
|---|---|---|---|---|
| 1 | 0.5988 | 0.5988 | 0.5764 | 0.5764 |
| 2 | 0.6324 | 0.6123 | 0.6120 | 0.5908 |
| 3 | 0.6473 | 0.6170 | 0.6289 | 0.5952 |
| 5 | 0.6633 | 0.6200 | 0.6464 | 0.5977 |
| 8 | 0.6783 | 0.6227 | 0.6633 | 0.6003 |
| 10 | 0.6848 | 0.6238 | 0.6712 | 0.6028 |
| 14 | 0.6957 | 0.6269 | 0.6843 | 0.6076 |
| **20** | **0.7066** | **0.6335** | **0.6975** | **0.6164** |

**The random baseline is flat to four decimals at every rung in every stratum**
(per-doubling slope −0.0001 / +0.0001 / +0.0001), which is the check that the subsampling
is unbiased. It passes.

### The rate per doubling, against the +0.0125 benchmark

| stratum | **oracle / doubling** | **selected / doubling** | **conversion** |
|---|---|---|---|
| **Type I (n=14)** | **+0.0428** [+0.0284, +0.0574] | **+0.0127** [−0.0004, +0.0280] | **29.6%** |
| Type II (n=73) | +0.0242 [+0.0203, +0.0283] | +0.0069 [+0.0037, +0.0104] | 28.5% |
| all (n=87) | +0.0272 [+0.0230, +0.0317] | +0.0078 [+0.0043, +0.0115] | 28.6% |
| *local rate in the top octave, 10 → 20* | *I +0.0501 · II +0.0218* | ***I +0.0333 · II +0.0097*** | |

**Three things, in order of importance.**

1. **Depth converts at the same fraction in both strata — 29.6% against 28.5%.** Type I
   gets a bigger absolute gain per doubling **only** because its oracle climbs 1.8× faster
   (+0.0428 vs +0.0242): its pool is further from saturation, exactly as its larger
   headroom (0.190 vs 0.107) predicts. Nothing about binding mode changes how depth is
   converted; it changes how much there is to convert. That is the same conclusion
   `FINDING_033` reached about the selector, arrived at from the generation side.
2. **The Type I selection rate is +0.0127 per doubling — the `FINDING_004` benchmark to
   three decimals** — and 1.8× the Type II rate. In the top octave, which is the region a
   20→40 purchase would actually extend, it is **+0.0333**, 3.4× Type II's +0.0097.
   So on this evidence `FINDING_033`'s recommendation is *quantitatively right about the
   rate*: depth is worth more on Type I than on Type II, by about a factor of two to three.
3. **The +0.0125 benchmark is selector-dependent and does not reproduce for the shipped
   selector.** `FINDING_004` measured it with the `FINDING_003` selector. Re-measured here
   on the same pool with `-z(xeng)`, the all-87 rate is **+0.0078**, not +0.0125. The
   like-for-like comparison is therefore Type I **+0.0127 against +0.0078**, a ratio of
   1.6. The headline number is unchanged in direction and should be quoted with its
   selector attached from now on.

At n = 14 the Type I selection rate's CI touches zero ([−0.0004, +0.0280]). The oracle
rate does not ([+0.0284, +0.0574]).

---

## 3. What the new poses actually are

420 complexes, 151 distinct poses, **10.79 per ligand** (8–12). And a venue fact that cost
56 of the 124 jobs to learn:

| replicate jobs per ligand | 2 | 4 | 6 | 8 | 12 | **16** | 20 | 24 | **30** |
|---|---|---|---|---|---|---|---|---|---|
| mean distinct poses | 1.86 | 2.57 | 3.50 | 5.00 | 7.21 | **10.79** | 10.79 | 10.79 | **10.79** |

**Replicates 16 through 29 — 196 complexes, 56 jobs — added exactly zero distinct
poses.** Checked by md5 on all fourteen ligands and not inferred from the count:
**14 of 14 ligands have all fourteen of replicates 16–29 byte-identical to replicate 15**,
196 of 196 complexes. This is `FINDING_009`/`015`/`031` for the fourth time, in a fourth
shape.

The obvious explanation is a server-side result cache, and **it was tested rather than
assumed**: four further jobs with a byte-identical payload, submitted as their own wave.
They returned a pose identical to each other to 0.0000 Å and **different from replicate 15
by 3.14 Å and from replicate 0 by 8.16 Å** — a *new* twelfth mode. So it is neither a
permanent payload cache nor a finite model repertoire. Within the first wave the
duplication is largely pairwise among adjacent replicates (MWS: 6=7, 8=9, 10=11, 12=13;
partially so on YNV and 08J), i.e. jobs submitted within a second or two of each other
tend to agree.

**Whatever the mechanism, the operational rule is the reportable part: distinct poses
track submission WAVES, not jobs, and replicate count is not a purchase order for pose
count.** Doubling the replicates from 15 to 30 bought nothing at all.

### The new arm's quality, stated before it is mixed into anything

| | new arm (OpenProtein boltz2, single-sequence) | existing pool (Modal Boltz-2 + MSA + heme bond) |
|---|---|---|
| random (pool mean) | **0.4025** | 0.4597 |
| oracle at ~10.8 distinct poses | **0.5994** | 0.5998 *(at depth 10)* |

**The new configuration's best pose is as good; its average pose is 0.057 worse.** Losing
the MSA widened the pool downward without lowering its ceiling. That is confound C1,
measured.

And the decisive control — **does depth convert differently inside the new configuration?**
Both arms subsampled on the same rungs 1/2/3/5/8, same 14 ligands:

| | oracle / doubling | selected / doubling | **conversion** |
|---|---|---|---|
| **new arm alone** | +0.0587 | +0.0134 | **22.8%** |
| **existing pool, same rungs** | +0.0414 | +0.0094 | **22.7%** |

**22.8% against 22.7%.** The depth law is a property of the problem, not of the
configuration. A homogeneous pool of the cheaper poses converts depth exactly as well.

---

## 4. The augmented pool — the primary endpoint

All 20 existing poses kept, `k` distinct new poses added, 256 draws, matched depth
`N = 8` (the minimum distinct count over the 14 ligands).

| pool | **oracle** | selected | random |
|---|---|---|---|
| **20** (existing only) | **0.6499** | **0.5270** | 0.4597 |
| 22 (+2 new) | 0.6527 | 0.5251 | 0.4545 |
| 24 (+4 new) | 0.6560 | 0.5234 | 0.4506 |
| **28 (+8 new)** | **0.6617** | **0.5205** | 0.4434 |
| *all available, mean depth 30.8* | *0.6670* | *0.5188* | *0.4403* |

| | value | 95% CI | sign |
|---|---|---|---|
| **Δ selected (20 → 28)** | **−0.0065** | **[−0.0276, +0.0135]** | 1 better / 3 worse / **10 unchanged** |
| Δ oracle (20 → 28) | **+0.0118** | [+0.0001, +0.0282] | |
| Δ random (20 → 28) | −0.0163 | | |
| Δ selected, all available (20 → 30.8) | **−0.0082** | | |
| Δ oracle, all available | **+0.0171** | | |
| poses the selector actually took from the new arm | **4 of 14** | | Wilcoxon p = 0.465 |

**The oracle rises and the score falls.** On ten of fourteen ligands the selector keeps an
old pose and nothing changes; on three of the four where it switches, it switches to a
worse one. This is CLAIM F of `CLAUDE.md` and FINDINGS 013 / 016 / 027 / 031 — **for the
fifth time.**

### The pre-registered verdict

| clause | required | measured | fires? |
|---|---|---|---|
| SHIPS (1) | Δ_sel ≥ +0.020 | **−0.0065** | no |
| SHIPS (2) | 95% CI excludes 0 | [−0.0276, +0.0135] | no |
| SHIPS (3) | Type I selection rate ≥ +0.0125 | **+0.0127** | **yes** |
| MEASURED — expansion | Δ_oracle ≥ +0.020 with Δ_sel short | +0.0118 (matched) / +0.0171 (all) | **no**, narrowly |
| **REFUTED** | **Δ_sel ≤ 0 and CI upper < +0.020** | **−0.0065, upper +0.0135** | **YES** |
| MEASURED — unbuyable | yield < 0.25 or < 5 distinct new poses | yield 0.50 at the probe, 8–12 distinct | no |

**REFUTED**, by the clause that excludes a useful effect rather than the one that merely
fails to find one — and **SHIPS (3) fires**, which is the whole shape of this finding: the
*rate* `FINDING_033` relied on is real and is the largest in the repo, and the *purchase*
available at the permitted venue destroys it.

---

## 5. Projection onto a Type I-rich test set

Expected selected LDDT-PLI at test-set Type I fraction `f`, per-stratum bootstrap over
ligands, 10,000 resamples:

| f | expected now | 95% CI | after the depth spend as executed | after a **matched-quality** doubling |
|---|---|---|---|---|
| 0.00 | 0.6335 | [0.5917, 0.6734] | 0.6335 | 0.6335 |
| 0.25 | 0.6069 | [0.5707, 0.6423] | 0.6053 (**−0.0016**) | 0.6101 to 0.6152 |
| **0.50** | **0.5803** | **[0.5385, 0.6211]** | **0.5770 (−0.0032)** | **0.5866 to 0.5970 (+0.006 to +0.017)** |
| 1.00 | 0.5270 | [0.4550, 0.6009] | 0.5205 (−0.0065) | 0.5397 to 0.5603 |

The matched-quality column applies the §2 curve — the whole-curve rate **+0.0127** and the
top-octave rate **+0.0333** — to one doubling of the Type I half only. **At a 50/50 test
set, a genuine 20 → 40 doubling on the Type I ligands is worth +0.006 to +0.017
LDDT-PLI**, against the selector's own whole contribution of +0.0395. Worth having,
smaller than the uncertainty on the test-set composition, and **not available at this
venue**.

---

## Verdict — **REFUTED** (the purchase) · **MEASURED** (the law)

1. **The depth law is real, and it is steepest exactly where `FINDING_033` said to spend.**
   Type I converts depth at **+0.0127 selected per doubling** (+0.0333 in the top octave)
   against Type II's +0.0069 (+0.0097) and the shipped selector's all-87 +0.0078. Its
   oracle climbs at **+0.0428 per doubling** with no saturation at 20.
2. **Depth is converted at the same 29% in both strata, and at the same 23% inside a
   completely different configuration.** Binding mode does not change how depth converts,
   only how much there is to convert. Three independent measurements now say the same
   thing about mode-dependence from three directions (033's selector test, 033's
   headroom-matched family test, this conversion ratio).
3. **The depth actually purchasable is not depth.** Losing the MSA — forced, because
   `boltz2` with an uploaded MSA fails server-side on OpenProtein — produced poses whose
   ceiling matches and whose average is 0.057 lower. Mixed into the pool they add
   **+0.0118 to +0.0171 of oracle and subtract 0.0065 to 0.0082 of score**. Fifth
   consecutive confirmation that the pool's ceiling and the submitted score are different
   quantities.
4. **Replicate jobs are not a purchase order for poses.** 30 replicates bought 10.8
   distinct poses; replicates 16–29 bought **zero**, at a cost of 56 jobs. Diversity tracks
   submission waves, not jobs, and a fresh wave with an identical payload produced a new
   pose that all four of its jobs shared exactly.

**What this does NOT license.** It does not say depth is worthless — the opposite, and the
rate is measured. It does not test depth at matched conditioning, which is the experiment
that would actually settle it and which needs Modal or a staged-MSA Explorer run. It does
not extend past 20→40; the curve's slope beyond 40 is unmeasured. And at n = 14 the Type I
selection-rate CI touches zero, so the rate itself is a point estimate with a wide interval.

---

## What to do at submission time — superseding `FINDING_033` item 3

`FINDING_033` items **1, 2, 4 and 5 stand unchanged** (label prediction-side; do not change
the selector; pre-announce a lower absolute score if the set is Type I-rich; watch QDY's
failure mode). Item 3 is replaced by:

3. **Do not buy Type I depth from a degraded configuration. Buy it at matched
   conditioning or not at all.** Additional poses only pay if they come from the same
   generator as the pool they join: measured here, a matched-quality doubling of the Type I
   half is worth **+0.006 to +0.017** at a 50/50 test set, and the MSA-less substitute is
   worth **−0.003**. The requirement is an MSA and, ideally, the heme bond — i.e. Boltz-2
   on Modal (over cap) or Explorer with a login-node-staged MSA (`CLAUDE.md`'s budgeted
   task), **not** OpenProtein, where the only `boltz2` configuration that runs is
   single-sequence.
3b. **If depth is bought, buy every replicate in ONE submission wave.** Distinct poses
   saturated at 16 replicates and the next 14 returned duplicates. Count distinct poses
   before paying for more, per configuration, every time — and re-count, because this is
   the fourth configuration in which the answer changed.
3c. **The prediction-side stratum is the one to report.** It is what exists at submission
   time, it differs from the crystal-side stratum by one ligand, and on it the shipped
   selector clears its own n=14 floor (+0.0673, null p95 +0.0431, p = 0.0063) — which the
   crystal-side stratum could not. Two ligands carry most of it; quote the leverage
   (drop-best-2 → +0.0234) alongside.

---

## Artefacts

| what | where |
|---|---|
| pre-registration | `docs/PREREG_type_i_depth.md` (`686487a`) |
| runner (plan / submit / collect / yield / score) | `scripts/cofold/type_i_depth.py` |
| analysis (controls / subsample / stratum / augment / project) | `scripts/structure/type_i_depth_analysis.py` |
| controls, incl. N2 and C7 | `data/processed/type_i_depth_controls.json` |
| the depth curve, all three strata | `data/processed/type_i_depth_subsample.json` |
| per-ligand, per-rung | `data/processed/type_i_depth_subsample_per_ligand.csv` |
| stratum at full depth, own null, leverage | `data/processed/type_i_depth_stratum.json` |
| the 420 new poses, scored | `data/processed/type_i_depth_poses.csv` |
| distinct-pose saturation curve | `data/processed/type_i_depth_yield_curve.json` |
| the identical-payload probe | `data/processed/type_i_depth_cache_probe_result.json` |
| augmented pool, primary endpoint, new-arm curve | `data/processed/type_i_depth_augment.json` |
| projection | `data/processed/type_i_depth_projection.json` |
| job ledger | `data/processed/type_i_depth_jobs.json` |
| the 424 mmCIFs (144 MB) | `C:\cyp_struct\type_i_depth\` — scratch, never `D:` |
