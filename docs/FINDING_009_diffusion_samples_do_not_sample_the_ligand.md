# FINDING 009 — `diffusion_samples` does not sample the ligand, and a second engine's confidence fails exactly like the first

**Date:** 2026-09-13
**Engine:** Protenix-v2 via OpenProtein.ai
**Data:** 20 collected CYP3A4 ligands × 20 requested samples

---

## 1. The bug that would have quietly wasted the whole campaign

I asked for `diffusion_samples=20` and got back 20 models per complex. They are not 20
poses. Across all 20 models of a single job:

| quantity | per-atom sd across the 20 models |
|---|---|
| **ligand** coordinates | **0.0000 Å** |
| **heme** coordinates | **0.0000 Å** |
| protein coordinates | 13.0 Å |

The ligand and the cofactor are byte-identical in every model; only the protein moves.
This is in the **raw mmCIF returned by the API**, not an artefact of splitting it — the
per-model ligand centroid is `[-4.192, 2.212, -2.101]` in model 0, 1 and 2 alike.

Two *separate jobs* for the same ligand, by contrast, differ by **8.63 Å ligand RMSD**.

So on this engine **pose diversity comes from replicate jobs, not from
`diffusion_samples`**. The campaign as first launched — 87 ligands × 20 samples — bought
87 distinct ligand poses, not 1,740, while looking in every log and manifest exactly like
a pool twenty times its real size.

The runner now takes `--replicates`, resumes on `(replicate, ligand)` so a re-run tops the
pool up rather than resubmitting it, and names files `{lig}__r{rep}s{k}.cif`.

### Why this was nearly invisible

Every surface check passed. 20 files, all with **distinct md5s**. LDDT-PLI varying widely
within a ligand, 0.04 to 0.79. A plausible oracle of 0.6664. The variation is real — it
is just the *protein* moving around a fixed ligand, which changes protein–ligand contacts
and therefore changes LDDT-PLI.

What actually exposed it was a number that was **too clean**: the Fe-to-nearest-ligand-atom
distance was identical to four decimal places across all 20 samples (sd 0.0000). Twenty
independent diffusion samples do not agree to 0.00005 Å. That is the same tell as the
pyrrole N–H filter that scored every nitrogen 1.00, and the third time this session that
an implausibly tidy number has been the thread worth pulling.

**Rule to carry forward:** when a new engine joins, verify that the samples differ *in the
thing being sampled* before scoring any of them. Distinct file hashes are not evidence of
distinct poses.

---

## 2. Protenix's confidence does not rank poses either

FINDING 001-E established that Boltz-2's `complex_ipde` ranks poses within a ligand at
chance (ρ = −0.092 / −0.033) and that picking the most confident pose is worse than
picking at random. The obvious hope was that this was a Boltz problem.

It is not. Protenix-v2 exposes seven confidence fields. Within-ligand Spearman against
true LDDT-PLI, over 16 ligands × 20 samples:

| field | mean ρ | median ρ | frac ρ > 0 | p (Wilcoxon) |
|---|---|---|---|---|
| ranking_score | +0.020 | −0.011 | 0.50 | 0.67 |
| ptm | −0.002 | −0.000 | 0.50 | 0.98 |
| iptm | +0.011 | −0.014 | 0.50 | 0.71 |
| plddt | −0.034 | +0.008 | 0.56 | 0.78 |
| gpde | −0.018 | −0.006 | 0.50 | 0.86 |
| has_clash, disorder | — | — | no within-ligand variance | — |

`frac ρ > 0` is **0.50** for four of five fields — indistinguishable from a coin.

Note the caveat that matters: because of §1 these 20 samples share one ligand pose, so
this measures whether confidence ranks *protein conformations* around a fixed ligand. It
is still the ranking the engine offers for choosing which structure to submit, and it is
still at chance. The replicate pool will let the same test run over genuinely different
ligand poses.

**Two independent architectures, trained separately, both fail to rank their own samples.**
This is no longer a quirk of one model — it is the central fact of the problem, and it is
what makes the physics scorer the whole project rather than a nice-to-have.

---

## 3. Protenix is a third correlated engine, not a decorrelated one

Per the FINDING 005 gate, on the 16 ligands collected so far:

| | Protenix-v2 | Chai-1 (FINDING 005) |
|---|---|---|
| oracle ρ vs Boltz | **+0.474** (p = 0.064) | +0.45 |
| own oracle | 0.6664 | — |
| Boltz oracle, same ligands | 0.7399 | — |
| union oracle | 0.7503 (**+0.0104**) | ~null |
| beats Boltz on | 6 of 16 ligands | — |

Chai +0.45, Protenix +0.474. Two independent engines, the same correlation with Boltz's
errors. The working hypothesis that architectural diversity buys decorrelation is now
0 for 2, and the ligands that are hard appear to be hard for everyone — which points the
remaining upside at scoring rather than at buying yet another engine.

Treat the ρ as provisional: n = 16, p = 0.064, and §1 means this pool is one ligand pose
per ligand. Re-run the gate on the replicate pool before making it a settled number.

---

## What changed as a result

- `--replicates` added; `--samples` demoted to a protein-conformation knob.
- Resume keyed on `(replicate, ligand)`; manifests accumulate across replicates.
- Two collector bugs fixed en route: `gemmi.read_structure_from_string` does not exist
  (it is `read_structure_string`), and a batch was being retired as done even when every
  save inside it failed — which discarded a whole job's poses while reporting success.
