# FINDING 014 — the QM scorer's tier-1 donor prior cannot select, and the reason redirects the plan

**Date:** 2026-09-14
**Check:** 800 unsteered Boltz poses, 40 CYP3A4 ligands
**Cost of the check:** about two minutes, against an xTB pipeline it would otherwise gate

---

## The pre-check the repo's own rules demand

`QM_SCORER_DESIGN.md` §2B makes the per-atom coordination prior the leverage of the whole
design: proton affinity, Fukui f⁻ and %V_bur at each candidate nitrogen, computed **once
per molecule** and then read against each pose's geometry, so 24,000 poses cost only 200 QM
jobs. That cost argument is sound. It is also irrelevant if the term has no within-ligand
variance, which is what FINDING 006's rule exists to catch: *check within-ligand CV before
testing any new feature.*

Which atom actually coordinates, across poses of the same ligand:

| distinct coordinating atoms per ligand | ligands |
|---|---|
| **1 (all poses identical)** | **34 of 40** |
| 3 | 2 |
| 5 | 3 |
| 6 | 1 |

Median 1.0. The donor **element** varies on only 6 of 40.

**So for 85% of ligands a per-atom donor prior is a constant**, and a constant cannot rank
a ligand's poses against each other. This is exactly the shape of `is_coordinated`, which
had within-ligand CV 0.023, was constant on 98% of ligands, and selected nothing.

## Why this was predictable in hindsight, and worth measuring anyway

It is the same saturation that has now defeated five separate ideas. Once the heme is
bonded to Cys442, Boltz reaches the coordination geometry reliably (thesis claim C,
falsified) and *picks the same donor atom every time*. The engine is not uncertain about
which nitrogen binds. It is uncertain about where the rest of the molecule goes — which is
FINDING 010's conclusion arrived at from a different direction.

## What this does and does not kill

**Does not kill tier 1 as chemistry.** The donor prior may well be right about which
nitrogen *should* coordinate; that is a real and checkable claim, and it would matter for
generation, for flagging poses that coordinate through a hindered nitrogen, and for the 6
of 40 ligands where the engines genuinely disagree about the donor.

**Does kill it as the main selector.** A term active on 15% of ligands, and only when it is
also correct, cannot be the leverage the design assigns it.

**Redirects to tier 2.** GFN2-xTB interaction energy and strain on a pocket cutout is
computed **per pose**, so it has within-ligand variance by construction — and it judges a
pose on its own merits rather than by agreement with siblings, which FINDING 013 showed is
the only thing that could collect the +0.0375 of oracle a union pool adds. Tier 2 is
CPU-bound and embarrassingly parallel, and the design already scopes it to the top K≈5
poses per ligand, so the cost is bounded.

**Run the same CV pre-check on tier 2 before building it.** It should pass by construction,
but "should" is what this finding just cost two minutes to disprove about tier 1.

---

## Addendum — tier 2 feasibility, and an honest prior against it

Tier 2 needs GFN2-xTB. Where it can actually run, checked rather than assumed:

| venue | status |
|---|---|
| local | `xtb` absent; `tblite` fails to install (no wheel); OpenMM and RDKit present |
| Explorer (NEU cluster) | connects fine, **no xtb module, no conda on PATH** - needs a from-scratch install |
| Modal | over its $25 cap and blocked; the $140 hackathon workspace token is still not on this box |

So tier 2 costs a real installation on a cluster before a single number comes out.

**And the prior from related terms is poor, which is the part worth stating before
spending the time.** The crude, classical versions of "score a pose by its physical
interaction with the pocket" are already measured, and they do not merely fail - two of
them are significantly *worse than random*:

| term | Δ vs random | p |
|---|---|---|
| `max_clash` (low) | **−0.0359** | 0.0024 |
| `n_contacts` (high) | −0.0122 | 0.50 |
| `s_fe_donor_angle` (high) | **−0.0358** | 0.0007 |
| `n_pocket_residues_touched` (high) | **+0.0108** | — works, and is in the incumbent |

The one member of this family that works is a *count of distinct residues touched*, not an
energy. Steric energy proxies actively mislead. GFN2-xTB would add electrostatics,
polarisation and real strain, which is a genuine difference from a clash count - but it is
a bet that better physics rescues a family where the cheap members are negative.

**Recommendation: do not start the xTB installation on the strength of FINDING 013 alone.**
The argument for tier 2 is real (it is the only known way to collect the +0.0375 of oracle
a union pool adds) but it is an argument from elimination, not from positive evidence that
an interaction energy discriminates here. With 10 days to the interim deadline, a validated
+0.0381 selector already shipping, and ~30 features dead, the expected value is better spent
on the P450 generalisation set - which is accumulating for free - than on a cluster install
with a negative prior.

**What would flip this:** a cheap classical interaction energy (OpenMM is already available)
showing within-ligand variance *and* any positive correlation on held-out ligands. That is
a bounded experiment and the right gate before the xTB work, rather than after it.

---

## The gate ran, and it fails: a classical interaction energy does not discriminate

Implemented the **AutoDock Vina** scoring terms - a real published function, not a proxy
invented here - over 900 poses / 45 ligands. Two gates, set in advance.

**GATE 1 passes.** Unlike the tier-1 donor prior, these terms genuinely vary within a
ligand, so they *can* rank its poses:

| term | median within-ligand CV |
|---|---|
| hbond | 0.473 |
| repulsion | 0.328 |
| gauss1 | 0.096 |
| full Vina score | 0.081 |
| hydrophobic | 0.075 |

**GATE 2 fails.** Random 0.5925, null 99th pct **+0.0248**:

| term | direction | gain | rho | p |
|---|---|---|---|---|
| **full Vina score** | low (better energy) | **+0.0205** | −0.083 | 0.028 |
| repulsion | low | +0.0034 | +0.062 | 0.33 |
| gauss1 | high | −0.0016 | +0.127 | 0.49 |
| hydrophobic | low | **−0.0380** | +0.157 | 0.999 |
| hbond | low | **−0.0461** | +0.087 | 1.000 |

The best term does **not** clear the null, and the individual physical terms are null or
strongly negative. The full score at least points the chemically correct way - lower
interaction energy is better - but at ρ = −0.083 and +0.0205 against a +0.0248 threshold it
is indistinguishable from the noise, and it is one of ten tests reported here.

### Decision

**The xTB installation is not justified.** The recommendation above was made on a prior;
it now rests on a measurement. A classical interaction energy has exactly the property
tier 1 lacked - within-ligand variance - and still cannot discriminate, which means the
limitation is not "the term is constant" but that **agreement with the pocket does not
separate good CYP3A4 poses from bad ones at this pool's resolution**.

**Honest caveat, because it cuts the other way.** This Vina implementation is approximate:
protein atom types are inferred from the first character of the atom name, there is no
torsional entropy term, and no desolvation. GFN2-xTB would add electrostatics,
polarisation and real strain, which is a genuine physical difference rather than a
refinement. So this is evidence against the *family*, not proof against xTB specifically.
What it removes is the justification for paying a cluster install up front: the cheap
member of the family should have shown *something*, and it showed +0.0205 on ten tries.

**Where the effort goes instead:** the P450 generalisation set, which accumulates for free
and has already produced the strongest result of the campaign.
