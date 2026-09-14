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
