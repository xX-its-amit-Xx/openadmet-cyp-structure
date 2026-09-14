# FINDING 013 — a union pool adds oracle headroom that selection cannot reach

**Date:** 2026-09-14
**Pools:** 1,740 Boltz + 1,029 Protenix-v2 poses, 87 CYP3A4 ligands
**Selector:** cross-engine agreement, leak-free — each pose scored only against poses from
the *other* engine, so nothing is ever compared with itself.

---

## The result

| pool | oracle | random | **selected** |
|---|---|---|---|
| Boltz only | 0.6975 | 0.5784 | **0.6112** |
| union (Boltz + Protenix) | **0.7350** | 0.5874 | **0.6095** |
| difference | **+0.0375** | +0.0090 | **−0.0017** |

Adding a second engine's poses to the pool buys **+0.0375 of oracle** and selection
captures **none of it** — it goes very slightly backwards.

## This is the PXR lesson, reproduced with a mechanism

The inherited campaign note says Protenix beat the pool on 6 of 8 holo cases yet
*regressed* the board 0.5551 → 0.5241 as a full pool member, while an 8-ligand tail swap
gained +0.017 and won. That was recorded as a rule of thumb. It is now reproduced under
controlled conditions with the reason visible:

**The extra poses are genuinely better — the oracle proves it — but nothing available can
identify them.** Cross-engine agreement is the strongest selector this project has, and it
still cannot pick out Protenix's wins. The union's own random baseline rises (+0.0090),
which confirms the added poses are decent on average; the selected score does not follow.

## Consequences

- **Do not put a second engine into the submission pool.** Keep Protenix as the
  *reference* for cross-engine agreement, where it is worth +0.0381, not as a pool member,
  where it is worth −0.0017.
- **The prize got bigger, not smaller.** Against the union oracle of 0.7350, current
  selection at 0.6095 leaves **0.126 LDDT-PLI per ligand** unclaimed — more than the 0.127
  the project started from, because generation improved while selection did not keep pace.
- **It reframes what "more sampling" is worth.** FINDING 004 measured the selector
  tracking the oracle at +0.0125 per doubling. That was within one engine. Across engines
  the coupling is zero here, so "add another engine's poses" and "sample more deeply from
  one engine" are not interchangeable ways of buying oracle.

## What would change this

A selector that works on the *union* would collect the +0.0375 immediately. The obvious
candidate is a term that judges a pose without reference to its siblings at all — a
physics score, which is what `docs/QM_SCORER_DESIGN.md` was always for. Consensus methods
are structurally unable to do this: they rank poses by agreement within a set, and the
union's value lies precisely in poses that *disagree* with the majority engine.
