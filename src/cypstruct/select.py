"""Pose selection: the cross-model z-hybrid, and the tail-rescue rule that beat it.

Ported from the PXR structure campaign, where this recipe finished 2nd of ~50 teams at
0.5640 LDDT-PLI. Two things carry over, and the second one is the non-obvious one.

**1. The z-hybrid.** Per engine, pick the best sample by that engine's *native* confidence
signal; z-score that signal across all ligands *within that engine*; then per ligand take
the engine with the highest z. Cross-model diversity works because engines fail on
different ligands; z-normalisation is what stops one engine's inflated raw scale from
winning everything.

**2. A new engine joins as a TAIL RESCUER, not as a pool member.** Measured in PXR:
Protenix beat the pool on 6 of 8 holo ground-truth cases, yet adding it as a full z-hybrid
member *regressed* the leaderboard from 0.5551 to 0.5241. Swapping it in for only the 8
ligands with the lowest pool confidence gained **+0.017** and produced the best submission
of the campaign. Bracket: N=4 -> 0.5578, **N=8 -> 0.5640**, N=12 -> 0.5629, N=20 -> 0.5587.

The mechanism is that a strong-but-different engine helps where the incumbent is unsure and
hurts where it is confident. Averaging that into one pool destroys the distinction.

Ties: the PXR implementation broke them positionally by argument order, silently. Here the
tie-break is explicit and deterministic.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# Per-engine confidence field and sign. Positive means higher-is-better AFTER the sign is
# applied, so every engine ends up on the same footing before z-scoring.
#
# The sign flips are not cosmetic: PDE and PAE are ERROR estimates, so a raw max() over
# them selects the *worst* pose. That inversion is an easy and completely silent bug.
CONFIDENCE_FIELDS: dict[str, tuple[str, float]] = {
    "boltz2":    ("complex_ipde", -1.0),      # interface predicted distance error
    "boltz2_aff": ("affinity_pred_value", 1.0),
    "af3":       ("ranking_score", 1.0),
    "openfold3": ("pae_pocket_mean", -1.0),   # mean protein-ligand PAE over the pocket
    "chai":      ("aggregate_score", 1.0),
    "protenix":  ("ranking_score", 1.0),
}


@dataclass
class Pose:
    ligand: str
    engine: str
    path: str
    raw_confidence: float          # already sign-corrected: higher is better
    seed: int | None = None
    sample: int | None = None
    z: float = math.nan
    extra: dict | None = None


def best_within_engine(poses: list[Pose]) -> dict[tuple[str, str], Pose]:
    """For each (ligand, engine), the single best sample by that engine's own confidence."""
    best: dict[tuple[str, str], Pose] = {}
    for p in poses:
        k = (p.ligand, p.engine)
        cur = best.get(k)
        if cur is None or p.raw_confidence > cur.raw_confidence:
            best[k] = p
    return best


def zscore_within_engine(best: dict[tuple[str, str], Pose]) -> dict[tuple[str, str], Pose]:
    """Z-score each engine's confidences over the ligands THAT ENGINE actually covers.

    Normalising over the engine's own coverage (rather than over the union of ligands) is
    deliberate: an engine that only produced poses for the easy half of the set would
    otherwise be compared on a scale built from a different population than it competes on.
    """
    by_engine: dict[str, list[Pose]] = {}
    for p in best.values():
        by_engine.setdefault(p.engine, []).append(p)
    for _eng, ps in by_engine.items():
        v = np.array([p.raw_confidence for p in ps], float)
        mu = float(v.mean())
        sd = float(v.std())
        if sd == 0.0 or not np.isfinite(sd):
            sd = 1.0                       # zero-variance guard; all poses tie at z=0
        for p in ps:
            p.z = (p.raw_confidence - mu) / sd
    return best


def z_hybrid(poses: list[Pose], engine_priority: list[str] | None = None
             ) -> tuple[dict[str, Pose], dict[str, int]]:
    """Cross-model selection. Returns ({ligand: chosen pose}, {engine: n_wins}).

    `engine_priority` makes tie-breaking explicit and reproducible. Without it, ties fall
    to whichever engine happened to be enumerated first, which makes the result depend on
    dict ordering rather than on anything meaningful.
    """
    best = zscore_within_engine(best_within_engine(poses))
    by_lig: dict[str, list[Pose]] = {}
    for p in best.values():
        by_lig.setdefault(p.ligand, []).append(p)

    prio = {e: i for i, e in enumerate(engine_priority or sorted(CONFIDENCE_FIELDS))}
    chosen, wins = {}, {}
    for lig, cands in by_lig.items():
        pick = min(cands, key=lambda p: (-p.z, prio.get(p.engine, 999), p.engine, p.path))
        chosen[lig] = pick
        wins[pick.engine] = wins.get(pick.engine, 0) + 1
    return chosen, wins


def pool_confidence(poses: list[Pose]) -> dict[str, float]:
    """Per ligand, the best z across engines — how sure the pool as a whole is.

    This is the ordering used to find the tail. It is a *relative* measure: a low value
    means every engine was unusually unsure about this ligand compared with the rest of
    the set, which is exactly the situation where a different engine is worth trying.
    """
    best = zscore_within_engine(best_within_engine(poses))
    out: dict[str, float] = {}
    for p in best.values():
        out[p.ligand] = max(out.get(p.ligand, -math.inf), p.z)
    return out


def tail_rescue(base: dict[str, Pose], rescue: dict[str, Pose],
                pool_conf: dict[str, float], n: int = 8
                ) -> tuple[dict[str, Pose], list[str]]:
    """Replace the `n` least-confident ligands' poses with a rescue engine's.

    This is the rule that actually won PXR, and the rule that a naive reading of
    "more models is better" gets exactly backwards. `n` is a real hyperparameter with an
    interior optimum — sweep it against held-out ground truth, do not assume 8 transfers.
    """
    eligible = [lig for lig in sorted(pool_conf, key=lambda x: pool_conf[x])
                if lig in rescue and lig in base]
    swap = eligible[:n]
    out = dict(base)
    for lig in swap:
        out[lig] = rescue[lig]
    return out, swap


def sweep_tail_rescue(base: dict[str, Pose], rescue: dict[str, Pose],
                      pool_conf: dict[str, float], truth: dict[str, float],
                      ns: tuple[int, ...] = (0, 2, 4, 6, 8, 12, 16, 20, 30)
                      ) -> list[dict]:
    """Score the rescue depth against ground truth. `truth` maps pose path -> LDDT-PLI.

    Returns one row per n so the interior optimum (or its absence) is visible rather than
    assumed. If the curve is flat, that is the finding — say so instead of picking a peak
    out of noise.
    """
    rows = []
    for n in ns:
        merged, swapped = tail_rescue(base, rescue, pool_conf, n)
        vals = [truth[p.path] for p in merged.values() if p.path in truth]
        if not vals:
            continue
        rows.append({"n": n, "n_swapped": len(swapped), "n_scored": len(vals),
                     "mean_lddt_pli": round(float(np.mean(vals)), 4),
                     "median_lddt_pli": round(float(np.median(vals)), 4)})
    return rows


def consensus_clusters(poses: list[Pose], coords: dict[str, np.ndarray],
                       cutoff: float = 2.0) -> dict[str, list[list[str]]]:
    """Cluster a ligand's poses by mutual RMSD. Returns {ligand: [[path, ...], ...]}.

    Relevant because OpenADMET reports that for several CYP3A4 ligands the cryoEM density
    fits **multiple mutually exclusive conformations**. When a pool splits into two tight,
    well-populated clusters, that is a signal about the system rather than noise to be
    averaged away — and picking the larger, better-scoring cluster beats picking a lone
    high-confidence outlier.
    """
    by_lig: dict[str, list[Pose]] = {}
    for p in poses:
        by_lig.setdefault(p.ligand, []).append(p)
    out: dict[str, list[list[str]]] = {}
    for lig, ps in by_lig.items():
        ps = [p for p in ps if p.path in coords]
        if not ps:
            continue
        n = len(ps)
        assigned = [-1] * n
        clusters: list[list[int]] = []
        for i in range(n):
            if assigned[i] >= 0:
                continue
            cid = len(clusters)
            clusters.append([i])
            assigned[i] = cid
            for j in range(i + 1, n):
                if assigned[j] >= 0:
                    continue
                a, b = coords[ps[i].path], coords[ps[j].path]
                if a.shape != b.shape:
                    continue
                if float(np.sqrt(((a - b) ** 2).sum(1).mean())) <= cutoff:
                    assigned[j] = cid
                    clusters[cid].append(j)
        out[lig] = [[ps[i].path for i in c] for c in
                    sorted(clusters, key=len, reverse=True)]
    return out
