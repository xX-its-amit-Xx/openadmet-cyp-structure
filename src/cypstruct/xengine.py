"""Cross-engine consensus selection — the first term measured to beat the incumbent.

FINDING 011: score a pose by how closely INDEPENDENT engines reproduce it. On 63 CYP3A4
ligands this scores +0.0284 against random (null 99th pct +0.0204, 0 of 3,000 draws
beating it) where the incumbent `0.5*z(contacts) - z(rmsd_to_siblings)` scores +0.0207.

Three properties that make it usable rather than merely significant:

* **No fitted parameters.** There is no weight to overfit and therefore no held-out fold
  to demand. Fitting a three-weight blend leave-one-cluster-out collapses it to +0.0149,
  inside the noise — so the shipped form is deliberately unweighted.
* **Permutation-invariant.** Engines order a ligand's atoms differently, and index-for-
  index comparison is the trap that once halved every LDDT-PLI. Chamfer distance needs no
  correspondence at all.
* **Superposition-free.** Each pose is placed in its OWN heme frame, so there is no
  protein alignment and no reference frame to get wrong.

**It needs ~4 independent poses per ligand, and it needs them to be genuinely
independent.** Below that the estimate is noise: with one reference pose the same feature
measured −0.0055. Verify with `reference_depth()` before trusting any selection, because
`diffusion_samples` does NOT produce independent poses on OpenProtein's Protenix
(FINDING 009) — replicate jobs do.

**protenix-v1 is DETERMINISTIC — buy exactly one replicate of it, ever.** Measured at
per-atom sd 0.0000 Å across replicates on 6 of 6 standard ligands and 14 of 14
organometallics. It contributes exactly one independent pose per ligand no matter how many
jobs you run, and that one pose is worth having (adding it to 8 Protenix-v2 poses moved
+0.0240 → +0.0310) — but every replicate after the first is pure waste. Protenix-v2 on the
same ligands is properly diverse (median sd 1.29 Å).

**Do NOT pass every engine you have.** Measured at matched depth, the best reference set
is the two Protenix checkpoints alone (+0.0380). Adding esmfold2 — a genuinely different
architecture, sampling properly, at depth 4 — drops it to +0.0178, and all four engines
give +0.0267. Reference quality beats reference variety: a pose you should disagree with
costs signal when it is averaged into the mean.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

__all__ = ["heme_frame", "in_heme_frame", "chamfer", "reference_poses",
           "reference_depth", "xeng_score", "select"]

_FNAME = re.compile(r"(.+)__r(\d+)s(\d+)\.cif$")


def heme_frame(heme_xyz, heme_elem, axial_sg):
    """(Fe, distal normal, in-plane x) from heme atoms — computable from a prediction.

    The plane comes from the four pyrrole nitrogens by SVD rather than a cross product of
    two of them: a ruffled porphyrin tilts any single pair by degrees. The normal is
    flipped away from the proximal thiolate so "distal" is a fact, not a sign convention.
    The x-axis uses the propionate oxygens — the only O atoms in a heme — because the
    porphyrin's 4-fold pseudo-symmetry leaves any nitrogen-derived axis arbitrary.
    """
    heme_xyz = np.asarray(heme_xyz, float)
    if len(heme_xyz) == 0:
        return None
    fe, ns, ox = None, [], []
    for p, e in zip(heme_xyz, heme_elem):
        e = str(e).upper()
        if e == "FE":
            fe = p
        elif e == "N":
            ns.append(p)
        elif e == "O":
            ox.append(p)
    if fe is None or len(ns) < 4 or not ox:
        return None
    ns = sorted(ns, key=lambda p: np.linalg.norm(p - fe))[:4]
    pm = np.array(ns) - np.array(ns).mean(axis=0)
    n = np.linalg.svd(pm)[2][-1]
    n /= np.linalg.norm(n)
    if axial_sg is not None and np.dot(n, np.asarray(axial_sg, float) - fe) > 0:
        n = -n
    v = np.mean(ox, axis=0) - fe
    v = v - n * (v @ n)
    if np.linalg.norm(v) < 1e-6:
        return None
    return fe, n, v / np.linalg.norm(v)


def in_heme_frame(complex_) -> np.ndarray | None:
    """Ligand heavy atoms expressed in the complex's own heme frame."""
    fr = heme_frame(complex_.heme_xyz, complex_.heme_elem, complex_.axial_sg)
    if fr is None or len(complex_.lig_xyz) == 0:
        return None
    fe, n, x = fr
    y = np.cross(n, x)
    d = np.asarray(complex_.lig_xyz, float) - fe
    return np.column_stack([d @ x, d @ y, d @ n])


def chamfer(a: np.ndarray, b: np.ndarray) -> float:
    """Symmetric nearest-neighbour distance. No atom correspondence required."""
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    return float(0.5 * (d.min(axis=1).mean() + d.min(axis=0).mean()))


def reference_poses(pool_dirs, loader) -> dict[str, list[np.ndarray]]:
    """One pose per (engine, replicate) — the independent opinions about each ligand.

    Keyed on the replicate index, never on file order. Replicate 0 of an OpenProtein
    Protenix job contributes 20 models that share ONE ligand conformation (FINDING 009),
    so a file-order slice would return N copies of a single pose and present them as N
    independent opinions.
    """
    ref: dict[str, dict[tuple[str, int], np.ndarray]] = {}
    for pool in pool_dirs:
        pool = Path(pool)
        engine = pool.name
        for f in sorted(pool.glob("*__r*.cif")):
            m = _FNAME.match(f.name)
            if m is None:
                continue
            lig, rep = m.group(1), int(m.group(2))
            key = (engine, rep)
            if key in ref.get(lig, {}):
                continue
            try:
                v = in_heme_frame(loader(f))
            except Exception:
                continue
            if v is not None:
                ref.setdefault(lig, {})[key] = v
    return {lig: _dedupe(list(d.values())) for lig, d in ref.items()}


def _dedupe(poses: list[np.ndarray], tol: float = 0.05) -> list[np.ndarray]:
    """Drop poses identical to one already kept. A replicate is not an opinion.

    Being a separate job is not sufficient for independence. RoseTTAFold-3 in
    single-sequence mode is *sometimes* deterministic - measured across replicates of one
    ligand at per-atom sd 0.64 and 1.59 A, and of another at exactly 0.0000 - so the
    replicate count overstates how many independent opinions there really are, and
    `xeng_score` would then average a duplicate in twice and weight it double.

    Counting distinct poses rather than distinct jobs is the general form of the FINDING
    009 lesson, and it costs one comparison per pair.
    """
    keep: list[np.ndarray] = []
    for p in poses:
        if not any(len(q) == len(p) and float(np.abs(q - p).max()) < tol for q in keep):
            keep.append(p)
    return keep


def reference_depth(ref: dict[str, list[np.ndarray]]) -> dict:
    """How many independent opinions per ligand, and is it enough to trust the score?"""
    depth = {k: len(v) for k, v in ref.items()}
    if not depth:
        return {"ligands": 0, "usable": False}
    vals = np.array(list(depth.values()))
    return {"ligands": len(depth), "min": int(vals.min()),
            "median": float(np.median(vals)), "max": int(vals.max()),
            "at_least_4": int((vals >= 4).sum()),
            "usable": bool((vals >= 4).mean() > 0.5),
            "note": ("below 4 independent poses the feature measured -0.0055; "
                     "verify replicates are not identical (FINDING 009)")}


def xeng_score(pose_xyz_in_frame: np.ndarray, refs: list[np.ndarray]) -> float:
    """Mean Chamfer distance to the independent opinions. LOWER is better."""
    if not refs:
        return float("nan")
    return float(np.mean([chamfer(pose_xyz_in_frame, r) for r in refs]))


def select(df, ligand_col="ligand", xeng_col="xeng", sibling_rmsd_col=None):
    """Pick one pose per ligand. Unweighted by design — fitting weights loses the gain.

    Default is `-z(xeng)` ALONE, which at full depth scores +0.0336 against the incumbent
    `0.5*zc - zm` at +0.0264. Passing `sibling_rmsd_col` adds `-z(sibling rmsd)` and makes
    it WORSE (+0.0183): once ~10 independent opinions from other engines are available,
    asking whether one engine agrees with itself adds nothing. The combination looked best
    at n=63 (+0.0305) and did not replicate - it is kept only as an option, not a default.
    Terms are z-scored WITHIN the ligand, since selection only compares poses of the same
    molecule.
    """
    import pandas as pd

    d = pd.DataFrame(df).copy()

    def zc(col):
        g = d.groupby(ligand_col)[col]
        return (d[col] - g.transform("mean")) / g.transform("std").replace(0, 1)

    score = -zc(xeng_col)
    if sibling_rmsd_col is not None:
        score = score - zc(sibling_rmsd_col)
    d["_score"] = score
    return d.loc[d.groupby(ligand_col)["_score"].idxmax()].drop(columns="_score")
