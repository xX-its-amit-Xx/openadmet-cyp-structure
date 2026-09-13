"""An orientation prior in the FULL heme frame, including the azimuth.

The cylindrical occupancy prior gained +0.0018 against a +0.0142 noise floor - a clean
null. The reason it failed is worth more than the attempt: `(r, z)` is exactly the
marginal that FINDING 008's addendum already showed is saturated. Every pose coordinates
the iron and sits in the pocket, so every pose has roughly the right radius and height.
Averaging over the azimuth threw away the only dimension in which the poses actually
differ.

So this restores it. The frame is completed with an in-plane x-axis taken from the heme's
propionate oxygens - the only O atoms in a heme, so the direction is fixed by chemistry
rather than by atom names, and the porphyrin's 4-fold pseudo-symmetry no longer leaves the
azimuthal zero arbitrary.

It also scores the right object. Per-atom histograms over 3 dimensions are hopelessly
sparse at 17.5k atoms; instead each pose is reduced to a handful of numbers describing
**where the bulk of the molecule sits and which way it points**, and those are scored
under a KDE fitted to the crystal poses:

    r, z, cos(phi), sin(phi)   centroid placement in the heme frame
    tilt                       angle between the ligand's principal axis and the normal
    r_gyr                      how extended the pose is

Held out exactly as before: the KDE is fitted on **non-CYP3A4** P450s and tested on CYP3A4
predictions, so nothing about a CYP3A4 ligand can enter its own score.

**Stated in advance, so this cannot be graded after the fact:** the bar is FINDING 007's -
a gain over random above the 99th percentile of a random-feature null measured on this
same pool. Anything below +0.020 is noise and gets logged as a negative.

    python scripts/structure/orientation_prior.py build
    python scripts/structure/orientation_prior.py test
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

UNI = DATA_PROCESSED / "p450_universe"
MODEL = UNI / "orientation_prior.npz"


def frame(heme_xyz, heme_elem, sg):
    """(Fe, normal, x-axis) from heme atoms alone - all computable from a prediction."""
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
    Pm = np.array(ns) - np.array(ns).mean(axis=0)
    n = np.linalg.svd(Pm)[2][-1]
    n /= np.linalg.norm(n)
    if sg is not None and np.dot(n, np.asarray(sg, float) - fe) > 0:
        n = -n
    v = np.mean(ox, axis=0) - fe
    v = v - n * (v @ n)
    if np.linalg.norm(v) < 1e-6:
        return None
    return fe, n, v / np.linalg.norm(v)


def descriptors(lig_xyz, fe, n, x) -> np.ndarray | None:
    """Six numbers describing where the bulk sits and which way it points."""
    L = np.asarray(lig_xyz, float)
    if len(L) < 3:
        return None
    y = np.cross(n, x)
    d = L - fe
    z = d @ n
    a, b = d @ x, d @ y
    r = np.sqrt(np.maximum(a * a + b * b, 0))

    cz, ca, cb = z.mean(), a.mean(), b.mean()
    cr = float(np.hypot(ca, cb))
    phi = float(np.arctan2(cb, ca))

    # principal axis of the pose, and how far it leans off the heme normal. Sign is
    # arbitrary for an eigenvector, so fold to [0, 90] deg.
    C = L - L.mean(axis=0)
    ax = np.linalg.svd(C)[2][0]
    tilt = float(np.degrees(np.arccos(abs(float(ax @ n)))))
    rg = float(np.sqrt((C * C).sum(1).mean()))
    return np.array([cr, float(cz), np.cos(phi), np.sin(phi), tilt, rg])


def cmd_build(exclude: str = "P08684") -> dict:
    atoms = pd.read_parquet(UNI / "p450_atoms.parquet")
    geo = pd.read_parquet(UNI / "p450_geometry.parquet")
    key = ["pdb", "chain", "lig", "seqid"]
    m = atoms.merge(geo[key + ["closest_fe"]], on=key)
    m = m[m.closest_fe <= 6.5]
    cset = pd.read_csv(UNI / "p450_cofold_set.csv")
    drop = set(cset[cset.uniprot == exclude].pdb)
    m = m[~m.pdb.isin(drop)]

    X, n_no_axis = [], 0
    for row in m.itertuples():
        if not row.xaxis:
            n_no_axis += 1
            continue
        xyz = np.frombuffer(row.xyz, np.float32).reshape(-1, 3).astype(float)
        fe = np.frombuffer(row.fe, np.float32).astype(float)
        n = np.frombuffer(row.normal, np.float32).astype(float)
        x = np.frombuffer(row.xaxis, np.float32).astype(float)
        v = descriptors(xyz, fe, n, x)
        if v is not None:
            X.append(v)
    X = np.array(X)
    np.savez(MODEL, X=X, excluded=exclude)
    return {"crystal_poses": len(X), "skipped_no_axis": n_no_axis,
            "excluded": exclude,
            "descriptor_means": [round(float(v), 2) for v in X.mean(0)]}


def _kde_logpdf(train: np.ndarray, query: np.ndarray) -> np.ndarray:
    from scipy.stats import gaussian_kde

    mu, sd = train.mean(0), train.std(0) + 1e-9
    k = gaussian_kde(((train - mu) / sd).T)
    return k.logpdf(((query - mu) / sd).T)


def cmd_test(pool_root: str, n_null: int = 2000) -> dict:
    from scipy import stats

    from cypstruct import pose as P

    d = np.load(MODEL)
    train = d["X"]
    sc = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    sc = sc[sc.arm == "unsteered"].copy()

    rows = []
    for dirp in sorted(Path(pool_root).glob("*__unsteered__*")):
        lig = dirp.name.split("__")[0]
        for f in sorted(dirp.glob("input_model_*.cif")):
            try:
                mo = P.load_structure(f)
                fr = frame(mo.heme_xyz, mo.heme_elem, mo.axial_sg)
                if fr is None or len(mo.lig_xyz) == 0:
                    continue
                v = descriptors(mo.lig_xyz, *fr)
                if v is not None:
                    rows.append({"ligand": lig, "sample": f.stem, "v": v})
            except Exception:
                continue
    if not rows:
        return {"error": "no poses read", "root": pool_root}

    Q = np.array([r["v"] for r in rows])
    lp = _kde_logpdf(train, Q)
    f = pd.DataFrame({"ligand": [r["ligand"] for r in rows],
                      "sample": [r["sample"] for r in rows], "orient": lp})
    m = sc.merge(f, on=["ligand", "sample"], how="inner").dropna(subset=["orient"])
    if m.empty:
        return {"error": "join produced nothing"}

    rs = []
    for lig, g in m.groupby("ligand"):
        if g.orient.nunique() >= 3:
            rs.append(stats.spearmanr(g.orient, g.lddt_pli).statistic)
    rs = np.array([r for r in rs if np.isfinite(r)])

    rng = np.random.default_rng(0)
    rand = float(np.mean([
        m.groupby("ligand").lddt_pli.apply(lambda s: s.sample(1, random_state=i).iloc[0])
        .mean() for i in range(200)]))
    oracle = float(m.groupby("ligand").lddt_pli.max().mean())
    selected = float(m.loc[m.groupby("ligand").orient.idxmax()].lddt_pli.mean())
    gain = selected - rand

    null = []
    for _ in range(n_null):
        mm = m.assign(s=rng.normal(size=len(m)))
        null.append(float(mm.loc[mm.groupby("ligand").s.idxmax()].lddt_pli.mean()) - rand)
    null = np.array(null)

    return {
        "poses": len(m), "ligands": int(m.ligand.nunique()),
        "prior_from": int(len(train)), "excluded": str(d["excluded"]),
        "within_ligand_rho_mean": round(float(rs.mean()), 4),
        "frac_rho_positive": round(float((rs > 0).mean()), 3),
        "wilcoxon_p": round(float(stats.wilcoxon(rs).pvalue), 5),
        "random": round(rand, 4), "oracle": round(oracle, 4),
        "selected": round(selected, 4), "gain_vs_random": round(gain, 4),
        "null_p95": round(float(np.percentile(null, 95)), 4),
        "null_p99": round(float(np.percentile(null, 99)), 4),
        "empirical_p": round(float((null >= gain).mean()), 4),
        "verdict": ("REAL - clears the 99th pct of the null"
                    if gain > np.percentile(null, 99) else
                    "inside the noise floor - log as a negative"),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "test"])
    ap.add_argument("--exclude", default="P08684")
    ap.add_argument("--pool", default="D:/cyp_scratch/val87b_unsteered")
    a = ap.parse_args()
    print(json.dumps(cmd_build(a.exclude) if a.cmd == "build"
                     else cmd_test(a.pool), indent=2))
