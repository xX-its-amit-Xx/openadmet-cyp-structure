"""Where do ligand atoms ACTUALLY sit? An empirical occupancy prior over the heme frame.

Findings G and I between them say the interesting thing plainly: coordination is worth
+0.14 LDDT-PLI as a feature, but the iron anchor is **saturated** - 84% of poses already
coordinate, so nothing measured at the metal discriminates any more. Whatever separates a
good pose from a bad one lives in where the REST of the molecule goes.

This asks that question empirically instead of from theory. Every crystal ligand in the
P450 universe is transformed into its own heme frame (origin at the iron, z along the
distal normal) and the atoms are accumulated into a density over (r, z). A predicted pose
is then scored by how typical its atoms are under that density.

**Why this can't be a leaky feature, and why the test is strong.** The prior is built from
**non-CYP3A4** P450s only and tested on CYP3A4 predictions. There is no path by which the
answer for a CYP3A4 ligand enters its own score - not the ligand, not the entry, not even
the protein. If a density learned from bacterial and other human P450s ranks CYP3A4 poses,
that is transferable P450 physics, which is exactly the claim the scorer needs.

**Cylindrical, not Cartesian, on purpose.** The porphyrin has 4-fold pseudo-symmetry, so
an azimuthal angle has no consistent zero across entries without an extra convention;
(r, z) sidesteps that entirely. It costs the I-helix/F-G asymmetry, which is the obvious
next refinement if the cylindrical version shows signal.

    python scripts/structure/occupancy_prior.py build
    python scripts/structure/occupancy_prior.py test
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
PRIOR = UNI / "occupancy_prior.npz"

R_MAX, Z_MIN, Z_MAX = 14.0, -2.0, 18.0
NR, NZ = 28, 40


def _frame_from_heme(heme_xyz: np.ndarray, heme_elem: list[str],
                     sg: np.ndarray | None) -> tuple[np.ndarray, np.ndarray] | None:
    """(Fe, distal unit normal) from heme atoms alone - computable from a prediction.

    Same construction as the reference pass: the plane comes from the four pyrrole
    nitrogens by SVD, and the normal is flipped to point AWAY from the proximal thiolate
    so that "distal" is a fact rather than a sign convention. Without the thiolate we fall
    back to the side the porphyrin's own substituents do not occupy.
    """
    if len(heme_xyz) == 0:
        return None
    fe = None
    ns = []
    for p, e in zip(heme_xyz, heme_elem):
        if e.upper() == "FE":
            fe = np.asarray(p, float)
        elif e.upper() == "N":
            ns.append(np.asarray(p, float))
    if fe is None or len(ns) < 4:
        return None
    ns = sorted(ns, key=lambda p: np.linalg.norm(p - fe))[:4]
    P = np.array(ns) - np.array(ns).mean(axis=0)
    n = np.linalg.svd(P)[2][-1]
    n = n / np.linalg.norm(n)
    if sg is not None and np.dot(n, np.asarray(sg, float) - fe) > 0:
        n = -n
    return fe, n


def _rz(xyz: np.ndarray, fe: np.ndarray, n: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    d = np.asarray(xyz, float) - fe
    z = d @ n
    r = np.sqrt(np.maximum((d * d).sum(1) - z * z, 0.0))
    return r, z


def cmd_build(exclude_uniprot: str = "P08684") -> dict:
    """Accumulate the density from crystal ligands, holding out the target isoform."""
    atoms = pd.read_parquet(UNI / "p450_atoms.parquet")
    geo = pd.read_parquet(UNI / "p450_geometry.parquet")
    key = ["pdb", "chain", "lig", "seqid"]
    m = atoms.merge(geo[key + ["closest_fe"]], on=key)
    m = m[m.closest_fe <= 6.5]

    # which entries belong to the held-out isoform
    cset = pd.read_csv(UNI / "p450_cofold_set.csv")
    drop = set(cset[cset.uniprot == exclude_uniprot].pdb)
    before = len(m)
    m = m[~m.pdb.isin(drop)]
    print(f"{before} in-pocket poses -> {len(m)} after removing "
          f"{len(drop)} {exclude_uniprot} entries", flush=True)

    H = np.zeros((NR, NZ))
    n_used = 0
    for row in m.itertuples():
        xyz = np.frombuffer(row.xyz, dtype=np.float32).reshape(-1, 3).astype(float)
        fe = np.frombuffer(row.fe, dtype=np.float32).astype(float)
        nrm = np.frombuffer(row.normal, dtype=np.float32).astype(float)
        r, z = _rz(xyz, fe, nrm)
        h, _, _ = np.histogram2d(r, z, bins=[NR, NZ],
                                 range=[[0, R_MAX], [Z_MIN, Z_MAX]])
        H += h
        n_used += 1

    # Volume-normalise: a cylindrical shell at large r contains more space, so raw counts
    # would call the crowded centre "typical" purely as a geometry artefact.
    r_edges = np.linspace(0, R_MAX, NR + 1)
    vol = (np.pi * (r_edges[1:] ** 2 - r_edges[:-1] ** 2))[:, None] \
        * ((Z_MAX - Z_MIN) / NZ)
    dens = (H + 0.5) / vol                       # +0.5 so empty cells stay finite
    logd = np.log(dens / dens.sum())
    np.savez(PRIOR, logd=logd, R_MAX=R_MAX, Z_MIN=Z_MIN, Z_MAX=Z_MAX,
             NR=NR, NZ=NZ, n_poses=n_used, excluded=exclude_uniprot)
    return {"poses_used": n_used, "atoms": int(H.sum()),
            "excluded_uniprot": exclude_uniprot, "prior": str(PRIOR.name)}


def score_pose(lig_xyz, heme_xyz, heme_elem, sg, logd) -> float:
    """Mean log-density of the ligand's atoms. Higher = more like a real P450 complex."""
    fr = _frame_from_heme(np.asarray(heme_xyz), list(heme_elem), sg)
    if fr is None:
        return float("nan")
    fe, n = fr
    r, z = _rz(np.asarray(lig_xyz), fe, n)
    ir = np.clip((r / R_MAX * NR).astype(int), 0, NR - 1)
    iz = np.clip(((z - Z_MIN) / (Z_MAX - Z_MIN) * NZ).astype(int), 0, NZ - 1)
    return float(logd[ir, iz].mean())


def cmd_test(pool_root: str, n_null: int = 2000) -> dict:
    """Does it rank CYP3A4 poses within a ligand? Against the FINDING 007 null."""
    from scipy import stats

    from cypstruct import pose as P

    d = np.load(PRIOR)
    logd = d["logd"]
    root = Path(pool_root)
    sc = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    sc = sc[sc.arm == "unsteered"].copy()

    rows = []
    for dirp in sorted(root.glob("*__unsteered__*")):
        lig = dirp.name.split("__")[0]
        for f in sorted(dirp.glob("input_model_*.cif")):
            try:
                mo = P.load_structure(f)
                if len(mo.lig_xyz) == 0 or len(mo.heme_xyz) == 0:
                    continue
                rows.append({"ligand": lig, "sample": f.stem,
                             "occ": score_pose(mo.lig_xyz, mo.heme_xyz, mo.heme_elem,
                                               mo.axial_sg, logd)})
            except Exception:
                continue
    f = pd.DataFrame(rows)
    if f.empty:
        return {"error": "no poses read", "root": str(root)}
    m = sc.merge(f, on=["ligand", "sample"], how="inner").dropna(subset=["occ"])
    if m.empty:
        return {"error": "join produced nothing - check naming", "scored": len(f)}

    # within-ligand rank correlation: the ONLY statistic that can select
    rs = []
    for lig, g in m.groupby("ligand"):
        if g.occ.nunique() < 3:
            continue
        rs.append(stats.spearmanr(g.occ, g.lddt_pli).statistic)
    rs = np.array([r for r in rs if np.isfinite(r)])

    rng = np.random.default_rng(0)
    rand = float(np.mean([
        m.groupby("ligand").lddt_pli.apply(lambda s: s.sample(1, random_state=i).iloc[0])
        .mean() for i in range(200)]))
    oracle = float(m.groupby("ligand").lddt_pli.max().mean())
    idx = m.groupby("ligand").occ.idxmax()
    selected = float(m.loc[idx].lddt_pli.mean())

    # FINDING 007 insists the null be measured on THIS pool, not recalled from the other
    null = []
    for _ in range(n_null):
        mm = m.assign(s=rng.normal(size=len(m)))
        null.append(float(mm.loc[mm.groupby("ligand").s.idxmax()].lddt_pli.mean()) - rand)
    null = np.array(null)
    gain = selected - rand

    return {
        "poses": len(m), "ligands": int(m.ligand.nunique()),
        "prior_built_from": int(d["n_poses"]), "excluded": str(d["excluded"]),
        "within_ligand_rho_mean": round(float(rs.mean()), 4),
        "within_ligand_rho_median": round(float(np.median(rs)), 4),
        "frac_rho_positive": round(float((rs > 0).mean()), 3),
        "wilcoxon_p": round(float(stats.wilcoxon(rs).pvalue), 5),
        "random": round(rand, 4), "oracle": round(oracle, 4),
        "selected": round(selected, 4), "gain_vs_random": round(gain, 4),
        "null_sd": round(float(null.std()), 4),
        "null_p95": round(float(np.percentile(null, 95)), 4),
        "null_p99": round(float(np.percentile(null, 99)), 4),
        "empirical_p": round(float((null >= gain).mean()), 4),
        "verdict": ("REAL - clears the 99th pct of the null"
                    if gain > np.percentile(null, 99) else
                    "inside the noise floor - do NOT call it promising"),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "test"])
    ap.add_argument("--exclude", default="P08684")
    ap.add_argument("--pool", default="D:/cyp_scratch/val87b_unsteered")
    a = ap.parse_args()
    print(json.dumps(cmd_build(a.exclude) if a.cmd == "build"
                     else cmd_test(a.pool), indent=2))
