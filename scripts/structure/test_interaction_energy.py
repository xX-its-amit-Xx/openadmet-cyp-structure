"""The gate on tier 2: does a classical interaction energy discriminate poses at all?

FINDING 013 says only a sibling-free term can collect the +0.0375 of oracle a union pool
adds, and FINDING 014 redirects that to tier-2 xTB. But xTB needs a from-scratch cluster
install, and the cheap members of its family are *negative* - `max_clash` -0.0359,
`n_contacts` -0.0122 - so the prior is poor. This is the bounded experiment that decides
whether the install is worth it.

It implements the **AutoDock Vina** scoring terms, which are a real published function
rather than a proxy invented here: two gaussians on surface distance, a repulsion term, a
hydrophobic term and a hydrogen-bond term. Every term is per pose and needs no siblings,
so it has within-ligand variance by construction - the thing tier 1 turned out to lack.

Two gates, in order, and the second only matters if the first passes:

1. **within-ligand CV** - a term constant inside a ligand cannot rank that ligand's poses,
   which is how `is_coordinated` (CV 0.023) and the tier-1 donor prior both died.
2. **held-out correlation and selection gain** against the random-feature null measured on
   this same pool.

    python scripts/structure/test_interaction_energy.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

# Vina's published weights and radii. Distances are SURFACE distances d = r - R_i - R_j.
W = {"gauss1": -0.0356, "gauss2": -0.00516, "repulsion": 0.840,
     "hydrophobic": -0.0351, "hbond": -0.587}
RADII = {"C": 1.9, "N": 1.8, "O": 1.7, "S": 2.0, "P": 2.1, "F": 1.5,
         "CL": 1.8, "BR": 2.0, "I": 2.2, "FE": 1.2}
HYDROPHOBIC = {"C", "F", "CL", "BR", "I"}
POLAR = {"N", "O"}


def vina_terms(lig_xyz, lig_el, prot_xyz, prot_el, cutoff: float = 8.0) -> dict:
    """The five Vina terms for one pose. Pure geometry, no siblings, no ground truth."""
    L = np.asarray(lig_xyz, float)
    Pp = np.asarray(prot_xyz, float)
    if len(L) == 0 or len(Pp) == 0:
        return {k: np.nan for k in W}
    lr = np.array([RADII.get(str(e).upper(), 1.8) for e in lig_el])
    pr = np.array([RADII.get(str(e).upper(), 1.8) for e in prot_el])
    lh = np.array([str(e).upper() in HYDROPHOBIC for e in lig_el])
    ph = np.array([str(e).upper() in HYDROPHOBIC for e in prot_el])
    lp = np.array([str(e).upper() in POLAR for e in lig_el])
    pp_ = np.array([str(e).upper() in POLAR for e in prot_el])

    r = np.linalg.norm(L[:, None, :] - Pp[None, :, :], axis=2)
    near = r < cutoff
    if not near.any():
        return {k: 0.0 for k in W}
    d = r - lr[:, None] - pr[None, :]

    g1 = np.exp(-((d / 0.5) ** 2))
    g2 = np.exp(-(((d - 3.0) / 2.0) ** 2))
    rep = np.where(d < 0, d * d, 0.0)
    # hydrophobic: 1 below 0.5 A surface separation, ramping to 0 at 1.5 A
    hyd = np.clip(1.5 - d, 0, 1) * (lh[:, None] & ph[None, :])
    # hbond: 1 below -0.7 A (interpenetrating polar pair), ramping to 0 at 0
    hb = np.clip(-d / 0.7, 0, 1) * (lp[:, None] & pp_[None, :])
    m = near
    return {"gauss1": float(g1[m].sum()), "gauss2": float(g2[m].sum()),
            "repulsion": float(rep[m].sum()), "hydrophobic": float(hyd[m].sum()),
            "hbond": float(hb[m].sum())}


def main() -> int:
    from scipy import stats

    from cypstruct import pose as P

    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", default="D:/cyp_scratch/val87b_unsteered")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    sc = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    sc = sc[sc.arm == "unsteered"]
    dirs = sorted(Path(a.pool).glob("*__unsteered__*"))
    if a.limit:
        dirs = dirs[:a.limit]

    rows = []
    for i, dirp in enumerate(dirs):
        lig = dirp.name.split("__")[0]
        for f in sorted(dirp.glob("input_model_*.cif")):
            try:
                mo = P.load_structure(f)
            except Exception:
                continue
            if len(mo.lig_xyz) == 0:
                continue
            pel = [k[2][0] if isinstance(k, tuple) else "C" for k in mo.prot_key]
            t = vina_terms(mo.lig_xyz, mo.lig_elem, mo.prot_xyz, pel)
            t["vina"] = sum(W[k] * t[k] for k in W)
            rows.append({"ligand": lig, "sample": f.stem, **t})
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(dirs)}", flush=True)

    f = pd.DataFrame(rows)
    m = sc.merge(f, on=["ligand", "sample"]).dropna(subset=["vina"])
    print(f"\n{len(m)} poses, {m.ligand.nunique()} ligands")

    print("\nGATE 1 - within-ligand coefficient of variation (constant terms cannot rank):")
    for col in ["vina", "gauss1", "repulsion", "hydrophobic", "hbond"]:
        g = m.groupby("ligand")[col]
        cv = (g.transform("std") / g.transform("mean").abs().replace(0, np.nan))
        cvm = float(cv.groupby(m.ligand).first().median())
        print(f"  {col:12s} median within-ligand CV {cvm:.4f}"
              f"{'   <- too flat to rank' if cvm < 0.03 else ''}")

    rng = np.random.default_rng(0)
    rand = float(np.mean([
        m.groupby("ligand").lddt_pli.apply(lambda s: s.sample(1, random_state=i).iloc[0])
        .mean() for i in range(200)]))
    null = np.array([
        float(m.assign(s=rng.normal(size=len(m)))
              .pipe(lambda t: t.loc[t.groupby("ligand").s.idxmax()]).lddt_pli.mean()) - rand
        for _ in range(2000)])
    p99 = float(np.percentile(null, 99))
    print(f"\nGATE 2 - selection, random {rand:.4f}, null p99 {p99:+.4f}")
    for col in ["vina", "gauss1", "repulsion", "hydrophobic", "hbond"]:
        for sign, lbl in ((-1, "low"), (1, "high")):
            v = m.assign(s=sign * m[col])
            sel = float(v.loc[v.groupby("ligand").s.idxmax()].lddt_pli.mean())
            g = sel - rand
            rs = np.array([stats.spearmanr(gg[col], gg.lddt_pli).statistic
                           for _l, gg in m.groupby("ligand") if gg[col].nunique() >= 3])
            rs = rs[np.isfinite(rs)]
            flag = "  CLEARS p99" if g > p99 else ""
            print(f"  {col:12s} {lbl:4s} gain {g:+.4f}  rho {rs.mean():+.3f}"
                  f"  p={float((null >= g).mean()):.4f}{flag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
