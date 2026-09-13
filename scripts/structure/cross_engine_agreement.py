"""Does a Boltz pose that Protenix independently reproduces score better?

FINDING 005 asked whether a second engine can REPLACE the first on the ligands it gets
wrong, and the answer was no - Chai's errors correlate with Boltz's at rho +0.45, Protenix
at +0.474, and tail rescue does not transfer. That is a question about which engine to
trust per ligand.

This is a different question. Consensus is the one thing that has ever worked here:
`mean_rmsd_to_others` contributes +0.0100 to the only selector that beats the noise floor,
and FINDING 010 explains why - it is a WITHIN-ligand, ligand-specific comparison, which is
the only kind that can select. But that consensus is computed among samples of one engine,
so it can only measure whether Boltz agrees with itself. An architecturally different
model agreeing is strictly more information than the same model agreeing twice.

**The honest worry, stated first.** Because the two engines' errors are correlated
(rho +0.474), they may agree on being wrong, and cross-engine agreement could be no better
than self-consistency. That is exactly what this measures rather than assumes.

**Permutation-invariant by construction.** Two engines order a ligand's atoms differently,
and index-for-index comparison is the trap that halved every LDDT-PLI once already. So
agreement is a symmetric Chamfer distance between atom clouds - each atom to its nearest
neighbour in the other pose - which needs no correspondence at all. Both poses are placed
in their OWN heme frame first, so this is superposition-free: no protein alignment, and
therefore no reference frame to get wrong.

    python scripts/structure/cross_engine_agreement.py
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "structure"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402
from orientation_prior import frame  # noqa: E402

BOLTZ_POOL = Path("D:/cyp_scratch/val87b_unsteered")
PROT_POOL = DATA_PROCESSED / "openprotein" / "op1" / "protenix_v2"


def in_frame(mo) -> np.ndarray | None:
    """Ligand atoms expressed in the complex's own heme frame."""
    fr = frame(mo.heme_xyz, mo.heme_elem, mo.axial_sg)
    if fr is None or len(mo.lig_xyz) == 0:
        return None
    fe, n, x = fr
    y = np.cross(n, x)
    d = np.asarray(mo.lig_xyz, float) - fe
    return np.column_stack([d @ x, d @ y, d @ n])


def chamfer(a: np.ndarray, b: np.ndarray) -> float:
    """Symmetric nearest-neighbour distance. No atom correspondence required."""
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    return float(0.5 * (d.min(axis=1).mean() + d.min(axis=0).mean()))


def main(n_null: int = 2000) -> dict:
    from scipy import stats

    from cypstruct import pose as P

    # Protenix side, keyed by ligand
    # ONE pose per replicate, keyed on the replicate index. Taking the first N *files*
    # would be the FINDING 009 trap all over again: replicate 0 alone contributes 20
    # models that share a single ligand conformation, so a file-order slice would hand
    # back eight copies of one pose and call it eight independent opinions.
    byrep: dict[str, dict[int, np.ndarray]] = {}
    for f in sorted(PROT_POOL.glob("*__r*.cif")):
        mm = re.match(r"(.+)__r(\d+)s(\d+)\.cif$", f.name)
        if mm is None:
            continue
        lig, rep = mm.group(1), int(mm.group(2))
        if rep in byrep.get(lig, {}):
            continue
        try:
            v = in_frame(P.load_structure(f))
        except Exception:
            continue
        if v is not None:
            byrep.setdefault(lig, {})[rep] = v
    prot = {lig: list(d.values()) for lig, d in byrep.items()}
    depth = {k: len(v) for k, v in prot.items()}
    print(f"protenix ligands: {len(prot)}, independent poses each: "
          f"min {min(depth.values())} median {int(np.median(list(depth.values())))} "
          f"max {max(depth.values())}", flush=True)
    # the pre-registered re-test needs >= 4 independent poses to mean anything
    prot = {k: v for k, v in prot.items() if len(v) >= 4}
    print(f"ligands with >= 4 independent poses: {len(prot)}", flush=True)

    rows = []
    for dirp in sorted(BOLTZ_POOL.glob("*__unsteered__*")):
        lig = dirp.name.split("__")[0]
        if lig not in prot:
            continue
        for f in sorted(dirp.glob("input_model_*.cif")):
            try:
                v = in_frame(P.load_structure(f))
            except Exception:
                continue
            if v is None:
                continue
            ds = [chamfer(v, q) for q in prot[lig]]
            rows.append({"ligand": lig, "sample": f.stem,
                         "xeng_min": float(np.min(ds)),
                         "xeng_mean": float(np.mean(ds))})
    f = pd.DataFrame(rows)
    if f.empty:
        return {"error": "no overlapping ligands between the two pools"}

    sc = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    sc = sc[sc.arm == "unsteered"]
    m = sc.merge(f, on=["ligand", "sample"], how="inner")
    if m.empty:
        return {"error": "join produced nothing"}

    rng = np.random.default_rng(0)
    rand = float(np.mean([
        m.groupby("ligand").lddt_pli.apply(lambda s: s.sample(1, random_state=i).iloc[0])
        .mean() for i in range(200)]))
    oracle = float(m.groupby("ligand").lddt_pli.max().mean())

    null = []
    for _ in range(n_null):
        mm = m.assign(s=rng.normal(size=len(m)))
        null.append(float(mm.loc[mm.groupby("ligand").s.idxmax()].lddt_pli.mean()) - rand)
    null = np.array(null)

    out = {"poses": len(m), "ligands": int(m.ligand.nunique()),
           "random": round(rand, 4), "oracle": round(oracle, 4),
           "null_p95": round(float(np.percentile(null, 95)), 4),
           "null_p99": round(float(np.percentile(null, 99)), 4)}

    for col in ("xeng_min", "xeng_mean"):
        # agreement is a DISTANCE, so the good pose is the one with the smallest value
        sel = float(m.loc[m.groupby("ligand")[col].idxmin()].lddt_pli.mean())
        rs = [stats.spearmanr(g[col], g.lddt_pli).statistic
              for _l, g in m.groupby("ligand") if g[col].nunique() >= 3]
        rs = np.array([r for r in rs if np.isfinite(r)])
        gain = sel - rand
        out[col] = {
            "selected": round(sel, 4), "gain_vs_random": round(gain, 4),
            "within_ligand_rho_mean": round(float(rs.mean()), 4),
            "frac_rho_negative": round(float((rs < 0).mean()), 3),
            "wilcoxon_p": round(float(stats.wilcoxon(rs).pvalue), 5),
            "empirical_p": round(float((null >= gain).mean()), 4),
            "verdict": ("REAL - clears the 99th pct"
                        if gain > np.percentile(null, 99) else "inside the noise floor"),
        }
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.parse_args()
    print(json.dumps(main(), indent=2))
