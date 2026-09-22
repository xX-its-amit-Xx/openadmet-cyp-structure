"""Test the one post-hoc direction this study turned up, on data it was not found in.

The a-priori-sign sweep put `fg_depth` - the ligand centroid's height above the porphyrin
plane, measured in the heme frame - at **-0.0376 with the physics sign**, i.e. the sign is
backwards and a pose whose ligand sits LOWER, closer to the heme, is the better pose.
Reversing a sign after seeing the answer is not a result, it is a hypothesis; the honest
move is to state the direction and test it where it was not found.

Two things are computed here, both from predictions only:

* on the **CYP3A4 run5 pool**, the same quantity plus two simpler cousins
  (`fe_centroid_dist`, `fe_min_dist`), to see whether "lower" is anything more than
  "closer to the iron", which FINDINGs 008/010 already declared saturated;
* on the **arm4 P450 holdout** - 85 pairs of OTHER P450 proteins with other ligands,
  five samples each, scored during the fine-tuning campaign - as an independent test of
  the stated direction.

FINDING 010 tested a population KDE over `(r, z)` in this frame and found it null. That is
not this: a density learns a peak, and the claim here is MONOTONE - the co-folder
systematically floats ligands too high, so less height is better at any height.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path("/scratch/shenoy.am/cyp-finetune")
EXP = Path("/scratch/shenoy.am/zexp")
sys.path.insert(0, str(ROOT / "src"))

import numpy as np                                            # noqa: E402
import pandas as pd                                           # noqa: E402
from scipy.stats import spearmanr, wilcoxon                    # noqa: E402

RNG = np.random.default_rng(23)


def depth_terms(cif) -> dict | None:
    from cypstruct.pose import load_structure
    from cypstruct.qmscore.geometry import porphyrin_frame
    cx = load_structure(cif)
    if len(cx.lig_xyz) == 0 or cx.fe is None:
        return None
    fe, normal = porphyrin_frame(cx.heme_xyz, cx.heme_atom, cx.axial_sg, cx.heme_elem)
    if fe is None:
        return None
    cen = cx.lig_xyz.mean(0)
    d = np.linalg.norm(cx.lig_xyz - fe, axis=1)
    out = {"fe_centroid_dist": float(np.linalg.norm(cen - fe)),
           "fe_min_dist": float(d.min())}
    out["fg_depth"] = float(np.dot(cen - fe, normal)) if normal is not None else np.nan
    return out


def evaluate(df: pd.DataFrame, col: str, sign: int, n_draws: int = 64) -> dict:
    groups = list(df.groupby("name"))
    rnd_lig = np.array([g["lddt_pli"].mean() for _, g in groups])
    per = np.zeros(len(groups))
    for _ in range(n_draws):
        for i, (_nm, g) in enumerate(groups):
            v = sign * g[col].astype(float).values
            v = np.where(np.isfinite(v), v, -np.inf)
            per[i] += g["lddt_pli"].values[RNG.choice(np.flatnonzero(v == v.max()))]
    per /= n_draws
    rs = []
    for _nm, g in groups:
        v = sign * g[col].astype(float).values
        ok = np.isfinite(v)
        if ok.sum() < 3 or np.nanstd(v[ok]) == 0:
            continue
        r = spearmanr(v[ok], g["lddt_pli"].values[ok]).statistic
        if r == r:
            rs.append(r)
    rs = np.array(rs)
    try:
        p = float(wilcoxon(per, rnd_lig).pvalue)
    except ValueError:
        p = float("nan")
    # permutation null on this pool: pick a pose at random
    draws = np.array([np.mean([RNG.choice(g["lddt_pli"].values) for _n, g in groups])
                      for _ in range(2000)])
    return {"ligands": len(groups),
            "oracle": round(float(np.mean([g["lddt_pli"].max() for _, g in groups])), 4),
            "random": round(float(rnd_lig.mean()), 4),
            "selected": round(float(per.mean()), 4),
            "gain": round(float(per.mean() - rnd_lig.mean()), 4),
            "rho": round(float(rs.mean()), 3),
            "frac_rho_pos": round(float((rs > 0).mean()), 3),
            "beats_random_frac": round(float((per > rnd_lig).mean()), 3),
            "p_wilcoxon": p,
            "null_p95": round(float(np.percentile(draws, 95) - rnd_lig.mean()), 4),
            "p_permutation": round(float((draws >= per.mean()).mean()), 4)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["cyp3a4", "holdout"], required=True)
    a = ap.parse_args()
    rows = []

    if a.set == "cyp3a4":
        for tag in ("run5_s42", "run5_s43", "run5_s44", "run5_s45"):
            preds = sorted((EXP / "out" / tag).glob("boltz_results_*/predictions"))[0]
            for d in sorted(preds.iterdir()):
                if not d.is_dir():
                    continue
                for cif in sorted(d.glob("*_model_*.cif")):
                    t = depth_terms(cif)
                    if t is None:
                        continue
                    t.update(name=d.name, tag=tag,
                             rank=int(cif.stem.rsplit("_", 1)[1]))
                    rows.append(t)
        df = pd.DataFrame(rows)
        truth = pd.concat([pd.read_csv(f) for f in
                           sorted(EXP.glob("truth_run5_s4*.csv"))], ignore_index=True)
        df = truth.merge(df, on=["name", "tag", "rank"])
    else:
        sc = json.loads((ROOT / "holdout_out" / "scores_arm4_mix_base.json").read_text())
        preds = sorted((ROOT / "holdout_out" / "arm4_mix_base").glob(
            "boltz_results_*/predictions"))[0]
        for d in sorted(preds.iterdir()):
            if not d.is_dir() or d.name not in sc:
                continue
            lp = sc[d.name]["lddt_pli"]
            for cif in sorted(d.glob("*_model_*.cif")):
                k = int(cif.stem.rsplit("_", 1)[1])
                if k >= len(lp):
                    continue
                t = depth_terms(cif)
                if t is None:
                    continue
                t.update(name=d.name, rank=k, lddt_pli=float(lp[k]))
                rows.append(t)
        df = pd.DataFrame(rows)

    out = {"set": a.set, "poses": int(len(df))}
    for col in ("fg_depth", "fe_centroid_dist", "fe_min_dist"):
        if col in df:
            out[f"{col} (lower is better)"] = evaluate(df, col, -1)
    dst = EXP / f"depth_replication_{a.set}.json"
    dst.write_text(json.dumps(out, indent=2))
    df.to_csv(EXP / f"depth_terms_{a.set}.csv", index=False)
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
