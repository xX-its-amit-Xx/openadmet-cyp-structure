"""The three controls that decide whether the physics ensemble's +0.025 is real and useful.

The ensemble beats its own matched-dimensionality null while **every one of its four term
groups is individually a null**. That pattern has two readings - a weak multivariate signal,
or a fit finding structure the null did not happen to find - and only controls separate them.

1. **Within-ligand shuffle.** Permute the feature rows inside each ligand, leaving the label
   column alone, and refit. This destroys the pose-to-feature pairing while preserving every
   marginal distribution, every between-ligand effect and the exact dimensionality. If the
   gain survives, it was never about the pose.
2. **Complementarity.** `CLAUDE.md`: a term ships only if its gain correlates with where the
   incumbent errs. A score that improves the ligands the incumbent already gets right is
   absorbed and worth nothing at inference.
3. **Seed stability.** Re-fit under different random states. The deploy rule requires the win
   to survive resampling.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

EXP = Path("/scratch/shenoy.am/zexp")
sys.path.insert(0, "/scratch/shenoy.am/cyp-finetune/src")
sys.path.insert(0, str(EXP))

import numpy as np                                            # noqa: E402
import pandas as pd                                           # noqa: E402
from scipy.stats import pearsonr, spearmanr, wilcoxon          # noqa: E402
from sklearn.ensemble import HistGradientBoostingRegressor     # noqa: E402

from cypstruct.qmscore.pocket import FEATURES                  # noqa: E402

RNG = np.random.default_rng(7)


def loo(X, y, names, seed=0):
    p = np.empty(len(y))
    for n in np.unique(names):
        m = names == n
        g = HistGradientBoostingRegressor(max_iter=200, random_state=seed).fit(X[~m], y[~m])
        p[m] = g.predict(X[m])
    return p


def per_ligand_select(groups, index, score, n=64):
    out = np.zeros(len(groups))
    for k in range(n):
        for i, (nm, g) in enumerate(groups):
            v = score[index[nm]]
            v = np.where(np.isfinite(v), v, -np.inf)
            out[i] += g["lddt_pli"].values[RNG.choice(np.flatnonzero(v == v.max()))]
    return out / n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--out", default=str(EXP / "pocket_controls.json"))
    a = ap.parse_args()

    truth = pd.concat([pd.read_csv(f) for f in
                       sorted(glob.glob(str(EXP / "truth_run5_s4*.csv")))],
                      ignore_index=True)
    pk = pd.concat([pd.read_csv(f) for f in
                    sorted(glob.glob(str(EXP / "pocket_terms_run5_s4*.csv")))],
                   ignore_index=True)
    d = truth.merge(pk, on=["name", "tag", "rank"]).merge(
        pd.read_csv(EXP / "xeng_4seed.csv"), on=["name", "tag", "rank"], how="left")
    d = d.sort_values(["name", "tag", "rank"]).reset_index(drop=True)

    groups = list(d.groupby("name"))
    index = {nm: g.index.to_numpy() for nm, g in groups}
    y = d["lddt_pli"].values
    names = d["name"].values
    X = d[FEATURES].astype(float).values

    oracle = np.array([g["lddt_pli"].max() for _, g in groups])
    rnd = np.array([g["lddt_pli"].mean() for _, g in groups])

    phys = per_ligand_select(groups, index, loo(X, y, names, seed=0))
    inc = per_ligand_select(groups, index, -d["xeng"].fillna(np.inf).values)

    rep: dict = {
        "oracle": round(float(oracle.mean()), 4),
        "random": round(float(rnd.mean()), 4),
        "physics": round(float(phys.mean()), 4),
        "incumbent": round(float(inc.mean()), 4),
    }

    # --- 1. within-ligand shuffle ----------------------------------------
    sh = []
    for r in range(a.reps):
        rg = np.random.default_rng(500 + r)
        Xs = X.copy()
        for nm in index:
            idx = index[nm]
            Xs[idx] = X[rg.permutation(idx)]
        sh.append(float(per_ligand_select(groups, index, loo(Xs, y, names)).mean()))
    rep["shuffle_within_ligand"] = {
        "reps": a.reps, "mean": round(float(np.mean(sh)), 4),
        "gains": [round(v - rnd.mean(), 4) for v in sh],
        "mean_gain": round(float(np.mean(sh) - rnd.mean()), 4),
    }

    # --- 2. complementarity with the incumbent ---------------------------
    phys_gain = phys - rnd                 # what physics buys, per ligand
    inc_short = oracle - inc               # where the incumbent still errs
    r_p = pearsonr(phys_gain, inc_short)
    r_s = spearmanr(phys_gain, inc_short)
    rep["complementarity"] = {
        "pearson_r": round(float(r_p.statistic), 3), "pearson_p": float(r_p.pvalue),
        "spearman_r": round(float(r_s.statistic), 3), "spearman_p": float(r_s.pvalue),
        "physics_beats_incumbent_frac": round(float((phys > inc).mean()), 3),
        "physics_vs_incumbent_wilcoxon_p": float(wilcoxon(phys, inc).pvalue),
        "both_selectors_agree_rho": round(float(np.corrcoef(phys, inc)[0, 1]), 3),
    }

    # tail rescue: hand the K ligands with the widest pool spread to physics
    spread = np.array([g["lddt_pli"].max() - g["lddt_pli"].min() for _, g in groups])
    for k in (10, 20, 30):
        pick = np.argsort(-spread)[:k]
        mixed = inc.copy()
        mixed[pick] = phys[pick]
        rep.setdefault("tail_rescue", {})[f"k={k}"] = round(float(mixed.mean()), 4)

    # --- 3. seed stability ------------------------------------------------
    seeds = [round(float(per_ligand_select(groups, index, loo(X, y, names, seed=s)).mean()), 4)
             for s in (0, 1, 2, 3, 4)]
    rep["seed_stability"] = {"selected": seeds,
                             "gains": [round(v - rnd.mean(), 4) for v in seeds],
                             "sd": round(float(np.std(seeds)), 4)}

    Path(a.out).write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
