"""Matched-dimensionality null for the physics ensemble, run wide enough to price a tail.

60 reps put the null's p95 at +0.0123 and its MAXIMUM at +0.0257 - above the observed
+0.0252. At that resolution the empirical p-value is 1/60 and the verdict is a coin toss
about a tail, which is exactly the situation FINDING 007 exists to prevent. This runs the
same pipeline on 34 Gaussian features many more times and writes every gain out.
"""
import argparse, glob, sys
from pathlib import Path
EXP = Path("/scratch/shenoy.am/zexp")
sys.path.insert(0, "/scratch/shenoy.am/cyp-finetune/src")
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

ap = argparse.ArgumentParser()
ap.add_argument("--start", type=int, required=True)
ap.add_argument("--reps", type=int, default=25)
ap.add_argument("--dim", type=int, default=34)
a = ap.parse_args()

truth = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(EXP/"truth_run5_s4*.csv")))], ignore_index=True)
pk = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(EXP/"pocket_terms_run5_s4*.csv")))], ignore_index=True)
d = truth.merge(pk, on=["name","tag","rank"]).sort_values(["name","tag","rank"]).reset_index(drop=True)
groups = list(d.groupby("name")); index = {nm: g.index.to_numpy() for nm,g in groups}
y = d["lddt_pli"].values; names = d["name"].values
rnd = float(np.mean([g["lddt_pli"].mean() for _,g in groups]))

out = []
for r in range(a.start, a.start + a.reps):
    rg = np.random.default_rng(90000 + r)
    X = rg.standard_normal((len(d), a.dim))
    p = np.empty(len(y))
    for n in np.unique(names):
        m = names == n
        p[m] = HistGradientBoostingRegressor(max_iter=200, random_state=0).fit(X[~m], y[~m]).predict(X[m])
    tot = 0.0
    for k in range(32):
        for nm, g in groups:
            v = p[index[nm]]
            tot += g["lddt_pli"].values[rg.choice(np.flatnonzero(v == v.max()))]
    out.append(tot / (32 * len(groups)) - rnd)
    print(r, round(out[-1], 5), flush=True)
np.save(EXP / f"nullgains_{a.start}.npy", np.array(out))
