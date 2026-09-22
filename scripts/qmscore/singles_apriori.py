"""Single-term selection with the sign fixed BY PHYSICS, not by a fit.

The first pass oriented each feature by its out-of-fold mean within-ligand rho. That rule
**degenerates on exactly the features it is meant to test**: when a term's true mean rho is
zero, the leave-one-out mean is dominated by the removal of the held-out ligand's own rho,
so the chosen sign is anti-correlated with it and the term reports a confident negative.
`n_hb_backbone` came back at rho -0.316 with 0 of 57 ligands positive and p = 9e-09 while
its RAW mean rho is +0.0005 on 51% of ligands. That is the artifact, not a result.

A physics term does not need its sign fitted. More satisfied H-bonds is better, shorter is
better, more stacking is better, less strain and less clash are better. Fixing the sign a
priori removes the artifact and makes the single-term test parameter-free.
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

EXP = Path("/scratch/shenoy.am/zexp")
sys.path.insert(0, "/scratch/shenoy.am/cyp-finetune/src")

import numpy as np                                            # noqa: E402
import pandas as pd                                           # noqa: E402
from scipy.stats import spearmanr, wilcoxon                    # noqa: E402

from cypstruct.qmscore.pocket import GROUPS                    # noqa: E402

RNG = np.random.default_rng(11)

# +1 = larger is a better pose, -1 = smaller is a better pose. Chemistry, not fitting.
SIGN = {
    "hb_ser119": +1, "hb_arg106": +1, "hb_arg212": +1, "hb_asp214": +1,
    "hb_thr224": +1, "anchor_score": +1, "n_anchor_hb": +1,
    "hb_backbone": +1, "n_hb_backbone": +1, "unsat_buried_polar": -1,
    "d_ser119": -1, "d_arg106": -1, "d_arg212": -1, "d_asp214": -1, "d_thr224": -1,
    "phe_n_contacts": +1, "phe_n_pd": +1, "phe_n_t": +1, "phe_stack_score": +1,
    "phe_best_score": +1, "phe_min_dist": -1, "phe_n_res": +1,
    "phe_heavy_contacts": +1, "f304_dist": -1,
    "fg_contacts": +1, "fg_sc_contacts": +1, "fg_frac": +1, "fg_min_dist": -1,
    "fg_n_res": +1, "fg_depth": +1,
    "strain": -1, "max_clash": -1, "n_clash": -1, "self_clash": -1,
}


def main() -> int:
    truth = pd.concat([pd.read_csv(f) for f in
                       sorted(glob.glob(str(EXP / "truth_run5_s4*.csv")))],
                      ignore_index=True)
    pk = pd.concat([pd.read_csv(f) for f in
                    sorted(glob.glob(str(EXP / "pocket_terms_run5_s4*.csv")))],
                   ignore_index=True)
    d = truth.merge(pk, on=["name", "tag", "rank"])
    d = d.sort_values(["name", "tag", "rank"]).reset_index(drop=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        d["n_contacts_all"] = d["fg_contacts"] / d["fg_frac"]

    groups = list(d.groupby("name"))
    index = {nm: g.index.to_numpy() for nm, g in groups}
    rnd_lig = np.array([g["lddt_pli"].mean() for _, g in groups])
    rnd = float(rnd_lig.mean())
    oracle = float(np.mean([g["lddt_pli"].max() for _, g in groups]))

    def evaluate(v: np.ndarray, n_draws: int = 64):
        per = np.zeros(len(groups))
        for _ in range(n_draws):
            for i, (nm, g) in enumerate(groups):
                s = v[index[nm]]
                s = np.where(np.isfinite(s), s, -np.inf)
                per[i] += g["lddt_pli"].values[RNG.choice(np.flatnonzero(s == s.max()))]
        per /= n_draws
        rs = []
        for nm, g in groups:
            s = v[index[nm]]
            ok = np.isfinite(s)
            if ok.sum() < 3 or np.nanstd(s[ok]) == 0:
                continue
            r = spearmanr(s[ok], g["lddt_pli"].values[ok]).statistic
            if r == r:
                rs.append(r)
        rs = np.array(rs)
        try:
            p = float(wilcoxon(per, rnd_lig).pvalue)
        except ValueError:
            p = float("nan")
        return {"selected": round(float(per.mean()), 4),
                "gain": round(float(per.mean() - rnd), 4),
                "rho": round(float(rs.mean()) if len(rs) else float("nan"), 3),
                "frac_rho_pos": round(float((rs > 0).mean()) if len(rs) else float("nan"), 3),
                "n_rho": int(len(rs)),
                "beats_random_frac": round(float((per > rnd_lig).mean()), 3),
                "p_wilcoxon": p}

    rows = []
    for f, s in SIGN.items():
        r = evaluate(s * d[f].astype(float).values)
        r["term"] = f
        r["sign"] = s
        r["group"] = [k for k, c in GROUPS.items() if f in c][0]
        rows.append(r)
    r = evaluate(d["n_contacts_all"].astype(float).values)
    r.update(term="CONTROL generic contact count", sign=+1, group="control")
    rows.append(r)

    out = {"oracle": round(oracle, 4), "random": round(rnd, 4), "results": rows}
    (EXP / "pocket_singles_apriori.json").write_text(json.dumps(out, indent=2))
    print(f"oracle {oracle:.4f}  random {rnd:.4f}\n")
    print(f"{'term':24s} {'sgn':>3s} {'sel':>7s} {'gain':>8s} {'rho':>7s} "
          f"{'r>0':>5s} {'n':>4s} {'beat':>5s} {'p':>8s}")
    for r in sorted(rows, key=lambda x: -x["gain"]):
        print(f"{r['term']:24s} {r['sign']:+3d} {r['selected']:7.4f} {r['gain']:+8.4f} "
              f"{r['rho']:+7.3f} {r['frac_rho_pos']:5.2f} {r['n_rho']:4d} "
              f"{r['beats_random_frac']:5.2f} {r['p_wilcoxon']:8.2g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
