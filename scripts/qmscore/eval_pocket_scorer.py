"""Honest evaluation of the CYP3A4 physics scorer on the run5 pool.

The protocol is the one this repo has settled on after four negatives, and every guard
here exists because its absence has already produced a wrong answer once:

* **leave-one-LIGAND-out**, never random rows. A random split puts the same ligand on
  both sides of the fold and inflates everything.
* **ties broken at random over 64 draws.** A constant or near-constant score otherwise
  always lands on `_model_0`, which is Boltz's own confidence order reported under a
  different name. That trap has fired here before.
* **the null is re-derived for THIS setup** with a matched-dimensionality random feature
  pushed through the identical fit-and-select pipeline. A borrowed noise floor is not a
  noise floor.
* **per-term numbers, not just the ensemble**, so a null term is visible as a null.
* **within-ligand Spearman**, because between-ligand structure is not selection.

    ./env/bin/python eval_pocket_scorer.py --pocket .../pocket_terms.csv
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

EXP = Path("/scratch/shenoy.am/zexp")
sys.path.insert(0, "/scratch/shenoy.am/cyp-finetune/src")

import numpy as np                                            # noqa: E402
import pandas as pd                                           # noqa: E402
from scipy.stats import spearmanr, wilcoxon                    # noqa: E402
from sklearn.ensemble import HistGradientBoostingRegressor     # noqa: E402

from cypstruct.qmscore.pocket import FEATURES, GROUPS          # noqa: E402

RNG = np.random.default_rng(20260922)
N_DRAWS = 64


# --------------------------------------------------------------------------


def loo_predict(X: np.ndarray, y: np.ndarray, names: np.ndarray) -> np.ndarray:
    """Leave-one-ligand-out predictions. One model per held-out ligand, 87 of them."""
    pred = np.empty(len(y))
    for n in np.unique(names):
        m = names == n
        g = HistGradientBoostingRegressor(max_iter=200, random_state=0)
        g.fit(X[~m], y[~m])
        pred[m] = g.predict(X[m])
    return pred


def loo_sign(v: np.ndarray, y: np.ndarray, names: np.ndarray) -> np.ndarray:
    """Orient ONE raw feature so that higher = better, deciding the sign out-of-fold.

    Single features need no fitting beyond a sign, and a sign chosen on the held-out
    ligand would be a leak all by itself - with 87 ligands and a binary choice it is
    worth roughly the whole effect being measured.
    """
    out = np.empty(len(v))
    uniq = np.unique(names)
    per: dict[str, float] = {}
    for t in uniq:
        mm = names == t
        a, b = v[mm], y[mm]
        ok = np.isfinite(a)
        if ok.sum() < 3 or np.nanstd(a[ok]) == 0:
            continue
        r = spearmanr(a[ok], b[ok]).statistic
        if r == r:
            per[t] = float(r)
    tot, cnt = sum(per.values()), len(per)
    for n in uniq:
        m = names == n
        held = per.get(n)
        s_tot = tot - (held or 0.0)
        s_cnt = cnt - (1 if held is not None else 0)
        s = 1.0 if (s_cnt > 0 and s_tot / s_cnt >= 0) else -1.0
        out[m] = s * v[m]
    return out


def select(groups: list[tuple[str, pd.DataFrame]], score: np.ndarray,
           index: dict, n_draws: int = N_DRAWS) -> tuple[float, np.ndarray]:
    """Mean LDDT-PLI of the argmax pose, ties broken at random, averaged over draws.

    Returns (mean over draws, per-ligand mean over draws) so that a paired test against
    random selection is possible.
    """
    per_lig = np.zeros((n_draws, len(groups)))
    for k in range(n_draws):
        for i, (nm, g) in enumerate(groups):
            idx = index[nm]
            v = score[idx]
            v = np.where(np.isfinite(v), v, -np.inf)
            cand = np.flatnonzero(v == v.max())
            per_lig[k, i] = g["lddt_pli"].values[RNG.choice(cand)]
    return float(per_lig.mean()), per_lig.mean(0)


def within_rho(groups, score, index, y_col="lddt_pli"):
    rs = []
    for nm, g in groups:
        v = score[index[nm]]
        ok = np.isfinite(v)
        if ok.sum() < 3 or np.nanstd(v[ok]) == 0:
            continue
        r = spearmanr(v[ok], g[y_col].values[ok]).statistic
        if r == r:
            rs.append(r)
    rs = np.array(rs)
    return float(rs.mean()), float((rs > 0).mean()), len(rs)


# --------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pocket", default=str(EXP / "pocket_terms_run5_s4*.csv"))
    ap.add_argument("--null-reps", type=int, default=60)
    ap.add_argument("--out", default=str(EXP / "pocket_eval.json"))
    a = ap.parse_args()

    truth = pd.concat([pd.read_csv(f) for f in
                       sorted(glob.glob(str(EXP / "truth_run5_s4*.csv")))],
                      ignore_index=True)
    pk = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(a.pocket))],
                   ignore_index=True)
    d = truth.merge(pk, on=["name", "tag", "rank"], how="inner")
    d = d.sort_values(["name", "tag", "rank"]).reset_index(drop=True)

    # the incumbent, for a same-pool comparison rather than a quoted one
    xe = EXP / "xeng_4seed.csv"
    has_x = xe.exists()
    if has_x:
        d = d.merge(pd.read_csv(xe), on=["name", "tag", "rank"], how="left")

    groups = list(d.groupby("name"))
    index = {nm: g.index.to_numpy() for nm, g in groups}
    y = d["lddt_pli"].values
    names = d["name"].values

    oracle = float(np.mean([g["lddt_pli"].max() for _, g in groups]))
    random_ = float(np.mean([g["lddt_pli"].mean() for _, g in groups]))
    per_lig_random = np.array([g["lddt_pli"].mean() for _, g in groups])

    report: dict = {
        "poses": int(len(d)), "ligands": int(d["name"].nunique()),
        "oracle": round(oracle, 4), "random": round(random_, 4),
        "numbering": {
            "median_offset": float(d["resnum_offset_y"].median())
            if "resnum_offset_y" in d else float(d["resnum_offset"].median()),
            "min_seq_identity": round(float(d["seq_identity"].min()), 3),
        },
    }

    # --- the random-SELECTION null (which pose you land on by luck) -------
    draws = np.array([np.mean([RNG.choice(g["lddt_pli"].values) for _, g in groups])
                      for _ in range(4000)])
    report["null_random_selection"] = {
        "p95": round(float(np.percentile(draws, 95) - random_), 4),
        "p99": round(float(np.percentile(draws, 99) - random_), 4),
    }

    # --- the matched-dimensionality FEATURE null -------------------------
    # A random feature block of exactly the same width, through the identical
    # LOO-fit-and-select pipeline. This is the number a fitted scorer has to clear;
    # it is strictly larger than the selection null above because the fit itself
    # can find structure in noise.
    dim = len(FEATURES)
    null_gains = []
    for rep in range(a.null_reps):
        rg = np.random.default_rng(1000 + rep)
        Xn = rg.standard_normal((len(d), dim))
        p = loo_predict(Xn, y, names)
        s, _ = select(groups, p, index, n_draws=16)
        null_gains.append(s - random_)
    null_gains = np.array(null_gains)
    null95 = float(np.percentile(null_gains, 95))
    report["null_matched_feature"] = {
        "dim": dim, "reps": int(a.null_reps),
        "mean": round(float(null_gains.mean()), 4),
        "sd": round(float(null_gains.std()), 4),
        "p95": round(null95, 4), "max": round(float(null_gains.max()), 4),
    }

    # --- incumbent on the same pool --------------------------------------
    rows = []

    def add(label: str, score: np.ndarray, kind: str):
        sel, per_lig = select(groups, score, index)
        rho, frac_pos, n = within_rho(groups, score, index)
        beats = float((per_lig > per_lig_random).mean())
        try:
            p = float(wilcoxon(per_lig, per_lig_random).pvalue)
        except ValueError:
            p = float("nan")
        rows.append({"term": label, "kind": kind, "selected": round(sel, 4),
                     "gain": round(sel - random_, 4), "rho": round(rho, 3),
                     "frac_rho_pos": round(frac_pos, 3), "n_rho": n,
                     "beats_random_frac": round(beats, 3), "p_wilcoxon": p})

    if has_x and d["xeng"].notna().any():
        add("INCUMBENT xeng", -d["xeng"].fillna(np.inf).values, "incumbent")

    # --- the control that decides whether any of this is CYP3A4-specific --
    # FINDING 003 already showed that a plain pocket-contact count selects (+0.0220).
    # If a residue-specific term does not beat a generic burial count, it is that
    # result wearing a Phe cluster for a hat. `fg_frac` is fg_contacts / all contacts,
    # so the generic count is recoverable exactly from what was written out.
    with np.errstate(divide="ignore", invalid="ignore"):
        n_all = d["fg_contacts"].values / d["fg_frac"].values
    d["n_contacts_all"] = n_all
    add("CONTROL generic contact count", loo_sign(
        np.nan_to_num(n_all, nan=float(np.nanmedian(n_all))), y, names), "control")
    add("CONTROL Boltz sample rank", -d["rank"].astype(float).values, "control")

    # --- every single feature on its own ---------------------------------
    for f in FEATURES:
        v = d[f].astype(float).values
        if np.nanstd(v) == 0 or np.isfinite(v).sum() < len(v) * 0.5:
            rows.append({"term": f, "kind": "single", "selected": None,
                         "gain": None, "rho": None, "frac_rho_pos": None,
                         "n_rho": 0, "beats_random_frac": None,
                         "p_wilcoxon": None, "note": "constant or mostly missing"})
            continue
        add(f, loo_sign(np.nan_to_num(v, nan=np.nanmedian(v)), y, names), "single")

    # --- each group on its own, and the ensemble minus each group --------
    def block(cols):
        # NaN is left in place on purpose: HistGradientBoosting handles missing
        # natively, and "this pose has no aromatic ring near the Phe roof" is
        # information, not a value to be imputed away.
        return d[cols].astype(float).values

    for gname, cols in GROUPS.items():
        add(f"GROUP {gname}", loo_predict(block(cols), y, names), "group")

    add("ENSEMBLE all", loo_predict(block(FEATURES), y, names), "ensemble")
    for gname, cols in GROUPS.items():
        rest = [c for c in FEATURES if c not in cols]
        add(f"ABLATE -{gname}", loo_predict(block(rest), y, names), "ablation")

    # physics + incumbent, on within-ligand ranks so no weight is fitted
    if has_x and d["xeng"].notna().any():
        ens = loo_predict(block(FEATURES), y, names)
        def rk(v, idx):
            s = pd.Series(v[idx]).rank()
            return (s / max(len(s), 1)).values
        comb = np.empty(len(d))
        for nm, g in groups:
            idx = index[nm]
            comb[idx] = rk(ens, idx) + rk(-d["xeng"].fillna(np.inf).values, idx)
        add("ENSEMBLE + incumbent (rank-avg)", comb, "combined")

    report["results"] = rows
    Path(a.out).write_text(json.dumps(report, indent=2))

    hdr = f"{'term':34s} {'sel':>7s} {'gain':>8s} {'rho':>7s} {'r>0':>6s} {'beat':>6s} {'p':>9s}"
    print(f"\nposes {len(d)}  ligands {d['name'].nunique()}  "
          f"oracle {oracle:.4f}  random {random_:.4f}")
    print(f"matched-dim feature null: mean {null_gains.mean():+.4f}  "
          f"p95 {null95:+.4f}  (dim {dim}, {a.null_reps} reps)")
    print(f"random-selection null p95 {report['null_random_selection']['p95']:+.4f}\n")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        if r["gain"] is None:
            print(f"{r['term']:34s} {'-':>7s} {'-':>8s}   {r.get('note','')}")
            continue
        print(f"{r['term']:34s} {r['selected']:7.4f} {r['gain']:+8.4f} "
              f"{r['rho']:+7.3f} {r['frac_rho_pos']:6.2f} "
              f"{r['beats_random_frac']:6.2f} {r['p_wilcoxon']:9.2g}")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
