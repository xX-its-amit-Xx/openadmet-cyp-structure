"""The depth curve on the predicted-Type-I stratum — oracle first, then selection.

Pre-registered in `docs/PREREG_type_i_depth.md`. Rungs, draw counts, strata and the
acceptance rule are fixed there; nothing here may move them.

Four stages:

    python scripts/structure/type_i_depth_analysis.py controls    # N2, C7, pool diversity
    python scripts/structure/type_i_depth_analysis.py subsample    # 4.1, free, exact
    python scripts/structure/type_i_depth_analysis.py augment      # 4.2, with new poses
    python scripts/structure/type_i_depth_analysis.py project      # 4.3

**Why `selected` is computed as argmin(xeng) within the drawn subset.** The shipped
selector is `-z(xeng)` z-scored WITHIN the ligand, and z-scoring is a monotone transform
of xeng inside a ligand, so the argmax of `-z(xeng)` is the argmin of `xeng`. That
identity is asserted against `cypstruct.xengine.select()` at full depth in `controls`
rather than assumed, because it is the one shortcut that makes 256 draws x 8 rungs x
4,000 null features affordable.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

# --- fixed by the pre-registration -----------------------------------------
COORD_MAX = 2.6
RUNGS = [1, 2, 3, 5, 8, 10, 14, 20]
N_DRAWS = 256          # prereg asks for >= 64
N_NULL = 4000
N_BOOT = 10000
SEED = 20260922
BENCH_SEL = 0.0125     # FINDING 004, selection, per doubling
BENCH_ORA = 0.024      # FINDING 004, oracle, per doubling

OUT_CTRL = DATA_PROCESSED / "type_i_depth_controls.json"
OUT_SUB = DATA_PROCESSED / "type_i_depth_subsample.json"
OUT_SUB_CSV = DATA_PROCESSED / "type_i_depth_subsample_per_ligand.csv"
OUT_AUG = DATA_PROCESSED / "type_i_depth_augment.json"
OUT_PROJ = DATA_PROCESSED / "type_i_depth_projection.json"


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

def load_pool() -> pd.DataFrame:
    P = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    P = P[P.arm == "unsteered"].copy()
    X = pd.read_csv(DATA_PROCESSED / "xeng_val87b.csv")
    d = P.merge(X, on=["ligand", "sample"], how="left", validate="one_to_one")
    assert d.xeng.notna().all(), "xeng did not join onto every pose"
    return d


def strata() -> dict[str, list[str]]:
    L = pd.read_csv(DATA_PROCESSED / "binding_mode_labels_cyp3a4.csv")
    ti = sorted(L.loc[L.pred_fe_donor_median > COORD_MAX, "id"])
    tii = sorted(L.loc[L.pred_fe_donor_median <= COORD_MAX, "id"])
    return {"type_I_pred": ti, "type_II_pred": tii, "all": sorted(L.id)}


def as_arrays(d: pd.DataFrame, ligands: list[str]):
    """-> (ligand order, lddt [n_lig, n_pose], xeng [n_lig, n_pose])."""
    ligs, L, XE = [], [], []
    for lig in ligands:
        g = d[d.ligand == lig]
        if len(g) == 0:
            continue
        ligs.append(lig)
        L.append(g.lddt_pli.to_numpy(float))
        XE.append(g.xeng.to_numpy(float))
    n = min(len(a) for a in L)
    assert all(len(a) == n for a in L), "ragged pool"
    return ligs, np.array(L), np.array(XE)


# ---------------------------------------------------------------------------
# the subsample engine
# ---------------------------------------------------------------------------

def curve(L: np.ndarray, XE: np.ndarray, rungs=RUNGS, n_draws=N_DRAWS, seed=SEED):
    """Per-ligand oracle / selected / random at each rung, averaged over draws.

    Random tie-breaking: ties in `xeng` are broken by a per-draw random key, so a pool
    with duplicate feature values cannot be silently ranked by file order.
    """
    rng = np.random.default_rng(seed)
    n_lig, n_pose = L.shape
    out = {}
    for k in rungs:
        ora = np.zeros((n_lig, n_draws))
        sel = np.zeros((n_lig, n_draws))
        ran = np.zeros((n_lig, n_draws))
        for t in range(n_draws):
            idx = np.argsort(rng.random((n_lig, n_pose)), axis=1)[:, :k]
            l_sub = np.take_along_axis(L, idx, axis=1)
            x_sub = np.take_along_axis(XE, idx, axis=1)
            tie = rng.random(x_sub.shape) * 1e-9
            pick = np.argmin(x_sub + tie, axis=1)
            ora[:, t] = l_sub.max(axis=1)
            sel[:, t] = l_sub[np.arange(n_lig), pick]
            ran[:, t] = l_sub.mean(axis=1)
        out[k] = {"oracle": ora.mean(axis=1), "selected": sel.mean(axis=1),
                  "random": ran.mean(axis=1)}
    return out


def slope_per_doubling(rungs, vals) -> float:
    x = np.log2(np.asarray(rungs, float))
    y = np.asarray(vals, float)
    return float(np.polyfit(x, y, 1)[0])


def boot_slope(per_lig: dict, rungs, key, n_boot=N_BOOT, seed=SEED):
    """Bootstrap the per-doubling slope over LIGANDS."""
    mat = np.array([per_lig[k][key] for k in rungs])        # [n_rung, n_lig]
    rng = np.random.default_rng(seed + 1)
    n_lig = mat.shape[1]
    point = slope_per_doubling(rungs, mat.mean(axis=1))
    draws = np.empty(n_boot)
    x = np.log2(np.asarray(rungs, float))
    xc = x - x.mean()
    denom = float((xc ** 2).sum())
    for b in range(n_boot):
        ii = rng.integers(0, n_lig, n_lig)
        y = mat[:, ii].mean(axis=1)
        draws[b] = float((xc * (y - y.mean())).sum() / denom)
    return point, float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def null_floor(L: np.ndarray, k: int, n_null=N_NULL, seed=SEED) -> dict:
    """FINDING 007's floor, recomputed INSIDE this stratum at this depth.

    A random feature stands in for xeng; the statistic is the selection gain over the
    exact random baseline. The floor is n-dependent, which is the whole point.
    """
    rng = np.random.default_rng(seed + 7)
    n_lig, n_pose = L.shape
    base = L.mean(axis=1).mean()
    gains = np.empty(n_null)
    for b in range(n_null):
        idx = np.argsort(rng.random((n_lig, n_pose)), axis=1)[:, :k]
        l_sub = np.take_along_axis(L, idx, axis=1)
        pick = np.argmin(rng.random(l_sub.shape), axis=1)   # a random feature
        gains[b] = l_sub[np.arange(n_lig), pick].mean() - base
    return {"p50": float(np.percentile(gains, 50)),
            "p95": float(np.percentile(gains, 95)),
            "p99": float(np.percentile(gains, 99))}


# ---------------------------------------------------------------------------
# stage 1 — controls
# ---------------------------------------------------------------------------

def controls() -> dict:
    from cypstruct import xengine as X

    d = load_pool()
    S = strata()
    res = {"n_rows": int(len(d)), "n_ligands": int(d.ligand.nunique()),
           "strata_sizes": {k: len(v) for k, v in S.items()}}

    # N2 - the shipped column reproduces
    sel = X.select(d)
    exact_random = d.groupby("ligand").lddt_pli.mean().mean()
    rhos, signs = [], []
    for lig, g in d.groupby("ligand"):
        if g.lddt_pli.std() > 0:
            r = stats.spearmanr(g.xeng, g.lddt_pli).statistic
            rhos.append(r)
            signs.append(r < 0)
    res["N2_shipped_column"] = {
        "selected": round(float(sel.lddt_pli.mean()), 4),
        "oracle": round(float(d.groupby("ligand").lddt_pli.max().mean()), 4),
        "random_exact": round(float(exact_random), 4),
        "gain": round(float(sel.lddt_pli.mean() - exact_random), 4),
        "within_ligand_rho": round(float(np.mean(rhos)), 4),
        "frac_correct_sign": round(float(np.mean(signs)), 4),
        "expected": {"selected": 0.6164, "oracle": 0.6975, "random_exact": 0.5769,
                     "gain": 0.0395, "within_ligand_rho": -0.2582,
                     "frac_correct_sign": 0.759},
    }
    res["N2_passes"] = bool(
        abs(sel.lddt_pli.mean() - 0.6164) < 5e-4
        and abs(d.groupby("ligand").lddt_pli.max().mean() - 0.6975) < 5e-4
        and abs(np.mean(rhos) + 0.2582) < 5e-4)

    # the argmin(xeng) shortcut equals select() exactly
    ligs, L, XE = as_arrays(d, S["all"])
    pick = np.argmin(XE, axis=1)
    short = float(L[np.arange(len(ligs)), pick].mean())
    res["argmin_shortcut_equals_select"] = {
        "shortcut": round(short, 6), "select": round(float(sel.lddt_pli.mean()), 6),
        "equal_to_1e-9": bool(abs(short - sel.lddt_pli.mean()) < 1e-9)}

    # C7 - the four exact zeros are ejections, not FINDING 021
    z = d[d.lddt_pli == 0]
    res["C7_exact_zeros"] = {
        "n": int(len(z)),
        "by_ligand": z.ligand.value_counts().to_dict(),
        "bisy_rmsd": [round(float(v), 1) for v in z.bisy_rmsd],
        "all_mapped": bool(z.mapped.all()) if len(z) else None,
        "other_poses_median_lddt": {
            lig: round(float(d[(d.ligand == lig) & (d.lddt_pli > 0)].lddt_pli.median()), 3)
            for lig in z.ligand.unique()},
    }

    # pool diversity - distinct poses, not jobs (README trap 1, applied to the POOL)
    g = d.groupby("ligand")
    res["pool_diversity"] = {
        "median_unique_xeng_per_ligand": float(g.xeng.nunique().median()),
        "min_unique_xeng_per_ligand": int(g.xeng.nunique().min()),
        "median_within_ligand_sd_lddt": round(float(g.lddt_pli.std().median()), 4),
        "n_ligands_with_zero_lddt_sd": int((g.lddt_pli.std() == 0).sum())}

    # stratum absolute numbers, reproducing FINDING 033 table 2
    for name in ("type_I_pred", "type_II_pred", "all"):
        ligs, L, XE = as_arrays(d, S[name])
        pick = np.argmin(XE, axis=1)
        res.setdefault("strata_at_depth_20", {})[name] = {
            "n": len(ligs),
            "random": round(float(L.mean(axis=1).mean()), 4),
            "oracle": round(float(L.max(axis=1).mean()), 4),
            "selected": round(float(L[np.arange(len(ligs)), pick].mean()), 4),
            "headroom": round(float(L.max(axis=1).mean() - L.mean(axis=1).mean()), 4),
        }
    OUT_CTRL.write_text(json.dumps(res, indent=2))
    return res


# ---------------------------------------------------------------------------
# stage 2 — the subsample depth curve (free, internally consistent)
# ---------------------------------------------------------------------------

def subsample() -> dict:
    d = load_pool()
    S = strata()
    res = {"rungs": RUNGS, "n_draws": N_DRAWS, "benchmarks":
           {"selection_per_doubling_FINDING_004": BENCH_SEL,
            "oracle_per_doubling_FINDING_004": BENCH_ORA}}
    rows = []
    for name in ("type_I_pred", "type_II_pred", "all"):
        ligs, L, XE = as_arrays(d, S[name])
        per = curve(L, XE)
        tab = {}
        for k in RUNGS:
            tab[k] = {"oracle": round(float(per[k]["oracle"].mean()), 4),
                      "selected": round(float(per[k]["selected"].mean()), 4),
                      "random_empirical": round(float(per[k]["random"].mean()), 4),
                      "random_exact": round(float(L.mean(axis=1).mean()), 4),
                      "gain_vs_random": round(float(per[k]["selected"].mean()
                                                    - L.mean(axis=1).mean()), 4),
                      "null_p95": round(null_floor(L, k)["p95"], 4)}
            for i, lig in enumerate(ligs):
                rows.append({"stratum": name, "ligand": lig, "depth": k,
                             "oracle": per[k]["oracle"][i],
                             "selected": per[k]["selected"][i],
                             "random": per[k]["random"][i]})
        so, so_lo, so_hi = boot_slope(per, RUNGS, "oracle")
        ss, ss_lo, ss_hi = boot_slope(per, RUNGS, "selected")
        sr, _, _ = boot_slope(per, RUNGS, "random")
        res[name] = {
            "n": len(ligs), "table": tab,
            "oracle_per_doubling": round(so, 5),
            "oracle_per_doubling_CI": [round(so_lo, 5), round(so_hi, 5)],
            "selected_per_doubling": round(ss, 5),
            "selected_per_doubling_CI": [round(ss_lo, 5), round(ss_hi, 5)],
            "random_per_doubling": round(sr, 6),
            "random_flat": bool(abs(sr) < 0.002),
        }
    # is the Type I selection rate different from Type II's?
    res["rate_difference_I_minus_II"] = round(
        res["type_I_pred"]["selected_per_doubling"]
        - res["type_II_pred"]["selected_per_doubling"], 5)
    pd.DataFrame(rows).to_csv(OUT_SUB_CSV, index=False)
    OUT_SUB.write_text(json.dumps(res, indent=2))
    return res


# ---------------------------------------------------------------------------
# stage 3 — the augmented curve (new poses)
# ---------------------------------------------------------------------------

def augment() -> dict:
    new_path = DATA_PROCESSED / "type_i_depth_poses.csv"
    if not new_path.exists():
        return {"error": "no new poses yet"}
    from cypstruct import xengine as X

    d = load_pool()
    S = strata()
    new = pd.read_csv(new_path)

    # N4 - distinct poses, not jobs. Dedupe on (ligand, rounded xeng+lddt) is not
    # enough; the frame-space dedupe lives in the runner's `yield` stage, and this
    # stage additionally drops rows identical in BOTH scored quantities.
    n_raw = len(new)
    new = new[new.lddt_pli.notna() & new.xeng.notna()].copy()
    n_scoreable = len(new)
    new["_k"] = (new.ligand + "|" + new.lddt_pli.round(6).astype(str)
                 + "|" + new.xeng.round(6).astype(str))
    new = new.drop_duplicates("_k").drop(columns="_k")
    n_distinct = len(new)

    ligs = S["type_I_pred"]
    depth_new = new.groupby("ligand").size()
    N = int(depth_new.min()) if len(depth_new) else 0
    res = {"n_new_rows_raw": n_raw, "n_new_scoreable": n_scoreable,
           "n_new_distinct": n_distinct,
           "new_per_ligand": depth_new.to_dict(),
           "N_matched_depth": N,
           "new_arm_quality": {
               "mean_lddt": round(float(new.lddt_pli.mean()), 4),
               "oracle": round(float(new.groupby("ligand").lddt_pli.max().mean()), 4),
               "random": round(float(new.groupby("ligand").lddt_pli.mean().mean()), 4),
           }}
    if N < 1:
        res["verdict_note"] = "no usable new depth"
        OUT_AUG.write_text(json.dumps(res, indent=2))
        return res

    old = d[d.ligand.isin(ligs)][["ligand", "lddt_pli", "xeng"]]
    new_k = new[new.ligand.isin(ligs)][["ligand", "lddt_pli", "xeng"]]
    common = sorted(set(old.ligand) & set(new_k.ligand))
    res["n_ligands_both"] = len(common)

    rungs = sorted({0, int(np.ceil(N / 4)), int(np.ceil(N / 2)), N})
    rng = np.random.default_rng(SEED + 3)
    tab = {}
    per_lig_delta = {}
    for k in rungs:
        ora, sel, ran = [], [], []
        per = {}
        for lig in common:
            lo = old[old.ligand == lig]
            ln = new_k[new_k.ligand == lig]
            o_l, o_x = lo.lddt_pli.to_numpy(), lo.xeng.to_numpy()
            n_l, n_x = ln.lddt_pli.to_numpy(), ln.xeng.to_numpy()
            oo, ss, rr = [], [], []
            for _t in range(N_DRAWS):
                j = rng.permutation(len(n_l))[:k]
                l_all = np.concatenate([o_l, n_l[j]])
                x_all = np.concatenate([o_x, n_x[j]])
                tie = rng.random(x_all.shape) * 1e-9
                oo.append(l_all.max())
                ss.append(l_all[np.argmin(x_all + tie)])
                rr.append(l_all.mean())
            ora.append(np.mean(oo)); sel.append(np.mean(ss)); ran.append(np.mean(rr))
            per[lig] = {"oracle": float(np.mean(oo)), "selected": float(np.mean(ss)),
                        "random": float(np.mean(rr))}
        tab[20 + k] = {"oracle": round(float(np.mean(ora)), 4),
                       "selected": round(float(np.mean(sel)), 4),
                       "random": round(float(np.mean(ran)), 4),
                       "n_new_added": k}
        per_lig_delta[20 + k] = per
    res["table"] = tab
    res["rungs"] = [20 + k for k in rungs]

    # primary endpoint: paired delta, deepest rung minus depth 20
    base, deep = per_lig_delta[20], per_lig_delta[20 + rungs[-1]]
    dsel = np.array([deep[l]["selected"] - base[l]["selected"] for l in common])
    dora = np.array([deep[l]["oracle"] - base[l]["oracle"] for l in common])
    dran = np.array([deep[l]["random"] - base[l]["random"] for l in common])
    rngb = np.random.default_rng(SEED + 5)

    def boot(v):
        dr = np.array([v[rngb.integers(0, len(v), len(v))].mean()
                       for _ in range(N_BOOT)])
        return [round(float(np.percentile(dr, 2.5)), 4),
                round(float(np.percentile(dr, 97.5)), 4)]

    res["primary"] = {
        "delta_selected": round(float(dsel.mean()), 4), "CI": boot(dsel),
        "n_pos": int((dsel > 0).sum()), "n_neg": int((dsel < 0).sum()),
        "wilcoxon_p": (float(stats.wilcoxon(dsel).pvalue)
                       if np.any(dsel != 0) else None)}
    res["secondary"] = {
        "delta_oracle": round(float(dora.mean()), 4), "CI": boot(dora),
        "delta_random": round(float(dran.mean()), 4)}
    res["per_ligand_delta_selected"] = {l: round(float(v), 4)
                                        for l, v in zip(common, dsel)}

    # unmatched "everything available" rung: each ligand contributes ALL its distinct
    # new poses. Reported alongside the matched rungs, never in place of them, because
    # ligands then sit at different depths.
    ora, sel, ran, dep = [], [], [], []
    for lig in common:
        lo = old[old.ligand == lig]
        ln = new_k[new_k.ligand == lig]
        l_all = np.concatenate([lo.lddt_pli.to_numpy(), ln.lddt_pli.to_numpy()])
        x_all = np.concatenate([lo.xeng.to_numpy(), ln.xeng.to_numpy()])
        ora.append(l_all.max()); ran.append(l_all.mean())
        sel.append(l_all[np.argmin(x_all)])
        dep.append(len(l_all))
    res["all_available"] = {
        "mean_depth": round(float(np.mean(dep)), 2),
        "oracle": round(float(np.mean(ora)), 4),
        "selected": round(float(np.mean(sel)), 4),
        "random": round(float(np.mean(ran)), 4),
        "delta_selected_vs_20": round(float(np.mean(sel)
                                            - np.mean([base[l]["selected"]
                                                       for l in common])), 4),
        "delta_oracle_vs_20": round(float(np.mean(ora)
                                          - np.mean([base[l]["oracle"]
                                                     for l in common])), 4),
        "n_new_poses_kept_by_selector": int(sum(
            1 for lig in common
            if np.argmin(np.concatenate([old[old.ligand == lig].xeng.to_numpy(),
                                         new_k[new_k.ligand == lig].xeng.to_numpy()]))
            >= len(old[old.ligand == lig]))),
    }

    # the NEW arm on its own, as a homogeneous pool: does depth convert at the same
    # ~29% inside a single configuration? This separates "the new configuration is
    # worse" from "mixing two configurations breaks selection".
    nrungs = [r for r in (1, 2, 3, 5, 8) if r <= N]
    nl, nx = [], []
    for lig in common:
        ln = new_k[new_k.ligand == lig]
        nl.append(ln.lddt_pli.to_numpy()[:N])
        nx.append(ln.xeng.to_numpy()[:N])
    NL, NX = np.array(nl), np.array(nx)
    ncur = curve(NL, NX, rungs=nrungs)
    res["new_arm_curve"] = {
        "rungs": nrungs, "matched_depth": N,
        "table": {k: {"oracle": round(float(ncur[k]["oracle"].mean()), 4),
                      "selected": round(float(ncur[k]["selected"].mean()), 4),
                      "random": round(float(NL.mean(axis=1).mean()), 4)}
                  for k in nrungs},
        "oracle_per_doubling": round(slope_per_doubling(
            nrungs, [ncur[k]["oracle"].mean() for k in nrungs]), 5),
        "selected_per_doubling": round(slope_per_doubling(
            nrungs, [ncur[k]["selected"].mean() for k in nrungs]), 5),
    }
    # and the OLD pool restricted to the same rungs, for a like-for-like comparison
    d2 = load_pool()
    _ligs, OL, OX = as_arrays(d2[d2.ligand.isin(common)], common)
    ocur = curve(OL, OX, rungs=nrungs)
    res["old_arm_curve_same_rungs"] = {
        "table": {k: {"oracle": round(float(ocur[k]["oracle"].mean()), 4),
                      "selected": round(float(ocur[k]["selected"].mean()), 4)}
                  for k in nrungs},
        "oracle_per_doubling": round(slope_per_doubling(
            nrungs, [ocur[k]["oracle"].mean() for k in nrungs]), 5),
        "selected_per_doubling": round(slope_per_doubling(
            nrungs, [ocur[k]["selected"].mean() for k in nrungs]), 5),
    }

    OUT_AUG.write_text(json.dumps(res, indent=2))
    return res


# ---------------------------------------------------------------------------
# stage 4 — projection onto a Type I-rich test set
# ---------------------------------------------------------------------------

def project() -> dict:
    d = load_pool()
    S = strata()
    sub = json.loads(OUT_SUB.read_text()) if OUT_SUB.exists() else None
    aug = json.loads(OUT_AUG.read_text()) if OUT_AUG.exists() else None

    ligsI, LI, XI = as_arrays(d, S["type_I_pred"])
    ligsII, LII, XII = as_arrays(d, S["type_II_pred"])
    sI = LI[np.arange(len(ligsI)), np.argmin(XI, axis=1)]
    sII = LII[np.arange(len(ligsII)), np.argmin(XII, axis=1)]

    gain_I = 0.0
    if aug and "primary" in aug:
        gain_I = aug["primary"]["delta_selected"]

    rng = np.random.default_rng(SEED + 9)
    out = {"S_I_now": round(float(sI.mean()), 4), "S_II_now": round(float(sII.mean()), 4),
           "type_I_gain_from_depth": gain_I, "n_I": len(ligsI), "n_II": len(ligsII)}
    for f in (0.0, 0.25, 0.5, 1.0):
        now = f * sI.mean() + (1 - f) * sII.mean()
        after = f * (sI.mean() + gain_I) + (1 - f) * sII.mean()
        draws = np.empty(N_BOOT)
        for b in range(N_BOOT):
            a = sI[rng.integers(0, len(sI), len(sI))].mean()
            c = sII[rng.integers(0, len(sII), len(sII))].mean()
            draws[b] = f * a + (1 - f) * c
        out[f"f={f}"] = {"expected_now": round(float(now), 4),
                         "CI_now": [round(float(np.percentile(draws, 2.5)), 4),
                                    round(float(np.percentile(draws, 97.5)), 4)],
                         "expected_after_depth": round(float(after), 4),
                         "delta": round(float(after - now), 4)}
    if sub:
        out["selected_per_doubling"] = {
            k: sub[k]["selected_per_doubling"]
            for k in ("type_I_pred", "type_II_pred", "all")}
    OUT_PROJ.write_text(json.dumps(out, indent=2))
    return out


def stratum() -> dict:
    """The prediction-side Type I stratum at full depth, with its own null and leverage.

    FINDING 033 stratified on the CRYSTAL label. The label we would actually have at
    submission time is the prediction-side one, and the two differ by a swap (D0R in,
    QDY out), so the stratum's own numbers are re-measured here rather than inherited.
    """
    d = load_pool()
    S = strata()
    rng = np.random.default_rng(SEED + 11)
    out = {}
    for name in ("type_I_pred", "type_II_pred"):
        ligs, L, XE = as_arrays(d, S[name])
        base = L.mean(axis=1)
        selp = L[np.arange(len(ligs)), np.argmin(XE, axis=1)]
        gains = selp - base
        n_lig, n_pose = L.shape
        null = np.empty(N_NULL)
        for b in range(N_NULL):
            pick = np.argmin(rng.random((n_lig, n_pose)), axis=1)
            null[b] = (L[np.arange(n_lig), pick] - base).mean()
        order = np.argsort(gains)[::-1]
        rhos = [stats.spearmanr(XE[i], L[i]).statistic for i in range(n_lig)]
        out[name] = {
            "n": n_lig,
            "gain": round(float(gains.mean()), 4),
            "null_p95": round(float(np.percentile(null, 95)), 4),
            "null_p99": round(float(np.percentile(null, 99)), 4),
            "empirical_p": round(float((null >= gains.mean()).mean()), 4),
            "drop_best_1": round(float(np.delete(gains, order[0]).mean()), 4),
            "drop_best_2": round(float(np.delete(gains, order[:2]).mean()), 4),
            "frac_ligands_positive": round(float((gains > 0).mean()), 4),
            "within_ligand_rho": round(float(np.mean(rhos)), 4),
            "frac_rho_negative": round(float(np.mean(np.array(rhos) < 0)), 4),
            "per_ligand_gain": {l: round(float(g), 4) for l, g in zip(ligs, gains)},
        }
    (DATA_PROCESSED / "type_i_depth_stratum.json").write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["controls", "subsample", "augment", "project",
                                    "stratum"])
    a = ap.parse_args()
    fn = {"controls": controls, "subsample": subsample, "stratum": stratum,
          "augment": augment, "project": project}[a.cmd]
    print(json.dumps(fn(), indent=2, default=str))
