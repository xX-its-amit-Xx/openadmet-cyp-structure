"""FINDING 035 — does depth at MATCHED conditioning convert, as the law predicts?

Pre-registered in `docs/PREREG_matched_depth.md`. Every stage below is named there and
nothing was added after a score was seen.

The arithmetic is deliberately `FINDING_034`'s, re-run against a different new arm, so the
two are comparable line for line: same stratum, same selector, same exact-expectation
random baseline, same 256-draw subsampling, same null recomputed INSIDE the n=14 stratum
because the `FINDING_007` floor is n-dependent (it is +0.043 here, triple the pooled
+0.0137).

    python scripts/structure/matched_depth_analysis.py controls
    python scripts/structure/matched_depth_analysis.py distinct
    python scripts/structure/matched_depth_analysis.py augment
    python scripts/structure/matched_depth_analysis.py all
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
sys.path.insert(0, str(REPO / "scripts" / "cofold"))
sys.path.insert(0, str(REPO / "scripts" / "explorer"))

from cypstruct.paths import DATA_PROCESSED, SCRATCH  # noqa: E402

COORD_MAX = 2.6
DEDUPE_TOL = 0.05          # cypstruct.xengine._dedupe default, A
DRAWS = 256
NULL_DRAWS = 4000
BOOT = 10_000
RNG = np.random.default_rng(20260922)

POSES_CSV = DATA_PROCESSED / "matched_depth_poses.csv"
OUT_CONTROLS = DATA_PROCESSED / "matched_depth_controls.json"
OUT_DISTINCT = DATA_PROCESSED / "matched_depth_distinct.json"
OUT_AUGMENT = DATA_PROCESSED / "matched_depth_augment.json"
OUT_PERLIG = DATA_PROCESSED / "matched_depth_per_ligand.csv"


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

def stratum() -> list[str]:
    L = pd.read_csv(DATA_PROCESSED / "binding_mode_labels_cyp3a4.csv")
    return sorted(L[L.pred_fe_donor_median > COORD_MAX].id.tolist())


def existing_pool() -> pd.DataFrame:
    """The 20-pose Modal pool, unsteered arm, with the SHIPPED xeng column attached."""
    p = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    x = pd.read_csv(DATA_PROCESSED / "xeng_val87b.csv")
    d = p[p.arm == "unsteered"].merge(x, on=["ligand", "sample"], how="left")
    assert d.xeng.notna().all(), "xeng did not merge onto every pool row"
    return d


def new_pool() -> pd.DataFrame:
    return pd.read_csv(POSES_CSV)


# ---------------------------------------------------------------------------
# the selector, and the three summary statistics
# ---------------------------------------------------------------------------

def _pick(sub: pd.DataFrame) -> float:
    """Selected LDDT-PLI: the pose with the lowest xeng.

    `cypstruct.xengine.select` maximises `-z(xeng)`; z-scoring is monotone within a
    ligand, so this is argmin(xeng) exactly. Asserted, not assumed, in control N2b.
    """
    return float(sub.lddt_pli.iloc[int(np.argmin(sub.xeng.values))])


def summarise(df: pd.DataFrame, ligs: list[str]) -> dict:
    orc, sel, rnd = [], [], []
    for lid in ligs:
        s = df[df.ligand == lid]
        orc.append(s.lddt_pli.max())
        sel.append(_pick(s))
        rnd.append(s.lddt_pli.mean())          # exact per-ligand expectation
    return {"oracle": float(np.mean(orc)), "selected": float(np.mean(sel)),
            "random": float(np.mean(rnd)), "gain": float(np.mean(sel) - np.mean(rnd)),
            "per_ligand_selected": sel, "per_ligand_oracle": orc,
            "per_ligand_random": rnd}


def within_ligand_rho(df: pd.DataFrame, ligs: list[str]) -> tuple[float, float]:
    from scipy.stats import spearmanr
    rs = []
    for lid in ligs:
        s = df[df.ligand == lid]
        if s.lddt_pli.nunique() < 2 or s.xeng.nunique() < 2:
            continue
        r = spearmanr(s.xeng, s.lddt_pli).statistic
        if np.isfinite(r):
            rs.append(r)
    rs = np.array(rs)
    return float(rs.mean()), float((rs < 0).mean())


# ---------------------------------------------------------------------------
# A1 -- controls
# ---------------------------------------------------------------------------

def controls() -> dict:
    from cypstruct import xengine as X

    pool = existing_pool()
    all87 = sorted(pool.ligand.unique())
    ligs = stratum()

    # N2 -- the shipped column reproduces, on all 87
    s87 = summarise(pool, all87)
    rho87, sign87 = within_ligand_rho(pool, all87)
    n2 = {"selected": round(s87["selected"], 4), "oracle": round(s87["oracle"], 4),
          "random": round(s87["random"], 4), "gain": round(s87["gain"], 4),
          "within_ligand_rho": round(rho87, 4), "correct_sign_pct": round(100 * sign87, 2),
          "expected": {"selected": 0.6164, "oracle": 0.6975, "random": 0.5769,
                       "gain": 0.0395, "rho": -0.2582, "correct_sign_pct": 75.86}}
    n2["matches_finding_011_033"] = bool(
        abs(n2["selected"] - 0.6164) < 5e-4 and abs(n2["oracle"] - 0.6975) < 5e-4
        and abs(n2["gain"] - 0.0395) < 5e-4)

    # N2b -- argmin(xeng) IS cypstruct.xengine.select(), asserted to 1e-9
    shipped = X.select(pool)
    mine = np.array([_pick(pool[pool.ligand == l]) for l in all87])
    ship = shipped.set_index("ligand").loc[all87, "lddt_pli"].values
    n2b = {"max_abs_diff": float(np.max(np.abs(mine - ship))),
           "equal_to_1e-9": bool(np.max(np.abs(mine - ship)) < 1e-9),
           "n": len(all87)}

    # baseline inside the stratum, depth 20
    base = summarise(pool[pool.ligand.isin(ligs)], ligs)
    rho14, sign14 = within_ligand_rho(pool[pool.ligand.isin(ligs)], ligs)

    # pool diversity of the EXISTING pool -- 20 samples must be 20 poses
    g = pool[pool.ligand.isin(ligs)].groupby("ligand")
    div = {"unique_xeng_per_ligand": sorted(set(g.xeng.nunique().tolist())),
           "median_within_ligand_lddt_sd": round(float(g.lddt_pli.std().median()), 4),
           "ligands_with_sd_zero": int((g.lddt_pli.std() == 0).sum())}

    out = {"N2_shipped_column": n2, "N2b_argmin_identity": n2b,
           "N5_filters": {"validation_ligands": int(pool.ligand.nunique()),
                          "pass_pred_fe_donor_median_gt_2.6": len(ligs),
                          "excluded_predicted_type_II": int(pool.ligand.nunique()) - len(ligs),
                          "stratum": ligs},
           "existing_pool_diversity_in_stratum": div,
           "baseline_depth20_stratum": {
               k: round(v, 4) for k, v in base.items() if not k.startswith("per_")},
           "baseline_rho": round(rho14, 4),
           "baseline_correct_sign_pct": round(100 * sign14, 2)}

    # N1 / C7 / N4 on the NEW arm, if it has been collected yet
    if POSES_CSV.exists():
        nw = new_pool()
        out["N1_numbering"] = {
            "n_poses": int(len(nw)),
            "renumber_offsets": {str(k): int(v)
                                 for k, v in nw.renumber_offset.value_counts().items()},
            "min_renumber_identity": float(nw.renumber_identity.min()),
            "poses_below_0.95_identity": int((nw.renumber_identity < 0.95).sum()),
            "mapped_true": int(nw.mapped.sum()), "mapped_false": int((~nw.mapped).sum()),
        }
        z = nw[nw.lddt_pli == 0]
        out["C7_exact_zero_lddt"] = {
            "n": int(len(z)),
            "ligands": z.ligand.value_counts().to_dict(),
            "bisy_rmsd_range": [round(float(z.bisy_rmsd.min()), 2),
                                round(float(z.bisy_rmsd.max()), 2)] if len(z) else None,
            "fe_donor_range": [round(float(z.fe_donor_dist.min()), 2),
                               round(float(z.fe_donor_dist.max()), 2)] if len(z) else None,
            "all_offset_zero": bool((z.renumber_offset == 0).all()) if len(z) else None,
            "all_mapped": bool(z.mapped.all()) if len(z) else None,
            "verdict": ("ejections, not FINDING_021 numbering" if len(z) and
                        (z.renumber_offset == 0).all() else
                        "none" if not len(z) else "CHECK"),
        }
        out["N4_md5_distinct"] = {
            "rows": int(len(nw)), "distinct_md5": int(nw.md5.nunique()),
            "per_ligand_distinct_md5": nw.groupby("ligand").md5.nunique().to_dict(),
        }
        out["A0_new_arm_quality"] = {
            "new_arm_pool_mean": round(float(nw.lddt_pli.mean()), 4),
            "existing_pool_mean_same_ligands": round(base["random"], 4),
            "delta": round(float(nw.lddt_pli.mean()) - base["random"], 4),
            "matched_within_0.02": bool(
                abs(float(nw.lddt_pli.mean()) - base["random"]) <= 0.02),
            "finding_034_openprotein_delta_for_reference": -0.057,
        }
    OUT_CONTROLS.write_text(json.dumps(out, indent=2, default=str))
    return out


# ---------------------------------------------------------------------------
# N4 -- DISTINCT poses, by geometry, not by job
# ---------------------------------------------------------------------------

def distinct() -> dict:
    """Dedupe the new poses on ligand coordinates in their own heme frame.

    md5 catches byte-identical files; this catches a pose that was written twice with
    different floating-point noise. `FINDING_034` found 39% of one arm's "replicates" were
    duplicates and 196 of 196 complexes byte-identical across fourteen replicate jobs, so
    the count is never taken from the job count.
    """
    from cypstruct import pose as P
    from cypstruct import xengine as X
    from boltz_depth import POSE_DIR

    flat = POSE_DIR / "stratum_flat"
    nw = new_pool()
    keep: dict[str, list[str]] = {}
    stats = {}
    for lid, grp in nw.groupby("ligand"):
        coords, names = [], []
        for s in grp["sample"]:
            f = flat / f"{lid}__{s}.cif"
            if not f.exists():
                continue
            try:
                hf = X.in_heme_frame(P.load_structure(f))
            except Exception:
                continue
            if hf is None:
                continue
            coords.append(hf)
            names.append(s)
        # exactly cypstruct.xengine._dedupe's rule, kept inline so the surviving NAMES
        # can be recorded alongside the surviving coordinates
        kept_xyz: list[np.ndarray] = []
        kept_names: list[str] = []
        for c, n in zip(coords, names):
            if not any(len(q) == len(c) and float(np.abs(q - c).max()) < DEDUPE_TOL
                       for q in kept_xyz):
                kept_xyz.append(c)
                kept_names.append(n)
        keep[lid] = kept_names
        stats[lid] = {"collected": len(names), "distinct": len(kept_names)}
    out = {"tol_angstrom": DEDUPE_TOL,
           "per_ligand": stats,
           "total_collected": int(sum(v["collected"] for v in stats.values())),
           "total_distinct": int(sum(v["distinct"] for v in stats.values())),
           "min_distinct": int(min(v["distinct"] for v in stats.values())) if stats else 0,
           "mean_distinct": round(float(np.mean([v["distinct"] for v in stats.values()])), 2)
           if stats else 0,
           "keep": keep}
    OUT_DISTINCT.write_text(json.dumps(out, indent=2))
    return {k: v for k, v in out.items() if k != "keep"}


# ---------------------------------------------------------------------------
# A2 -- the noise floor, recomputed INSIDE this stratum
# ---------------------------------------------------------------------------

def noise_floor(df: pd.DataFrame, ligs: list[str], draws: int = NULL_DRAWS) -> dict:
    """A RANDOM feature's gain over random selection, at this n and this depth."""
    lig_scores = [df[df.ligand == l].lddt_pli.values for l in ligs]
    means = np.array([v.mean() for v in lig_scores])
    gains = np.empty(draws)
    for i in range(draws):
        pick = np.array([v[RNG.integers(len(v))] for v in lig_scores])
        gains[i] = pick.mean() - means.mean()
    return {"n": len(ligs), "draws": draws,
            "p95": round(float(np.percentile(gains, 95)), 4),
            "p99": round(float(np.percentile(gains, 99)), 4)}


# ---------------------------------------------------------------------------
# A3 -- the primary endpoint
# ---------------------------------------------------------------------------

def _augment_once(pool_e: pd.DataFrame, pool_n: pd.DataFrame, ligs: list[str],
                  keep: dict[str, list[str]], k: int, draws: int) -> dict:
    """All existing poses kept, k distinct new poses added. Per-ligand means over draws."""
    sel = np.zeros((draws, len(ligs)))
    orc = np.zeros((draws, len(ligs)))
    rnd = np.zeros((draws, len(ligs)))
    took_new = np.zeros(len(ligs))
    for j, lid in enumerate(ligs):
        e = pool_e[pool_e.ligand == lid][["lddt_pli", "xeng"]].values
        avail = keep.get(lid, [])
        nsub = pool_n[(pool_n.ligand == lid) & (pool_n["sample"].isin(avail))]
        n = nsub[["lddt_pli", "xeng"]].values
        for d in range(draws):
            if k > 0 and len(n):
                idx = RNG.choice(len(n), size=min(k, len(n)), replace=False)
                cand = np.vstack([e, n[idx]])
                is_new = np.concatenate([np.zeros(len(e)), np.ones(len(idx))])
            else:
                cand, is_new = e, np.zeros(len(e))
            w = int(np.argmin(cand[:, 1]))
            sel[d, j] = cand[w, 0]
            orc[d, j] = cand[:, 0].max()
            rnd[d, j] = cand[:, 0].mean()
            took_new[j] += is_new[w]
    return {"selected": float(sel.mean()), "oracle": float(orc.mean()),
            "random": float(rnd.mean()),
            "per_ligand_selected": sel.mean(axis=0).tolist(),
            "per_ligand_oracle": orc.mean(axis=0).tolist(),
            "per_ligand_random": rnd.mean(axis=0).tolist(),
            "frac_draws_took_new": (took_new / draws).tolist()}


def augment() -> dict:
    from scipy.stats import wilcoxon

    pool_e = existing_pool()
    pool_n = new_pool()
    ligs = stratum()
    pool_e = pool_e[pool_e.ligand.isin(ligs)]
    dist = json.loads(OUT_DISTINCT.read_text())
    keep = dist["keep"]
    kmax = min(len(v) for v in keep.values())

    rungs = sorted({0, 1, 2, 4, 8, max(1, kmax // 2), kmax})
    rungs = [r for r in rungs if r <= kmax]
    curve = {}
    for k in rungs:
        curve[k] = _augment_once(pool_e, pool_n, ligs, keep, k, DRAWS)

    b0, bk = curve[0], curve[kmax]
    d_sel = np.array(bk["per_ligand_selected"]) - np.array(b0["per_ligand_selected"])
    d_orc = np.array(bk["per_ligand_oracle"]) - np.array(b0["per_ligand_oracle"])

    def boot(v):
        m = np.array([RNG.choice(v, len(v), replace=True).mean() for _ in range(BOOT)])
        return [round(float(np.percentile(m, 2.5)), 4),
                round(float(np.percentile(m, 97.5)), 4)]

    ties = int(np.sum(np.abs(d_sel) < 1e-9))
    better = int(np.sum(d_sel > 1e-9))
    worse = int(np.sum(d_sel < -1e-9))
    try:
        w = wilcoxon(d_sel) if (better + worse) else None
        wp = float(w.pvalue) if w is not None else float("nan")
    except ValueError:
        wp = float("nan")

    # the rate, over the actual doubling that was bought
    depth0, depthk = 20, 20 + kmax
    doublings = np.log2(depthk / depth0)
    rate_sel = (bk["selected"] - b0["selected"]) / doublings if doublings > 0 else np.nan
    rate_orc = (bk["oracle"] - b0["oracle"]) / doublings if doublings > 0 else np.nan

    # filter the new arm to its DISTINCT poses PER LIGAND. A global `isin` over the union
    # of kept sample names would be wrong: sample names ("s101_m3") repeat across ligands,
    # so a pose deduped out of ligand B would be readmitted because ligand A kept the same
    # index. That would quietly double-count a duplicate in rho and in the noise floor.
    keep_pairs = {(lid, s) for lid, v in keep.items() for s in v}
    nsel = pool_n[[(r.ligand, r.sample) in keep_pairs
                   for r in pool_n.itertuples()]]
    aug_df = pd.concat([pool_e, nsel], ignore_index=True)
    rho_a, sign_a = within_ligand_rho(aug_df, ligs)

    out = {
        "stratum_n": len(ligs), "ligands": ligs,
        "k_matched": kmax, "depth_before": depth0, "depth_after": depthk,
        "rungs": {str(k): {kk: round(vv, 4) for kk, vv in v.items()
                           if not kk.startswith(("per_", "frac_"))}
                  for k, v in curve.items()},
        "primary": {
            "delta_selected": round(float(d_sel.mean()), 4),
            "delta_selected_CI95": boot(d_sel),
            "delta_oracle": round(float(d_orc.mean()), 4),
            "delta_oracle_CI95": boot(d_orc),
            "delta_random": round(bk["random"] - b0["random"], 4),
            "wilcoxon_p": round(wp, 4) if np.isfinite(wp) else None,
            "ligands_better": better, "ligands_worse": worse, "ligands_tied": ties,
            "ligands_where_selector_ever_took_a_new_pose":
                int(np.sum(np.array(bk["frac_draws_took_new"]) > 0)),
            "mean_frac_draws_took_new": round(
                float(np.mean(bk["frac_draws_took_new"])), 4),
        },
        "rate_per_doubling": {
            "selected": round(float(rate_sel), 4),
            "oracle": round(float(rate_orc), 4),
            "conversion_pct": round(100 * float(rate_sel / rate_orc), 1)
            if rate_orc else None,
            "finding_034_type_I_selected": 0.0127,
            "finding_034_type_I_oracle": 0.0428,
            "finding_004_benchmark": 0.0125,
        },
        "augmented_within_ligand_rho": round(rho_a, 4),
        "augmented_correct_sign_pct": round(100 * sign_a, 2),
        "noise_floor_in_stratum": noise_floor(pool_e, ligs),
        "noise_floor_augmented": noise_floor(aug_df, ligs),
    }

    pd.DataFrame({
        "ligand": ligs,
        "sel_before": b0["per_ligand_selected"], "sel_after": bk["per_ligand_selected"],
        "d_sel": d_sel,
        "orc_before": b0["per_ligand_oracle"], "orc_after": bk["per_ligand_oracle"],
        "d_orc": d_orc,
        "rnd_before": b0["per_ligand_random"], "rnd_after": bk["per_ligand_random"],
        "frac_draws_took_new": bk["frac_draws_took_new"],
    }).to_csv(OUT_PERLIG, index=False)

    OUT_AUGMENT.write_text(json.dumps(out, indent=2, default=str))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["controls", "distinct", "augment", "all"])
    a = ap.parse_args()
    if a.cmd in ("controls", "all"):
        print(json.dumps(controls(), indent=2, default=str))
    if a.cmd in ("distinct", "all"):
        print(json.dumps(distinct(), indent=2, default=str))
    if a.cmd in ("augment", "all"):
        print(json.dumps(augment(), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
