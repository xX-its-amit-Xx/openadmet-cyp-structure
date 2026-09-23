"""FINDING 037 — is the right pose a MINORITY REPORT the consensus down-ranks?

Pre-registered in `docs/PREREG_minority_report.md`, committed at `6f547ef` before a single
endpoint existed.

`FINDING_036` closed tie-breaking by its term oracle and left one question: **60% of the
0.0811 oracle gap sits outside the top 5**, so the question is *"why is the right pose
ranked 8th?"*. 036 also supplied the mechanism to test first — the co-folders SHARE
CYP3A4's 30 degree orientation error (`FINDING_024`), so the consensus is displaced from
the truth, and the best pose may be an outlier the consensus actively down-ranks.

Everything is on disk. Zero new inference, zero downloads, CPU only.

    python scripts/structure/minority_report.py controls   # C-XENG C-SEL C-N2 C-NUM ...
    python scripts/structure/minority_report.py rank       # STEP 1 -- where the oracle sits
    python scripts/structure/minority_report.py cluster    # STEPS 2 + 4 -- minority + oracle
    python scripts/structure/minority_report.py confident  # STEP 3 -- confidently wrong
    python scripts/structure/minority_report.py all
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

TIEBREAK_CACHE = Path("C:/Temp/cyp_tiebreak")        # frames + pose-reference D, verified
CACHE = Path("C:/Temp/cyp_minority")                 # this experiment's own intermediates
FROZEN = DATA_PROCESSED / "reference_set_cyp3a4.npz"

# --- fixed a priori in the pre-registration, section 4 -----------------------------
CUTS_CHAMFER = (0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 2.50, 3.00)
CUTS_ROT = (15.0, 30.0, 45.0, 60.0, 90.0, 120.0)
CUT_MEDIAN = 1.25                                    # the "median cut", for step 3 only
STAKE_HIGH = 0.05                                    # FINDING 036's stake threshold
WRONG_CUT = 0.50                                     # step 3, below this pool's random
NULL_DRAWS = 4000
BOOT = 10_000
SEED = 20260923

OUT_CONTROLS = DATA_PROCESSED / "minority_controls.json"
OUT_RANK = DATA_PROCESSED / "minority_rank.json"
OUT_CLUST = DATA_PROCESSED / "minority_clusters.json"
OUT_CONF = DATA_PROCESSED / "minority_confident.json"
OUT_PERLIG = DATA_PROCESSED / "minority_per_ligand.csv"


# ---------------------------------------------------------------------------
# loading -- the caches FINDING 036 built and verified, re-verified here
# ---------------------------------------------------------------------------

def load_frames() -> dict[str, tuple[list[str], np.ndarray]]:
    d = np.load(TIEBREAK_CACHE / "frames_A.npz", allow_pickle=False)
    ligs = sorted({k.split("__", 1)[1] for k in d.files if k.startswith("xyz__")})
    return {l: ([str(x) for x in d[f"nm__{l}"]], d[f"xyz__{l}"]) for l in ligs}


def load_dmat() -> dict[str, tuple[list[str], np.ndarray]]:
    p = TIEBREAK_CACHE / "dmat_A.npz"
    if not p.exists():
        raise SystemExit(f"missing {p}; run tiebreak_analysis.py frames first")
    d = np.load(p, allow_pickle=False)
    ligs = sorted({k.split("__", 1)[1] for k in d.files if k.startswith("D__")})
    return {l: ([str(x) for x in d[f"nm__{l}"]], d[f"D__{l}"]) for l in ligs}


def truth() -> pd.DataFrame:
    p = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    return p[p.arm == "unsteered"].copy()


# ---------------------------------------------------------------------------
# the two pre-registered distances
# ---------------------------------------------------------------------------

def _chamfer(a: np.ndarray, b: np.ndarray) -> float:
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    return float(0.5 * (d.min(axis=1).mean() + d.min(axis=0).mean()))


def _rot_angle(a: np.ndarray, b: np.ndarray) -> float:
    """FINDING 024's rotation measure: Kabsch after centroid alignment, degrees.

    Identity atom correspondence, which is legal ONLY within one ligand's pool -- all 20
    poses come from one generator with one atom order. C-ATOM checks that.
    """
    ac = a - a.mean(axis=0)
    bc = b - b.mean(axis=0)
    H = ac.T @ bc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return float(np.degrees(np.arccos(np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0))))


def pose_pose(xyz: np.ndarray, kind: str) -> np.ndarray:
    n = len(xyz)
    M = np.zeros((n, n))
    f = _chamfer if kind == "chamfer" else _rot_angle
    for i in range(n):
        for j in range(i + 1, n):
            M[i, j] = M[j, i] = f(np.asarray(xyz[i], float), np.asarray(xyz[j], float))
    return M


def pp_cache() -> dict[str, dict[str, np.ndarray]]:
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / "posepose.npz"
    fr = load_frames()
    if p.exists():
        d = np.load(p, allow_pickle=False)
        return {l: {"chamfer": d[f"C__{l}"], "rot": d[f"R__{l}"]} for l in sorted(fr)}
    out, payload = {}, {}
    for i, (lig, (_, xyz)) in enumerate(sorted(fr.items()), 1):
        C = pose_pose(xyz, "chamfer")
        R = pose_pose(xyz, "rot")
        out[lig] = {"chamfer": C, "rot": R}
        payload[f"C__{lig}"] = C.astype(np.float32)
        payload[f"R__{lig}"] = R.astype(np.float32)
        if i % 20 == 0:
            print(f"  pose-pose {i}/{len(fr)}", flush=True)
    np.savez_compressed(p, **payload)
    return out


# ---------------------------------------------------------------------------
# assemble -- feature rows aligned to truth rows, exactly as FINDING 036 did
# ---------------------------------------------------------------------------

def assemble() -> dict:
    m = load_dmat()
    fr = load_frames()
    t = truth().set_index(["ligand", "sample"])
    out = {}
    dropped_ref = dropped_truth = 0
    for lig, (names, D) in m.items():
        assert fr[lig][0] == names
        try:
            lddt = np.array([float(t.loc[(lig, nm), "lddt_pli"]) for nm in names])
        except KeyError:
            dropped_truth += 1
            continue
        out[lig] = {"names": names, "D": D, "xeng": D.mean(axis=1), "lddt": lddt,
                    "xyz": fr[lig][1]}
    for lig in fr:
        if lig not in m:
            dropped_ref += 1
    return out, {"dropped_no_reference": dropped_ref, "dropped_no_truth": dropped_truth}


def order_of(rec) -> np.ndarray:
    return np.argsort(rec["xeng"], kind="stable")


# ---------------------------------------------------------------------------
# CONTROLS
# ---------------------------------------------------------------------------

def controls() -> dict:
    from scipy.stats import spearmanr
    from cypstruct import xengine as X

    out: dict = {}
    t = truth()
    ship = pd.read_csv(DATA_PROCESSED / "xeng_val87b.csv")
    m = load_dmat()

    rows = [{"ligand": l, "sample": nm, "xeng_rebuilt": D[i].mean()}
            for l, (names, D) in m.items() for i, nm in enumerate(names)]
    reb = pd.DataFrame(rows)
    j = ship.merge(reb, on=["ligand", "sample"], how="outer", indicator=True)
    md = float(np.nanmax(np.abs(j.xeng - j.xeng_rebuilt)))
    out["C_XENG"] = {"shipped_rows": int(len(ship)), "rebuilt_rows": int(len(reb)),
                     "merge_both": int((j._merge == "both").sum()),
                     "merge_left_only": int((j._merge == "left_only").sum()),
                     "merge_right_only": int((j._merge == "right_only").sum()),
                     "max_abs_diff": md, "passes_1e-6": bool(md < 1e-6)}

    d = t.merge(ship, on=["ligand", "sample"], how="left")
    assert d.xeng.notna().all()
    all87 = sorted(d.ligand.unique())
    shipped = X.select(d).set_index("ligand").loc[all87, "lddt_pli"].values
    mine = np.array([float(g.lddt_pli.iloc[int(np.argmin(g.xeng.values))])
                     for _, g in d.groupby("ligand")])
    out["C_SEL"] = {"n": len(all87),
                    "max_abs_diff": float(np.max(np.abs(mine - shipped))),
                    "equal_to_1e-9": bool(np.max(np.abs(mine - shipped)) < 1e-9)}

    sel = float(mine.mean())
    orc = float(d.groupby("ligand").lddt_pli.max().mean())
    rnd = float(d.groupby("ligand").lddt_pli.mean().mean())
    rhos = np.array([spearmanr(g.xeng.values, g.lddt_pli.values).statistic
                     for _, g in d.groupby("ligand")])
    rhos = rhos[np.isfinite(rhos)]
    out["C_N2"] = {"selected": round(sel, 4), "oracle": round(orc, 4),
                   "random": round(rnd, 4), "gain": round(sel - rnd, 4),
                   "within_ligand_rho": round(float(rhos.mean()), 4),
                   "correct_sign_pct": round(100 * float((rhos < 0).mean()), 2),
                   "expected": {"selected": 0.6164, "oracle": 0.6975, "random": 0.5769,
                                "gain": 0.0395, "rho": -0.2582,
                                "correct_sign_pct": 75.86}}
    out["C_N2"]["matches"] = bool(abs(sel - 0.6164) < 5e-4 and abs(orc - 0.6975) < 5e-4
                                  and abs(sel - rnd - 0.0395) < 5e-4)

    z = t[t.lddt_pli == 0]
    out["C_NUM"] = {"rows": int(len(t)), "mapped_true": int(t.mapped.sum()),
                    "mapped_false": int((~t.mapped).sum()),
                    "exact_zero_lddt": int(len(z)),
                    "zero_ligands": z.ligand.value_counts().to_dict(),
                    "zero_bisy_rmsd": [round(float(v), 2) for v in z.bisy_rmsd.values],
                    "verdict": ("ejections, not FINDING_021 numbering"
                                if len(z) and (z.bisy_rmsd > 10).all() else "CHECK")}

    ref = X.load_reference(FROZEN)
    dep = np.array([len(v) for v in ref.values()])
    out["C_REF"] = {"ligands": len(ref), "min": int(dep.min()),
                    "median": float(np.median(dep)), "max": int(dep.max()),
                    "all_at_least_4": bool((dep >= 4).all())}

    fr = load_frames()
    homog = {l: len({np.asarray(v).shape[0] for v in xyz}) == 1
             for l, (_, xyz) in fr.items()}
    out["C_ATOM"] = {"ligands": len(homog),
                     "homogeneous_atom_count": int(sum(homog.values())),
                     "heterogeneous": [l for l, v in homog.items() if not v],
                     "passes": bool(all(homog.values()))}

    data, drops = assemble()
    npose = np.array([len(v["lddt"]) for v in data.values()])
    out["C_FILT"] = {"validation_ligands": int(t.ligand.nunique()),
                     "pool_A_poses": int(len(t)),
                     "pool_A_poses_with_features": int(len(reb)),
                     "ligands_assembled": len(data),
                     "poses_assembled": int(npose.sum()),
                     "poses_per_ligand_min": int(npose.min()),
                     "poses_per_ligand_max": int(npose.max()),
                     **drops}

    OUT_CONTROLS.write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str))
    return out


# ---------------------------------------------------------------------------
# STEP 1 -- where does the oracle pose sit in the xeng ranking?
# ---------------------------------------------------------------------------

def _per_ligand(data: dict) -> pd.DataFrame:
    rows = []
    for lig, rec in sorted(data.items()):
        x, y = rec["xeng"], rec["lddt"]
        o = order_of(rec)
        best = y.max()
        # ties at the maximum -> SMALLEST rank. Conservative against this hypothesis.
        tied = np.flatnonzero(y >= best - 1e-9)
        rank = int(min(np.flatnonzero(np.isin(o, tied))) + 1)
        rows.append({"ligand": lig, "n_poses": len(x),
                     "rank_oracle": rank,
                     "n_tied_at_oracle": int(len(tied)),
                     "sel": float(y[o[0]]), "oracle": float(best),
                     "rand": float(y.mean()),
                     "stake_all": float(best - y[o[0]]),
                     "margin_12": float(x[o[1]] - x[o[0]]),
                     "xeng_mean": float(x.mean()),
                     "n_heavy": int(rec["xyz"].shape[1])})
    return pd.DataFrame(rows)


def rank() -> dict:
    from scipy.stats import wilcoxon, spearmanr
    data, _ = assemble()
    df = _per_ligand(data)
    rng = np.random.default_rng(SEED)

    def dist(sub: pd.DataFrame) -> dict:
        n = len(sub)
        r = sub.rank_oracle.values
        return {"n_ligands": int(n),
                "mean_rank": round(float(r.mean()), 3),
                "median_rank": round(float(np.median(r)), 1),
                "rank_1": f"{int((r == 1).sum())} / {n}",
                "rank_1_pct": round(100 * float((r == 1).mean()), 1),
                "top_3": f"{int((r <= 3).sum())} / {n}",
                "top_3_pct": round(100 * float((r <= 3).mean()), 1),
                "top_5": f"{int((r <= 5).sum())} / {n}",
                "top_5_pct": round(100 * float((r <= 5).mean()), 1),
                "bottom_half_gt10": f"{int((r > 10).sum())} / {n}",
                "bottom_half_pct": round(100 * float((r > 10).mean()), 1),
                "dead_last_20": f"{int((r == 20).sum())} / {n}",
                "dead_last_pct": round(100 * float((r == 20).mean()), 1),
                "histogram_rank_1_to_20": [int((r == k).sum()) for k in range(1, 21)],
                "mean_stake": round(float(sub.stake_all.mean()), 4),
                "wilcoxon_vs_uniform_median_10.5": float(
                    wilcoxon(r - 10.5).pvalue) if n > 5 else None}

    out = {"pool": "A -- 87 ligands x 20 unsteered Boltz-2 poses",
           "incumbent": round(float(df.sel.mean()), 4),
           "oracle": round(float(df.oracle.mean()), 4),
           "random": round(float(df["rand"].mean()), 4),
           "oracle_gap": round(float(df.oracle.mean() - df.sel.mean()), 4),
           "uniform_null_mean_rank": 10.5,
           "ALL": dist(df)}
    out["ALL"]["stake_weighted_mean_rank"] = round(
        float(np.average(df.rank_oracle, weights=np.maximum(df.stake_all, 1e-12))), 3)

    hi = df[df.stake_all > STAKE_HIGH]
    lo = df[df.stake_all <= STAKE_HIGH]
    out["HIGH_stake_gt_0.05"] = dist(hi)
    out["LOW_stake_le_0.05"] = dist(lo)
    out["stake_terciles"] = {}
    tq = pd.qcut(df.stake_all, 3, labels=["T1 least", "T2", "T3 most"], duplicates="drop")
    for g, s in df.groupby(tq, observed=True):
        out["stake_terciles"][str(g)] = {
            "mean_stake": round(float(s.stake_all.mean()), 4), **dist(s)}

    r = spearmanr(df.rank_oracle, df.stake_all)
    out["rho_rank_vs_stake"] = {"rho": round(float(r.statistic), 4), "p": float(r.pvalue)}

    # how much of the 0.0811 gap is reachable inside the top-k -- FINDING 036 reproduced
    reach = {}
    for k in (1, 2, 3, 5, 10):
        v = []
        for lig, rec in sorted(data.items()):
            o = order_of(rec)[:k]
            v.append(float(rec["lddt"][o].max() - rec["lddt"][o[0]]))
        reach[f"k{k}"] = {"ORACLE_topk": round(float(np.mean(v)), 4),
                          "share_of_gap": round(float(np.mean(v)) /
                                                float(df.oracle.mean() - df.sel.mean()), 4)}
    out["reachable_inside_topk"] = reach
    del rng
    OUT_RANK.write_text(json.dumps(out, indent=2, default=str))
    df.to_csv(OUT_PERLIG, index=False)
    print(json.dumps(out, indent=2, default=str))
    return out


# ---------------------------------------------------------------------------
# clustering
# ---------------------------------------------------------------------------

def _labels(M: np.ndarray, cut: float) -> np.ndarray:
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    n = len(M)
    if n < 2:
        return np.ones(n, int)
    Z = linkage(squareform(M, checks=False), method="average")
    return fcluster(Z, t=cut, criterion="distance")


def _oracle_cluster(lab: np.ndarray, x: np.ndarray, y: np.ndarray,
                    o: np.ndarray) -> tuple[float, int, int]:
    """Best achievable when the ORACLE picks the cluster and `xeng` picks inside it.

    Poses in ascending `xeng` order; the first occurrence of each label is that cluster's
    prediction-side representative.
    """
    lo = lab[o]
    _, first = np.unique(lo, return_index=True)
    reps = o[first]
    best = int(reps[int(np.argmax(y[reps]))])
    return float(y[best]), int(len(reps)), best


def _matched_null(lab: np.ndarray, x: np.ndarray, y: np.ndarray, o: np.ndarray,
                  rng, draws: int) -> np.ndarray:
    """Shuffle cluster labels within the ligand, preserving the size multiset."""
    out = np.empty(draws)
    n = len(lab)
    for d in range(draws):
        p = rng.permutation(n)
        lo = lab[p][o]
        _, first = np.unique(lo, return_index=True)
        out[d] = y[o[first]].max()
    return out


def cluster() -> dict:
    from scipy.stats import binomtest, spearmanr, wilcoxon
    data, _ = assemble()
    pp = pp_cache()
    df = _per_ligand(data)
    hi = set(df[df.stake_all > STAKE_HIGH].ligand)
    ligs = sorted(data)
    inc = np.array([data[l]["lddt"][order_of(data[l])[0]] for l in ligs])
    rnd = float(np.mean([data[l]["lddt"].mean() for l in ligs]))
    orc = float(np.mean([data[l]["lddt"].max() for l in ligs]))

    # ---- BAR 1: the FINDING 007 random-feature floor, recomputed on this pool ----
    rng = np.random.default_rng(SEED)
    Y = [data[l]["lddt"] for l in ligs]
    g = np.array([np.mean([float(rng.choice(y)) for y in Y]) - rnd
                  for _ in range(NULL_DRAWS)])
    floor = {"p95": round(float(np.percentile(g, 95)), 4),
             "p99": round(float(np.percentile(g, 99)), 4), "draws": NULL_DRAWS,
             "FINDING_007_quoted_p95": 0.0137}

    res = {"n_ligands": len(ligs), "incumbent": round(float(inc.mean()), 4),
           "random": round(rnd, 4), "POOL_ORACLE": round(orc, 4),
           "oracle_gap_O_C": round(orc - float(inc.mean()), 4),
           "F87_random_feature_floor": floor,
           "cut_grids": {"chamfer_A": list(CUTS_CHAMFER), "rotation_deg": list(CUTS_ROT)},
           "by_cut": {}}

    per_cut_rows = []
    for kind, cuts, unit in (("chamfer", CUTS_CHAMFER, "A"), ("rot", CUTS_ROT, "deg")):
        for cut in cuts:
            key = f"{kind}@{cut:g}{unit}"
            sizes_all, sizes_orc, sizes_inc = [], [], []
            orc_in_largest, inc_in_largest, chance = [], [], []
            rhos, oa, switched, ncl = [], [], [], []
            nullmeans = []
            rng2 = np.random.default_rng(SEED + int(cut * 100) + (0 if kind == "chamfer"
                                                                  else 7919))
            for l in ligs:
                rec = data[l]
                M = pp[l][kind]
                lab = _labels(M, cut)
                x, y = rec["xeng"], rec["lddt"]
                o = order_of(rec)
                n = len(y)
                size = np.array([int((lab == c).sum()) for c in lab])
                big = size.max()
                i_inc = int(o[0])
                best = y.max()
                tied = np.flatnonzero(y >= best - 1e-9)
                # the oracle pose: among ties at the max, the one xeng ranks BEST
                i_orc = int(tied[int(np.argmin(np.array([list(o).index(t)
                                                         for t in tied])))])
                sizes_all.append(size.mean())
                sizes_orc.append(size[i_orc])
                sizes_inc.append(size[i_inc])
                orc_in_largest.append(size[i_orc] == big)
                inc_in_largest.append(size[i_inc] == big)
                chance.append(big / n)
                ncl.append(len(np.unique(lab)))
                if len(np.unique(size)) > 1:
                    rk = np.empty(n)
                    rk[o] = np.arange(1, n + 1)
                    rr = spearmanr(size, rk).statistic
                    if np.isfinite(rr):
                        rhos.append(float(rr))
                v, k_, best_i = _oracle_cluster(lab, x, y, o)
                oa.append(v)
                switched.append(best_i != i_inc)
                nullmeans.append(_matched_null(lab, x, y, o, rng2, NULL_DRAWS))
                per_cut_rows.append({"cut": key, "ligand": l,
                                     "n_clusters": int(len(np.unique(lab))),
                                     "size_oracle_pose": int(size[i_orc]),
                                     "size_incumbent_pose": int(size[i_inc]),
                                     "largest_cluster": int(big),
                                     "oracle_in_largest": bool(size[i_orc] == big),
                                     "O_A_pose_lddt": v,
                                     "incumbent_lddt": float(y[i_inc]),
                                     "high_stake": l in hi})
            oa = np.array(oa)
            sizes_orc = np.array(sizes_orc, float)
            sizes_inc = np.array(sizes_inc, float)
            rhos = np.array(rhos)
            d_oa = oa - inc
            boot = np.array([d_oa[rng2.integers(0, len(d_oa), len(d_oa))].mean()
                             for _ in range(BOOT)])
            NA = np.stack(nullmeans).mean(axis=0) - inc.mean()   # NULL_DRAWS matched draws
            d_sz = sizes_orc - sizes_inc
            bt_l = binomtest(int(np.sum(orc_in_largest)), len(ligs), 0.5)
            hi_idx = np.array([l in hi for l in ligs])

            entry = {
                "n_clusters": {"mean": round(float(np.mean(ncl)), 2),
                               "median": float(np.median(ncl)),
                               "min": int(np.min(ncl)), "max": int(np.max(ncl)),
                               "ligands_single_cluster": int(np.sum(np.array(ncl) == 1))},
                "E2a_oracle_in_largest": {
                    "observed": round(float(np.mean(orc_in_largest)), 4),
                    "matched_chance": round(float(np.mean(chance)), 4),
                    "delta_vs_chance": round(float(np.mean(orc_in_largest)
                                                   - np.mean(chance)), 4),
                    "count": f"{int(np.sum(orc_in_largest))} / {len(ligs)}"},
                "E2b_incumbent_in_largest": {
                    "observed": round(float(np.mean(inc_in_largest)), 4),
                    "matched_chance": round(float(np.mean(chance)), 4),
                    "delta_vs_chance": round(float(np.mean(inc_in_largest)
                                                   - np.mean(chance)), 4),
                    "count": f"{int(np.sum(inc_in_largest))} / {len(ligs)}"},
                "E2c_rho_size_vs_xengrank": {
                    "mean_rho": round(float(rhos.mean()), 4) if len(rhos) else None,
                    "n_ligands_defined": int(len(rhos)),
                    "pct_negative_P2_direction": round(
                        100 * float((rhos < 0).mean()), 2) if len(rhos) else None,
                    "binomial_p": float(binomtest(int((rhos < 0).sum()), len(rhos),
                                                  0.5).pvalue) if len(rhos) else None},
                "E2d_size_oracle_minus_incumbent": {
                    "mean": round(float(d_sz.mean()), 4),
                    "mean_size_oracle_pose": round(float(sizes_orc.mean()), 3),
                    "mean_size_incumbent_pose": round(float(sizes_inc.mean()), 3),
                    "wilcoxon_p": float(wilcoxon(d_sz).pvalue)
                    if np.any(d_sz != 0) else 1.0},
                "binom_oracle_in_largest_vs_half": float(bt_l.pvalue),
                "STEP4_term_oracle": {
                    "O_A_selected": round(float(oa.mean()), 4),
                    "O_A_gain_over_incumbent": round(float(d_oa.mean()), 4),
                    "ci95": [round(float(np.percentile(boot, 2.5)), 4),
                             round(float(np.percentile(boot, 97.5)), 4)],
                    "ligands_switched": int(np.sum(switched)),
                    "share_of_O_C_gap": round(float(d_oa.mean()) / (orc - inc.mean()), 4),
                    "NA_matched_cluster_shuffle_null": {
                        "mean": round(float(NA.mean()), 4),
                        "p95": round(float(np.percentile(NA, 95)), 4),
                        "draws": NULL_DRAWS},
                    "beats_F87_p95": bool(d_oa.mean() > floor["p95"]),
                    "beats_NA_p95": bool(d_oa.mean() > float(np.percentile(NA, 95)))},
            }
            # conditioned on stake -- E2e
            if hi_idx.sum() >= 5:
                entry["E2e_HIGH_stake"] = {
                    "n": int(hi_idx.sum()),
                    "oracle_in_largest": round(float(np.array(orc_in_largest)[hi_idx]
                                                     .mean()), 4),
                    "matched_chance": round(float(np.array(chance)[hi_idx].mean()), 4),
                    "mean_size_oracle_pose": round(float(sizes_orc[hi_idx].mean()), 3),
                    "mean_size_incumbent_pose": round(float(sizes_inc[hi_idx].mean()), 3),
                    "O_A_gain_over_incumbent": round(float(d_oa[hi_idx].mean()), 4)}
            res["by_cut"][key] = entry

    ch = [res["by_cut"][f"chamfer@{c:g}A"]["STEP4_term_oracle"]["O_A_gain_over_incumbent"]
          for c in CUTS_CHAMFER]
    ro = [res["by_cut"][f"rot@{c:g}deg"]["STEP4_term_oracle"]["O_A_gain_over_incumbent"]
          for c in CUTS_ROT]
    res["STEP4_SUMMARY"] = {
        "O_A_chamfer_by_cut": dict(zip([f"{c:g}" for c in CUTS_CHAMFER],
                                       [round(v, 4) for v in ch])),
        "O_A_rotation_by_cut": dict(zip([f"{c:g}" for c in CUTS_ROT],
                                        [round(v, 4) for v in ro])),
        "O_A_max_over_chamfer_grid": round(float(np.max(ch)), 4),
        "O_A_median_over_chamfer_grid": round(float(np.median(ch)), 4),
        "O_A_max_over_rotation_grid": round(float(np.max(ro)), 4),
        "O_C_full_pose_oracle": round(orc - float(inc.mean()), 4),
        "F87_p95": floor["p95"],
        "cuts_beating_F87": [k for k, v in res["by_cut"].items()
                             if v["STEP4_term_oracle"]["beats_F87_p95"]],
        "cuts_beating_NA": [k for k, v in res["by_cut"].items()
                            if v["STEP4_term_oracle"]["beats_NA_p95"]],
    }
    mx = float(np.max(ch))
    md = float(np.median(ch))
    nbeat = sum(1 for c in CUTS_CHAMFER
                if res["by_cut"][f"chamfer@{c:g}A"]["STEP4_term_oracle"]["beats_NA_p95"])
    if mx < floor["p95"]:
        v = ("CLOSED / REFUTED -- the max of the term oracle over the whole cut grid is "
             "below the recomputed random-feature floor. No selector is built.")
    elif md >= 2 * floor["p95"] and nbeat >= len(CUTS_CHAMFER) / 2:
        v = ("LICENSED -- median-over-grid ceiling is at least twice the floor AND beats "
             "the matched cluster-shuffle null at half the cuts. A follow-up "
             "pre-registration may build a rule. Nothing is built in this sitting.")
    else:
        v = ("MEASURED -- the ceiling clears the floor at some cuts but fails the "
             "pre-registered LICENSED clause. Report and stop. No selector.")
    res["VERDICT"] = {"rule": v, "O_A_max": round(mx, 4),
                      "O_A_median": round(md, 4), "F87_p95": floor["p95"],
                      "cuts_beating_matched_null": nbeat, "of_cuts": len(CUTS_CHAMFER)}

    OUT_CLUST.write_text(json.dumps(res, indent=2, default=str))
    pd.DataFrame(per_cut_rows).to_csv(
        DATA_PROCESSED / "minority_cluster_per_ligand.csv", index=False)
    print(json.dumps(res["STEP4_SUMMARY"], indent=2, default=str))
    print(json.dumps(res["VERDICT"], indent=2, default=str))
    return res


# ---------------------------------------------------------------------------
# STEP 3 -- when is the consensus CONFIDENTLY wrong?
# ---------------------------------------------------------------------------

def confident() -> dict:
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    data, _ = assemble()
    pp = pp_cache()
    df = _per_ligand(data)
    lab = pd.read_csv(DATA_PROCESSED / "binding_mode_labels_cyp3a4.csv")
    lab = lab.set_index("id")

    rows = []
    for lig, rec in sorted(data.items()):
        M = pp[lig]["chamfer"]
        n = len(M)
        spread = float(M[np.triu_indices(n, 1)].mean())
        ncl = int(len(np.unique(_labels(M, CUT_MEDIAN))))
        r = df[df.ligand == lig].iloc[0]
        rows.append({"ligand": lig, "spread": spread, "n_clusters": ncl,
                     "margin_12": float(r.margin_12), "xeng_mean": float(r.xeng_mean),
                     "n_heavy": int(r.n_heavy), "sel": float(r.sel),
                     "oracle": float(r.oracle), "stake_all": float(r.stake_all),
                     "pred_fe_donor_median": float(lab.loc[lig, "pred_fe_donor_median"])
                     if lig in lab.index else np.nan,
                     "pred_frac_coordinated": float(lab.loc[lig, "pred_frac_coordinated"])
                     if lig in lab.index else np.nan})
    c = pd.DataFrame(rows)
    t1 = c.spread.quantile(1 / 3)
    c["tight"] = c.spread <= t1
    c["wrong"] = c.sel < WRONG_CUT

    tab = {}
    for tg in (True, False):
        for wr in (True, False):
            k = f"{'TIGHT' if tg else 'WIDE'}_{'WRONG' if wr else 'RIGHT'}"
            s = c[(c.tight == tg) & (c.wrong == wr)]
            tab[k] = {"n": int(len(s)),
                      "mean_selected": round(float(s.sel.mean()), 4) if len(s) else None,
                      "mean_oracle": round(float(s.oracle.mean()), 4) if len(s) else None,
                      "mean_spread": round(float(s.spread.mean()), 3) if len(s) else None,
                      "ligands": sorted(s.ligand.tolist())}

    covs = ["spread", "n_clusters", "margin_12", "xeng_mean", "n_heavy",
            "pred_fe_donor_median", "pred_frac_coordinated"]
    tw = c[c.tight & c.wrong]
    tr = c[c.tight & ~c.wrong]
    tests, ps = {}, []
    for cv in covs:
        a = tw[cv].dropna().values
        b = tr[cv].dropna().values
        if len(a) < 3 or len(b) < 3:
            tests[cv] = {"n_wrong": len(a), "n_right": len(b), "p": None,
                         "note": "too few"}
            ps.append(1.0)
            continue
        u = mannwhitneyu(a, b)
        tests[cv] = {"n_wrong": len(a), "n_right": len(b),
                     "median_tight_wrong": round(float(np.median(a)), 4),
                     "median_tight_right": round(float(np.median(b)), 4),
                     "p_raw": float(u.pvalue)}
        ps.append(float(u.pvalue))
    holm = multipletests(ps, method="holm")[1]
    for cv, p in zip(covs, holm):
        tests[cv]["p_holm"] = float(p)
        tests[cv]["survives_holm_0.05"] = bool(p < 0.05)

    out = {"definitions": {"TIGHT": f"bottom tercile of mean pairwise Chamfer "
                                    f"(<= {t1:.3f} A)",
                           "WRONG": f"incumbent-selected LDDT-PLI < {WRONG_CUT}",
                           "n_clusters at": f"chamfer cut {CUT_MEDIAN} A"},
           "n_ligands": int(len(c)),
           "tight_tercile_cut_A": round(float(t1), 4),
           "rate_wrong_overall": round(float(c.wrong.mean()), 4),
           "rate_wrong_given_tight": round(float(c[c.tight].wrong.mean()), 4),
           "rate_wrong_given_wide": round(float(c[~c.tight].wrong.mean()), 4),
           "two_by_two": tab,
           "separators_tight_wrong_vs_tight_right": tests,
           "any_survives_holm": bool(any(v.get("survives_holm_0.05") for v in
                                         tests.values()))}
    OUT_CONF.write_text(json.dumps(out, indent=2, default=str))
    c.to_csv(DATA_PROCESSED / "minority_confident_per_ligand.csv", index=False)
    print(json.dumps(out, indent=2, default=str))
    return out


# ---------------------------------------------------------------------------
# ADDENDUM -- declared POST-HOC in FINDING 037 section 5. Three de-confounding
# checks that the pre-registration did not name. Reported as diagnostics, and
# NOT eligible to change the pre-registered verdict.
# ---------------------------------------------------------------------------

def addendum() -> dict:
    from scipy.stats import spearmanr
    data, _ = assemble()
    pp = pp_cache()
    df = _per_ligand(data)
    ligs = sorted(data)
    rng = np.random.default_rng(SEED + 1)
    out: dict = {}

    # A1 -- rank-1 IMPLIES stake 0, so rho(rank, stake) is tautological at the bottom.
    # Recompute excluding the ligands the incumbent already solves.
    nz = df[df.stake_all > 1e-9]
    r_all = spearmanr(df.rank_oracle, df.stake_all)
    r_nz = spearmanr(nz.rank_oracle, nz.stake_all)
    out["A1_rho_rank_vs_stake"] = {
        "all_87": {"rho": round(float(r_all.statistic), 4), "p": float(r_all.pvalue),
                   "n": int(len(df))},
        "stake_gt_0_only": {"rho": round(float(r_nz.statistic), 4),
                            "p": float(r_nz.pvalue), "n": int(len(nz)),
                            "note": "rank_oracle == 1 forces stake == 0; this removes it"},
    }

    # A2 -- a proper null for E2a: replace the ORACLE pose by a RANDOM pose of the same
    # ligand and recompute the in-largest-cluster rate. Poisson-binomial, by simulation.
    out["A2_E2a_null"] = {}
    for kind, cuts, unit in (("chamfer", CUTS_CHAMFER, "A"), ("rot", CUTS_ROT, "deg")):
        for cut in cuts:
            inlarge_o, inlarge_i, memb = [], [], []
            hi = []
            for l in ligs:
                rec = data[l]
                lab = _labels(pp[l][kind], cut)
                y, o = rec["lddt"], order_of(rec)
                size = np.array([int((lab == c).sum()) for c in lab])
                big = size.max()
                tied = np.flatnonzero(y >= y.max() - 1e-9)
                i_orc = int(tied[int(np.argmin([list(o).index(t) for t in tied]))])
                inlarge_o.append(size[i_orc] == big)
                inlarge_i.append(size[int(o[0])] == big)
                memb.append(size == big)
                hi.append(float(df[df.ligand == l].stake_all.iloc[0]) > STAKE_HIGH)
            memb = np.array(memb)
            hi = np.array(hi)
            draws = np.array([np.mean([m[rng.integers(0, len(m))] for m in memb])
                              for _ in range(NULL_DRAWS)])
            obs = float(np.mean(inlarge_o))
            dr_hi = np.array([np.mean([m[rng.integers(0, len(m))] for m in memb[hi]])
                              for _ in range(NULL_DRAWS)])
            out["A2_E2a_null"][f"{kind}@{cut:g}{unit}"] = {
                "observed_oracle": round(obs, 4),
                "null_mean": round(float(draws.mean()), 4),
                "null_p05": round(float(np.percentile(draws, 5)), 4),
                "p_one_sided_below": float(np.mean(draws <= obs)),
                "observed_incumbent": round(float(np.mean(inlarge_i)), 4),
                "HIGH_stake": {
                    "n": int(hi.sum()),
                    "observed_oracle": round(float(np.mean(np.array(inlarge_o)[hi])), 4),
                    "null_mean": round(float(dr_hi.mean()), 4),
                    "p_one_sided_below": float(
                        np.mean(dr_hi <= np.mean(np.array(inlarge_o)[hi]))),
                    "minority_rate_observed": round(
                        1 - float(np.mean(np.array(inlarge_o)[hi])), 4),
                    "minority_rate_chance": round(1 - float(dr_hi.mean()), 4)},
            }

    # A3 -- the orientation-specific excess: O-A minus the matched cluster-shuffle null.
    # O-A is a best-of-k oracle and rises with k by construction; this is the part that
    # is attributable to clustering by ORIENTATION rather than into same-sized groups.
    c = json.loads(OUT_CLUST.read_text())
    ex = {}
    for k, v in c["by_cut"].items():
        t = v["STEP4_term_oracle"]
        ex[k] = {"O_A": t["O_A_gain_over_incumbent"],
                 "NA_mean": t["NA_matched_cluster_shuffle_null"]["mean"],
                 "excess": round(t["O_A_gain_over_incumbent"]
                                 - t["NA_matched_cluster_shuffle_null"]["mean"], 4),
                 "n_clusters_mean": v["n_clusters"]["mean"]}
    ch = [v["excess"] for k, v in ex.items() if k.startswith("chamfer")]
    ro = [v["excess"] for k, v in ex.items() if k.startswith("rot")]
    out["A3_orientation_specific_excess"] = {
        "per_cut": ex,
        "chamfer_median": round(float(np.median(ch)), 4),
        "chamfer_max": round(float(np.max(ch)), 4),
        "rotation_median": round(float(np.median(ro)), 4),
        "rotation_max": round(float(np.max(ro)), 4),
        "F87_p95": c["F87_random_feature_floor"]["p95"],
        "reading": ("O-A converges on the FULL pose oracle as clusters approach "
                    "singletons, where 'take the minority cluster' carries no "
                    "information at all. The excess over a size-matched shuffle is the "
                    "part attributable to orientation."),
    }

    # A4 -- does the MAJORITY VOTE pay? Compare the incumbent (which is by construction
    # one cluster's representative) against the MEAN over clusters of that cluster's
    # xeng-representative, i.e. picking an orientation mode at random. If the two are
    # equal the majority vote buys nothing; if the incumbent wins, it pays even though
    # the oracle pose is a minority report.
    from scipy.stats import wilcoxon as _wx
    out["A4_majority_vote_pays"] = {}
    for kind, cuts, unit in (("chamfer", CUTS_CHAMFER, "A"), ("rot", CUTS_ROT, "deg")):
        for cut in cuts:
            inc_v, rndc_v, worst_v, nk = [], [], [], []
            for l in ligs:
                rec = data[l]
                lab = _labels(pp[l][kind], cut)
                y, o = rec["lddt"], order_of(rec)
                lo = lab[o]
                _, first = np.unique(lo, return_index=True)
                reps = o[first]
                inc_v.append(float(y[int(o[0])]))
                rndc_v.append(float(y[reps].mean()))
                worst_v.append(float(y[reps].min()))
                nk.append(len(reps))
            inc_v = np.array(inc_v); rndc_v = np.array(rndc_v); worst_v = np.array(worst_v)
            d = inc_v - rndc_v
            out["A4_majority_vote_pays"][f"{kind}@{cut:g}{unit}"] = {
                "mean_clusters": round(float(np.mean(nk)), 2),
                "incumbent": round(float(inc_v.mean()), 4),
                "random_cluster_representative": round(float(rndc_v.mean()), 4),
                "worst_cluster_representative": round(float(worst_v.mean()), 4),
                "incumbent_minus_random_cluster": round(float(d.mean()), 4),
                "wilcoxon_p": float(_wx(d).pvalue) if np.any(d != 0) else 1.0}

    (DATA_PROCESSED / "minority_addendum.json").write_text(
        json.dumps(out, indent=2, default=str))
    print(json.dumps({"A4": out["A4_majority_vote_pays"],
                      "A1": out["A1_rho_rank_vs_stake"],
                      "A3": {k: v for k, v in out["A3_orientation_specific_excess"].items()
                             if k != "per_cut"}}, indent=2, default=str))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["controls", "rank", "cluster", "confident",
                                      "addendum", "all"])
    a = ap.parse_args()
    if a.stage in ("controls", "all"):
        controls()
    if a.stage in ("rank", "all"):
        rank()
    if a.stage in ("cluster", "all"):
        cluster()
    if a.stage in ("confident", "all"):
        confident()
    if a.stage in ("addendum", "all"):
        addendum()
    return 0


if __name__ == "__main__":
    os.environ.setdefault("PYTHONWARNINGS", "ignore")
    raise SystemExit(main())
