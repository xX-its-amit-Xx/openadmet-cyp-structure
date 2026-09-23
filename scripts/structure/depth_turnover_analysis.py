"""FINDING 040 - does the SELECTED depth curve turn over? Second matched-depth stratum.

Pre-registered in `docs/PREREG_depth_turnover.md`, committed at `c4a18b8` before a single
Boltz job was submitted. No ligand, rung, threshold or acceptance clause was moved after.

`FINDING_039` s3b recomputed `FINDING_035`'s depth curve in CLOSED FORM and found the
selected score peaks at pool 33 and falls by 40 while the oracle stays monotone -- and
that the whole effect is ONE ligand (08J), with ten of fourteen unable to move at all. It
asked for the same exact curve on a second matched-depth stratum of >= 30 ligands.

This is that, on n = 73: the prediction-side Type II majority, the exact complement of
035's Type-I stratum.

Everything here is EXACT. `argmin` over a uniformly random subset has a combinatorial
distribution, so no rung is subsampled and no rung carries Monte-Carlo error. The one
thing that is sampled is the noise floor, at 2e6 draws with its SE reported.

    python scripts/structure/depth_turnover_analysis.py controls
    python scripts/structure/depth_turnover_analysis.py distinct
    python scripts/structure/depth_turnover_analysis.py curves
    python scripts/structure/depth_turnover_analysis.py all
"""
from __future__ import annotations

import argparse
import json
import sys
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "cofold"))
sys.path.insert(0, str(REPO / "scripts" / "explorer"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

COORD_MAX = 2.6
DEDUPE_TOL = 0.05           # cypstruct.xengine._dedupe default, A
TAGS = ["dt2a", "dt2b", "dt2c", "dt2d"]
DEPTH_E = 20                # existing Modal arm
DEPTH_N = 20                # new Explorer arm
FLOOR_BLOCKS, FLOOR_BLOCK = 20, 100_000      # 2e6 draws, SE from block spread
BOOT = 10_000
RNG = np.random.default_rng(20260923)

RUNGS = [5, 10, 15, 20, 25, 28, 30, 33, 35, 38, 40]
PERLIG_RUNGS = [20, 30, 33, 40]

OUT_POSES = DATA_PROCESSED / "depth_turnover_poses.csv"
OUT_CONTROLS = DATA_PROCESSED / "depth_turnover_controls.json"
OUT_DISTINCT = DATA_PROCESSED / "depth_turnover_distinct.json"
OUT_CURVES = DATA_PROCESSED / "depth_turnover_curves.json"
OUT_PERLIG = DATA_PROCESSED / "depth_turnover_per_ligand.csv"


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

def stratum() -> list[str]:
    """The pre-registered rule, applied. Prediction-side; nothing from a crystal."""
    L = pd.read_csv(DATA_PROCESSED / "binding_mode_labels_cyp3a4.csv")
    return sorted(L[L.pred_fe_donor_median <= COORD_MAX].id.tolist())


def type_i_stratum() -> list[str]:
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
    """The four scored Explorer batches, concatenated with their tag kept (confound K6)."""
    frames = []
    for t in TAGS:
        f = DATA_PROCESSED / f"matched_depth_{t}_poses.csv"
        if not f.exists():
            continue
        d = pd.read_csv(f)
        d["batch"] = t
        frames.append(d)
    if not frames:
        raise SystemExit("no scored batches yet -- run boltz_depth.py score --tag dt2*")
    out = pd.concat(frames, ignore_index=True)
    out.to_csv(OUT_POSES, index=False)
    return out


def _pick(sub: pd.DataFrame) -> float:
    return float(sub.lddt_pli.iloc[int(np.argmin(sub.xeng.values))])


def summarise(df: pd.DataFrame, ligs: list[str]) -> dict:
    orc, sel, rnd = [], [], []
    for lid in ligs:
        s = df[df.ligand == lid]
        orc.append(s.lddt_pli.max())
        sel.append(_pick(s))
        rnd.append(s.lddt_pli.mean())
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
# controls
# ---------------------------------------------------------------------------

def controls() -> dict:
    from cypstruct import xengine as X

    pool = existing_pool()
    all87 = sorted(pool.ligand.unique())
    ligs = stratum()

    # N2 -- the shipped board reproduces, on all 87
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

    # N2b -- argmin(xeng) IS cypstruct.xengine.select(), to 1e-9
    shipped = X.select(pool)
    mine = np.array([_pick(pool[pool.ligand == lid]) for lid in all87])
    ship = shipped.set_index("ligand").loc[all87, "lddt_pli"].values
    n2b = {"max_abs_diff": float(np.max(np.abs(mine - ship))),
           "equal_to_1e-9": bool(np.max(np.abs(mine - ship)) < 1e-9), "n": len(all87)}

    # REF -- reference depth, the FINDING 011 precondition
    refset = X.load_reference(DATA_PROCESSED / "reference_set_cyp3a4.npz")
    depths = {lid: len(refset.get(lid, [])) for lid in ligs}
    ref = {"min": int(min(depths.values())), "median": float(np.median(list(depths.values()))),
           "below_4": [k for k, v in depths.items() if v < 4],
           "histogram": {str(k): int(v) for k, v in
                         pd.Series(list(depths.values())).value_counts().sort_index().items()}}

    # N5 -- the filter, firing, with counts, and the zero overlap asserted
    t1 = type_i_stratum()
    n5 = {"validation_ligands": int(pool.ligand.nunique()),
          "pass_pred_fe_donor_median_le_2.6": len(ligs),
          "excluded_as_finding_035_type_I_stratum": len(t1),
          "overlap_with_finding_035_stratum": sorted(set(ligs) & set(t1)),
          "stratum": ligs}
    assert not (set(ligs) & set(t1)), "the two strata overlap -- the rule is broken"

    base = summarise(pool[pool.ligand.isin(ligs)], ligs)
    rho73, sign73 = within_ligand_rho(pool[pool.ligand.isin(ligs)], ligs)
    g = pool[pool.ligand.isin(ligs)].groupby("ligand")
    div = {"unique_xeng_per_ligand": sorted(set(g.xeng.nunique().tolist())),
           "median_within_ligand_lddt_sd": round(float(g.lddt_pli.std().median()), 4),
           "ligands_with_sd_zero": int((g.lddt_pli.std() == 0).sum())}

    out = {"N2_shipped_column": n2, "N2b_argmin_identity": n2b,
           "N5_filters": n5, "REF_reference_depth": ref,
           "existing_pool_diversity_in_stratum": div,
           "baseline_depth20_stratum": {k: round(v, 4) for k, v in base.items()
                                        if not k.startswith("per_")},
           "baseline_rho": round(rho73, 4),
           "baseline_correct_sign_pct": round(100 * sign73, 2)}

    if any((DATA_PROCESSED / f"matched_depth_{t}_poses.csv").exists() for t in TAGS):
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
            "n": int(len(z)), "ligands": z.ligand.value_counts().to_dict(),
            "bisy_rmsd_range": [round(float(z.bisy_rmsd.min()), 2),
                                round(float(z.bisy_rmsd.max()), 2)] if len(z) else None,
            "fe_donor_range": [round(float(z.fe_donor_dist.min()), 2),
                               round(float(z.fe_donor_dist.max()), 2)] if len(z) else None,
            "all_offset_zero": bool((z.renumber_offset == 0).all()) if len(z) else None,
            "all_mapped": bool(z.mapped.all()) if len(z) else None,
        }
        out["N4_md5_distinct"] = {
            "rows": int(len(nw)), "distinct_md5": int(nw.md5.nunique()),
            "min_per_ligand_distinct_md5": int(nw.groupby("ligand").md5.nunique().min()),
            "ligands": int(nw.ligand.nunique()),
            "poses_per_ligand": sorted(set(nw.groupby("ligand").size().tolist())),
        }
        # A0 -- THE GATE. Reported before any curve.
        out["A0_conditioning_match"] = {
            "new_arm_pool_mean": round(float(nw.lddt_pli.mean()), 4),
            "existing_pool_mean_same_ligands": round(base["random"], 4),
            "delta": round(float(nw.lddt_pli.mean()) - base["random"], 4),
            "matched_within_0.02": bool(
                abs(float(nw.lddt_pli.mean()) - base["random"]) <= 0.02),
            "finding_035_delta_for_reference": -0.0019,
            "finding_034_openprotein_delta_for_reference": -0.057,
        }
        ex = pool[pool.ligand.isin(ligs)]
        out["A0_secondary_binding_character"] = {
            "new_frac_coordinated": round(float(nw.is_coordinated.mean()), 4),
            "existing_frac_coordinated": round(float(ex.is_coordinated.mean()), 4),
            "new_median_fe_donor": round(float(nw.fe_donor_dist.median()), 3),
            "existing_median_fe_donor": round(float(ex.fe_donor_dist.median()), 3),
        }
        # K6 -- per-batch means, so a job-level effect cannot hide inside the arm mean
        out["K6_per_batch"] = {
            t: {"ligands": int(d.ligand.nunique()), "poses": int(len(d)),
                "pool_mean": round(float(d.lddt_pli.mean()), 4),
                "existing_mean_same_ligands": round(
                    float(pool[pool.ligand.isin(d.ligand.unique())].lddt_pli.mean()), 4)}
            for t, d in nw.groupby("batch")}
    OUT_CONTROLS.write_text(json.dumps(out, indent=2, default=str))
    return out


# ---------------------------------------------------------------------------
# N4 -- DISTINCT poses, by geometry in the heme frame, not by job
# ---------------------------------------------------------------------------

def distinct() -> dict:
    from cypstruct import pose as P
    from cypstruct import xengine as X
    from boltz_depth import POSE_DIR

    nw = new_pool()
    keep: dict[str, list[str]] = {}
    stats = {}
    for lid, grp in nw.groupby("ligand"):
        tag = grp.batch.iloc[0]
        flat = POSE_DIR / f"{tag}_flat"
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
        kept_xyz: list[np.ndarray] = []
        kept_names: list[str] = []
        for c, n in zip(coords, names):
            if not any(len(q) == len(c) and float(np.abs(q - c).max()) < DEDUPE_TOL
                       for q in kept_xyz):
                kept_xyz.append(c)
                kept_names.append(n)
        keep[lid] = kept_names
        stats[lid] = {"collected": len(names), "distinct": len(kept_names)}
    d = [v["distinct"] for v in stats.values()]
    out = {"tol_angstrom": DEDUPE_TOL, "per_ligand": stats,
           "total_collected": int(sum(v["collected"] for v in stats.values())),
           "total_distinct": int(sum(d)),
           "min_distinct": int(min(d)) if d else 0,
           "mean_distinct": round(float(np.mean(d)), 3) if d else 0,
           "ligands_with_20_distinct": int(sum(x == 20 for x in d)),
           "ligands_with_ge_16_distinct": int(sum(x >= 16 for x in d)),
           "keep": keep}
    OUT_DISTINCT.write_text(json.dumps(out, indent=2))
    return {k: v for k, v in out.items() if k != "keep"}


# ---------------------------------------------------------------------------
# the exact curves
# ---------------------------------------------------------------------------

def _exact_union_selected(u: np.ndarray, d: int) -> float:
    """E[ lddt of argmin(xeng) ] over a uniform d-subset of this ligand's M poses.

    Pose of rank j (1-based, ascending xeng) wins iff it is drawn and none of the j-1
    strictly better ones are: P = C(M-j, d-1) / C(M, d).
    """
    u = u[np.argsort(u[:, 0], kind="stable")]
    M = len(u)
    den = comb(M, d)
    tot = 0.0
    for j in range(1, M + 1):
        rest = M - j
        if rest >= d - 1:
            tot += comb(rest, d - 1) / den * u[j - 1, 1]
    return tot


def _exact_union_oracle(u: np.ndarray, d: int) -> float:
    y = np.sort(u[:, 1])
    M = len(y)
    den = comb(M, d)
    return float(sum(comb(j - 1, d - 1) / den * y[j - 1]
                     for j in range(d, M + 1)))


def _exact_augment_selected(e: np.ndarray, n: np.ndarray, k: int) -> float:
    """All of `e` kept, k of `n` added. Exact."""
    i = int(np.argmin(e[:, 0]))
    ex, ey = e[i, 0], e[i, 1]
    if k == 0:
        return float(ey)
    w = n[n[:, 0] < ex]
    w = w[np.argsort(w[:, 0], kind="stable")]
    m, M = len(w), len(n)
    den = comb(M, k)
    tot = 0.0
    for j in range(1, m + 1):
        rest = M - j
        if rest >= k - 1:
            tot += comb(rest, k - 1) / den * w[j - 1, 1]
    pnone = comb(M - m, k) / den if M - m >= k else 0.0
    return float(tot + pnone * ey)


def _exact_augment_oracle(e: np.ndarray, n: np.ndarray, k: int) -> float:
    emax = float(e[:, 1].max())
    if k == 0:
        return emax
    v = n[n[:, 1] > emax][:, 1]
    v = np.sort(v)[::-1]                     # descending
    m, M = len(v), len(n)
    den = comb(M, k)
    tot = 0.0
    for j in range(1, m + 1):
        rest = M - j
        if rest >= k - 1:
            tot += comb(rest, k - 1) / den * v[j - 1]
    pnone = comb(M - m, k) / den if M - m >= k else 0.0
    return float(tot + pnone * emax)


def noise_floor(mat: np.ndarray, blocks: int = FLOOR_BLOCKS,
                block: int = FLOOR_BLOCK) -> dict:
    """p95/p99 of a random feature's gain, on THIS population, with the SE of the p95.

    `mat` is (n_ligands, depth) of LDDT-PLI. A random feature's argmin over a ligand's
    poses is a uniform draw from them (FINDING 039 s1a), so the gain is
    mean(drawn) - mean(pool means).
    """
    nl, dep = mat.shape
    base = mat.mean(axis=1).mean()
    rng = np.random.default_rng(20260923)
    p95s, p99s, allg = [], [], []
    for _ in range(blocks):
        idx = rng.integers(0, dep, size=(block, nl))
        g = mat[np.arange(nl)[None, :], idx].mean(axis=1) - base
        p95s.append(np.percentile(g, 95))
        p99s.append(np.percentile(g, 99))
        allg.append(g)
    allg = np.concatenate(allg)
    return {"n_ligands": int(nl), "depth": int(dep), "draws": int(blocks * block),
            "p95": round(float(np.percentile(allg, 95)), 5),
            "p99": round(float(np.percentile(allg, 99)), 5),
            "se_of_p95": round(float(np.std(p95s, ddof=1) / np.sqrt(blocks)), 6),
            "sd_of_gain": round(float(allg.std(ddof=1)), 5)}


def curves() -> dict:
    from scipy.stats import wilcoxon

    pool = existing_pool()
    nw = new_pool()
    ligs_all = stratum()
    dist = json.loads(OUT_DISTINCT.read_text())
    keep = dist["keep"]

    # PRE-REGISTERED yield rule: primary curve on ligands with exactly 20 distinct new
    # poses; a ragged depth axis is not a depth axis.
    ligs = [lid for lid in ligs_all if len(keep.get(lid, [])) == DEPTH_N]
    dropped = {lid: len(keep.get(lid, [])) for lid in ligs_all
               if len(keep.get(lid, [])) != DEPTH_N}

    keep_pairs = {(lid, s) for lid, v in keep.items() for s in v}
    nsel = nw[[(r.ligand, r.sample) in keep_pairs for r in nw.itertuples()]]

    E = {lid: pool[pool.ligand == lid][["xeng", "lddt_pli"]].values for lid in ligs}
    N = {lid: nsel[nsel.ligand == lid][["xeng", "lddt_pli"]].values for lid in ligs}
    for lid in ligs:
        assert len(E[lid]) == DEPTH_E and len(N[lid]) == DEPTH_N, (
            lid, len(E[lid]), len(N[lid]))
    U = {lid: np.vstack([E[lid], N[lid]]) for lid in ligs}
    M = DEPTH_E + DEPTH_N

    # --- the bound on what depth can buy at all --------------------------------
    movable = []
    for lid in ligs:
        inc = E[lid][:, 0].min()
        movable.append(bool((N[lid][:, 0] < inc).any()))
    movable = np.array(movable)
    n_better = np.array([int((N[lid][:, 0] < E[lid][:, 0].min()).sum()) for lid in ligs])

    # --- PRIMARY: the exact union curve, every rung ------------------------------
    sel_u = {d: float(np.mean([_exact_union_selected(U[lid], d) for lid in ligs]))
             for d in range(1, M + 1)}
    orc_u = {d: float(np.mean([_exact_union_oracle(U[lid], d) for lid in ligs]))
             for d in range(1, M + 1)}
    rnd_u = float(np.mean([U[lid][:, 1].mean() for lid in ligs]))

    seq = np.array([sel_u[d] for d in range(1, M + 1)])
    dstar = int(np.argmax(seq)) + 1
    per_lig_star = np.array([_exact_union_selected(U[lid], dstar) for lid in ligs])
    per_lig_40 = np.array([_exact_union_selected(U[lid], M) for lid in ligs])
    delta = per_lig_star - per_lig_40
    d_peak = float(delta.mean())

    # --- acceptance T1..T4 -------------------------------------------------------
    t1 = bool(dstar < M)
    t2 = bool(d_peak >= 0.0050)
    S = delta.sum()
    loo = (S - delta) / (len(delta) - 1)
    t3 = bool(np.all(loo > 0)) if t1 else False
    bs = np.array([RNG.choice(delta, len(delta), replace=True).mean()
                   for _ in range(BOOT)])
    ci = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
    try:
        wp = float(wilcoxon(delta).pvalue) if np.any(np.abs(delta) > 1e-12) else float("nan")
    except ValueError:
        wp = float("nan")
    t4 = bool((ci[0] > 0 or ci[1] < 0) and np.isfinite(wp) and wp < 0.05)

    if not t1:
        verdict = "NO TURNOVER (monotone)"
    elif not t3:
        verdict = "TURNOVER IS ONE LIGAND"
    elif t1 and t2 and t3 and t4:
        verdict = "TURNOVER REAL"
    else:
        verdict = "TURNOVER UNDECIDED"

    # --- SECONDARY: the augment (purchase) curve, 20 -> 40 -----------------------
    sel_a = {20 + k: float(np.mean([_exact_augment_selected(E[lid], N[lid], k)
                                    for lid in ligs])) for k in range(0, DEPTH_N + 1)}
    orc_a = {20 + k: float(np.mean([_exact_augment_oracle(E[lid], N[lid], k)
                                    for lid in ligs])) for k in range(0, DEPTH_N + 1)}
    rnd_a = {20 + k: float(np.mean([(E[lid][:, 1].sum() + k * N[lid][:, 1].mean())
                                    / (20 + k) for lid in ligs]))
             for k in range(0, DEPTH_N + 1)}
    aseq = np.array([sel_a[20 + k] for k in range(0, DEPTH_N + 1)])
    a_star = 20 + int(np.argmax(aseq))

    pl_a0 = np.array([_exact_augment_selected(E[lid], N[lid], 0) for lid in ligs])
    pl_a20 = np.array([_exact_augment_selected(E[lid], N[lid], DEPTH_N) for lid in ligs])
    pl_o0 = np.array([_exact_augment_oracle(E[lid], N[lid], 0) for lid in ligs])
    pl_o20 = np.array([_exact_augment_oracle(E[lid], N[lid], DEPTH_N) for lid in ligs])
    dsel = pl_a20 - pl_a0
    dorc = pl_o20 - pl_o0

    def boot(v):
        m = np.array([RNG.choice(v, len(v), replace=True).mean() for _ in range(BOOT)])
        return [round(float(np.percentile(m, 2.5)), 4),
                round(float(np.percentile(m, 97.5)), 4)]

    try:
        wp_a = float(wilcoxon(dsel).pvalue) if np.any(np.abs(dsel) > 1e-12) else float("nan")
    except ValueError:
        wp_a = float("nan")

    # --- rates per doubling, on the unbiased union curve -------------------------
    oct_pairs = [(2.5, 5), (5, 10), (10, 20), (20, 40)]
    rates = {}
    for a, b in oct_pairs:
        ai, bi = int(round(a)), int(round(b))
        if ai < 1:
            continue
        rates[f"{ai}->{bi}"] = {
            "selected": round(sel_u[bi] - sel_u[ai], 4),
            "oracle": round(orc_u[bi] - orc_u[ai], 4)}

    # where the marginal selected rate per doubling falls below FINDING 004's +0.0125
    rec_depth = None
    for a, b in [(5, 10), (10, 20), (20, 40)]:
        if sel_u[b] - sel_u[a] < 0.0125:
            rec_depth = a
            break

    # --- nulls, on THIS population, 2e6 draws ------------------------------------
    mat20 = np.array([E[lid][:, 1] for lid in ligs])
    mat40 = np.array([U[lid][:, 1] for lid in ligs])
    floor20 = noise_floor(mat20)
    floor40 = noise_floor(mat40)

    aug_df = pd.concat([pool[pool.ligand.isin(ligs)], nsel[nsel.ligand.isin(ligs)]],
                       ignore_index=True)
    rho_a, sign_a = within_ligand_rho(aug_df, ligs)
    rho_b, sign_b = within_ligand_rho(pool[pool.ligand.isin(ligs)], ligs)
    gain40 = sel_u[M] - rnd_u
    gain20 = float(np.mean([_exact_augment_selected(E[lid], N[lid], 0)
                            for lid in ligs])) - float(np.mean([E[lid][:, 1].mean()
                                                                for lid in ligs]))

    out = {
        "stratum_n": len(ligs), "stratum_n_planned": len(ligs_all),
        "ligands": ligs, "dropped_for_yield": dropped,
        "movable": {
            "n_with_any_new_pose_beating_incumbent_xeng": int(movable.sum()),
            "fraction": round(float(movable.mean()), 4),
            "n_cannot_move": int((~movable).sum()),
            "median_new_poses_beating_incumbent": float(np.median(n_better)),
            "finding_039_comparison": "10 of 14 could not move (28.6% movable)",
        },
        "union_curve_exact": {
            str(d): {"oracle": round(orc_u[d], 5), "selected": round(sel_u[d], 5),
                     "random": round(rnd_u, 5)}
            for d in RUNGS},
        "union_curve_full_selected": {str(d): round(sel_u[d], 5)
                                      for d in range(1, M + 1)},
        "union_curve_full_oracle": {str(d): round(orc_u[d], 5) for d in range(1, M + 1)},
        "union_selected_argmax_depth": dstar,
        "union_selected_monotone": bool(np.all(np.diff(seq) >= -1e-12)),
        "union_oracle_monotone": bool(
            np.all(np.diff([orc_u[d] for d in range(1, M + 1)]) >= -1e-12)),
        "turnover": {
            "d_star": dstar,
            "selected_at_d_star": round(float(sel_u[dstar]), 5),
            "selected_at_40": round(float(sel_u[M]), 5),
            "delta_peak": round(d_peak, 5),
            "delta_peak_CI95": [round(ci[0], 5), round(ci[1], 5)],
            "wilcoxon_p": round(wp, 6) if np.isfinite(wp) else None,
            "n_ligands_moving": int(np.sum(np.abs(delta) > 1e-12)),
            "loo_min_delta": round(float(loo.min()), 5),
            "loo_max_delta": round(float(loo.max()), 5),
            "loo_all_positive": t3,
            "T1_interior_peak": t1, "T2_magnitude_ge_0.0050": t2,
            "T3_survives_leave_one_ligand_out": t3, "T4_CI_excludes_0_and_p<0.05": t4,
            "VERDICT": verdict,
        },
        "augment_curve_exact": {
            str(20 + k): {"oracle": round(orc_a[20 + k], 5),
                          "selected": round(sel_a[20 + k], 5),
                          "random": round(rnd_a[20 + k], 5)}
            for k in range(0, DEPTH_N + 1)},
        "augment_selected_argmax_pool": a_star,
        "augment_primary_20_to_40": {
            "delta_selected": round(float(dsel.mean()), 5),
            "delta_selected_CI95": boot(dsel),
            "delta_oracle": round(float(dorc.mean()), 5),
            "delta_oracle_CI95": boot(dorc),
            "wilcoxon_p": round(wp_a, 6) if np.isfinite(wp_a) else None,
            "ligands_better": int((dsel > 1e-12).sum()),
            "ligands_worse": int((dsel < -1e-12).sum()),
            "ligands_tied": int((np.abs(dsel) <= 1e-12).sum()),
            "finding_035_delta_selected": 0.0116,
            "finding_035_ligands_tied": "10 of 14",
        },
        "rate_per_doubling_union": rates,
        "depth_at_which_selected_rate_falls_below_0.0125": rec_depth,
        "noise_floor_depth20_this_population": floor20,
        "noise_floor_depth40_this_population": floor40,
        "gain_vs_floor": {
            "gain_depth20": round(gain20, 5), "floor20_p95": floor20["p95"],
            "clears20": bool(gain20 > floor20["p95"]),
            "gain_depth40": round(gain40, 5), "floor40_p95": floor40["p95"],
            "clears40": bool(gain40 > floor40["p95"]),
        },
        "within_ligand_rho": {
            "depth20": round(rho_b, 4), "depth20_correct_sign_pct": round(100 * sign_b, 2),
            "depth40": round(rho_a, 4), "depth40_correct_sign_pct": round(100 * sign_a, 2),
        },
    }

    rows = {"ligand": ligs,
            "movable": movable,
            "n_new_beating_incumbent": n_better,
            "delta_peak_minus_40": delta,
            "aug_sel_20": pl_a0, "aug_sel_40": pl_a20, "aug_d_sel": dsel,
            "aug_orc_20": pl_o0, "aug_orc_40": pl_o20, "aug_d_orc": dorc}
    for d in sorted(set(PERLIG_RUNGS + [dstar])):
        rows[f"union_sel_d{d}"] = [_exact_union_selected(U[lid], d) for lid in ligs]
        rows[f"union_orc_d{d}"] = [_exact_union_oracle(U[lid], d) for lid in ligs]
    pd.DataFrame(rows).to_csv(OUT_PERLIG, index=False)

    OUT_CURVES.write_text(json.dumps(out, indent=2, default=str))
    return {k: v for k, v in out.items()
            if k not in ("ligands", "union_curve_full_selected",
                         "union_curve_full_oracle", "augment_curve_exact")}


# ---------------------------------------------------------------------------
# X1 -- the closed-form functions reproduce FINDING 039 s3b, exactly
# ---------------------------------------------------------------------------

def x1() -> dict:
    """Replay FINDING 039 s3b's exact curve with THIS file's estimators.

    The closed forms below are new code on a new stratum, so before any of their numbers
    are believed they are made to reproduce, to four decimals, every rung FINDING 039
    published on the 14-ligand Type-I stratum -- and FINDING 035's two oracle endpoints.
    A curve that agrees with a Monte-Carlo estimate to 0.0004 and with a published exact
    curve to 0.0000 is not a new implementation of a different quantity.
    """
    pool = existing_pool()
    S = type_i_stratum()
    new = pd.read_csv(DATA_PROCESSED / "matched_depth_poses.csv")
    dist = json.loads((DATA_PROCESSED / "matched_depth_distinct.json").read_text())
    keep = {(l, s) for l, v in dist["keep"].items() for s in v}
    new = new[[(r.ligand, r.sample) in keep for r in new.itertuples()]]
    E = {l: pool[pool.ligand == l][["xeng", "lddt_pli"]].values for l in S}
    N = {l: new[new.ligand == l][["xeng", "lddt_pli"]].values for l in S}
    U = {l: np.vstack([E[l], N[l]]) for l in S}
    aug = lambda k: float(np.mean([_exact_augment_selected(E[l], N[l], k) for l in S]))
    uni = lambda d: float(np.mean([_exact_union_selected(U[l], d) for l in S]))
    pub = {20: 0.5270, 24: 0.5349, 28: 0.5400, 30: 0.5415, 33: 0.5424,
           36: 0.5418, 38: 0.5406, 40: 0.5386}
    got = {str(k): round(aug(k - 20), 4) for k in pub}
    seq = [uni(d) for d in range(1, 41)]
    out = {
        "finding_039_augment_rungs_published": {str(k): v for k, v in pub.items()},
        "reproduced_here": got,
        "max_abs_deviation": round(max(abs(got[str(k)] - v) for k, v in pub.items()), 5),
        "union_argmax_depth": 1 + int(np.argmax(seq)),
        "union_peak": round(float(max(seq)), 4),
        "union_depth40": round(float(seq[39]), 4),
        "finding_039_union": {"argmax_depth": 34, "peak": 0.5410, "depth40": 0.5386},
        "oracle_augment_k0": round(float(np.mean(
            [_exact_augment_oracle(E[l], N[l], 0) for l in S])), 4),
        "oracle_augment_k20": round(float(np.mean(
            [_exact_augment_oracle(E[l], N[l], 20) for l in S])), 4),
        "finding_035_oracle": {"pool20": 0.6499, "pool40": 0.6595},
    }
    out["PASSES"] = bool(out["max_abs_deviation"] < 5e-4
                         and out["union_argmax_depth"] == 34)
    (DATA_PROCESSED / "depth_turnover_x1.json").write_text(json.dumps(out, indent=2))
    return out


# ---------------------------------------------------------------------------
# decomposition -- per ligand, and the combined 87-ligand curve
# ---------------------------------------------------------------------------

def decomp() -> dict:
    """Who moves, by how much, and what happens if the biggest movers are removed.

    FINDING 039's whole point was that a 14-ligand curve was ONE ligand. The same
    leverage arithmetic is run here, at every rung that matters, so the answer cannot
    hide the same way twice. The combined 87-ligand curve at the bottom is NOT
    pre-registered -- it is a descriptive aggregation of this stratum with FINDING 035's,
    labelled as such.
    """
    pool = existing_pool()
    nw = new_pool()
    dist = json.loads(OUT_DISTINCT.read_text())
    keep_pairs = {(lid, s) for lid, v in dist["keep"].items() for s in v}
    nsel = nw[[(r.ligand, r.sample) in keep_pairs for r in nw.itertuples()]]
    ligs = stratum()
    E = {l: pool[pool.ligand == l][["xeng", "lddt_pli"]].values for l in ligs}
    N = {l: nsel[nsel.ligand == l][["xeng", "lddt_pli"]].values for l in ligs}
    U = {l: np.vstack([E[l], N[l]]) for l in ligs}

    def sel(d, L, UU):
        return float(np.mean([_exact_union_selected(UU[l], d) for l in L]))

    # per-ligand 30 -> 40, the exact comparison FINDING 039 s3b made
    rows = []
    for l in ligs:
        s30 = _exact_union_selected(U[l], 30)
        s40 = _exact_union_selected(U[l], 40)
        a20 = _exact_augment_selected(E[l], N[l], 0)
        a40 = _exact_augment_selected(E[l], N[l], 20)
        rows.append({"ligand": l, "union_sel_30": s30, "union_sel_40": s40,
                     "d_30_minus_40": s30 - s40,
                     "aug_sel_20": a20, "aug_sel_40": a40, "aug_delta": a40 - a20,
                     "n_new_beating_incumbent": int((N[l][:, 0] < E[l][:, 0].min()).sum())})
    df = pd.DataFrame(rows).sort_values("d_30_minus_40", ascending=False)
    df.to_csv(DATA_PROCESSED / "depth_turnover_decomp.csv", index=False)

    # leverage: drop the k ligands with the largest |augment delta|
    aug = df.set_index("ligand").aug_delta
    order = aug.abs().sort_values(ascending=False).index.tolist()
    lev = {}
    for k in (0, 1, 2, 3, 5):
        rem = [l for l in ligs if l not in order[:k]]
        lev[f"drop_top_{k}"] = {
            "n": len(rem),
            "delta_selected_20_to_40": round(float(aug.loc[rem].mean()), 5),
            "union_sel_20": round(sel(20, rem, U), 5),
            "union_sel_40": round(sel(40, rem, U), 5),
            "union_delta": round(sel(40, rem, U) - sel(20, rem, U), 5),
            "dropped": order[:k]}

    # can ANY single ligand invert the 30->40 comparison, as 08J did at n=14?
    tot = float(df.d_30_minus_40.sum())
    worst = df.iloc[0]
    out = {
        "top5_pushing_toward_a_turnover_30_minus_40": df.head(5).round(5).to_dict("records"),
        "top5_pushing_against_30_minus_40": df.tail(5).round(5).to_dict("records"),
        "n_ligands_with_nonzero_30_minus_40": int((df.d_30_minus_40.abs() > 1e-12).sum()),
        "mean_30_minus_40": round(tot / len(ligs), 5),
        "largest_single_ligand_30_minus_40": round(float(worst.d_30_minus_40), 5),
        "finding_039_largest_was_08J": 0.1721,
        "leverage_on_the_20_to_40_purchase": lev,
    }

    # combined 87 = this stratum + FINDING 035's, both at matched depth 40
    old = pd.read_csv(DATA_PROCESSED / "matched_depth_poses.csv")
    odist = json.loads((DATA_PROCESSED / "matched_depth_distinct.json").read_text())
    okeep = {(l, s) for l, v in odist["keep"].items() for s in v}
    old = old[[(r.ligand, r.sample) in okeep for r in old.itertuples()]]
    ligs87 = ligs + type_i_stratum()
    U87 = dict(U)
    for l in type_i_stratum():
        U87[l] = np.vstack([pool[pool.ligand == l][["xeng", "lddt_pli"]].values,
                            old[old.ligand == l][["xeng", "lddt_pli"]].values])
    seq87 = [sel(d, ligs87, U87) for d in range(1, 41)]
    out["combined_87_union_selected"] = {
        str(d): round(seq87[d - 1], 5) for d in RUNGS}
    out["combined_87_oracle"] = {
        str(d): round(float(np.mean([_exact_union_oracle(U87[l], d) for l in ligs87])), 5)
        for d in RUNGS}
    out["combined_87_argmax_depth"] = 1 + int(np.argmax(seq87))
    out["combined_87_monotone"] = bool(np.all(np.diff(seq87) >= -1e-12))
    out["combined_87_note"] = ("NOT pre-registered -- a descriptive aggregation of this "
                               "stratum with FINDING 035's, at matched depth 40 on all 87")
    (DATA_PROCESSED / "depth_turnover_decomp.json").write_text(
        json.dumps(out, indent=2, default=str))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["x1", "controls", "distinct", "curves", "decomp", "all"])
    a = ap.parse_args()
    if a.cmd in ("x1", "all"):
        print(json.dumps(x1(), indent=2, default=str))
    if a.cmd in ("controls", "all"):
        print(json.dumps(controls(), indent=2, default=str))
    if a.cmd in ("distinct", "all"):
        print(json.dumps(distinct(), indent=2, default=str))
    if a.cmd in ("curves", "all"):
        print(json.dumps(curves(), indent=2, default=str))
    if a.cmd in ("decomp", "all"):
        print(json.dumps(decomp(), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
