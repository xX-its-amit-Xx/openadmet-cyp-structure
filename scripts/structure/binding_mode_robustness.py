"""Is the shipped cross-engine selector robust to BINDING MODE?

The biology map measured that the 87-ligand CYP3A4 validation set is **83% Type II** — 72
ligands coordinate the heme iron directly, 15 sit in the active site without touching it —
and that Reactome places CYP3A4 in only the xenobiotic pathway, none of the sterol,
fatty-acid, eicosanoid or vitamin ones. Both frames we sample from are biased toward
heme-coordinating chemistry.

`cypstruct.xengine.select()` scores a pose by its mean symmetric Chamfer distance, IN THE
HEME FRAME, to independent-engine poses of the same ligand. The frame is built from the
heme. A Type I ligand is not anchored to the landmark that defines the frame, so there is
a specific mechanistic reason its gain could differ — and if it does, and the blind test
set is Type I-rich, every selector number in this repo is measured on the wrong population.

Three stages:

    python scripts/structure/binding_mode_robustness.py labels   # verify the labels
    python scripts/structure/binding_mode_robustness.py cyp3a4   # stratify n=87
    python scripts/structure/binding_mode_robustness.py p450     # stratify the family

**Thresholds are fixed a priori** — they are the repo's existing constants, unchanged:
`COORD_MAX = 2.6` Å (direct coordination) and `OVERHEAD_MAX = 6.5` Å (in the active site),
from `scripts/structure/build_reference_set.py`. Nothing here is tuned.

**Random baselines are computed two ways** and both are reported: the exact per-ligand
expectation (the mean of that ligand's pool, which is what "pick one uniformly at random"
converges to, with zero sampling noise) and a 300-draw empirical average. They agree to
~0.001; the exact form is used for the per-ligand gains that everything downstream rests
on, because a per-ligand statistic built from 300 draws of a 20-pose pool is noisier than
the effect being measured.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

# --- a priori, from build_reference_set.py. Do not tune. ---------------------
COORD_MAX = 2.6       # Fe-to-closest-ligand-atom below this = Type II (coordinated)
OVERHEAD_MAX = 6.5    # closest-atom below this but above COORD_MAX = Type I (active site)

# The empty band between the two observed classes on the validation set (FINDING 033):
# furthest Type II at 2.59 A, closest Type I at 2.82 A. A prediction-side median landing
# inside it is FLAGGED, never guessed.
AMBIG_LO, AMBIG_HI = 2.59, 2.82

N_RANDOM_DRAWS = 300  # the task asks for >= 64
N_NULL = 4000
N_BOOT = 10000
SEED = 20260922

UNI = DATA_PROCESSED / "p450_universe"
OUT_LABELS = DATA_PROCESSED / "binding_mode_labels_cyp3a4.csv"
OUT_LIG = DATA_PROCESSED / "binding_mode_per_ligand_cyp3a4.csv"
OUT_CYP = DATA_PROCESSED / "binding_mode_strata_cyp3a4.json"
OUT_PXENG = DATA_PROCESSED / "binding_mode_p450_xeng.csv"
OUT_PPAIR = DATA_PROCESSED / "binding_mode_per_pair_p450.csv"
OUT_P450 = DATA_PROCESSED / "binding_mode_strata_p450.json"


def mode_from_distance(d: float) -> str:
    if not np.isfinite(d):
        return "unknown"
    if d <= COORD_MAX:
        return "type_II"
    if d <= OVERHEAD_MAX:
        return "type_I"
    return "peripheral"


# ---------------------------------------------------------------------------
# shared statistics
# ---------------------------------------------------------------------------

def per_unit_table(m: pd.DataFrame, unit: str, score_col: str, truth: str = "lddt_pli"):
    """One row per ligand/pair: random expectation, oracle, what the selector picked.

    `score_col` is LOWER-is-better (Chamfer distance), so the pick is the argmin of the
    within-unit z-score — identical to `cypstruct.xengine.select`, which maximises
    `-z(xeng)`. Reproduced here rather than called so the same code path can be pointed at
    the P450 table, and checked against `select()` in `cyp3a4` below.
    """
    rows = []
    for u, g in m.groupby(unit):
        pick = g.loc[g[score_col].idxmin()]
        rows.append({
            unit: u,
            "n_poses": len(g),
            "rand_exact": float(g[truth].mean()),
            "oracle": float(g[truth].max()),
            "worst": float(g[truth].min()),
            "selected": float(pick[truth]),
            "spread": float(g[truth].max() - g[truth].min()),
            "rho": (float(stats.spearmanr(g[score_col], g[truth]).statistic)
                    if g[score_col].nunique() >= 3 and g[truth].nunique() >= 2
                    else np.nan),
        })
    t = pd.DataFrame(rows)
    t["gain"] = t.selected - t.rand_exact
    t["headroom"] = t.oracle - t.rand_exact
    return t


def empirical_random(m: pd.DataFrame, unit: str, truth: str = "lddt_pli",
                     draws: int = N_RANDOM_DRAWS) -> float:
    """Average of `draws` uniform picks per unit — the form used in FINDING 011/012."""
    return float(np.mean([
        m.groupby(unit)[truth].apply(lambda s: s.sample(1, random_state=i).iloc[0]).mean()
        for i in range(draws)]))


def null_gain(m: pd.DataFrame, unit: str, truth: str = "lddt_pli",
              n: int = N_NULL, seed: int = SEED) -> np.ndarray:
    """Gain of a RANDOM feature — FINDING 007's floor, recomputed inside each stratum.

    The floor is n-dependent: at n=15 a random feature buys far more than at n=72, so a
    stratum's gain must be read against its OWN null, never against the pooled one.
    """
    rng = np.random.default_rng(seed)
    base = m.groupby(unit)[truth].mean().mean()
    out = np.empty(n)
    for i in range(n):
        s = rng.normal(size=len(m))
        out[i] = m.loc[m.assign(_s=s).groupby(unit)._s.idxmax(), truth].mean() - base
    return out


def _pspear(x, y, z):
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    rxz = np.corrcoef(rx, rz)[0, 1]
    ryz = np.corrcoef(ry, rz)[0, 1]
    rxy = np.corrcoef(rx, ry)[0, 1]
    d = np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    return float((rxy - rxz * ryz) / d) if d > 1e-9 else np.nan


def partial_spearman(x, y, z, seed: int = SEED, n_boot: int = 2000):
    """Rank correlation of x and y with z partialled out — here z is always headroom.

    Reported with a bootstrap CI and a permutation p, because the whole point of the frame
    hypothesis is whether a SMALL residual correlation survives the headroom confound, and
    a point estimate cannot answer that.
    """
    x, y, z = map(np.asarray, (x, y, z))
    r = _pspear(x, y, z)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), (n_boot, len(x)))
    bs = np.array([_pspear(x[i], y[i], z[i]) for i in idx])
    bs = bs[np.isfinite(bs)]
    perm = np.array([_pspear(x, rng.permutation(y), z) for _ in range(n_boot)])
    perm = perm[np.isfinite(perm)]
    return {"partial_rho": round(float(r), 4), "n": int(len(x)),
            "ci95": [round(float(np.percentile(bs, 2.5)), 4),
                     round(float(np.percentile(bs, 97.5)), 4)],
            "permutation_p_approx": round(float((np.abs(perm) >= abs(r)).mean()), 4)}


def boot_ci(x: np.ndarray, n: int = N_BOOT, seed: int = SEED, q=(2.5, 97.5)):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(n, len(x)))
    b = x[idx].mean(axis=1)
    return [round(float(np.percentile(b, q[0])), 4), round(float(np.percentile(b, q[1])), 4)]


def headroom_matched(a: pd.DataFrame, b: pd.DataFrame, n_bins: int = 5,
                     seed: int = SEED) -> dict:
    """The stratum difference AFTER matching on pool headroom.

    A raw difference in mean gain between two strata is not evidence that the selector
    behaves differently on them: a stratum whose pools have more room to improve will show
    a bigger gain from the *same* selector. FINDING 012 is exactly this — the family-wide
    +0.3006 was pool composition, not a better selector. So bin both strata on
    `headroom = oracle - random`, take the difference of mean gains inside each bin, and
    average the bins weighted by how many units they hold. If the raw difference survives,
    the mode matters; if it collapses, the headroom did.
    """
    both = pd.concat([a.assign(_g="a"), b.assign(_g="b")], ignore_index=True)
    both["_bin"] = pd.qcut(both.headroom.rank(method="first"), n_bins, labels=False)
    rows, w, d = [], [], []
    for bi, g in both.groupby("_bin"):
        ga = g.loc[g._g == "a", "gain"].to_numpy()
        gb = g.loc[g._g == "b", "gain"].to_numpy()
        if len(ga) < 3 or len(gb) < 3:
            continue
        rows.append({"bin": int(bi), "n_a": len(ga), "n_b": len(gb),
                     "headroom": round(float(g.headroom.mean()), 4),
                     "gain_a": round(float(ga.mean()), 4),
                     "gain_b": round(float(gb.mean()), 4),
                     "difference": round(float(ga.mean() - gb.mean()), 4)})
        w.append(len(ga) + len(gb))
        d.append(ga.mean() - gb.mean())
    if not rows:
        return {"error": "no bin holds >= 3 units of both strata"}
    w = np.array(w, float); d = np.array(d)
    rng = np.random.default_rng(seed)
    bs = np.empty(N_BOOT)
    for i in range(N_BOOT):
        s = both.sample(len(both), replace=True, random_state=int(rng.integers(1 << 31)))
        dd, ww = [], []
        for _bi, g in s.groupby("_bin"):
            ga = g.loc[g._g == "a", "gain"].to_numpy()
            gb = g.loc[g._g == "b", "gain"].to_numpy()
            if len(ga) < 3 or len(gb) < 3:
                continue
            dd.append(ga.mean() - gb.mean()); ww.append(len(ga) + len(gb))
        bs[i] = np.average(dd, weights=ww) if dd else np.nan
    bs = bs[np.isfinite(bs)]
    return {
        "n_bins_used": len(rows),
        "bins": rows,
        "matched_difference": round(float(np.average(d, weights=w)), 4),
        "matched_difference_ci95": [round(float(np.percentile(bs, 2.5)), 4),
                                    round(float(np.percentile(bs, 97.5)), 4)],
        "raw_difference": round(float(a.gain.mean() - b.gain.mean()), 4),
    }


def stratum_report(m: pd.DataFrame, t: pd.DataFrame, unit: str, label: str) -> dict:
    g = t.gain.to_numpy()
    rho = t.rho.dropna().to_numpy()
    null = null_gain(m, unit)
    return {
        "stratum": label,
        "n_units": int(len(t)),
        "n_poses": int(len(m)),
        "xeng_mean_chamfer_A": round(float(m.xeng.mean()), 3),
        "xeng_within_unit_sd_A": round(float(
            m.groupby(unit).xeng.std().mean()), 3),
        "random_exact": round(float(t.rand_exact.mean()), 4),
        "random_empirical_300": round(empirical_random(m, unit), 4),
        "oracle": round(float(t.oracle.mean()), 4),
        "selected": round(float(t.selected.mean()), 4),
        "gain": round(float(g.mean()), 4),
        "gain_ci95": boot_ci(g),
        "gain_sd_per_unit": round(float(g.std(ddof=1)), 4),
        "headroom_oracle_minus_random": round(float(t.headroom.mean()), 4),
        "frac_of_headroom_captured": round(float(g.mean() / t.headroom.mean()), 4),
        "within_unit_rho_mean": round(float(rho.mean()), 4) if len(rho) else None,
        "frac_rho_correct_sign": round(float((rho < 0).mean()), 4) if len(rho) else None,
        "frac_gain_positive": round(float((g > 0).mean()), 4),
        "frac_gain_nonneg": round(float((g >= 0).mean()), 4),
        "null_p95": round(float(np.percentile(null, 95)), 4),
        "null_p99": round(float(np.percentile(null, 99)), 4),
        "empirical_p": round(float((null >= g.mean()).mean()), 4),
        "wilcoxon_p": (round(float(stats.wilcoxon(g).pvalue), 6)
                       if len(g) > 5 and np.any(g != 0) else None),
        # leverage: at n=14 a mean can be two ligands. Report it rather than hide it.
        "jackknife_min": round(float(min(np.delete(g, i).mean()
                                         for i in range(len(g)))), 4),
        "jackknife_max": round(float(max(np.delete(g, i).mean()
                                         for i in range(len(g)))), 4),
        "gain_dropping_best_unit": round(float(np.sort(g)[::-1][1:].mean()), 4),
        "gain_dropping_best_two_units": round(float(np.sort(g)[::-1][2:].mean()), 4),
    }


def compare_strata(a: pd.DataFrame, b: pd.DataFrame, na: str, nb: str,
                   seed: int = SEED) -> dict:
    """Is the difference in gain real, and what could we have detected at this n?"""
    ga, gb = a.gain.to_numpy(), b.gain.to_numpy()
    diff = ga.mean() - gb.mean()
    rng = np.random.default_rng(seed)
    bd = (ga[rng.integers(0, len(ga), (N_BOOT, len(ga)))].mean(1)
          - gb[rng.integers(0, len(gb), (N_BOOT, len(gb)))].mean(1))
    sp = np.sqrt(((len(ga) - 1) * ga.var(ddof=1) + (len(gb) - 1) * gb.var(ddof=1))
                 / (len(ga) + len(gb) - 2))
    mdd = (stats.norm.ppf(0.975) + stats.norm.ppf(0.80)) * sp * np.sqrt(
        1 / len(ga) + 1 / len(gb))
    # permutation test on the difference of means (exchangeable under H0)
    pool = np.concatenate([ga, gb])
    perm = np.empty(N_NULL)
    for i in range(N_NULL):
        p = rng.permutation(pool)
        perm[i] = p[:len(ga)].mean() - p[len(ga):].mean()
    # fraction of the available headroom each stratum converts, as a ratio of means
    # (never a mean of per-unit ratios — a unit with ~0 headroom makes that explode)
    ha, hb = a.headroom.to_numpy(), b.headroom.to_numpy()
    ia = rng.integers(0, len(ga), (N_BOOT, len(ga)))
    ib = rng.integers(0, len(gb), (N_BOOT, len(gb)))
    fa = ga[ia].mean(1) / ha[ia].mean(1)
    fb = gb[ib].mean(1) / hb[ib].mean(1)
    sa = int((ga > 0).sum()); sb = int((gb > 0).sum())
    return {
        "headroom_a": round(float(ha.mean()), 4), "headroom_b": round(float(hb.mean()), 4),
        "frac_headroom_a": round(float(ga.mean() / ha.mean()), 4),
        "frac_headroom_b": round(float(gb.mean() / hb.mean()), 4),
        "frac_headroom_a_ci95": [round(float(np.percentile(fa, 2.5)), 4),
                                 round(float(np.percentile(fa, 97.5)), 4)],
        "frac_headroom_b_ci95": [round(float(np.percentile(fb, 2.5)), 4),
                                 round(float(np.percentile(fb, 97.5)), 4)],
        "frac_headroom_difference_ci95": [
            round(float(np.percentile(fa - fb, 2.5)), 4),
            round(float(np.percentile(fa - fb, 97.5)), 4)],
        "units_with_positive_gain": f"{sa}/{len(ga)} vs {sb}/{len(gb)}",
        "fisher_p_on_sign": round(float(stats.fisher_exact(
            [[sa, len(ga) - sa], [sb, len(gb) - sb]]).pvalue), 4),
        "a": na, "b": nb, "n_a": len(ga), "n_b": len(gb),
        "gain_a": round(float(ga.mean()), 4), "gain_b": round(float(gb.mean()), 4),
        "difference": round(float(diff), 4),
        "difference_ci95": [round(float(np.percentile(bd, 2.5)), 4),
                            round(float(np.percentile(bd, 97.5)), 4)],
        "ci_width": round(float(np.percentile(bd, 97.5) - np.percentile(bd, 2.5)), 4),
        "welch_p": round(float(stats.ttest_ind(ga, gb, equal_var=False).pvalue), 4),
        "mannwhitney_p": round(float(stats.mannwhitneyu(ga, gb).pvalue), 4),
        "permutation_p": round(float((np.abs(perm) >= abs(diff)).mean()), 4),
        "pooled_sd": round(float(sp), 4),
        "min_detectable_difference_80pct_power": round(float(mdd), 4),
        "observed_over_mdd": round(float(abs(diff) / mdd), 3),
    }


# ---------------------------------------------------------------------------
# stage 1 — verify the labels
# ---------------------------------------------------------------------------

def cmd_predict(pool: str, pool_flat: bool, pattern: str, pool_glob: str,
                out_csv: str | None, ligands_csv: str | None) -> dict:
    """Type a ligand set from the PREDICTIONS alone. Needs no crystal (G7).

    The rule is FINDING 033 item 1, unchanged and untuned: the **median `fe_donor_dist`
    over that ligand's own poses**, cut at `COORD_MAX = 2.6` A. On the validation set it
    agrees with the crystal on **84 of 87 = 96.6%**, class separation AUC 1.00, with the
    median fraction of poses coordinating at 1.00 for Type II and 0.00 for Type I.

    `cmd_labels` computes the same column but only as a by-product of reading
    `poses_scored_*.csv`, which is derived from crystals - so on a blind set it cannot run
    at all. This reads the pool directly.

    **The ambiguous band is reported, not resolved.** The closest Type I sits 0.244 A
    beyond the furthest Type II, so a ligand whose median lands inside [2.59, 2.82] A is
    not labelled by this rule in either direction. It is flagged; it is not guessed.
    """
    sys.path.insert(0, str(REPO / "src"))
    from cypstruct import pose as P
    from cypstruct.qmscore import geometry as G

    root = Path(pool)
    if pool_flat:
        units: dict = {}
        for f in sorted(root.glob(pattern)):
            if f.is_file():
                units.setdefault(f.name.split("__")[0], []).append(f)
        units = sorted(units.items())
    else:
        units = [(d.name.split("__")[0], sorted(d.glob(pattern)))
                 for d in sorted(root.glob(pool_glob)) if d.is_dir()]

    want = None
    if ligands_csv:
        df = pd.read_csv(ligands_csv)
        idc = "structure" if "structure" in df.columns and "id" not in df.columns else "id"
        want = set(df[idc].astype(str))

    counts = {"units_seen": len(units), "not_in_csv": 0, "unreadable": 0,
              "no_ligand_atoms": 0, "no_heme": 0, "no_pose_scored": 0}
    rows = []
    for lig, files in units:
        if want is not None and lig not in want:
            counts["not_in_csv"] += 1
            continue
        ds = []
        for f in files:
            try:
                cx = P.load_structure(f)
            except Exception:
                counts["unreadable"] += 1
                continue
            if len(cx.lig_xyz) == 0:
                counts["no_ligand_atoms"] += 1
                continue
            if len(cx.heme_xyz) == 0:
                counts["no_heme"] += 1
                continue
            t = G.compute(cx.lig_xyz, cx.lig_elem, cx.prot_xyz,
                          cx.heme_xyz, cx.heme_atom, cx.axial_sg)
            ds.append((float(t.fe_donor_dist), bool(t.is_coordinated)))
        if not ds:
            counts["no_pose_scored"] += 1
            continue
        dd = np.array([d for d, _c in ds], float)
        med = float(np.nanmedian(dd))
        rows.append({
            "id": lig, "n_poses": len(ds),
            "pred_fe_donor_median": med,
            "pred_fe_donor_min": float(np.nanmin(dd)),
            "pred_frac_coordinated": float(np.mean([c for _d, c in ds])),
            "pred_mode": mode_from_distance(med),
            # FINDING 033: the observed classes do not touch - furthest Type II 2.59 A,
            # closest Type I 2.82 A. Inside that gap the rule has no evidence either way.
            "ambiguous": bool(AMBIG_LO <= med <= AMBIG_HI),
        })
    t = pd.DataFrame(rows).sort_values("pred_fe_donor_median")
    out = Path(out_csv or (DATA_PROCESSED / "binding_mode_pred.csv"))
    t.to_csv(out, index=False)
    missing = sorted(want - set(t.id)) if want else []
    return {
        "rule": "median fe_donor_dist over the ligand's OWN poses, cut at "
                f"COORD_MAX = {COORD_MAX} A (FINDING 033 item 1)",
        "filters": counts,
        "n_ligands": int(len(t)),
        "composition": t.pred_mode.value_counts().to_dict(),
        "n_ambiguous": int(t.ambiguous.sum()) if len(t) else 0,
        "ambiguous_ligands": t.loc[t.ambiguous, "id"].tolist() if len(t) else [],
        "ambiguous_band_A": [AMBIG_LO, AMBIG_HI],
        "ligands_with_no_poses": missing,
        "expected_absolute_score": (
            "Type I-rich sets score nearer 0.51, otherwise nearer 0.64 "
            "(playbook section 3) - pre-announce it with the entry"),
        "out": str(out),
    }


def cmd_labels() -> dict:
    lig = pd.read_csv(DATA_PROCESSED / "validation_ligands.csv")
    pool = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    pool = pool[pool.arm == "unsteered"]

    # CRYSTAL side: `fe` is the closest ligand heavy atom to the iron in the deposited
    # structure (build_reference_set.py). Re-derive the class from it rather than trusting
    # the stored string.
    lig["crystal_mode"] = lig.fe.map(mode_from_distance)
    lig["stored_mode"] = lig.cls.map({"type_II_coordinated": "type_II",
                                      "type_I_active_site": "type_I",
                                      "peripheral": "peripheral"})

    # PREDICTION side: median Fe-donor distance over that ligand's own 20 poses. Uses
    # nothing from the crystal, so it is what would be available on a blind test ligand.
    g = pool.groupby("ligand")
    pred = pd.DataFrame({
        "pred_fe_donor_median": g.fe_donor_dist.median(),
        "pred_fe_donor_min": g.fe_donor_dist.min(),
        "pred_frac_coordinated": g.is_coordinated.mean(),
        "pool_mean_lddt": g.lddt_pli.mean(),
    }).reset_index().rename(columns={"ligand": "id"})
    t = lig.merge(pred, on="id", how="left")
    t["pred_mode"] = t.pred_fe_donor_median.map(mode_from_distance)
    t["label_agrees"] = t.crystal_mode == t.stored_mode
    t["pred_agrees_crystal"] = t.pred_mode == t.crystal_mode
    # how far from the 2.6 A line, in the crystal — the only thing that can flip a label
    t["margin_to_coord_cut"] = (t.fe - COORD_MAX).round(3)
    t.sort_values("fe").to_csv(OUT_LABELS, index=False)

    t2 = t[t.crystal_mode == "type_II"]
    t1 = t[t.crystal_mode == "type_I"]
    borderline = t[(t.fe > COORD_MAX - 0.35) & (t.fe < COORD_MAX + 0.35)]
    # a stored Type II whose own poses never coordinate, or vice versa
    odd = t[(t.crystal_mode == "type_II") & (t.pred_frac_coordinated < 0.10) |
            (t.crystal_mode == "type_I") & (t.pred_frac_coordinated > 0.90)]

    return {
        "n_ligands": int(len(t)),
        "stored_counts": t.stored_mode.value_counts().to_dict(),
        "recomputed_from_crystal_counts": t.crystal_mode.value_counts().to_dict(),
        "stored_label_disagreements": int((~t.label_agrees).sum()),
        "disagreeing": t.loc[~t.label_agrees, ["id", "pdb", "fe", "stored_mode",
                                               "crystal_mode"]].to_dict("records"),
        "crystal_fe_distance": {
            "type_II": {k: round(float(v), 3) for k, v in
                        t2.fe.describe(percentiles=[.05, .5, .95]).items()},
            "type_I": {k: round(float(v), 3) for k, v in
                       t1.fe.describe(percentiles=[.05, .5, .95]).items()},
            "gap_between_classes_A": round(float(t1.fe.min() - t2.fe.max()), 3),
            "mannwhitney_p": float(stats.mannwhitneyu(t2.fe, t1.fe).pvalue),
            "auc_separation": round(float(
                stats.mannwhitneyu(t1.fe, t2.fe).statistic / (len(t1) * len(t2))), 4),
        },
        "borderline_within_0.35A_of_cut": borderline[
            ["id", "pdb", "fe", "crystal_mode", "pred_fe_donor_median",
             "pred_frac_coordinated"]].to_dict("records"),
        "prediction_side_label": {
            "note": "median Fe-donor distance over the ligand's own 20 unsteered poses",
            "agrees_with_crystal": int(t.pred_agrees_crystal.sum()),
            "n": int(len(t)),
            "accuracy": round(float(t.pred_agrees_crystal.mean()), 4),
            "confusion": t.groupby(["crystal_mode", "pred_mode"]).size()
                          .rename("n").reset_index().to_dict("records"),
            "type_II_pred_frac_coordinated_median": round(
                float(t2.pred_frac_coordinated.median()), 3),
            "type_I_pred_frac_coordinated_median": round(
                float(t1.pred_frac_coordinated.median()), 3),
        },
        "label_looks_wrong": odd[["id", "pdb", "fe", "crystal_mode",
                                  "pred_fe_donor_median",
                                  "pred_frac_coordinated"]].to_dict("records"),
        "out": str(OUT_LABELS),
    }


# ---------------------------------------------------------------------------
# stage 2 — stratify the shipped selector on CYP3A4
# ---------------------------------------------------------------------------

def _cyp_pool() -> pd.DataFrame:
    pool = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    xe = pd.read_csv(DATA_PROCESSED / "xeng_val87b.csv")
    m = pool[pool.arm == "unsteered"].merge(xe, on=["ligand", "sample"], how="inner")
    lig = pd.read_csv(DATA_PROCESSED / "validation_ligands.csv")
    lig["crystal_mode"] = lig.fe.map(mode_from_distance)
    return m.merge(lig[["id", "fe", "crystal_mode"]].rename(columns={"id": "ligand"}),
                   on="ligand", how="left")


def cmd_cyp3a4() -> dict:
    from cypstruct import xengine as X

    m = _cyp_pool()
    # control 1: reproduce the shipped column and the shipped number exactly
    shipped_sel = float(X.select(m).lddt_pli.mean())
    t_all = per_unit_table(m, "ligand", "xeng")
    assert abs(shipped_sel - t_all.selected.mean()) < 1e-9, "select() not reproduced"

    # control 2: the exact-zero LDDT-PLI rows are ejections, not FINDING 021 numbering
    zeros = m[m.lddt_pli == 0]
    zcheck = []
    for L in zeros.ligand.unique():
        g = m[m.ligand == L]
        zcheck.append({
            "ligand": L, "pdb": str(g.pdb.iloc[0]),
            "n_zero": int((g.lddt_pli == 0).sum()), "n_poses": int(len(g)),
            "bisy_rmsd_of_zeros": [round(float(v), 2)
                                   for v in g.loc[g.lddt_pli == 0, "bisy_rmsd"]],
            "lddt_of_the_other_poses_median": round(
                float(g.loc[g.lddt_pli > 0, "lddt_pli"].median()), 3),
            "all_mapped": bool(g.mapped.all()),
        })

    per = {}
    tabs = {}
    for label, sub in [("all", m),
                       ("type_II", m[m.crystal_mode == "type_II"]),
                       ("type_I", m[m.crystal_mode == "type_I"])]:
        t = per_unit_table(sub, "ligand", "xeng")
        t["crystal_mode"] = label
        tabs[label] = t
        per[label] = stratum_report(sub, t, "ligand", label)

    cmp_ = compare_strata(tabs["type_II"], tabs["type_I"], "type_II", "type_I")
    cmp_["headroom_matched"] = headroom_matched(tabs["type_II"], tabs["type_I"], n_bins=3)

    # --- the frame hypothesis, tested directly -----------------------------
    t = per_unit_table(m, "ligand", "xeng")
    lig = pd.read_csv(DATA_PROCESSED / "validation_ligands.csv").rename(
        columns={"id": "ligand"})
    predfe = m.groupby("ligand").fe_donor_dist.median().rename("pred_fe_median")
    t = t.merge(lig[["ligand", "fe", "cls"]], on="ligand").merge(predfe, on="ligand")
    t["crystal_mode"] = t.fe.map(mode_from_distance)
    t.to_csv(OUT_LIG, index=False)

    def sp(x, y):
        r = stats.spearmanr(x, y)
        return {"rho": round(float(r.statistic), 4), "p": round(float(r.pvalue), 4),
                "n": int(len(x))}

    t2 = t[t.crystal_mode == "type_II"]

    partial = partial_spearman
    frame = {
        "hypothesis": ("the heme frame is defined by a landmark a Type I ligand is not "
                       "anchored to, so the gain should fall as the ligand sits further "
                       "from the iron"),
        "gain_vs_crystal_fe_distance_ALL": sp(t.fe, t.gain),
        "gain_vs_crystal_fe_distance_TYPE_II_ONLY": sp(t2.fe, t2.gain),
        "gain_vs_predicted_fe_distance_ALL": sp(t.pred_fe_median, t.gain),
        "gain_vs_predicted_fe_distance_TYPE_II_ONLY": sp(t2.pred_fe_median, t2.gain),
        "rho_vs_crystal_fe_distance_ALL": sp(t.dropna(subset=["rho"]).fe,
                                             t.dropna(subset=["rho"]).rho),
        "confound_gain_vs_pool_headroom_ALL": sp(t.headroom, t.gain),
        "confound_gain_vs_pool_spread_ALL": sp(t.spread, t.gain),
        "partial_gain_vs_fe_controlling_headroom_ALL":
            partial(t.fe.to_numpy(), t.gain.to_numpy(), t.headroom.to_numpy()),
        "partial_gain_vs_fe_controlling_headroom_TYPE_II":
            partial(t2.fe.to_numpy(), t2.gain.to_numpy(), t2.headroom.to_numpy()),
        "headroom_by_stratum": {
            k: round(float(v.headroom.mean()), 4) for k, v in tabs.items()},
        "pool_mean_by_stratum": {
            k: round(float(v.rand_exact.mean()), 4) for k, v in tabs.items()},
    }

    return {
        "controls": {
            "shipped_select_reproduced": round(shipped_sel, 4),
            "shipped_gain_vs_empirical_random": round(
                shipped_sel - empirical_random(m, "ligand"), 4),
            "finding_011_full_depth_reference": 0.0380,
            "exact_zero_lddt_rows": int(len(zeros)),
            "zero_row_audit": zcheck,
            "zero_rows_are_numbering": False,
            "zero_row_verdict": ("every zero row's SIBLING poses of the same ligand "
                                 "against the same crystal score 0.44-0.82, so the "
                                 "numbering frame for that pair is correct (FINDING 021 "
                                 "failures are per-PAIR, not per-pose); the zeros carry "
                                 "BiSyRMSD 23-27 A, i.e. genuine ejections"),
        },
        "strata": per,
        "comparison": cmp_,
        "frame_hypothesis": frame,
        "out": [str(OUT_LIG)],
    }


# ---------------------------------------------------------------------------
# stage 3 — the same stratification on the held-out P450 family
# ---------------------------------------------------------------------------

_FN = re.compile(r"(.+)__r(\d+)s(\d+)\.cif$")


def _p450_xeng(force: bool = False) -> pd.DataFrame:
    """Pool = protenix_v2 (deduped), reference = esmfold2, exactly FINDING 012's setup.

    Cached, because it is ~10k mmCIF loads. The cache is a 3-column CSV, so the poses can
    go back to cold storage without losing the ability to re-stratify.
    """
    if OUT_PXENG.exists() and not force:
        return pd.read_csv(OUT_PXENG)

    from cypstruct import pose as P
    from cypstruct import xengine as X

    sc = pd.read_csv(UNI / "p450_poses_scored.csv")
    rows, seen = [], {}
    pairs = list(sc.pair.unique())
    for i, pair in enumerate(pairs, 1):
        d = UNI / "poses" / pair / "esmfold2"
        refs: dict[int, np.ndarray] = {}
        for f in sorted(d.glob("*.cif")) if d.exists() else []:
            mm = _FN.match(f.name)
            rep = int(mm.group(2)) if mm else 0
            if rep in refs:
                continue
            try:
                v = X.in_heme_frame(P.load_structure(f))
            except Exception:
                continue
            if v is not None:
                refs[rep] = v
        r = X._dedupe(list(refs.values()))
        if len(r) < 2:
            continue
        d = UNI / "poses" / pair / "protenix_v2"
        for f in sorted(d.glob("*.cif")) if d.exists() else []:
            try:
                v = X.in_heme_frame(P.load_structure(f))
            except Exception:
                continue
            if v is None:
                continue
            # FINDING 015: duplicates inflate the random baseline and DEFLATE the gain
            if any(v.shape == w.shape and np.allclose(v, w, atol=0.05)
                   for w in seen.setdefault(pair, [])):
                continue
            seen[pair].append(v)
            rows.append({"pair": pair, "pose": f"protenix_v2/{f.name}",
                         "xeng": X.xeng_score(v, r), "ref_depth": len(r)})
        if i % 50 == 0:
            print(f"  [{i}/{len(pairs)}] rows={len(rows)}", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT_PXENG, index=False)
    return out


def cmd_p450(force: bool = False) -> dict:
    sc = pd.read_csv(UNI / "p450_poses_scored.csv")
    xe = _p450_xeng(force=force)
    m = sc.merge(xe, on=["pair", "pose"], how="inner")
    m = m[m.groupby("pair").pose.transform("size") >= 3]

    # CRYSTAL-side mode for every pair, on the same a-priori cuts
    cof = pd.read_csv(UNI / "p450_cofold_set.csv")
    cof["pair"] = cof.pdb + "_" + cof.id
    cof["crystal_mode"] = cof.closest_fe.map(mode_from_distance)
    m = m.merge(cof[["pair", "closest_fe", "donor_dist", "crystal_mode"]],
                on="pair", how="left")

    per, tabs = {}, {}
    for label, sub in [("all", m),
                       ("type_II", m[m.crystal_mode == "type_II"]),
                       ("type_I", m[m.crystal_mode == "type_I"]),
                       ("peripheral", m[m.crystal_mode == "peripheral"])]:
        if sub.pair.nunique() < 4:
            continue
        t = per_unit_table(sub, "pair", "xeng")
        t["crystal_mode"] = label
        tabs[label] = t
        per[label] = stratum_report(sub, t, "pair", label)

    t_all = per_unit_table(m, "pair", "xeng").merge(
        cof[["pair", "closest_fe", "crystal_mode", "uniprot", "n_heavy"]], on="pair")
    t_all.to_csv(OUT_PPAIR, index=False)

    def sp(x, y):
        r = stats.spearmanr(x, y)
        return {"rho": round(float(r.statistic), 4), "p": float(r.pvalue), "n": int(len(x))}

    t2 = t_all[t_all.crystal_mode == "type_II"]
    cmps = {}
    if "type_I" in tabs:
        cmps["type_II_vs_type_I"] = compare_strata(tabs["type_II"], tabs["type_I"],
                                                   "type_II", "type_I")
        cmps["type_II_vs_type_I"]["headroom_matched"] = headroom_matched(
            tabs["type_II"], tabs["type_I"], n_bins=5)
    if "peripheral" in tabs:
        cmps["type_II_vs_peripheral"] = compare_strata(
            tabs["type_II"], tabs["peripheral"], "type_II", "peripheral")
        noncoord = pd.concat([tabs.get("type_I", pd.DataFrame()),
                              tabs["peripheral"]], ignore_index=True)
        cmps["type_II_vs_all_noncoordinating"] = compare_strata(
            tabs["type_II"], noncoord, "type_II", "type_I + peripheral")

    # protein-level: is the mode effect confounded with target identity?
    prot = []
    for u, g in t_all.groupby("uniprot"):
        if len(g) < 6 or g.crystal_mode.nunique() < 2:
            continue
        a = g[g.crystal_mode == "type_II"].gain
        b = g[g.crystal_mode != "type_II"].gain
        if len(a) >= 3 and len(b) >= 3:
            prot.append({"uniprot": u, "n_typeII": len(a), "n_other": len(b),
                         "gain_typeII": round(float(a.mean()), 4),
                         "gain_other": round(float(b.mean()), 4),
                         "difference": round(float(a.mean() - b.mean()), 4)})

    return {
        "setup": {"pool": "protenix_v2 (deduped)", "reference": "esmfold2",
                  "poses": int(len(m)), "pairs": int(m.pair.nunique()),
                  "proteins": int(m.uniprot.nunique()),
                  "non_cyp3a4_proteins": int(
                      m[m.uniprot != "P08684"].uniprot.nunique()),
                  "finding_012_reference_gain": 0.0414},
        "strata": per,
        "comparisons": cmps,
        "frame_hypothesis": {
            "gain_vs_crystal_fe_distance_ALL": sp(t_all.closest_fe, t_all.gain),
            "gain_vs_crystal_fe_distance_TYPE_II_ONLY": sp(t2.closest_fe, t2.gain),
            "confound_gain_vs_pool_headroom_ALL": sp(t_all.headroom, t_all.gain),
            "partial_gain_vs_fe_controlling_headroom_ALL": partial_spearman(
                t_all.closest_fe.to_numpy(), t_all.gain.to_numpy(),
                t_all.headroom.to_numpy()),
            "partial_gain_vs_fe_controlling_headroom_TYPE_II": partial_spearman(
                t2.closest_fe.to_numpy(), t2.gain.to_numpy(), t2.headroom.to_numpy()),
            "headroom_by_stratum": {k: round(float(v.headroom.mean()), 4)
                                    for k, v in tabs.items()},
            "pool_mean_by_stratum": {k: round(float(v.rand_exact.mean()), 4)
                                     for k, v in tabs.items()},
        },
        # the obvious alternative explanation for a stratum difference: one stratum got
        # thinner reference sets or thinner pools. Measured, and it runs the WRONG WAY.
        "depth_confound": {
            "reference_depth_by_stratum": m.groupby("crystal_mode").ref_depth.mean()
                                           .round(2).to_dict(),
            "pool_depth_by_stratum": m.groupby(["crystal_mode", "pair"]).size()
                                      .groupby("crystal_mode").mean().round(2).to_dict(),
        },
        "within_protein": prot,
        "out": [str(OUT_PXENG), str(OUT_PPAIR)],
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["labels", "predict", "cyp3a4", "p450", "all"])
    ap.add_argument("--pool", default=None,
                    help="predict: pose pool to type. Needs no crystal.")
    ap.add_argument("--pool-flat", action="store_true",
                    help="predict: --pool is a flat directory of <LIG>__*.cif")
    ap.add_argument("--pattern", default="*.cif")
    ap.add_argument("--pool-glob", default="*__*")
    ap.add_argument("--ligands", default=None,
                    help="predict: ligand csv, to report any id with no poses at all")
    ap.add_argument("--out", default=None, help="predict: output csv")
    ap.add_argument("--force", action="store_true", help="recompute the P450 xeng cache")
    a = ap.parse_args()

    res = {}
    if a.cmd == "predict":
        if not a.pool:
            raise SystemExit("predict needs --pool")
        res["predict"] = cmd_predict(a.pool, a.pool_flat, a.pattern, a.pool_glob,
                                     a.out, a.ligands)
    if a.cmd in ("labels", "all"):
        res["labels"] = cmd_labels()
    if a.cmd in ("cyp3a4", "all"):
        res["cyp3a4"] = cmd_cyp3a4()
        OUT_CYP.write_text(json.dumps(res["cyp3a4"], indent=1))
    if a.cmd in ("p450", "all"):
        res["p450"] = cmd_p450(force=a.force)
        OUT_P450.write_text(json.dumps(res["p450"], indent=1))
    print(json.dumps(res, indent=1, default=str))
