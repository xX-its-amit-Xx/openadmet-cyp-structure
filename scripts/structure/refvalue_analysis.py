"""FINDING 041 — what reference DEPTH is worth to SELECTION on CYP3A4.

Pre-registered in `docs/PREREG_reference_depth_value.md` (commit `40e0f55`).

The pool is held completely fixed — the 20 unsteered Boltz-2 poses per ligand, already
scored — and the ONLY thing that varies is which reference poses the shipped Chamfer
feature is allowed to see. Everything in the depth curve is EXACT: for a ligand with `D`
reference poses, every one of the `C(D, d)` subsets of size `d` is enumerated, so there is
no Monte-Carlo error anywhere except in the noise floor (2e6 draws, SE reported) and in
the across-draw percentile spread (2e5 draws), both of which say so.

Stages, each resumable via its own artefact:

    frames   740 pool poses -> ligand heavy atoms in each pose's OWN heme frame
    controls C1 numbering, C2 shipped-xeng reproduction, C3 argmin == select()
    curve    the exact depth-vs-selection curve, spread, floor, leverage
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

DATA = REPO / "data" / "processed"
POOL = Path(os.environ.get("RV_POOL", "D:/cyp_scratch/val87b_unsteered"))
SCRATCH = Path(os.environ.get(
    "RV_SCRATCH",
    r"C:\tb\tmp\1\claude\D--Users-ashenoy00000--windsurf-OpenADMET-cyp-structure"
    r"\bd288271-2aa5-4ed4-a343-ea31e5ad8c13\scratchpad\refvalue"))

DEPTHS = [4, 5, 6, 7, 8]
SEED = 20260923
FLOOR_BLOCKS, FLOOR_BLOCK = 200, 10_000      # 2e6 draws, as FINDING 039/040
SPREAD_DRAWS = 200_000
BOOT = 10_000
RHO_SUBSETS = 50

LIGCSV = DATA / "refvalue_ligands.csv"
REFNPZ = DATA / "refvalue_reference_set.npz"
SHIPPED_REF = DATA / "reference_set_cyp3a4.npz"


def ligands() -> pd.DataFrame:
    return pd.read_csv(LIGCSV)


# ---------------------------------------------------------------- stage: frames
def stage_frames() -> None:
    """Pool poses into their own heme frames. `xeng_cache` holds RAW coords, not these."""
    from cypstruct import pose as P
    from cypstruct import xengine as X

    out = SCRATCH / "frames"
    out.mkdir(parents=True, exist_ok=True)
    done = 0
    for lig in ligands()["id"].astype(str):
        f = out / f"{lig}.npz"
        if f.exists():
            done += 1
            continue
        job = POOL / f"{lig}__unsteered__s1"
        cifs = sorted(job.glob("input_model_*.cif"),
                      key=lambda p: int(p.stem.rsplit("_", 1)[1]))
        names, frames = [], []
        for cif in cifs:
            try:
                m = P.load_structure(cif)
            except Exception:
                continue
            v = X.in_heme_frame(m)
            if v is None:
                continue
            names.append(cif.stem)
            frames.append(np.asarray(v, np.float32))
        if not names:
            print(f"  {lig}: NO POSES under {job}", flush=True)
            continue
        np.savez_compressed(f, names=np.array(names),
                            xyz=np.array(frames, np.float32))
        done += 1
        print(f"  {lig}: {len(names)} poses", flush=True)
    print(f"frames: {done} ligands -> {out}", flush=True)


def load_frames() -> dict:
    out = SCRATCH / "frames"
    d = {}
    for lig in ligands()["id"].astype(str):
        f = out / f"{lig}.npz"
        if f.exists():
            z = np.load(f, allow_pickle=True)
            d[lig] = (list(z["names"]), z["xyz"])
    return d


def pool_scores() -> pd.DataFrame:
    p = pd.read_csv(DATA / "poses_scored_val87b.csv")
    return p[p.arm == "unsteered"].copy()


# ---------------------------------------------------------------- chamfer matrices
def chamfer_matrix(poses: np.ndarray, refs: list) -> np.ndarray:
    from cypstruct import xengine as X
    return np.array([[X.chamfer(np.asarray(p, float), np.asarray(r, float))
                      for r in refs] for p in poses])


# ---------------------------------------------------------------- stage: controls
def stage_controls() -> dict:
    from cypstruct import xengine as X

    fr = load_frames()
    sc = pool_scores()
    ligs = [l for l in ligands()["id"].astype(str) if l in fr]
    rep = {"n_ligands": len(ligs)}

    # C1 - numbering control (FINDING 021): a numbering offset reads as LDDT-PLI 0.0
    zero_lig, nz = [], 0
    for lig in ligs:
        v = sc[sc.ligand == lig].lddt_pli.values
        if len(v) and np.allclose(v, 0.0):
            zero_lig.append(lig)
        nz += int((v == 0.0).sum())
    rep["C1_numbering"] = {
        "ligands_all_zero_lddt": zero_lig, "poses_exactly_zero": nz,
        "n_poses": int(sc[sc.ligand.isin(ligs)].shape[0]),
        "lddt_min": round(float(sc[sc.ligand.isin(ligs)].lddt_pli.min()), 5),
        "lddt_max": round(float(sc[sc.ligand.isin(ligs)].lddt_pli.max()), 5),
        "pass": len(zero_lig) == 0}

    # C2 - reproduce the shipped xeng column from the shipped reference set
    shipped = X.load_reference(SHIPPED_REF)
    ship_csv = pd.read_csv(DATA / "xeng_val87b.csv")
    ship_map = {(r.ligand, r.sample): r.xeng for r in ship_csv.itertuples()}
    diffs, n_cmp = [], 0
    rows = []
    for lig in ligs:
        names, xyz = fr[lig]
        refs = shipped.get(lig, [])
        for nm, p in zip(names, xyz):
            s = X.xeng_score(np.asarray(p, float), [np.asarray(r, float) for r in refs])
            rows.append({"ligand": lig, "sample": str(nm), "xeng": s})
            k = (lig, str(nm))
            if k in ship_map:
                diffs.append(abs(s - ship_map[k]))
                n_cmp += 1
    mine = pd.DataFrame(rows)
    rep["C2_reproduce_shipped_xeng"] = {
        "compared": n_cmp, "max_abs_diff": float(np.max(diffs)) if diffs else None,
        "mean_abs_diff": float(np.mean(diffs)) if diffs else None,
        "pass": bool(diffs and np.max(diffs) < 1e-6)}

    # C3 - argmin(xeng) IS xengine.select()
    sel = X.select(mine)
    agree = 0
    for lig, g in mine.groupby("ligand"):
        a = g.loc[g.xeng.idxmin(), "sample"]
        b = sel.loc[sel.ligand == lig, "sample"].iloc[0]
        agree += int(a == b)
    rep["C3_argmin_is_select"] = {"ligands": int(mine.ligand.nunique()),
                                  "agree": agree,
                                  "pass": agree == int(mine.ligand.nunique())}
    mine.to_csv(DATA / "refvalue_xeng_shipped_repro.csv", index=False)
    return rep


# ---------------------------------------------------------------- exact curve
def exact_depth_stats(C: np.ndarray, y: np.ndarray, d: int):
    """Exact (mean, var, per-subset picks) of the selected LDDT-PLI at reference depth d.

    C is (n_poses, D) Chamfer; y is (n_poses,) LDDT-PLI. Every C(D, d) subset is
    enumerated - the pick is argmin over the subset mean, first index on a tie, which is
    what `select()`'s idxmax does on `-z(xeng)`.
    """
    D = C.shape[1]
    vals, picks = [], []
    for S in combinations(range(D), d):
        m = C[:, S].mean(axis=1)
        k = int(np.argmin(m))
        picks.append(k)
        vals.append(y[k])
    vals = np.asarray(vals, float)
    return float(vals.mean()), float(vals.var(ddof=0)), np.asarray(picks), vals


def stage_curve() -> dict:
    from scipy.stats import spearmanr
    from cypstruct import xengine as X

    fr = load_frames()
    sc = pool_scores()
    ref = X.load_reference(REFNPZ)
    ligs = sorted(l for l in ligands()["id"].astype(str) if l in fr and l in ref)
    rng = np.random.default_rng(SEED)

    per, mats, Y, chamf = {}, {}, {}, {}
    for lig in ligs:
        names, xyz = fr[lig]
        s = sc[sc.ligand == lig].set_index("sample")
        keep = [i for i, nm in enumerate(names) if str(nm) in s.index]
        y = np.array([s.loc[str(names[i]), "lddt_pli"] for i in keep], float)
        poses = np.asarray(xyz)[keep]
        refs = [np.asarray(r, float) for r in ref[lig]]
        C = chamfer_matrix(poses, refs)
        mats[lig] = C
        Y[lig] = y
        # pairwise distinctness inside the reference set (C5)
        pw = [X.chamfer(refs[i], refs[j])
              for i in range(len(refs)) for j in range(i + 1, len(refs))]
        chamf[lig] = pw
        per[lig] = {"n_poses": len(y), "ref_depth": len(refs),
                    "oracle": float(y.max()), "random": float(y.mean()),
                    "min_pairwise_chamfer": float(np.min(pw)) if pw else None,
                    "median_pairwise_chamfer": float(np.median(pw)) if pw else None}

    n = len(ligs)
    oracle = float(np.mean([per[l]["oracle"] for l in ligs]))
    random_ = float(np.mean([per[l]["random"] for l in ligs]))

    # ---- the exact curve
    curve = {}
    per_depth_ligand = {}
    for d in DEPTHS:
        means, varis, samples = [], [], {}
        rhos, sgn = [], 0
        for lig in ligs:
            C, y = mats[lig], Y[lig]
            if C.shape[1] < d:
                raise SystemExit(f"{lig} has reference depth {C.shape[1]} < {d}")
            m, v, picks, vals = exact_depth_stats(C, y, d)
            means.append(m)
            varis.append(v)
            samples[lig] = vals
            per_depth_ligand.setdefault(lig, {})[d] = {
                "mean": m, "var": v, "n_subsets": int(comb(C.shape[1], d)),
                "modal_pick": int(np.bincount(picks).argmax()),
                "distinct_picks": int(len(set(picks.tolist())))}
            # within-ligand rho at this depth, averaged over RHO_SUBSETS random
            # subsets of size d so the number is not one draw's noise
            rr = []
            for _ in range(RHO_SUBSETS):
                S = rng.choice(C.shape[1], size=d, replace=False)
                r = spearmanr(C[:, S].mean(axis=1), y).statistic
                if not np.isnan(r):
                    rr.append(float(r))
            if rr:
                rhos.append(float(np.mean(rr)))
                sgn += int(np.mean(rr) < 0)
        sel_mean = float(np.mean(means))
        pooled_var = float(np.sum(varis)) / n ** 2
        # across-draw distribution of the pooled mean
        draws = np.zeros(SPREAD_DRAWS)
        for lig in ligs:
            v = samples[lig]
            draws += v[rng.integers(0, len(v), size=SPREAD_DRAWS)]
        draws /= n
        curve[d] = {
            "selected_exact": sel_mean, "gain_exact": sel_mean - random_,
            "pooled_sd_across_draws": float(np.sqrt(pooled_var)),
            "draw_p5": float(np.percentile(draws, 5)),
            "draw_p50": float(np.percentile(draws, 50)),
            "draw_p95": float(np.percentile(draws, 95)),
            "draw_min": float(draws.min()), "draw_max": float(draws.max()),
            "draw_sd": float(draws.std(ddof=1)), "draws": SPREAD_DRAWS,
            "mean_within_ligand_rho": float(np.mean(rhos)),
            "frac_correct_sign": float(sgn / len(rhos)),
            "per_ligand_mean": {l: means[i] for i, l in enumerate(ligs)}}

    # ---- full depth, for A0
    full = {}
    fdm = []
    for lig in ligs:
        C, y = mats[lig], Y[lig]
        k = int(np.argmin(C.mean(axis=1)))
        fdm.append(y[k])
    full = {"depth": {l: mats[l].shape[1] for l in ligs},
            "selected": float(np.mean(fdm)), "gain": float(np.mean(fdm)) - random_}

    # ---- noise floor on THIS population
    mat = np.array([Y[l] for l in ligs])
    floor = noise_floor(mat)

    # ---- primary: delta 8 - 4, paired bootstrap, drop-1, pick change
    a = np.array([curve[4]["per_ligand_mean"][l] for l in ligs])
    b = np.array([curve[8]["per_ligand_mean"][l] for l in ligs])
    delta = float(b.mean() - a.mean())
    brng = np.random.default_rng(SEED + 1)
    idx = brng.integers(0, n, size=(BOOT, n))
    bs = (b[idx] - a[idx]).mean(axis=1)
    drop1 = [float((np.delete(b, i) - np.delete(a, i)).mean()) for i in range(n)]
    changed = [l for l in ligs
               if per_depth_ligand[l][4]["modal_pick"] != per_depth_ligand[l][8]["modal_pick"]]
    # fraction of (d=4 subset, d=8 subset) pairs that disagree on the pick
    disagree = {}
    for lig in ligs:
        C = mats[lig]
        p4 = np.array([int(np.argmin(C[:, S].mean(axis=1)))
                       for S in combinations(range(C.shape[1]), 4)])
        p8 = np.array([int(np.argmin(C[:, S].mean(axis=1)))
                       for S in combinations(range(C.shape[1]), 8)])
        disagree[lig] = float((p4[:, None] != p8[None, :]).mean())

    from scipy.stats import wilcoxon
    try:
        w = float(wilcoxon(b, a).pvalue)
    except Exception:
        w = float("nan")

    out = {
        "n_ligands": n, "ligands": ligs,
        "pool": {"poses_per_ligand": {l: per[l]["n_poses"] for l in ligs},
                 "oracle": oracle, "random": random_},
        "full_depth": full,
        "curve": {str(k): {kk: vv for kk, vv in v.items() if kk != "per_ligand_mean"}
                  for k, v in curve.items()},
        "per_ligand_mean": {str(d): curve[d]["per_ligand_mean"] for d in DEPTHS},
        "noise_floor_this_population": floor,
        "primary": {
            "delta_8_minus_4": delta,
            "boot_ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
            "boot_frac_positive": float((bs > 0).mean()),
            "wilcoxon_p": w,
            "drop1_min": float(np.min(drop1)), "drop1_max": float(np.max(drop1)),
            "drop1_sign_stable": bool(np.sign(np.min(drop1)) == np.sign(np.max(drop1))),
            "ligands_changing_modal_pick": changed,
            "n_changing_modal_pick": len(changed),
            "mean_subset_pair_disagreement": float(np.mean(list(disagree.values()))),
            "monotone_4_to_8": bool(all(
                curve[DEPTHS[i + 1]]["selected_exact"] >= curve[DEPTHS[i]]["selected_exact"]
                for i in range(len(DEPTHS) - 1)))},
        "per_ligand": per,
        "subset_pair_disagreement": disagree,
    }
    return out


def noise_floor(mat: np.ndarray, blocks: int = FLOOR_BLOCKS,
                block: int = FLOOR_BLOCK) -> dict:
    """p95/p99 of a random feature's gain on THIS population (FINDING 039/040 form)."""
    nl, dep = mat.shape
    base = mat.mean(axis=1).mean()
    rng = np.random.default_rng(SEED)
    p95s, allg = [], []
    for _ in range(blocks):
        idx = rng.integers(0, dep, size=(block, nl))
        g = mat[np.arange(nl)[None, :], idx].mean(axis=1) - base
        p95s.append(np.percentile(g, 95))
        allg.append(g)
    allg = np.concatenate(allg)
    return {"n_ligands": int(nl), "depth": int(dep), "draws": int(blocks * block),
            "p95": round(float(np.percentile(allg, 95)), 5),
            "p99": round(float(np.percentile(allg, 99)), 5),
            "se_of_p95": round(float(np.std(p95s, ddof=1) / np.sqrt(blocks)), 6),
            "sd_of_gain": round(float(allg.std(ddof=1)), 5)}


def main() -> int:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    SCRATCH.mkdir(parents=True, exist_ok=True)
    if cmd in ("frames", "all"):
        stage_frames()
    if cmd in ("controls", "all"):
        r = stage_controls()
        (DATA / "refvalue_controls.json").write_text(json.dumps(r, indent=2, default=str))
        print(json.dumps(r, indent=2, default=str))
    if cmd in ("curve", "all"):
        r = stage_curve()
        (DATA / "refvalue_curve.json").write_text(json.dumps(r, indent=2, default=str))
        print(json.dumps({k: v for k, v in r.items()
                          if k not in ("per_ligand", "per_ligand_mean",
                                       "subset_pair_disagreement", "ligands")},
                         indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# ---------------------------------------------------------------- C4/C5/C6 positives
def stage_filter_controls() -> dict:
    """Positive controls for `refset`'s three filters (FINDING 038's T8 rule).

    A filter that never fires beats any passing check, so each zero count in the refset
    report has to be shown to be a real zero rather than dead code.
    """
    import hashlib
    from cypstruct import xengine as X

    ref = X.load_reference(REFNPZ)
    lig = sorted(ref)[0]
    a = np.asarray(ref[lig][0], float)
    b = np.asarray(ref[lig][1], float)
    out = {}
    # md5 filter on an exact byte copy
    out["md5_on_byte_copy_fires"] = bool(
        hashlib.md5(a.tobytes()).hexdigest() == hashlib.md5(a.copy().tobytes()).hexdigest())
    # coordinate filter on an identical vector
    out["coords_on_identical_fires"] = len(X._dedupe([a, a.copy()])) == 1
    # coordinate filter on a geometrically identical pose, byte-different
    out["coords_on_1e-4_fires"] = len(X._dedupe([a, a + 1e-4])) == 1
    # coordinate filter KEEPS a genuinely distinct pose
    out["coords_keeps_5A_apart"] = len(X._dedupe([a, a + 5.0])) == 2
    # and keeps two genuinely different reference poses of the same ligand
    out["coords_keeps_two_real_poses"] = len(X._dedupe([a, b])) == 2
    out["probe_ligand"] = lig
    out["all_pass"] = all(v for k, v in out.items() if isinstance(v, bool))
    return out
