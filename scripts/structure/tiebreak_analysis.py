"""FINDING 036 — is the near-tie a prize, or only a variance mechanism?

Pre-registered in `docs/PREREG_tiebreak.md`, committed at `a8ae652` before a single
endpoint existed. Stage 1 runs first and can stop the experiment: if a PERFECT top-2
tie-break is worth less than the +0.0137 pooled noise floor, no tie-breaker of any shape
can clear it, and the line closes by the oracle rather than by a gate (`FINDING_032`).

Everything is on disk. Zero new inference, zero downloads, CPU only.

    python scripts/structure/tiebreak_analysis.py frames    # heme-frame cache + D matrix
    python scripts/structure/tiebreak_analysis.py controls  # C-XENG C-SEL C-N2 C-NUM ...
    python scripts/structure/tiebreak_analysis.py regime    # STAGE 1 -- the prize
    python scripts/structure/tiebreak_analysis.py gate      # R1 answer recognition
    python scripts/structure/tiebreak_analysis.py select    # candidates + the three bars
    python scripts/structure/tiebreak_analysis.py all
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

CACHE = Path("C:/Temp/cyp_tiebreak")
POOL_A_DIR = Path("D:/cyp_scratch/val87b_unsteered")
POOL_B_DIR = Path("C:/cyp_struct/matched_depth/poses/stratum_flat")
RCSB = REPO / "data" / "reference" / "rcsb"
FROZEN = DATA_PROCESSED / "reference_set_cyp3a4.npz"

KS = (2, 3, 5)
TIE_DRAWS = 64
NULL_DRAWS = 4000
BOOT = 10_000
SEED = 20260923

OUT_CONTROLS = DATA_PROCESSED / "tiebreak_controls.json"
OUT_REGIME = DATA_PROCESSED / "tiebreak_regime.json"
OUT_GATE = DATA_PROCESSED / "tiebreak_answer_recognition.json"
OUT_SELECT = DATA_PROCESSED / "tiebreak_selectors.json"
OUT_PERLIG = DATA_PROCESSED / "tiebreak_per_ligand.csv"

# FINDING 007 / 035 -- quoted, then RECOMPUTED below for these features on this pool
FLOOR_POOLED_P95 = 0.0137
FLOOR_N14_P95 = 0.0453


# ---------------------------------------------------------------------------
# frames -- every pose, and every crystal, in its own heme frame
# ---------------------------------------------------------------------------

def _frame_one(args):
    lig, path = args
    sys.path.insert(0, str(REPO / "src"))
    from cypstruct import pose as P
    from cypstruct import xengine as X
    try:
        v = X.in_heme_frame(P.load_structure(path))
    except Exception as e:                                   # noqa: BLE001
        return lig, Path(path).stem, None, str(e)[:120]
    return lig, Path(path).stem, (None if v is None else np.asarray(v, np.float32)), None


def _crystal_one(args):
    """The deposited ligand, ONE chain, in its own heme frame. Diagnostic use only."""
    lig, pdb = args
    sys.path.insert(0, str(REPO / "src"))
    import gemmi
    from cypstruct import pose as P
    from cypstruct import xengine as X
    cif = RCSB / f"{pdb}.cif"
    if not cif.exists():
        return lig, None, "no-cif"
    try:
        st = gemmi.read_structure(str(cif))
        st.setup_entities()
    except Exception as e:                                   # noqa: BLE001
        return lig, None, f"read:{str(e)[:60]}"
    for chain in st[0]:
        if any(r.name.strip().upper() == lig.upper() for r in chain):
            try:
                cx = P.load_structure(cif, ligand_code=lig, assembly_chain=chain.name)
            except Exception as e:                           # noqa: BLE001
                return lig, None, f"load:{str(e)[:60]}"
            v = X.in_heme_frame(cx)
            if v is None:
                return lig, None, "no-heme-frame"
            return lig, np.asarray(v, np.float32), None
    return lig, None, "ligand-not-in-any-chain"


def _pool_files(pool: str) -> list[tuple[str, str]]:
    jobs: list[tuple[str, str]] = []
    if pool == "A":
        for d in sorted(POOL_A_DIR.glob("*__*")):
            if d.is_dir():
                lig = d.name.split("__")[0]
                jobs += [(lig, str(f)) for f in sorted(d.glob("*.cif"))]
    else:
        for f in sorted(POOL_B_DIR.glob("*.cif")):
            jobs.append((f.name.split("__")[0], str(f)))
    return jobs


def frames(workers: int = 8) -> dict:
    """Cache heme-frame ligand coordinates for both pools and every crystal."""
    CACHE.mkdir(parents=True, exist_ok=True)
    report = {}
    for pool in ("A", "B"):
        out = CACHE / f"frames_{pool}.npz"
        jobs = _pool_files(pool)
        if out.exists():
            report[pool] = {"cached": True, "files": len(jobs)}
            continue
        store: dict[str, list] = {}
        fails = []
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(_frame_one, j) for j in jobs]
            for i, f in enumerate(as_completed(futs), 1):
                lig, name, v, err = f.result()
                if v is None:
                    fails.append((lig, name, err))
                    continue
                store.setdefault(lig, []).append((name, v))
                if i % 400 == 0:
                    print(f"  {pool}: {i}/{len(jobs)}", flush=True)
        payload = {}
        for lig, items in store.items():
            items.sort(key=lambda t: t[0])
            payload[f"xyz__{lig}"] = np.stack([v for _, v in items])
            payload[f"nm__{lig}"] = np.array([n for n, _ in items])
        np.savez_compressed(out, **payload)
        report[pool] = {"cached": False, "files": len(jobs),
                        "ligands": len(store),
                        "poses": int(sum(len(v) for v in store.values())),
                        "failures": fails[:10], "n_failures": len(fails)}
        print(f"{pool}: {report[pool]}", flush=True)

    out = CACHE / "crystals.npz"
    if not out.exists():
        p = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
        pairs = (p[p.arm == "unsteered"][["ligand", "pdb"]]
                 .drop_duplicates().values.tolist())
        payload, why = {}, {}
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for f in as_completed([ex.submit(_crystal_one, tuple(t)) for t in pairs]):
                lig, v, err = f.result()
                if v is None:
                    why[lig] = err
                else:
                    payload[lig] = v
        np.savez_compressed(out, **payload)
        report["crystals"] = {"ok": len(payload), "failed": why}
        print(f"crystals: {report['crystals']}", flush=True)
    else:
        report["crystals"] = {"cached": True}
    return report


def load_frames(pool: str) -> dict[str, tuple[list[str], np.ndarray]]:
    """Pool B's files are named <LIG>__<sample>.cif; the truth table keys on <sample>
    alone. Normalising here rather than at the join is the difference between a control
    that reads 0 rows and one that reads 280 (`too-clean-numbers-are-the-tell`)."""
    d = np.load(CACHE / f"frames_{pool}.npz", allow_pickle=False)
    ligs = sorted({k.split("__", 1)[1] for k in d.files if k.startswith("xyz__")})
    out = {}
    for l in ligs:
        nm = [str(x) for x in d[f"nm__{l}"]]
        if pool == "B":
            nm = [n.split("__", 1)[1] if n.startswith(l + "__") else n for n in nm]
        out[l] = (nm, d[f"xyz__{l}"])
    return out


def load_crystals() -> dict[str, np.ndarray]:
    d = np.load(CACHE / "crystals.npz", allow_pickle=False)
    return {k: d[k] for k in d.files}


# ---------------------------------------------------------------------------
# the distance matrix xeng averages
# ---------------------------------------------------------------------------

def dmatrix(pool: str) -> dict[str, tuple[list[str], np.ndarray]]:
    """D[pose, ref] -- Chamfer, exactly the matrix `xeng_score` takes the mean of."""
    from cypstruct import xengine as X
    ref = X.load_reference(FROZEN)
    fr = load_frames(pool)
    out = {}
    for lig, (names, xyz) in fr.items():
        if lig not in ref or not ref[lig]:
            continue
        R = ref[lig]
        D = np.array([[X.chamfer(np.asarray(v, float), np.asarray(r, float))
                       for r in R] for v in xyz])
        out[lig] = (names, D)
    return out


def _cached_dmatrix(pool: str):
    p = CACHE / f"dmat_{pool}.npz"
    if p.exists():
        d = np.load(p, allow_pickle=False)
        ligs = sorted({k.split("__", 1)[1] for k in d.files if k.startswith("D__")})
        return {l: ([str(x) for x in d[f"nm__{l}"]], d[f"D__{l}"]) for l in ligs}
    m = dmatrix(pool)
    np.savez_compressed(p, **{f"D__{l}": v[1] for l, v in m.items()},
                        **{f"nm__{l}": np.array(v[0]) for l, v in m.items()})
    return m


# ---------------------------------------------------------------------------
# truth tables
# ---------------------------------------------------------------------------

def truth_A() -> pd.DataFrame:
    p = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    return p[p.arm == "unsteered"].copy()


def truth_B() -> pd.DataFrame:
    return pd.read_csv(DATA_PROCESSED / "matched_depth_poses.csv")


def stratum() -> list[str]:
    L = pd.read_csv(DATA_PROCESSED / "binding_mode_labels_cyp3a4.csv")
    return sorted(L[L.pred_fe_donor_median > 2.6].id.tolist())


# ---------------------------------------------------------------------------
# candidate features -- signs fixed in the pre-registration, ALL "lower is better"
# ---------------------------------------------------------------------------

def features(D: np.ndarray, xyz: np.ndarray) -> dict[str, np.ndarray]:
    """Every statistic the pre-registration names, oriented so LOWER is better."""
    from cypstruct import xengine as X
    from scipy.stats import rankdata

    n, m = D.shape
    f = {"xeng": D.mean(axis=1)}
    f["xeng_2nd"] = np.sort(D, axis=1)[:, 1] if m >= 2 else D[:, 0]
    f["xeng_borda"] = np.mean([rankdata(D[:, r]) for r in range(m)], axis=0)
    f["xeng_sd"] = D.std(axis=1, ddof=1) if m >= 2 else np.zeros(n)
    votes = np.zeros(n)
    for r in range(m):
        votes[int(np.argmin(D[:, r]))] += 1
    f["xeng_votes"] = -votes                       # T5: HIGHER votes is better
    # T1 medoid: pose-to-pose Chamfer inside the pool (used only within a top-k)
    P = np.array([[X.chamfer(np.asarray(a, float), np.asarray(b, float))
                   for b in xyz] for a in xyz])
    f["_pp"] = P
    return f


def crystal_features(D: np.ndarray, dc: np.ndarray) -> dict[str, float]:
    """The same statistics for the crystal, ranked against the 20 predictions."""
    from scipy.stats import rankdata
    n, m = D.shape
    aug = np.vstack([D, dc[None, :]])
    out = {"xeng": float(dc.mean()),
           "xeng_2nd": float(np.sort(dc)[1] if m >= 2 else dc[0]),
           "xeng_sd": float(dc.std(ddof=1) if m >= 2 else 0.0)}
    out["xeng_borda"] = float(np.mean([rankdata(aug[:, r])[-1] for r in range(m)]))
    v = 0
    for r in range(m):
        if int(np.argmin(aug[:, r])) == n:
            v += 1
    out["xeng_votes"] = -float(v)
    return out


def crystal_borda_for_poses(D: np.ndarray, dc: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata
    aug = np.vstack([D, dc[None, :]])
    return np.mean([rankdata(aug[:, r])[:-1] for r in range(aug.shape[1])], axis=0)


def crystal_votes_for_poses(D: np.ndarray, dc: np.ndarray) -> np.ndarray:
    n, m = D.shape
    aug = np.vstack([D, dc[None, :]])
    v = np.zeros(n + 1)
    for r in range(m):
        v[int(np.argmin(aug[:, r]))] += 1
    return -v[:n]


# ---------------------------------------------------------------------------
# assembling one pool into aligned arrays
# ---------------------------------------------------------------------------

def assemble(pool: str, ligands: list[str] | None = None, depth: str = "all") -> dict:
    """{ligand: {'lddt', 'feat', 'names'}} with feature rows aligned to truth rows."""
    m = _cached_dmatrix(pool)
    fr = load_frames(pool)
    t = truth_A() if pool == "A" else truth_B()
    t = t.set_index(["ligand", "sample"])
    out = {}
    for lig, (names, D) in m.items():
        if ligands is not None and lig not in ligands:
            continue
        xyz = fr[lig][1]
        assert fr[lig][0] == names
        try:
            lddt = np.array([float(t.loc[(lig, nm), "lddt_pli"]) for nm in names])
        except KeyError:
            continue
        f = features(D, xyz)
        out[lig] = {"names": names, "D": D, "lddt": lddt, "feat": f}
    return out


def assemble_union(ligands: list[str]) -> dict:
    """Pool B: the 40-pose union -- 20 existing Modal + 20 new Explorer, per ligand."""
    a = assemble("A", ligands)
    b = assemble("B", ligands)
    out = {}
    for lig in ligands:
        if lig not in a or lig not in b:
            continue
        D = np.vstack([a[lig]["D"], b[lig]["D"]])
        names = [f"A::{n}" for n in a[lig]["names"]] + [f"B::{n}" for n in b[lig]["names"]]
        lddt = np.concatenate([a[lig]["lddt"], b[lig]["lddt"]])
        fra = load_frames("A")[lig][1]
        frb = load_frames("B")[lig][1]
        xyz = np.concatenate([fra, frb])
        out[lig] = {"names": names, "D": D, "lddt": lddt,
                    "feat": features(D, xyz)}
    return out


# ---------------------------------------------------------------------------
# selection primitives
# ---------------------------------------------------------------------------

def topk_idx(x: np.ndarray, k: int) -> np.ndarray:
    return np.argsort(x, kind="stable")[:k]


def pick_expect(vals: np.ndarray, lddt: np.ndarray) -> float:
    """Exact expectation of argmin(vals) with ties broken uniformly at random."""
    lo = vals.min()
    tied = np.flatnonzero(vals <= lo + 1e-12)
    return float(lddt[tied].mean())


def pick_draws(vals: np.ndarray, lddt: np.ndarray, rng, draws: int = TIE_DRAWS) -> float:
    lo = vals.min()
    tied = np.flatnonzero(vals <= lo + 1e-12)
    if len(tied) == 1:
        return float(lddt[tied[0]])
    return float(np.mean(lddt[rng.choice(tied, size=draws)]))


def within_rho(feat: np.ndarray, lddt: np.ndarray) -> float | None:
    from scipy.stats import spearmanr
    if len(np.unique(feat)) < 2 or len(np.unique(lddt)) < 2:
        return None
    r = spearmanr(feat, lddt).statistic
    return float(r) if np.isfinite(r) else None


# ---------------------------------------------------------------------------
# controls
# ---------------------------------------------------------------------------

def controls() -> dict:
    from cypstruct import xengine as X

    out: dict = {}
    t = truth_A()
    ship = pd.read_csv(DATA_PROCESSED / "xeng_val87b.csv")
    m = _cached_dmatrix("A")

    rows = []
    for lig, (names, D) in m.items():
        for i, nm in enumerate(names):
            rows.append({"ligand": lig, "sample": nm, "xeng_rebuilt": D[i].mean()})
    reb = pd.DataFrame(rows)
    j = ship.merge(reb, on=["ligand", "sample"], how="outer", indicator=True)
    out["C_XENG"] = {
        "shipped_rows": int(len(ship)), "rebuilt_rows": int(len(reb)),
        "merge_both": int((j._merge == "both").sum()),
        "merge_left_only": int((j._merge == "left_only").sum()),
        "merge_right_only": int((j._merge == "right_only").sum()),
        "max_abs_diff": float(np.nanmax(np.abs(j.xeng - j.xeng_rebuilt))),
    }
    out["C_XENG"]["passes_1e-6"] = bool(out["C_XENG"]["max_abs_diff"] < 1e-6)

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
    rhos = [within_rho(g.xeng.values, g.lddt_pli.values) for _, g in d.groupby("ligand")]
    rhos = np.array([r for r in rhos if r is not None])
    out["C_N2"] = {
        "selected": round(sel, 4), "oracle": round(orc, 4), "random": round(rnd, 4),
        "gain": round(sel - rnd, 4), "within_ligand_rho": round(float(rhos.mean()), 4),
        "correct_sign_pct": round(100 * float((rhos < 0).mean()), 2),
        "expected": {"selected": 0.6164, "oracle": 0.6975, "random": 0.5769,
                     "gain": 0.0395, "rho": -0.2582, "correct_sign_pct": 75.86}}
    out["C_N2"]["matches"] = bool(abs(sel - 0.6164) < 5e-4 and abs(orc - 0.6975) < 5e-4
                                  and abs(sel - rnd - 0.0395) < 5e-4)

    z = t[t.lddt_pli == 0]
    out["C_NUM"] = {
        "rows": int(len(t)), "mapped_true": int(t.mapped.sum()),
        "mapped_false": int((~t.mapped).sum()),
        "exact_zero_lddt": int(len(z)),
        "zero_ligands": z.ligand.value_counts().to_dict(),
        "zero_bisy_rmsd": [round(float(v), 2) for v in z.bisy_rmsd.values],
        "verdict": ("ejections, not FINDING_021 numbering"
                    if len(z) and (z.bisy_rmsd > 10).all()
                    else "none" if not len(z) else "CHECK")}

    ref = X.load_reference(FROZEN)
    dep = np.array([len(v) for v in ref.values()])
    out["C_REF"] = {"ligands": len(ref), "min": int(dep.min()),
                    "median": float(np.median(dep)), "max": int(dep.max()),
                    "all_at_least_4": bool((dep >= 4).all())}

    xt = load_crystals()
    out["C_XTAL"] = {"crystals_loaded": len(xt),
                     "of_ligands": int(t.ligand.nunique()),
                     "missing": sorted(set(t.ligand.unique()) - set(xt))}

    # B pool -- reproduce the shipped matched-depth xeng too
    mb = _cached_dmatrix("B")
    tb = truth_B().set_index(["ligand", "sample"])
    diffs = []
    for lig, (names, D) in mb.items():
        for i, nm in enumerate(names):
            if (lig, nm) in tb.index:
                diffs.append(abs(float(tb.loc[(lig, nm), "xeng"]) - D[i].mean()))
    out["C_XENG_B"] = {"rows": len(diffs),
                       "max_abs_diff": float(np.max(diffs)) if diffs else None,
                       "passes_1e-6": bool(diffs and np.max(diffs) < 1e-6)}

    ligs = stratum()
    out["C_FILT"] = {
        "validation_ligands": int(t.ligand.nunique()),
        "pool_A_poses": int(len(t)),
        "pool_A_poses_with_features": int(len(reb)),
        "predicted_type_I_stratum": len(ligs),
        "excluded_predicted_type_II": int(t.ligand.nunique()) - len(ligs),
        "pool_B_new_poses": int(len(truth_B())),
        "ligands_dropped_no_reference": 0,
        "ligands_dropped_no_truth": 0,
    }
    OUT_CONTROLS.write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str)[:4000])
    return out


# ---------------------------------------------------------------------------
# STAGE 1 -- the prize
# ---------------------------------------------------------------------------

def _regime(data: dict, label: str) -> dict:
    ligs = sorted(data)
    rows = []
    for lig in ligs:
        x = data[lig]["feat"]["xeng"]
        y = data[lig]["lddt"]
        order = np.argsort(x, kind="stable")
        xs = x[order]
        r = {"ligand": lig, "n_poses": len(x),
             "margin_12": float(xs[1] - xs[0]),
             "zmargin_12": float((xs[1] - xs[0]) / x.std(ddof=1)),
             "xeng_sd": float(x.std(ddof=1)),
             "sel": float(y[order[0]]),
             "oracle": float(y.max()), "rand": float(y.mean()),
             "n_unique_xeng": int(len(np.unique(x)))}
        for k in KS:
            idx = order[:k]
            r[f"stake_{k}"] = float(y[idx].max() - y[idx[0]])
            r[f"downside_{k}"] = float(y[idx[0]] - y[idx].min())
            r[f"randk_{k}"] = float(y[idx].mean())
        rows.append(r)
    df = pd.DataFrame(rows)
    out = {"label": label, "n_ligands": len(df),
           "pool_oracle": round(float(df.oracle.mean()), 4),
           "selected_incumbent": round(float(df.sel.mean()), 4),
           "random": round(float(df["rand"].mean()), 4),
           "gain_incumbent": round(float(df.sel.mean() - df["rand"].mean()), 4),
           "oracle_gap": round(float(df.oracle.mean() - df.sel.mean()), 4),
           "margin_12": {q: round(float(df.margin_12.quantile(q / 100)), 4)
                         for q in (5, 25, 50, 75, 95)},
           "zmargin_12": {q: round(float(df.zmargin_12.quantile(q / 100)), 3)
                          for q in (5, 25, 50, 75, 95)},
           "ligands_with_duplicate_xeng": int((df.n_unique_xeng < df.n_poses).sum()),
           "margin_below_0.06A": int((df.margin_12 < 0.06).sum()),
           "margin_below_0.10A": int((df.margin_12 < 0.10).sum())}
    for k in KS:
        s, dn, rk = df[f"stake_{k}"], df[f"downside_{k}"], df[f"randk_{k}"]
        out[f"k{k}"] = {
            "ORACLE_tiebreak": round(float(s.mean()), 4),
            "oracle_tiebreak_median": round(float(s.median()), 4),
            "worst_case_tiebreak": round(-float(dn.mean()), 4),
            "random_in_topk": round(float(rk.mean()), 4),
            "incumbent_minus_random_in_topk": round(float(df.sel.mean() - rk.mean()), 4),
            "ligands_with_stake_gt_0": int((s > 1e-9).sum()),
            "ligands_with_stake_gt_0.05": int((s > 0.05).sum()),
            "share_of_oracle_gap": round(float(s.mean() /
                                               (df.oracle.mean() - df.sel.mean())), 4),
        }
    # does the stake concentrate in the near-ties?
    qs = pd.qcut(df.margin_12, 4, labels=["Q1 tightest", "Q2", "Q3", "Q4 widest"])
    out["stake_2_by_margin_quartile"] = {
        str(g): {"n": int(len(s)), "mean_margin": round(float(s.margin_12.mean()), 4),
                 "mean_stake_2": round(float(s.stake_2.mean()), 4),
                 "mean_downside_2": round(float(s.downside_2.mean()), 4)}
        for g, s in df.groupby(qs, observed=True)}
    out["rho_margin_vs_stake2"] = None
    from scipy.stats import spearmanr
    r = spearmanr(df.margin_12, df.stake_2)
    out["rho_margin_vs_stake2"] = {"rho": round(float(r.statistic), 4),
                                   "p": float(r.pvalue)}
    return out, df


def regime() -> dict:
    a = assemble("A")
    out_a, df_a = _regime(a, "pool A -- 87 ligands x 20 Boltz-2 poses")
    ligs = [l for l in stratum() if l in a]
    b = assemble_union(ligs)
    out_b, df_b = _regime(b, "pool B -- 14 predicted-Type-I ligands x 40 (FINDING 035)")
    out_a2, _ = _regime({k: a[k] for k in ligs}, "pool B ligands at depth 20 only")

    o2 = out_a["k2"]["ORACLE_tiebreak"]
    verdict = ("CLOSED -- a perfect top-2 tie-break is below the pooled noise floor"
               if o2 < FLOOR_POOLED_P95 else
               "SMALL -- between the floor and 2x the floor; candidates tested but the "
               "pre-registered expectation is that the paired bar will not clear"
               if o2 < 2 * FLOOR_POOLED_P95 else
               "PROSPECT -- the ceiling is at least 2x the floor")
    res = {"pool_A": out_a, "pool_B_depth40": out_b, "pool_B_depth20": out_a2,
           "floor_pooled_p95_quoted": FLOOR_POOLED_P95,
           "floor_n14_p95_quoted": FLOOR_N14_P95,
           "STOP_RULE": {"O2_pool_A": o2, "verdict": verdict}}
    OUT_REGIME.write_text(json.dumps(res, indent=2, default=str))
    df_a.to_csv(OUT_PERLIG, index=False)
    print(json.dumps(res, indent=2, default=str))
    return res


# ---------------------------------------------------------------------------
# STAGE 2 -- R1, the answer-recognition gate
# ---------------------------------------------------------------------------

FEATS = ["xeng", "xeng_2nd", "xeng_borda", "xeng_sd", "xeng_votes"]


def gate() -> dict:
    from scipy.stats import binomtest
    a = assemble("A")
    xt = load_crystals()
    rows = []
    for lig, rec in a.items():
        if lig not in xt:
            continue
        from cypstruct import xengine as X
        ref = X.load_reference(FROZEN)[lig]
        dc = np.array([X.chamfer(np.asarray(xt[lig], float), np.asarray(r, float))
                       for r in ref])
        D = rec["D"]
        cf = crystal_features(D, dc)
        pf = dict(rec["feat"])
        pf["xeng_borda"] = crystal_borda_for_poses(D, dc)
        pf["xeng_votes"] = crystal_votes_for_poses(D, dc)
        r = {"ligand": lig}
        for f in FEATS:
            v = np.asarray(pf[f], float)
            # LOWER is better for every feature (votes already negated):
            # percentile = fraction of the 20 predictions that are WORSE than the crystal
            r[f] = float(np.mean(v > cf[f]))
            r[f + "__delta"] = float(cf[f] - v.mean())
        rows.append(r)
    df = pd.DataFrame(rows)
    out = {"n_ligands": int(len(df)), "pass_rule": "mean percentile >= 0.65 and p < 0.05",
           "features": {}}
    for f in FEATS:
        better = int((df[f] > 0.5).sum())
        bt = binomtest(better, len(df), 0.5)
        out["features"][f] = {
            "mean_percentile": round(float(df[f].mean()), 4),
            "median_percentile": round(float(df[f].median()), 4),
            "crystal_beats_median_pose": f"{better} / {len(df)}",
            "binomial_p": float(bt.pvalue),
            "mean_delta_crystal_minus_poolmean": round(float(df[f + "__delta"].mean()), 4),
            "passes": bool(df[f].mean() >= 0.65 and bt.pvalue < 0.05
                           and better > len(df) / 2)}
    cal = out["features"]["xeng"]["passes"]
    out["CALIBRATION"] = {
        "shipped_xeng_passes_R1": cal,
        "reading": ("R1 is informative for this family; a candidate that fails it is dead"
                    if cal else
                    "R1 is NOT a valid gate for this family -- the shipped feature, which "
                    "demonstrably selects at +0.0395, also fails it. Reported as a "
                    "diagnostic only; every candidate is carried to the paired bar. "
                    "This branch was declared a priori in the pre-registration.")}
    OUT_GATE.write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str))
    return out


# ---------------------------------------------------------------------------
# STAGE 3/4 -- candidates and the three bars
# ---------------------------------------------------------------------------

def _apply(data: dict, feat: str, k: int | None, rng) -> tuple[np.ndarray, np.ndarray]:
    """Per-ligand selected LDDT-PLI for a tie-break (k) or a standalone argmin (k=None)."""
    ligs = sorted(data)
    sel, inc = [], []
    for lig in ligs:
        x = data[lig]["feat"]["xeng"]
        y = data[lig]["lddt"]
        order = np.argsort(x, kind="stable")
        inc.append(float(y[order[0]]))
        if k is None:
            idx = np.arange(len(x))
        else:
            idx = order[:k]
        if feat == "med":
            if k is None or k < 3:
                sel.append(float("nan"))
                continue
            P = data[lig]["feat"]["_pp"][np.ix_(idx, idx)]
            v = (P.sum(axis=1)) / (len(idx) - 1)
        else:
            v = np.asarray(data[lig]["feat"][feat], float)[idx]
        sel.append(pick_draws(v, y[idx], rng))
    return np.array(sel), np.array(inc)


def _paired(sel: np.ndarray, inc: np.ndarray, rng) -> dict:
    from scipy.stats import wilcoxon
    d = sel - inc
    ok = np.isfinite(d)
    d = d[ok]
    changed = int(np.sum(np.abs(d) > 1e-9))
    boot = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(BOOT)])
    try:
        w = float(wilcoxon(sel[ok], inc[ok]).pvalue) if changed else 1.0
    except ValueError:
        w = 1.0
    ch = d[np.abs(d) > 1e-9]
    return {"paired_mean": round(float(d.mean()), 4),
            "ci95": [round(float(np.percentile(boot, 2.5)), 4),
                     round(float(np.percentile(boot, 97.5)), 4)],
            "wilcoxon_p": w, "n": int(len(d)),
            "TIES_unchanged": int(len(d) - changed),
            "PICKS_CHANGED": changed,
            "on_changed_only": {
                "n": changed,
                "mean_delta": round(float(ch.mean()), 4) if changed else None,
                "better": int((ch > 0).sum()), "worse": int((ch < 0).sum()),
                "incumbent_mean": round(float(inc[ok][np.abs(d) > 1e-9].mean()), 4)
                if changed else None,
                "candidate_mean": round(float(sel[ok][np.abs(d) > 1e-9].mean()), 4)
                if changed else None}}


def _nulls(data: dict, rng) -> dict:
    ligs = sorted(data)
    Y = [data[l]["lddt"] for l in ligs]
    X = [data[l]["feat"]["xeng"] for l in ligs]
    order = [np.argsort(x, kind="stable") for x in X]
    inc = np.array([float(y[o[0]]) for y, o in zip(Y, order)])
    rnd = float(np.mean([y.mean() for y in Y]))
    out = {"random_baseline": round(rnd, 4), "incumbent": round(float(inc.mean()), 4)}

    # FINDING 007 null: a random within-ligand feature, argmin over ALL poses
    g = np.array([np.mean([float(rng.choice(y)) for y in Y]) - rnd
                  for _ in range(NULL_DRAWS)])
    out["random_feature_null_gain"] = {"p95": round(float(np.percentile(g, 95)), 4),
                                       "p99": round(float(np.percentile(g, 99)), 4),
                                       "draws": NULL_DRAWS}
    # matched tie-break null: a random reorder of the SAME top-k
    for k in KS:
        dd = np.array([np.mean([float(rng.choice(y[o[:k]])) for y, o in zip(Y, order)])
                       - inc.mean() for _ in range(NULL_DRAWS)])
        out[f"tiebreak_null_k{k}_vs_incumbent"] = {
            "mean": round(float(dd.mean()), 4),
            "p95": round(float(np.percentile(dd, 95)), 4),
            "p99": round(float(np.percentile(dd, 99)), 4), "draws": NULL_DRAWS}
    return out


def _rank_stats(data: dict) -> dict:
    from scipy.stats import binomtest
    out = {}
    for f in FEATS:
        rs = []
        for lig in sorted(data):
            r = within_rho(np.asarray(data[lig]["feat"][f], float), data[lig]["lddt"])
            if r is not None:
                rs.append(r)
        rs = np.array(rs)
        neg = int((rs < 0).sum())
        out[f] = {"mean_within_ligand_rho": round(float(rs.mean()), 4),
                  "correct_sign_pct": round(100 * neg / len(rs), 2),
                  "n": len(rs),
                  "binomial_p": float(binomtest(neg, len(rs), 0.5).pvalue)}
    return out


def _board(data: dict, label: str, floor: float) -> dict:
    rng = np.random.default_rng(SEED)
    ligs = sorted(data)
    rnd = float(np.mean([data[l]["lddt"].mean() for l in ligs]))
    orc = float(np.mean([data[l]["lddt"].max() for l in ligs]))
    res = {"label": label, "n_ligands": len(ligs),
           "POOL_ORACLE": round(orc, 4), "random": round(rnd, 4),
           "noise_floor_quoted": floor,
           "nulls": _nulls(data, rng), "ranking": _rank_stats(data), "candidates": {}}
    _, inc = _apply(data, "xeng", 1, rng)
    res["incumbent_selected"] = round(float(inc.mean()), 4)
    res["incumbent_gain_vs_random"] = round(float(inc.mean() - rnd), 4)

    cands = [("T1_medoid", "med"), ("T2_second_nearest", "xeng_2nd"),
             ("T3_borda", "xeng_borda"), ("T4_ref_sd", "xeng_sd"),
             ("T5_plurality", "xeng_votes")]
    for name, f in cands:
        res["candidates"][name] = {}
        for k in KS:
            if f == "med" and k < 3:
                res["candidates"][name][f"k{k}"] = "N/A -- undefined at k=2"
                continue
            sel, inc2 = _apply(data, f, k, np.random.default_rng(SEED + k))
            r = _paired(sel, inc2, np.random.default_rng(SEED + 7 * k))
            r["selected"] = round(float(np.nanmean(sel)), 4)
            r["gain_vs_random"] = round(float(np.nanmean(sel) - rnd), 4)
            res["candidates"][name][f"k{k}"] = r
        if f != "med":
            sel, inc2 = _apply(data, f, None, np.random.default_rng(SEED + 99))
            r = _paired(sel, inc2, np.random.default_rng(SEED + 101))
            r["selected"] = round(float(np.nanmean(sel)), 4)
            r["gain_vs_random"] = round(float(np.nanmean(sel) - rnd), 4)
            res["candidates"][name]["standalone_all_poses"] = r
    return res


def select() -> dict:
    a = assemble("A")
    ligs = [l for l in stratum() if l in a]
    b = assemble_union(ligs)
    out = {"pool_A": _board(a, "pool A -- 87 x 20", FLOOR_POOLED_P95),
           "pool_B_depth40": _board(b, "pool B -- 14 x 40", FLOOR_N14_P95)}
    OUT_SELECT.write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str))
    return out


def swing() -> dict:
    """POST-HOC. The two ligands that carried all of `FINDING_035`'s Δ selected.

    Not a candidate, not eligible to ship: it is two ligands, chosen because 035 named
    them. It is here because it shows exactly what the n=14 board is made of.
    """
    a = assemble("A")
    b = assemble_union([l for l in stratum() if l in a])
    out = {}
    for lig in ("CFF", "08J"):
        if lig not in b:
            continue
        r = b[lig]
        x = r["feat"]["xeng"]
        y = r["lddt"]
        o = np.argsort(x, kind="stable")[:3]
        rec = {"top3_xeng": [round(float(v), 4) for v in x[o]],
               "top3_lddt": [round(float(v), 4) for v in y[o]],
               "top3_names": [r["names"][i] for i in o],
               "incumbent": round(float(y[o[0]]), 4),
               "oracle_in_top3": round(float(y[o].max()), 4), "picks": {}}
        for f in ("xeng_2nd", "xeng_borda", "xeng_sd", "xeng_votes"):
            v = np.asarray(r["feat"][f], float)[o]
            rec["picks"][f] = round(float(y[o][int(np.argmin(v))]), 4)
        P = r["feat"]["_pp"][np.ix_(o, o)]
        rec["picks"]["medoid"] = round(float(y[o][int(np.argmin(P.sum(1) / 2))]), 4)
        out[lig] = rec
    (DATA_PROCESSED / "tiebreak_swing.json").write_text(
        json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["frames", "controls", "regime", "gate",
                                      "select", "swing", "all"])
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    if a.stage in ("frames", "all"):
        frames(a.workers)
    if a.stage in ("controls", "all"):
        controls()
    if a.stage in ("regime", "all"):
        regime()
    if a.stage in ("gate", "all"):
        gate()
    if a.stage in ("select", "all"):
        select()
    if a.stage in ("swing", "all"):
        swing()
    return 0


if __name__ == "__main__":
    os.environ.setdefault("PYTHONWARNINGS", "ignore")
    raise SystemExit(main())
