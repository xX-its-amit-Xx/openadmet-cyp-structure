"""CONTRADICTION 2 - FINDING 025 arm4 holdout: does +0.0408 reproduce, and does the
   catastrophe mechanism survive the FINDING 021 numbering fix?"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon

REPO = Path(r"D:/Users/ashenoy00000/.windsurf/OpenADMET-cyp-structure")
DP = REPO / "data" / "processed"
OUT = {}

df = pd.read_csv(r"C:/Temp/cyp_corrections/depth_terms_holdout.csv")
print("FILTER: depth_terms_holdout.csv rows =", len(df), "cols", df.columns.tolist())
sc = json.loads((DP / "finetune" / "scores_arm4_mix_base.json").read_text())
print("FILTER: scores_arm4_mix_base.json pairs =", len(sc))

# ---- control: the lddt column in the feature file IS the post-fix score ----
mism = 0
for _, r in df.iterrows():
    v = sc[r["name"]]["lddt_pli"][int(r["rank"])]
    if abs(v - r["lddt_pli"]) > 1e-9:
        mism += 1
print("CONTROL C-LDDT: feature-file lddt vs post-fix scores json, mismatches =",
      mism, "/", len(df))
off = np.array([sc[k]["resnum_offset"] for k in sc])
ident = np.array([sc[k]["seq_identity_at_offset"] for k in sc])
print("CONTROL C-NUM (FINDING 021): pairs with offset != 0 =", int((off != 0).sum()),
      "; min seq identity at chosen offset =", round(float(ident.min()), 3),
      "; pairs below 0.80 identity =", int((ident < 0.80).sum()))


def evaluate(d, col, sign=-1, n_draws=64, rng=None):
    rng = rng if rng is not None else np.random.default_rng(23)
    groups = list(d.groupby("name"))
    rnd_lig = np.array([g["lddt_pli"].mean() for _, g in groups])
    per = np.zeros(len(groups))
    for _ in range(n_draws):
        for i, (_nm, g) in enumerate(groups):
            v = sign * g[col].astype(float).values
            v = np.where(np.isfinite(v), v, -np.inf)
            per[i] += g["lddt_pli"].values[rng.choice(np.flatnonzero(v == v.max()))]
    per /= n_draws
    rs = []
    for _nm, g in groups:
        v = sign * g[col].astype(float).values
        ok = np.isfinite(v)
        if ok.sum() < 3 or np.nanstd(v[ok]) == 0:
            continue
        r = spearmanr(v[ok], g["lddt_pli"].values[ok]).statistic
        if r == r:
            rs.append(r)
    rs = np.array(rs)
    draws = np.array([np.mean([rng.choice(g["lddt_pli"].values) for _n, g in groups])
                      for _ in range(2000)])
    return dict(ligands=len(groups),
                oracle=round(float(np.mean([g["lddt_pli"].max() for _, g in groups])), 4),
                random=round(float(rnd_lig.mean()), 4),
                selected=round(float(per.mean()), 4),
                gain=round(float(per.mean() - rnd_lig.mean()), 4),
                rho=round(float(rs.mean()), 3),
                frac_rho_pos=round(float((rs > 0).mean()), 3),
                beats_random_frac=round(float((per > rnd_lig).mean()), 3),
                p_wilcoxon=float(wilcoxon(per, rnd_lig).pvalue),
                null_p95=round(float(np.percentile(draws, 95) - rnd_lig.mean()), 4),
                p_permutation=round(float((draws >= per.mean()).mean()), 4),
                per_pair_gain={nm: float(p - r) for (nm, _), p, r
                               in zip(groups, per, rnd_lig)})


pub_all = json.loads((DP / "depth_replication_holdout.json").read_text())
rep = {}
for col in ("fg_depth", "fe_centroid_dist", "fe_min_dist"):
    r = evaluate(df, col, -1, rng=np.random.default_rng(23))
    rep[col] = r
    pub = pub_all[col + " (lower is better)"]
    print("REPRO %-17s gain %+.4f (published %+.4f)  rho %.3f (pub %.3f)  "
          "p_perm %.4f (pub %.4f)  random %.4f (pub %.4f)"
          % (col, r["gain"], pub["gain"], r["rho"], pub["rho"],
             r["p_permutation"], pub["p_permutation"], r["random"], pub["random"]))
OUT["reproduction"] = {k: {kk: vv for kk, vv in v.items() if kk != "per_pair_gain"}
                       for k, v in rep.items()}

# ---- the catastrophe census, post-fix, on THIS holdout ----
groups = {nm: g.sort_values("rank") for nm, g in df.groupby("name")}
L = {nm: g["lddt_pli"].values for nm, g in groups.items()}
pose_all = np.concatenate(list(L.values()))
cat_pose = pose_all < 0.1
pairs_any = [nm for nm, v in L.items() if (v < 0.1).any()]
pairs_all = [nm for nm, v in L.items() if (v < 0.1).all()]
pairs_mean = [nm for nm, v in L.items() if v.mean() < 0.1]
print()
print("CATASTROPHE CENSUS (arm4 holdout, post-FINDING-021):")
print("  pairs = %d, poses = %d" % (len(L), len(pose_all)))
print("  poses < 0.1             : %d/%d = %.1f%%"
      % (cat_pose.sum(), len(pose_all), 100 * cat_pose.mean()))
print("  pairs with ANY pose<0.1 : %d/%d" % (len(pairs_any), len(L)))
print("  pairs with ALL poses<0.1: %d/%d   (FINDING 021 definition)"
      % (len(pairs_all), len(L)))
print("  pairs with mean < 0.1   : %d/%d" % (len(pairs_mean), len(L)))
print("  FINDING 025 asserted 39/85; FINDING 021 re-measured 1/84 on the bonded variant")
OUT["catastrophe_census"] = dict(
    pairs=len(L), poses=int(len(pose_all)),
    poses_below_0_1=int(cat_pose.sum()),
    pairs_any_pose_below_0_1=len(pairs_any),
    pairs_all_poses_below_0_1=len(pairs_all),
    pairs_mean_below_0_1=len(pairs_mean),
    asserted_in_finding_025=39, denominator_in_finding_025=85)

# ---- DECISIVE TEST: the catastrophe-avoidance ceiling ----
def ceiling(thresh, draws=20000, seed=7):
    rng = np.random.default_rng(seed)
    base = np.mean([v.mean() for v in L.values()])
    g = np.empty(draws)
    for i in range(draws):
        picks = [rng.choice(v[v >= thresh]) if (v >= thresh).any() else rng.choice(v)
                 for v in L.values()]
        g[i] = np.mean(picks) - base
    return float(g.mean())


print()
print("CEILING of a PERFECT catastrophe-avoider (random among survivors):")
for t in (0.1, 0.2, 0.3, 0.5):
    c = ceiling(t)
    print("  avoid every pose < %.1f : %+.4f" % (t, c))
    OUT.setdefault("avoidance_ceiling", {})[str(t)] = round(c, 4)
print("  observed fe_centroid_dist gain: %+.4f" % rep["fe_centroid_dist"]["gain"])

# ---- decomposition of the observed gain ----
g = rep["fe_centroid_dist"]["per_pair_gain"]
names = list(L)
gv = np.array([g[nm] for nm in names])
has_cat = np.array([nm in pairs_any for nm in names])
print()
print("DECOMPOSITION of the fe_centroid_dist gain by pair:")
print("  pairs containing a catastrophic pose: n=%d  mean %+.4f  contributes %+.4f of the total"
      % (has_cat.sum(), gv[has_cat].mean() if has_cat.any() else 0,
         gv[has_cat].sum() / len(gv)))
print("  pairs with NO catastrophic pose     : n=%d  mean %+.4f  contributes %+.4f of the total"
      % ((~has_cat).sum(), gv[~has_cat].mean(), gv[~has_cat].sum() / len(gv)))
OUT["decomposition"] = dict(
    total_gain=round(float(gv.mean()), 4),
    n_pairs_with_catastrophe=int(has_cat.sum()),
    contribution_from_catastrophe_pairs=round(float(gv[has_cat].sum() / len(gv)), 4),
    n_pairs_clean=int((~has_cat).sum()),
    contribution_from_clean_pairs=round(float(gv[~has_cat].sum() / len(gv)), 4),
    mean_gain_catastrophe_pairs=round(float(gv[has_cat].mean()), 4) if has_cat.any() else None,
    mean_gain_clean_pairs=round(float(gv[~has_cat].mean()), 4))

df2 = df[df.lddt_pli >= 0.1]
print("FILTER: drop poses < 0.1 -> %d/%d rows, %d/%d pairs retained"
      % (len(df2), len(df), df2.name.nunique(), df.name.nunique()))
r2 = evaluate(df2, "fe_centroid_dist", -1, rng=np.random.default_rng(23))
print("  fe_centroid_dist on the CATASTROPHE-FREE holdout: gain %+.4f  rho %.3f  "
      "p_perm %.4f  random %.4f  oracle %.4f"
      % (r2["gain"], r2["rho"], r2["p_permutation"], r2["random"], r2["oracle"]))
OUT["catastrophe_free_rerun"] = {k: v for k, v in r2.items() if k != "per_pair_gain"}

# ---- the alternative account: headroom and feature spread ----
cy = pd.read_csv(r"C:/Temp/cyp_corrections/depth_terms_cyp3a4.csv")
pubc = json.loads((DP / "depth_replication_cyp3a4.json").read_text())[
    "fe_centroid_dist (lower is better)"]
hh = rep["fe_centroid_dist"]["oracle"] - rep["fe_centroid_dist"]["random"]
hc = pubc["oracle"] - pubc["random"]
print()
print("HEADROOM (oracle - random), the 'pool uncertainty' account:")
print("  arm4 holdout : %+.4f  (fe_centroid_dist gain %+.4f)" % (hh, rep["fe_centroid_dist"]["gain"]))
print("  CYP3A4 pool  : %+.4f  (fe_centroid_dist gain %+.4f)" % (hc, pubc["gain"]))
OUT["headroom"] = {"arm4_holdout": round(hh, 4), "cyp3a4_pool": round(hc, 4),
                   "cyp3a4_gain": pubc["gain"],
                   "arm4_gain": rep["fe_centroid_dist"]["gain"]}


def featstats(d, col="fe_centroid_dist"):
    within = [g[col].std() for _, g in d.groupby("name")]
    rr = [g[col].max() - g[col].min() for _, g in d.groupby("name")]
    return [round(float(np.mean(within)), 3), round(float(np.mean(rr)), 3)]


print("FEATURE SPREAD within a pair [sd, range] in A:")
print("  arm4 holdout:", featstats(df), "  CYP3A4 pool:", featstats(cy))
OUT["feature_spread"] = {"arm4_holdout": featstats(df), "cyp3a4_pool": featstats(cy)}

(DP / "corrections_c2_arm4.json").write_text(json.dumps(OUT, indent=2))
print("\nwrote data/processed/corrections_c2_arm4.json")
