"""CONTRADICTION 1 - the authoritative n=14 noise floor, and why five values were quoted.

The statistic in every one of the five quotations is the SAME estimator:

    gain = mean_over_ligands( LDDT-PLI of one pose drawn uniformly at random from that
                              ligand's pool )  -  mean_over_ligands( that pool's mean )

which is what "a random feature's argmin" reduces to, and the floor is its p95 / p99 over
independent draws.  What differs between the quotations is the POPULATION it is evaluated
on, and the Monte-Carlo stream.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(r"D:/Users/ashenoy00000/.windsurf/OpenADMET-cyp-structure")
sys.path.insert(0, str(REPO / "src"))
DP = REPO / "data" / "processed"
COORD_MAX = 2.6
OUT = {}


# --------------------------------------------------------------------------- pools
def existing_pool():
    p = pd.read_csv(DP / "poses_scored_val87b.csv")
    x = pd.read_csv(DP / "xeng_val87b.csv")
    d = p[p.arm == "unsteered"].merge(x, on=["ligand", "sample"], how="left")
    assert d.xeng.notna().all(), "xeng did not merge onto every pool row"
    return d


def summarise(df, ligs):
    orc, sel, rnd = [], [], []
    for lid in ligs:
        s = df[df.ligand == lid]
        orc.append(s.lddt_pli.max())
        sel.append(float(s.lddt_pli.iloc[int(np.argmin(s.xeng.values))]))
        rnd.append(s.lddt_pli.mean())
    return dict(n=len(ligs), oracle=round(float(np.mean(orc)), 4),
                selected=round(float(np.mean(sel)), 4),
                random=round(float(np.mean(rnd)), 4),
                gain=round(float(np.mean(sel) - np.mean(rnd)), 4))


pool = existing_pool()
raw = pd.read_csv(DP / "poses_scored_val87b.csv")
all87 = sorted(pool.ligand.unique())
print("FILTER poses_scored_val87b: %d rows -> arm=='unsteered' %d rows, %d ligands, "
      "depth/ligand %s" % (len(raw), len(pool), len(all87),
                           sorted(pool.groupby("ligand").size().unique())))

lab = pd.read_csv(DP / "binding_mode_labels_cyp3a4.csv")
print("FILTER binding_mode_labels_cyp3a4: %d ligands" % len(lab))
S_PRED = sorted(lab[lab.pred_fe_donor_median > COORD_MAX].id.tolist())
S_XTAL = sorted(lab[lab.crystal_mode == "type_I"].id.tolist())
print("FILTER pred_fe_donor_median > %.1f  -> prediction-side Type I, n=%d" % (COORD_MAX, len(S_PRED)))
print("FILTER crystal_mode == 'type_I'    -> crystal-side   Type I, n=%d" % len(S_XTAL))
print("   membership difference: pred-only %s ; crystal-only %s"
      % (sorted(set(S_PRED) - set(S_XTAL)), sorted(set(S_XTAL) - set(S_PRED))))

# ---- CONTROL C-XENG: the shipped board reproduces ----
b = summarise(pool, all87)
ok = (b["selected"] == 0.6164 and b["oracle"] == 0.6975 and b["random"] == 0.5769
      and b["gain"] == 0.0395)
print("CONTROL C-XENG (all 87): %s  -> reproduces shipped board: %s" % (b, ok))
assert ok, "shipped board did not reproduce; stop"
OUT["control_C_XENG"] = dict(b, reproduces=bool(ok))

# ---- CONTROL C-NUM: FINDING 021's numbering control on this pool ----
if "mapped" in pool.columns:
    print("CONTROL C-NUM: pool rows with mapped==True: %d/%d ; zero-score rows: %d"
          % (int(pool.mapped.sum()) if pool.mapped.dtype != object else -1,
             len(pool), int((pool.lddt_pli == 0).sum())))
OUT["control_C_NUM"] = {"pool_rows_scoring_exactly_zero": int((pool.lddt_pli == 0).sum()),
                        "note": "FINDING 020's signature was 28/85 exact zeros; this pool has none"}

# ---- the augmented depth-40 pool (FINDING 035 / 036 pool B) ----
new = pd.read_csv(DP / "matched_depth_poses.csv")
dist = json.loads((DP / "matched_depth_distinct.json").read_text())
keep_pairs = {(lid, s) for lid, v in dist["keep"].items() for s in v}
nsel = new[[(r.ligand, r.sample) in keep_pairs for r in new.itertuples()]]
print("FILTER matched_depth_poses: %d rows -> distinct-kept %d rows (%d/ligand)"
      % (len(new), len(nsel), len(nsel) // 14))
pool_e14 = pool[pool.ligand.isin(S_PRED)]
aug = pd.concat([pool_e14, nsel], ignore_index=True)
print("FILTER augmented pool B: %d rows, depth/ligand %s"
      % (len(aug), sorted(aug.groupby("ligand").size().unique())))

for nm, d, ligs in (("pool B depth 20", pool_e14, S_PRED),
                    ("pool B depth 40", aug, S_PRED),
                    ("crystal-side Type I, depth 20", pool[pool.ligand.isin(S_XTAL)], S_XTAL)):
    print("CONTROL board %-30s %s" % (nm, summarise(d, ligs)))
    OUT.setdefault("control_boards", {})[nm] = summarise(d, ligs)


# --------------------------------------------------------------------------- the null
def floor(df, ligs, draws, seed):
    """p95/p99 of the random-pose-selection gain. Vectorised; identical statistic to
    FINDING 033/034/035/036."""
    rng = np.random.default_rng(seed)
    cols = [df[df.ligand == l].lddt_pli.values for l in ligs]
    depth = {len(c) for c in cols}
    assert len(depth) == 1, "ragged depth %s - the vectorised form assumes rectangular" % depth
    M = np.stack(cols)                       # (n_lig, n_pose)
    base = M.mean(axis=1).mean()
    n_lig, n_pose = M.shape
    g = np.empty(draws)
    CH = 200_000
    done = 0
    while done < draws:
        b = min(CH, draws - done)
        idx = rng.integers(0, n_pose, size=(b, n_lig))
        g[done:done + b] = M[np.arange(n_lig)[None, :], idx].mean(axis=1) - base
        done += b
    return g


def report(name, df, ligs, quoted, draws=2_000_000, seed=20260923):
    g = floor(df, ligs, draws, seed)
    p95, p99 = np.percentile(g, [95, 99])
    # bootstrap SE of the p95 estimate AT 4,000 DRAWS -- the spread to be explained
    rng = np.random.default_rng(11)
    sub = np.array([np.percentile(rng.choice(g, 4000, replace=True), 95) for _ in range(2000)])
    se4k = float(sub.std(ddof=1))
    # and the SE of the authoritative estimate itself
    subA = np.array([np.percentile(rng.choice(g, draws, replace=True), 95) for _ in range(200)])
    seA = float(subA.std(ddof=1))
    print("\n%s  (n=%d, depth=%d, %d draws)" % (name, len(ligs), len(df) // len(ligs), draws))
    print("   AUTHORITATIVE  p95 = %+.5f  (SE %.5f)   p99 = %+.5f" % (p95, seA, p99))
    print("   sd of the gain = %.5f ; a 4,000-draw p95 has SE %.5f, "
          "so a +-%.4f spread between two 4,000-draw runs is 1 sd of the estimator"
          % (g.std(ddof=1), se4k, se4k * np.sqrt(2)))
    for q in quoted:
        z = (q["value"] - p95) / se4k
        print("     quoted %+.4f  in %-12s  -> %+.2f sd of a 4,000-draw estimate"
              % (q["value"], q["where"], z))
    OUT.setdefault("floors", {})[name] = dict(
        n=len(ligs), depth=len(df) // len(ligs), draws=draws,
        p95=round(float(p95), 5), p99=round(float(p99), 5),
        se_of_authoritative_p95=round(seA, 5),
        se_of_a_4000_draw_p95=round(se4k, 5),
        sd_of_gain=round(float(g.std(ddof=1)), 5),
        quoted=quoted, ligands=ligs)
    return p95, se4k


print("\n" + "=" * 78)
print("THE THREE POPULATIONS THE FIVE NUMBERS BELONG TO")
print("=" * 78)

p_d1, se1 = report(
    "D1 - crystal-side Type I (FINDING 033), pool A depth 20",
    pool[pool.ligand.isin(S_XTAL)], S_XTAL,
    [{"value": 0.0433, "where": "FINDING 033"}])

p_d2, se2 = report(
    "D2 - prediction-side Type I, pool A depth 20",
    pool_e14, S_PRED,
    [{"value": 0.0431, "where": "FINDING 034"},
     {"value": 0.0432, "where": "034 curve"},
     {"value": 0.0453, "where": "FINDING 035"}])

p_d3, se3 = report(
    "D3 - prediction-side Type I, augmented pool B depth 40",
    aug, S_PRED,
    [{"value": 0.0456, "where": "FINDING 035"},
     {"value": 0.0435, "where": "FINDING 036"}])

# --------------------------------------------------------- the sign that was flipped
print("\n" + "=" * 78)
print("THE CONCLUSION THAT TURNS ON THIS: FINDING 036 SS4c")
print("=" * 78)
reg = json.loads((DP / "tiebreak_regime.json").read_text())
ceil40 = reg["pool_B_depth40"]["k2"]["ORACLE_tiebreak"]
ceil20 = reg["pool_B_depth20"]["k2"]["ORACLE_tiebreak"]
print("  pool B depth 40, PERFECT top-2 tie-break ceiling = %+.4f" % ceil40)
print("  authoritative D3 floor                          = %+.5f  "
      "(its own SE 0.00004; a 4,000-draw estimate of it has SE %.5f)" % (p_d3, se3))
print("  -> ceiling %s floor by %+.5f, i.e. %.2f sd of a 4,000-draw floor estimate"
      % ("ABOVE" if ceil40 > p_d3 else "BELOW", ceil40 - p_d3, abs(ceil40 - p_d3) / se3))
# exact probability that a single random-pose draw beats the observed ceiling
g3 = floor(aug, S_PRED, 2_000_000, 20260923)
pval = float((g3 >= ceil40).mean())
print("  exact one-sided p of the ceiling against its own null: %.4f (2e6 draws)" % pval)
print("  pool B depth 20 ceiling %+.4f vs D2 floor %+.5f -> %s"
      % (ceil20, p_d2, "ABOVE" if ceil20 > p_d2 else "BELOW"))
g2 = floor(pool_e14, S_PRED, 2_000_000, 20260923)
print("  exact one-sided p of the depth-20 ceiling: %.4f" % float((g2 >= ceil20).mean()))
OUT["sign_test"] = {
    "pool_B_depth40_perfect_top2_ceiling": ceil40,
    "authoritative_floor_D3_p95": round(float(p_d3), 5),
    "ceiling_minus_floor": round(float(ceil40 - p_d3), 5),
    "exact_p_of_ceiling_vs_its_own_null": round(pval, 5),
    "pool_B_depth20_ceiling": ceil20,
    "authoritative_floor_D2_p95": round(float(p_d2), 5),
    "exact_p_depth20": round(float((g2 >= ceil20).mean()), 5)}

# other candidates in 036 that were judged against a n=14 floor
print("\n  other FINDING 036 numbers judged against an n=14 floor:")
for label, val, flo in (("036 depth-40 k=3 ceiling", reg["pool_B_depth40"]["k3"]["ORACLE_tiebreak"], p_d3),
                        ("036 depth-40 k=5 ceiling", reg["pool_B_depth40"]["k5"]["ORACLE_tiebreak"], p_d3),
                        ("036 best pool-B candidate +0.0392", 0.0392, p_d3)):
    print("    %-38s %+.4f vs %+.5f -> %s" % (label, val, flo, "ABOVE" if val > flo else "BELOW"))

# --------------------------------------------------- the four gain-vs-floor verdicts
print("\n" + "=" * 78)
print("EXACT p-VALUES FOR THE PUBLISHED FLOOR-CLEARING VERDICTS (2e6 draws each)")
print("=" * 78)
pv = {}
for nm, dfx, ligs, gain, pubp in (
        ("FINDING 033 crystal-side Type I, depth 20", pool[pool.ligand.isin(S_XTAL)], S_XTAL, 0.0411, 0.0595),
        ("FINDING 034 prediction-side Type I, depth 20", pool_e14, S_PRED, 0.0673, 0.0063),
        ("FINDING 035 prediction-side Type I, depth 40", aug, S_PRED, 0.0799, None)):
    g = floor(dfx, ligs, 2_000_000, 20260923)
    ep = float((g >= gain).mean())
    print("  %-44s gain %+.4f  exact p = %.4f  (published %s)"
          % (nm, gain, ep, pubp))
    pv[nm] = {"gain": gain, "exact_p": round(ep, 5), "published_p": pubp,
              "p95": round(float(np.percentile(g, 95)), 5),
              "p99": round(float(np.percentile(g, 99)), 5)}
pv["FINDING 036 pool-B perfect top-2 ceiling, depth 40"] = {
    "gain": ceil40, "exact_p": round(pval, 5), "published_p": None,
    "p95": round(float(p_d3), 5)}
(DP / "corrections_c1_pvalues.json").write_text(json.dumps(pv, indent=2))

(DP / "corrections_c1_floor.json").write_text(json.dumps(OUT, indent=2))
print("\nwrote data/processed/corrections_c1_floor.json and corrections_c1_pvalues.json")
