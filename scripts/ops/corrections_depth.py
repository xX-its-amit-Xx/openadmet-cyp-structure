"""Secondary (a) - FINDING 035's two 20-pose baselines, and the pool-30 > pool-40 peak.

Both curves in 035 are Monte-Carlo averages over 256 subsample draws. Both are exactly
computable in closed form, because argmin over a random subset has a combinatorial
distribution. Computing them exactly removes the Monte-Carlo error entirely and settles
whether the non-monotonicity is real.
"""
import json
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(r"D:/Users/ashenoy00000/.windsurf/OpenADMET-cyp-structure")
DP = REPO / "data" / "processed"
OUT = {}

p = pd.read_csv(DP / "poses_scored_val87b.csv")
x = pd.read_csv(DP / "xeng_val87b.csv")
pool = p[p.arm == "unsteered"].merge(x, on=["ligand", "sample"], how="left")
lab = pd.read_csv(DP / "binding_mode_labels_cyp3a4.csv")
S = sorted(lab[lab.pred_fe_donor_median > 2.6].id.tolist())
new = pd.read_csv(DP / "matched_depth_poses.csv")
dist = json.loads((DP / "matched_depth_distinct.json").read_text())
keep = {(l, s) for l, v in dist["keep"].items() for s in v}
new = new[[(r.ligand, r.sample) in keep for r in new.itertuples()]]
print("FILTER: stratum n=%d; existing %d rows; new distinct %d rows"
      % (len(S), len(pool[pool.ligand.isin(S)]), len(new)))

E = {l: pool[pool.ligand == l][["xeng", "lddt_pli"]].values for l in S}
N = {l: new[new.ligand == l][["xeng", "lddt_pli"]].values for l in S}
for l in S:
    assert len(E[l]) == 20 and len(N[l]) == 20, (l, len(E[l]), len(N[l]))


def exact_augment(k):
    """E[selected] when all 20 existing poses are kept and k of the 20 new are added."""
    tot = 0.0
    for l in S:
        e, n = E[l], N[l]
        i = int(np.argmin(e[:, 0]))
        estar_x, estar_y = e[i, 0], e[i, 1]
        if k == 0:
            tot += estar_y
            continue
        w = n[n[:, 0] < estar_x]
        w = w[np.argsort(w[:, 0], kind="stable")]
        m, M = len(w), len(n)
        denom = comb(M, k)
        exp = 0.0
        pnone = comb(M - m, k) / denom if M - m >= k else 0.0
        for j in range(1, m + 1):
            rest = M - j
            pj = comb(rest, k - 1) / denom if rest >= k - 1 else 0.0
            exp += pj * w[j - 1, 1]
        tot += exp + pnone * estar_y
    return tot / len(S)


def exact_union(d):
    """E[selected] when d poses are drawn from the 40-pose union."""
    tot = 0.0
    for l in S:
        u = np.vstack([E[l], N[l]])
        u = u[np.argsort(u[:, 0], kind="stable")]
        M = len(u)
        denom = comb(M, d)
        exp = 0.0
        for j in range(1, M + 1):
            rest = M - j
            pj = comb(rest, d - 1) / denom if rest >= d - 1 else 0.0
            exp += pj * u[j - 1, 1]
        tot += exp
    return tot / len(S)


def exact_union_oracle_random(d):
    orc = rnd = 0.0
    for l in S:
        u = np.vstack([E[l], N[l]])
        y = np.sort(u[:, 1])
        M = len(y)
        denom = comb(M, d)
        orc += sum(comb(j - 1, d - 1) / denom * y[j - 1] for j in range(d, M + 1))
        rnd += y.mean()
    return orc / len(S), rnd / len(S)


print("\nA. THE PRIMARY (augment) CURVE, exact vs 035's 256-draw Monte Carlo")
pub_aug = {20: 0.5270, 22: 0.5316, 24: 0.5354, 28: 0.5393, 30: 0.5411, 40: 0.5386}
rows = []
for k in range(0, 21):
    v = exact_augment(k)
    pubv = pub_aug.get(20 + k)
    rows.append((20 + k, v, pubv))
    if pubv is not None or k in (10, 12, 14, 16, 18):
        print("   pool %2d  exact %.4f%s" % (20 + k, v,
              "   (035 published %.4f, diff %+.4f)" % (pubv, v - pubv) if pubv else ""))
OUT["exact_augment_curve"] = {str(a): round(b, 5) for a, b, _ in rows}
vals = np.array([b for _, b, _ in rows])
print("   exact curve monotone non-decreasing in k: %s ; argmax at pool %d"
      % (bool(np.all(np.diff(vals) >= -1e-12)), 20 + int(np.argmax(vals))))
print("   exact pool 30 = %.4f, exact pool 40 = %.4f, difference %+.4f"
      % (vals[10], vals[20], vals[10] - vals[20]))
OUT["augment_peak"] = {"exact_pool30": round(float(vals[10]), 5),
                       "exact_pool40": round(float(vals[20]), 5),
                       "exact_diff": round(float(vals[10] - vals[20]), 5),
                       "published_pool30": 0.5411, "published_pool40": 0.5386,
                       "published_diff": 0.0025,
                       "exact_curve_is_monotone": bool(np.all(np.diff(vals) >= -1e-12))}

# Monte-Carlo error of the published pool-30 number (256 draws)
rng = np.random.default_rng(20260922)
draws = np.empty(4000)
for i in range(4000):
    tot = 0.0
    for l in S:
        idx = rng.choice(20, 10, replace=False)
        cand = np.vstack([E[l], N[l][idx]])
        tot += cand[int(np.argmin(cand[:, 0])), 1]
    draws[i] = tot / len(S)
se256 = draws.std(ddof=1) / np.sqrt(256)
print("   sd of one draw = %.4f -> SE of a 256-draw mean = %.4f" % (draws.std(ddof=1), se256))
print("   published 0.5411 is %+.2f SE from the exact value %.4f"
      % ((0.5411 - vals[10]) / se256, vals[10]))
OUT["augment_peak"]["se_of_256_draw_mean"] = round(float(se256), 5)
OUT["augment_peak"]["z_of_published_pool30"] = round(float((0.5411 - vals[10]) / se256), 2)

print("\nB. THE UNION SUBSAMPLE CURVE, exact")
pub_un = {1: 0.4597, 2: 0.4800, 5: 0.4906, 10: 0.5064, 20: 0.5291, 30: 0.5393, 40: 0.5386}
uvals = {}
for d in (1, 2, 5, 10, 20, 25, 30, 35, 38, 39, 40):
    v = exact_union(d)
    uvals[d] = v
    pubv = pub_un.get(d)
    print("   depth %2d  exact selected %.4f%s" % (d, v,
          "   (035 published %.4f, diff %+.4f)" % (pubv, v - pubv) if pubv else ""))
OUT["exact_union_curve"] = {str(k): round(v, 5) for k, v in uvals.items()}
seq = [exact_union(d) for d in range(1, 41)]
print("   exact union curve monotone non-decreasing: %s ; argmax at depth %d"
      % (bool(np.all(np.diff(seq) >= -1e-12)), 1 + int(np.argmax(seq))))
OUT["union_curve_monotone"] = bool(np.all(np.diff(seq) >= -1e-12))
OUT["union_argmax_depth"] = 1 + int(np.argmax(seq))

o20, r20 = exact_union_oracle_random(20)
o40, r40 = exact_union_oracle_random(40)
print("   exact union oracle: depth20 %.4f depth40 %.4f (035: 0.6264 / 0.6595)" % (o20, o40))
print("   exact union random: depth20 %.4f depth40 %.4f (035: 0.4587 / 0.4587)" % (r20, r40))

print("\nC. THE TWO 20-POSE BASELINES ARE DIFFERENT ESTIMANDS")
print("   0.5270 = argmin(xeng) over THE 20 Modal poses that exist (deterministic, no draw)")
print("   0.5291 = E[argmin(xeng) over 20 poses drawn from the 40-pose UNION] (exact %.4f)"
      % uvals[20])
print("   the second pool is half Explorer poses, so it is a different population entirely")
det20 = np.mean([E[l][int(np.argmin(E[l][:, 0])), 1] for l in S])
print("   deterministic 20-pose baseline recomputed here: %.4f" % det20)
OUT["two_baselines"] = {"deterministic_existing_20": round(float(det20), 5),
                        "expected_20_of_union_40_exact": round(float(uvals[20]), 5),
                        "published_deterministic": 0.5270,
                        "published_union_subsample": 0.5291}

(DP / "corrections_c3_depth.json").write_text(json.dumps(OUT, indent=2))
print("\nwrote data/processed/corrections_c3_depth.json")
