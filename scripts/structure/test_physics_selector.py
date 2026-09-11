"""Can the tier-0 physics terms select a better pose than chance?

The bar is unusually low and that is the point. FINDING 001 measured Boltz's own
confidence selecting *worse* than random within a ligand, so the incumbent to beat is
not a strong learned selector — it is a coin flip.

What is tested, in increasing order of how much it could fool us:

1. **Single terms**, used directly as a ranking. No fitting, no parameters, nothing to
   overfit. If a raw geometric quantity beats random, that is close to unimpeachable.
2. **A fitted ranker** over all tier-0 terms, validated by **leave-one-scaffold-cluster-out**.
   Never random CV: the 87 ligands cluster hard by chemical series (the ritonavir analogs
   alone contribute several near-identical entries), and a random split puts analogs on
   both sides and reports a large meaningless gain.

Reported for every selector: mean LDDT-PLI of the pose it picks, against the random
baseline and the oracle, plus the fraction of ligands where it beats random. The oracle
is printed first, always, because a selection number without its ceiling cannot be read.

    python scripts/structure/test_physics_selector.py --tag val87b
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

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

# Tier-0 terms and the direction that SHOULD be better, from the physics.
# Sign is asserted from chemistry up front rather than fitted, so that a term "working"
# with the wrong sign is visible as the red flag it is.
SINGLE_TERMS: dict[str, int] = {
    "is_coordinated": +1,        # coordinating the iron is the native mode for most ligands
    "n_contacts": +1,            # more pocket engagement
    "max_clash": -1,             # steric overlap is bad
    "frac_proximal": -1,         # ligand density on the buried heme face is impossible
    "s_fe_donor_angle": +1,      # trans to the thiolate; crystal median 171 deg
    "fe_donor_dist": -1,         # nearer the iron, within reason
}


def scaffold_clusters(ligands: pd.DataFrame) -> dict[str, int]:
    """Cluster ligands by Murcko scaffold, then by Tanimoto at 0.6.

    Cluster-held-out splitting is not a refinement here, it is the difference between a
    real number and a fake one. These 87 ligands are not 87 independent examples.
    """
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    ids, fps, scaf = [], [], {}
    for r in ligands.itertuples():
        m = Chem.MolFromSmiles(r.smiles)
        if m is None:
            continue
        ids.append(r.id)
        fps.append(gen.GetFingerprint(m))
        try:
            scaf[r.id] = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
        except Exception:
            scaf[r.id] = ""

    cluster: dict[str, int] = {}
    nxt = 0
    for i, lid in enumerate(ids):
        if lid in cluster:
            continue
        cluster[lid] = nxt
        for j in range(i + 1, len(ids)):
            other = ids[j]
            if other in cluster:
                continue
            same_scaf = scaf.get(lid) and scaf.get(lid) == scaf.get(other)
            tani = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            if same_scaf or tani >= 0.6:
                cluster[other] = nxt
        nxt += 1
    return cluster


def evaluate(df: pd.DataFrame, score: pd.Series, label: str) -> dict:
    """Pick the top-scoring pose per ligand and report what it is worth."""
    d = df.assign(_s=score).dropna(subset=["_s", "lddt_pli"])
    picked, randmean, oracle = [], [], []
    for _lig, g in d.groupby("ligand"):
        if len(g) < 2:
            continue
        picked.append(g.loc[g._s.idxmax(), "lddt_pli"])
        randmean.append(g.lddt_pli.mean())
        oracle.append(g.lddt_pli.max())
    if not picked:
        return {"selector": label, "n": 0}
    picked, randmean, oracle = map(np.array, (picked, randmean, oracle))
    beat = float((picked > randmean).mean())
    # paired test against the per-ligand random expectation
    from scipy import stats
    t = stats.wilcoxon(picked, randmean) if len(picked) > 10 else None
    return {
        "selector": label, "n": int(len(picked)),
        "selected": round(float(picked.mean()), 4),
        "random": round(float(randmean.mean()), 4),
        "oracle": round(float(oracle.mean()), 4),
        "delta_vs_random": round(float((picked - randmean).mean()), 4),
        "frac_ligands_beating_random": round(beat, 3),
        "wilcoxon_p": (round(float(t.pvalue), 5) if t is not None else None),
        "oracle_captured": round(float((picked.mean() - randmean.mean())
                                       / max(1e-9, oracle.mean() - randmean.mean())), 3),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="val87b")
    ap.add_argument("--arm", default="unsteered",
                    help="steered arm is deprecated; see FINDING 001")
    a = ap.parse_args()

    df = pd.read_csv(DATA_PROCESSED / f"poses_scored_{a.tag}.csv")
    df = df[df.arm == a.arm].copy()
    ligands = pd.read_csv(DATA_PROCESSED / "validation_ligands.csv")
    print(f"arm={a.arm}  poses={len(df)}  ligands={df.ligand.nunique()}", flush=True)

    results = [evaluate(df, df.confidence, "boltz_confidence (incumbent)")]

    for term, sign in SINGLE_TERMS.items():
        if term not in df.columns:
            continue
        col = df[term].astype(float) * sign
        results.append(evaluate(df, col, f"single: {term} ({'+' if sign > 0 else '-'})"))

    # --- fitted ranker, leave-one-cluster-out ------------------------------
    clusters = scaffold_clusters(ligands)
    df["cluster"] = df.ligand.map(clusters)
    feats = [c for c in ("is_coordinated", "n_contacts", "max_clash", "frac_proximal",
                         "s_fe_donor_angle", "fe_donor_dist", "confidence")
             if c in df.columns]
    work = df.dropna(subset=feats + ["lddt_pli", "cluster"]).copy()
    work["is_coordinated"] = work["is_coordinated"].astype(float)

    n_clusters = int(work.cluster.nunique())
    print(f"scaffold clusters: {n_clusters} over {work.ligand.nunique()} ligands", flush=True)

    try:
        import lightgbm as lgb

        preds = np.full(len(work), np.nan)
        for cl in sorted(work.cluster.unique()):
            tr = work[work.cluster != cl]
            te_idx = np.where((work.cluster == cl).values)[0]
            if len(tr) < 200 or not len(te_idx):
                continue
            grp = tr.groupby("ligand", sort=False).size().values
            # rank within ligand: we only ever need the argmax, never a calibrated value
            rk = lgb.LGBMRanker(n_estimators=250, learning_rate=0.05, num_leaves=15,
                                min_child_samples=30, verbose=-1, random_state=0)
            # LambdaRank needs integer relevance; bucket LDDT-PLI into 5 grades per ligand
            y = (tr.groupby("ligand").lddt_pli
                   .transform(lambda s: pd.qcut(s.rank(method="first"), 5,
                                                labels=False, duplicates="drop"))
                   .fillna(0).astype(int))
            rk.fit(tr[feats], y, group=grp)
            preds[te_idx] = rk.predict(work.iloc[te_idx][feats])
        results.append(evaluate(work, pd.Series(preds, index=work.index),
                                "fitted ranker (leave-one-cluster-out)"))
    except ImportError:
        print("lightgbm unavailable; skipping the fitted ranker", flush=True)

    print(f"\n{'selector':44s} {'sel':>7s} {'rand':>7s} {'oracle':>7s} "
          f"{'delta':>7s} {'beat%':>6s} {'p':>8s} {'capt':>6s}")
    print("-" * 100)
    for r in sorted(results, key=lambda x: -(x.get("delta_vs_random") or -9)):
        if not r.get("n"):
            continue
        print(f"{r['selector']:44s} {r['selected']:7.4f} {r['random']:7.4f} "
              f"{r['oracle']:7.4f} {r['delta_vs_random']:+7.4f} "
              f"{r['frac_ligands_beating_random']*100:5.1f}% "
              f"{str(r['wilcoxon_p']):>8s} {r['oracle_captured']:6.3f}")

    out = DATA_PROCESSED / f"physics_selector_{a.tag}_{a.arm}.json"
    out.write_text(json.dumps({"arm": a.arm, "n_clusters": n_clusters,
                               "results": results}, indent=1))
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
