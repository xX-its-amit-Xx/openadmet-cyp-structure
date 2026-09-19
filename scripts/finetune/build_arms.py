"""Build the four fine-tuning datasets, with splits designed not to lie.

Amit's four arms, mapped onto data we already have with crystal ground truth:

  1. PROMISCUOUS  - ligands that bind many different P450s, and the proteins that bind
                    many different ligands. Chemical+biological breadth.
  2. FAMILY       - the whole PF00067 P450 superfamily. 491 pairs / 87 proteins.
  3. CYP3A4 ONLY  - the target itself. Narrowest, most on-distribution.
  4. MIX          - stratified: all of CYP3A4, plus family coverage weighted toward
                    proteins whose pockets resemble CYP3A4, plus the promiscuous tail.

**Why fine-tuning is the one intervention our own findings do not argue against.** Every
diversity experiment here added oracle the selector could not reach: a union of engines
gave +0.0375 oracle / -0.0017 selected, a sampler sweep +0.0355 / -0.0038. Those add poses
to the *tail*. A fine-tune instead moves the whole distribution - it raises the MEAN, and a
better mean needs no selector at all. That is a different mechanism, not another helping of
the one that keeps failing.

**The split is where a fake win would come from.** The 367 ligands are not 367 independent
examples: 62 of 83 coordinated CYP3A4 ligands present a pyridine donor, and the
ritonavir-analog series alone contributes several near-identical entries. A random split
puts analogs on both sides and reports a large meaningless gain.

So we hold out whole **Murcko-scaffold clusters**, and then check the held-out set's
nearest-neighbour similarity against FINDING 017's measured challenge distribution
(median 0.587). A split that is far more novel than the real task is not conservative - it
is answering a different question and will reject a method that would have worked.

    python scripts/finetune/build_arms.py
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "structure"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

UNI = DATA_PROCESSED / "p450_universe"
OUT = DATA_PROCESSED / "finetune_arms"
CYP3A4_UNIPROT = "P08684"
CHALLENGE_MEDIAN_NN = 0.587       # FINDING 017


def scaffold(smiles: str) -> str:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if m is None:
        return ""
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
    except Exception:
        return ""


def load_pairs() -> pd.DataFrame:
    sc = pd.read_csv(UNI / "p450_poses_scored.csv")
    cs = pd.read_csv(UNI / "p450_cofold_set.csv")
    cs["pair"] = cs.pdb + "_" + cs.id
    df = cs[cs.pair.isin(set(sc.pair))].drop_duplicates("pair").copy()
    # carry the measured pose quality through, so an arm can be filtered on it later
    best = sc.groupby("pair").lddt_pli.max().rename("best_lddt")
    df = df.merge(best, left_on="pair", right_index=True, how="left")
    df["scaffold"] = [scaffold(s) for s in df.smiles]
    return df


def promiscuity(df: pd.DataFrame) -> pd.DataFrame:
    """How many distinct proteins does each ligand bind, and vice versa."""
    lig_targets = df.groupby("id").uniprot.nunique().rename("n_targets")
    prot_ligs = df.groupby("uniprot").id.nunique().rename("n_ligands")
    return df.merge(lig_targets, left_on="id", right_index=True) \
             .merge(prot_ligs, left_on="uniprot", right_index=True)


def holdout_by_scaffold(df: pd.DataFrame, frac: float, seed: int) -> tuple[list, list]:
    """Hold out WHOLE scaffold clusters, never individual ligands.

    Splitting on ligands would put ritonavir analogs on both sides. Splitting on scaffolds
    is the coarsest unit that keeps a chemical series intact, so a held-out series is
    genuinely unseen rather than half-memorised.
    """
    rng = np.random.default_rng(seed)
    groups = defaultdict(list)
    for r in df.itertuples():
        groups[r.scaffold or f"__singleton_{r.pair}"].append(r.pair)
    keys = list(groups)
    rng.shuffle(keys)
    want = int(round(len(df) * frac))
    test, n = [], 0
    for k in keys:
        if n >= want:
            break
        test.extend(groups[k])
        n += len(groups[k])
    test_set = set(test)
    train = [p for p in df.pair if p not in test_set]
    return train, sorted(test_set)


def novelty_check(df: pd.DataFrame, train: list, test: list) -> dict:
    """Is the held-out chemistry as novel as the real task, or harder for no reason?"""
    from similarity_matched_splits import fps, max_sim_to
    tr_smi = df[df.pair.isin(set(train))].smiles.tolist()
    te_smi = df[df.pair.isin(set(test))].smiles.tolist()
    tr_fp, _, DS = fps(tr_smi)
    te_fp, _, _ = fps(te_smi)
    if not tr_fp or not te_fp:
        return {}
    nn = max_sim_to(te_fp, tr_fp, DS)
    return {"median_nn": round(float(np.median(nn)), 4),
            "mean_nn": round(float(nn.mean()), 4),
            "vs_challenge_median": round(float(np.median(nn) - CHALLENGE_MEDIAN_NN), 4)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--holdout", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    df = promiscuity(load_pairs())
    print(f"{len(df)} pairs | {df.id.nunique()} ligands | {df.uniprot.nunique()} proteins "
          f"| {df.scaffold.nunique()} Murcko scaffolds\n")

    cyp3a4 = df[df.uniprot == CYP3A4_UNIPROT]
    # promiscuous: ligands seen on >1 protein, or proteins with a broad ligand set.
    # Both halves matter - a promiscuous LIGAND teaches cross-target binding modes, a
    # promiscuous PROTEIN teaches one pocket accommodating diverse chemistry.
    promisc = df[(df.n_targets >= 2) | (df.n_ligands >= df.n_ligands.quantile(0.75))]

    arms = {
        "arm1_promiscuous": promisc,
        "arm2_family": df,
        "arm3_cyp3a4_only": cyp3a4,
    }
    # arm 4: everything CYP3A4, plus the promiscuous tail, plus family coverage sampled to
    # spread over proteins rather than over rows - otherwise the few heavily-deposited
    # proteins dominate and "family" collapses to "a handful of well-studied enzymes".
    rng = np.random.default_rng(a.seed)
    per_protein = (df[~df.pair.isin(set(cyp3a4.pair) | set(promisc.pair))]
                   .groupby("uniprot", group_keys=False)
                   .apply(lambda g: g.sample(min(len(g), 4), random_state=a.seed)))
    arms["arm4_mix"] = pd.concat([cyp3a4, promisc, per_protein]).drop_duplicates("pair")

    OUT.mkdir(parents=True, exist_ok=True)
    report = {}
    for name, sub in arms.items():
        if len(sub) < 20:
            print(f"{name}: only {len(sub)} pairs - skipped")
            continue
        train, test = holdout_by_scaffold(sub, a.holdout, a.seed)
        nov = novelty_check(sub, train, test)
        rec = {
            "pairs": int(len(sub)),
            "ligands": int(sub.id.nunique()),
            "proteins": int(sub.uniprot.nunique()),
            "scaffolds": int(sub.scaffold.nunique()),
            "train": len(train), "test": len(test),
            "held_out_scaffold_clusters": int(
                sub[sub.pair.isin(set(test))].scaffold.nunique()),
            "novelty": nov,
            "mean_best_lddt": round(float(sub.best_lddt.mean()), 4),
        }
        report[name] = rec
        sub.assign(split=np.where(sub.pair.isin(set(test)), "test", "train")) \
           .to_csv(OUT / f"{name}.csv", index=False)
        print(f"{name:22s} {len(sub):4d} pairs  {sub.uniprot.nunique():3d} proteins  "
              f"{rec['held_out_scaffold_clusters']:3d} held-out clusters  "
              f"NN {nov.get('median_nn','?')} (challenge {CHALLENGE_MEDIAN_NN})")

    (OUT / "arms_report.json").write_text(json.dumps(report, indent=1))
    print(f"\nwritten to {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
