"""Build P450 validation splits whose chemical novelty MATCHES the real challenge split.

The problem this fixes. Every generalisation number in this repo comes from
leave-one-TARGET-out on the P450 superfamily: hold out a protein, select on the rest. That
is a good test of whether a method transfers across proteins. It is **not** a test of
whether it transfers across the chemical gap the challenge will actually impose, because
the challenge holds out *compounds*, and how hard that is depends entirely on how novel the
test chemistry is relative to train.

If the challenge test set is analog-dense — most test compounds have a close neighbour in
train — our leave-one-target-out number is pessimistic. If it is scaffold-novel, our number
is optimistic and drop day will disappoint. Neither is knowable without measuring, and
nothing in this repo has measured it.

So: measure the real train→test novelty distribution (max Tanimoto from each test compound
to the whole train set), then construct P450 ligand splits that reproduce that distribution.
A method validated on a matched split is being asked the same question the challenge asks.

This is step 1 and 2 of the method-hopping ladder in docs/IDEAS.md: optimise on an
analogous target with matched similarity structure, then hop.

    python scripts/structure/similarity_matched_splits.py measure   # the real gap
    python scripts/structure/similarity_matched_splits.py build     # matched P450 splits
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

RAW = Path(r"D:\Users\ashenoy00000\.windsurf\OpenADMET-cyp-challenge\data\raw")
UNI = DATA_PROCESSED / "p450_universe"
OUT = DATA_PROCESSED / "similarity_splits"


def fps(smiles: list[str]):
    """Morgan/ECFP4 bit vectors. The standard choice, so the number is comparable to
    anything else anyone reports on this challenge."""
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    out, keep = [], []
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if m is None:
            continue
        out.append(gen.GetFingerprint(m))
        keep.append(i)
    return out, keep, DataStructs


def max_sim_to(query_fps, ref_fps, DataStructs) -> np.ndarray:
    """For each query, its similarity to the NEAREST reference compound.

    Nearest-neighbour rather than mean: what makes a test compound easy is having one close
    analogue in train, not being vaguely similar to everything. A mean would wash that out
    and report every compound as equally novel.
    """
    out = np.empty(len(query_fps))
    for i, q in enumerate(query_fps):
        sims = DataStructs.BulkTanimotoSimilarity(q, ref_fps)
        out[i] = max(sims) if sims else 0.0
    return out


def cmd_measure() -> dict:
    tr = pd.read_csv(RAW / "cyp-challenge-TRAIN_inhibition.csv")
    te = pd.read_csv(RAW / "cyp-challenge-TEST-BLINDED.csv")
    tr_fp, _, DS = fps(list(tr.SMILES))
    te_fp, _, _ = fps(list(te.SMILES))
    print(f"train {len(tr_fp)} fingerprints | test {len(te_fp)}")
    nn = max_sim_to(te_fp, tr_fp, DS)
    q = {f"p{p}": round(float(np.percentile(nn, p)), 4) for p in (5, 25, 50, 75, 95)}
    rep = {
        "n_train": len(tr_fp), "n_test": len(te_fp),
        "nn_tanimoto": q,
        "mean": round(float(nn.mean()), 4),
        "frac_analog_dense_ge_0.7": round(float((nn >= 0.7).mean()), 4),
        "frac_novel_lt_0.4": round(float((nn < 0.4).mean()), 4),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    np.save(OUT / "challenge_nn_sim.npy", nn)
    (OUT / "challenge_gap.json").write_text(json.dumps(rep, indent=1))
    return rep


def cmd_build(n_splits: int, seed: int) -> dict:
    """Construct P450 ligand splits reproducing the challenge's novelty distribution.

    Greedy: draw a target novelty from the challenge distribution, then pick the held-out
    ligand whose nearest-neighbour similarity to the current training pool is closest to it.
    Rebuilding the neighbour distances after each pick is what makes the resulting split
    actually match rather than approximately match.
    """
    target = np.load(OUT / "challenge_nn_sim.npy")
    sc = pd.read_csv(UNI / "p450_poses_scored.csv")
    cs = pd.read_csv(UNI / "p450_cofold_set.csv")
    cs["pair"] = cs.pdb + "_" + cs.id
    lig = (cs[cs.pair.isin(set(sc.pair))][["id", "smiles"]]
           .drop_duplicates("id").reset_index(drop=True))
    fp, keep, DS = fps(list(lig.smiles))
    lig = lig.iloc[keep].reset_index(drop=True)
    print(f"{len(lig)} distinct P450 ligands with scored poses")

    rng = np.random.default_rng(seed)
    splits = []
    n_test = max(10, int(round(len(lig) * 0.2)))
    for s in range(n_splits):
        wanted = rng.choice(target, size=n_test, replace=True)
        wanted.sort()
        avail = list(range(len(lig)))
        test_idx: list[int] = []
        # start from the whole set as "train", move ligands out one at a time
        for w in wanted:
            train_idx = [i for i in avail if i not in set(test_idx)]
            if len(train_idx) < 20:
                break
            cand = [i for i in train_idx]
            sims = np.array([
                max(DS.BulkTanimotoSimilarity(fp[i], [fp[j] for j in train_idx if j != i]) or [0.0])
                for i in cand])
            pick = cand[int(np.argmin(np.abs(sims - w)))]
            test_idx.append(pick)
        got = []
        tr_set = [i for i in range(len(lig)) if i not in set(test_idx)]
        for i in test_idx:
            got.append(max(DS.BulkTanimotoSimilarity(fp[i], [fp[j] for j in tr_set]) or [0.0]))
        got = np.array(got)
        splits.append({
            "split": s,
            "n_test": len(test_idx),
            "test_ligands": [str(lig.id.iloc[i]) for i in test_idx],
            "achieved_mean_nn": round(float(got.mean()), 4),
            "achieved_p50": round(float(np.median(got)), 4),
        })
        print(f"  split {s}: {len(test_idx)} held out, "
              f"mean NN {got.mean():.4f} (target {wanted.mean():.4f})")
    (OUT / "p450_matched_splits.json").write_text(json.dumps(splits, indent=1))
    return {"splits": len(splits), "n_ligands": len(lig),
            "target_mean": round(float(target.mean()), 4),
            "achieved_mean": round(float(np.mean([s["achieved_mean_nn"] for s in splits])), 4),
            "file": str((OUT / "p450_matched_splits.json").relative_to(REPO))}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["measure", "build"])
    ap.add_argument("--splits", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out = cmd_measure() if a.cmd == "measure" else cmd_build(a.splits, a.seed)
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
