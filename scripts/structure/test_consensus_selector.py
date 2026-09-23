"""Does cross-sample consensus select a better pose than chance?

Motivation. FINDING 001 showed Boltz's confidence ranks poses worse than random, and the
tier-0 geometry test showed no single geometric term beats random at significance either.
Two terms were significantly *worse* than random - picking the most perfectly trans
Fe-donor angle, or the least clashing pose, actively selects a worse structure.

That pattern is informative rather than just disappointing. It says the coordination
anchor is **necessary but not sufficient**: once the donor nitrogen is pinned at 2.2 A,
the molecule can still rotate about the Fe-donor axis, and LDDT-PLI scores the contacts
of the *whole* ligand. A pose can have textbook metal geometry and have everything else
pointing the wrong way. Optimising the anchor alone therefore selects for idealised
geometry rather than for a correct pose.

Consensus does not have that failure mode. It asks a different question - not "is this
pose good by some model of goodness" but "do independent samples agree on it" - and it
is the standard fallback when a model's own confidence is uninformative.

Three consensus selectors, all parameter-free:

  - **medoid**: the pose with the lowest mean RMSD to all other samples of that ligand.
  - **largest-cluster medoid**: cluster samples at 2 A, take the medoid of the biggest
    cluster. Differs from the plain medoid when the pool is genuinely bimodal - which
    OpenADMET says is real here, reporting density consistent with multiple mutually
    exclusive conformations.
  - **cluster size of the pose's own cluster**: a per-pose feature, for the ranker.

Also records, per ligand, the pool's **spread** (mean pairwise RMSD) and **cluster count**
as an uncertainty estimate, and checks whether spread predicts which ligands are hard.

    python scripts/structure/test_consensus_selector.py --tag val87b --arm unsteered
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

# `modal` is imported lazily, inside the branch that needs it (G3). At module scope it
# made this file unimportable on any box without the package - and on drop day Modal is
# over its spend cap, so the Modal branch is the one that will NOT be taken.
from cypstruct import pose as P  # noqa: E402
from cypstruct.paths import DATA_PROCESSED, free_gb, safe_workers  # noqa: E402


def scratch_root() -> Path:
    d = Path("D:/cyp_scratch")
    try:
        d.mkdir(parents=True, exist_ok=True)
        if free_gb(d) >= 5.0:
            return d
    except OSError:
        pass
    c = Path("C:/cyp_scratch")
    c.mkdir(parents=True, exist_ok=True)
    return c


def _wanted(csv: str) -> set:
    d = pd.read_csv(csv)
    idc = "structure" if "structure" in d.columns and "id" not in d.columns else "id"
    return {str(v) for v in d[idc]}


def _local_units(pool: Path) -> list:
    dirs = [d.name for d in sorted(pool.glob("*")) if d.is_dir() and any(d.glob("*.cif"))]
    if dirs:
        return dirs
    return sorted({f.name.split("__")[0] for f in pool.glob("*.cif")})


def _local_files(pool: Path, unit: str) -> list:
    d = pool / unit
    if d.is_dir():
        return sorted(d.glob("*.cif"))
    lig = unit.split("__")[0]
    d = pool / lig
    if d.is_dir():
        return sorted(d.glob("*.cif"))
    return sorted(pool.glob(f"{lig}__*.cif"))


def symmetric_rmsd(xa: np.ndarray, ea: list[str], xb: np.ndarray, eb: list[str]) -> float:
    """Element-constrained optimal assignment RMSD between two copies of one molecule.

    A lower bound on the true symmetry-corrected RMSD, which is fine here: we only need
    a consistent pairwise distance to cluster on, not an absolute value to report.
    """
    if len(xa) != len(xb):
        return float("nan")
    from scipy.optimize import linear_sum_assignment

    d = np.linalg.norm(xa[:, None, :] - xb[None, :, :], axis=2)
    cost = d.copy()
    cost[np.array(ea)[:, None] != np.array(eb)[None, :]] = d.max() + 1e3
    r, c = linear_sum_assignment(cost)
    return float(np.sqrt((d[r, c] ** 2).mean()))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="val87b")
    ap.add_argument("--arm", default="unsteered")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--cluster-cutoff", type=float, default=2.0)
    ap.add_argument("--pool-dir", default=None,
                    help="read poses from a LOCAL directory instead of the Modal volume")
    ap.add_argument("--ligands", default=None,
                    help="ligand csv; with --pool-dir it replaces poses_scored_<tag>")
    ap.add_argument("--no-eval", action="store_true",
                    help="write consensus_features_<tag>_<arm>.csv and stop; the "
                         "evaluation needs crystals")
    a = ap.parse_args()

    sc_path = DATA_PROCESSED / f"poses_scored_{a.tag}.csv"
    blind = a.no_eval or not sc_path.exists()
    local = Path(a.pool_dir) if a.pool_dir else None
    if blind:
        truth = {}
        if a.ligands:
            ligands = sorted(_wanted(a.ligands))
        elif local is not None:
            ligands = sorted({u.split("__")[0] for u in _local_units(local)})
        else:
            raise SystemExit("blind run needs --pool-dir and/or --ligands")
    else:
        scored = pd.read_csv(sc_path)
        scored = scored[scored.arm == a.arm]
        truth = {(r.ligand, r.sample): r.lddt_pli for r in scored.itertuples()}
        ligands = sorted(scored.ligand.unique())
    print(f"arm={a.arm}  ligands={len(ligands)}  blind={blind}", flush=True)

    vol = None
    if local is None:
        import modal
        vol = modal.Volume.from_name("cyp-pool")
    root = Path(tempfile.mkdtemp(prefix="cypcons_", dir=str(scratch_root())))
    nw = safe_workers(a.workers, ram_per_worker_gb=1.0)
    print(f"using {nw} worker(s)", flush=True)
    lock = threading.Lock()

    per_ligand: dict[str, dict] = {}
    feature_rows: list[dict] = []

    def do_ligand(lig: str) -> tuple[str, dict, list[dict]]:
        job = f"{lig}__{a.arm}__s1"
        work = root / job
        work.mkdir(parents=True, exist_ok=True)
        try:
            if vol is None:
                files = _local_files(local, job)
            else:
                for e in vol.iterdir(f"/{a.tag}/{job}"):
                    fn = e.path.split("/")[-1]
                    if fn.endswith(".cif"):
                        (work / fn).write_bytes(b"".join(vol.read_file(e.path)))
                files = sorted(work.glob("*.cif"))

            names, xyz, elems = [], [], None
            for cif in files:
                try:
                    m = P.load_structure(cif)
                except Exception:
                    continue
                if len(m.lig_xyz) == 0:
                    continue
                names.append(cif.stem)
                xyz.append(m.lig_xyz)
                elems = [e.upper() for e in m.lig_elem]
            if len(names) < 3:
                return lig, {}, []

            n = len(names)
            D = np.zeros((n, n))
            for i in range(n):
                for j in range(i + 1, n):
                    D[i, j] = D[j, i] = symmetric_rmsd(xyz[i], elems, xyz[j], elems)

            mean_to_others = D.sum(1) / max(1, n - 1)
            medoid = int(np.argmin(mean_to_others))

            # single-linkage clustering at the cutoff
            assigned = [-1] * n
            clusters: list[list[int]] = []
            for i in range(n):
                if assigned[i] >= 0:
                    continue
                cid = len(clusters)
                stack, members = [i], []
                while stack:
                    k = stack.pop()
                    if assigned[k] >= 0:
                        continue
                    assigned[k] = cid
                    members.append(k)
                    stack.extend(int(j) for j in np.where(D[k] <= a.cluster_cutoff)[0]
                                 if assigned[j] < 0)
                clusters.append(members)
            biggest = max(clusters, key=len)
            sub = D[np.ix_(biggest, biggest)]
            big_medoid = biggest[int(np.argmin(sub.sum(1)))]

            info = {
                "n_samples": n,
                "medoid_sample": names[medoid],
                "largest_cluster_medoid_sample": names[big_medoid],
                "n_clusters": len(clusters),
                "largest_cluster_frac": round(len(biggest) / n, 3),
                "mean_pairwise_rmsd": round(float(D[np.triu_indices(n, 1)].mean()), 3),
            }
            rows = [{"ligand": lig, "sample": names[i],
                     "mean_rmsd_to_others": float(mean_to_others[i]),
                     "cluster_size": len(clusters[assigned[i]]),
                     "in_largest_cluster": assigned[i] == assigned[big_medoid]}
                    for i in range(n)]
            return lig, info, rows
        except Exception:
            return lig, {}, []
        finally:
            shutil.rmtree(work, ignore_errors=True)

    try:
        with ThreadPoolExecutor(max_workers=nw) as ex:
            futs = {ex.submit(do_ligand, l): l for l in ligands}
            for k, fut in enumerate(as_completed(futs), 1):
                lig, info, rows = fut.result()
                if info:
                    with lock:
                        per_ligand[lig] = info
                        feature_rows.extend(rows)
                if k % 20 == 0:
                    print(f"  {k}/{len(ligands)} ligands", flush=True)
    finally:
        shutil.rmtree(root, ignore_errors=True)

    feats = pd.DataFrame(feature_rows)
    out_feats = DATA_PROCESSED / f"consensus_features_{a.tag}_{a.arm}.csv"
    feats.to_csv(out_feats, index=False)

    if blind:
        # The features ARE the deliverable on a blind set: build_submission's fallback
        # rule reads this csv plus orientation_features, and neither needs a crystal.
        # Everything below scores against one and would be vacuous.
        print(f"\n{len(feats)} poses over "
              f"{feats.ligand.nunique() if len(feats) else 0} ligands -> "
              f"{out_feats.name}; evaluation skipped (blind).", flush=True)
        return

    # --- evaluate ----------------------------------------------------------
    def ev(pick_fn, label: str) -> dict:
        picked, rand, oracle = [], [], []
        for lig, info in per_ligand.items():
            g = scored[scored.ligand == lig]
            if len(g) < 2:
                continue
            sample = pick_fn(lig, info)
            t = truth.get((lig, sample))
            if t is None or t != t:
                continue
            picked.append(t)
            rand.append(g.lddt_pli.mean())
            oracle.append(g.lddt_pli.max())
        if not picked:
            return {"selector": label, "n": 0}
        from scipy import stats
        picked, rand, oracle = map(np.array, (picked, rand, oracle))
        w = stats.wilcoxon(picked, rand) if len(picked) > 10 else None
        return {"selector": label, "n": int(len(picked)),
                "selected": round(float(picked.mean()), 4),
                "random": round(float(rand.mean()), 4),
                "oracle": round(float(oracle.mean()), 4),
                "delta_vs_random": round(float((picked - rand).mean()), 4),
                "frac_beating_random": round(float((picked > rand).mean()), 3),
                "wilcoxon_p": (round(float(w.pvalue), 5) if w is not None else None),
                "oracle_captured": round(float((picked.mean() - rand.mean())
                                               / max(1e-9, oracle.mean() - rand.mean())), 3)}

    results = [
        ev(lambda l, i: i["medoid_sample"], "consensus: medoid"),
        ev(lambda l, i: i["largest_cluster_medoid_sample"],
           "consensus: largest-cluster medoid"),
    ]

    print(f"\n{'selector':40s} {'sel':>7s} {'rand':>7s} {'oracle':>7s} {'delta':>8s} "
          f"{'beat%':>6s} {'p':>9s} {'capt':>6s}")
    print("-" * 96)
    for r in results:
        if not r.get("n"):
            continue
        print(f"{r['selector']:40s} {r['selected']:7.4f} {r['random']:7.4f} "
              f"{r['oracle']:7.4f} {r['delta_vs_random']:+8.4f} "
              f"{r['frac_beating_random']*100:5.1f}% {str(r['wilcoxon_p']):>9s} "
              f"{r['oracle_captured']:6.3f}")

    # does pool spread predict which ligands are hard?
    spread = pd.DataFrame([{"ligand": l, **i} for l, i in per_ligand.items()])
    best = scored.groupby("ligand").lddt_pli.max().rename("oracle")
    j = spread.set_index("ligand").join(best).dropna()
    from scipy import stats as st
    rho_spread = st.spearmanr(j.mean_pairwise_rmsd, j.oracle)
    rho_frac = st.spearmanr(j.largest_cluster_frac, j.oracle)
    print(f"\npool spread vs achievable quality: rho={rho_spread.statistic:+.3f} "
          f"(p={rho_spread.pvalue:.4f})")
    print(f"largest-cluster fraction vs quality: rho={rho_frac.statistic:+.3f} "
          f"(p={rho_frac.pvalue:.4f})")
    print(f"ligands whose pool splits into >1 cluster: "
          f"{int((j.n_clusters > 1).sum())}/{len(j)}")

    out = DATA_PROCESSED / f"consensus_selector_{a.tag}_{a.arm}.json"
    out.write_text(json.dumps({
        "arm": a.arm, "results": results,
        "spread_vs_oracle_rho": round(float(rho_spread.statistic), 4),
        "spread_vs_oracle_p": round(float(rho_spread.pvalue), 5),
        "largest_cluster_frac_vs_oracle_rho": round(float(rho_frac.statistic), 4),
        "n_multimodal_pools": int((j.n_clusters > 1).sum()),
        "per_ligand": per_ligand}, indent=1))
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
