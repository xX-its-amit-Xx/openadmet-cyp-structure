"""Does cross-ENGINE agreement select a better pose than within-engine consensus?

This is the signal FINDING 002 named as most promising and FINDING 003 left untested.
The logic: two independently trained co-folders making the same mistake is less likely
than one making it alone, so a pose that a *different* engine also produces is more
likely to be right. It is a different question from "do this engine's own samples agree",
which is already half of the FINDING 003 selector and only worth +0.010 alone.

Features computed per pose (all within-ligand, none touching the reference structure):

  - `xeng_min_rmsd`  — distance from this pose to the NEAREST pose from the other engine.
  - `xeng_mean_rmsd` — distance to that engine's poses on average.
  - `xeng_support`   — how many of the other engine's poses sit within 2 A.

Two passes are needed because the two pools live in different Modal workspaces, and a
process can only hold one `MODAL_PROFILE`. Each pass caches ligand coordinates to a local
npz, so the analysis itself runs offline and is repeatable without touching Modal.

    python scripts/structure/cross_engine.py cache --tag val87b --profile ashenoydiscovery
    python scripts/structure/cross_engine.py cache --tag chai87  --profile xx-its-amit-xx
    python scripts/structure/cross_engine.py test  --a val87b --b chai87
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED, free_gb, safe_workers  # noqa: E402

CACHE = DATA_PROCESSED / "xeng_cache"


def scratch_root() -> Path:
    d = Path("D:/cyp_scratch")
    d.mkdir(parents=True, exist_ok=True)
    return d if free_gb(d) >= 5.0 else Path("C:/cyp_scratch")


def cache_pool(tag: str, profile: str | None, workers: int) -> None:
    """Pull each pose's ligand coordinates + elements and store one npz per ligand."""
    if profile:
        os.environ["MODAL_PROFILE"] = profile
    import modal

    from cypstruct import pose as P

    out = CACHE / tag
    out.mkdir(parents=True, exist_ok=True)
    vol = modal.Volume.from_name("cyp-pool")
    entries = list(vol.iterdir(f"/{tag}", recursive=True))
    jobs = sorted({e.path.strip("/").split("/")[1] for e in entries
                   if len(e.path.strip("/").split("/")) > 1
                   and "__" in e.path.strip("/").split("/")[1]})
    todo = [j for j in jobs if not (out / f"{j.split('__')[0]}.npz").exists()]
    print(f"{tag}: {len(jobs)} jobs, {len(todo)} to cache", flush=True)
    if not todo:
        return

    root = Path(tempfile.mkdtemp(prefix=f"xeng_{tag}_", dir=str(scratch_root())))
    nw = safe_workers(workers, ram_per_worker_gb=1.0)
    lock = threading.Lock()

    def do(job: str) -> str:
        lig = job.split("__")[0]
        w = root / job
        w.mkdir(parents=True, exist_ok=True)
        try:
            for e in vol.iterdir(f"/{tag}/{job}"):
                fn = e.path.split("/")[-1]
                if fn.endswith(".cif"):
                    (w / fn).write_bytes(b"".join(vol.read_file(e.path)))
            names, coords, elems = [], [], None
            for cif in sorted(w.glob("*.cif")):
                try:
                    m = P.load_structure(cif)
                except Exception:
                    continue
                if len(m.lig_xyz) == 0:
                    continue
                names.append(cif.stem)
                coords.append(m.lig_xyz)
                elems = [e.upper() for e in m.lig_elem]
            if names:
                with lock:
                    np.savez_compressed(out / f"{lig}.npz",
                                        names=np.array(names),
                                        xyz=np.array(coords, dtype=np.float32),
                                        elems=np.array(elems))
            return lig
        finally:
            shutil.rmtree(w, ignore_errors=True)

    try:
        with ThreadPoolExecutor(max_workers=nw) as ex:
            futs = {ex.submit(do, j): j for j in todo}
            for k, f in enumerate(as_completed(futs), 1):
                f.result()
                if k % 20 == 0:
                    print(f"  cached {k}/{len(todo)}", flush=True)
    finally:
        shutil.rmtree(root, ignore_errors=True)
    print(f"-> {out}", flush=True)


def sym_rmsd(xa: np.ndarray, xb: np.ndarray, el: np.ndarray) -> float:
    """Element-constrained assignment RMSD. A lower bound, consistent across pairs."""
    from scipy.optimize import linear_sum_assignment

    d = np.linalg.norm(xa[:, None, :] - xb[None, :, :], axis=2)
    cost = d.copy()
    cost[el[:, None] != el[None, :]] = d.max() + 1e3
    r, c = linear_sum_assignment(cost)
    return float(np.sqrt((d[r, c] ** 2).mean()))


def test(tag_a: str, tag_b: str) -> None:
    from scipy import stats

    sa = pd.read_csv(DATA_PROCESSED / f"poses_scored_{tag_a}.csv")
    sa = sa[sa.arm == "unsteered"]
    ligs_a = {p.stem for p in (CACHE / tag_a).glob("*.npz")}
    ligs_b = {p.stem for p in (CACHE / tag_b).glob("*.npz")}
    shared = sorted(ligs_a & ligs_b)
    print(f"{tag_a}: {len(ligs_a)} ligands | {tag_b}: {len(ligs_b)} | shared: {len(shared)}",
          flush=True)

    rows = []
    for lig in shared:
        A = np.load(CACHE / tag_a / f"{lig}.npz", allow_pickle=True)
        B = np.load(CACHE / tag_b / f"{lig}.npz", allow_pickle=True)
        xa, xb = A["xyz"], B["xyz"]
        ea, eb = A["elems"], B["elems"]
        if xa.shape[1] != xb.shape[1] or sorted(ea.tolist()) != sorted(eb.tolist()):
            continue      # different molecule parsed; skip rather than compare nonsense
        for i, nm in enumerate(A["names"]):
            ds = [sym_rmsd(xa[i], xb[j], ea) for j in range(len(xb))]
            rows.append({"ligand": lig, "sample": str(nm),
                         "xeng_min_rmsd": float(np.min(ds)),
                         "xeng_mean_rmsd": float(np.mean(ds)),
                         "xeng_support": int(sum(1 for d in ds if d <= 2.0))})
    feats = pd.DataFrame(rows)
    feats.to_csv(DATA_PROCESSED / f"xeng_features_{tag_a}_vs_{tag_b}.csv", index=False)
    df = sa.merge(feats, on=["ligand", "sample"], how="inner")
    print(f"merged {len(df)} poses over {df.ligand.nunique()} ligands", flush=True)

    def ev(score, label):
        d = df.assign(_s=score).dropna(subset=["_s", "lddt_pli"])
        pick, rnd, orc = [], [], []
        for _l, g in d.groupby("ligand"):
            if len(g) < 2:
                continue
            pick.append(g.loc[g._s.idxmax(), "lddt_pli"])
            rnd.append(g.lddt_pli.mean())
            orc.append(g.lddt_pli.max())
        pick, rnd, orc = map(np.array, (pick, rnd, orc))
        w = stats.wilcoxon(pick, rnd)
        print(f"  {label:44s} sel={pick.mean():.4f} rand={rnd.mean():.4f} "
              f"orc={orc.mean():.4f} delta={(pick-rnd).mean():+.4f} "
              f"beat={100*(pick>rnd).mean():4.1f}% p={w.pvalue:.4f} n={len(pick)}")
        return (pick - rnd).mean(), w.pvalue

    zw = lambda c: ((df[c] - df.groupby("ligand")[c].transform("mean"))
                    / (df.groupby("ligand")[c].transform("std") + 1e-9))

    print("\n=== cross-engine agreement as a selector ===")
    ev(-zw("xeng_min_rmsd"), "xeng: nearest other-engine pose (closer better)")
    ev(-zw("xeng_mean_rmsd"), "xeng: mean distance to other engine")
    ev(zw("xeng_support"), "xeng: count of other-engine poses within 2 A")

    print("\n=== against, and combined with, the FINDING 003 selector ===")
    orient = DATA_PROCESSED / f"orientation_features_{tag_a}_unsteered.csv"
    cons = DATA_PROCESSED / f"consensus_features_{tag_a}_unsteered.csv"
    if orient.exists() and cons.exists():
        o = pd.read_csv(orient)
        c = pd.read_csv(cons)
        df2 = df.merge(o, on=["ligand", "sample"], how="left").merge(
            c, on=["ligand", "sample"], how="left").dropna(
            subset=["n_pocket_residues_touched", "mean_rmsd_to_others"])
        globals()["df"] = df2
        zw2 = lambda col: ((df2[col] - df2.groupby("ligand")[col].transform("mean"))
                           / (df2.groupby("ligand")[col].transform("std") + 1e-9))
        base = 0.5 * zw2("n_pocket_residues_touched") - zw2("mean_rmsd_to_others")
        df = df2
        ev(base, "FINDING 003 baseline (contacts + medoid)")
        for w_x in (0.25, 0.5, 1.0):
            ev(base - w_x * zw2("xeng_min_rmsd"), f"  + {w_x} * cross-engine agreement")

    out = DATA_PROCESSED / f"cross_engine_{tag_a}_vs_{tag_b}.json"
    out.write_text(json.dumps({"a": tag_a, "b": tag_b, "n_shared": len(shared),
                               "n_poses": int(len(df))}, indent=1))
    print(f"\n-> {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["cache", "test"])
    ap.add_argument("--tag")
    ap.add_argument("--profile")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--a")
    ap.add_argument("--b")
    a = ap.parse_args()
    if a.cmd == "cache":
        cache_pool(a.tag, a.profile, a.workers)
    else:
        test(a.a, a.b)
