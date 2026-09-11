"""Orientation-aware pose features: where the rest of the molecule sits.

FINDING 002 diagnosed why every feature tried so far failed. All of them are
**anchor-local** — they describe the ligand atom nearest the iron. But the iron anchor is
nearly constant across a ligand's samples (1,424 of 1,440 type II poses coordinate), so
those features have no within-ligand variance to rank with, however large their
between-ligand effect. And the addendum showed the anchor can be textbook-perfect while
the molecule is rotated the wrong way about it.

So these features describe the **rest** of the molecule:

1. **Pocket contact fingerprint.** A binary vector over the CYP3A4 pocket residues: which
   ones does this pose touch? This is exactly "where the substituents sit", expressed in
   a frame that needs no alignment to anything.

2. **Contact-space consensus.** Per ligand, the mean Jaccard similarity of a pose's
   fingerprint to every other sample's. The coordinate-space medoid was the best selector
   so far (+0.0100, p=0.087); this is the same idea moved into contact space, where a pose
   that is rotated about the anchor looks *different* even though its anchor atom has
   barely moved.

3. **Rotation about the Fe-donor axis.** For coordinated poses, the azimuth of the
   ligand's centre of mass about the Fe→donor vector. This is the degree of freedom that
   the anchor leaves completely unconstrained, and therefore the one that should carry the
   discriminating information.

4. **Anchor-distal geometry.** Radius of gyration about the donor atom, the distance from
   the donor to the ligand centroid, and the buried fraction of the half of the molecule
   furthest from the iron.

Everything is reported and tested **within ligand**, per the standing rule: a feature
that is constant across a ligand's samples is not a selector at all.

    python scripts/structure/orientation_features.py --tag val87b --arm unsteered
"""
from __future__ import annotations

import argparse
import json
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

import modal  # noqa: E402

from cypstruct import pose as P  # noqa: E402
from cypstruct.paths import DATA_PROCESSED, free_gb, safe_workers  # noqa: E402
from cypstruct.qmscore import geometry as G  # noqa: E402
from cypstruct.targets import CYP3A4_POCKET  # noqa: E402

POCKET = sorted({r for v in CYP3A4_POCKET.values() for r in v})
CONTACT_CUTOFF = 4.5


def scratch_root() -> Path:
    d = Path("D:/cyp_scratch")
    d.mkdir(parents=True, exist_ok=True)
    if free_gb(d) >= 5.0:
        return d
    c = Path("C:/cyp_scratch")
    c.mkdir(parents=True, exist_ok=True)
    return c


def pose_features(model: P.Complex) -> dict:
    """Orientation features for one pose. No reference structure is used anywhere here.

    That matters: anything derived from the answer is a leak. Every quantity below is
    computable from the prediction alone at inference time.
    """
    out: dict = {}
    lig = model.lig_xyz
    if len(lig) == 0:
        return out

    # --- 1. pocket contact fingerprint ------------------------------------
    fp = np.zeros(len(POCKET), dtype=np.int8)
    if len(model.prot_xyz):
        d = np.linalg.norm(model.prot_xyz[:, None, :] - lig[None, :, :], axis=2).min(1)
        touched = {model.prot_key[i][1] for i in np.where(d <= CONTACT_CUTOFF)[0]}
        for k, r in enumerate(POCKET):
            fp[k] = 1 if r in touched else 0
    out["_fp"] = fp
    out["n_pocket_residues_touched"] = int(fp.sum())

    # --- 3/4. anchor-relative geometry ------------------------------------
    fe, normal = G.porphyrin_frame(np.asarray(model.heme_xyz, float), model.heme_atom,
                                   model.axial_sg, model.heme_elem)
    centroid = lig.mean(0)
    if fe is not None:
        d_fe = np.linalg.norm(lig - fe, axis=1)
        donor_idx = [i for i, e in enumerate(model.lig_elem) if e.upper() in ("N", "O", "S")]
        anchor = int(min(donor_idx, key=lambda i: d_fe[i])) if donor_idx else int(np.argmin(d_fe))
        a = lig[anchor]
        out["donor_to_centroid"] = float(np.linalg.norm(centroid - a))
        out["rg_about_donor"] = float(np.sqrt(((lig - a) ** 2).sum(1).mean()))
        out["centroid_fe_dist"] = float(np.linalg.norm(centroid - fe))

        # azimuth of the centroid about the Fe->donor axis. The anchor leaves this free.
        axis = a - fe
        na = np.linalg.norm(axis)
        if na > 1e-6 and normal is not None:
            axis = axis / na
            # build an orthonormal frame perpendicular to the axis, using the porphyrin
            # normal as the reference direction so the angle is comparable across poses
            ref = normal - np.dot(normal, axis) * axis
            if np.linalg.norm(ref) > 1e-6:
                ref /= np.linalg.norm(ref)
                perp = np.cross(axis, ref)
                v = centroid - a
                v = v - np.dot(v, axis) * axis
                if np.linalg.norm(v) > 1e-6:
                    out["azimuth_about_axis"] = float(
                        np.degrees(np.arctan2(np.dot(v, perp), np.dot(v, ref))))

        # burial of the half of the molecule FURTHEST from the iron
        far = lig[d_fe > np.median(d_fe)]
        if len(far) and len(model.prot_xyz):
            dd = np.linalg.norm(model.prot_xyz[:, None, :] - far[None, :, :], axis=2)
            out["distal_half_contacts"] = int((dd < CONTACT_CUTOFF).sum())
            out["distal_half_buried_frac"] = float(((dd < CONTACT_CUTOFF).sum(0) >= 6).mean())
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="val87b")
    ap.add_argument("--arm", default="unsteered")
    ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()

    cache = DATA_PROCESSED / f"orientation_features_{a.tag}_{a.arm}.csv"
    if cache.exists():
        print(f"using cached {cache}", flush=True)
        feats = pd.read_csv(cache)
    else:
        scored = pd.read_csv(DATA_PROCESSED / f"poses_scored_{a.tag}.csv")
        scored = scored[scored.arm == a.arm]
        jobs = sorted({f"{r.ligand}__{a.arm}__s1" for r in scored.itertuples()})
        vol = modal.Volume.from_name("cyp-pool")
        root = Path(tempfile.mkdtemp(prefix="cyporient_", dir=str(scratch_root())))
        nw = safe_workers(a.workers, ram_per_worker_gb=1.0)
        print(f"{len(jobs)} jobs, {nw} worker(s)", flush=True)
        lock = threading.Lock()
        rows: list[dict] = []

        def do(job: str) -> list[dict]:
            lig_id = job.split("__")[0]
            w = root / job
            w.mkdir(parents=True, exist_ok=True)
            local: list[dict] = []
            try:
                for e in vol.iterdir(f"/{a.tag}/{job}"):
                    fn = e.path.split("/")[-1]
                    if fn.endswith(".cif"):
                        (w / fn).write_bytes(b"".join(vol.read_file(e.path)))
                fps, recs = [], []
                for cif in sorted(w.glob("*.cif")):
                    try:
                        m = P.load_structure(cif)
                    except Exception:
                        continue
                    f = pose_features(m)
                    if not f:
                        continue
                    fps.append(f.pop("_fp"))
                    recs.append({"ligand": lig_id, "sample": cif.stem, **f})
                # contact-space consensus: mean Jaccard to the other samples
                if len(fps) >= 2:
                    A = np.array(fps, dtype=float)
                    inter = A @ A.T
                    tot = A.sum(1)
                    union = tot[:, None] + tot[None, :] - inter
                    jac = np.divide(inter, np.maximum(union, 1e-9))
                    np.fill_diagonal(jac, np.nan)
                    mean_jac = np.nanmean(jac, axis=1)
                    for i, r in enumerate(recs):
                        r["contact_consensus"] = float(mean_jac[i])
                local = recs
            except Exception:
                local = []
            finally:
                shutil.rmtree(w, ignore_errors=True)
            return local

        try:
            with ThreadPoolExecutor(max_workers=nw) as ex:
                futs = {ex.submit(do, j): j for j in jobs}
                for k, fut in enumerate(as_completed(futs), 1):
                    r = fut.result()
                    with lock:
                        rows.extend(r)
                    if k % 20 == 0:
                        print(f"  {k}/{len(jobs)} ligands, {len(rows)} poses", flush=True)
        finally:
            shutil.rmtree(root, ignore_errors=True)

        feats = pd.DataFrame(rows)
        feats.to_csv(cache, index=False)
        print(f"-> {cache}", flush=True)

    # --- evaluate, WITHIN LIGAND ------------------------------------------
    scored = pd.read_csv(DATA_PROCESSED / f"poses_scored_{a.tag}.csv")
    scored = scored[scored.arm == a.arm]
    df = scored.merge(feats, on=["ligand", "sample"], how="inner")
    print(f"\nmerged {len(df)} poses over {df.ligand.nunique()} ligands", flush=True)

    from scipy import stats

    cand = [c for c in ("contact_consensus", "n_pocket_residues_touched",
                        "distal_half_contacts", "distal_half_buried_frac",
                        "donor_to_centroid", "rg_about_donor", "centroid_fe_dist")
            if c in df.columns]

    print(f"\n{'feature':32s} {'sign':>5s} {'within-rho':>11s} {'sel':>7s} {'rand':>7s} "
          f"{'delta':>8s} {'beat%':>6s} {'p':>9s}")
    print("-" * 96)
    results = []
    for c in cand:
        # within-ligand Spearman: the only correlation that means anything for selection
        rhos = []
        for _l, g in df.groupby("ligand"):
            gg = g.dropna(subset=[c, "lddt_pli"])
            if len(gg) >= 5 and gg[c].nunique() > 1:
                rhos.append(stats.spearmanr(gg[c], gg.lddt_pli).statistic)
        rhos = [r for r in rhos if r == r]
        mrho = float(np.mean(rhos)) if rhos else float("nan")
        sign = 1 if mrho >= 0 else -1
        pick, rnd = [], []
        for _l, g in df.groupby("ligand"):
            gg = g.dropna(subset=[c, "lddt_pli"])
            if len(gg) < 2:
                continue
            s = gg[c] * sign
            pick.append(gg.loc[s.idxmax(), "lddt_pli"])
            rnd.append(gg.lddt_pli.mean())
        if not pick:
            continue
        pick, rnd = np.array(pick), np.array(rnd)
        w = stats.wilcoxon(pick, rnd) if len(pick) > 10 else None
        pv = float(w.pvalue) if w else float("nan")
        results.append({"feature": c, "sign": sign, "within_rho": round(mrho, 4),
                        "selected": round(float(pick.mean()), 4),
                        "random": round(float(rnd.mean()), 4),
                        "delta": round(float((pick - rnd).mean()), 4),
                        "beat_frac": round(float((pick > rnd).mean()), 3),
                        "p": round(pv, 5), "n": int(len(pick))})
        print(f"{c:32s} {sign:+5d} {mrho:+11.4f} {pick.mean():7.4f} {rnd.mean():7.4f} "
              f"{(pick-rnd).mean():+8.4f} {100*(pick>rnd).mean():5.1f}% {pv:9.4f}")

    out = DATA_PROCESSED / f"orientation_selector_{a.tag}_{a.arm}.json"
    out.write_text(json.dumps({"arm": a.arm, "n_poses": int(len(df)),
                               "results": results}, indent=1))
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
