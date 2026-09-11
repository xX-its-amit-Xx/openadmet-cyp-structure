"""Pull a Modal pose pool one job at a time, score it, and discard the files.

Streaming rather than bulk-download because the numbers do not fit: 168 jobs x 20
samples is 3,360 mmCIFs of a 503-residue protein plus heme plus ligand, on the order of
1.7 GB, against a local disk that has repeatedly sat near zero. Peak usage here is one
job, about 10 MB.

Reports the **pool oracle before any selection number**. A selected score without its
ceiling is uninterpretable: 0.45 against a 0.48 oracle means the selector is nearly
perfect and the generator is the bottleneck; 0.45 against 0.72 means the reverse. Those
demand opposite responses and look identical if only one is printed.

    python scripts/structure/collect_and_score.py --tag val87b
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

import modal  # noqa: E402

from cypstruct import pose as P  # noqa: E402
from cypstruct.paths import DATA_PROCESSED, resource_headroom, safe_workers  # noqa: E402
from cypstruct.qmscore import geometry as G  # noqa: E402
from cypstruct.select import Pose, z_hybrid  # noqa: E402
from cypstruct.targets import fetch_cif  # noqa: E402

_REF_CACHE: dict[str, P.Complex | None] = {}


def reference_for(pdb_id: str, ligand_code: str) -> P.Complex | None:
    """Deposited structure for one ligand, restricted to the chain that contains it.

    Chain restriction is mandatory, not tidiness. CYP3A4 deposits as monomer, dimer and
    (in the new cryoEM work) a symmetric trimer; parsing the whole file would compare one
    predicted copy against several superimposed reference copies.
    """
    key = f"{pdb_id}:{ligand_code}"
    if key in _REF_CACHE:
        return _REF_CACHE[key]
    out = None
    try:
        import gemmi

        cif = fetch_cif(pdb_id)
        st = gemmi.read_structure(str(cif))
        st.setup_entities()
        for chain in st[0]:
            if any(r.name.strip().upper() == ligand_code.upper() for r in chain):
                out = P.load_structure(cif, ligand_code=ligand_code,
                                       assembly_chain=chain.name)
                break
    except Exception:
        out = None
    _REF_CACHE[key] = out
    return out


def confidence_from(path: Path) -> float | None:
    """Boltz confidence, sign-corrected so that higher is better.

    `complex_ipde` is an ERROR estimate. A raw max over it selects the WORST pose, and
    does so silently because the numbers all look reasonable.
    """
    cand = path.parent / f"confidence_{path.stem}.json"
    if not cand.exists():
        return None
    try:
        d = json.loads(cand.read_text())
    except json.JSONDecodeError:
        return None
    if "complex_ipde" in d:
        return -float(d["complex_ipde"])
    for k in ("confidence_score", "ptm", "iptm", "complex_plddt"):
        if k in d:
            return float(d[k])
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--ligands", default=str(DATA_PROCESSED / "validation_ligands.csv"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()

    ligands = pd.read_csv(a.ligands)
    meta = {r.id: r for r in ligands.itertuples()}
    vol = modal.Volume.from_name("cyp-pool")

    job_ids = sorted({e.path.rstrip("/").split("/")[-1] for e in vol.iterdir(f"/{a.tag}")})
    job_ids = [j for j in job_ids if "__" in j]
    if a.limit:
        job_ids = job_ids[: a.limit]
    print(f"{len(job_ids)} jobs in /{a.tag}", flush=True)

    rows: list[dict] = []
    poses_by_arm: dict[str, list[Pose]] = defaultdict(list)
    tmp_root = Path(tempfile.mkdtemp(prefix="cypscore_"))

    # Parallelism is capped by measured headroom, not by core count. This box has no GPU
    # and has hit zero free disk more than once; saturating it degrades every other tool
    # in the session. Each worker holds one job's files (~10 MB) at a time.
    nw = safe_workers(a.workers, ram_per_worker_gb=1.5)
    print(f"headroom: {resource_headroom()} -> {nw} worker(s)", flush=True)

    lock = threading.Lock()

    def score_job(job: str) -> tuple[list[dict], list[Pose]]:
        lig_id, arm = job.split("__")[0], job.split("__")[1]
        m = meta.get(lig_id)
        if m is None:
            return [], []
        with lock:                       # reference parse + cache is not thread-safe
            ref = reference_for(m.pdb, lig_id)
        if ref is None or len(ref.lig_xyz) == 0:
            return [], []

        work = tmp_root / job
        work.mkdir(parents=True, exist_ok=True)
        out_rows: list[dict] = []
        out_poses: list[Pose] = []
        try:
            for e in vol.iterdir(f"/{a.tag}/{job}"):
                fn = e.path.split("/")[-1]
                if fn.endswith(".cif") or fn.startswith("confidence_"):
                    (work / fn).write_bytes(b"".join(vol.read_file(e.path)))

            for cif in sorted(work.glob("*.cif")):
                try:
                    model = P.load_structure(cif)
                except Exception:
                    continue
                if len(model.lig_xyz) == 0:
                    continue
                perm = P.best_ligand_mapping(m.smiles, model, ref)
                lddt = P.lddt_pli(model, ref, lig_perm=perm)
                rmsd = P.bisy_rmsd(model, ref, lig_perm=perm)
                geo = G.compute(model.lig_xyz, model.lig_elem, model.prot_xyz,
                                model.heme_xyz, model.heme_atom, model.axial_sg)
                valid, why = G.hard_validity(geo)
                conf = confidence_from(cif)
                out_rows.append(dict(
                    ligand=lig_id, arm=arm, pdb=m.pdb, cls=m.cls, sample=cif.stem,
                    lddt_pli=lddt, bisy_rmsd=rmsd, confidence=conf,
                    valid=valid, invalid_reason=why,
                    fe_donor_dist=geo.fe_donor_dist,
                    s_fe_donor_angle=geo.s_fe_donor_angle,
                    is_coordinated=geo.is_coordinated,
                    frac_proximal=geo.frac_atoms_proximal,
                    max_clash=geo.max_clash, n_contacts=geo.n_protein_contacts,
                    mapped=perm is not None))
                if conf is not None and lddt == lddt:
                    out_poses.append(Pose(ligand=lig_id, engine="boltz2",
                                          path=f"{job}/{cif.name}", raw_confidence=conf))
        except Exception:
            pass
        finally:
            shutil.rmtree(work, ignore_errors=True)
        return out_rows, out_poses

    try:
        with ThreadPoolExecutor(max_workers=nw) as ex:
            futs = {ex.submit(score_job, j): j for j in job_ids}
            for n, fut in enumerate(as_completed(futs), 1):
                try:
                    r, ps = fut.result()
                except Exception as exc:
                    print(f"  [job failed] {futs[fut]}: {type(exc).__name__}: {exc}",
                          flush=True)
                    continue
                rows.extend(r)
                arm_of_job = futs[fut].split("__")[1]
                for pz in ps:
                    poses_by_arm[arm_of_job].append(pz)
                if n % 20 == 0:
                    print(f"  scored {n}/{len(job_ids)} jobs, {len(rows)} poses",
                          flush=True)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("no scoreable poses")
    out_csv = DATA_PROCESSED / f"poses_scored_{a.tag}.csv"
    df.to_csv(out_csv, index=False)
    truth = {f"{r.ligand}__{r.arm}__s1/{r.sample}.cif": r.lddt_pli
             for r in df.itertuples() if r.lddt_pli == r.lddt_pli}

    summary: dict = {"tag": a.tag, "n_poses": int(len(df)),
                     "n_ligands": int(df.ligand.nunique()), "arms": {}}

    for arm, sub in df.groupby("arm"):
        oracle = sub.groupby("ligand").lddt_pli.max()
        chosen, wins = z_hybrid(poses_by_arm.get(arm, []))
        sel = [truth[p.path] for p in chosen.values() if p.path in truth]
        rng = np.random.default_rng(0)
        rand = [float(rng.choice(g.lddt_pli.dropna().values))
                for _l, g in sub.groupby("ligand") if g.lddt_pli.notna().any()]
        best_rmsd = sub.groupby("ligand").bisy_rmsd.min()
        summary["arms"][arm] = {
            "n_poses": int(len(sub)), "n_ligands": int(sub.ligand.nunique()),
            "pool_oracle_lddt_pli": round(float(oracle.mean()), 4),
            "zhybrid_selected_lddt_pli": (round(float(np.mean(sel)), 4) if sel else None),
            "random_pose_lddt_pli": (round(float(np.mean(rand)), 4) if rand else None),
            "selection_gap_to_oracle": (round(float(oracle.mean() - np.mean(sel)), 4)
                                        if sel else None),
            "mean_bisy_rmsd": round(float(sub.bisy_rmsd.mean(skipna=True)), 3),
            "best_per_ligand_bisy_rmsd": round(float(best_rmsd.mean()), 3),
            "frac_poses_under_2A": round(float((sub.bisy_rmsd < 2.0).mean()), 4),
            "frac_ligands_with_a_sub2A_pose": round(float((best_rmsd < 2.0).mean()), 4),
            "frac_poses_valid": round(float(sub.valid.mean()), 4),
            "frac_poses_coordinated": round(float(sub.is_coordinated.mean()), 4),
            "frac_mapped": round(float(sub.mapped.mean()), 4),
            "fe_donor_dist_p": {f"p{q}": (round(float(np.percentile(
                sub.fe_donor_dist.dropna(), q)), 3) if sub.fe_donor_dist.notna().any()
                else None) for q in (5, 25, 50, 75, 95)},
            "n_inside_coordination_window": int(
                ((sub.fe_donor_dist >= G.COORD_LO) & (sub.fe_donor_dist <= G.COORD_HI)).sum()),
            "engine_wins": wins,
        }

    if {"steered", "unsteered"} <= set(summary["arms"]):
        s, u = summary["arms"]["steered"], summary["arms"]["unsteered"]
        summary["steering_effect"] = {
            "oracle_delta": round(s["pool_oracle_lddt_pli"] - u["pool_oracle_lddt_pli"], 4),
            "selected_delta": (round(s["zhybrid_selected_lddt_pli"]
                                     - u["zhybrid_selected_lddt_pli"], 4)
                               if s["zhybrid_selected_lddt_pli"] is not None
                               and u["zhybrid_selected_lddt_pli"] is not None else None),
            "coordination_rate_delta": round(s["frac_poses_coordinated"]
                                             - u["frac_poses_coordinated"], 4),
            "sub2A_ligand_rate_delta": round(s["frac_ligands_with_a_sub2A_pose"]
                                             - u["frac_ligands_with_a_sub2A_pose"], 4),
            "caveat": ("Both arms carry the Cys442-SG to heme-FE bond, so 'unsteered' "
                       "means heme-anchored but ligand-unconstrained, not unconstrained."),
        }

    out = DATA_PROCESSED / "pool_scorecard.json"
    prev = json.loads(out.read_text()) if out.exists() else []
    if isinstance(prev, dict):
        prev = [prev]
    prev.append(summary)
    out.write_text(json.dumps(prev, indent=2))
    summary["per_pose_csv"] = str(out_csv)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
