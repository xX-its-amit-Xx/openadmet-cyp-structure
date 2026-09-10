"""Score a pose pool against crystal ground truth, and report the oracle first.

This is the measurement loop the whole project runs on. Every claim about steering,
selection, or the physics scorer has to come out of here.

**The pool oracle is reported before anything else, and that is deliberate.** The oracle
is the best LDDT-PLI achievable by *any* pose in the pool — the ceiling a perfect selector
would reach. A selection result quoted without it is uninterpretable: 0.45 selected against
a 0.48 oracle means the selector is nearly perfect and the generator is the problem, while
0.45 against a 0.72 oracle means the opposite. Those need opposite responses, and the two
numbers look identical if you only print one.

Outputs `data/processed/pool_scorecard.json`:
  - per arm (steered / unsteered): oracle, z-hybrid selection, random-pose baseline
  - the realised Fe-donor distance distribution, which is how we tell whether the
    steering constraint was actually honoured or silently ignored
  - per-ligand rows, so failures can be looked at individually

Usage:
    python scripts/structure/score_pool.py --pool <dir> --tag val99
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cypstruct import pose as P  # noqa: E402
from cypstruct.paths import DATA_PROCESSED, REFERENCE  # noqa: E402
from cypstruct.qmscore import geometry as G  # noqa: E402
from cypstruct.select import Pose, pool_confidence, z_hybrid  # noqa: E402
from cypstruct.targets import fetch_cif  # noqa: E402


def load_reference(pdb_id: str, ligand_code: str) -> P.Complex | None:
    """The deposited structure, restricted to ONE chain.

    Chain restriction is not optional. CYP3A4 deposits as monomer, dimer and (in the new
    cryoEM work) a symmetric trimer; without it we would compare a single predicted copy
    against three superimposed reference copies and every distance would be meaningless.
    """
    try:
        cif = fetch_cif(pdb_id)
    except Exception:
        return None
    import gemmi
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    for chain in st[0]:
        if any(r.name.strip().upper() == ligand_code.upper() for r in chain):
            try:
                return P.load_structure(cif, ligand_code=ligand_code,
                                        assembly_chain=chain.name)
            except Exception:
                return None
    return None


def parse_confidence(path: Path) -> float | None:
    """Boltz confidence for one sample, sign-corrected so higher is better.

    `complex_ipde` is an ERROR estimate. Taking a raw max over it selects the worst pose,
    and does so silently — the numbers all look plausible.
    """
    for cand in (path.parent / f"confidence_{path.stem}.json",
                 path.with_suffix(".json"),
                 path.parent / f"{path.stem}_confidence.json"):
        if cand.exists():
            try:
                d = json.loads(cand.read_text())
            except json.JSONDecodeError:
                continue
            if "complex_ipde" in d:
                return -float(d["complex_ipde"])
            for k in ("confidence_score", "ptm", "iptm", "complex_plddt"):
                if k in d:
                    return float(d[k])
    return None


def score_pool(pool_dir: Path, ligands: pd.DataFrame) -> dict:
    ref_cache: dict[str, P.Complex | None] = {}
    rows: list[dict] = []
    poses_by_arm: dict[str, list[Pose]] = defaultdict(list)

    lig_meta = {r.id: r for r in ligands.itertuples()}

    for job_dir in sorted(p for p in pool_dir.iterdir() if p.is_dir()):
        name = job_dir.name
        parts = name.split("__")
        lig_id, arm = parts[0], (parts[1] if len(parts) > 1 else "unknown")
        meta = lig_meta.get(lig_id)
        if meta is None:
            continue
        key = f"{meta.pdb}:{lig_id}"
        if key not in ref_cache:
            ref_cache[key] = load_reference(meta.pdb, lig_id)
        ref = ref_cache[key]
        if ref is None or len(ref.lig_xyz) == 0:
            continue

        for cif in sorted(job_dir.glob("*.cif")):
            try:
                model = P.load_structure(cif)
            except Exception:
                continue
            if len(model.lig_xyz) == 0:
                continue
            perm = P.best_ligand_mapping(meta.smiles, model, ref)
            lddt = P.lddt_pli(model, ref, lig_perm=perm)
            rmsd = P.bisy_rmsd(model, ref, lig_perm=perm)
            geo = G.compute(model.lig_xyz, model.lig_elem, model.prot_xyz,
                            model.heme_xyz, model.heme_atom, model.axial_sg)
            valid, why = G.hard_validity(geo)
            conf = parse_confidence(cif)
            rows.append(dict(ligand=lig_id, arm=arm, path=str(cif), pdb=meta.pdb,
                             cls=meta.cls, lddt_pli=lddt, bisy_rmsd=rmsd,
                             confidence=conf, valid=valid, invalid_reason=why,
                             fe_donor_dist=geo.fe_donor_dist,
                             s_fe_donor_angle=geo.s_fe_donor_angle,
                             is_coordinated=geo.is_coordinated,
                             frac_proximal=geo.frac_atoms_proximal))
            if conf is not None and lddt == lddt:
                poses_by_arm[arm].append(
                    Pose(ligand=lig_id, engine="boltz2", path=str(cif),
                         raw_confidence=conf))

    df = pd.DataFrame(rows)
    if df.empty:
        return {"error": "no scoreable poses found", "pool_dir": str(pool_dir)}
    truth = {r.path: r.lddt_pli for r in df.itertuples() if r.lddt_pli == r.lddt_pli}

    summary: dict = {"pool_dir": str(pool_dir), "n_poses": int(len(df)),
                     "n_ligands": int(df.ligand.nunique()), "arms": {}}

    for arm, sub in df.groupby("arm"):
        per_lig_oracle = sub.groupby("ligand").lddt_pli.max()
        chosen, wins = z_hybrid(poses_by_arm.get(arm, []))
        sel = [truth[p.path] for p in chosen.values() if p.path in truth]
        rng = np.random.default_rng(0)
        rand = [float(rng.choice(g.lddt_pli.dropna().values))
                for _l, g in sub.groupby("ligand") if g.lddt_pli.notna().any()]
        coord = sub[sub.is_coordinated]
        summary["arms"][arm] = {
            "n_poses": int(len(sub)),
            "n_ligands": int(sub.ligand.nunique()),
            "pool_oracle_mean_lddt_pli": round(float(per_lig_oracle.mean()), 4),
            "zhybrid_selected_mean_lddt_pli": (round(float(np.mean(sel)), 4) if sel else None),
            "random_pose_mean_lddt_pli": (round(float(np.mean(rand)), 4) if rand else None),
            "selection_gap_to_oracle": (
                round(float(per_lig_oracle.mean() - np.mean(sel)), 4) if sel else None),
            "frac_poses_valid": round(float(sub.valid.mean()), 4),
            "frac_poses_coordinated": round(float(sub.is_coordinated.mean()), 4),
            # The check that says whether steering was honoured at all. If the steered arm
            # does not shift this distribution toward 2.0-2.4 A, the constraint is being
            # ignored and no conclusion about selection is warranted.
            "fe_donor_dist_percentiles": {
                f"p{q}": (round(float(np.percentile(sub.fe_donor_dist.dropna(), q)), 3)
                          if sub.fe_donor_dist.notna().any() else None)
                for q in (5, 25, 50, 75, 95)},
            "n_inside_coordination_window": int(
                ((sub.fe_donor_dist >= G.COORD_LO) & (sub.fe_donor_dist <= G.COORD_HI)).sum()),
            "mean_bisy_rmsd": round(float(sub.bisy_rmsd.mean(skipna=True)), 3),
            "frac_under_2A": round(float((sub.bisy_rmsd < 2.0).mean()), 4),
            "engine_wins": wins,
            "lddt_pli_on_coordinated_poses": (
                round(float(coord.lddt_pli.mean()), 4) if len(coord) else None),
        }

    if len(summary["arms"]) == 2 and "steered" in summary["arms"] and "unsteered" in summary["arms"]:
        s, u = summary["arms"]["steered"], summary["arms"]["unsteered"]
        summary["steering_effect"] = {
            "oracle_delta": round(s["pool_oracle_mean_lddt_pli"]
                                  - u["pool_oracle_mean_lddt_pli"], 4),
            "selected_delta": (round(s["zhybrid_selected_mean_lddt_pli"]
                                     - u["zhybrid_selected_mean_lddt_pli"], 4)
                               if s["zhybrid_selected_mean_lddt_pli"] is not None
                               and u["zhybrid_selected_mean_lddt_pli"] is not None else None),
            "coordination_rate_delta": round(s["frac_poses_coordinated"]
                                             - u["frac_poses_coordinated"], 4),
            "note": ("If coordination_rate_delta is ~0 the constraint was not honoured; "
                     "no claim about selection follows until that is fixed."),
        }

    out_csv = DATA_PROCESSED / "pool_poses_scored.csv"
    df.to_csv(out_csv, index=False)
    summary["per_pose_csv"] = str(out_csv)
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True, help="directory of <ligand>__<arm>__<seed>/ dirs")
    ap.add_argument("--ligands", default=str(DATA_PROCESSED / "validation_ligands.csv"))
    a = ap.parse_args()
    ligands = pd.read_csv(a.ligands)
    summary = score_pool(Path(a.pool), ligands)
    out = DATA_PROCESSED / "pool_scorecard.json"
    prev = json.loads(out.read_text()) if out.exists() else []
    if isinstance(prev, dict):
        prev = [prev]
    prev.append(summary)
    out.write_text(json.dumps(prev, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
