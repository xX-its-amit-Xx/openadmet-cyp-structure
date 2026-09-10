"""Harvest every deposited CYP3A4 structure and measure the active-site geometry.

Purpose: calibrate the physics scorer's windows on DATA rather than on my reading of
the literature. Everything the scorer asserts about "type II coordinates at 2.0-2.3 A,
trans to the thiolate" should be a measured distribution over the ~120 public CYP3A4
entries, not a number quoted from a review.

It also produces the training substrate for the learned ranker: a list of
(pdb_id, chain, ligand_code, SMILES-able geometry, binding class).

Two traps this script is written around, both of which produce confident nonsense:

1. **Multiple copies.** CYP3A4 deposits as monomer, dimer, and (in the new cryoEM work)
   a symmetric trimer. Parsing a whole file gives N copies of the heme and the ligand;
   naive distance stats then mix intra- and inter-copy distances. We iterate per chain
   and pair each ligand to the heme *in its own chain*.

2. **The largest HET group is not always the ligand.** CYP3A4 has a well-documented
   peripheral / access-channel site well away from the iron, and entries often contain
   two copies of the same compound. Every HET group is measured and classified by its
   own distance to iron; nothing is silently designated "the ligand".

Writes: data/processed/cyp3a4_reference_geometry.parquet  (+ a summary JSON)
CIFs are cached to cold storage (O:), never to the near-full local disks.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cypstruct.paths import DATA_PROCESSED, REFERENCE  # noqa: E402
from cypstruct.targets import (  # noqa: E402
    HEME_ALIASES,
    IGNORE_HET,
    fetch_cif,
    rcsb_holo_structures,
)

# Binding-class boundaries. Deliberately wide; the point of this script is to
# replace them with measured percentiles.
COORD_MAX = 2.6      # Fe-donor distance below this = direct coordination (type II)
OVERHEAD_MAX = 6.5   # closest-atom distance below this = in the active site (type I)


def heme_frame(heme_xyz: np.ndarray, heme_atom: list[str]) -> tuple[np.ndarray | None, np.ndarray | None]:
    """(Fe position, unit normal of the porphyrin plane) from the four pyrrole N.

    The normal defines the distal/proximal axis. Without it there is no way to tell a
    physically possible pose from one placed on the buried face of the heme.
    """
    fe = None
    ns = []
    for p, a in zip(heme_xyz, heme_atom):
        au = a.upper()
        if au == "FE":
            fe = np.asarray(p, float)
        elif au in ("NA", "NB", "NC", "ND"):
            ns.append(p)
    if fe is None or len(ns) < 3:
        return fe, None
    N = np.asarray(ns, float) - np.asarray(ns, float).mean(0)
    # plane normal = singular vector with the smallest singular value
    _u, _s, vt = np.linalg.svd(N)
    n = vt[-1]
    return fe, n / np.linalg.norm(n)


def analyse_entry(pdb_id: str) -> list[dict]:
    """One row per (chain, HET group). Returns [] if the entry has no heme iron."""
    import gemmi

    path = fetch_cif(pdb_id)
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if len(st) == 0:
        return []
    md = st[0]
    resolution = float(st.resolution) if st.resolution else float("nan")
    method = (st.raw_remarks and "") or ""
    try:
        method = st.get_info("_exptl.method") or ""
    except Exception:
        method = ""

    rows = []
    for chain in md:
        cid = chain.name
        heme_xyz, heme_atom, hets, cys_sg = [], [], [], []
        for res in chain:
            rname = res.name.strip().upper()
            if rname in HEME_ALIASES:
                for at in res:
                    heme_xyz.append([at.pos.x, at.pos.y, at.pos.z])
                    heme_atom.append(at.name.strip())
            elif rname == "CYS":
                for at in res:
                    if at.name.strip() == "SG":
                        cys_sg.append((res.seqid.num, np.array([at.pos.x, at.pos.y, at.pos.z])))
            elif rname not in IGNORE_HET:
                info = gemmi.find_tabulated_residue(rname)
                if info and info.is_amino_acid():
                    continue
                at_xyz = [[a.pos.x, a.pos.y, a.pos.z] for a in res]
                at_el = [a.element.name for a in res]
                at_nm = [a.name.strip() for a in res]
                if at_xyz:
                    hets.append((rname, res.seqid.num, np.asarray(at_xyz, float), at_el, at_nm))

        if not heme_xyz:
            continue
        fe, normal = heme_frame(np.asarray(heme_xyz, float), heme_atom)
        if fe is None:
            continue

        # proximal thiolate: the Cys SG nearest the iron
        sg, sg_res, d_sg = None, None, float("nan")
        if cys_sg:
            dists = [float(np.linalg.norm(p - fe)) for _n, p in cys_sg]
            j = int(np.argmin(dists))
            if dists[j] < 4.0:
                sg_res, sg, d_sg = cys_sg[j][0], cys_sg[j][1], dists[j]

        # orient the normal so +normal points to the DISTAL face (away from the thiolate)
        if sg is not None and normal is not None:
            if np.dot(sg - fe, normal) > 0:
                normal = -normal

        if not hets:
            rows.append(dict(pdb_id=pdb_id, chain=cid, resolution=resolution, method=method,
                             ligand="", n_atoms=0, min_fe_dist=np.nan, donor_element="",
                             s_fe_l_angle=np.nan, distal_side=np.nan, fe_out_of_plane=np.nan,
                             binding_class="apo", axial_cys=sg_res, fe_sg_dist=d_sg))
            continue

        for rname, _seq, xyz, els, _nms in hets:
            d = np.linalg.norm(xyz - fe, axis=1)
            k = int(np.argmin(d))
            dmin = float(d[k])
            donor = els[k]
            # angle S(Cys)-Fe-Ligand_atom: 180 deg means trans to the thiolate
            ang = float("nan")
            if sg is not None:
                v1, v2 = sg - fe, xyz[k] - fe
                c = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
                ang = float(np.degrees(np.arccos(np.clip(c, -1, 1))))
            # which side of the porphyrin plane is the closest ligand atom on
            side = float(np.dot(xyz[k] - fe, normal)) if normal is not None else float("nan")
            # iron displacement out of the pyrrole-N plane, signed toward distal
            oop = float("nan")
            if normal is not None:
                npos = np.asarray([p for p, a in zip(heme_xyz, heme_atom)
                                   if a.upper() in ("NA", "NB", "NC", "ND")], float)
                if len(npos) >= 3:
                    oop = float(np.dot(fe - npos.mean(0), normal))

            if dmin <= COORD_MAX and donor in ("N", "O", "S", "C"):
                cls = "type_II_coordinated"
            elif dmin <= OVERHEAD_MAX:
                cls = "type_I_active_site"
            else:
                cls = "peripheral"
            rows.append(dict(pdb_id=pdb_id, chain=cid, resolution=resolution, method=method,
                             ligand=rname, n_atoms=int(len(xyz)), min_fe_dist=dmin,
                             donor_element=donor, s_fe_l_angle=ang, distal_side=side,
                             fe_out_of_plane=oop, binding_class=cls,
                             axial_cys=sg_res, fe_sg_dist=d_sg))
    return rows


def main(limit: int | None = None) -> None:
    ids = rcsb_holo_structures("cyp3a4")
    if limit:
        ids = ids[:limit]
    print(f"CYP3A4 entries from RCSB: {len(ids)}", flush=True)

    all_rows, failed = [], []
    for i, pid in enumerate(ids, 1):
        try:
            r = analyse_entry(pid)
            all_rows.extend(r)
        except Exception as exc:  # keep going; report the tail honestly at the end
            failed.append((pid, f"{type(exc).__name__}: {exc}"))
        if i % 10 == 0:
            print(f"  {i}/{len(ids)}  rows={len(all_rows)}  failed={len(failed)}", flush=True)

    df = pd.DataFrame(all_rows)
    out = DATA_PROCESSED / "cyp3a4_reference_geometry.parquet"
    df.to_parquet(out, index=False)

    lig = df[df.ligand != ""]
    t2 = lig[lig.binding_class == "type_II_coordinated"]
    t1 = lig[lig.binding_class == "type_I_active_site"]

    def pct(s, qs=(5, 25, 50, 75, 95)):
        s = s.dropna()
        return {f"p{q}": round(float(np.percentile(s, q)), 3) for q in qs} if len(s) else {}

    summary = {
        "generated": "2026-09-09",
        "n_entries_queried": len(ids),
        "n_entries_parsed": int(df.pdb_id.nunique()),
        "n_failed": len(failed),
        "failed": failed[:20],
        "n_chain_ligand_rows": int(len(df)),
        "binding_class_counts": lig.binding_class.value_counts().to_dict(),
        "n_apo_chains": int((df.binding_class == "apo").sum()),
        "fe_sg_dist": pct(df.fe_sg_dist),
        "type_II": {
            "n": int(len(t2)),
            "n_unique_ligands": int(t2.ligand.nunique()),
            "donor_elements": t2.donor_element.value_counts().to_dict(),
            "fe_donor_dist": pct(t2.min_fe_dist),
            "s_fe_ligand_angle": pct(t2.s_fe_l_angle),
            "fe_out_of_plane": pct(t2.fe_out_of_plane),
        },
        "type_I": {
            "n": int(len(t1)),
            "n_unique_ligands": int(t1.ligand.nunique()),
            "min_fe_dist": pct(t1.min_fe_dist),
            "s_fe_ligand_angle": pct(t1.s_fe_l_angle),
        },
        "n_peripheral": int((lig.binding_class == "peripheral").sum()),
        "distal_side_negative_fraction": (
            round(float((lig.distal_side < 0).mean()), 4) if len(lig) else None),
        "parquet": str(out),
    }
    (DATA_PROCESSED / "cyp3a4_reference_geometry.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main(limit=int(sys.argv[1]) if len(sys.argv) > 1 else None)
