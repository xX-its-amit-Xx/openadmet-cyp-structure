"""WHICH side chains are wrong on CYP3A4, which ones exclude the true pose, and how far
they would have to move to admit it.

FINDING 027 established that the crystal ligand, dropped unmoved into the co-folded
protein, makes a median closest heavy-atom contact of 1.42 A and clashes below 2.2 A on
71.2% of poses, while 0 of 87 deposited complexes clash in their own protein. The
backbone is right (FINDING 024: pocket CA 0.73 A) and the side chains are moulded around
whatever orientation the model chose. This script asks the four follow-up questions with
numbers, and it is a DIAGNOSIS: nothing here is or may become a selection feature, because
every quantity below is computed against the crystal.

    python scripts/structure/side_chain_diagnosis.py measure --workers 12   # stage 1
    python scripts/structure/side_chain_diagnosis.py analyse                # Q1-Q3
    python scripts/structure/side_chain_diagnosis.py repack  --workers 12   # Q4

Everything pre-declared, before any number was read, is collected in `PREREG` below and
echoed into the output JSON, so the record shows the thresholds were fixed first. The
project's own failure mode is picking a sign or a cutoff after seeing a scored result
(`loo-sign-selection-fakes-negatives`): the expected signs for every correlation in Q3 are
written down in `PREREG` too.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

POOL = Path(os.environ.get("OE_POOL", "D:/cyp_scratch/val87b_unsteered"))
SCRATCH = Path(os.environ.get(
    "SCD_SCRATCH",
    r"C:\tb\tmp\1\claude\D--Users-ashenoy00000--windsurf-OpenADMET-cyp-structure"
    r"\bd288271-2aa5-4ed4-a343-ea31e5ad8c13\scratchpad\side_chain_diagnosis"))

PREREG = {
    "contact_radius_A": 6.0,        # a residue "lines the pocket" for a pair at this cutoff
    "align_radius_A": 8.0,          # pocket-CA superposition set, identical to FINDING 027
    "clash_cut_A": 2.2,             # ligand-protein, FINDING 027, validated on 87 crystals
    "self_clash_cut_A": 2.6,        # repacked side chain vs the REST of the protein
    "chi_off_deg": 40.0,            # "chi1 is wrong" threshold for the Q1 headline
    "core_pocket_min_pairs": 20,    # a residue enters the per-residue table at >= 20 of 87
    "chi_grid_deg": 10.0,           # Q4 scan resolution on chi1 and chi2
    "chi3_grid_deg": 20.0,          # Q4 escalation resolution on chi3, only when 1x2 fails
    "expected_sign": {
        "rho(sidechain_error, lddt_pli)": "negative",
        "rho(sidechain_error, bisy_rmsd)": "positive",
    },
    "note": "signs fixed before any correlation was computed; no per-fold sign choice",
}

# --------------------------------------------------------------------------
# side-chain chemistry
# --------------------------------------------------------------------------

SC_ATOMS = {
    "ALA": ["CB"],
    "ARG": ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["CB", "CG", "OD1", "ND2"],
    "ASP": ["CB", "CG", "OD1", "OD2"],
    "CYS": ["CB", "SG"],
    "GLN": ["CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["CB", "CG", "CD", "OE1", "OE2"],
    "GLY": [],
    "HIS": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["CB", "CG1", "CG2", "CD1"],
    "LEU": ["CB", "CG", "CD1", "CD2"],
    "LYS": ["CB", "CG", "CD", "CE", "NZ"],
    "MET": ["CB", "CG", "SD", "CE"],
    "PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["CB", "CG", "CD"],
    "SER": ["CB", "OG"],
    "THR": ["CB", "OG1", "CG2"],
    "TRP": ["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL": ["CB", "CG1", "CG2"],
}

CHI1_G = {"ARG": "CG", "ASN": "CG", "ASP": "CG", "CYS": "SG", "GLN": "CG", "GLU": "CG",
          "HIS": "CG", "ILE": "CG1", "LEU": "CG", "LYS": "CG", "MET": "CG", "PHE": "CG",
          "PRO": "CG", "SER": "OG", "THR": "OG1", "TRP": "CG", "TYR": "CG", "VAL": "CG1"}
CHI2_D = {"ARG": "CD", "ASN": "OD1", "ASP": "OD1", "GLN": "CD", "GLU": "CD", "HIS": "ND1",
          "ILE": "CD1", "LEU": "CD1", "LYS": "CD", "MET": "SD", "PHE": "CD1", "PRO": "CD",
          "TRP": "CD1", "TYR": "CD1"}
CHI3_E = {"ARG": "NE", "GLN": "OE1", "GLU": "OE1", "LYS": "CE", "MET": "CE"}

# chi angles that are degenerate modulo 180 because the two distal branches are
# chemically identical. NOTE the brief said "Phe/Tyr/Asp/Glu chi2"; that is right for
# Phe/Tyr/Asp and WRONG for Glu, whose degenerate torsion is chi3 (CB-CG-CD-OE1), not
# chi2 (CA-CB-CG-CD). Implemented correctly, and Glu chi2 is treated as non-degenerate.
CHI_MOD180 = {("PHE", 2), ("TYR", 2), ("ASP", 2), ("GLU", 3), ("ARG", 5)}

# atom-label swaps that leave the molecule unchanged; used only for heavy-atom RMSD
SWAPS = {
    "ASP": [("OD1", "OD2")],
    "GLU": [("OE1", "OE2")],
    "PHE": [("CD1", "CD2"), ("CE1", "CE2")],
    "TYR": [("CD1", "CD2"), ("CE1", "CE2")],
    "ARG": [("NH1", "NH2")],
    "VAL": [("CG1", "CG2")],
    "LEU": [("CD1", "CD2")],
}

# branch pairs whose LABELS are a prochiral convention rather than chemistry. Both
# structures are canonicalised by the SAME geometric rule below, which removes naming
# swaps while preserving the genuine gauche+/trans/gauche- distinction that a plain
# "minimise over the swap" would destroy.
PROCHIRAL = {"VAL": ("CB", "CA", "CG1", "CG2"), "LEU": ("CG", "CB", "CD1", "CD2")}

# rotation axes and the atoms distal to each
ROT = {
    "ARG": [("CA", "CB", ["CG", "CD", "NE", "CZ", "NH1", "NH2"]),
            ("CB", "CG", ["CD", "NE", "CZ", "NH1", "NH2"]),
            ("CG", "CD", ["NE", "CZ", "NH1", "NH2"])],
    "LYS": [("CA", "CB", ["CG", "CD", "CE", "NZ"]), ("CB", "CG", ["CD", "CE", "NZ"]),
            ("CG", "CD", ["CE", "NZ"])],
    "MET": [("CA", "CB", ["CG", "SD", "CE"]), ("CB", "CG", ["SD", "CE"]),
            ("CG", "SD", ["CE"])],
    "GLU": [("CA", "CB", ["CG", "CD", "OE1", "OE2"]), ("CB", "CG", ["CD", "OE1", "OE2"]),
            ("CG", "CD", ["OE1", "OE2"])],
    "GLN": [("CA", "CB", ["CG", "CD", "OE1", "NE2"]), ("CB", "CG", ["CD", "OE1", "NE2"]),
            ("CG", "CD", ["OE1", "NE2"])],
    "ASP": [("CA", "CB", ["CG", "OD1", "OD2"]), ("CB", "CG", ["OD1", "OD2"])],
    "ASN": [("CA", "CB", ["CG", "OD1", "ND2"]), ("CB", "CG", ["OD1", "ND2"])],
    "PHE": [("CA", "CB", ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]),
            ("CB", "CG", ["CD1", "CD2", "CE1", "CE2", "CZ"])],
    "TYR": [("CA", "CB", ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]),
            ("CB", "CG", ["CD1", "CD2", "CE1", "CE2", "CZ", "OH"])],
    "TRP": [("CA", "CB", ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]),
            ("CB", "CG", ["CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"])],
    "HIS": [("CA", "CB", ["CG", "ND1", "CD2", "CE1", "NE2"]),
            ("CB", "CG", ["ND1", "CD2", "CE1", "NE2"])],
    "LEU": [("CA", "CB", ["CG", "CD1", "CD2"]), ("CB", "CG", ["CD1", "CD2"])],
    "ILE": [("CA", "CB", ["CG1", "CG2", "CD1"]), ("CB", "CG1", ["CD1"])],
    "VAL": [("CA", "CB", ["CG1", "CG2"])],
    "THR": [("CA", "CB", ["OG1", "CG2"])],
    "SER": [("CA", "CB", ["OG"])],
    "CYS": [("CA", "CB", ["SG"])],
}
BACKBONE = {"N", "CA", "C", "O", "OXT"}

# The ligand-free alignment frame and the plastic span, both declared in
# `cypstruct.targets` long before this experiment and imported rather than re-chosen.
from cypstruct.targets import CYP3A4_FG_LOOP as FG          # noqa: E402
from cypstruct.targets import CYP3A4_RIGID_CORE             # noqa: E402

CORE = set(CYP3A4_RIGID_CORE)


def dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-12)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def angdiff(a: float, b: float, mod180: bool = False) -> float:
    """Signed circular difference in degrees, wrapped to (-180,180] or (-90,90]."""
    if not np.isfinite(a) or not np.isfinite(b):
        return float("nan")
    p = 180.0 if mod180 else 360.0
    return float(((a - b + p / 2) % p) - p / 2)


def circ_sd(angles) -> float:
    """Circular standard deviation in degrees; nan if fewer than 2 finite values."""
    a = np.asarray([x for x in angles if np.isfinite(x)], float)
    if len(a) < 2:
        return float("nan")
    r = np.hypot(np.cos(np.radians(a)).mean(), np.sin(np.radians(a)).mean())
    r = min(max(r, 1e-9), 1.0)
    return float(np.degrees(np.sqrt(-2.0 * np.log(r))))


def circ_mean(angles) -> float:
    a = np.asarray([x for x in angles if np.isfinite(x)], float)
    if not len(a):
        return float("nan")
    return float(np.degrees(np.arctan2(np.sin(np.radians(a)).mean(),
                                       np.cos(np.radians(a)).mean())))


def residue_atoms(cx, resnum: int) -> dict[str, np.ndarray]:
    return {a: cx.prot_xyz[i] for i, (_c, r, a) in enumerate(cx.prot_key) if r == resnum}


def canonicalise_prochiral(at: dict[str, np.ndarray], rname: str) -> bool:
    """Relabel VAL CG1/CG2 and LEU CD1/CD2 by a fixed geometric sign, in place.

    The two branches are chemically identical, so which one is called `1` is a
    convention. Applying ONE rule to both the crystal and the prediction removes a
    spurious ~120 deg chi1 error for Val (and ~120 deg chi2 for Leu) without collapsing
    the three genuine rotamers into one, which is what "minimise over the swap" would do.
    Returns True when a swap was applied, so the count can be reported.
    """
    spec = PROCHIRAL.get(rname)
    if not spec:
        return False
    cen, base, a1, a2 = spec
    if not all(k in at for k in spec):
        return False
    s = np.dot(np.cross(at[a1] - at[cen], at[a2] - at[cen]), at[base] - at[cen])
    if s < 0:
        at[a1], at[a2] = at[a2].copy(), at[a1].copy()
        return True
    return False


def chis(at: dict[str, np.ndarray], rname: str) -> tuple[float, float]:
    c1 = c2 = float("nan")
    g = CHI1_G.get(rname)
    if g and all(k in at for k in ("N", "CA", "CB", g)):
        c1 = dihedral(at["N"], at["CA"], at["CB"], at[g])
    d = CHI2_D.get(rname)
    if g and d and all(k in at for k in ("CA", "CB", g, d)):
        c2 = dihedral(at["CA"], at["CB"], at[g], at[d])
    return c1, c2


def sc_rmsd(a: dict, b: dict, rname: str, beyond_cb: bool = True) -> float:
    """Heavy-atom RMSD of one side chain, minimised over label swaps."""
    names = [n for n in SC_ATOMS.get(rname, []) if not beyond_cb or n != "CB"]
    names = [n for n in names if n in a and n in b]
    if not names:
        return float("nan")
    best = None
    for k in range(1 << len(SWAPS.get(rname, []))):
        bb = dict(b)
        for j, (x, y) in enumerate(SWAPS.get(rname, [])):
            if (k >> j) & 1 and x in bb and y in bb:
                bb[x], bb[y] = bb[y], bb[x]
        v = np.sqrt(np.mean([np.sum((a[n] - bb[n]) ** 2) for n in names]))
        best = v if best is None else min(best, v)
    return float(best)


def rot_about(pts: np.ndarray, p0: np.ndarray, p1: np.ndarray,
              deg: np.ndarray) -> np.ndarray:
    """Rotate (n,3) `pts` about the axis p0->p1 by each angle in `deg`. -> (len(deg),n,3)"""
    k = p1 - p0
    k = k / (np.linalg.norm(k) + 1e-12)
    v = pts - p0
    th = np.radians(deg)[:, None, None]
    kv = np.cross(np.broadcast_to(k, v.shape), v)
    kkv = (v @ k)[:, None] * k
    return p0 + v * np.cos(th) + kv * np.sin(th) + kkv * (1 - np.cos(th))


# --------------------------------------------------------------------------
# shared loading (same helpers FINDING 027 used, so the frame is the same one)
# --------------------------------------------------------------------------

def load_reference(pdb_id: str, ligand_code: str):
    import gemmi
    from cypstruct import pose as P
    from cypstruct.targets import fetch_cif
    cif = fetch_cif(pdb_id)
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    for chain in st[0]:
        if any(r.name.strip().upper() == ligand_code.upper() for r in chain):
            return P.load_structure(cif, ligand_code=ligand_code,
                                    assembly_chain=chain.name), cif, chain.name
    return None, None, None


def renumber_to_reference(pred, ref):
    """FINDING 021: residue-NAME-agreement offset. Never array position."""
    pm = {num: name for (_c, num), name in pred.prot_res.items()}
    rm = {num: name for (_c, num), name in ref.prot_res.items()}
    if len(pm) < 30 or len(rm) < 30:
        return pred, 0, 0.0
    best, best_id = 0, -1.0
    for off in range(int(min(rm) - max(pm)), int(max(rm) - min(pm)) + 1):
        shared = [n for n in pm if n + off in rm]
        if len(shared) < 30:
            continue
        ident = sum(pm[n] == rm[n + off] for n in shared) / len(shared)
        if ident > best_id:
            best, best_id = off, ident
    if best_id < 0:
        return pred, 0, 0.0
    if best != 0:
        pred.prot_key = [(c, n + best, a) for (c, n, a) in pred.prot_key]
        pred.prot_res = {(c, n + best): v for (c, n), v in pred.prot_res.items()}
    return pred, best, best_id


def crystal_quality(cif_path, chain: str, resnums: set[int]) -> dict:
    """Occupancy / B / completeness of the deposited SIDE CHAINS, plus the resolution.

    A residue that is unmodelled or at zero occupancy cannot support a conclusion about
    its rotamer, so the count is reported rather than assumed to be zero.
    """
    import gemmi
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    res = float("nan")
    try:
        res = float(st.resolution)
    except Exception:
        pass
    out = {}
    for ch in st[0]:
        if ch.name != chain:
            continue
        for r in ch:
            n = r.seqid.num
            if n not in resnums:
                continue
            rn = r.name.strip().upper()
            exp = SC_ATOMS.get(rn)
            if exp is None:
                continue
            seen, occ, b, alt = set(), [], [], False
            for at in r:
                if at.element == gemmi.Element("H"):
                    continue
                nm = at.name.strip()
                if at.altloc not in ("", "\0", " "):
                    alt = True
                if nm in exp:
                    seen.add(nm)
                    occ.append(float(at.occ))
                    b.append(float(at.b_iso))
            out[n] = dict(resname=rn, n_sc=len(seen), n_sc_expected=len(exp),
                          complete=bool(len(exp) == 0 or len(seen) == len(exp)),
                          occ_min=float(min(occ)) if occ else float("nan"),
                          b_mean=float(np.mean(b)) if b else float("nan"),
                          altloc=alt)
    return {"resolution": res, "residues": out, "method": str(st.spacegroup_hm)}


# --------------------------------------------------------------------------
# stage 1 — one ligand
# --------------------------------------------------------------------------

def measure_ligand(args) -> dict:
    lig_id, pdb, smiles = args
    from scipy.spatial import cKDTree
    from cypstruct import pose as P

    out = SCRATCH / f"{lig_id}.json"
    if out.exists():
        return {"ligand": lig_id, "status": "cached"}

    ref, cif_path, chain = load_reference(pdb, lig_id)
    if ref is None or len(ref.lig_xyz) == 0:
        return {"ligand": lig_id, "status": "no-reference"}
    job = POOL / f"{lig_id}__unsteered__s1"
    cifs = sorted(job.glob("input_model_*.cif"),
                  key=lambda p: int(p.stem.rsplit("_", 1)[1]))
    if not cifs:
        return {"ligand": lig_id, "status": "no-poses"}

    contact = P.pocket_residues_from_structure(ref, radius=PREREG["contact_radius_A"])
    alignset = P.pocket_residues_from_structure(ref, radius=PREREG["align_radius_A"])
    qual = crystal_quality(cif_path, chain, set(alignset))

    # crystal side chains, canonicalised
    ref_at, ref_chi, n_swap_ref = {}, {}, 0
    for rn in alignset:
        rname = ref.prot_res.get((ref.lig_chain, rn)) or \
            next((v for (c, n), v in ref.prot_res.items() if n == rn), None)
        if rname is None:
            continue
        at = residue_atoms(ref, rn)
        n_swap_ref += int(canonicalise_prochiral(at, rname))
        ref_at[rn] = (rname, at)
        ref_chi[rn] = chis(at, rname)

    # crystal ligand distance to each crystal residue (the CONTROL for Q2)
    ctrl_rows = []
    hm = np.asarray(ref.heme_xyz, float)
    he = [e.upper() for e in ref.heme_elem]
    hn = hm[[i for i, e in enumerate(he) if e != "FE"]] if len(hm) else hm
    envr = np.vstack([ref.prot_xyz, hn]) if len(hn) else ref.prot_xyz
    d_native = cKDTree(envr).query(ref.lig_xyz, k=1)[0].min()
    dr = np.linalg.norm(ref.prot_xyz[:, None, :] - ref.lig_xyz[None, :, :], axis=2).min(1)
    j = int(np.argmin(dr))
    native_worst = dict(resnum=int(ref.prot_key[j][1]), atom=ref.prot_key[j][2],
                        dist=float(dr[j]))

    per_res_ref = {}
    for rn in alignset:
        idx = [i for i, (_c, r, a) in enumerate(ref.prot_key) if r == rn]
        if not idx:
            continue
        per_res_ref[rn] = float(dr[idx].min())

    poses = []
    n_swap_mod, n_unmapped = 0, 0
    for cif in cifs:
        m = P.load_structure(cif)
        m, off, ident = renumber_to_reference(m, ref)
        perm = P.best_ligand_mapping(smiles, m, ref)
        mapped = perm is not None
        if perm is None:
            n_unmapped += 1
            continue
        # crystal -> model, on pocket CA matched by residue NUMBER.
        # This is FINDING 027's frame, reproduced exactly so Q2 is comparable to it.
        mn = {k: v for (_c, k), v in m.ca().items()}
        rn_ca = {k: v for (_c, k), v in ref.ca().items()}
        sh = sorted(set(mn) & set(rn_ca) & set(alignset))
        Rk, tk, fit_pocket = P.kabsch(np.array([rn_ca[k] for k in sh]),
                                      np.array([mn[k] for k in sh]))
        Ltrue = ref.lig_xyz[perm] @ Rk.T + tk

        # The LIGAND-FREE frame: CYP3A4_RIGID_CORE, declared in targets.py long before
        # this experiment as "the heme-proximal half, which does not move between holo
        # forms". Used for every backbone/side-chain error number, because a pocket-CA
        # fit is dragged by whichever pocket residues moved.
        shc = sorted(set(mn) & set(rn_ca) & CORE)
        Rm, tm, fit_core = P.kabsch(np.array([mn[k] for k in shc]),
                                    np.array([rn_ca[k] for k in shc]))
        Rc, tc, _ = P.kabsch(np.array([rn_ca[k] for k in shc]),
                             np.array([mn[k] for k in shc]))
        Ltrue_core = ref.lig_xyz[perm] @ Rc.T + tc
        ca_dev = {k: float(np.linalg.norm(mn[k] @ Rm.T + tm - rn_ca[k]))
                  for k in set(mn) & set(rn_ca)}

        def _rms(v):
            v = [x for x in v if np.isfinite(x)]
            return float(np.sqrt(np.mean(np.square(v)))) if v else float("nan")

        inpk = [k for k in ca_dev if k in set(alignset)]
        ca_stats = dict(
            fit_core=float(fit_core), fit_pocket=float(fit_pocket),
            ca_core=_rms([ca_dev[k] for k in ca_dev if k in CORE]),
            ca_global=_rms(list(ca_dev.values())),
            ca_global_med=float(np.median(list(ca_dev.values()))),
            ca_pocket=_rms([ca_dev[k] for k in inpk]),
            ca_pocket_med=float(np.median([ca_dev[k] for k in inpk])),
            ca_pocket_fg=_rms([ca_dev[k] for k in inpk if FG[0] <= k <= FG[1]]),
            ca_pocket_nonfg=_rms([ca_dev[k] for k in inpk if not FG[0] <= k <= FG[1]]),
            n_pocket_fg=int(sum(1 for k in inpk if FG[0] <= k <= FG[1])),
            n_pocket=len(inpk))

        hm = np.asarray(m.heme_xyz, float)
        he = [e.upper() for e in m.heme_elem]
        hn = hm[[i for i, e in enumerate(he) if e != "FE"]] if len(hm) else hm
        env = np.vstack([m.prot_xyz, hn]) if len(hn) else m.prot_xyz
        env_lbl = [(k[1], k[2]) for k in m.prot_key] + [(-1, "HEME")] * len(hn)
        t = cKDTree(env)
        dq, iq = t.query(Ltrue, k=1)
        kbest = int(np.argmin(dq))
        worst = dict(dist=float(dq[kbest]), resnum=int(env_lbl[iq[kbest]][0]),
                     atom=env_lbl[iq[kbest]][1])
        worst["resname"] = next((v for (c, r), v in m.prot_res.items()
                                 if r == worst["resnum"]), "HEME")
        worst["sidechain"] = worst["atom"] not in BACKBONE and worst["atom"] != "HEME"
        worst["in_fg"] = bool(FG[0] <= worst["resnum"] <= FG[1])
        clash_core = float(t.query(Ltrue_core, k=1)[0].min())
        # per-residue contact of the crystal ligand inside the model protein
        dm = np.linalg.norm(m.prot_xyz[:, None, :] - Ltrue[None, :, :], axis=2).min(1)
        pred_min = float(cKDTree(env).query(np.asarray(m.lig_xyz, float), k=1)[0].min())

        pr = {}
        for rnum in alignset:
            idx = [i for i, (_c, r, a) in enumerate(m.prot_key) if r == rnum]
            if not idx:
                continue
            sc = [i for i in idx if m.prot_key[i][2] not in BACKBONE]
            pr[rnum] = (float(dm[idx].min()),
                        float(dm[sc].min()) if sc else float("nan"))

        # side-chain geometry, this pose
        chi_rows = []
        for rnum, (rname, rat) in ref_at.items():
            mat = residue_atoms(m, rnum)
            if not mat:
                continue
            mrname = next((v for (c, r), v in m.prot_res.items() if r == rnum), None)
            if mrname != rname:
                continue
            n_swap_mod += int(canonicalise_prochiral(mat, rname))
            c1, c2 = chis(mat, rname)
            r1, r2 = ref_chi[rnum]
            d1 = angdiff(c1, r1, (rname, 1) in CHI_MOD180)
            d2 = angdiff(c2, r2, (rname, 2) in CHI_MOD180)
            # LOCAL frame: superpose this residue's own N,CA,C, then side chain only
            loc = float("nan")
            if all(k in mat for k in ("N", "CA", "C")) and \
               all(k in rat for k in ("N", "CA", "C")):
                A = np.array([mat[k] for k in ("N", "CA", "C")])
                B = np.array([rat[k] for k in ("N", "CA", "C")])
                Rl, tl, _ = P.kabsch(A, B)
                loc = sc_rmsd({k: v @ Rl.T + tl for k, v in mat.items()}, rat, rname)
            glob = sc_rmsd({k: v @ Rm.T + tm for k, v in mat.items()}, rat, rname)
            ca_err = ca_dev.get(int(rnum), float("nan"))
            chi_rows.append(dict(resnum=int(rnum), resname=rname,
                                 chi1_mod=c1, chi2_mod=c2, chi1_ref=r1, chi2_ref=r2,
                                 dchi1=d1, dchi2=d2, sc_rmsd_local=loc,
                                 sc_rmsd_coreframe=glob, ca_err=ca_err,
                                 in_fg=bool(FG[0] <= rnum <= FG[1]),
                                 dist_true_lig=pr.get(int(rnum), (np.nan, np.nan))[0],
                                 dist_true_lig_sc=pr.get(int(rnum),
                                                         (np.nan, np.nan))[1]))
        poses.append(dict(sample=cif.stem, offset=int(off), identity=float(ident),
                          mapped=bool(mapped), worst=worst, clash_pred=pred_min,
                          clash_true_core_frame=clash_core,
                          n_align=len(sh), ca=ca_stats, chi=chi_rows))

    rec = dict(ligand=lig_id, pdb=pdb, chain=chain,
               contact=[int(x) for x in contact], alignset=[int(x) for x in alignset],
               crystal=dict(resolution=qual["resolution"], residues={
                   str(k): v for k, v in qual["residues"].items()},
                   native_min_contact=float(d_native), native_worst=native_worst,
                   per_res_dist={str(k): v for k, v in per_res_ref.items()},
                   chi={str(k): [ref_chi[k][0], ref_chi[k][1]] for k in ref_chi},
                   resname={str(k): ref_at[k][0] for k in ref_at}),
               prochiral_swaps_ref=int(n_swap_ref), prochiral_swaps_model=int(n_swap_mod),
               n_unmapped=int(n_unmapped), n_cifs=len(cifs), poses=poses)
    out.write_text(json.dumps(rec))
    return {"ligand": lig_id, "status": "ok", "poses": len(poses)}


def cmd_measure(workers: int, limit: int | None) -> None:
    import pandas as pd
    SCRATCH.mkdir(parents=True, exist_ok=True)
    lig = pd.read_csv(REPO / "data" / "processed" / "validation_ligands.csv")
    have = {p.name.split("__")[0] for p in POOL.iterdir() if p.is_dir()}
    jobs = [(r.id, r.pdb, r.smiles) for r in lig.itertuples() if r.id in have]
    if limit:
        jobs = jobs[:limit]
    todo = [j for j in jobs if not (SCRATCH / f"{j[0]}.json").exists()]
    print(f"{len(jobs)} ligands, {len(todo)} to do, {workers} workers", flush=True)
    if not todo:
        return
    if workers <= 1:
        for j in todo:
            print(measure_ligand(j), flush=True)
        return
    import multiprocessing as mp
    with mp.Pool(workers) as pool:
        for res in pool.imap_unordered(measure_ligand, todo):
            print(res, flush=True)


# --------------------------------------------------------------------------
# stage 2 — Q1, Q2, Q3
# --------------------------------------------------------------------------

def _load_all():
    return [json.loads(f.read_text()) for f in sorted(SCRATCH.glob("*.json"))
            if not f.name.startswith("repack_")]


def cmd_analyse() -> None:
    import pandas as pd
    from scipy import stats

    recs = _load_all()
    chi_rows, pose_rows, cry_rows = [], [], []
    for rec in recs:
        lig, contact = rec["ligand"], set(rec["contact"])
        for p in rec["poses"]:
            w = p["worst"]
            pose_rows.append(dict(ligand=lig, sample=p["sample"], offset=p["offset"],
                                  identity=p["identity"], clash_pred=p["clash_pred"],
                                  clash_core=p["clash_true_core_frame"],
                                  worst_dist=w["dist"], worst_resnum=w["resnum"],
                                  worst_resname=w["resname"], worst_atom=w["atom"],
                                  worst_sidechain=w["sidechain"], worst_fg=w["in_fg"],
                                  **p["ca"]))
            for c in p["chi"]:
                c2 = dict(c, ligand=lig, sample=p["sample"],
                          contact=c["resnum"] in contact)
                chi_rows.append(c2)
        q = rec["crystal"]["residues"]
        for k, v in q.items():
            cry_rows.append(dict(ligand=lig, resnum=int(k), contact=int(k) in contact,
                                 resolution=rec["crystal"]["resolution"], **v))

    chi = pd.DataFrame(chi_rows)
    pos = pd.DataFrame(pose_rows)
    cry = pd.DataFrame(cry_rows)

    # ---- control block: FINDING 021 and the crystal-quality audit ---------
    n_pairs = len(recs)
    cq = cry[cry.contact]
    controls = {
        "n_pairs": n_pairs,
        "n_poses": int(len(pos)),
        "C_offsets": {str(k): int(v) for k, v in pos.offset.value_counts().items()},
        "C_min_residue_name_identity": float(pos.identity.min()),
        "C_poses_dropped_unmapped_ligand": int(sum(r["n_unmapped"] for r in recs)),
        "C_cifs_available": int(sum(r["n_cifs"] for r in recs)),
        "C_prochiral_swaps_crystal": int(sum(r["prochiral_swaps_ref"] for r in recs)),
        "C_prochiral_swaps_model": int(sum(r["prochiral_swaps_model"] for r in recs)),
        "crystal_quality": {
            "n_contact_residue_observations": int(len(cq)),
            "incomplete_side_chains": int((~cq.complete).sum()),
            "frac_incomplete": float((~cq.complete).mean()),
            "occupancy_below_1": int((cq.occ_min < 1.0).sum()),
            "occupancy_zero": int((cq.occ_min <= 0.0).sum()),
            "altloc_present": int(cq.altloc.sum()),
            "b_mean_median": float(cq.b_mean.median()),
            "b_mean_p90": float(cq.b_mean.quantile(0.90)),
            "resolution_median": float(cry.groupby("ligand").resolution.first().median()),
            "resolution_worst": float(cry.groupby("ligand").resolution.first().max()),
        },
    }
    # residues whose side chain is missing in >= 25% of the crystals that contact them
    miss = (cq.groupby("resnum")
              .agg(n=("complete", "size"), frac_incomplete=("complete",
                                                            lambda s: 1 - s.mean()),
                   b=("b_mean", "median"))
              .query("n >= @PREREG['core_pocket_min_pairs']")
              .sort_values("frac_incomplete", ascending=False))
    controls["crystal_quality"]["worst_modelled_core_residues"] = [
        {"resnum": int(i), "n": int(r.n), "frac_incomplete": round(float(r.frac_incomplete), 3),
         "median_B": round(float(r.b), 1)} for i, r in miss.head(8).iterrows()]

    # a residue enters the tables only where the crystal side chain is COMPLETE
    ok = cry[cry.complete][["ligand", "resnum"]].assign(usable=True)
    chi = chi.merge(ok, on=["ligand", "resnum"], how="left")
    chi["usable"] = chi.usable.fillna(False).astype(bool)
    controls["chi_rows_total"] = int(len(chi))
    controls["chi_rows_dropped_incomplete_crystal"] = int((~chi.usable).sum())
    chi = chi[chi.usable].copy()
    controls["chi_rows_used"] = int(len(chi))
    controls["chi1_nan_rows"] = int(chi.dchi1.isna().sum())

    # per-pair collapse: 20 samples -> one number per (ligand, residue)
    def cmean(s):
        return circ_mean(s.values)

    pair = (chi.groupby(["ligand", "resnum", "resname", "contact", "in_fg"])
               .agg(dchi1=("dchi1", lambda s: circ_mean(s.values)),
                    dchi2=("dchi2", lambda s: circ_mean(s.values)),
                    chi1_mod=("chi1_mod", lambda s: circ_mean(s.values)),
                    chi2_mod=("chi2_mod", lambda s: circ_mean(s.values)),
                    chi1_ref=("chi1_ref", "first"), chi2_ref=("chi2_ref", "first"),
                    sc_local=("sc_rmsd_local", "median"),
                    sc_pocket=("sc_rmsd_coreframe", "median"),
                    ca_err=("ca_err", "median"),
                    wl_sd_chi1=("chi1_mod", lambda s: circ_sd(s.values)),
                    dist_true=("dist_true_lig", "median"),
                    dist_true_sc=("dist_true_lig_sc", "median"))
               .reset_index())
    pair["adchi1"] = pair.dchi1.abs()
    pair["adchi2"] = pair.dchi2.abs()

    # ---- Q0: where is the BACKBONE wrong? The frame matters, so measure it -----
    bb = {
        "note": "all CA errors after superposing on CYP3A4_RIGID_CORE, a ligand-free "
                "set declared in targets.py before this experiment",
        "core_fit_rmsd_median_A": float(pos.fit_core.median()),
        "pocket_fit_rmsd_median_A": float(pos.fit_pocket.median()),
        "ca_core_median_A": float(pos.ca_core.median()),
        "ca_global_rms_median_A": float(pos.ca_global.median()),
        "ca_global_per_residue_median_A": float(pos.ca_global_med.median()),
        "ca_pocket_rms_median_A": float(pos.ca_pocket.median()),
        "ca_pocket_per_residue_median_A": float(pos.ca_pocket_med.median()),
        "ca_pocket_nonFG_rms_median_A": float(pos.ca_pocket_nonfg.median()),
        "ca_pocket_FG_rms_median_A": float(pos.ca_pocket_fg.median()),
        "ca_pocket_FG_rms_p90_A": float(pos.ca_pocket_fg.quantile(0.90)),
        "frac_poses_pocket_FG_over_2A": float((pos.ca_pocket_fg > 2.0).mean()),
        "frac_poses_pocket_nonFG_over_2A": float((pos.ca_pocket_nonfg > 2.0).mean()),
        "median_pocket_residues": float(pos.n_pocket.median()),
        "median_pocket_residues_in_FG": float(pos.n_pocket_fg.median()),
    }

    # ---------------- Q1 ---------------------------------------------------
    cp = pair[pair.contact]
    g = cp.groupby(["resnum", "resname"])
    q1 = g.agg(n_pairs=("ligand", "nunique"),
               med_dchi1=("adchi1", "median"),
               q1_dchi1=("adchi1", lambda s: s.quantile(0.25)),
               q3_dchi1=("adchi1", lambda s: s.quantile(0.75)),
               frac_chi1_off=("adchi1", lambda s: float((s > PREREG["chi_off_deg"]).mean())),
               med_dchi2=("adchi2", "median"),
               frac_chi2_off=("adchi2", lambda s: float((s > PREREG["chi_off_deg"]).mean())),
               med_sc_local=("sc_local", "median"),
               med_sc_core=("sc_pocket", "median"),
               med_ca_err=("ca_err", "median"),
               in_fg=("in_fg", "first")).reset_index()
    q1 = q1[q1.n_pairs >= PREREG["core_pocket_min_pairs"]].copy()
    q1 = q1.sort_values("frac_chi1_off", ascending=False)
    q1.to_csv(REPO / "data" / "processed" / "side_chain_q1_per_residue.csv", index=False)

    pooled = {
        "core_pocket_residues": int(len(q1)),
        "median_abs_dchi1_deg": float(cp.adchi1.median()),
        "iqr_abs_dchi1_deg": [float(cp.adchi1.quantile(0.25)),
                              float(cp.adchi1.quantile(0.75))],
        "frac_chi1_off_gt40": float((cp.adchi1 > PREREG["chi_off_deg"]).mean()),
        "median_abs_dchi2_deg": float(cp.adchi2.median()),
        "frac_chi2_off_gt40": float((cp.adchi2 > PREREG["chi_off_deg"]).mean()),
        "median_sc_rmsd_local_A": float(cp.sc_local.median()),
        "median_sc_rmsd_coreframe_A": float(cp.sc_pocket.median()),
        "median_ca_err_A": float(cp.ca_err.median()),
        "FG_split": {
            "median_abs_dchi1_FG": float(cp[cp.in_fg].adchi1.median()),
            "median_abs_dchi1_nonFG": float(cp[~cp.in_fg].adchi1.median()),
            "frac_chi1_off_FG": float((cp[cp.in_fg].adchi1 > PREREG["chi_off_deg"]).mean()),
            "frac_chi1_off_nonFG": float(
                (cp[~cp.in_fg].adchi1 > PREREG["chi_off_deg"]).mean()),
            "median_sc_rmsd_local_FG": float(cp[cp.in_fg].sc_local.median()),
            "median_sc_rmsd_local_nonFG": float(cp[~cp.in_fg].sc_local.median()),
            "n_FG": int(cp.in_fg.sum()), "n_nonFG": int((~cp.in_fg).sum()),
            "mannwhitney_p": float(stats.mannwhitneyu(
                cp[cp.in_fg].adchi1.dropna(), cp[~cp.in_fg].adchi1.dropna()).pvalue)},
    }
    # concentration: how much of the total chi1 error sits in the worst k residues
    order = q1.sort_values("frac_chi1_off", ascending=False)
    tot = float(order.frac_chi1_off.sum())
    pooled["frac_chi1_off_share_top5"] = float(order.frac_chi1_off.head(5).sum() / tot)
    pooled["frac_chi1_off_min"] = float(order.frac_chi1_off.min())
    pooled["frac_chi1_off_max"] = float(order.frac_chi1_off.max())
    pooled["frac_chi1_off_median_across_residues"] = float(order.frac_chi1_off.median())
    # per-PAIR: how many of its contact residues are off by >40 deg
    perpair = (cp.assign(off=cp.adchi1 > PREREG["chi_off_deg"])
                 .groupby("ligand").agg(n=("off", "size"), k=("off", "sum")))
    pooled["per_pair_contact_residues_median"] = float(perpair.n.median())
    pooled["per_pair_chi1_off_count_median"] = float(perpair.k.median())
    pooled["per_pair_chi1_off_frac_median"] = float((perpair.k / perpair.n).median())

    # ---------------- Q2 ---------------------------------------------------
    cut = PREREG["clash_cut_A"]
    natives = [r["crystal"]["native_min_contact"] for r in recs]
    q2h = {
        "control_crystal_ligand_in_own_crystal": {
            "n": len(natives), "median": float(np.median(natives)),
            "min": float(np.min(natives)),
            "frac_below_cut": float(np.mean(np.array(natives) < cut))},
        "crystal_ligand_in_model_protein": {
            "n": int(len(pos)), "median": float(pos.worst_dist.median()),
            "frac_below_cut": float((pos.worst_dist < cut).mean())},
        "predicted_pose_in_own_model_protein": {
            "median": float(pos.clash_pred.median()),
            "frac_below_cut": float((pos.clash_pred < cut).mean())},
        "crystal_ligand_in_model_protein_CORE_frame": {
            "median": float(pos.clash_core.median()),
            "frac_below_cut": float((pos.clash_core < cut).mean())},
        "worst_contact_is_a_side_chain_atom": float(pos.worst_sidechain.mean()),
        "worst_contact_is_the_heme": float((pos.worst_resname == "HEME").mean()),
        "worst_contact_is_in_the_FG_span": float(pos.worst_fg.mean()),
    }
    clashing = pos[pos.worst_dist < cut]
    wa = (clashing.groupby(["worst_resnum", "worst_resname"])
                  .agg(n_poses=("ligand", "size"), n_ligands=("ligand", "nunique"),
                       med_dist=("worst_dist", "median")).reset_index()
                  .sort_values("n_poses", ascending=False))
    wa["frac_of_clashing_poses"] = wa.n_poses / max(len(clashing), 1)
    wa.to_csv(REPO / "data" / "processed" / "side_chain_q2_blockers.csv", index=False)

    # blocking rate per residue: fraction of poses where THIS residue is under the cut
    blk = (chi.assign(blocks=chi.dist_true_lig < cut,
                      blocks_sc=chi.dist_true_lig_sc < cut)
              .groupby(["resnum", "resname"])
              .agg(n_obs=("blocks", "size"),
                   n_lig=("ligand", "nunique"),
                   frac_blocking=("blocks", "mean"),
                   frac_blocking_sc=("blocks_sc", "mean"),
                   med_dist=("dist_true_lig", "median")).reset_index())
    blk = blk[blk.n_lig >= PREREG["core_pocket_min_pairs"]]
    x = q1.merge(blk, on=["resnum", "resname"], how="inner")
    xtab_rho = float(stats.spearmanr(x.frac_chi1_off, x.frac_blocking).statistic)
    xtab_p = float(stats.spearmanr(x.frac_chi1_off, x.frac_blocking).pvalue)
    x = x.sort_values("frac_blocking", ascending=False)
    x.to_csv(REPO / "data" / "processed" / "side_chain_q1q2_crosstab.csv", index=False)
    # No threshold is chosen here on purpose: both rankings are printed in full with
    # BOTH columns attached, so "wrong" and "blocking" can be read against each other
    # without a cutoff that could be tuned after the fact. The full table is in
    # side_chain_q1q2_crosstab.csv.
    def _pair_rows(frame):
        return [{"res": f"{r.resname}{r.resnum}", "fg": bool(r.in_fg),
                 "med_dchi1": round(float(r.med_dchi1), 1)
                 if np.isfinite(r.med_dchi1) else None,
                 "frac_chi1_off": round(float(r.frac_chi1_off), 3),
                 "frac_blocking": round(float(r.frac_blocking), 4),
                 "med_sc_rmsd_local": (round(float(r.med_sc_local), 2)
                                       if np.isfinite(r.med_sc_local) else None),
                 "med_ca_err": round(float(r.med_ca_err), 2)}
                for _, r in frame.iterrows()]

    q2h["crosstab"] = {
        "n_residues": int(len(x)),
        "rho_frac_chi1_off_vs_frac_blocking": xtab_rho, "p": xtab_p,
        "top8_by_chi1_error": _pair_rows(
            x.sort_values("frac_chi1_off", ascending=False).head(8)),
        "top8_by_blocking": _pair_rows(x.head(8)),
    }

    # ---------------- Q3 ---------------------------------------------------
    # per-residue: sd of chi1 ACROSS ligands, model vs crystal
    lig_model = (chi.groupby(["ligand", "resnum", "resname"])
                    .agg(chi1=("chi1_mod", lambda s: circ_mean(s.values)),
                         wl=("chi1_mod", lambda s: circ_sd(s.values))).reset_index())
    lig_cry = pair[["ligand", "resnum", "resname", "chi1_ref"]]
    plast = (lig_model.groupby(["resnum", "resname"])
                      .agg(n_lig=("ligand", "nunique"),
                           sd_model=("chi1", lambda s: circ_sd(s.values)),
                           sd_within_ligand=("wl", "median")).reset_index())
    plc = (lig_cry.groupby(["resnum", "resname"])
                  .agg(sd_crystal=("chi1_ref", lambda s: circ_sd(s.values))).reset_index())
    plast = plast.merge(plc, on=["resnum", "resname"], how="inner")
    plast = plast[plast.n_lig >= PREREG["core_pocket_min_pairs"]].copy()
    n_before = len(plast)
    # residues with no chi1 at all (Ala, Gly) carry no rotamer and are dropped here
    plast = plast.dropna(subset=["sd_model", "sd_crystal"]).copy()
    plast["ratio"] = plast.sd_model / plast.sd_crystal
    plast = plast.sort_values("sd_crystal", ascending=False)
    plast.to_csv(REPO / "data" / "processed" / "side_chain_q3_plasticity.csv", index=False)

    q3 = {
        "n_residues": int(len(plast)),
        "n_residues_dropped_no_chi1": int(n_before - len(plast)),
        "pooled_sd_chi1_across_ligands_model_deg": float(plast.sd_model.median()),
        "pooled_sd_chi1_across_ligands_crystal_deg": float(plast.sd_crystal.median()),
        "pooled_sd_chi1_within_ligand_model_deg": float(plast.sd_within_ligand.median()),
        "ratio_model_over_crystal_median": float(plast.ratio.median()),
        "residues_model_more_variable": int((plast.ratio > 1).sum()),
        "residues_crystal_more_variable": int((plast.ratio < 1).sum()),
        "wilcoxon_model_vs_crystal_p": float(
            stats.wilcoxon(plast.sd_model, plast.sd_crystal).pvalue),
    }

    # correlation of pocket side-chain error with the pose metrics
    ship = pd.read_csv(REPO / "data" / "processed" / "poses_scored_val87b.csv")
    ship = ship[ship.arm == "unsteered"][["ligand", "sample", "lddt_pli", "bisy_rmsd"]]
    feat = (chi[chi.contact].groupby(["ligand", "sample"])
                .agg(mean_adchi1=("dchi1", lambda s: float(np.nanmean(np.abs(s)))),
                     mean_sc_local=("sc_rmsd_local", "mean"),
                     mean_sc_pocket=("sc_rmsd_coreframe", "mean"),
                     mean_ca=("ca_err", "mean")).reset_index())
    mm = feat.merge(ship, on=["ligand", "sample"], how="inner")
    q3["n_poses_scored"] = int(len(mm))
    q3["exact_zero_lddt_rows"] = int((mm.lddt_pli == 0).sum())

    def corr_block(col):
        bt = {k: float(stats.spearmanr(mm.groupby("ligand")[col].mean(),
                                       mm.groupby("ligand")[k].mean()).statistic)
              for k in ("lddt_pli", "bisy_rmsd")}
        wl = {}
        for k in ("lddt_pli", "bisy_rmsd"):
            rs = []
            for _l, gg in mm.groupby("ligand"):
                if gg[col].nunique() < 3 or gg[k].nunique() < 3:
                    continue
                rs.append(float(stats.spearmanr(gg[col], gg[k]).statistic))
            rs = np.array([r for r in rs if np.isfinite(r)])
            exp_neg = PREREG["expected_sign"][f"rho(sidechain_error, {k})"] == "negative"
            wl[k] = {"n_ligands": int(len(rs)), "mean_rho": float(rs.mean()),
                     "frac_expected_sign": float((rs < 0).mean() if exp_neg
                                                 else (rs > 0).mean()),
                     "wilcoxon_p": float(stats.wilcoxon(rs).pvalue) if len(rs) > 5 else 1.0}
        return {"between_ligand": bt, "within_ligand": wl}

    q3["correlations"] = {c: corr_block(c) for c in
                          ("mean_adchi1", "mean_sc_local", "mean_sc_pocket", "mean_ca")}
    ctl = REPO / "data" / "processed" / "side_chain_plasticity_control.json"
    if ctl.exists():
        q3["cartesian_control"] = json.loads(ctl.read_text())

    summary = {"prereg": PREREG, "controls": controls,
               "Q0_backbone": bb,
               "Q1_which_residues_are_wrong": {
                   "pooled": pooled,
                   "top_by_frac_chi1_off": [
                       {"res": f"{r.resname}{r.resnum}", "n": int(r.n_pairs),
                        "med_dchi1": round(float(r.med_dchi1), 1),
                        "iqr": [round(float(r.q1_dchi1), 1), round(float(r.q3_dchi1), 1)],
                        "frac_off": round(float(r.frac_chi1_off), 3),
                        "med_sc_rmsd_local": round(float(r.med_sc_local), 2),
                        "med_ca_err": round(float(r.med_ca_err), 2),
                        "fg": bool(r.in_fg)}
                       for _, r in q1.head(18).iterrows()],
                   "bottom_by_frac_chi1_off": [
                       {"res": f"{r.resname}{r.resnum}", "n": int(r.n_pairs),
                        "frac_off": round(float(r.frac_chi1_off), 3),
                        "med_sc_rmsd_local": round(float(r.med_sc_local), 2)}
                       for _, r in q1.tail(8).iterrows()]},
               "Q2_which_residues_block": q2h,
               "Q2_top_blockers": [
                   {"res": f"{r.worst_resname}{r.worst_resnum}",
                    "n_poses_worst": int(r.n_poses), "n_ligands": int(r.n_ligands),
                    "frac_of_clashing_poses": round(float(r.frac_of_clashing_poses), 3),
                    "med_dist": round(float(r.med_dist), 2)}
                   for _, r in wa.head(15).iterrows()],
               "Q3_rigid_or_wrongly_adaptive": q3}
    out = REPO / "data" / "processed" / "side_chain_diagnosis.json"
    out.write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


# --------------------------------------------------------------------------
# stage 2b — the Cartesian control on Q3
# --------------------------------------------------------------------------

def cmd_plasticity() -> None:
    """Does the model's pocket move from ligand to ligand AT ALL, in Cartesian space?

    The chi-angle answer can be doubted on symmetry handling, so this repeats it with no
    torsions in it: superpose residue r of structure i onto residue r of structure j by
    that residue's OWN N/CA/C, then take the side-chain heavy-atom RMSD. Frame-free, and
    the same operation for the 87 predictions and for the 87 crystals.
    """
    import pandas as pd
    from cypstruct import pose as P

    q1 = pd.read_csv(REPO / "data" / "processed" / "side_chain_q1_per_residue.csv")
    core_pocket = sorted(q1.resnum.tolist())
    lig = pd.read_csv(REPO / "data" / "processed" / "validation_ligands.csv")
    have = {p.name.split("__")[0] for p in POOL.iterdir() if p.is_dir()}

    def collect(cx, rnames):
        out = {}
        for r in core_pocket:
            at = residue_atoms(cx, r)
            rn = next((v for (_c, k), v in cx.prot_res.items() if k == r), None)
            if rn is None or not all(k in at for k in ("N", "CA", "C")):
                continue
            canonicalise_prochiral(at, rn)
            names = [a for a in SC_ATOMS.get(rn, []) if a != "CB" and a in at]
            if not names:
                continue
            out[r] = (rn, at)
            rnames.setdefault(r, set()).add(rn)
        return out

    models, crystals, rnames = [], [], {}
    for r in lig.itertuples():
        if r.id not in have:
            continue
        ref, _c, _ch = load_reference(r.pdb, r.id)
        if ref is None:
            continue
        m = P.load_structure(POOL / f"{r.id}__unsteered__s1" / "input_model_0.cif")
        m, _o, _i = renumber_to_reference(m, ref)
        models.append(collect(m, rnames))
        crystals.append(collect(ref, rnames))

    def spread(sets):
        acc = {}
        for r in core_pocket:
            present = [s[r] for s in sets if r in s]
            vals = []
            for i in range(len(present)):
                for j in range(i + 1, len(present)):
                    (na, a), (nb, b) = present[i], present[j]
                    if na != nb:
                        continue
                    A = np.array([a[k] for k in ("N", "CA", "C")])
                    B = np.array([b[k] for k in ("N", "CA", "C")])
                    R, t, _ = P.kabsch(A, B)
                    vals.append(sc_rmsd({k: v @ R.T + t for k, v in a.items()}, b, na))
            vals = [v for v in vals if np.isfinite(v)]
            if len(vals) >= 50:
                acc[r] = (float(np.median(vals)), len(present))
        return acc

    sm, sc = spread(models), spread(crystals)
    common = sorted(set(sm) & set(sc))
    rows = [{"resnum": r, "resname": sorted(rnames[r])[0],
             "n_model": sm[r][1], "n_crystal": sc[r][1],
             "median_pairwise_sc_rmsd_model": round(sm[r][0], 3),
             "median_pairwise_sc_rmsd_crystal": round(sc[r][0], 3)}
            for r in common]
    rows.sort(key=lambda x: -x["median_pairwise_sc_rmsd_crystal"])
    mv = np.array([r["median_pairwise_sc_rmsd_model"] for r in rows])
    cv = np.array([r["median_pairwise_sc_rmsd_crystal"] for r in rows])
    from scipy import stats
    res = {
        "what": "median pairwise side-chain heavy-atom RMSD (A) between two different "
                "ligands' structures, after superposing that residue's own N/CA/C",
        "n_residues": len(rows), "n_structures": len(models),
        "model_median_A": float(np.median(mv)), "crystal_median_A": float(np.median(cv)),
        "model_p90_A": float(np.percentile(mv, 90)),
        "crystal_p90_A": float(np.percentile(cv, 90)),
        "residues_model_more_variable": int((mv > cv).sum()),
        "residues_crystal_more_variable": int((cv > mv).sum()),
        "wilcoxon_p": float(stats.wilcoxon(mv, cv).pvalue),
        "per_residue": rows,
    }
    (REPO / "data" / "processed" / "side_chain_plasticity_control.json").write_text(
        json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "per_residue"}, indent=1))
    print(json.dumps(rows[:12], indent=1))


# --------------------------------------------------------------------------
# stage 2c — the template ceiling over the LESION, 210-216 (addendum)
# --------------------------------------------------------------------------

LESION = list(range(210, 217))     # the span Q4 named; fixed before this stage was run


def _entry_chain(path, want_chains):
    """One chain per deposited entry: the listed auth chain with the most residues.

    Two chains of the same crystal are not two opinions, so an entry contributes once.
    """
    import gemmi
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    best, bestn = None, -1
    for ch in st[0]:
        if want_chains and ch.name not in want_chains:
            continue
        n = sum(1 for r in ch
                if (gemmi.find_tabulated_residue(r.name.strip().upper()) or None)
                and gemmi.find_tabulated_residue(r.name.strip().upper()).is_amino_acid())
        if n > bestn:
            best, bestn = ch, n
    if best is None or bestn < 200:
        return None
    ca, occ, bf = {}, {}, {}
    for r in best:
        info = gemmi.find_tabulated_residue(r.name.strip().upper())
        if not (info and info.is_amino_acid()):
            continue
        for at in r:
            if at.name.strip() == "CA":
                ca[r.seqid.num] = np.array([at.pos.x, at.pos.y, at.pos.z])
                occ[r.seqid.num] = float(at.occ)
                bf[r.seqid.num] = float(at.b_iso)
    bbcb = {}
    for r in best:
        info = gemmi.find_tabulated_residue(r.name.strip().upper())
        if not (info and info.is_amino_acid()) or r.seqid.num not in LESION:
            continue
        for at in r:
            nm = at.name.strip()
            if nm in ("N", "CA", "C", "O", "CB"):
                bbcb[(r.seqid.num, nm)] = np.array([at.pos.x, at.pos.y, at.pos.z])
    res = float("nan")
    try:
        res = float(st.resolution)
    except Exception:
        pass
    return dict(ca=ca, occ=occ, b=bf, bbcb=bbcb, chain=best.name,
                resolution=res, spacegroup=str(st.spacegroup_hm))


def cmd_templates() -> None:
    """The measurement FINDING 028 named as deciding whether holo templates come back.

    A blind template can only ever be as good as one deposited crystal is at predicting
    another. So: over residues 210-216 ALONE, superpose two deposited CYP3A4 entries on
    the ligand-free rigid core and measure the CA deviation in the span. If that spread
    is not clearly below the model's own error there, a template cannot help in the
    lesion and FINDING 024's dismissal stands unqualified.

    Decision rule, fixed before any number was read (see the addendum's preamble):
    LICENSED only if the holo-holo DIFFERENT-ligand spread is clearly below the model's
    median AND p90 error over the same span.
    """
    import pandas as pd
    from cypstruct import pose as P
    from cypstruct.targets import IGNORE_HET

    meta = json.loads((REPO / "data" / "processed" / "p450_universe"
                       / "entry_meta.json").read_text())
    rcsb = REPO / "data" / "reference" / "rcsb"

    entries = {}
    for pid, v in meta.items():
        accs = {a.get("reference_database_accession")
                for pe in (v.get("polymer_entities") or [])
                for a in (pe.get("rcsb_polymer_entity_align") or [])}
        if "P08684" not in accs:
            continue
        cif = rcsb / f"{pid}.cif"
        if not cif.exists():
            continue
        chains = [c for pe in (v.get("polymer_entities") or [])
                  for c in (pe.get("rcsb_polymer_entity_container_identifiers", {})
                            .get("auth_asym_ids") or [])]
        ligs = set()
        for ne in (v.get("nonpolymer_entities") or []):
            cid = (ne.get("nonpolymer_comp", {}).get("chem_comp", {}) or {}).get("id")
            if cid and cid.upper() not in IGNORE_HET and cid.upper() != "HEM":
                ligs.add(cid.upper())
        e = _entry_chain(cif, set(chains))
        if e is None:
            continue
        e["ligands"] = sorted(ligs)
        e["apo"] = not ligs
        entries[pid] = e
    ids = sorted(entries)

    # ---- FIRST: is the span even modelled? -------------------------------
    cov = {p: sum(1 for r in LESION if r in entries[p]["ca"]) for p in ids}
    full = [p for p in ids if cov[p] == len(LESION)]
    part = [p for p in ids if 4 <= cov[p] < len(LESION)]
    none_ = [p for p in ids if cov[p] == 0]
    modelled = {
        "n_entries": len(ids),
        "span": f"{LESION[0]}-{LESION[-1]}",
        "fully_modelled": len(full), "frac_fully_modelled": len(full) / len(ids),
        "partially_modelled_4_to_6": len(part),
        "completely_absent": len(none_),
        "per_residue_modelled_frac": {
            str(r): round(float(np.mean([r in entries[p]["ca"] for p in ids])), 3)
            for r in LESION},
        "apo_entries": int(sum(entries[p]["apo"] for p in ids)),
        "holo_entries": int(sum(not entries[p]["apo"] for p in ids)),
        "apo_fully_modelled": int(sum(entries[p]["apo"] for p in full)),
        "median_CA_B_in_span": float(np.median(
            [entries[p]["b"][r] for p in full for r in LESION])),
        "min_CA_occupancy_in_span": float(min(
            entries[p]["occ"][r] for p in full for r in LESION)),
        "distinct_space_groups": len({entries[p]["spacegroup"] for p in ids}),
        "space_group_counts": {k: int(v) for k, v in
                               pd.Series([entries[p]["spacegroup"] for p in ids])
                               .value_counts().items()},
        "resolution_median": float(np.nanmedian([entries[p]["resolution"] for p in ids])),
        "resolution_range": [float(np.nanmin([entries[p]["resolution"] for p in ids])),
                             float(np.nanmax([entries[p]["resolution"] for p in ids]))],
        "distinct_ligand_codes": len({tuple(entries[p]["ligands"]) for p in ids}),
    }

    # ---- the crystal-to-crystal spread ------------------------------------
    rows = []
    for i in range(len(full)):
        for j in range(i + 1, len(full)):
            a, b = entries[full[i]], entries[full[j]]
            sh = sorted(set(a["ca"]) & set(b["ca"]) & CORE)
            if len(sh) < 50:
                continue
            R, t, fit = P.kabsch(np.array([a["ca"][k] for k in sh]),
                                 np.array([b["ca"][k] for k in sh]))
            d = np.array([np.linalg.norm(a["ca"][r] @ R.T + t - b["ca"][r])
                          for r in LESION])
            if a["apo"] and b["apo"]:
                s = "apo-apo"
            elif a["apo"] or b["apo"]:
                s = "apo-holo"
            elif set(a["ligands"]) & set(b["ligands"]):
                s = "holo-holo same ligand"
            else:
                s = "holo-holo different ligand"
            rows.append(dict(a=full[i], b=full[j], stratum=s, core_fit=fit,
                             rms=float(np.sqrt((d ** 2).mean())), max=float(d.max()),
                             same_sg=a["spacegroup"] == b["spacegroup"],
                             same_res=abs((a["resolution"] or 0) -
                                          (b["resolution"] or 0)) < 0.05))
    cc = pd.DataFrame(rows)
    cc.to_csv(REPO / "data" / "processed" / "template_ceiling_pairs.csv", index=False)

    def block(frame):
        if not len(frame):
            return None
        q = frame["rms"]
        return {"n_pairs": int(len(frame)), "median_A": float(q.median()),
                "p25_A": float(q.quantile(.25)), "p75_A": float(q.quantile(.75)),
                "p90_A": float(q.quantile(.90)), "max_A": float(q.max()),
                "frac_over_1A": float((q > 1.0).mean()),
                "frac_over_2A": float((q > 2.0).mean()),
                "median_core_fit_A": float(frame.core_fit.median())}

    spread = {s: block(cc[cc.stratum == s]) for s in sorted(cc.stratum.unique())}
    spread["ALL pairs"] = block(cc)
    hh = cc[cc.stratum == "holo-holo different ligand"]
    spread["holo-holo different ligand, CROSS space group only"] = block(
        hh[~hh.same_sg])

    # Is the spread one broad distribution or several discrete loop states? Descriptive,
    # not a decision: an entry's median deviation to every OTHER entry, then split at
    # 1 A. The cut is descriptive and the whole distribution is in the CSV.
    med = {e: float(pd.concat([cc[cc.a == e]["rms"], cc[cc.b == e]["rms"]]).median())
           for e in full}
    major = sorted([e for e, v in med.items() if v < 1.0])
    minor = sorted([e for e, v in med.items() if v >= 3.0])
    within = cc[cc.a.isin(major) & cc.b.isin(major)]
    states = {
        "entries_with_span": len(full),
        "majority_cluster_n": len(major),
        "within_majority_median_A": float(within["rms"].median()),
        "within_majority_p90_A": float(within["rms"].quantile(.90)),
        "minority_conformers_n": len(minor), "minority_conformers": minor,
        "intermediate_n": len(full) - len(major) - len(minor),
        "note": "a blind template cannot know which state the query is in, and "
                "selecting inside the majority cluster is not available blind",
    }

    # ---- the model's own error over the SAME span, same frame -------------
    lig = pd.read_csv(REPO / "data" / "processed" / "validation_ligands.csv")
    have = {p.name.split("__")[0] for p in POOL.iterdir() if p.is_dir()}
    mrows, nspan_full = [], []
    for r in lig.itertuples():
        if r.id not in have:
            continue
        ref, _c, _ch = load_reference(r.pdb, r.id)
        if ref is None:
            continue
        rn = {k: v for (_c2, k), v in ref.ca().items()}
        span = [x for x in LESION if x in rn]
        nspan_full.append(len(span) == len(LESION))
        if len(span) < 4:
            mrows.append(dict(ligand=r.id, pdb=r.pdb, n_span=len(span), rms=np.nan))
            continue
        for cif in sorted((POOL / f"{r.id}__unsteered__s1").glob("input_model_*.cif")):
            m = P.load_structure(cif)
            m, _o, _i = renumber_to_reference(m, ref)
            mn = {k: v for (_c2, k), v in m.ca().items()}
            shc = sorted(set(mn) & set(rn) & CORE)
            R, t, _f = P.kabsch(np.array([mn[k] for k in shc]),
                                np.array([rn[k] for k in shc]))
            d = np.array([np.linalg.norm(mn[x] @ R.T + t - rn[x])
                          for x in span if x in mn])
            mrows.append(dict(ligand=r.id, pdb=r.pdb, sample=cif.stem,
                              n_span=len(d), rms=float(np.sqrt((d ** 2).mean()))))
    md = pd.DataFrame(mrows)

    mq = md.rms.dropna()
    model = {"n_poses": int(len(mq)),
             "n_validation_crystals": int(len(nspan_full)),
             "n_validation_crystals_span_FULLY_modelled": int(sum(nspan_full)),
             "n_pairs_with_span_modelled": int(md.dropna(subset=["rms"]).ligand.nunique()),
             "median_A": float(mq.median()), "p90_A": float(mq.quantile(.90)),
             "frac_over_1A": float((mq > 1.0).mean()),
             "frac_over_2A": float((mq > 2.0).mean())}

    # ---- the seven unrescuable ligands ------------------------------------
    hard7 = ["1RD", "5AW", "A1A4T", "ERY", "MWY", "QEP", "X7P"]
    seven = []
    for lg in hard7:
        r = lig[lig.id == lg]
        if not len(r):
            continue
        r = r.iloc[0]
        ref, _c, _ch = load_reference(r.pdb, lg)
        rn = {k: v for (_c2, k), v in ref.ca().items()}
        L = np.asarray(ref.lig_xyz, float)
        best, clears = None, 0
        tried = 0
        for p in full:
            if p.upper() == str(r.pdb).upper():
                continue
            e = entries[p]
            sh = sorted(set(e["ca"]) & set(rn) & CORE)
            if len(sh) < 50:
                continue
            R, t, _f = P.kabsch(np.array([e["ca"][k] for k in sh]),
                                np.array([rn[k] for k in sh]))
            pts = np.array([v @ R.T + t for v in e["bbcb"].values()])
            if not len(pts):
                continue
            tried += 1
            dmin = float(np.linalg.norm(pts[:, None, :] - L[None, :, :],
                                        axis=2).min())
            clears += int(dmin >= PREREG["clash_cut_A"])
            if best is None or dmin > best[1]:
                best = (p, dmin)
        # the control: the query's OWN crystal, same atoms
        own = np.array([v for k, v in
                        ((kk, vv) for kk, vv in _entry_chain(
                            REPO / "data" / "reference" / "rcsb" / f"{r.pdb}.cif",
                            None)["bbcb"].items())])
        own_min = (float(np.linalg.norm(own[:, None, :] - L[None, :, :], axis=2).min())
                   if len(own) else None)
        seven.append(dict(ligand=lg, pdb=str(r.pdb),
                          own_span_fully_modelled=bool(str(r.pdb).upper() in med),
                          own_median_dev_to_other_entries_A=(
                              round(med[str(r.pdb).upper()], 2)
                              if str(r.pdb).upper() in med else None),
                          donors_tested=tried,
                          donors_clearing=clears,
                          frac_clearing=(clears / tried) if tried else None,
                          best_donor=(best[0] if best else None),
                          best_donor_min_contact_A=(round(best[1], 2) if best else None),
                          own_crystal_min_contact_A=(round(own_min, 2) if own_min
                                                     else None)))

    out = {
        "what": "template ceiling over the lesion: how well one deposited CYP3A4 "
                "crystal predicts another's CA over residues 210-216, after "
                "superposition on CYP3A4_RIGID_CORE",
        "decision_rule_fixed_in_advance":
            "LICENSED only if the holo-holo DIFFERENT-ligand spread is clearly below "
            "the model's median AND p90 error over the same span; comparable or worse "
            "means templates are refuted in the lesion too",
        "modelled_first": modelled,
        "crystal_to_crystal_spread": spread,
        "conformational_states": states,
        "model_error_same_span_same_frame": model,
        "seven_unrescuable": seven,
    }
    (REPO / "data" / "processed" / "template_ceiling_lesion.json").write_text(
        json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


# --------------------------------------------------------------------------
# stage 3 — Q4, the repack scan
# --------------------------------------------------------------------------

def repack_ligand(args) -> dict:
    lig_id, pdb, smiles = args
    from scipy.spatial import cKDTree
    from cypstruct import pose as P

    cut = PREREG["clash_cut_A"]
    selfcut = PREREG["self_clash_cut_A"]
    g1 = np.arange(-180.0, 180.0, PREREG["chi_grid_deg"])
    g3 = np.arange(-180.0, 180.0, PREREG["chi3_grid_deg"])

    ref, _cif, _ch = load_reference(pdb, lig_id)
    if ref is None:
        return {"ligand": lig_id, "status": "no-reference"}
    alignset = P.pocket_residues_from_structure(ref, radius=PREREG["align_radius_A"])
    job = POOL / f"{lig_id}__unsteered__s1"
    cifs = sorted(job.glob("input_model_*.cif"),
                  key=lambda p: int(p.stem.rsplit("_", 1)[1]))

    out_poses = []
    for cif in cifs:
        m = P.load_structure(cif)
        m, _o, _i = renumber_to_reference(m, ref)
        perm = P.best_ligand_mapping(smiles, m, ref)
        if perm is None:
            continue
        mn = {k: v for (_c, k), v in m.ca().items()}
        rn = {k: v for (_c, k), v in ref.ca().items()}
        sh = sorted(set(mn) & set(rn) & set(alignset))
        Rk, tk, _ = P.kabsch(np.array([rn[k] for k in sh]),
                             np.array([mn[k] for k in sh]))
        L = ref.lig_xyz[perm] @ Rk.T + tk

        hm = np.asarray(m.heme_xyz, float)
        he = [e.upper() for e in m.heme_elem]
        hn = hm[[i for i, e in enumerate(he) if e != "FE"]] if len(hm) else hm

        # index the model's protein by residue
        byres: dict[int, dict[str, int]] = {}
        for i, (_c, r, a) in enumerate(m.prot_key):
            byres.setdefault(r, {})[a] = i
        resname = {r: v for (_c, r), v in m.prot_res.items()}

        # residues whose atoms sit under the cutoff against the crystal ligand
        dall = np.linalg.norm(m.prot_xyz[:, None, :] - L[None, :, :], axis=2).min(1)
        d_heme = (np.linalg.norm(hn[:, None, :] - L[None, :, :], axis=2).min()
                  if len(hn) else np.inf)
        clashing = sorted({m.prot_key[i][1] for i in np.where(dall < cut)[0]})
        bb_clash = sorted({m.prot_key[i][1] for i in np.where(dall < cut)[0]
                           if m.prot_key[i][2] in BACKBONE})
        movable = [r for r in clashing
                   if resname.get(r) in ROT and any(
                       a in byres[r] for a in ("CA", "CB"))]

        # environment that does NOT move: everything except the movable residues
        mov_idx = set()
        for r in movable:
            for a, i in byres[r].items():
                if a not in BACKBONE:
                    mov_idx.add(i)
        fixed = np.array([i for i in range(len(m.prot_xyz)) if i not in mov_idx])
        env_fixed = np.vstack([m.prot_xyz[fixed], hn]) if len(hn) else m.prot_xyz[fixed]
        fixed_res = np.array([m.prot_key[i][1] for i in fixed] + [-1] * len(hn))
        t_lig = cKDTree(L)

        per_res = []
        new_pos = {}
        for r in movable:
            rname = resname[r]
            at = {a: m.prot_xyz[i].copy() for a, i in byres[r].items()}
            axes = [ax for ax in ROT[rname]
                    if ax[0] in at and ax[1] in at and any(x in at for x in ax[2])]
            if not axes:
                per_res.append(dict(resnum=int(r), resname=rname, cleared=False,
                                    reason="no-axis", pert=None))
                continue
            # neighbours excluded from the self-clash test: own residue and i +/- 1
            keep = ~np.isin(fixed_res, [r - 1, r, r + 1])
            tree_env = cKDTree(env_fixed[keep])

            names = [a for a in SC_ATOMS[rname] if a in at]
            base = np.array([at[a] for a in names])

            def scan(grids):
                """Cartesian grid over the listed torsions. Returns the minimum-magnitude
                solution that clears the ligand AND stays self-consistent, plus the
                minimum-magnitude solution that clears the ligand while IGNORING the
                self-clash test — the pair separates "no room in the pocket" from "the
                rest of the protein is in the way"."""
                best = (np.inf, base, [])
                free = np.inf
                from itertools import product
                for deltas in product(*grids):
                    pts = base.copy()
                    for k, dg in enumerate(deltas):
                        if dg == 0:
                            continue
                        a0, a1, distal = axes[k]
                        sel = [j for j, nm in enumerate(names) if nm in distal]
                        if not sel:
                            continue
                        pts[sel] = rot_about(pts[sel], at[a0], at[a1],
                                             np.array([dg]))[0]
                    if t_lig.query(pts, k=1)[0].min() < cut:
                        continue
                    mag = float(max(abs(d) for d in deltas)) if deltas else 0.0
                    free = min(free, mag)
                    if tree_env.query(pts, k=1)[0].min() < selfcut:
                        continue
                    if mag < best[0]:
                        best = (mag, pts.copy(), list(deltas))
                return best, free

            grids = [g1 for _ in axes[:2]]
            res, free = scan(grids)
            if not np.isfinite(res[0]) and len(axes) >= 3:
                res, free2 = scan([g1, g1, g3])   # chi3 only when chi1xchi2 fails
                free = min(free, free2)
            cbd = (float(np.linalg.norm(at["CB"] - L, axis=1).min())
                   if "CB" in at else float("nan"))
            if np.isfinite(res[0]):
                new_pos[r] = (names, res[1])
                per_res.append(dict(resnum=int(r), resname=rname, cleared=True,
                                    pert=float(res[0]), cb_dist=cbd,
                                    deltas=[float(d) for d in res[2]]))
            else:
                per_res.append(dict(
                    resnum=int(r), resname=rname, cleared=False, pert=None, cb_dist=cbd,
                    reason=("CB-itself-clashes" if cbd < cut else
                            "self-clash-blocks" if np.isfinite(free) else
                            "no-clearing-rotamer"),
                    free_pert=(None if not np.isfinite(free) else float(free))))

        # concerted check: apply every clearing rotamer at once
        moved = m.prot_xyz.copy()
        for r, (names, pts) in new_pos.items():
            for a, p in zip(names, pts):
                moved[byres[r][a]] = p
        env2 = np.vstack([moved, hn]) if len(hn) else moved
        final_min = float(cKDTree(env2).query(L, k=1)[0].min())
        # do the repacked side chains crash into each other?
        inter = np.inf
        keys = list(new_pos)
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                a, b = new_pos[keys[i]][1], new_pos[keys[j]][1]
                if abs(keys[i] - keys[j]) <= 1:
                    continue
                inter = min(inter, float(np.linalg.norm(a[:, None] - b[None], axis=2).min()))

        out_poses.append(dict(
            sample=cif.stem, n_clashing=len(clashing), n_movable=len(movable),
            n_backbone_clashing=len(bb_clash), backbone_res=bb_clash,
            heme_min=float(d_heme), start_min=float(min(dall.min(), d_heme)),
            per_res=per_res,
            n_cleared=int(sum(p["cleared"] for p in per_res)),
            n_moved=int(sum(p["cleared"] and p["pert"] > 0 for p in per_res)),
            n_uncleared=int(sum(not p["cleared"] for p in per_res)),
            max_pert=float(max([p["pert"] for p in per_res if p["cleared"]], default=0.0)),
            sum_pert=float(sum(p["pert"] for p in per_res if p["cleared"])),
            final_min=final_min,
            solved=bool(final_min >= cut),
            reasons={k: int(sum(1 for p in per_res if p.get("reason") == k))
                     for k in ("CB-itself-clashes", "self-clash-blocks",
                               "no-clearing-rotamer", "no-axis")},
            inter_sidechain_min=float(inter) if np.isfinite(inter) else None))
    rec = {"ligand": lig_id, "pdb": pdb, "poses": out_poses}
    (SCRATCH / f"repack_{lig_id}.json").write_text(json.dumps(rec))
    return {"ligand": lig_id, "status": "ok", "poses": len(out_poses),
            "solved": int(sum(p["solved"] for p in out_poses))}


def cmd_repack(workers: int) -> None:
    import pandas as pd
    from scipy.spatial import cKDTree
    from cypstruct import pose as P
    from cypstruct.targets import CYP3A4_FG_LOOP

    SCRATCH.mkdir(parents=True, exist_ok=True)
    exp = pd.read_csv(REPO / "data" / "processed" / "orientation_expansion_per_ligand.csv")
    lig = pd.read_csv(REPO / "data" / "processed" / "validation_ligands.csv")
    hard = exp[(exp.orig_min_rmsd >= 2) & (exp.grid_min_rmsd >= 2)].ligand.tolist()
    fail = exp[exp.orig_min_rmsd >= 2].ligand.tolist()
    soft = [x for x in fail if x not in set(hard)]
    solved_pairs = exp[exp.orig_min_rmsd < 2].ligand.tolist()
    # the whole 87 are scanned, so the 46 pairs the model ALREADY solves act as the
    # control: on those the crystal ligand should need little or no repacking at all.
    jobs = [(r.id, r.pdb, r.smiles) for r in lig.itertuples()
            if r.id in set(exp.ligand)]
    todo = [j for j in jobs if not (SCRATCH / f"repack_{j[0]}.json").exists()]
    print(f"rotation-unreachable failures: {len(hard)} of {len(fail)} failures; "
          f"{len(todo)} of {len(jobs)} to do", flush=True)
    if todo:
        if workers <= 1:
            for j in todo:
                print(repack_ligand(j), flush=True)
        else:
            import multiprocessing as mp
            with mp.Pool(workers) as pool:
                for r in pool.imap_unordered(repack_ligand, todo):
                    print(r, flush=True)

    allrecs = [json.loads((SCRATCH / f"repack_{l}.json").read_text())
               for l in exp.ligand if (SCRATCH / f"repack_{l}.json").exists()]
    strata = {}
    for name, keys in (("rotation_unreachable_failures", hard),
                       ("rotation_rescuable_failures", soft),
                       ("already_solved_control", solved_pairs)):
        rr = [dict(ligand=r["ligand"], **{k: v for k, v in p.items()
                                          if k not in ("per_res", "backbone_res",
                                                       "reasons")})
              for r in allrecs if r["ligand"] in set(keys) for p in r["poses"]]
        if not rr:
            continue
        dd = pd.DataFrame(rr)
        strata[name] = {
            "n_ligands": int(dd.ligand.nunique()), "n_poses": int(len(dd)),
            "median_start_contact_A": float(dd.start_min.median()),
            "frac_poses_clashing_before": float((dd.start_min < PREREG["clash_cut_A"]).mean()),
            "median_clashing_residues": float(dd.n_clashing.median()),
            "frac_poses_with_backbone_clash": float((dd.n_backbone_clashing > 0).mean()),
            "frac_poses_solved_by_chi_scan": float(dd.solved.mean()),
            "median_max_chi_perturbation_deg_solved": (
                float(dd[dd.solved].max_pert.median()) if dd.solved.any() else None),
            "median_residues_moved_solved": (
                float(dd[dd.solved].n_moved.median()) if dd.solved.any() else None)}

    recs = [r for r in allrecs if r["ligand"] in set(hard)]
    rows, reasons = [], {}
    immov: dict[str, dict] = {}
    for rec in recs:
        for p in rec["poses"]:
            for k, v in p["reasons"].items():
                reasons[k] = reasons.get(k, 0) + v
            for pr in p["per_res"]:
                if pr.get("reason") == "CB-itself-clashes":
                    key = f"{pr['resname']}{pr['resnum']}"
                    e = immov.setdefault(key, {"poses": 0, "ligands": set(),
                                               "cb_dist": []})
                    e["poses"] += 1
                    e["ligands"].add(rec["ligand"])
                    e["cb_dist"].append(pr["cb_dist"])
            rows.append(dict(ligand=rec["ligand"],
                             **{k: v for k, v in p.items()
                                if k not in ("per_res", "backbone_res", "reasons")},
                             n_backbone=p["n_backbone_clashing"]))
    d = pd.DataFrame(rows)
    d.to_csv(REPO / "data" / "processed" / "side_chain_q4_repack.csv", index=False)

    perlig = d.groupby("ligand").agg(
        n_poses=("solved", "size"), n_solved=("solved", "sum"),
        best_max_pert=("max_pert", "min"), med_max_pert=("max_pert", "median"),
        med_moved=("n_moved", "median"), med_clashing=("n_clashing", "median"),
        med_bb=("n_backbone", "median"), med_start=("start_min", "median"),
        med_final=("final_min", "median")).reset_index()
    perlig["any_solved"] = perlig.n_solved > 0
    solved = d[d.solved]

    res = {
        "prereg": {k: PREREG[k] for k in ("clash_cut_A", "self_clash_cut_A",
                                          "chi_grid_deg", "chi3_grid_deg")},
        "scope": {"n_ligands": int(len(perlig)), "n_poses": int(len(d)),
                  "definition": "failures (no sub-2A pose in 20 samples) that the whole "
                                "1024-orientation grid also cannot rescue"},
        "strata": strata,
        "clashing_residues_per_pose": {
            "median": float(d.n_clashing.median()),
            "iqr": [float(d.n_clashing.quantile(.25)), float(d.n_clashing.quantile(.75))],
            "max": int(d.n_clashing.max()),
            "median_movable": float(d.n_movable.median()),
            "median_backbone_involved": float(d.n_backbone.median()),
            "frac_poses_with_backbone_clash": float((d.n_backbone > 0).mean())},
        "outcome": {
            "poses_solved_by_chi_scan": int(d.solved.sum()), "poses": int(len(d)),
            "frac_poses_solved": float(d.solved.mean()),
            "ligands_with_any_solved_pose": int(perlig.any_solved.sum()),
            "ligands": int(len(perlig)),
            "why_a_residue_could_not_be_cleared": reasons},
        "magnitude_on_solved_poses": {
            "n": int(len(solved)),
            "median_residues_moved": float(solved.n_moved.median()),
            "residues_moved_iqr": [float(solved.n_moved.quantile(.25)),
                                   float(solved.n_moved.quantile(.75))],
            "median_max_chi_perturbation_deg": float(solved.max_pert.median()),
            "max_chi_perturbation_iqr": [float(solved.max_pert.quantile(.25)),
                                         float(solved.max_pert.quantile(.75))],
            "median_sum_chi_perturbation_deg": float(solved.sum_pert.median()),
            "frac_single_residue": float((solved.n_moved <= 1).mean()),
            "frac_two_or_fewer": float((solved.n_moved <= 2).mean()),
            "frac_zero_residues_moved": float((solved.n_moved == 0).mean()),
            "frac_max_pert_le_30deg": float((solved.max_pert <= 30).mean()),
            "frac_max_pert_le_60deg": float((solved.max_pert <= 60).mean())},
        "per_ligand_best_case": [
            {"ligand": r.ligand, "solved_poses": int(r.n_solved),
             "best_max_pert_deg": (None if not np.isfinite(r.best_max_pert)
                                   else round(float(r.best_max_pert), 1)),
             "median_clashing_residues": float(r.med_clashing),
             "median_backbone_clashes": float(r.med_bb),
             "median_start_contact": round(float(r.med_start), 2),
             "median_final_contact": round(float(r.med_final), 2)}
            for _, r in perlig.sort_values("n_solved").iterrows()],
    }

    res["immovable_blockers"] = sorted(
        [{"res": k, "poses": v["poses"], "n_ligands": len(v["ligands"]),
          "median_cb_dist_A": round(float(np.median(v["cb_dist"])), 2),
          "in_fg_span": bool(FG[0] <= int(k[3:]) <= FG[1])}
         for k, v in immov.items()], key=lambda x: -x["poses"])

    # what is wrong where the scan fails
    unsolved_lig = perlig[~perlig.any_solved].ligand.tolist()
    res["unrepackable"] = {"n_ligands": len(unsolved_lig), "ligands": unsolved_lig}
    extra = []
    for lg in unsolved_lig:
        r = lig[lig.id == lg].iloc[0]
        ref, _c, _ch = load_reference(r.pdb, lg)
        alignset = P.pocket_residues_from_structure(ref, radius=PREREG["align_radius_A"])
        cifs = sorted((POOL / f"{lg}__unsteered__s1").glob("input_model_*.cif"))
        m = P.load_structure(cifs[0])
        m, _o, _i = renumber_to_reference(m, ref)
        perm = P.best_ligand_mapping(r.smiles, m, ref)
        mn = {k: v for (_c, k), v in m.ca().items()}
        rn = {k: v for (_c, k), v in ref.ca().items()}
        sh = sorted(set(mn) & set(rn) & set(alignset))
        Rk, tk, _ = P.kabsch(np.array([rn[k] for k in sh]),
                             np.array([mn[k] for k in sh]))
        L = ref.lig_xyz[perm] @ Rk.T + tk
        pocket_ca = float(np.sqrt(np.mean([
            np.sum((mn[k] - (rn[k] @ Rk.T + tk)) ** 2) for k in sh])))
        fg = [k for k in set(mn) & set(rn)
              if CYP3A4_FG_LOOP[0] <= k <= CYP3A4_FG_LOOP[1]]
        fg_ca = (float(np.sqrt(np.mean([np.sum((mn[k] - (rn[k] @ Rk.T + tk)) ** 2)
                                        for k in fg]))) if fg else None)
        bbi = [i for i, (_c, _r, a) in enumerate(m.prot_key) if a in BACKBONE]
        dbb = np.linalg.norm(m.prot_xyz[bbi][:, None, :] - L[None, :, :], axis=2).min(1)
        j = int(np.argmin(dbb))
        extra.append(dict(ligand=lg, pocket_ca_rmsd=round(pocket_ca, 2),
                          fg_loop_ca_rmsd=None if fg_ca is None else round(fg_ca, 2),
                          n_fg_residues=len(fg),
                          min_backbone_contact=round(float(dbb.min()), 2),
                          worst_backbone=f"{m.prot_key[bbi[j]][1]}:{m.prot_key[bbi[j]][2]}"))
    res["unrepackable_diagnosis"] = extra

    # Is the 2.6 A self-clash criterion reasonable? Measure what the model's OWN pocket
    # side chains already satisfy. If they routinely sit below 2.6 A the criterion is too
    # strict and would manufacture failures; `self-clash-blocks` above says how often it
    # actually bound.
    from scipy.spatial import cKDTree
    nat = []
    for lg in hard:
        r = lig[lig.id == lg].iloc[0]
        ref, _c, _ch = load_reference(r.pdb, lg)
        pk = set(P.pocket_residues_from_structure(ref, radius=PREREG["align_radius_A"]))
        m = P.load_structure(POOL / f"{lg}__unsteered__s1" / "input_model_0.cif")
        m, _o, _i = renumber_to_reference(m, ref)
        for rr in pk:
            sel = [i for i, (_c2, rn2, a) in enumerate(m.prot_key)
                   if rn2 == rr and a not in BACKBONE and a != "CB"]
            oth = [i for i, (_c2, rn2, _a) in enumerate(m.prot_key)
                   if abs(rn2 - rr) > 1]
            if sel and oth:
                nat.append(float(cKDTree(m.prot_xyz[oth]).query(
                    m.prot_xyz[sel], k=1)[0].min()))
    res["self_clash_control"] = {
        "n_native_pocket_side_chains": len(nat),
        "median_min_distance_to_rest_of_protein_A": float(np.median(nat)),
        "p5_A": float(np.percentile(nat, 5)),
        "frac_below_2.6A": float(np.mean(np.array(nat) < PREREG["self_clash_cut_A"])),
    }
    out = REPO / "data" / "processed" / "side_chain_q4_repack.json"
    out.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["measure", "analyse", "plasticity",
                                      "templates", "repack"])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    if a.stage == "measure":
        cmd_measure(a.workers, a.limit)
    elif a.stage == "analyse":
        cmd_analyse()
    elif a.stage == "plasticity":
        cmd_plasticity()
    elif a.stage == "templates":
        cmd_templates()
    else:
        cmd_repack(a.workers)


if __name__ == "__main__":
    main()
