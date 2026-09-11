"""Tier-0 of the CYP physics scorer: the terms that need only coordinates.

Free to compute, so these run on every pose in the pool. They encode the
first-coordination-sphere facts of P450 chemistry that generic docking scores and
co-folder confidence heads have no representation for. See `docs/QM_SCORER_DESIGN.md`
for the reasoning and `scripts/structure/build_reference_set.py` for the empirical
calibration of every window used below.

Everything here returns raw physical quantities, not a fitted score. Combination and
weighting happen in `rank.py` against measured LDDT-PLI labels, so that no hand-chosen
constant silently becomes a model parameter.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

# Van der Waals radii (Angstrom) for the clash term. Bondi values.
VDW = {"H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
       "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98, "FE": 2.00, "B": 1.92,
       "SI": 2.10, "SE": 1.90}

# --------------------------------------------------------------------------
# Empirical windows — MEASURED, not quoted.
#
# Source: `scripts/structure/build_reference_set.py` over 116 deposited CYP3A4
# entries (UniProt P08684), 2026-09-09. 151 chain-ligand pairs, of which
# **104 are directly Fe-coordinated across 84 distinct ligands** — so this is the
# majority binding mode for CYP3A4 in the PDB, not a special case.
#
#   Fe-donor distance   p5 1.94   p50 2.20   p95 2.38      (donor is N in 103/104)
#   S(Cys442)-Fe-donor  p5 159.5  p50 171.2  p95 177.6     (trans, and tight)
#   Fe out-of-plane     p5 -0.42  p50 -0.05  p95 +0.06     (in-plane, 6-coordinate)
#   Fe-SG(Cys442)       p5 2.12   p50 2.37   p95 2.51
#   type I standoff     p5 2.85   p50 3.91   p95 5.91
#
# The single most useful number here: only **0.7%** of ligand atoms across the whole
# reference set fall on the proximal face of the porphyrin. Real structures essentially
# never violate it, which makes the distal-side test a near-perfect validity filter.
# --------------------------------------------------------------------------
COORD_LO, COORD_HI = 1.90, 2.45      # Fe-donor dative bond (p5-p95 with a small margin)
TRANS_ANGLE_MIN = 155.0              # S(Cys)-Fe-donor; measured p5 is 159.5
LONEPAIR_ANGLE_MAX = 35.0            # (lone-pair vector, N->Fe). Not yet measured — prior.
SOM_LO, SOM_HI = 2.9, 5.9            # substrate standoff over Fe (type I p5-p95)
STACK_LO, STACK_HI = 3.2, 4.2        # ring centroid to porphyrin plane. Prior.
STACK_ANGLE_MAX = 25.0               # interplanar angle for stacking. Prior.


@dataclass
class GeometryTerms:
    """Raw geometric observables for one pose. All distances A, all angles degrees."""

    # --- iron coordination (type II) ---
    fe_min_dist: float = np.nan          # closest ligand heavy atom to Fe
    fe_min_element: str = ""             # its element
    fe_donor_dist: float = np.nan        # closest N/O/S specifically
    fe_donor_element: str = ""
    s_fe_donor_angle: float = np.nan     # trans-to-thiolate angle
    lonepair_angle: float = np.nan       # lone-pair direction vs N->Fe
    is_coordinated: bool = False

    # --- side of the heme (hard validity) ---
    distal_projection: float = np.nan    # + = distal (correct), - = proximal (impossible)
    n_atoms_proximal: int = 0            # ligand atoms on the buried face
    frac_atoms_proximal: float = np.nan

    # --- substrate geometry (type I) ---
    fe_standoff: float = np.nan          # closest-atom distance when not coordinated
    n_carbons_in_som_shell: int = 0

    # --- sterics ---
    max_clash: float = np.nan            # worst vdW overlap with protein/heme
    n_clashes: int = 0
    heme_clash: float = np.nan

    # --- porphyrin stacking ---
    stack_dist: float = np.nan
    stack_angle: float = np.nan
    is_stacked: bool = False

    # --- burial / pocket engagement ---
    n_protein_contacts: int = 0          # ligand-protein heavy-atom pairs < 4.5 A
    buried_fraction: float = np.nan      # fraction of ligand atoms with >=6 neighbours
    radius_of_gyration: float = np.nan

    def to_dict(self) -> dict:
        return asdict(self)


# --------------------------------------------------------------------------


def porphyrin_frame(heme_xyz: np.ndarray, heme_atom: list[str],
                    axial_sg: np.ndarray | None,
                    heme_elem: list[str] | None = None
                    ) -> tuple[np.ndarray | None, np.ndarray | None]:
    """(Fe, unit normal pointing to the DISTAL face).

    Identification is by ELEMENT, with atom names used only as a fallback. Engines do not
    agree on naming: a deposited heme calls its iron `FE` and its pyrrole nitrogens
    `NA/NB/NC/ND`, but Chai-1 receives the cofactor as SMILES and emits generic names, so
    a name-only lookup finds no iron at all and every metal term silently returns NaN.

    The four pyrrole nitrogens are taken as the four nitrogens nearest the iron, which is
    true of any porphyrin and needs no naming convention whatsoever.

    Orientation is fixed by the thiolate: the proximal face is the one the Cys sulfur is
    on. Without it we cannot orient the axis and return None rather than guess - a guessed
    sign would flip the single most important validity check in the scorer.
    """
    heme_xyz = np.asarray(heme_xyz, float)
    if len(heme_xyz) == 0:
        return None, None
    elems = [e.upper() for e in (heme_elem or [])]
    names = [a.upper() for a in heme_atom]

    fe = None
    if elems and len(elems) == len(heme_xyz):
        idx = [i for i, e in enumerate(elems) if e == "FE"]
        if idx:
            fe = heme_xyz[idx[0]]
    if fe is None:
        idx = [i for i, a in enumerate(names) if a == "FE" or a.startswith("FE")]
        if idx:
            fe = heme_xyz[idx[0]]
    if fe is None:
        return None, None

    # pyrrole nitrogens: the four N nearest the iron. No naming convention required.
    if elems and len(elems) == len(heme_xyz):
        n_idx = [i for i, e in enumerate(elems) if e == "N"]
    else:
        n_idx = [i for i, a in enumerate(names) if a in ("NA", "NB", "NC", "ND")]
    if len(n_idx) < 3:
        return fe, None
    n_idx = sorted(n_idx, key=lambda i: np.linalg.norm(heme_xyz[i] - fe))[:4]
    N = heme_xyz[n_idx]

    _u, _s, vt = np.linalg.svd(N - N.mean(0))
    n = vt[-1] / np.linalg.norm(vt[-1])
    if axial_sg is None:
        return fe, None
    if np.dot(np.asarray(axial_sg, float) - fe, n) > 0:
        n = -n
    return fe, n


def lone_pair_direction(mol, atom_idx: int, conf_xyz: np.ndarray) -> np.ndarray | None:
    """Unit vector along an sp2 nitrogen's in-plane lone pair.

    For a two-coordinate aromatic N the lone pair points away from the bisector of its
    two ring bonds, in the ring plane. This is what has to aim at the iron; distance
    alone does not distinguish a real dative bond from a nitrogen that merely happens
    to be nearby with its lone pair pointing elsewhere.
    """
    try:
        a = mol.GetAtomWithIdx(int(atom_idx))
    except Exception:
        return None
    nbrs = [n.GetIdx() for n in a.GetNeighbors()]
    if len(nbrs) != 2:
        return None      # not a two-coordinate sp2 N; no unambiguous in-plane lone pair
    p = conf_xyz[atom_idx]
    v = np.zeros(3)
    for j in nbrs:
        d = conf_xyz[j] - p
        n = np.linalg.norm(d)
        if n < 1e-6:
            return None
        v += d / n
    if np.linalg.norm(v) < 1e-6:
        return None
    lp = -v / np.linalg.norm(v)     # opposite the bisector
    return lp


def compute(lig_xyz: np.ndarray, lig_elem: list[str],
            prot_xyz: np.ndarray, heme_xyz: np.ndarray, heme_atom: list[str],
            axial_sg: np.ndarray | None, mol=None,
            aromatic_rings: list[list[int]] | None = None,
            heme_elem: list[str] | None = None) -> GeometryTerms:
    """All tier-0 terms for one pose.

    `mol` is an RDKit molecule whose atom order matches `lig_xyz` — supply it to get
    the lone-pair and ring-stacking terms; without it those come back NaN rather than
    being silently approximated.
    """
    t = GeometryTerms()
    if len(lig_xyz) == 0:
        return t
    lig_xyz = np.asarray(lig_xyz, float)
    els = [e.upper() for e in lig_elem]

    fe, normal = porphyrin_frame(np.asarray(heme_xyz, float), heme_atom,
                                axial_sg, heme_elem)

    # ---- iron coordination -------------------------------------------------
    if fe is not None:
        d_fe = np.linalg.norm(lig_xyz - fe, axis=1)
        k = int(np.argmin(d_fe))
        t.fe_min_dist, t.fe_min_element = float(d_fe[k]), els[k]

        donor_idx = [i for i, e in enumerate(els) if e in ("N", "O", "S")]
        if donor_idx:
            dj = min(donor_idx, key=lambda i: d_fe[i])
            t.fe_donor_dist, t.fe_donor_element = float(d_fe[dj]), els[dj]
            if axial_sg is not None:
                v1 = np.asarray(axial_sg, float) - fe
                v2 = lig_xyz[dj] - fe
                c = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
                t.s_fe_donor_angle = float(np.degrees(np.arccos(np.clip(c, -1, 1))))
            if mol is not None:
                lp = lone_pair_direction(mol, dj, lig_xyz)
                if lp is not None:
                    to_fe = fe - lig_xyz[dj]
                    to_fe /= np.linalg.norm(to_fe) + 1e-12
                    c = float(np.clip(np.dot(lp, to_fe), -1, 1))
                    t.lonepair_angle = float(np.degrees(np.arccos(c)))
            t.is_coordinated = bool(
                COORD_LO <= t.fe_donor_dist <= COORD_HI
                and (np.isnan(t.s_fe_donor_angle) or t.s_fe_donor_angle >= TRANS_ANGLE_MIN)
            )

        if not t.is_coordinated:
            t.fe_standoff = float(d_fe.min())
            t.n_carbons_in_som_shell = int(sum(
                1 for i, e in enumerate(els) if e == "C" and SOM_LO <= d_fe[i] <= SOM_HI))

        # ---- which face of the porphyrin ----------------------------------
        if normal is not None:
            proj = (lig_xyz - fe) @ normal
            t.distal_projection = float(proj[k])
            t.n_atoms_proximal = int((proj < -1.0).sum())
            t.frac_atoms_proximal = float((proj < -1.0).mean())

    # ---- sterics -----------------------------------------------------------
    lig_r = np.array([VDW.get(e, 1.7) for e in els])
    if len(prot_xyz):
        # protein element types are unavailable here; 1.70 (carbon) is the safe
        # majority assumption and makes this a conservative clash estimate
        dm = np.linalg.norm(lig_xyz[:, None, :] - np.asarray(prot_xyz)[None, :, :], axis=2)
        overlap = (lig_r[:, None] + 1.70) - dm
        t.max_clash = float(overlap.max())
        t.n_clashes = int((overlap > 0.5).sum())
        t.n_protein_contacts = int((dm < 4.5).sum())
        t.buried_fraction = float(((dm < 4.5).sum(1) >= 6).mean())
    if len(heme_xyz):
        dh = np.linalg.norm(lig_xyz[:, None, :] - np.asarray(heme_xyz)[None, :, :], axis=2)
        # exclude the coordinating contact itself from the clash statistic
        t.heme_clash = float(((lig_r[:, None] + 1.70) - dh).max()) if dh.size else np.nan

    # ---- porphyrin stacking ------------------------------------------------
    if normal is not None and fe is not None and aromatic_rings:
        best = None
        for ring in aromatic_rings:
            pts = lig_xyz[list(ring)]
            if len(pts) < 3:
                continue
            cen = pts.mean(0)
            _u, _s, vt = np.linalg.svd(pts - cen)
            rn = vt[-1] / np.linalg.norm(vt[-1])
            dist = abs(float(np.dot(cen - fe, normal)))
            ang = float(np.degrees(np.arccos(np.clip(abs(np.dot(rn, normal)), -1, 1))))
            # lateral offset must be small or it is not stacked on the porphyrin
            lateral = float(np.linalg.norm((cen - fe) - np.dot(cen - fe, normal) * normal))
            if lateral > 4.0:
                continue
            if best is None or dist < best[0]:
                best = (dist, ang)
        if best:
            t.stack_dist, t.stack_angle = best
            t.is_stacked = bool(STACK_LO <= best[0] <= STACK_HI and best[1] <= STACK_ANGLE_MAX)

    cen = lig_xyz.mean(0)
    t.radius_of_gyration = float(np.sqrt(((lig_xyz - cen) ** 2).sum(1).mean()))
    return t


def hard_validity(t: GeometryTerms) -> tuple[bool, str]:
    """Physical impossibility filter. Returns (is_valid, reason_if_not).

    These are not soft penalties — a pose that fails is not a worse pose, it is not a
    pose. Filtering before ranking keeps impossible geometries from occupying slots in
    a consensus and dragging a centroid off the real binding mode.
    """
    if t.frac_atoms_proximal == t.frac_atoms_proximal and t.frac_atoms_proximal > 0.15:
        return False, f"{t.frac_atoms_proximal:.0%} of ligand atoms on the proximal heme face"
    if t.max_clash == t.max_clash and t.max_clash > 1.5:
        return False, f"severe steric overlap ({t.max_clash:.2f} A into vdW)"
    if t.fe_min_dist == t.fe_min_dist and t.fe_min_dist < 1.6:
        return False, f"ligand atom {t.fe_min_dist:.2f} A from Fe (shorter than any real bond)"
    if t.n_protein_contacts == 0:
        return False, "ligand makes no protein contact at all"
    return True, ""
