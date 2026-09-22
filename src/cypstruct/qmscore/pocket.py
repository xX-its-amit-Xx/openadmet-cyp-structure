"""CYP3A4-specific pose terms ABOVE the iron: the part of the scorer that is not saturated.

`geometry.py` implements the first-coordination-sphere terms (design sections A, B, C).
Those are **inert** and are deliberately not rebuilt here — five independent measurements
say the iron anchor is already reproduced by unbonded co-folding (FINDING 008, 019, 021).
What is left, and what this module implements, is everything the ligand does with the part
of the pocket that actually differs between a right pose and a wrong one:

* **section E, residue-specific polar anchors** — Ser119, Arg106, Arg212, Asp214, plus
  Thr224 and backbone amides from `docs/worldmodel/CYP3A4_EVOLUTION.md` §4.4/§5. Scored at
  those positions specifically, never as a generic H-bond count.
* **the phenylalanine cluster** — Phe57/108/213/215/219/220/241/304, the hydrophobic roof
  over the distal face. Absent from the design document entirely; it is the most
  CYP3A4-specific recognition feature there is, and the evolutionary decomposition names
  aromatic stacking as the primary recognition mode.
* **the F/G roof** — residues 202-260, the region OpenADMET names as the remodelling
  hotspot and the reason co-folding underperforms. Engagement is measured as contact
  counts and depth, never as an RMSD to anything (that would be leaky).
* **section D, strain and clash** — included so the PXR replication is logged. The PXR
  campaign built an MMFF strain validator, tested it and killed it; this is expected to
  gate null and is here to be *measured* doing so, not to be made to work.

Every quantity is computed from the PREDICTION ALONE. Nothing here reads the reference
structure, including the residue numbering, which is resolved against the canonical
UniProt sequence rather than against the crystal (see `resolve_offset`).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np

# --------------------------------------------------------------------------
# The residues. UniProt P08684 (full-length CYP3A4) numbering throughout.
# --------------------------------------------------------------------------

# Section E of the design, plus the two the design predates.
#   Ser119  - midazolam imidazole N, ritonavir; the canonical CYP3A4 H-bond (5TE8).
#   Arg106  - DHEA-S sulfate in one chain of 8GK3; B-C loop rim.
#   Arg212  - confirmed by mutagenesis as the active-site caffeine anchor (8SO1/8SO2).
#   Asp214  - design section E.
#   Thr224  - intrachannel anchor, mutagenesis-confirmed (8SO1/8SO2); fluorol.
POLAR_ANCHORS: dict[int, str] = {106: "ARG", 119: "SER", 212: "ARG",
                                 214: "ASP", 224: "THR"}

# Side-chain polar atoms per anchor type, and whether the side chain can donate/accept.
_ANCHOR_ATOMS: dict[str, tuple[tuple[str, ...], bool, bool]] = {
    "ARG": (("NE", "NH1", "NH2"), True, False),
    "SER": (("OG",), True, True),
    "THR": (("OG1",), True, True),
    "ASP": (("OD1", "OD2"), False, True),
    "ASN": (("OD1", "ND2"), True, True),
    "GLU": (("OE1", "OE2"), False, True),
    "LYS": (("NZ",), True, False),
    "TYR": (("OH",), True, True),
}

# The hydrophobic roof. Williams 2004 and every CYP3A4 structure since.
# Phe304 is the gatekeeper - different rotamers for different ligand sizes.
PHE_CLUSTER: tuple[int, ...] = (57, 108, 213, 215, 219, 220, 241, 304)
PHE_GATEKEEPER = 304
_PHE_RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")

# F-helix -> F' -> G' -> G. The plastic region; where the evolutionary argument says
# discrimination lives.
FG_SPAN: tuple[int, int] = (202, 260)

# H-bond geometry. A well centred on 2.9 A, zero outside 2.4-3.6 A.
_HB_OPT, _HB_HALF = 2.9, 0.7

# Aromatic ring geometry for pi contacts.
_RING_MAX = 6.5          # centroid-centroid beyond this is not a stacking contact
_PD_ANGLE = 35.0         # parallel-displaced: interplanar angle below this
_T_ANGLE = 55.0          # T-shaped / edge-to-face: interplanar angle above this

VDW = {"H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
       "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98, "FE": 2.00, "B": 1.92,
       "SI": 2.10, "SE": 1.90}

_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


# --------------------------------------------------------------------------
# numbering
# --------------------------------------------------------------------------


def resolve_offset(prot_res: dict[tuple[str, int], str], canon_seq: str,
                   span: int = 120) -> tuple[int, float]:
    """Offset `k` such that `model_resnum + k` is the canonical (UniProt) number.

    Co-folders number 1..N from the input construct; the CYP3A4 literature, the PDB and
    every residue name in this module use full-length P08684 numbering, which for the
    usual crystallography construct differs by about +22. Getting this wrong does not
    fail loudly - it silently scores the ligand against the wrong seven residues.

    The scan maximises residue-NAME agreement against the CANONICAL SEQUENCE, not against
    the reference structure. That matters for two reasons: it is available on a blind
    target, where no crystal exists, and it cannot manufacture agreement, because a wrong
    offset scores near-chance identity (~5%) on a 480-residue chain.

    Returns `(offset, identity)`; callers should refuse anything under ~0.8.
    """
    obs = {r: _THREE_TO_ONE.get(n.upper(), "X") for (_c, r), n in prot_res.items()}
    if not obs:
        return 0, 0.0
    best_k, best_n = 0, -1
    for k in range(-span, span + 1):
        agree = 0
        for r, aa in obs.items():
            j = r + k - 1
            if 0 <= j < len(canon_seq) and canon_seq[j] == aa:
                agree += 1
        if agree > best_n:
            best_k, best_n = k, agree
    return best_k, best_n / len(obs)


# --------------------------------------------------------------------------
# small geometric helpers
# --------------------------------------------------------------------------


def _hb_well(d: float) -> float:
    """1.0 at an ideal heavy-atom H-bond distance, 0 outside the window."""
    if not np.isfinite(d):
        return 0.0
    return float(max(0.0, 1.0 - abs(d - _HB_OPT) / _HB_HALF))


def _plane(pts: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """(centroid, unit normal, rms deviation) of the best-fit plane."""
    c = pts.mean(0)
    _u, s, vt = np.linalg.svd(pts - c)
    n = vt[-1] / (np.linalg.norm(vt[-1]) + 1e-12)
    rms = float(np.sqrt((((pts - c) @ n) ** 2).mean()))
    return c, n, rms


def _elem_from_name(atom_name: str) -> str:
    """Element from a PDB atom name. Protein atoms only, so first character suffices."""
    a = atom_name.strip().upper()
    return a[0] if a else "C"


def planar_rings(xyz: np.ndarray, elems: list[str],
                 tol: float = 0.45, flat: float = 0.15) -> list[list[int]]:
    """Planar 5- and 6-membered rings, perceived from coordinates alone.

    Deliberately NOT read off an RDKit molecule built from the manifest SMILES: that
    would tie the term to an input the blind pipeline supplies separately, and bond
    perception from the pose's own coordinates is what a scorer sees at inference time.
    Planarity is the aromaticity test - a chair cyclohexane fails it, a phenyl passes.
    """
    n = len(xyz)
    if n < 5:
        return []
    cov = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "F": 0.57, "CL": 1.02,
           "BR": 1.20, "I": 1.39, "P": 1.07, "B": 0.84, "SE": 1.20, "SI": 1.11}
    r = np.array([cov.get(e.upper(), 0.77) for e in elems])
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    adj = (d < (r[:, None] + r[None, :] + tol)) & (d > 0.1)
    try:
        import networkx as nx
    except ImportError:            # pragma: no cover
        return []
    g = nx.from_numpy_array(adj.astype(int))
    g.remove_edges_from(nx.selfloop_edges(g))
    out, seen = [], set()
    for cyc in nx.minimum_cycle_basis(g):
        if len(cyc) not in (5, 6):
            continue
        key = tuple(sorted(cyc))
        if key in seen:
            continue
        pts = xyz[list(cyc)]
        _c, _n, rms = _plane(pts)
        if rms <= flat:
            seen.add(key)
            out.append(list(cyc))
    return out


def _ring_pair(c1, n1, c2, n2) -> tuple[float, float, float]:
    """(centroid distance, interplanar angle in [0,90], lateral offset)."""
    v = c2 - c1
    dist = float(np.linalg.norm(v))
    ang = float(np.degrees(np.arccos(np.clip(abs(float(np.dot(n1, n2))), -1, 1))))
    # offset measured on ring 1's axis, the conventional definition
    axial = float(abs(np.dot(v, n1)))
    off = float(np.sqrt(max(dist ** 2 - axial ** 2, 0.0)))
    return dist, ang, off


# --------------------------------------------------------------------------
# the terms
# --------------------------------------------------------------------------


@dataclass
class PocketTerms:
    """Everything this module measures for one pose. Distances A, angles degrees."""

    # --- E. residue-specific polar anchors ---
    hb_ser119: float = 0.0
    hb_arg106: float = 0.0
    hb_arg212: float = 0.0
    hb_asp214: float = 0.0
    hb_thr224: float = 0.0
    d_ser119: float = np.nan
    d_arg106: float = np.nan
    d_arg212: float = np.nan
    d_asp214: float = np.nan
    d_thr224: float = np.nan
    anchor_score: float = 0.0            # sum of the five wells
    n_anchor_hb: int = 0                 # how many are satisfied at all
    hb_backbone: float = 0.0             # best backbone N/O H-bond well
    n_hb_backbone: int = 0
    unsat_buried_polar: int = 0          # buried ligand N/O with no protein partner

    # --- the phenylalanine cluster ---
    phe_n_rings: int = 0                 # ligand planar rings found
    phe_n_contacts: int = 0              # (lig ring, Phe ring) pairs inside 6.5 A
    phe_n_pd: int = 0                    # parallel-displaced
    phe_n_t: int = 0                     # T-shaped / edge-to-face
    phe_stack_score: float = 0.0         # summed smooth geometric score
    phe_best_score: float = 0.0
    phe_min_dist: float = np.nan         # best ring-centroid distance in the cluster
    phe_n_res: int = 0                   # distinct cluster Phe engaged
    phe_heavy_contacts: int = 0          # heavy-atom pairs < 4.5 A to cluster side chains
    f304_dist: float = np.nan            # gatekeeper standoff

    # --- the F/G roof ---
    fg_contacts: int = 0
    fg_sc_contacts: int = 0
    fg_frac: float = np.nan              # share of all protein contacts that are F/G
    fg_min_dist: float = np.nan
    fg_n_res: int = 0
    fg_depth: float = np.nan             # ligand centroid height above the heme plane

    # --- D. strain and clash ---
    max_clash: float = np.nan            # worst vdW overlap with protein
    n_clash: int = 0
    self_clash: float = np.nan           # worst 1-4+ intramolecular overlap
    strain: float = np.nan               # MMFF94s local strain, kcal/mol
    strain_ok: bool = False

    # --- bookkeeping, never a feature ---
    resnum_offset: int = 0
    seq_identity: float = 0.0
    n_anchors_present: int = 0
    n_phe_present: int = 0

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class _Res:
    name: str
    atoms: dict[str, np.ndarray] = field(default_factory=dict)


def _index_residues(prot_xyz, prot_key, prot_res, offset) -> dict[int, _Res]:
    """{canonical resnum: _Res} for the single protein chain of a prediction."""
    out: dict[int, _Res] = {}
    for i, (c, r, a) in enumerate(prot_key):
        rn = r + offset
        res = out.get(rn)
        if res is None:
            res = out[rn] = _Res(name=prot_res.get((c, r), "UNK").upper())
        res.atoms[a.strip().upper()] = prot_xyz[i]
    return out


def compute(lig_xyz: np.ndarray, lig_elem: list[str],
            prot_xyz: np.ndarray, prot_key: list[tuple[str, int, str]],
            prot_res: dict[tuple[str, int], str],
            canon_seq: str,
            fe: np.ndarray | None = None,
            heme_xyz: np.ndarray | None = None,
            heme_atom: list[str] | None = None,
            heme_elem: list[str] | None = None,
            axial_sg: np.ndarray | None = None,
            smiles: str | None = None,
            offset: int | None = None) -> PocketTerms:
    """All above-the-iron terms for one pose.

    `canon_seq` is the canonical UniProt sequence of the target; it fixes the residue
    numbering and is the only input that is not the prediction itself. Supplying
    `offset` skips the scan (useful when a whole pool shares one construct).
    """
    t = PocketTerms()
    lig_xyz = np.asarray(lig_xyz, float)
    prot_xyz = np.asarray(prot_xyz, float)
    if len(lig_xyz) == 0 or len(prot_xyz) == 0:
        return t

    if offset is None:
        offset, ident = resolve_offset(prot_res, canon_seq)
    else:
        _k, ident = resolve_offset(prot_res, canon_seq)
    t.resnum_offset, t.seq_identity = int(offset), float(ident)

    res = _index_residues(prot_xyz, prot_key, prot_res, offset)
    els = [e.upper() for e in lig_elem]
    lig_polar = [i for i, e in enumerate(els) if e in ("N", "O", "S")]

    # ---------------- E. residue-specific polar anchors -------------------
    wells = {}
    for num, expect in POLAR_ANCHORS.items():
        r = res.get(num)
        if r is None or r.name != expect:
            continue                        # not this protein, or a mutant; stay 0/NaN
        t.n_anchors_present += 1
        names, _don, _acc = _ANCHOR_ATOMS[expect]
        pts = [r.atoms[a] for a in names if a in r.atoms]
        if not pts or not lig_polar:
            continue
        P = np.asarray(pts)
        d = np.linalg.norm(lig_xyz[lig_polar][:, None, :] - P[None, :, :], axis=2)
        dm = float(d.min())
        wells[num] = _hb_well(dm)
        setattr(t, f"d_{expect.lower()[:3]}{num}", dm)
    # named fields, in the design's own order
    for num, attr in ((106, "arg106"), (119, "ser119"), (212, "arg212"),
                      (214, "asp214"), (224, "thr224")):
        setattr(t, f"hb_{attr}", float(wells.get(num, 0.0)))
    t.anchor_score = float(sum(wells.values()))
    t.n_anchor_hb = int(sum(1 for v in wells.values() if v > 0))

    # backbone amides and carbonyls - the DHEA-S/Gly481 motif, kept separate from
    # the named side chains so a generic H-bond cannot be mistaken for a specific one
    bb_idx = [i for i, (_c, _r, a) in enumerate(prot_key)
              if a.strip().upper() in ("N", "O")]
    if bb_idx and lig_polar:
        d = np.linalg.norm(lig_xyz[lig_polar][:, None, :]
                           - prot_xyz[bb_idx][None, :, :], axis=2)
        t.hb_backbone = _hb_well(float(d.min()))
        t.n_hb_backbone = int((np.abs(d - _HB_OPT) < _HB_HALF).any(1).sum())

    # buried but unsatisfied ligand polar atoms - the classic desolvation penalty,
    # here restricted to atoms the pose has actually buried
    prot_polar = [i for i, (_c, _r, a) in enumerate(prot_key)
                  if _elem_from_name(a) in ("N", "O", "S")]
    if lig_polar:
        dall = np.linalg.norm(lig_xyz[:, None, :] - prot_xyz[None, :, :], axis=2)
        nbr = (dall < 4.5).sum(1)
        if prot_polar:
            dpol = dall[:, prot_polar]
            for i in lig_polar:
                if nbr[i] >= 8 and dpol[i].min() > 3.6:
                    t.unsat_buried_polar += 1
    else:
        dall = np.linalg.norm(lig_xyz[:, None, :] - prot_xyz[None, :, :], axis=2)

    # ---------------- the phenylalanine cluster ---------------------------
    rings = planar_rings(lig_xyz, els)
    t.phe_n_rings = len(rings)
    lig_rings = []
    for ring in rings:
        c, n, _rms = _plane(lig_xyz[list(ring)])
        lig_rings.append((c, n))

    phe_sc_idx: list[int] = []
    best_d, engaged = np.inf, set()
    for num in PHE_CLUSTER:
        r = res.get(num)
        if r is None or r.name != "PHE":
            continue
        t.n_phe_present += 1
        pts = [r.atoms[a] for a in _PHE_RING if a in r.atoms]
        if len(pts) < 5:
            continue
        pc, pn, _rms = _plane(np.asarray(pts))
        phe_sc_idx.extend(
            i for i, (_c, rr, a) in enumerate(prot_key)
            if rr + offset == num and a.strip().upper() not in ("N", "CA", "C", "O"))
        for lc, ln in lig_rings:
            dist, ang, off = _ring_pair(lc, ln, pc, pn)
            if dist > _RING_MAX:
                continue
            t.phe_n_contacts += 1
            engaged.add(num)
            best_d = min(best_d, dist)
            s = 0.0
            if ang <= _PD_ANGLE and dist <= 5.5 and off <= 3.0:
                t.phe_n_pd += 1
                s = (1 - abs(dist - 3.9) / 1.6) * (1 - ang / _PD_ANGLE)
            elif ang >= _T_ANGLE and 4.2 <= dist <= 6.2:
                t.phe_n_t += 1
                s = (1 - abs(dist - 5.2) / 1.0) * ((ang - _T_ANGLE) / (90 - _T_ANGLE))
            s = float(max(s, 0.0))
            t.phe_stack_score += s
            t.phe_best_score = max(t.phe_best_score, s)
        if num == PHE_GATEKEEPER:
            sc = [r.atoms[a] for a in _PHE_RING + ("CB",) if a in r.atoms]
            if sc:
                t.f304_dist = float(np.linalg.norm(
                    lig_xyz[:, None, :] - np.asarray(sc)[None, :, :], axis=2).min())
    t.phe_min_dist = float(best_d) if np.isfinite(best_d) else np.nan
    t.phe_n_res = len(engaged)
    if phe_sc_idx:
        t.phe_heavy_contacts = int((dall[:, phe_sc_idx] < 4.5).sum())

    # ---------------- the F/G roof ---------------------------------------
    lo, hi = FG_SPAN
    fg_idx = [i for i, (_c, r, _a) in enumerate(prot_key) if lo <= r + offset <= hi]
    fg_sc = [i for i, (_c, r, a) in enumerate(prot_key)
             if lo <= r + offset <= hi and a.strip().upper() not in ("N", "CA", "C", "O")]
    n_contacts_all = int((dall < 4.5).sum())
    if fg_idx:
        dfg = dall[:, fg_idx]
        t.fg_contacts = int((dfg < 4.5).sum())
        t.fg_min_dist = float(dfg.min())
        t.fg_n_res = len({prot_key[fg_idx[j]][1]
                          for j in np.unique(np.nonzero(dfg < 4.5)[1])})
        t.fg_frac = float(t.fg_contacts / n_contacts_all) if n_contacts_all else np.nan
    if fg_sc:
        t.fg_sc_contacts = int((dall[:, fg_sc] < 4.5).sum())

    # depth: how far the ligand centroid sits above the porphyrin, i.e. how far up
    # toward the roof it has been placed. Distal is positive.
    if heme_xyz is not None and len(heme_xyz):
        try:
            from .geometry import porphyrin_frame
            f2, normal = porphyrin_frame(np.asarray(heme_xyz, float),
                                         list(heme_atom or []), axial_sg,
                                         list(heme_elem or []))
            if f2 is not None and normal is not None:
                t.fg_depth = float(np.dot(lig_xyz.mean(0) - f2, normal))
        except Exception:                                   # noqa: BLE001
            pass

    # ---------------- D. strain and clash --------------------------------
    lig_r = np.array([VDW.get(e, 1.7) for e in els])
    prot_r = np.array([VDW.get(_elem_from_name(a), 1.7) for (_c, _r, a) in prot_key])
    overlap = (lig_r[:, None] + prot_r[None, :]) - dall
    t.max_clash = float(overlap.max())
    t.n_clash = int((overlap > 0.5).sum())

    dl = np.linalg.norm(lig_xyz[:, None, :] - lig_xyz[None, :, :], axis=2)
    np.fill_diagonal(dl, np.inf)
    far = dl > 2.9            # crude 1-4 exclusion: bonded/angle pairs are closer
    t.self_clash = float(np.max(
        np.where(far, (lig_r[:, None] + lig_r[None, :]) - dl, -np.inf)))

    if smiles:
        s = mmff_strain(smiles, lig_xyz, els)
        if s is not None:
            t.strain = float(s)
            t.strain_ok = bool(s < 10.0)
    return t


# --------------------------------------------------------------------------
# strain — design section D, expected to gate null
# --------------------------------------------------------------------------


def mmff_strain(smiles: str, lig_xyz: np.ndarray, lig_elem: list[str]) -> float | None:
    """Local MMFF94s strain: E(pose, hydrogens relaxed) - E(nearest local minimum).

    Local rather than global on purpose. A global strain needs a conformer search and
    turns a pose property into a molecule property with noise attached; the quantity
    that is supposed to discriminate is how far the *posed* torsions sit from the
    nearest relaxed geometry.

    The PXR campaign's verdict on the term it replicates: rejected, with the note that
    DFT would have ordered the energies identically. Returns None rather than a
    guess whenever the molecule cannot be matched to the pose.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem
    except ImportError:                                     # pragma: no cover
        return None
    RDLogger.DisableLog("rdApp.*")
    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        return None
    ref = Chem.RemoveHs(ref)
    if ref.GetNumAtoms() != len(lig_xyz):
        return None

    perm = _match_smiles_to_pose(ref, lig_xyz, lig_elem)
    if perm is None:
        return None

    mol = Chem.RWMol(ref)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i_mol, i_pose in perm.items():
        x, y, z = lig_xyz[i_pose]
        conf.SetAtomPosition(i_mol, (float(x), float(y), float(z)))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    m = mol.GetMol()
    try:
        Chem.SanitizeMol(m)
        m = Chem.AddHs(m, addCoords=True)
        props = AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94s")
        if props is None:
            return None
        # relax hydrogens only: the pose carries no H, so their placement is ours
        ff = AllChem.MMFFGetMoleculeForceField(m, props)
        for a in m.GetAtoms():
            if a.GetAtomicNum() > 1:
                ff.AddFixedPoint(a.GetIdx())
        ff.Minimize(maxIts=400)
        e_pose = ff.CalcEnergy()

        ff2 = AllChem.MMFFGetMoleculeForceField(m, props)
        ff2.Minimize(maxIts=2000)
        e_min = ff2.CalcEnergy()
    except Exception:                                       # noqa: BLE001
        return None
    return float(e_pose - e_min)


def _match_smiles_to_pose(ref, lig_xyz, lig_elem) -> dict[int, int] | None:
    """{smiles atom index -> pose atom index} by element-labelled graph isomorphism."""
    try:
        import networkx as nx
        from networkx.algorithms.isomorphism import GraphMatcher, categorical_node_match
    except ImportError:                                     # pragma: no cover
        return None
    cov = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "F": 0.57, "CL": 1.02,
           "BR": 1.20, "I": 1.39, "P": 1.07, "B": 0.84, "SE": 1.20, "SI": 1.11}
    els = [e.upper() for e in lig_elem]
    r = np.array([cov.get(e, 0.77) for e in els])
    xyz = np.asarray(lig_xyz, float)
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    adj = (d < (r[:, None] + r[None, :] + 0.45)) & (d > 0.1)
    gp = nx.from_numpy_array(adj.astype(int))
    gp.remove_edges_from(nx.selfloop_edges(gp))
    for i, e in enumerate(els):
        gp.nodes[i]["el"] = e

    gs = nx.Graph()
    for a in ref.GetAtoms():
        gs.add_node(a.GetIdx(), el=a.GetSymbol().upper())
    for b in ref.GetBonds():
        gs.add_edge(b.GetBeginAtomIdx(), b.GetEndAtomIdx())

    gm = GraphMatcher(gs, gp, node_match=categorical_node_match("el", ""))
    try:
        it = gm.isomorphisms_iter()
        return dict(next(it))
    except StopIteration:
        return None
    except Exception:                                       # noqa: BLE001
        return None


# --------------------------------------------------------------------------
# the feature block the evaluation uses
# --------------------------------------------------------------------------

#: Feature groups, so that a null term is visible as a null rather than hidden in a sum.
GROUPS: dict[str, list[str]] = {
    "polar": ["hb_ser119", "hb_arg106", "hb_arg212", "hb_asp214", "hb_thr224",
              "anchor_score", "n_anchor_hb", "hb_backbone", "n_hb_backbone",
              "unsat_buried_polar", "d_ser119", "d_arg106", "d_arg212",
              "d_asp214", "d_thr224"],
    "phe": ["phe_n_contacts", "phe_n_pd", "phe_n_t", "phe_stack_score",
            "phe_best_score", "phe_min_dist", "phe_n_res", "phe_heavy_contacts",
            "f304_dist"],
    "fg": ["fg_contacts", "fg_sc_contacts", "fg_frac", "fg_min_dist", "fg_n_res",
           "fg_depth"],
    "strain": ["strain", "max_clash", "n_clash", "self_clash"],
}

FEATURES: list[str] = [f for g in GROUPS.values() for f in g]
