"""Pose geometry: structure loading, alignment, and the challenge metrics.

The reason this module exists rather than a pile of ad-hoc parsers: in the PXR
campaign the single most expensive mistake was trusting cross-model metrics that
were computed in an unvalidated reference frame. Some engine exports landed ~20 A
off the crystal and every downstream agreement number was inflated about twofold
before anyone noticed. So all frame handling is centralised here, alignment is
always by explicit residue identity (never by array position), and `lddt_pli` is
superposition-free by construction so it cannot be fooled by a bad fit.

Metrics implemented:
  - `lddt_pli`  — the OpenADMET structure-track primary metric. Superposition-free
    local distance difference test over protein-ligand contacts.
  - `bisy_rmsd` — binding-site-superposed, symmetry-corrected ligand RMSD, the
    secondary metric on the PXR board.
  - `lddt_lp`   — the protein-side local distance test restricted to pocket atoms.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .targets import HEME_ALIASES, IGNORE_HET

# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------


@dataclass
class Complex:
    """A protein-ligand-heme complex in one flat, index-stable representation."""

    name: str
    # protein heavy atoms
    prot_xyz: np.ndarray                     # (N, 3)
    prot_key: list[tuple[str, int, str]]     # (chain, resnum, atom_name), parallel to prot_xyz
    prot_res: dict[tuple[str, int], str] = field(default_factory=dict)  # -> resname
    # the query ligand
    lig_xyz: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    lig_elem: list[str] = field(default_factory=list)
    lig_name: str = ""
    lig_chain: str = ""
    # heme
    fe: np.ndarray | None = None             # (3,)
    heme_xyz: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    heme_atom: list[str] = field(default_factory=list)
    heme_elem: list[str] = field(default_factory=list)   # elements, for name-free parsing
    # proximal thiolate sulfur (Cys442 SG in CYP3A4)
    axial_sg: np.ndarray | None = None

    def ca(self) -> dict[tuple[str, int], np.ndarray]:
        """{(chain, resnum): CA coordinate} — the alignment handle."""
        return {(c, r): self.prot_xyz[i]
                for i, (c, r, a) in enumerate(self.prot_key) if a == "CA"}

    def transformed(self, R: np.ndarray, t: np.ndarray) -> "Complex":
        """Rigid-body copy. Applies to every coordinate set so nothing desynchronises."""
        def tx(x):
            return x @ R.T + t if x is not None and len(x) else x
        return Complex(
            name=self.name, prot_xyz=tx(self.prot_xyz), prot_key=list(self.prot_key),
            prot_res=dict(self.prot_res), lig_xyz=tx(self.lig_xyz),
            lig_elem=list(self.lig_elem), lig_name=self.lig_name, lig_chain=self.lig_chain,
            fe=(self.fe @ R.T + t) if self.fe is not None else None,
            heme_xyz=tx(self.heme_xyz), heme_atom=list(self.heme_atom),
            heme_elem=list(self.heme_elem),
            axial_sg=(self.axial_sg @ R.T + t) if self.axial_sg is not None else None,
        )


def _looks_like_heme(res) -> bool:
    """Recognise an iron-porphyrin by composition, regardless of what it is called.

    A heme is an iron atom inside a ~43-atom macrocycle with four nitrogens. That is a
    far more reliable signature than the residue name, which varies by engine and by
    depositor: `HEM` from Boltz and most crystals, `HEC` in some CYP2C9 entries, and a
    generic `LIG2` from Chai-1, which has no CCD input and receives the cofactor as SMILES.

    Getting this wrong is not a small error. An unrecognised heme is also the largest HET
    group in the file, so it gets selected as the QUERY LIGAND and every distance, RMSD
    and LDDT-PLI for that structure silently describes the wrong molecule.
    """
    n_fe = n_n = n_at = 0
    for at in res:
        n_at += 1
        el = at.element.name.upper()
        if el == "FE":
            n_fe += 1
        elif el == "N":
            n_n += 1
    return n_fe == 1 and n_n >= 4 and n_at >= 30



def load_structure(path: str | Path, ligand_code: str | None = None,
                   ligand_chain: str | None = None, model: int = 0,
                   assembly_chain: str | None = None) -> Complex:
    """Read an mmCIF/PDB into a `Complex`.

    `ligand_code` pins which HET group is the query ligand. When omitted we take the
    largest non-heme, non-solvent HET group — correct for challenge outputs (one
    ligand) and for most holo depositions, but *pin it explicitly* for anything that
    matters, because "largest HET group" quietly picks a cryoprotectant on the rare
    entry where the real ligand is small and a PEG chain is long.

    `assembly_chain` restricts to one chain. The CYP3A4 cryoEM construct is a
    symmetric trimer, so a whole-file parse yields three copies of everything; taking
    all three would make every distance metric meaningless.
    """
    import gemmi

    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if len(st) == 0:
        raise ValueError(f"no models in {path}")
    md = st[model]

    prot_xyz, prot_key, prot_res = [], [], {}
    heme_xyz, heme_atom, heme_elem, fe = [], [], [], None
    axial_sg = None
    het: dict[tuple[str, str, int], list] = {}

    for chain in md:
        cid = chain.name
        if assembly_chain and cid != assembly_chain:
            continue
        for res in chain:
            rname = res.name.strip().upper()
            info = gemmi.find_tabulated_residue(rname)
            is_aa = bool(info and info.is_amino_acid())
            if is_aa:
                prot_res[(cid, res.seqid.num)] = rname
                for at in res:
                    if at.element == gemmi.Element("H"):
                        continue
                    prot_xyz.append([at.pos.x, at.pos.y, at.pos.z])
                    prot_key.append((cid, res.seqid.num, at.name.strip()))
                    if rname == "CYS" and at.name.strip() == "SG":
                        # keep the LAST one only if we have no better candidate; the
                        # true axial SG is disambiguated later by distance to Fe
                        if axial_sg is None:
                            axial_sg = np.array([at.pos.x, at.pos.y, at.pos.z])
            elif rname in HEME_ALIASES or _looks_like_heme(res):
                # Name-based detection is not sufficient. Different engines name the same
                # cofactor differently: Boltz emits the CCD code `HEM`, but Chai-1 - which
                # takes the heme as a SMILES ligand because it has no CCD input - emits
                # generic `LIG2`. With only the alias list, the heme went unrecognised AND
                # was then picked up as the query ligand, because it is the largest HET
                # group present. Every metric for that structure compared the wrong
                # molecule and looked merely bad rather than wrong.
                for at in res:
                    p = [at.pos.x, at.pos.y, at.pos.z]
                    heme_xyz.append(p)
                    heme_atom.append(at.name.strip())
                    heme_elem.append(at.element.name.upper())
                    if at.element.name.upper() == "FE":
                        fe = np.array(p, float)
            elif rname not in IGNORE_HET:
                key = (cid, rname, res.seqid.num)
                het.setdefault(key, [])
                for at in res:
                    if at.element == gemmi.Element("H"):
                        continue
                    het[key].append(([at.pos.x, at.pos.y, at.pos.z], at.element.name))

    # pick the query ligand
    lig_xyz, lig_elem, lig_name, lig_chain = np.zeros((0, 3)), [], "", ""
    cands = list(het.items())
    if ligand_code:
        cands = [kv for kv in cands if kv[0][1] == ligand_code.upper()]
    if ligand_chain:
        cands = [kv for kv in cands if kv[0][0] == ligand_chain]
    if cands:
        (cid, rname, _), atoms = max(cands, key=lambda kv: len(kv[1]))
        lig_xyz = np.array([a[0] for a in atoms], float)
        lig_elem = [a[1] for a in atoms]
        lig_name, lig_chain = rname, cid

    prot_xyz = np.array(prot_xyz, float) if prot_xyz else np.zeros((0, 3))
    heme_xyz = np.array(heme_xyz, float) if heme_xyz else np.zeros((0, 3))

    cx = Complex(name=Path(path).stem, prot_xyz=prot_xyz, prot_key=prot_key,
                 prot_res=prot_res, lig_xyz=lig_xyz, lig_elem=lig_elem,
                 lig_name=lig_name, lig_chain=lig_chain, fe=fe,
                 heme_xyz=heme_xyz, heme_atom=heme_atom, heme_elem=heme_elem)
    # resolve the real axial thiolate: the Cys SG closest to Fe (should be ~2.3 A)
    if fe is not None:
        sgs = [(i, k) for i, k in enumerate(prot_key)
               if k[2] == "SG" and prot_res.get((k[0], k[1])) == "CYS"]
        if sgs:
            d = [np.linalg.norm(prot_xyz[i] - fe) for i, _ in sgs]
            j = int(np.argmin(d))
            if d[j] < 4.0:
                cx.axial_sg = prot_xyz[sgs[j][0]]
    else:
        cx.axial_sg = axial_sg
    return cx


# --------------------------------------------------------------------------
# alignment
# --------------------------------------------------------------------------


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Optimal rigid transform mapping P onto Q. Returns (R, t, rmsd)."""
    if len(P) != len(Q) or len(P) < 3:
        raise ValueError(f"kabsch needs >=3 paired points, got {len(P)} and {len(Q)}")
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = Qc - R @ Pc
    rmsd = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))
    return R, t, rmsd


def align_by_residue(mobile: Complex, ref: Complex,
                     resnums: list[int] | None = None,
                     atom: str = "CA") -> tuple[Complex, float, int]:
    """Superpose `mobile` onto `ref` using residues matched BY NUMBER, not by order.

    Returns (aligned_copy, rmsd_of_the_fit, n_atoms_used).

    Matching by residue number is the whole point. Positional matching breaks the
    moment two structures have different chain contents, different modelled ranges,
    or a gap — and it breaks silently, producing a plausible-looking wrong frame.
    """
    mca, rca = mobile.ca(), ref.ca()
    if atom != "CA":
        def sel(cx):
            return {(c, r): cx.prot_xyz[i]
                    for i, (c, r, a) in enumerate(cx.prot_key) if a == atom}
        mca, rca = sel(mobile), sel(ref)

    mnum = {r: v for (_c, r), v in mca.items()}
    rnum = {r: v for (_c, r), v in rca.items()}
    shared = sorted(set(mnum) & set(rnum))
    if resnums is not None:
        keep = set(resnums)
        shared = [r for r in shared if r in keep]
    if len(shared) < 3:
        raise ValueError(f"only {len(shared)} shared residues for alignment")
    P = np.array([mnum[r] for r in shared])
    Q = np.array([rnum[r] for r in shared])
    R, t, rmsd = kabsch(P, Q)
    return mobile.transformed(R, t), rmsd, len(shared)


def pocket_residues_from_structure(cx: Complex, radius: float = 6.0) -> list[int]:
    """Residues with any heavy atom within `radius` of the ligand. Empirical, not a prior."""
    if len(cx.lig_xyz) == 0:
        return []
    d = np.linalg.norm(cx.prot_xyz[:, None, :] - cx.lig_xyz[None, :, :], axis=2).min(1)
    return sorted({cx.prot_key[i][1] for i in np.where(d <= radius)[0]})


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------

_LDDT_THRESHOLDS = (0.5, 1.0, 2.0, 4.0)


def lddt_pli(model: Complex, ref: Complex, inclusion_radius: float = 6.0,
             lig_perm: np.ndarray | None = None,
             thresholds: tuple[float, ...] = _LDDT_THRESHOLDS) -> float:
    """Local Distance Difference Test over protein-ligand contacts (the primary metric).

    Superposition-free: it compares the *set of interatomic distances* between ligand
    atoms and nearby protein atoms, so no alignment can flatter or penalise it. For
    every (ligand atom, protein atom) pair inside `inclusion_radius` in the REFERENCE,
    the pair is "preserved" at threshold T if |d_model - d_ref| < T. The score is the
    fraction preserved, averaged over the four thresholds.

    `lig_perm` maps model ligand atom i -> reference ligand atom lig_perm[i]. Pass the
    symmetry-corrected mapping from `best_ligand_mapping`; without it, a molecule with
    a symmetric substituent is scored against an arbitrary atom labelling and loses
    points for being correct.
    """
    if len(model.lig_xyz) == 0 or len(ref.lig_xyz) == 0:
        return float("nan")

    # pair protein atoms by (chain-insensitive) residue number + atom name
    ref_idx = {(r, a): i for i, (_c, r, a) in enumerate(ref.prot_key)}
    mod_idx = {(r, a): i for i, (_c, r, a) in enumerate(model.prot_key)}
    shared = sorted(set(ref_idx) & set(mod_idx))
    if not shared:
        return float("nan")
    ri = np.array([ref_idx[k] for k in shared])
    mi = np.array([mod_idx[k] for k in shared])

    if lig_perm is None:
        n = min(len(model.lig_xyz), len(ref.lig_xyz))
        lig_perm = np.arange(n)
    lm = model.lig_xyz[: len(lig_perm)]
    lr = ref.lig_xyz[lig_perm]

    d_ref = np.linalg.norm(lr[:, None, :] - ref.prot_xyz[ri][None, :, :], axis=2)
    d_mod = np.linalg.norm(lm[:, None, :] - model.prot_xyz[mi][None, :, :], axis=2)

    mask = d_ref < inclusion_radius
    if not mask.any():
        return float("nan")
    diff = np.abs(d_mod[mask] - d_ref[mask])
    return float(np.mean([(diff < t).mean() for t in thresholds]))


def lddt_lp(model: Complex, ref: Complex, inclusion_radius: float = 15.0,
            pocket: list[int] | None = None,
            thresholds: tuple[float, ...] = _LDDT_THRESHOLDS) -> float:
    """Protein-side LDDT restricted to pocket residues (the 'lddt_lp' column)."""
    if pocket is None:
        pocket = pocket_residues_from_structure(ref)
    keep = set(pocket)
    ref_idx = {(r, a): i for i, (_c, r, a) in enumerate(ref.prot_key) if r in keep}
    mod_idx = {(r, a): i for i, (_c, r, a) in enumerate(model.prot_key) if r in keep}
    shared = sorted(set(ref_idx) & set(mod_idx))
    if len(shared) < 4:
        return float("nan")
    ri = np.array([ref_idx[k] for k in shared])
    mi = np.array([mod_idx[k] for k in shared])
    dr = np.linalg.norm(ref.prot_xyz[ri][:, None] - ref.prot_xyz[ri][None, :], axis=2)
    dm = np.linalg.norm(model.prot_xyz[mi][:, None] - model.prot_xyz[mi][None, :], axis=2)
    iu = np.triu_indices(len(ri), k=1)
    m = dr[iu] < inclusion_radius
    if not m.any():
        return float("nan")
    diff = np.abs(dm[iu][m] - dr[iu][m])
    return float(np.mean([(diff < t).mean() for t in thresholds]))


# Covalent radii (Angstrom) for perceiving ligand connectivity from coordinates.
_COV_R = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
          "CL": 1.02, "BR": 1.20, "I": 1.39, "B": 0.84, "SI": 1.11, "SE": 1.20, "FE": 1.32}


def _adjacency(xyz: np.ndarray, elems: list[str], tol: float = 0.45) -> np.ndarray:
    """Bond adjacency perceived from interatomic distances and covalent radii."""
    r = np.array([_COV_R.get(e.upper(), 0.77) for e in elems])
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    cut = r[:, None] + r[None, :] + tol
    adj = (d < cut) & (d > 0.1)
    return adj


def best_ligand_mapping(model_smiles: str, model: Complex, ref: Complex,
                        align_resnums: list[int] | None = None) -> np.ndarray | None:
    """Map model ligand atom i to reference ligand atom `perm[i]`, symmetry-aware.

    **This must not assume the two structures list atoms in the same order.** They do
    not: a Boltz mmCIF for ligand 1RD begins C,C,C,C,N,C while the deposited crystal
    begins O,C,O,C,C,S. An index-for-index comparison of those gives 11.4 A where the
    true answer is far smaller, and it does so silently, poisoning every LDDT-PLI and
    RMSD downstream. An earlier version of this function made exactly that assumption.

    Method: perceive each ligand's bond graph from its own coordinates, find graph
    isomorphisms with element labels (VF2), and among them take the mapping that
    minimises ligand RMSD after superposing on the binding site. Isomorphism enumeration
    is what makes this symmetry-aware for free — a symmetric substituent yields several
    valid mappings and we take the best, rather than penalising a correct pose for an
    arbitrary atom labelling.

    Falls back to an element-constrained Hungarian assignment when no isomorphism is
    found (a co-folder occasionally emits distorted geometry that breaks bond
    perception). That fallback is a LOWER bound on RMSD, not the true value, so it is
    reported separately by callers that care.

    Returns None when the two ligands differ in heavy-atom composition, which is itself
    informative: the engine did not build the molecule we asked for.
    """
    if len(model.lig_xyz) == 0 or len(ref.lig_xyz) == 0:
        return None
    if len(model.lig_xyz) != len(ref.lig_xyz):
        return None
    me = [e.upper() for e in model.lig_elem]
    re_ = [e.upper() for e in ref.lig_elem]
    if sorted(me) != sorted(re_):
        return None

    # superpose on the binding site so that "best mapping" is measured in a shared frame
    resn = align_resnums or pocket_residues_from_structure(ref, radius=8.0)
    try:
        aligned, _fit, _n = align_by_residue(model, ref, resnums=resn)
        mxyz = aligned.lig_xyz
    except ValueError:
        mxyz = model.lig_xyz

    best, best_rmsd = None, np.inf
    try:
        import networkx as nx
        from networkx.algorithms.isomorphism import GraphMatcher, categorical_node_match

        gm_a = nx.from_numpy_array(_adjacency(model.lig_xyz, me))
        gm_b = nx.from_numpy_array(_adjacency(ref.lig_xyz, re_))
        nx.set_node_attributes(gm_a, {i: e for i, e in enumerate(me)}, "el")
        nx.set_node_attributes(gm_b, {i: e for i, e in enumerate(re_)}, "el")
        matcher = GraphMatcher(gm_a, gm_b, node_match=categorical_node_match("el", ""))
        for k, iso in enumerate(matcher.isomorphisms_iter()):
            if k > 20000:      # pathological symmetry; the best so far is good enough
                break
            p = np.array([iso[i] for i in range(len(me))])
            r = float(np.sqrt(((mxyz - ref.lig_xyz[p]) ** 2).sum(1).mean()))
            if r < best_rmsd:
                best, best_rmsd = p, r
    except Exception:
        best = None

    if best is not None:
        return best

    # fallback: element-constrained optimal assignment (a lower bound on true RMSD)
    try:
        from scipy.optimize import linear_sum_assignment

        d = np.linalg.norm(mxyz[:, None, :] - ref.lig_xyz[None, :, :], axis=2)
        cost = d.copy()
        cost[np.array(me)[:, None] != np.array(re_)[None, :]] = d.max() + 1e3
        _r, c = linear_sum_assignment(cost)
        return np.asarray(c)
    except Exception:
        return None


def bisy_rmsd(model: Complex, ref: Complex, align_resnums: list[int] | None = None,
              lig_perm: np.ndarray | None = None) -> float:
    """Binding-site-superposed, symmetry-corrected ligand RMSD.

    Superpose on the *binding site* (not the whole protein — a global fit lets a
    well-placed ligand in a shifted domain look bad and vice versa), then measure
    ligand RMSD under the best symmetry mapping.
    """
    if len(model.lig_xyz) == 0 or len(ref.lig_xyz) == 0:
        return float("nan")
    resn = align_resnums or pocket_residues_from_structure(ref, radius=8.0)
    try:
        aligned, _fit, _n = align_by_residue(model, ref, resnums=resn)
    except ValueError:
        return float("nan")
    if lig_perm is None:
        n = min(len(aligned.lig_xyz), len(ref.lig_xyz))
        lig_perm = np.arange(n)
    lm = aligned.lig_xyz[: len(lig_perm)]
    lr = ref.lig_xyz[lig_perm]
    return float(np.sqrt(((lm - lr) ** 2).sum(1).mean()))
