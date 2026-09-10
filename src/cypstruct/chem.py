"""Ligand chemistry: standardisation, and identifying which atom talks to the iron.

The important function here is `coordinating_atoms`. It ranks a molecule's candidate
heme-coordinating nitrogens, and that ranking is what lets us *steer* a co-folder.

Why this matters more than it looks. OpenADMET's own analysis of co-folding on ADMET
targets reports that on CYP3A4 the heme is placed correctly and the backbone is
accurate — the failure is **ligand orientation**, including near-180-degree flips, with
Boltz-2 reaching only 57.8% of poses under 2 A BiSyRMSD. Separately, our own parse of
116 deposited CYP3A4 entries finds **104 of 143 ligand-bound chains coordinate the iron
directly, through nitrogen in 103 of 104 cases**.

Those two facts compose into the central lever of this project: for the majority binding
mode, *which atom points at the iron completely determines the orientation*. So naming
the right donor atom turns an orientation problem into a constraint we can hand the
model, and turns pose selection into a check we can actually run.
"""
from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

# Ring systems whose sp2 N is a competent Fe(III) donor, best first. The ordering
# reflects the classical type II ligand series; it is a PRIOR used to break ties and
# is superseded by the computed proton affinity once tier-1 QM is available.
_DONOR_PATTERNS: list[tuple[str, str, float]] = [
    # (name, SMARTS matching the donor atom, prior donor strength 0-1)
    ("imidazole_N3",   "[nX2;r5]1[cX3][nX3][cX3][cX3]1", 1.00),
    ("triazole_N",     "[nX2;r5]1[nX2,nX3][nX2,nX3][cX3][cX3]1", 0.90),
    ("imidazole_generic", "[nX2;r5]", 0.85),
    ("pyridine_N",     "[nX2;r6]", 0.75),
    ("pyrimidine_N",   "[nX2;r6][cX3][nX2;r6]", 0.65),
    ("thiazole_N",     "[nX2;r5][cX3][sX2]", 0.60),
    ("oxazole_N",      "[nX2;r5][cX3][oX2]", 0.45),
    ("aliphatic_amine", "[NX3;H0,H1,H2;!$(N[C,S]=[O,S,N]);!$(N=*)]", 0.30),
    ("nitrile_N",      "[NX1]#[CX2]", 0.25),
]


@dataclass
class DonorSite:
    atom_idx: int
    element: str
    pattern: str
    prior: float
    n_neighbors: int
    steric_neighbors: int   # heavy atoms within 2 bonds of the donor (ortho crowding)
    aromatic: bool
    score: float            # prior discounted by steric crowding

    def __repr__(self) -> str:
        return (f"DonorSite(idx={self.atom_idx}, {self.element}, {self.pattern}, "
                f"score={self.score:.2f})")


_normalizer = rdMolStandardize.Normalizer()
_uncharger = rdMolStandardize.Uncharger()
_chooser = rdMolStandardize.LargestFragmentChooser()


def standardize(smiles: str) -> str | None:
    """Canonical, largest-fragment, neutralised SMILES. None if unparseable."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    try:
        m = _chooser.choose(_normalizer.normalize(m))
        m = _uncharger.uncharge(m)
        Chem.SanitizeMol(m)
    except Exception:
        return None
    return Chem.MolToSmiles(m)


def inchikey(smiles: str) -> str | None:
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(m) if m else None


def coordinating_atoms(smiles_or_mol, top_k: int = 3) -> list[DonorSite]:
    """Rank candidate heme-Fe coordinating atoms, best first.

    Two filters do most of the work, and both encode real coordination chemistry:

    - **Two-coordinate aromatic nitrogen only.** A pyrrole-type NH (three connections,
      lone pair in the pi system) has no available lone pair and cannot coordinate.
      Counting every aromatic N would nominate exactly the wrong atom on indoles,
      benzimidazoles and purines.
    - **Ortho crowding discount.** A 2,6-disubstituted pyridine is a poor ligand for
      steric reasons no electronic descriptor captures. Neighbours within two bonds
      of the donor discount the prior.

    Returns [] when nothing plausible is present — which is itself informative: that
    molecule is a type I substrate candidate and should be steered by site-of-metabolism
    geometry instead of by coordination.
    """
    mol = smiles_or_mol
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
    if mol is None:
        return []

    seen: dict[int, DonorSite] = {}
    for name, smarts, prior in _DONOR_PATTERNS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            for idx in match:
                a = mol.GetAtomWithIdx(idx)
                if a.GetSymbol() != "N":
                    continue
                # available lone pair: aromatic N must be two-coordinate;
                # aliphatic N must not be an amide and must have a free pair
                heavy_nbrs = [n for n in a.GetNeighbors()]
                if a.GetIsAromatic() and len(heavy_nbrs) != 2:
                    continue
                if a.GetFormalCharge() > 0:
                    continue          # protonated / quaternary N cannot donate
                if idx in seen and seen[idx].prior >= prior:
                    continue
                # ortho crowding: heavy atoms two bonds out, excluding the ring path
                two_bond = set()
                for n1 in heavy_nbrs:
                    for n2 in n1.GetNeighbors():
                        if n2.GetIdx() != idx:
                            two_bond.add(n2.GetIdx())
                crowd = sum(1 for j in two_bond
                            if mol.GetAtomWithIdx(j).GetDegree() > 2)
                score = prior * (0.75 ** max(0, crowd - 1))
                seen[idx] = DonorSite(atom_idx=idx, element="N", pattern=name,
                                      prior=prior, n_neighbors=len(heavy_nbrs),
                                      steric_neighbors=crowd,
                                      aromatic=a.GetIsAromatic(), score=score)
    return sorted(seen.values(), key=lambda d: -d.score)[:top_k]


def boltz_atom_name(smiles: str, atom_idx: int) -> str | None:
    """The name Boltz gives heavy atom `atom_idx` of a SMILES ligand.

    The rule is `SYMBOL + str(CanonicalRankAtoms(mol)[i] + 1)`, and the load-bearing
    detail is that the rank is computed **after hydrogens are added**, over all atoms.
    So the name is not a per-element counter and it is not the heavy-atom rank either:
    for the ritonavir analog 1RD the coordinating nitrogen is **N41**, where the
    heavy-atom-only rank would say N19 and a naive per-element counter would say N1.

    This was settled by running Boltz's own `parse_boltz_schema` on both candidates on
    CPU (`modal_boltz.probe_yaml`): N19 raised `KeyError: ('L', 0, 'N19')` and N41 parsed.
    Do not "simplify" this back to ranking the heavy-atom graph.

    A wrong name fails loudly at parse time rather than silently constraining the wrong
    atom, which is the one mercy here - but on the GPU path it costs a whole batch to
    find out, so `probe_yaml` is the cheap way to re-check after a Boltz upgrade.
    """
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or atom_idx >= mol.GetNumAtoms():
        return None
    # AddHs preserves heavy-atom indices, so atom_idx still addresses the same atom.
    with_h = Chem.AddHs(mol)
    rank = list(AllChem.CanonicalRankAtoms(with_h))
    name = with_h.GetAtomWithIdx(int(atom_idx)).GetSymbol().upper() + str(rank[atom_idx] + 1)
    return name if len(name) <= 4 else None



def ligand_class(smiles: str) -> str:
    """'type_II' if a competent Fe donor exists, else 'type_I_candidate'.

    A prediction, not a measurement — the deposited-structure classes in
    `data/processed/cyp3a4_reference_geometry.parquet` are the ground truth to score
    this against, and that check is worth running before trusting the steering.
    """
    d = coordinating_atoms(smiles, top_k=1)
    if d and d[0].score >= 0.55:
        return "type_II"
    if d:
        return "ambiguous"
    return "type_I_candidate"


def properties(smiles: str) -> dict:
    """Cheap descriptors used for stratifying splits and sanity-checking the set."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {}
    return {
        "mw": round(Descriptors.MolWt(m), 2),
        "clogp": round(Descriptors.MolLogP(m), 2),
        "tpsa": round(Descriptors.TPSA(m), 2),
        "n_heavy": m.GetNumHeavyAtoms(),
        "n_rot": rdMolDescriptors.CalcNumRotatableBonds(m),
        "n_rings": rdMolDescriptors.CalcNumRings(m),
        "n_arom_n": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N" and a.GetIsAromatic()),
        "formal_charge": Chem.GetFormalCharge(m),
    }


def aromatic_rings(mol) -> list[list[int]]:
    """Aromatic ring atom-index lists, for the porphyrin-stacking term."""
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
    if mol is None:
        return []
    out = []
    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            out.append(list(ring))
    return out


if __name__ == "__main__":
    # Known CYP3A4 ligands with known binding modes — a self-check of the classifier.
    tests = {
        "ketoconazole": "CC(=O)N1CCN(CC1)c1ccc(OC[C@@H]2CO[C@](Cn3ccnc3)(O2)c2ccc(Cl)cc2Cl)cc1",
        "ritonavir": "CC(C)c1nc(CN(C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](Cc2ccccc2)C[C@H](O)"
                     "[C@H](Cc2ccccc2)NC(=O)OCc2cncs2)cs1",
        "midazolam": "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
        "testosterone": "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O",
        "erythromycin_frag": "CC[C@H]1OC(=O)[C@H](C)[C@@H](O)[C@H](C)[C@@H](O)CC(C)(O)[C@@H](C)C1=O",
        "clotrimazole": "Clc1ccccc1C(n1ccnc1)(c1ccccc1)c1ccccc1",
    }
    for name, smi in tests.items():
        cls = ligand_class(smi)
        sites = coordinating_atoms(smi)
        print(f"{name:20s} {cls:18s} {sites}")
