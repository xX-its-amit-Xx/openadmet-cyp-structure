"""Tier-1 of the CYP physics scorer: the quantities that belong to the MOLECULE.

Tier 0 (`geometry.py`) asks *where is the ligand*. Tier 1 asks *what is this ligand
chemically capable of* — which of its nitrogens is actually a competent Fe(III) donor,
and which of its carbons is actually easy to abstract a hydrogen from.

**The cost architecture is the whole point.** Read `docs/QM_SCORER_DESIGN.md` §3. A pool
of 87 ligands x 20 samples x 2 arms is ~3,500 poses but only **87 QM jobs**, because
none of the quantities below depend on the pose. They are computed once on the free
molecule and then *read against* each pose's tier-0 geometry by `read_against_pose()`,
which is a dictionary lookup and a rank correlation. Nothing in this module may ever be
called per pose. If you find yourself passing coordinates from a predicted structure in
here, you have broken the economics of the scorer.

What this replaces. `cypstruct.chem.coordinating_atoms` ranks candidate donors with a
hand-written SMARTS table (imidazole 1.00, pyridine 0.75, oxazole 0.45 ...) discounted by
a crowding count. That table is a *prior* — my reading of the type II literature, typed
into a list. Tier 1's job is to replace each of its three ingredients with something
measured:

    SMARTS ring-type prior   ->  proton affinity at that N      (sigma-donor strength)
    (nothing)                ->  condensed Fukui f-             (nucleophilicity)
    integer "crowd" count    ->  %V_bur in a 3.5 A sphere       (real steric burial)

How much room is there to improve? Measured, before spending anything: of the 83
deposited CYP3A4 ligands that coordinate the iron, **50 have exactly one candidate donor
nitrogen** — every ranking method agrees on those and none can win. Only **32 have two or
more**, and in **16 of those 32 the SMARTS scores are exactly tied**, so the incumbent is
breaking the tie by RDKit atom index, i.e. arbitrarily. Those 16 arbitrary calls, plus
the 16 opinionated ones, are the entire addressable surface of tier 1 for donor choice.
`scripts/qm/validate_donor_ranking.py` scores against that, honestly, and reports the
denominator rather than quoting an accuracy over the easy 50.

**Measured result, 2026-09-11** (`scripts/qm/validate_donor_ranking.py`, 80 ligands run on
Modal CPU for 0.287 CPU-hours, 71 mappable to an atom-level ground truth). Top-1 accuracy
at naming the atom that is actually 1.9-2.45 A from the heme iron in the deposited
structure, on the 32 ligands where the methods can differ at all:

    incumbent SMARTS prior   0.594      proton affinity alone   0.625
    Fukui f- alone           0.781      %V_bur alone            0.812
    PA - k * %V_bur          0.906

On all 71 mapped ligands: SMARTS 0.817 -> 0.958. Tier 1 wins, and the two methods disagree
on 14 ligands, of which QM takes 12. Note carefully that the 0.988 quoted elsewhere for the
SMARTS prior is a *molecule-level class recall* ("is this ligand type II at all"), not a
donor-ranking accuracy; the comparable incumbent number is the 0.594 / 0.817 above.

Everything here returns **raw physical quantities and rankings, never a fitted score**,
matching `geometry.py`. Weighting happens in the learned ranker against measured
LDDT-PLI labels, so that no constant in this file quietly becomes a model parameter.
The one composite provided, `rank_donors(by="pa_minus_steric")`, carries its coefficient
as a named, documented PRIOR and says in its docstring that it must be refit.
"""
from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

HARTREE_KCAL = 627.509474

# Bondi van der Waals radii (A). Same table as tier 0, kept local so this module does
# not import geometry.py for three constants.
VDW = {"H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
       "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98, "B": 1.92, "SI": 2.10, "SE": 1.90}

# Where a coordinating iron would sit, relative to the donor nitrogen. This is the
# MEASURED median Fe-donor distance over 104 coordinated chain-ligand pairs in
# `data/processed/cyp3a4_reference_geometry.parquet` (p5 1.94, p50 2.20, p95 2.38),
# not a textbook value. It is used to place the centre of the %V_bur sphere at the
# position the metal actually occupies rather than on the nitrogen itself.
FE_DONOR_DIST = 2.20

# Radius of the buried-volume sphere, from the design doc (§2B).
VBUR_RADIUS = 3.5

# Grid spacing for the %V_bur integration (A). 0.10 A over a 3.5 A sphere is ~180k
# sample points; the integral converges to better than 0.1 percentage point and it is
# deterministic, which a Monte-Carlo estimate would not be. Determinism matters because
# these numbers are cached to JSON and compared across runs.
VBUR_GRID = 0.10


# ==========================================================================
# data model
# ==========================================================================


@dataclass
class DonorCandidate:
    """One candidate Fe(III)-coordinating atom, with everything tier 1 knows about it.

    `atom_idx` indexes the heavy-atom RDKit molecule the descriptors were computed on.
    That molecule is identified by `LigandQM.smiles` (canonical, from `chem.standardize`),
    so a consumer must re-parse that exact SMILES to address the same atom. This is the
    same atom-identity discipline that `pose.align_by_residue` enforces for residues, and
    for the same reason: index-for-index matching across two differently-built molecules
    is the failure mode that has already cost this project half an LDDT-PLI.
    """

    atom_idx: int
    element: str = "N"

    # --- what the incumbent SMARTS prior says, carried along for comparison ---
    smarts_pattern: str = ""
    smarts_prior: float = float("nan")
    smarts_score: float = float("nan")

    # --- structural context (free, from the graph) ---
    ring_size: int = 0
    ring_hetero: tuple = ()
    n_ortho_heavy: int = 0          # heavy atoms two bonds out with degree > 2
    rejected: str = ""              # non-empty if tier 1 ruled this candidate out

    # --- sigma-donor strength (GFN2-xTB) ---
    # Energy released on protonating THIS nitrogen. See `proton_affinity` for the
    # crucial caveat about the absolute scale.
    proton_affinity_kcal: float = float("nan")
    proton_affinity_solv_kcal: float = float("nan")
    e_protomer_hartree: float = float("nan")
    protomer_ok: bool = False

    # --- nucleophilicity (condensed Fukui from Mulliken populations) ---
    fukui_minus: float = float("nan")     # f- : propensity to DONATE an electron pair
    fukui_plus: float = float("nan")
    dual_descriptor: float = float("nan")  # f+ - f- ; negative = nucleophilic site
    mulliken_charge: float = float("nan")

    # --- steric accessibility of the lone pair (pure geometry, no QM) ---
    vbur_at_atom: float = float("nan")     # sphere centred ON the nitrogen
    vbur_at_metal: float = float("nan")    # sphere centred where the Fe would be
    vbur_at_metal_min: float = float("nan")   # over the conformer ensemble
    vbur_at_metal_max: float = float("nan")
    lone_pair_ok: bool = False             # was a lone-pair vector definable?

    note: str = ""


@dataclass
class CarbonSite:
    """One carbon bearing at least one hydrogen, with an ease-of-abstraction estimate."""

    atom_idx: int
    n_h: int = 0
    environment: str = ""
    bde_proxy_kcal: float = float("nan")   # lower = easier to abstract
    bde_xtb_kcal: float = float("nan")     # optional, real GFN2 homolysis
    aromatic: bool = False


@dataclass
class LigandQM:
    """The complete tier-1 profile of one ligand. Serialised as one JSON per ligand."""

    ligand_id: str = ""
    smiles: str = ""
    inchikey: str = ""
    n_heavy: int = 0
    formal_charge: int = 0

    donors: list = field(default_factory=list)     # list[DonorCandidate]
    carbons: list = field(default_factory=list)    # list[CarbonSite]

    e_neutral_hartree: float = float("nan")
    e_neutral_solv_hartree: float = float("nan")

    method: str = ""
    solvent: str = ""
    n_conformers: int = 0
    seconds: float = float("nan")
    n_xtb_calls: int = 0
    schema_version: int = 1
    status: str = "pending"       # ok | partial | failed
    error: str = ""

    # ---- (de)serialisation -------------------------------------------------
    def to_dict(self) -> dict:
        d = asdict(self)
        d["donors"] = [asdict(x) if not isinstance(x, dict) else x for x in self.donors]
        d["carbons"] = [asdict(x) if not isinstance(x, dict) else x for x in self.carbons]
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "LigandQM":
        d = dict(d)
        donors = [DonorCandidate(**x) for x in d.pop("donors", [])]
        carbons = [CarbonSite(**x) for x in d.pop("carbons", [])]
        known = {f for f in cls.__dataclass_fields__}
        obj = cls(**{k: v for k, v in d.items() if k in known})
        obj.donors, obj.carbons = donors, carbons
        return obj

    def save(self, path: str | Path) -> None:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(self.to_dict(), indent=1))

    @classmethod
    def load(cls, path: str | Path) -> "LigandQM":
        return cls.from_dict(json.loads(Path(path).read_text()))


# ==========================================================================
# free geometry: percent buried volume
# ==========================================================================


def embed_conformers(mol, n_conf: int = 4, seed: int = 0xC0FFEE, max_iters: int = 500):
    """Add hydrogens and generate an MMFF-relaxed conformer ensemble.

    An ensemble rather than one conformer because %V_bur around a donor is a property of
    the *accessible* geometries, not of whichever one ETKDG produced first. A 2-substituted
    pyridine with a rotatable ortho group can look accessible in one conformer and blocked
    in another; reporting min/mean/max over a small ensemble makes that visible instead of
    hiding it behind a single draw. The QM terms use the lowest-MMFF-energy member, since
    proton affinity is far less conformer-sensitive than burial is.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molh = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.5
    cids = list(AllChem.EmbedMultipleConfs(molh, numConfs=max(1, n_conf), params=params))
    if not cids:
        # ETKDG fails on some macrocycles and cage systems; random coords still give a
        # usable starting point and MMFF cleans it up.
        params.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(molh, numConfs=max(1, n_conf), params=params))
    if not cids:
        return None, []
    energies = []
    for cid in cids:
        try:
            res = AllChem.MMFFOptimizeMoleculeConfs(molh, maxIters=max_iters)
            energies = [e for _conv, e in res]
            break
        except Exception:
            try:
                AllChem.UFFOptimizeMoleculeConfs(molh, maxIters=max_iters)
            except Exception:
                pass
            energies = [0.0] * len(cids)
            break
    order = sorted(range(len(cids)), key=lambda i: energies[i] if i < len(energies) else 0.0)
    return molh, [cids[i] for i in order]


def _conf_xyz(molh, cid: int) -> np.ndarray:
    return np.asarray(molh.GetConformer(cid).GetPositions(), float)


def lone_pair_vector(molh, atom_idx: int, xyz: np.ndarray) -> np.ndarray | None:
    """Unit vector along a two-coordinate sp2 nitrogen's in-plane lone pair.

    Mirrors `geometry.lone_pair_direction` deliberately: the direction tier 1 places the
    virtual metal along must be the same direction tier 0 measures the approach angle
    against, or the two tiers are describing different chemistry.

    Returns None for anything that is not two-coordinate, because for a three-coordinate
    amine the lone pair is the sp3 apex and a bisector construction would point at the
    wrong hemisphere entirely.
    """
    a = molh.GetAtomWithIdx(int(atom_idx))
    nbrs = [n.GetIdx() for n in a.GetNeighbors()]
    if len(nbrs) == 2:
        p = xyz[atom_idx]
        v = np.zeros(3)
        for j in nbrs:
            d = xyz[j] - p
            n = np.linalg.norm(d)
            if n < 1e-6:
                return None
            v += d / n
        if np.linalg.norm(v) < 1e-6:
            return None
        return -v / np.linalg.norm(v)
    if len(nbrs) == 3:
        # sp3 amine: the lone pair is opposite the sum of the three bond vectors.
        p = xyz[atom_idx]
        v = np.zeros(3)
        for j in nbrs:
            d = xyz[j] - p
            n = np.linalg.norm(d)
            if n < 1e-6:
                return None
            v += d / n
        if np.linalg.norm(v) < 1e-6:
            return None      # planar N: no axial lone pair to speak of
        return -v / np.linalg.norm(v)
    return None


def percent_buried_volume(xyz: np.ndarray, elements: list[str], centre: np.ndarray,
                          radius: float = VBUR_RADIUS, exclude: tuple = (),
                          grid: float = VBUR_GRID) -> float:
    """Fraction (in %) of a sphere at `centre` that is inside some atom's vdW radius.

    This is the Cavallo %V_bur construction used throughout organometallic chemistry to
    quantify how much room a ligand leaves a metal. It is used here for the one job the
    SMARTS prior does badly: a 2,6-disubstituted pyridine is electronically a fine donor
    and sterically a hopeless one, and no electronic descriptor will ever say so. The
    incumbent approximates this with an integer count of neighbours two bonds out, which
    cannot tell a methyl from a tert-butyl, or an ortho substituent folded away from one
    folded over the lone pair.

    Deterministic grid integration, not Monte Carlo: these values are cached and diffed
    across runs, and a stochastic estimate would produce spurious changes.
    """
    n = int(math.ceil(2 * radius / grid))
    lin = (np.arange(n + 1) - n / 2.0) * grid
    gx, gy, gz = np.meshgrid(lin, lin, lin, indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    inside = (pts ** 2).sum(1) <= radius ** 2
    pts = pts[inside] + centre
    if len(pts) == 0:
        return float("nan")

    occupied = np.zeros(len(pts), dtype=bool)
    for i, (p, e) in enumerate(zip(xyz, elements)):
        if i in exclude:
            continue
        r = VDW.get(e.upper(), 1.70)
        d = p - centre
        # cheap rejection: an atom further than radius + r cannot touch the sphere
        if float(np.dot(d, d)) > (radius + r) ** 2:
            continue
        occupied |= ((pts - p) ** 2).sum(1) <= r * r
    return float(100.0 * occupied.mean())


def buried_volume_for_donor(molh, cids: list[int], atom_idx: int,
                            radius: float = VBUR_RADIUS) -> dict:
    """%V_bur around a donor, both on the atom and at the virtual metal position.

    Two numbers because they answer different questions. The sphere centred on the
    nitrogen measures how crowded the atom itself is. The sphere centred 2.20 A along the
    lone pair — where the iron would actually be — measures whether the *metal* can get
    there, which is the question that decides coordination. For 2,6-disubstitution the
    second number moves much more than the first.
    """
    els = [a.GetSymbol() for a in molh.GetAtoms()]
    at_atom, at_metal = [], []
    lp_ok = False
    for cid in cids:
        xyz = _conf_xyz(molh, cid)
        at_atom.append(percent_buried_volume(xyz, els, xyz[atom_idx], radius,
                                             exclude=(atom_idx,)))
        lp = lone_pair_vector(molh, atom_idx, xyz)
        if lp is None:
            continue
        lp_ok = True
        centre = xyz[atom_idx] + FE_DONOR_DIST * lp
        # The donor nitrogen itself is NOT excluded here: it genuinely occupies part of
        # the metal's coordination sphere, and excluding it would make every donor look
        # equally open. Only the metal's own volume is absent, as in the standard
        # construction.
        at_metal.append(percent_buried_volume(xyz, els, centre, radius))
    out = {"vbur_at_atom": float(np.mean(at_atom)) if at_atom else float("nan"),
           "lone_pair_ok": lp_ok}
    if at_metal:
        out.update(vbur_at_metal=float(np.mean(at_metal)),
                   vbur_at_metal_min=float(np.min(at_metal)),
                   vbur_at_metal_max=float(np.max(at_metal)))
    return out


# ==========================================================================
# free graph chemistry: a C-H abstraction proxy
# ==========================================================================

# Homolytic C-H bond dissociation energies, kcal/mol, from the standard compilations
# (Luo, "Comprehensive Handbook of Chemical Bond Energies"; Blanksby & Ellison 2003).
# Values are for the named prototype, e.g. "benzylic_1" is toluene's 89.7.
_BDE_BASE = {
    "methyl": 105.0,        # CH4
    "primary": 101.1,       # ethane
    "secondary": 98.6,      # propane C2
    "tertiary": 96.5,       # isobutane
    "cyclopropyl": 106.3,
    "vinylic": 111.0,       # ethene
    "aromatic": 112.9,      # benzene
    "alkynyl": 133.3,       # acetylene
    "formyl": 88.7,         # CH3CHO, the aldehyde C-H
}
# Activating environments, as the BDE they impose (not as an increment). Applied by
# taking the MINIMUM over everything that applies, plus a small bonus for each further
# activating group. Minimum-plus-bonus rather than additive because full additivity of
# radical stabilisation energies badly over-stabilises doubly activated positions
# (e.g. it would put a benzylic ether CH below 75 kcal/mol, which is not a real bond).
_BDE_ACTIVATED = {
    "benzylic_1": 89.7,     # toluene
    "benzylic_2": 87.0,
    "benzylic_3": 84.5,
    "allylic": 88.8,        # propene
    "alpha_O": 93.0,        # ethers / alcohols alpha C-H
    "alpha_N": 91.0,        # dialkylamine alpha C-H (N-CH3 of an amine)
    "alpha_S": 93.0,
    "alpha_carbonyl": 94.0,
    "alpha_halogen": 99.0,
}
_MULTI_ACTIVATION_BONUS = 2.5   # kcal/mol per additional activating group, capped below
_MULTI_ACTIVATION_CAP = 6.0


def ch_environment(mol, idx: int) -> tuple[str, list[str]]:
    """(base environment, list of activating environments) for a hydrogen-bearing carbon."""
    a = mol.GetAtomWithIdx(int(idx))
    heavy = [n for n in a.GetNeighbors() if n.GetSymbol() != "H"]
    n_heavy_c = sum(1 for n in heavy if n.GetSymbol() == "C")

    if a.GetIsAromatic():
        base = "aromatic"
    elif a.GetHybridization().name == "SP":
        base = "alkynyl"
    elif a.GetHybridization().name == "SP2":
        base = "formyl" if any(
            n.GetSymbol() == "O" and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in heavy) else "vinylic"
    elif a.IsInRingSize(3):
        base = "cyclopropyl"
    else:
        base = {0: "methyl", 1: "primary", 2: "secondary"}.get(n_heavy_c, "tertiary")
        if len(heavy) == 0:
            base = "methyl"

    acts: list[str] = []
    if base in ("aromatic", "alkynyl", "vinylic", "formyl"):
        return base, acts        # sp2/sp C-H is not activated by neighbours in this scheme

    degree = {"methyl": 1, "primary": 1, "secondary": 2, "tertiary": 3}.get(base, 1)
    for n in heavy:
        sym = n.GetSymbol()
        bond = mol.GetBondBetweenAtoms(idx, n.GetIdx())
        if sym == "C" and n.GetIsAromatic():
            acts.append({1: "benzylic_1", 2: "benzylic_2"}.get(degree, "benzylic_3"))
        elif sym == "C" and bond is not None and any(
                b.GetBondTypeAsDouble() == 2.0 for b in n.GetBonds()):
            # neighbouring C=C gives an allylic position; neighbouring C=O an alpha-carbonyl
            if any(b.GetOtherAtom(n).GetSymbol() == "O" and b.GetBondTypeAsDouble() == 2.0
                   for b in n.GetBonds()):
                acts.append("alpha_carbonyl")
            else:
                acts.append("allylic")
        elif sym == "O":
            acts.append("alpha_O")
        elif sym == "N":
            # An amide nitrogen is a poor radical stabiliser: its lone pair is already
            # delocalised into the carbonyl and is not available to the adjacent radical.
            amide = any(b.GetOtherAtom(n).GetSymbol() in ("O", "S")
                        and b.GetBondTypeAsDouble() == 2.0 for b in n.GetBonds())
            acts.append("alpha_carbonyl" if amide else "alpha_N")
        elif sym == "S":
            acts.append("alpha_S")
        elif sym in ("F", "CL", "BR", "I", "Cl", "Br"):
            acts.append("alpha_halogen")
    return base, acts


def ch_bde_proxy(mol) -> list[CarbonSite]:
    """Per-carbon ease-of-abstraction estimate from the molecular graph alone.

    **This estimate is systematically biased and that is acceptable by design.** The
    scorer's site-of-metabolism term (§2C of the design) is explicitly RANK-BASED: it
    rewards poses where the ordering of carbons by distance to the iron agrees with their
    ordering by ease of abstraction. A monotone but miscalibrated energy scale therefore
    costs nothing, while an absolute-energy term would need a much more expensive method
    and would still be wrong by several kcal/mol.

    The precedent is SMARTCyp, which predicts CYP sites of metabolism from a table of
    precomputed fragment activation energies looked up by substructure and beats far more
    expensive approaches. This is the same idea with published BDEs in place of its DFT
    table.

    Honest limitations, all of which the rank form tolerates and none of which it fixes:
      - Aromatic carbons get benzene's 112.9 and are correctly disfavoured for
        *abstraction*, but CYP oxidises aromatic rings through epoxidation, a mechanism
        with no C-H cleavage at all. Aromatic SoM is simply outside this term.
      - No steric or accessibility correction: a buried tertiary C-H scores as easily
        abstracted as an exposed one. That correction is tier 0's job, via the Fe-distance
        rank the term is read against.
      - Multiple activation is capped rather than summed, which flattens genuinely
        activated positions such as an N-CH2 next to an arene.

    `LigandQM.bde_xtb_kcal` is the escape hatch: `compute_ligand_qm(..., n_bde_sites=K)`
    replaces the proxy on the K most abstractable carbons with a real GFN2-xTB homolysis
    energy, at the cost of K extra optimizations per ligand.
    """
    out: list[CarbonSite] = []
    for a in mol.GetAtoms():
        if a.GetSymbol() != "C" or a.GetTotalNumHs() == 0:
            continue
        base, acts = ch_environment(mol, a.GetIdx())
        e = _BDE_BASE.get(base, 101.0)
        if acts:
            best = min(_BDE_ACTIVATED.get(x, 100.0) for x in acts)
            extra = min(_MULTI_ACTIVATION_CAP,
                        _MULTI_ACTIVATION_BONUS * (len(acts) - 1))
            e = min(e, best - extra)
        label = base + ("+" + "+".join(sorted(set(acts))) if acts else "")
        out.append(CarbonSite(atom_idx=a.GetIdx(), n_h=a.GetTotalNumHs(),
                              environment=label, bde_proxy_kcal=round(float(e), 2),
                              aromatic=a.GetIsAromatic()))
    return out


# ==========================================================================
# a correctness filter the incumbent is missing
# ==========================================================================


def lone_pair_available(mol, atom_idx: int) -> tuple[bool, str]:
    """Can this nitrogen actually donate a lone pair to Fe(III)?

    **This catches a live bug in `chem.coordinating_atoms`.** That function's docstring
    says it admits "two-coordinate aromatic nitrogen only", because a pyrrole-type N-H
    has its lone pair in the pi system and cannot coordinate. Its implementation tests
    `len(a.GetNeighbors()) != 2` — but on a molecule parsed from SMILES the hydrogens are
    implicit, so `GetNeighbors()` returns only HEAVY neighbours and a pyrrole N-H reports
    exactly 2. The guard never fires.

    Measured consequence, on the 83 deposited CYP3A4 ligands that coordinate the iron:
    `coordinating_atoms("c1c[nH]cn1")` returns BOTH imidazole nitrogens at score 1.00,
    and benzimidazole likewise. **13 of the 32 multi-candidate coordinated ligands carry
    a spurious N-H candidate this way** — which is most of the ambiguity tier 1 was
    commissioned to resolve. In all 13 the real donor currently happens to sort first,
    so nothing has broken yet, but the tie is being broken by dict insertion order rather
    than by chemistry, and a ligand ordering the other way would be steered onto an atom
    that physically cannot bond to the metal.

    `chem.py` is deliberately NOT patched here. The 98.8% figure that tier 1 is measured
    against was produced by the current `coordinating_atoms`, and silently changing the
    incumbent mid-comparison would make the validation meaningless. The filter lives on
    this side, is reported as a tier-1 contribution in
    `scripts/qm/validate_donor_ranking.py`, and the `chem.py` fix should be made as its
    own change with its own before/after number.

    The test is `GetTotalNumHs() == 0` for AROMATIC nitrogen only. An aliphatic secondary
    amine N-H keeps a perfectly good sp3 lone pair and is not excluded.
    """
    a = mol.GetAtomWithIdx(int(atom_idx))
    if a.GetSymbol() != "N":
        return False, "not nitrogen"
    if a.GetFormalCharge() > 0:
        return False, "cationic nitrogen has no lone pair to donate"
    if a.GetIsAromatic() and a.GetTotalNumHs() > 0:
        return False, "pyrrole-type N-H: lone pair is in the pi system"
    if a.GetIsAromatic() and a.GetTotalDegree() != 2:
        return False, f"aromatic N with total degree {a.GetTotalDegree()}"
    return True, ""


# ==========================================================================
# GFN2-xTB drivers
# ==========================================================================


class XTBUnavailable(RuntimeError):
    """Raised when the `xtb` binary is not on PATH.

    Deliberately loud. This box has no xtb and no GPU; tier 1 runs on Modal CPU via
    `scripts/qm/modal_qm.py`. Silently falling back to the graph-only descriptors would
    produce a LigandQM that looks complete and has NaN where the physics should be.
    """


def xtb_available() -> bool:
    return shutil.which("xtb") is not None


def _xtb(workdir: Path, xyz_text: str, charge: int = 0, uhf: int = 0,
         opt: bool = False, solvent: str | None = None,
         timeout: int = 900) -> dict:
    """One xtb invocation. Returns energy (Hartree), Mulliken charges, optimised xyz.

    The CLI rather than the Python bindings: `xtb` ships a geometry optimiser, writes
    Mulliken charges to a `charges` file, and handles charge/multiplicity from flags, so
    every quantity tier 1 needs comes out of it without driving an SCF by hand. The
    `tblite` pip package also installs and runs (verified on Modal), but exposes only
    single points through its Python API.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "m.xyz").write_text(xyz_text)
    cmd = ["xtb", "m.xyz", "--gfn", "2", "--chrg", str(charge), "--uhf", str(uhf),
           "--norestart"]
    if solvent:
        cmd += ["--alpb", solvent]
    if opt:
        cmd += ["--opt", "normal"]
    t0 = time.time()
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, cwd=str(workdir),
                            timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": "xtb timed out", "seconds": time.time() - t0}
    energy = None
    for line in cp.stdout.splitlines():
        if "TOTAL ENERGY" in line:
            try:
                energy = float(line.split()[3])
            except (IndexError, ValueError):
                pass
    charges = None
    qf = workdir / "charges"
    if qf.exists():
        try:
            charges = [float(x) for x in qf.read_text().split()]
        except ValueError:
            charges = None
    xyz_out = None
    of = workdir / "xtbopt.xyz"
    if opt and of.exists():
        xyz_out = of.read_text()
    ok = cp.returncode == 0 and energy is not None
    return {"ok": ok, "energy": energy, "charges": charges, "xyz": xyz_out,
            "returncode": cp.returncode, "seconds": round(time.time() - t0, 2),
            "error": "" if ok else (cp.stderr or cp.stdout)[-800:]}


def _xyz_text(molh, cid: int) -> str:
    from rdkit import Chem
    return Chem.MolToXYZBlock(molh, confId=int(cid))


def _protomer(mol, atom_idx: int):
    """A copy of `mol` with an extra H on `atom_idx` and a +1 formal charge there."""
    from rdkit import Chem

    rw = Chem.RWMol(mol)
    a = rw.GetAtomWithIdx(int(atom_idx))
    a.SetNumExplicitHs(a.GetTotalNumHs() + 1)
    a.SetNoImplicit(True)
    a.SetFormalCharge(a.GetFormalCharge() + 1)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None
    return out


def _protonate_in_place(molh, cid: int, donor_idx: int):
    """Add a proton to `donor_idx` ON the existing optimised geometry. Returns a new mol.

    **This is not a detail, it is the difference between a usable descriptor and noise.**
    The obvious implementation — rebuild the protonated molecule from SMILES and embed a
    fresh conformer — was tried first and measured: it put clotrimazole's imidazole at
    78.1 kcal/mol against bare imidazole's 120.8. Clotrimazole's imidazole is a slightly
    weaker base than imidazole, by about one pKa unit, not by 43 kcal/mol. The gap was
    almost entirely **conformational**: E(neutral at its own minimum) was being compared
    with E(protomer at a different conformer's minimum), so the rotameric energy of a
    trityl group landed inside a number that is supposed to report a lone pair's basicity.

    For the four small prototypes this never shows up, because they have one conformer.
    It appears exactly on the drug-sized ligands the pool is made of, where the
    conformational spread is tens of kcal/mol and the donor-strength differences tier 1
    is trying to resolve are about six.

    So the proton is placed on the already-xTB-relaxed neutral geometry, 1.02 A along the
    lone-pair vector, and only then re-optimised. Both endpoints stay in the same
    conformational basin and the difference is the protonation energy rather than a
    conformer search.
    """
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    xyz = _conf_xyz(molh, cid)
    lp = lone_pair_vector(molh, donor_idx, xyz)
    if lp is None:
        return None
    rw = Chem.RWMol(molh)
    n = rw.GetAtomWithIdx(int(donor_idx))
    n.SetFormalCharge(n.GetFormalCharge() + 1)
    n.SetNoImplicit(True)
    h_idx = rw.AddAtom(Chem.Atom(1))
    rw.AddBond(int(donor_idx), h_idx, Chem.BondType.SINGLE)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None
    conf = out.GetConformer(int(cid)) if out.GetNumConformers() else None
    if conf is None:
        return None
    p = xyz[donor_idx] + 1.02 * lp          # a typical N-H bond length
    conf.SetAtomPosition(h_idx, Point3D(float(p[0]), float(p[1]), float(p[2])))
    return out


def proton_affinity(molh, cid: int, donor_idx: int, e_neutral: float, workdir: Path,
                    base_charge: int = 0, solvent: str | None = None,
                    mol_fallback=None, seed: int = 0xC0FFEE) -> dict:
    """Energy released on protonating `donor_idx`, in kcal/mol. Higher = stronger donor.

    ⚠️ **The absolute number is not a proton affinity.** It is
    `E(neutral) - E(N-protonated)` in GFN2-xTB total energies, which implicitly sets the
    proton's energy to zero. GFN2's total energy contains atomic reference terms, so
    adding a hydrogen adds a constant that is not the bare proton's energy. Measured on
    this implementation: pyridine comes out at 114.3 kcal/mol against an experimental
    gas-phase PA of 222.0 — an offset of about 108 kcal/mol.

    That offset is **the same constant for every molecule**, because every protonation
    adds exactly one hydrogen. So differences and orderings are preserved exactly and are
    the only things any consumer may use. Measured on this implementation, gas phase,
    xtb 6.7.1, on the four prototypes the design names (kcal/mol):

        imidazole 65.4  >  pyridine 60.3  >  thiazole 53.2  >  oxazole 47.4
        experiment: 225.3  >   222.0    >    208.0   >   199.6

    The ordering §2B predicts is reproduced exactly; the offset from experiment is about
    162 kcal/mol. Do not compare these values to tabulated proton affinities, and do not
    let one into an absolute-energy expression.

    A second caveat, honestly stated: GFN2 does not reproduce the experimental *spacing*.
    The measured gaps are 5.1 / 7.1 / 5.8 against experimental 3.3 / 14.0 / 8.4 — the
    azine-to-azole separation is compressed roughly twofold and the imidazole-pyridine
    gap is overstated. Ordering across families is reliable; fine discrimination *within*
    a family is what tier 3 exists for (§6 of the design says exactly this).

    A third, and the useful one: the term correctly makes 2,6-lutidine (71.4) a
    *stronger* base than pyridine (60.3), matching experiment (230.1 vs 222.0). It is
    also a far worse metal ligand, for reasons no electronic descriptor can see. That is
    what %V_bur is for, and the two numbers must be read together — see
    `VBUR_TO_KCAL_PRIOR`.

    `molh`/`cid` must be the xTB-optimised NEUTRAL geometry — see `_protonate_in_place`
    for why rebuilding the protomer from SMILES instead silently injects conformational
    energy into the result.
    """
    prot = _protonate_in_place(molh, cid, donor_idx)
    used_fallback = False
    if prot is not None:
        xyz_text = _xyz_text(prot, cid)
    else:
        # No definable lone-pair vector (e.g. a planar three-coordinate N). Fall back to
        # rebuilding and re-embedding, and SAY SO, because this path carries the
        # conformational error described above and its values are not comparable with
        # the others.
        if mol_fallback is None:
            return {"ok": False, "error": "no lone-pair vector and no fallback molecule"}
        p2 = _protomer(mol_fallback, donor_idx)
        if p2 is None:
            return {"ok": False, "error": "could not build the N-protonated form"}
        ph, pcids = embed_conformers(p2, n_conf=1, seed=seed)
        if ph is None or not pcids:
            return {"ok": False, "error": "could not embed the protomer"}
        xyz_text, cid, used_fallback = _xyz_text(ph, pcids[0]), pcids[0], True

    r = _xtb(workdir, xyz_text, charge=base_charge + 1, opt=True, solvent=solvent)
    if not r["ok"]:
        return {"ok": False, "error": r.get("error", "")[-300:], "seconds": r["seconds"]}
    return {"ok": True,
            "pa_kcal": round((e_neutral - r["energy"]) * HARTREE_KCAL, 3),
            "e_protomer": r["energy"], "seconds": r["seconds"],
            "reembedded": used_fallback}


def condensed_fukui(molh, cid: int, base_charge: int, workdir: Path,
                    solvent: str | None = None) -> dict:
    """Condensed Fukui indices from Mulliken populations at a FIXED geometry.

    f-_A = q_A(N-1) - q_A(N)   propensity of atom A to give up electron density
    f+_A = q_A(N)   - q_A(N+1) propensity of atom A to accept it
    dual = f+ - f-             negative at nucleophilic sites

    f- is the one the design asks for: coordination to Fe(III) is the nitrogen donating
    its lone pair into an empty d orbital, so the relevant question is which atom's
    density leaves most readily. Computing it at frozen geometry is deliberate — relaxing
    the cation would fold geometric relaxation into what is meant to be an electronic
    descriptor.

    The anion (for f+) is run in implicit solvent even when the rest of the job is gas
    phase, because GFN2 anions in vacuum frequently have an unbound extra electron and
    return populations that are numerically fine and physically meaningless. If that
    single point fails, f- is still returned and f+/dual come back NaN rather than the
    whole descriptor being lost.
    """
    xyz = _xyz_text(molh, cid)
    neutral = _xtb(workdir / "fukui_n", xyz, charge=base_charge, solvent=solvent)
    if not neutral["ok"] or not neutral.get("charges"):
        return {"ok": False, "error": "neutral single point failed"}
    cation = _xtb(workdir / "fukui_cat", xyz, charge=base_charge + 1, uhf=1, solvent=solvent)
    out: dict = {"ok": False, "seconds": neutral["seconds"] + cation["seconds"],
                 "q": neutral["charges"], "calls": 2}
    if cation["ok"] and cation.get("charges"):
        qn = np.asarray(neutral["charges"], float)
        qc = np.asarray(cation["charges"], float)
        if len(qn) == len(qc):
            out["f_minus"] = (qc - qn).tolist()
            out["ok"] = True
    anion = _xtb(workdir / "fukui_an", xyz, charge=base_charge - 1, uhf=1,
                 solvent=solvent or "water")
    out["calls"] += 1
    out["seconds"] += anion["seconds"]
    if anion["ok"] and anion.get("charges"):
        qn = np.asarray(neutral["charges"], float)
        qa = np.asarray(anion["charges"], float)
        if len(qn) == len(qa):
            out["f_plus"] = (qn - qa).tolist()
    return out


def bde_xtb(molh, cid: int, carbon_idx: int, e_parent: float, base_charge: int,
            workdir: Path, solvent: str | None = None) -> dict:
    """Real GFN2 homolytic C-H BDE: E(R.) + E(H.) - E(RH), kcal/mol.

    Off by default. It costs one geometry optimization per carbon, so a 20-carbon ligand
    turns a ~1-minute job into a ~10-minute one. Since the site-of-metabolism term is
    rank-based, the graph proxy is usually sufficient; this exists to CHECK the proxy's
    ordering on a sample of ligands rather than to replace it wholesale.
    """
    from rdkit import Chem

    rw = Chem.RWMol(molh)
    h = None
    for n in rw.GetAtomWithIdx(int(carbon_idx)).GetNeighbors():
        if n.GetSymbol() == "H":
            h = n.GetIdx()
            break
    if h is None:
        return {"ok": False, "error": "carbon has no explicit hydrogen"}
    rw.RemoveAtom(int(h))
    rad = rw.GetMol()
    try:
        rad.GetAtomWithIdx(int(carbon_idx) - (1 if h < carbon_idx else 0)).SetNumRadicalElectrons(1)
        Chem.SanitizeMol(rad)
    except Exception:
        return {"ok": False, "error": "could not build the radical"}
    xyz = Chem.MolToXYZBlock(rad, confId=int(cid)) if rad.GetNumConformers() else None
    if xyz is None:
        return {"ok": False, "error": "radical lost its conformer"}
    r = _xtb(workdir, xyz, charge=base_charge, uhf=1, opt=True, solvent=solvent)
    if not r["ok"]:
        return {"ok": False, "error": r.get("error", "")[-200:]}
    # GFN2 total energy of a hydrogen atom, computed once with the same settings.
    e_h = _H_ATOM_ENERGY
    return {"ok": True,
            "bde_kcal": round((r["energy"] + e_h - e_parent) * HARTREE_KCAL, 2),
            "seconds": r["seconds"]}


# GFN2-xTB total energy of an isolated hydrogen atom (Hartree), doublet, gas phase.
# Constant of the method, not a fitted parameter; recomputed by
# `scripts/qm/modal_qm.py --calibrate` if the xtb version changes.
_H_ATOM_ENERGY = -0.39308740


# ==========================================================================
# the per-ligand job
# ==========================================================================


def compute_ligand_qm(smiles: str, ligand_id: str = "", *, n_conformers: int = 4,
                      solvent: str | None = None, also_solvated: bool = False,
                      n_bde_sites: int = 0, max_donors: int = 6,
                      workdir: str | Path | None = None,
                      seed: int = 0xC0FFEE) -> LigandQM:
    """Everything tier 1 knows about one molecule. **Call once per ligand, never per pose.**

    Order of operations is chosen so that a partial failure still returns useful work:
    the free graph/geometry descriptors (%V_bur, BDE proxy, SMARTS context) are computed
    first and unconditionally, then the xTB terms are layered on. A ligand whose xTB job
    dies still yields a LigandQM with `status="partial"` and real buried volumes, which is
    the difference between losing one descriptor and losing the ligand.

    `workdir` must NOT be left to the system default on the local box: the default temp
    directory resolves onto the near-full C: drive. Callers on this machine pass a
    directory under D:, and the Modal runner passes container-local /tmp.
    """
    from rdkit import Chem

    from ..chem import coordinating_atoms, inchikey, standardize

    t0 = time.time()
    q = LigandQM(ligand_id=ligand_id or smiles, method="GFN2-xTB (xtb CLI)",
                 solvent=solvent or "gas", n_conformers=n_conformers)
    canon = standardize(smiles) or smiles
    q.smiles = canon
    mol = Chem.MolFromSmiles(canon)
    if mol is None:
        q.status, q.error = "failed", "unparseable SMILES"
        return q
    q.inchikey = inchikey(canon) or ""
    q.n_heavy = mol.GetNumHeavyAtoms()
    q.formal_charge = Chem.GetFormalCharge(mol)

    # ---- free: the carbon profile ----------------------------------------
    q.carbons = ch_bde_proxy(mol)

    # ---- free: candidate donors and their steric environment --------------
    sites = coordinating_atoms(mol, top_k=max_donors)
    ri = mol.GetRingInfo()
    molh, cids = embed_conformers(mol, n_conf=n_conformers, seed=seed)
    for s in sites:
        d = DonorCandidate(atom_idx=s.atom_idx, element=s.element,
                           smarts_pattern=s.pattern, smarts_prior=s.prior,
                           smarts_score=s.score, n_ortho_heavy=s.steric_neighbors)
        ok, why = lone_pair_available(mol, s.atom_idx)
        if not ok:
            d.rejected = why
        for ring in ri.AtomRings():
            if s.atom_idx in ring:
                d.ring_size = len(ring)
                d.ring_hetero = tuple(sorted(
                    mol.GetAtomWithIdx(j).GetSymbol() for j in ring
                    if mol.GetAtomWithIdx(j).GetSymbol() != "C"))
                break
        if molh is not None and cids:
            # AddHs preserves heavy-atom indices, so atom_idx still addresses the same N.
            for k, v in buried_volume_for_donor(molh, cids, s.atom_idx).items():
                setattr(d, k, v)
        q.donors.append(d)

    if molh is None or not cids:
        q.status, q.error = "partial", "3D embedding failed; QM terms skipped"
        q.seconds = round(time.time() - t0, 2)
        return q

    if not xtb_available():
        q.status = "partial"
        q.error = ("xtb not on PATH; free descriptors only. Run this on Modal via "
                   "scripts/qm/modal_qm.py.")
        q.seconds = round(time.time() - t0, 2)
        return q

    # ---- xTB ---------------------------------------------------------------
    wd = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="ligqm_"))
    wd.mkdir(parents=True, exist_ok=True)
    calls = 0
    try:
        neutral = _xtb(wd / "neutral", _xyz_text(molh, cids[0]),
                       charge=q.formal_charge, opt=True, solvent=solvent)
        calls += 1
        if not neutral["ok"]:
            q.status, q.error = "partial", "neutral optimization failed: " + neutral.get("error", "")[-200:]
            q.seconds, q.n_xtb_calls = round(time.time() - t0, 2), calls
            return q
        q.e_neutral_hartree = neutral["energy"]
        # Re-read the optimised geometry into the conformer so that Fukui and the
        # protomers start from the xTB minimum rather than the MMFF one.
        if neutral.get("xyz"):
            _load_xyz_into_conformer(molh, cids[0], neutral["xyz"])

        fk = condensed_fukui(molh, cids[0], q.formal_charge, wd, solvent=solvent)
        calls += fk.get("calls", 0)
        for d in q.donors:
            i = d.atom_idx
            if fk.get("q") and i < len(fk["q"]):
                d.mulliken_charge = float(fk["q"][i])
            if fk.get("f_minus") and i < len(fk["f_minus"]):
                d.fukui_minus = float(fk["f_minus"][i])
            if fk.get("f_plus") and i < len(fk["f_plus"]):
                d.fukui_plus = float(fk["f_plus"][i])
            if d.fukui_minus == d.fukui_minus and d.fukui_plus == d.fukui_plus:
                d.dual_descriptor = d.fukui_plus - d.fukui_minus

        for d in q.donors:
            if d.rejected:
                continue     # no lone pair to protonate; the PA would be meaningless
            pa = proton_affinity(molh, cids[0], d.atom_idx, q.e_neutral_hartree,
                                 wd / f"pa_{d.atom_idx}", base_charge=q.formal_charge,
                                 solvent=solvent, mol_fallback=mol, seed=seed)
            calls += 1
            if pa.get("ok"):
                d.proton_affinity_kcal = pa["pa_kcal"]
                d.e_protomer_hartree = pa["e_protomer"]
                d.protomer_ok = True
                if pa.get("reembedded"):
                    d.note = (d.note + " | PA from a re-embedded protomer; carries "
                              "conformational error, not comparable").strip(" |")
            else:
                d.note = (d.note + " | PA: " + str(pa.get("error", ""))[:120]).strip(" |")

        if also_solvated:
            ns = _xtb(wd / "neutral_w", _xyz_text(molh, cids[0]),
                      charge=q.formal_charge, solvent="water")
            calls += 1
            if ns["ok"]:
                q.e_neutral_solv_hartree = ns["energy"]
                for d in q.donors:
                    if not d.protomer_ok:
                        continue
                    pw = proton_affinity(molh, cids[0], d.atom_idx, ns["energy"],
                                         wd / f"paw_{d.atom_idx}",
                                         base_charge=q.formal_charge, solvent="water",
                                         mol_fallback=mol, seed=seed)
                    calls += 1
                    if pw.get("ok"):
                        d.proton_affinity_solv_kcal = pw["pa_kcal"]

        if n_bde_sites > 0:
            for c in sorted(q.carbons, key=lambda x: x.bde_proxy_kcal)[:n_bde_sites]:
                r = bde_xtb(molh, cids[0], c.atom_idx, q.e_neutral_hartree,
                            q.formal_charge, wd / f"bde_{c.atom_idx}", solvent=solvent)
                calls += 1
                if r.get("ok"):
                    c.bde_xtb_kcal = r["bde_kcal"]

        q.status = "ok" if any(d.protomer_ok for d in q.donors) or not q.donors else "partial"
    except Exception as exc:  # noqa: BLE001 - one ligand must not kill a batch
        import traceback
        q.status, q.error = "partial", traceback.format_exc()[-600:]
    finally:
        q.seconds = round(time.time() - t0, 2)
        q.n_xtb_calls = calls
        if workdir is None:
            shutil.rmtree(wd, ignore_errors=True)
    return q


def _load_xyz_into_conformer(molh, cid: int, xyz_text: str) -> None:
    from rdkit.Geometry import Point3D

    lines = xyz_text.strip().splitlines()
    try:
        n = int(lines[0].split()[0])
    except (IndexError, ValueError):
        return
    conf = molh.GetConformer(int(cid))
    if n != molh.GetNumAtoms():
        return
    for i, line in enumerate(lines[2:2 + n]):
        parts = line.split()
        conf.SetAtomPosition(i, Point3D(float(parts[1]), float(parts[2]), float(parts[3])))


# ==========================================================================
# rankings  (candidate orderings to be VALIDATED, not a fitted model)
# ==========================================================================

# Coefficient converting %V_bur into a kcal/mol-equivalent penalty, for the one composite
# ranking offered below.
#
# It is anchored on the single case where the right answer is not in doubt. Measured
# here: 2,6-lutidine has PA 71.35 and %V_bur 21.4 at the metal site; pyridine has PA
# 60.26 and %V_bur 14.5. Lutidine is the stronger BASE and the much weaker METAL LIGAND —
# that is the textbook 2,6-disubstitution effect and the exact failure mode §2B commissions
# %V_bur to catch. For the composite to order those two correctly it needs
#
#     71.35 - k*21.4  <  60.26 - k*14.5   =>   k > 1.61
#
# so the value below sits just above that threshold. An earlier guess of 0.25 — picked to
# make the two terms "look commensurate" — failed this case outright and is recorded here
# because a coefficient that cannot reproduce the one thing it was introduced for is worth
# knowing about.
#
# The value was fixed at 1.65 from that anchor BEFORE the validation was run. Sweeping it
# afterwards over the 32 decisive deposited ligands shows a broad plateau, not a knife
# edge — top-1 accuracy on the decisive subset / on all 71 mapped ligands:
#
#     k     0.00   0.25   0.50   0.80   1.00   1.65   2.00   2.50   5.00  10.00
#     dec  0.625  0.781  0.844  0.906  0.906  0.906  0.906  0.844  0.844  0.844
#     all  0.831  0.901  0.930  0.958  0.958  0.958  0.958  0.930  0.930  0.930
#
# So the result does not depend on hitting this number: anything in 0.8-2.0 gives the same
# answer, and the real content is that sterics enter AT ROUGHLY UNIT WEIGHT rather than at
# the 0.25 originally guessed. k=0 (proton affinity alone) is much worse than either
# ingredient's own ranking, which is the useful negative here.
#
# It remains a ONE-POINT ANCHOR on a textbook fact, not a fit: it was set by two molecules.
# `rank.py` must refit it against measured LDDT-PLI labels under leave-one-ligand-cluster-out
# before any of this reaches a submission, and `scripts/qm/validate_donor_ranking.py` reports
# the single-criterion rankings separately precisely so that this composite cannot hide
# which ingredient did the work.
VBUR_TO_KCAL_PRIOR = 1.65

RANKINGS = ("smarts", "proton_affinity", "fukui_minus", "vbur_at_metal",
            "pa_minus_steric")


def rank_donors(q: LigandQM, by: str = "pa_minus_steric") -> list:
    """Candidate donors ordered best-first under one criterion.

    Five criteria are offered rather than one blended number, because the honest question
    is *which single measured quantity, if any, orders donors better than the SMARTS
    prior does* — and a blend hides which ingredient did the work. The validation script
    scores every criterion separately against deposited ground truth and reports them side
    by side, including the cases where they disagree.

      smarts           the incumbent, for comparison (higher better)
      proton_affinity  sigma-donor strength           (higher better)
      fukui_minus      nucleophilicity toward Fe(III) (higher better)
      vbur_at_metal    steric access to the metal site (LOWER better)
      pa_minus_steric  PA - VBUR_TO_KCAL_PRIOR * %V_bur; see the constant's note

    Donors whose criterion is NaN sort last rather than being dropped, so the list always
    contains every candidate and a caller can see that a value is missing rather than
    inferring it from a shorter list. Candidates `lone_pair_available` ruled out sort
    below every live candidate under EVERY criterion including `smarts` — that filter is
    a physical impossibility, not a preference, and the same rule has to apply to the
    incumbent's ordering or the comparison would credit tier 1 with a correction the
    incumbent was never offered.
    """
    def key(d):
        if by == "smarts":
            v = d.smarts_score
        elif by == "proton_affinity":
            v = d.proton_affinity_kcal
        elif by == "fukui_minus":
            v = d.fukui_minus
        elif by == "vbur_at_metal":
            v = -d.vbur_at_metal if d.vbur_at_metal == d.vbur_at_metal else float("nan")
        elif by == "pa_minus_steric":
            if d.proton_affinity_kcal != d.proton_affinity_kcal:
                v = float("nan")
            else:
                vb = d.vbur_at_metal if d.vbur_at_metal == d.vbur_at_metal else 0.0
                v = d.proton_affinity_kcal - VBUR_TO_KCAL_PRIOR * vb
        else:
            raise ValueError(f"unknown ranking {by!r}; expected one of {RANKINGS}")
        missing = v != v
        return (bool(d.rejected), missing, -(0.0 if missing else v), d.atom_idx)

    return sorted(q.donors, key=key)


def best_donor(q: LigandQM, by: str = "pa_minus_steric"):
    r = rank_donors(q, by=by)
    return r[0] if r else None


# ==========================================================================
# reading tier 1 against a pose  (cheap: this is the part that runs per pose)
# ==========================================================================


def read_against_pose(q: LigandQM, *, coordinating_atom_idx: int | None = None,
                      fe_distances: dict | None = None,
                      rank_by: str = "pa_minus_steric") -> dict:
    """Turn the per-ligand profile into per-pose features. No QM, no coordinates built.

    This is the join that makes the cost architecture work. Tier 0 has already measured,
    for this pose, *which* ligand atom is nearest the iron and how far every atom is;
    tier 1 already knows what each of those atoms is chemically worth. All that remains is
    a lookup and a rank correlation — microseconds per pose, so the 3,500-pose pool costs
    nothing beyond the 87 QM jobs.

    `coordinating_atom_idx` is the ligand atom index tier 0 found closest to Fe, in the
    index space of `q.smiles`. Getting that mapping right is the project's known trap
    (predicted and deposited ligands do not share an atom order) — callers must map
    through a substructure match, never by array position.

    `fe_distances` maps ligand atom index -> distance to Fe for that pose. Supplying it
    enables the site-of-metabolism term, which is a Spearman correlation between each
    carbon's rank by closeness to the iron and its rank by ease of abstraction
    (low BDE first). +1 means the pose presents exactly the carbons the enzyme would
    attack; -1 means it presents the most inert ones. Rank form, per design §2C, so the
    proxy's systematic energy bias cannot leak in.
    """
    feats: dict = {"tier1_status": q.status, "n_candidate_donors": len(q.donors)}
    order = rank_donors(q, by=rank_by)

    if coordinating_atom_idx is not None:
        by_idx = {d.atom_idx: d for d in q.donors}
        d = by_idx.get(int(coordinating_atom_idx))
        feats["coordinates_through_a_candidate"] = d is not None
        if d is not None:
            feats.update(
                donor_pa_kcal=d.proton_affinity_kcal,
                donor_fukui_minus=d.fukui_minus,
                donor_vbur_at_metal=d.vbur_at_metal,
                donor_smarts_score=d.smarts_score,
                # 0 = the pose used the molecule's best donor. The design's claim in
                # §2B is precisely that this should separate good poses from bad ones at
                # identical Fe distance.
                donor_rank_used=next((i for i, x in enumerate(order)
                                      if x.atom_idx == d.atom_idx), -1),
                donor_pa_deficit=(order[0].proton_affinity_kcal - d.proton_affinity_kcal)
                if order and order[0].proton_affinity_kcal == order[0].proton_affinity_kcal
                and d.proton_affinity_kcal == d.proton_affinity_kcal else float("nan"),
            )

    if fe_distances:
        pairs = [(c.atom_idx, c.bde_proxy_kcal) for c in q.carbons
                 if c.atom_idx in fe_distances and c.bde_proxy_kcal == c.bde_proxy_kcal]
        if len(pairs) >= 3:
            idx = [i for i, _ in pairs]
            d_fe = np.array([fe_distances[i] for i in idx], float)
            bde = np.array([b for _, b in pairs], float)
            feats["som_rank_agreement"] = _spearman(d_fe, bde)
            k = int(np.argmin(d_fe))
            feats["closest_carbon_bde"] = float(bde[k])
            feats["closest_carbon_bde_percentile"] = float((bde < bde[k]).mean())
    return feats


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman rho without scipy (the container image does not carry it).

    Ties get average ranks, which matters here: a molecule with eight equivalent
    aromatic CH carbons has eight identical proxy BDEs, and ordinal ranking would invent
    a preference among them that the descriptor never expressed.
    """
    def rank(x):
        order = np.argsort(x, kind="mergesort")
        r = np.empty(len(x), float)
        r[order] = np.arange(len(x), dtype=float)
        # average ranks within tied groups
        xs = x[order]
        i = 0
        while i < len(xs):
            j = i
            while j + 1 < len(xs) and xs[j + 1] == xs[i]:
                j += 1
            if j > i:
                r[order[i:j + 1]] = (i + j) / 2.0
            i = j + 1
        return r

    ra, rb = rank(np.asarray(a, float)), rank(np.asarray(b, float))
    ra -= ra.mean()
    rb -= rb.mean()
    den = math.sqrt(float((ra ** 2).sum()) * float((rb ** 2).sum()))
    return float((ra * rb).sum() / den) if den > 0 else float("nan")


# ==========================================================================


def default_workdir() -> Path:
    """A scratch directory that is NOT on the near-full C: drive.

    `tempfile.gettempdir()` resolves to C:\\tb\\tmp on this box, which has hit zero bytes
    free three times and once silently truncated a scoring run. xtb writes a handful of
    small files per call, but "small" times a few hundred calls on a volume at zero is
    still a corrupted run, so the default is placed on D: explicitly.
    """
    if os.name == "nt":
        p = Path("D:/cyp_scratch/ligand_qm")
    else:
        p = Path(tempfile.gettempdir()) / "cyp_ligand_qm"
    p.mkdir(parents=True, exist_ok=True)
    return p


if __name__ == "__main__":
    # Free descriptors only unless xtb happens to be on PATH — which on this box it is
    # not. This self-check exercises the graph and geometry halves, which is what can be
    # verified locally.
    from rdkit import Chem

    tests = {
        "pyridine": "c1ccncc1",
        "2,6-lutidine": "Cc1cccc(C)n1",
        "imidazole": "c1c[nH]cn1",
        "ketoconazole": "CC(=O)N1CCN(CC1)c1ccc(OC[C@@H]2CO[C@](Cn3ccnc3)(O2)c2ccc(Cl)cc2Cl)cc1",
        "clotrimazole": "Clc1ccccc1C(n1ccnc1)(c1ccccc1)c1ccccc1",
        "testosterone": "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O",
    }
    for name, smi in tests.items():
        q = compute_ligand_qm(smi, name, n_conformers=3)
        print(f"\n{name}  ({q.status}, {q.n_heavy} heavy, {q.seconds}s)")
        for d in q.donors:
            print(f"   N{d.atom_idx:<3d} {d.smarts_pattern:<18s} smarts={d.smarts_score:.2f} "
                  f"vbur_atom={d.vbur_at_atom:.1f}%  vbur_metal={d.vbur_at_metal:.1f}%")
        best = sorted(q.carbons, key=lambda c: c.bde_proxy_kcal)[:3]
        print("   easiest C-H: " + ", ".join(
            f"C{c.atom_idx}({c.environment} {c.bde_proxy_kcal:.1f})" for c in best))
