"""Extract compact, ALIGNED PDB files for the dashboard's 3D viewer.

The viewer overlays a predicted pose on the crystal truth, so the two must be in the same
reference frame or the comparison is theatre. The PXR campaign lost real time to exactly
this: some engine exports landed ~20 A off the crystal and inflated every agreement number
about twofold before anyone noticed. `align_by_residue` matches by residue NUMBER, never
by array position, which is what makes the overlay mean something.

Output is trimmed on purpose - backbone + heme + ligand + pocket side chains. A full
mmCIF is ~800 KB and most of it is solvent and alternate conformers the viewer will never
draw.

    python scripts/viz/extract_for_viewer.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

OUT = REPO / "data" / "processed" / "viewer"
BACKBONE = {"N", "CA", "C", "O"}


def write_pdb(path: Path, atoms: list[dict]) -> int:
    lines = []
    for i, a in enumerate(atoms, start=1):
        x, y, z = a["xyz"]
        name = a["name"]
        # PDB atom-name column rules: 4-char field, element right-justified in cols 13-14
        nm = f"{name:<4}" if len(name) >= 4 else f" {name:<3}"
        res3 = a['res'][:3]          # PDB residue column is 3 wide; "LIG1" shifts everything
        lines.append(
            f"{a.get('rec','ATOM'):<6}{i:>5} {nm}{'':1}{res3:>3} {a.get('chain','A')}"
            f"{a['resi']:>4}{'':4}{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}"
            f"{'':10}{a.get('elem',''):>2}"
        )
    lines.append("END")
    path.write_text("\n".join(lines))
    return len(atoms)


# Match the ligand STRUCTURALLY, not by a name list. The prediction calls it "LIG1",
# which an exact-match set of {"LIG", ...} silently missed - and a 4-character residue
# name also overflows the PDB residue column, shifting every field after it. Anything
# hetero that is not the heme and not solvent is the ligand.
NON_LIGAND_HET = {"HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "ACT", "CL", "NA", "MG"}


def is_ligand(res_name: str, heme_aliases) -> bool:
    n = res_name.strip().upper()
    return n not in heme_aliases and n not in NON_LIGAND_HET


def one_pocket_copy(atoms: list[dict], tag: str) -> list[dict]:
    """Keep ONE ligand copy - the one in the active site - and name it consistently.

    8SO1 deposits 42 caffeine atoms: three copies of a 14-atom molecule. Only one sits in
    the pocket; the others are peripheral surface sites. Drawing all three would suggest
    the crystal has three binding poses, which is the same `(pdb, chain, lig)`-is-not-a-
    unique-key trap that once manufactured 125 duplicate rows in this repo.

    The prediction calls its ligand LIG on a blank chain; the crystal calls it CFF. They
    are the same molecule and the viewer should not have to know that, so both are
    normalised to CFF on chain L.
    """
    fe = next((np.array(a["xyz"]) for a in atoms
               if a["res"] == "HEM" and a["name"].upper() == "FE"), None)
    from cypstruct.targets import HEME_ALIASES
    het = [a for a in atoms if a.get("rec") == "HETATM"]
    lig = [a for a in het if is_ligand(a["res"], HEME_ALIASES)]
    rest = [a for a in atoms if a not in lig]
    if not lig:
        return atoms
    # group copies by (chain, residue number) and keep the group nearest the iron
    groups: dict[tuple, list[dict]] = {}
    for a in lig:
        groups.setdefault((a.get("chain", " "), a["resi"]), []).append(a)
    if fe is not None and len(groups) > 1:
        def dist(g):
            return min(float(np.linalg.norm(np.array(a["xyz"]) - fe)) for a in g)
        key = min(groups, key=lambda k: dist(groups[k]))
        dropped = sum(len(v) for k, v in groups.items() if k != key)
        print(f"    {tag}: {len(groups)} ligand copies, keeping the one "
              f"{dist(groups[key]):.1f} A from Fe, dropping {dropped} peripheral atoms")
        lig = groups[key]
    for a in lig:
        a["res"], a["chain"] = "CFF", "L"
    return rest + lig


def kabsch_by_resnum(mobile, ref) -> tuple[np.ndarray, np.ndarray]:
    """Rigid transform putting `mobile` in `ref`'s frame, matched by residue NUMBER.

    `pose.align_by_residue` returns an aligned Complex rather than the transform, and the
    viewer needs the transform so it can move the whole gemmi structure (side chains,
    heme, ligand) and not just the arrays the Complex carries. Matching on residue number
    rather than array position is the non-negotiable part: the crystal is author-numbered
    (26..496 here) and the prediction is numbered from 1, so index-for-index pairing would
    silently align the wrong residues and still return a plausible-looking RMSD.
    """
    def ca_map(m):
        out = {}
        for k, xyz in zip(m.prot_key, m.prot_xyz):
            resi, name = (k[1], k[2]) if isinstance(k, tuple) and len(k) > 2 else (None, None)
            if name == "CA" and resi is not None:
                out[resi] = np.asarray(xyz, float)
        return out

    a, b = ca_map(mobile), ca_map(ref)
    shared = sorted(set(a) & set(b))
    if len(shared) < 30:
        raise ValueError(f"only {len(shared)} shared CA residues - refusing to align")
    X = np.array([a[r] for r in shared])
    Y = np.array([b[r] for r in shared])
    xc, yc = X.mean(0), Y.mean(0)
    U, S, Vt = np.linalg.svd((X - xc).T @ (Y - yc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    rms = float(np.sqrt(((((R @ (X - xc).T).T + yc) - Y) ** 2).sum(1).mean()))
    print(f"    aligned on {len(shared)} CA residues, backbone RMSD {rms:.2f} A")
    return R, yc - R @ xc


def collect(struct, keep_resis: set[int] | None, tag: str) -> list[dict]:
    """Pull backbone + hetero (HEM/ligand) + pocket side chains out of a gemmi structure."""
    import gemmi
    from cypstruct.targets import HEME_ALIASES

    atoms = []
    model = struct[0]
    for chain in model:
        for res in chain:
            nm = res.name.strip().upper()
            is_het = nm in HEME_ALIASES or res.het_flag == "H"
            if nm in ("HOH", "WAT", "DOD"):
                continue
            for at in res:
                if at.element == gemmi.Element("H"):
                    continue
                keep = (
                    is_het
                    or at.name in BACKBONE
                    or (keep_resis and res.seqid.num in keep_resis)
                )
                if not keep:
                    continue
                atoms.append({
                    "rec": "HETATM" if is_het else "ATOM",
                    "name": at.name, "res": nm, "chain": chain.name[:1] or "A",
                    "resi": res.seqid.num,
                    "xyz": (at.pos.x, at.pos.y, at.pos.z),
                    "elem": at.element.name.upper(),
                })
    return atoms


def main() -> int:
    import gemmi

    from cypstruct import pose as P
    from cypstruct import targets as T

    OUT.mkdir(parents=True, exist_ok=True)
    crystal_path = REPO / "data" / "reference" / "rcsb" / "8SO1.cif"
    pool = Path("D:/cyp_scratch/val87b_unsteered/CFF__unsteered__s1")
    jobs = {
        "crystal": crystal_path,
        "good": pool / "input_model_0.cif",     # LDDT-PLI 0.908
        "bad": pool / "input_model_17.cif",     # LDDT-PLI 0.320, same ligand, same pool
    }
    pocket = set(T.pocket_atoms()) if hasattr(T, "pocket_atoms") else set()

    ref_mo = P.load_structure(crystal_path)
    manifest = {}
    for tag, path in jobs.items():
        st = gemmi.read_structure(str(path))
        st.remove_alternative_conformations()
        st.remove_ligands_and_waters() if False else None
        if tag != "crystal":
            # put the prediction in the crystal's frame; residue-number matching only
            mo = P.load_structure(path)
            try:
                R, t = kabsch_by_resnum(mo, ref_mo)
                for model in st:
                    for chain in model:
                        for res in chain:
                            for at in res:
                                v = np.array([at.pos.x, at.pos.y, at.pos.z])
                                w = R @ v + t
                                at.pos = gemmi.Position(*w)
            except Exception as exc:
                print(f"  {tag}: alignment failed ({type(exc).__name__}) - "
                      "NOT writing, an unaligned overlay would be misleading")
                continue
        atoms = collect(st, pocket, tag)
        atoms = one_pocket_copy(atoms, tag)
        n = write_pdb(OUT / f"cff_{tag}.pdb", atoms)
        manifest[tag] = {"atoms": n, "src": str(path)}
        print(f"  {tag}: {n} atoms -> cff_{tag}.pdb")

    import json
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
