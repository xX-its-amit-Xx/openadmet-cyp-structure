"""Generate the interactive structure dashboard.

Written as a generator rather than a hand-authored page because ~500 KB of real atomic
coordinates have to be embedded exactly, and because every number on the page should come
from the measured artifacts rather than from memory. If a figure here disagrees with the
repo, the repo wins and this script is the bug.

    python scripts/viz/build_dashboard.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

VIEW = REPO / "data" / "processed" / "viewer"
# The Vercel dashboard is canonical. This page is generated straight into the shared site
# so it ships with every deploy instead of living as a separate artifact that drifts.
SITE = Path(r"D:\Users\ashenoy00000\.windsurf\OpenADMET-cyp-challenge\site")
OUT = SITE / "structures.html"


def ligand_frames() -> dict:
    """Crystal and predicted ligand coordinates in a SHARED atom order.

    The morph between truth and prediction is only chemically meaningful if atom i in one
    frame is the same atom in the other. Predicted and crystal ligands list their atoms
    differently - this repo has a standing note about it, because index-for-index pairing
    halves LDDT-PLI and looks merely disappointing rather than wrong. `best_ligand_mapping`
    resolves it, including molecular symmetry.
    """
    from cypstruct import pose as P
    smi = json.loads((REPO / "data/reference/ligand_smiles.json").read_text())["CFF"]
    ref = P.load_structure(REPO / "data/reference/rcsb/8SO1.cif")
    pool = Path("D:/cyp_scratch/val87b_unsteered/CFF__unsteered__s1")
    out = {}
    for tag, fn in (("good", "input_model_0.cif"), ("bad", "input_model_17.cif")):
        mo = P.load_structure(pool / fn)
        perm = P.best_ligand_mapping(smi, mo, ref)
        out[tag] = {
            "lddt": round(float(P.lddt_pli(mo, ref, lig_perm=perm)), 3),
            "rmsd": round(float(P.bisy_rmsd(mo, ref, lig_perm=perm)), 2),
        }
    return out


def morph_frames() -> dict:
    """Crystal and predicted ligand coordinates in ONE shared atom order, plus a
    ligand-only PDB written in that same order.

    The viewer interpolates between truth and prediction, which is only meaningful if
    atom i is the same atom in both frames. `best_ligand_mapping` gives that
    correspondence (`perm[i]` = reference atom for model atom i), symmetry-aware. The
    ligand-only PDB is emitted here rather than parsed out of the big file so the JS
    array order matches by construction instead of by assumption - the ligands carry no
    atom names to re-key on, so order is the only handle and it must not be guessed.
    """
    import json as _j
    sys.path.insert(0, str(REPO / "scripts" / "viz"))
    from extract_for_viewer import kabsch_by_resnum

    from cypstruct import pose as P
    smi = _j.loads((REPO / "data/reference/ligand_smiles.json").read_text())["CFF"]
    ref = P.load_structure(REPO / "data/reference/rcsb/8SO1.cif")
    pool = Path("D:/cyp_scratch/val87b_unsteered/CFF__unsteered__s1")
    out = {}
    for tag, fn in (("good", "input_model_0.cif"), ("bad", "input_model_17.cif")):
        mo = P.load_structure(pool / fn)
        perm = P.best_ligand_mapping(smi, mo, ref)
        R, t = kabsch_by_resnum(mo, ref)
        model_xyz = (np.asarray(mo.lig_xyz) @ R.T) + t        # into the crystal frame
        ref_xyz = np.asarray(ref.lig_xyz)
        if perm is None:
            print(f"    {tag}: no atom mapping - morph disabled for this pose")
            continue
        paired = ref_xyz[np.asarray(perm, int)]               # crystal atom for model atom i
        elems = list(mo.lig_elem)
        lines = []
        for i, (e, xyz) in enumerate(zip(elems, model_xyz), start=1):
            nm = f"{e}{i}"
            lines.append(
                f"HETATM{i:>5} {nm:<4} CFF L   1    "
                f"{xyz[0]:>8.3f}{xyz[1]:>8.3f}{xyz[2]:>8.3f}  1.00  0.00"
                f"{'':10}{e:>2}")
        lines.append("END")
        out[tag] = {
            "pdb": chr(10).join(lines),
            "to": [[round(float(v), 3) for v in p_] for p_ in model_xyz],
            "from": [[round(float(v), 3) for v in p_] for p_ in paired],
        }
        drift = float(np.linalg.norm(model_xyz - paired, axis=1).mean())
        print(f"    {tag}: morph over {len(elems)} atoms, mean atom drift {drift:.2f} A")
    return out


def pdb_text(tag: str) -> str:
    return (VIEW / f"cff_{tag}.pdb").read_text()


def lig_xyz(tag: str) -> list[list[float]]:
    out = []
    for ln in pdb_text(tag).splitlines():
        if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() == "CFF":
            out.append([round(float(ln[30:38]), 3), round(float(ln[38:46]), 3),
                        round(float(ln[46:54])), ])
    return out


def main() -> int:
    from cypstruct import targets as T

    scores = ligand_frames()
    pockets = sorted(set(T.pocket_atoms()))
    data = {
        "pdb": {t: pdb_text(t) for t in ("crystal", "good", "bad")},
        "morph": morph_frames(),
        "scores": scores,
        "pocket": pockets,
        "axialCys": T.CYP_TARGETS["cyp3a4"]["axial_cys"],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tpl = (Path(__file__).parent / "structures_template.html").read_text(encoding="utf-8")
    html = tpl.replace("/*__DATA__*/", json.dumps(data))
    OUT.write_text(html, encoding="utf-8")
    kb = OUT.stat().st_size / 1024
    print(f"wrote {OUT}  ({kb:.0f} KB)")
    print(f"  scores: {scores}")
    print(f"  pocket residues: {len(pockets)}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
