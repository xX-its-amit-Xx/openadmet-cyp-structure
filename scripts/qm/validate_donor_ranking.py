"""Does tier-1 QM pick the nitrogen that actually coordinates, better than SMARTS does?

This is the honest gate on tier 1, and it is set up to be able to return "no".

**What the incumbent's 98.8% actually measures, and why tier 1 cannot be compared to it
directly.** `scripts/structure/validate_donor_prediction.py` reports recall 0.988 = 82/83:
of the 83 deposited CYP3A4 ligands observed coordinating the iron, 82 were *classified*
`type_II`, i.e. the predictor offered any competent donor at all. That is a **molecule-level
class recall**. It says nothing about whether the top-ranked atom is the right atom, because
the classifier never looked. Quoting 98.8% as a donor-ranking accuracy would be a category
error, and beating or missing it with an atom-level number would be comparing two different
quantities.

So this script measures the atom-level question for BOTH methods on the same footing:

    given the deposited structure, which ligand atom is 1.9-2.45 A from the heme iron,
    and does the method rank that atom first?

**The denominator is reported honestly, and it is small.** Measured before any compute was
spent: of the 83 coordinated ligands, **50 have exactly one candidate donor nitrogen**.
Every method — SMARTS, proton affinity, Fukui, %V_bur — returns the same single atom for
those, so they are 50 free correct answers that inflate any accuracy quoted over the full
set. Only **32 ligands have two or more candidates**, and those 32 are the entire
addressable surface. Headline numbers here are quoted on that decisive subset, with the
full-set number alongside so neither can be cherry-picked.

Correctness is **symmetry-aware**. Two atoms related by molecular symmetry (the two
nitrogens of a symmetric bis-imidazole, say) are the same answer chemically, and RDKit's
canonical ranking with `breakTies=False` gives symmetry-equivalent atoms an equal rank. A
prediction counts as correct when it shares a rank with the observed donor, so the score
is not decided by which of two indistinguishable atoms the deposition happened to name.

Writes data/processed/tier1_donor_ranking_validation.json.

Usage:
    python scripts/qm/validate_donor_ranking.py ground-truth     # local, free
    python scripts/qm/validate_donor_ranking.py score --qm-dir data/processed/qm/tier1
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.chem import coordinating_atoms, standardize  # noqa: E402
from cypstruct.paths import DATA_PROCESSED, REFERENCE  # noqa: E402
from cypstruct.qmscore.ligand_qm import (  # noqa: E402
    RANKINGS,
    LigandQM,
    lone_pair_available,
    rank_donors,
)

GT_PATH = DATA_PROCESSED / "cyp3a4_donor_ground_truth.json"
OUT_PATH = DATA_PROCESSED / "tier1_donor_ranking_validation.json"

# The tier-0 coordination window, from geometry.py. Kept identical on purpose: the
# ground truth for "this atom coordinates" must be the same test the scorer applies.
COORD_LO, COORD_HI = 1.90, 2.45


# ==========================================================================
# ground truth: which ATOM, not just which molecule
# ==========================================================================


def build_ground_truth() -> dict:
    """For each coordinated ligand, the PDB atom name of the atom bonded to the iron.

    Read from the locally cached mmCIFs rather than from the parquet, because the parquet
    records only the donor ELEMENT (`donor_element`), which cannot distinguish the two
    nitrogens of an imidazole — exactly the discrimination under test.

    Per-chain, never per-file: CYP3A4 deposits as monomer, dimer and trimer, and pairing a
    ligand with a heme from a different copy produces a confident wrong answer. Same trap
    `build_reference_set.py` is written around.
    """
    import gemmi

    from cypstruct.targets import HEME_ALIASES, IGNORE_HET

    df = pd.read_parquet(DATA_PROCESSED / "cyp3a4_reference_geometry.parquet")
    coord = df[df.binding_class == "type_II_coordinated"]
    print(f"{len(coord)} coordinated chain-ligand rows, "
          f"{coord.ligand.nunique()} unique ligands")

    per_lig: dict[str, Counter] = {}
    dists: dict[str, list] = {}
    for pdb_id, sub in coord.groupby("pdb_id"):
        path = REFERENCE / "rcsb" / f"{pdb_id}.cif"
        if not path.exists():
            continue
        try:
            st = gemmi.read_structure(str(path))
            st.setup_entities()
            st.remove_alternative_conformations()
            st.remove_hydrogens()
        except Exception as exc:
            print(f"  [skip] {pdb_id}: {exc}")
            continue
        if len(st) == 0:
            continue
        want = {(r.chain, r.ligand) for _i, r in sub.iterrows()}
        for chain in st[0]:
            fe = None
            hets = []
            for res in chain:
                rn = res.name.strip().upper()
                if rn in HEME_ALIASES:
                    for at in res:
                        if at.name.strip().upper() == "FE":
                            fe = np.array([at.pos.x, at.pos.y, at.pos.z])
                elif rn not in IGNORE_HET:
                    info = gemmi.find_tabulated_residue(rn)
                    if info and info.is_amino_acid():
                        continue
                    hets.append(res)
            if fe is None:
                continue
            for res in hets:
                rn = res.name.strip().upper()
                if (chain.name, rn) not in want:
                    continue
                best = None
                for at in res:
                    d = float(np.linalg.norm(
                        np.array([at.pos.x, at.pos.y, at.pos.z]) - fe))
                    if best is None or d < best[0]:
                        best = (d, at.name.strip(), at.element.name)
                if best and COORD_LO <= best[0] <= COORD_HI:
                    per_lig.setdefault(rn, Counter())[(best[1], best[2])] += 1
                    dists.setdefault(rn, []).append(best[0])

    gt = {}
    for code, names in per_lig.items():
        (atom_name, element), n = names.most_common(1)[0]
        gt[code] = {"atom_name": atom_name, "element": element, "n_obs": int(n),
                    "n_chains": int(sum(names.values())),
                    "agrees_across_chains": len(names) == 1,
                    "mean_fe_dist": round(float(np.mean(dists[code])), 3)}
    GT_PATH.write_text(json.dumps(gt, indent=1))
    print(f"ground truth for {len(gt)} ligands -> {GT_PATH}")
    disagree = [k for k, v in gt.items() if not v["agrees_across_chains"]]
    if disagree:
        print(f"  {len(disagree)} ligands name different donors in different chains "
              f"(majority vote used): {disagree[:10]}")
    return gt


def map_atom_name_to_index(code: str, smiles: str, atom_name: str) -> int | None:
    """PDB atom name -> atom index in the standardised SMILES the descriptors used.

    The two do not share an index space, and assuming they do is the trap that has
    already cost this project half an LDDT-PLI (see `chem.boltz_atom_name`). The route
    taken here is the safe one: read the deposited residue's own connectivity, impose the
    SMILES bond orders on it with `AssignBondOrdersFromTemplate`, then recover the
    correspondence with a substructure match. Returns None — never a guess — when any
    step fails, and the caller counts those as unmapped rather than as wrong.
    """
    import gemmi
    from rdkit import Chem
    from rdkit.Chem import AllChem

    template = Chem.MolFromSmiles(smiles)
    if template is None:
        return None

    # Find any deposited copy of this ligand to read connectivity from.
    df = pd.read_parquet(DATA_PROCESSED / "cyp3a4_reference_geometry.parquet")
    rows = df[(df.ligand == code) & (df.binding_class == "type_II_coordinated")]
    for pdb_id in rows.pdb_id.unique():
        path = REFERENCE / "rcsb" / f"{pdb_id}.cif"
        if not path.exists():
            continue
        try:
            st = gemmi.read_structure(str(path))
            st.setup_entities()
            st.remove_alternative_conformations()
            st.remove_hydrogens()
        except Exception:
            continue
        for chain in st[0]:
            for res in chain:
                if res.name.strip().upper() != code:
                    continue
                lines, names = [], []
                for i, at in enumerate(res):
                    el = at.element.name
                    # The residue-name field is columns 18-20, exactly THREE characters.
                    # The PDB's newer five-character CCD codes (A1A06, A1ASO, ... — which
                    # is what the 2024-25 CYP3A4 depositions use) overflow it and shift
                    # every subsequent column right by two, so the coordinates parse as
                    # garbage and the molecule silently fails to match its own template.
                    # This cost 18 of the 27 ligands that first reported as unmappable.
                    # The name is cosmetic here, so truncating is safe; what must not
                    # move is the column alignment.
                    lines.append(
                        f"HETATM{i + 1:5d} {at.name.strip()[:4]:<4s} {code[:3]:>3s} A   1    "
                        f"{at.pos.x:8.3f}{at.pos.y:8.3f}{at.pos.z:8.3f}  1.00  0.00"
                        f"          {el.upper():>2s}")
                    names.append(at.name.strip())
                block = "\n".join(lines) + "\nEND\n"
                pdbmol = Chem.MolFromPDBBlock(block, sanitize=False,
                                              removeHs=True, proximityBonding=True)
                if pdbmol is None or pdbmol.GetNumAtoms() != template.GetNumAtoms():
                    continue
                try:
                    fixed = AllChem.AssignBondOrdersFromTemplate(template, pdbmol)
                except Exception:
                    continue
                # match[i_template] = i_pdb
                match = fixed.GetSubstructMatch(template)
                if not match:
                    continue
                try:
                    pdb_idx = names.index(atom_name)
                except ValueError:
                    continue
                for t_idx, p_idx in enumerate(match):
                    if p_idx == pdb_idx:
                        return t_idx
    return None


# ==========================================================================
# scoring
# ==========================================================================


def symmetry_classes(smiles: str) -> list[int]:
    """Canonical ranks with ties UNBROKEN: equal rank == symmetry-equivalent atom."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    return list(AllChem.CanonicalRankAtoms(mol, breakTies=False)) if mol else []


def score(qm_dir: str | Path) -> dict:
    """Atom-level top-1 accuracy for every ranking criterion, on the same ligands."""
    qm_dir = Path(qm_dir)
    if not GT_PATH.exists():
        raise SystemExit(f"missing {GT_PATH}; run `ground-truth` first")
    gt = json.loads(GT_PATH.read_text())
    smi_map = json.loads((REFERENCE / "ligand_smiles.json").read_text())

    # `modal_qm.plan` deduplicates jobs on the standardised SMILES, so a ligand that is
    # chemically identical to another PDB code has no file of its own. Index what we have
    # by SMILES as well as by code, or those ligands would be miscounted as "QM missing"
    # when in fact they were computed under a synonym.
    by_smiles = {}
    for p in sorted(qm_dir.glob("*.json")) if qm_dir.exists() else []:
        try:
            rec = LigandQM.load(p)
        except Exception:
            continue
        by_smiles.setdefault(rec.smiles, rec)

    per_ligand, unmapped, missing_qm = [], [], []
    for code, g in sorted(gt.items()):
        raw = smi_map.get(code, "")
        smi = standardize(raw) if raw else ""
        if not smi:
            unmapped.append({"ligand": code, "why": "no SMILES"})
            continue
        true_idx = map_atom_name_to_index(code, smi, g["atom_name"])
        if true_idx is None:
            unmapped.append({"ligand": code, "why": f"could not map {g['atom_name']}"})
            continue
        sym = symmetry_classes(smi)

        qpath = qm_dir / f"{code}.json"
        q = LigandQM.load(qpath) if qpath.exists() else by_smiles.get(smi)
        if q is None:
            missing_qm.append(code)

        # The incumbent, exactly as it stands today (no tier-1 filter applied).
        smarts_sites = coordinating_atoms(smi, top_k=99)
        rec = {
            "ligand": code, "true_atom_name": g["atom_name"],
            "true_idx": true_idx, "fe_dist": g["mean_fe_dist"],
            "n_candidates": len(smarts_sites),
            "decisive": len(smarts_sites) >= 2,
            "smarts_top": smarts_sites[0].atom_idx if smarts_sites else None,
            "smarts_offered_true": any(s.atom_idx == true_idx for s in smarts_sites),
        }
        rec["smarts_correct"] = bool(
            smarts_sites and sym and sym[smarts_sites[0].atom_idx] == sym[true_idx])

        # The free half of tier 1: the lone-pair filter alone, no QM at all.
        live = [s for s in smarts_sites if lone_pair_available(
            __import__("rdkit").Chem.MolFromSmiles(smi), s.atom_idx)[0]]
        rec["n_candidates_after_filter"] = len(live)
        rec["filter_correct"] = bool(
            live and sym and sym[live[0].atom_idx] == sym[true_idx])

        if q is not None:
            rec["qm_status"] = q.status
            for crit in RANKINGS:
                order = rank_donors(q, by=crit)
                top = order[0].atom_idx if order else None
                rec[f"{crit}_top"] = top
                rec[f"{crit}_correct"] = bool(
                    top is not None and sym and sym[top] == sym[true_idx])
        per_ligand.append(rec)

    def acc(rows, key):
        rows = [r for r in rows if key in r]
        if not rows:
            return None
        return {"n": len(rows), "correct": sum(bool(r[key]) for r in rows),
                "accuracy": round(sum(bool(r[key]) for r in rows) / len(rows), 4)}

    decisive = [r for r in per_ligand if r["decisive"]]
    with_qm = [r for r in per_ligand if "qm_status" in r]
    decisive_qm = [r for r in decisive if "qm_status" in r]

    summary = {
        "n_ligands_with_ground_truth": len(gt),
        "n_scored": len(per_ligand),
        "n_unmapped": len(unmapped),
        "n_missing_qm": len(missing_qm),
        "note_on_the_988_baseline": (
            "98.8% (82/83) in donor_prediction_validation.json is MOLECULE-LEVEL class "
            "recall - did the predictor call the ligand type II at all. It is not a "
            "donor-ranking accuracy and is not directly comparable with anything below. "
            "The comparable incumbent number is smarts_top1 on the same ligands."),
        "all_ligands": {
            "smarts_top1": acc(per_ligand, "smarts_correct"),
            "lonepair_filter_top1": acc(per_ligand, "filter_correct"),
            **{f"{c}_top1": acc(with_qm, f"{c}_correct") for c in RANKINGS},
        },
        "decisive_subset_multi_candidate_only": {
            "smarts_top1": acc(decisive, "smarts_correct"),
            "lonepair_filter_top1": acc(decisive, "filter_correct"),
            **{f"{c}_top1": acc(decisive_qm, f"{c}_correct") for c in RANKINGS},
        },
        "single_candidate_ligands": sum(1 for r in per_ligand if not r["decisive"]),
        "unmapped": unmapped,
        "missing_qm": missing_qm,
    }

    # Where the methods disagree is the only place any of this can matter.
    dis = [r for r in decisive_qm
           if r.get("smarts_correct") != r.get("pa_minus_steric_correct")]
    summary["disagreements_smarts_vs_pa_minus_steric"] = [
        {k: r.get(k) for k in ("ligand", "true_atom_name", "true_idx", "smarts_top",
                               "pa_minus_steric_top", "smarts_correct",
                               "pa_minus_steric_correct")} for r in dis]

    OUT_PATH.write_text(json.dumps({"summary": summary, "per_ligand": per_ligand},
                                   indent=1))
    return summary


def verdict(summary: dict) -> str:
    """State plainly whether tier 1 earned its place. A negative is a result."""
    d = summary["decisive_subset_multi_candidate_only"]
    base = d.get("smarts_top1") or {}
    lines = [f"decisive subset: n={base.get('n')} multi-candidate coordinated ligands",
             f"  incumbent SMARTS top-1 : {base.get('accuracy')}"]
    best, best_acc = None, -1.0
    for c in RANKINGS:
        a = d.get(f"{c}_top1")
        if not a:
            continue
        lines.append(f"  {c:<18s}top-1 : {a['accuracy']}  ({a['correct']}/{a['n']})")
        if c != "smarts" and a["accuracy"] > best_acc:
            best, best_acc = c, a["accuracy"]
    f = d.get("lonepair_filter_top1")
    if f:
        lines.append(f"  lone-pair filter alone : {f['accuracy']}  (free, no QM)")
    if best is not None and base.get("accuracy") is not None:
        delta = best_acc - base["accuracy"]
        lines.append("")
        if delta > 0:
            lines.append(f"TIER 1 WINS on the decisive subset: {best} beats SMARTS by "
                         f"{delta:+.3f} ({best_acc} vs {base['accuracy']}).")
        elif delta == 0:
            lines.append(f"TIER 1 TIES the incumbent ({best_acc}). No gain to report. "
                         "Per the repo's honest-gating rule this is a NEGATIVE: the QM "
                         "is absorbed by the existing prior and adds nothing at "
                         "inference. Log it and do not ship the term.")
        else:
            lines.append(f"TIER 1 LOSES to the incumbent by {delta:+.3f}. Negative "
                         "result - report it, do not bury it.")
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["ground-truth", "score"])
    ap.add_argument("--qm-dir", default=str(DATA_PROCESSED / "qm" / "tier1"))
    a = ap.parse_args()
    if a.cmd == "ground-truth":
        build_ground_truth()
    else:
        s = score(a.qm_dir)
        print(json.dumps(s["all_ligands"], indent=1))
        print()
        print(verdict(s))
        print(f"\nwritten -> {OUT_PATH}")
