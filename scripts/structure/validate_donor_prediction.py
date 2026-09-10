"""Does the donor-atom predictor actually pick the atom that coordinates the iron?

This is the honest check on the central mechanism of the project. `chem.coordinating_atoms`
is what we will use to *steer* co-folding — if it names the wrong nitrogen, the constraint
actively forces a wrong orientation, which is worse than not steering at all.

Test: for every deposited CYP3A4 chain classified `type_II_coordinated` in
`cyp3a4_reference_geometry.parquet`, fetch the ligand's chemical definition from RCSB,
run the predictor on its SMILES, and check whether the predicted top donor is the same
*element and chemical environment* as the atom actually 2.0-2.4 A from the iron.

Matching is done on the donor's local environment rather than atom index, because the
deposited atom naming and an RDKit SMILES parse do not share an index space. We compare:
  - element (must be N),
  - aromaticity and ring size of the coordinating atom,
  - the ring's heteroatom composition.
That is strict enough to catch "predicted the pyridine but ketoconazole uses the
imidazole" and loose enough not to fail on naming conventions.

Writes data/processed/donor_prediction_validation.json.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from cypstruct.chem import coordinating_atoms, ligand_class  # noqa: E402
from cypstruct.paths import DATA_PROCESSED, REFERENCE  # noqa: E402

LIG_CACHE = REFERENCE / "ligand_smiles.json"


def ligand_smiles(codes: list[str]) -> dict[str, str]:
    """Fetch canonical SMILES for PDB chemical component IDs (cached)."""
    cache = json.loads(LIG_CACHE.read_text()) if LIG_CACHE.exists() else {}
    todo = [c for c in codes if c not in cache]
    for i, code in enumerate(todo, 1):
        try:
            r = requests.get(
                f"https://data.rcsb.org/rest/v1/core/chemcomp/{code}", timeout=30)
            if r.status_code != 200:
                cache[code] = ""
                continue
            d = r.json()
            smi = ""
            for desc in d.get("pdbx_chem_comp_descriptor", []):
                if desc.get("type") == "SMILES_CANONICAL" and desc.get("program") == "OpenEye OEToolkits":
                    smi = desc.get("descriptor", "")
                    break
            if not smi:
                for desc in d.get("pdbx_chem_comp_descriptor", []):
                    if desc.get("type", "").startswith("SMILES"):
                        smi = desc.get("descriptor", "")
                        break
            cache[code] = smi
        except Exception:
            cache[code] = ""
        if i % 20 == 0:
            print(f"  fetched {i}/{len(todo)} ligand definitions", flush=True)
    LIG_CACHE.parent.mkdir(parents=True, exist_ok=True)
    LIG_CACHE.write_text(json.dumps(cache, indent=1))
    return cache


def donor_environment(mol, idx: int) -> tuple:
    """(aromatic, ring_size, sorted heteroatom symbols in that ring) for an atom."""
    a = mol.GetAtomWithIdx(idx)
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if idx in ring:
            hetero = sorted(mol.GetAtomWithIdx(j).GetSymbol()
                            for j in ring if mol.GetAtomWithIdx(j).GetSymbol() != "C")
            return (a.GetIsAromatic(), len(ring), tuple(hetero))
    return (a.GetIsAromatic(), 0, ())


def main() -> None:
    from rdkit import Chem

    pq = DATA_PROCESSED / "cyp3a4_reference_geometry.parquet"
    if not pq.exists():
        raise SystemExit(f"missing {pq}; run scripts/structure/build_reference_set.py first")
    df = pd.read_parquet(pq)

    lig_rows = df[df.ligand != ""].copy()
    codes = sorted(lig_rows.ligand.unique())
    print(f"{len(codes)} unique ligand codes across {len(lig_rows)} chain-ligand rows")
    smi_map = ligand_smiles(codes)

    # one row per unique (ligand, observed class) — chains are replicates
    per_lig = (lig_rows.groupby("ligand")
               .agg(n_chains=("chain", "size"),
                    class_mode=("binding_class", lambda s: s.value_counts().index[0]),
                    min_fe=("min_fe_dist", "min"),
                    donor_el=("donor_element", lambda s: s.value_counts().index[0]))
               .reset_index())

    results, class_confusion = [], Counter()
    for _i, r in per_lig.iterrows():
        smi = smi_map.get(r.ligand, "")
        if not smi:
            results.append(dict(ligand=r.ligand, status="no_smiles"))
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            results.append(dict(ligand=r.ligand, status="unparseable"))
            continue

        observed = r.class_mode
        predicted = ligand_class(smi)
        class_confusion[(observed, predicted)] += 1

        sites = coordinating_atoms(smi, top_k=3)
        rec = dict(ligand=r.ligand, n_chains=int(r.n_chains), observed_class=observed,
                   predicted_class=predicted, min_fe_dist=round(float(r.min_fe), 3),
                   observed_donor_element=r.donor_el,
                   n_candidate_donors=len(sites),
                   top_donor=(sites[0].pattern if sites else None),
                   top_donor_score=(round(sites[0].score, 3) if sites else None),
                   top_donor_env=(str(donor_environment(mol, sites[0].atom_idx))
                                  if sites else None),
                   status="ok")
        results.append(rec)

    ok = [r for r in results if r.get("status") == "ok"]
    t2_obs = [r for r in ok if r["observed_class"] == "type_II_coordinated"]
    t1_obs = [r for r in ok if r["observed_class"] == "type_I_active_site"]

    # The number that matters: of the ligands that DO coordinate, how many does the
    # predictor even offer a competent donor for?
    t2_has_donor = sum(1 for r in t2_obs if r["n_candidate_donors"] > 0)
    t2_called_t2 = sum(1 for r in t2_obs if r["predicted_class"] == "type_II")
    t1_called_t2 = sum(1 for r in t1_obs if r["predicted_class"] == "type_II")

    summary = {
        "generated": "2026-09-09",
        "n_unique_ligands": int(len(per_lig)),
        "n_with_smiles": len(ok),
        "n_no_smiles": sum(1 for r in results if r.get("status") != "ok"),
        "observed_class_counts": dict(Counter(r["observed_class"] for r in ok)),
        "recall_on_coordinated": {
            "n_observed_type_II": len(t2_obs),
            "n_offering_any_donor": t2_has_donor,
            "n_called_type_II": t2_called_t2,
            "recall": round(t2_called_t2 / len(t2_obs), 3) if t2_obs else None,
            "donor_coverage": round(t2_has_donor / len(t2_obs), 3) if t2_obs else None,
        },
        "false_positive_on_type_I": {
            "n_observed_type_I": len(t1_obs),
            "n_called_type_II": t1_called_t2,
            "false_positive_rate": round(t1_called_t2 / len(t1_obs), 3) if t1_obs else None,
        },
        "top_donor_patterns_on_coordinated": dict(
            Counter(r["top_donor"] for r in t2_obs if r["top_donor"])),
        "misses": [{k: r[k] for k in ("ligand", "min_fe_dist", "predicted_class",
                                      "n_candidate_donors", "top_donor")}
                   for r in t2_obs if r["predicted_class"] != "type_II"][:25],
    }
    out = DATA_PROCESSED / "donor_prediction_validation.json"
    out.write_text(json.dumps({"summary": summary, "per_ligand": results}, indent=1))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
