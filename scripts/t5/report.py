"""Emit the T5 acceptance-gate report in the schema the BU spec asks for.

T5 wants three things before a production run is authorised:

  1. pose validity checks
  2. completion rate  (its words: "critical, as co-folders fail entirely on flexible
     pockets rather than producing poor poses")
  3. whether derived features reduce downstream RMSE

This script delivers 1 and 2 and the manifest. Item 3 is T6's and needs the activity
model, so it is reported as NOT ANSWERED rather than guessed at - a gate answered with two
of three criteria is not a passed gate.

On validity: PoseBusters is not installed here, and for this target that is not the loss it
sounds like. `qmscore.geometry.hard_validity` encodes physical impossibility specific to
CYPs - ligand density on the PROXIMAL heme face (the iron's sixth site is the only one
free; the fifth is Cys442's thiolate), steric overlap, sub-bond-length Fe contacts, and
poses touching no protein at all. PoseBusters checks none of those, because it does not
know what a heme is. Both are reported where they overlap; neither is claimed as the other.

Output schema, exactly as T5 specifies:
    inchikey_block, isoform, method, structure_path, confidence,
    predicted_affinity, completed

    python scripts/t5/report.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "cofold"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

OUT = DATA_PROCESSED / "t5_bakeoff"
METHOD = "protenix_v2@openprotein"


def main() -> int:
    from cypstruct import pose as P
    from cypstruct.qmscore import geometry as G

    st = json.loads((OUT / "state.json").read_text())
    compounds = pd.read_csv(st["compounds"]).set_index("inchikey_block")
    rows = []
    for job_id, rec in st["folds"].items():
        ik, iso = rec["inchikey_block"], rec["isoform"]
        completed = rec.get("done") is True
        path, conf, valid, why = "", None, None, ""
        geo_fields = {}
        if completed:
            files = rec.get("files") or []
            if files:
                p = OUT / "poses" / iso / files[0]
                path = str(p.relative_to(REPO)) if p.exists() else ""
                try:
                    mo = P.load_structure(p)
                    if len(mo.lig_xyz):
                        g = G.compute(mo.lig_xyz, mo.lig_elem, mo.prot_xyz,
                                      mo.heme_xyz, mo.heme_atom, mo.axial_sg)
                        valid, why = G.hard_validity(g)
                        geo_fields = {
                            "fe_min_dist": round(float(g.fe_min_dist), 3)
                            if g.fe_min_dist == g.fe_min_dist else None,
                            "frac_proximal": round(float(g.frac_atoms_proximal), 4)
                            if g.frac_atoms_proximal == g.frac_atoms_proximal else None,
                            "max_clash": round(float(g.max_clash), 3)
                            if g.max_clash == g.max_clash else None,
                            "n_contacts": int(g.n_protein_contacts),
                            "is_coordinated": bool(g.is_coordinated),
                        }
                    else:
                        valid, why = False, "no ligand atoms parsed"
                except Exception as exc:
                    valid, why = False, f"{type(exc).__name__}: {exc}"[:120]
        rows.append({
            "inchikey_block": ik,
            "isoform": iso,
            "method": METHOD,
            "structure_path": path,
            "confidence": conf,          # this engine returns none per-pose here
            "predicted_affinity": None,  # not produced by co-folding; T6/activity model
            "completed": completed,
            "valid": valid,
            "invalid_reason": why,
            **geo_fields,
        })

    df = pd.DataFrame(rows)
    man = OUT / "t5_manifest.csv"
    df.to_csv(man, index=False)

    done = df[df.completed]
    rep = {
        "compounds": int(compounds.shape[0]),
        "isoforms": sorted(df.isoform.unique()),
        "method": METHOD,
        "attempted": int(len(df)),
        "completed": int(df.completed.sum()),
        "completion_rate": round(float(df.completed.mean()), 4),
        "by_isoform": {
            iso: {
                "attempted": int(len(g)),
                "completion_rate": round(float(g.completed.mean()), 4),
                "valid_rate": (round(float(g[g.completed].valid.mean()), 4)
                               if g.completed.any() else None),
                "coordinating_rate": (round(float(g[g.completed].is_coordinated.mean()), 4)
                                      if "is_coordinated" in g and g.completed.any() else None),
            } for iso, g in df.groupby("isoform")
        },
        "validity_failures": (done[~done.valid.astype(bool)].invalid_reason
                              .str.replace(r"[\d.]+", "N", regex=True)
                              .value_counts().head(6).to_dict() if len(done) else {}),
        "gate_item_3_downstream_rmse": "NOT ANSWERED - requires the activity model (T6)",
        "manifest": str(man.relative_to(REPO)),
    }
    (OUT / "t5_gate_report.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
