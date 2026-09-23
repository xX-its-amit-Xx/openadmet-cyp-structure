"""Run `openprotein_cofold.refset` for the refvalue tag and SAVE the full report.

The CLI prints the report without `per_ligand`; FINDING 041 needs every filter count,
including the zeros, and the per-ligand depth table, on disk.
"""
import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "opc", REPO / "scripts" / "cofold" / "openprotein_cofold.py")
opc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(opc)

DATA = REPO / "data" / "processed"
md = int(sys.argv[1]) if len(sys.argv) > 1 else 8
allow = "--allow-thin" in sys.argv
r = opc.refset("refvalue", ["protenix_v2", "protenix"],
               DATA / "refvalue_reference_set.npz", min_depth=md,
               ligands_csv=str(DATA / "refvalue_ligands.csv"), allow_thin=allow)
(DATA / "refvalue_refset_report.json").write_text(json.dumps(r, indent=2, default=str))
print(json.dumps({k: v for k, v in r.items() if k != "per_ligand"}, indent=2, default=str))
