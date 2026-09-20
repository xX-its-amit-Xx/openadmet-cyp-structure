"""Write Boltz inference YAMLs for one arm's held-out pairs.

Runs ON EXPLORER. The gate for the fine-tune is held-out LDDT-PLI against the crystal,
measured by `cypstruct.pose` - not any number Lightning logs. That means predicting the
test pairs twice from identical inputs, once with the base checkpoint and once with the
fine-tuned one, and comparing the pairs.

The YAMLs are written ONCE and used for both arms of that comparison. Generating them
separately per checkpoint would let an input difference masquerade as a fine-tuning
effect, which is the single easiest way to fake a win here.

MSAs are local a3m paths, not `--use_msa_server`: Explorer's GPU nodes have no direct
internet, which cost the PXR campaign twelve days.

**No bond constraint in this version.** An explicit Cys-SG -> heme-FE bond is what makes
Boltz reach crystallographic Fe-donor geometry (FINDING 001 addendum), but the ligating
cysteine's index differs per target and is not in the arm CSVs. Its absence lowers the
ABSOLUTE scores for both checkpoints equally, so the paired delta stays valid; the
absolute numbers from this script should not be compared to pool scores that had it.

    ./env/bin/python build_holdout_yaml.py --arm arm4_mix
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path("/scratch/shenoy.am/cyp-finetune")


def main() -> int:
    import pandas as pd
    from boltz.data.types import Manifest

    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="arm4_mix")
    a = ap.parse_args()

    arm_dir = ROOT / "records" / a.arm
    # manifest_test.json is post-validation, so a pair whose structure does not featurize
    # is already gone. Reading the CSV's split column instead would silently re-admit it.
    test_ids = {r.id for r in Manifest.load(arm_dir / "manifest_test.json").records}

    df = pd.read_csv(ROOT / "finetune_arms" / f"{a.arm}.csv")
    df = df[(df["split"] == "test") & (df["pdb"].isin(test_ids))]

    out = ROOT / "holdout_yaml" / a.arm
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob("*.yaml"):
        old.unlink()

    written, no_msa = 0, 0
    for r in df.itertuples():
        a3m = ROOT / "msa" / f"{r.target_key}.a3m"
        if not a3m.exists():
            no_msa += 1
            continue
        # SMILES go in single quotes: they routinely contain ':' and '#', which bare YAML
        # reads as a mapping or a comment and silently truncates the molecule.
        smiles = str(r.smiles).replace("'", "")
        yaml = (
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: {r.sequence}\n"
            f"      msa: {a3m}\n"
            "  - ligand:\n"
            "      id: B\n"
            "      ccd: HEM\n"
            "  - ligand:\n"
            "      id: C\n"
            f"      smiles: '{smiles}'\n"
        )
        (out / f"{r.pdb}_{r.id}.yaml").write_text(yaml)
        written += 1

    print(json.dumps({
        "arm": a.arm, "test_pdbs": len(test_ids), "pairs_written": written,
        "skipped_no_msa": no_msa, "out": str(out),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
