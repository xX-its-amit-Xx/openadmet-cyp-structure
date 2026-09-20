"""Score one checkpoint's held-out predictions against the crystals.

Runs ON EXPLORER, beside both the predictions and the reference mmCIFs, so nothing large
crosses the network. It imports `cypstruct.pose` from a copy of src/ pushed alongside -
the same code that scored every pool number in this repo, because a second implementation
of LDDT-PLI would make the fine-tune's number incomparable to everything it must beat.

Three things it does NOT take shortcuts on:

1. **Atom mapping.** `best_ligand_mapping` is called for every pose. Boltz and the PDB
   list ligand atoms in different orders; index-for-index gives 11.4 A where the truth is
   far smaller, silently, and halves LDDT-PLI.
2. **The reference chain.** Crystals here are often multi-chain and the arm CSV names
   which one the measured ligand sits in. Parsing all of them puts several copies of the
   protein in one Complex and makes every distance meaningless.
3. **Both per-sample and best-of-N.** Reporting only best-of-5 measures the pool's oracle,
   not what a submission would ship. Sample 0 is what you get with no selection at all.
   Both go in the output; the comparison script decides which to gate on.

    ./env/bin/python score_holdout.py --arm arm4_mix --tag base
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

ROOT = Path("/scratch/shenoy.am/cyp-finetune")
sys.path.insert(0, str(ROOT / "src"))


def main() -> int:  # noqa: PLR0915
    import pandas as pd

    from cypstruct.pose import best_ligand_mapping, bisy_rmsd, lddt_pli, load_structure

    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="arm4_mix")
    ap.add_argument("--tag", default="base", help="base | ft")
    a = ap.parse_args()

    preds = ROOT / "holdout_out" / f"{a.arm}_{a.tag}" / f"boltz_results_{a.arm}" / "predictions"
    if not preds.exists():
        msg = f"no predictions at {preds}"
        raise SystemExit(msg)

    df = pd.read_csv(ROOT / "finetune_arms" / f"{a.arm}.csv")
    rows = {f"{r.pdb}_{r.id}": r for r in df.itertuples()}

    out: dict[str, dict] = {}
    done = failed = 0
    reasons: dict[str, int] = {}
    for d in sorted(preds.iterdir()):
        if not d.is_dir() or d.name not in rows:
            continue
        r = rows[d.name]
        try:
            ref = load_structure(ROOT / "rcsb" / f"{r.pdb}.cif",
                                 ligand_code=str(r.id), assembly_chain=str(r.chain))
        except Exception as exc:  # noqa: BLE001
            failed += 1
            k = f"ref {type(exc).__name__}: {str(exc)[:60]}"
            reasons[k] = reasons.get(k, 0) + 1
            continue

        lddts, rmsds = [], []
        for cif in sorted(d.glob("*_model_*.cif")):
            try:
                mdl = load_structure(cif)
                perm = best_ligand_mapping(str(r.smiles), mdl, ref)
                lddts.append(float(lddt_pli(mdl, ref, lig_perm=perm)))
                rmsds.append(float(bisy_rmsd(mdl, ref, lig_perm=perm)))
            except Exception as exc:  # noqa: BLE001
                k = f"model {type(exc).__name__}: {str(exc)[:60]}"
                reasons[k] = reasons.get(k, 0) + 1
                if len(reasons) <= 2:
                    traceback.print_exc(limit=2)
        if not lddts:
            failed += 1
            continue
        out[d.name] = {
            "pdb": r.pdb, "ligand": str(r.id), "n_samples": len(lddts),
            "lddt_pli": lddts, "bisy_rmsd": rmsds,
            # sample_0 is "no selection"; best is the pool oracle. Never quote one alone.
            "lddt_sample0": lddts[0], "lddt_best": max(lddts),
            "rmsd_sample0": rmsds[0], "rmsd_best": min(rmsds),
        }
        done += 1
        if done % 20 == 0:
            print(f"  scored {done}", flush=True)

    dst = ROOT / "holdout_out" / f"scores_{a.arm}_{a.tag}.json"
    dst.write_text(json.dumps(out, indent=1))

    import statistics as st
    vals0 = [v["lddt_sample0"] for v in out.values()]
    valsb = [v["lddt_best"] for v in out.values()]
    print(json.dumps({
        "arm": a.arm, "tag": a.tag, "pairs_scored": done, "pairs_failed": failed,
        "mean_lddt_sample0": round(st.mean(vals0), 4) if vals0 else None,
        "mean_lddt_best_of_n": round(st.mean(valsb), 4) if valsb else None,
        "reasons": dict(sorted(reasons.items(), key=lambda kv: -kv[1])[:4]),
        "out": str(dst),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
