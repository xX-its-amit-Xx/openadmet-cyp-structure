"""Did the heme bond actually change the iron coordination?

Runs ON EXPLORER. The bonded holdout run came back a null (-0.0081 LDDT-PLI, p=0.086,
catastrophes 39 -> 38). That has two completely different explanations and they demand
opposite next steps:

  (a) the constraint never took effect - a broken experiment, worth fixing; or
  (b) it worked, the poses coordinate the iron properly, and LDDT-PLI did not care -
      in which case coordination is not the bottleneck and the whole anchor family of
      ideas is dead.

Measuring the Fe-to-ligand distance in both prediction sets separates them in one pass.
Reference distances from the crystals are printed alongside as the target: this repo's
deposited-structure survey puts Fe-donor at 1.94 / 2.20 / 2.38 A (p5/p50/p95).

    ./env/bin/python check_coordination.py --arm arm4_mix --tags base base_bonded
"""
from __future__ import annotations

import argparse
import json
import statistics as st
import sys
from pathlib import Path

ROOT = Path("/scratch/shenoy.am/cyp-finetune")
sys.path.insert(0, str(ROOT / "src"))


def fe_distances(cx) -> tuple[float, float]:
    """(nearest ligand atom to Fe, nearest ligand NITROGEN to Fe). inf when absent."""
    import numpy as np

    if cx.fe is None or len(cx.lig_xyz) == 0:
        return float("inf"), float("inf")
    d = np.linalg.norm(cx.lig_xyz - cx.fe, axis=1)
    any_d = float(d.min())
    ns = [i for i, e in enumerate(cx.lig_elem) if e.upper() == "N"]
    n_d = float(d[ns].min()) if ns else float("inf")
    return any_d, n_d


def summarise(name: str, vals: list[float]) -> dict:
    v = [x for x in vals if x != float("inf")]
    if not v:
        return {"set": name, "n": 0}
    v.sort()
    return {
        "set": name, "n": len(v),
        "p5": round(v[int(0.05 * len(v))], 3),
        "median": round(st.median(v), 3),
        "p95": round(v[min(int(0.95 * len(v)), len(v) - 1)], 3),
        # 1.90-2.45 A is the window this repo calibrated on deposited P450 structures.
        "frac_in_window": round(sum(1.90 <= x <= 2.45 for x in v) / len(v), 3),
        "frac_under_3A": round(sum(x < 3.0 for x in v) / len(v), 3),
    }


def main() -> int:
    import pandas as pd

    from cypstruct.pose import load_structure

    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="arm4_mix")
    ap.add_argument("--tags", nargs="+", default=["base", "base_bonded"])
    a = ap.parse_args()

    df = pd.read_csv(ROOT / "finetune_arms" / f"{a.arm}.csv")
    rows = {f"{r.pdb}_{r.id}": r for r in df.itertuples()}

    report = []
    crystal_any, crystal_n = [], []
    seen_crystal = set()

    for tag in a.tags:
        root = ROOT / "holdout_out" / f"{a.arm}_{tag}"
        cands = sorted(root.glob("boltz_results_*/predictions"))
        if not cands:
            print(f"skip {tag}: no predictions")
            continue
        anys, ns = [], []
        for d in sorted(cands[0].iterdir()):
            if not d.is_dir() or d.name not in rows:
                continue
            r = rows[d.name]
            # Only sample 0, so this measures what a submission would ship rather than
            # the best of five - the same convention as the LDDT-PLI comparison.
            cif = sorted(d.glob("*_model_0.cif"))
            if not cif:
                continue
            try:
                cx = load_structure(cif[0])
            except Exception:  # noqa: BLE001, S112
                continue
            x, y = fe_distances(cx)
            anys.append(x)
            ns.append(y)

            if d.name not in seen_crystal:
                seen_crystal.add(d.name)
                try:
                    ref = load_structure(ROOT / "rcsb" / f"{r.pdb}.cif",
                                         ligand_code=str(r.id), assembly_chain=str(r.chain))
                    cx_, cn_ = fe_distances(ref)
                    crystal_any.append(cx_)
                    crystal_n.append(cn_)
                except Exception:  # noqa: BLE001, S112
                    pass
        report.append(summarise(f"{tag} (Fe-any)", anys))
        report.append(summarise(f"{tag} (Fe-N)", ns))

    report.append(summarise("CRYSTAL (Fe-any)", crystal_any))
    report.append(summarise("CRYSTAL (Fe-N)", crystal_n))
    print(json.dumps(report, indent=1))
    (ROOT / "holdout_out" / f"coordination_{a.arm}.json").write_text(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
