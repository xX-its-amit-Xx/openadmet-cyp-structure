"""Render FINDING 041's tables from the artefacts. No computation lives here."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data" / "processed"


def main() -> int:
    c = json.loads((DATA / "refvalue_curve.json").read_text())
    ctl = json.loads((DATA / "refvalue_controls.json").read_text())
    ship = json.loads((DATA / "refvalue_shipped_ref_arm.json").read_text())
    rs = json.loads((DATA / "refvalue_refset_report.json").read_text())
    fl = c["noise_floor_this_population"]
    pool = c["pool"]

    print(f"n = {c['n_ligands']}  oracle = {pool['oracle']:.5f}  "
          f"random = {pool['random']:.5f}")
    print(f"floor p95 = {fl['p95']} (SE {fl['se_of_p95']}, {fl['draws']:,} draws), "
          f"p99 = {fl['p99']}")
    print()
    print("| depth | subsets/ligand | selected (exact) | gain over random | clears floor |"
          " sd across draws | p5 | p95 | mean within-ligand rho | correct sign |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for d in ("4", "5", "6", "7", "8"):
        v = c["curve"][d]
        print(f"| {d} | - | {v['selected_exact']:.5f} | {v['gain_exact']:+.5f} | "
              f"{'yes' if v['gain_exact'] > fl['p95'] else 'NO'} | "
              f"{v['draw_sd']:.5f} | {v['draw_p5']:.5f} | {v['draw_p95']:.5f} | "
              f"{v['mean_within_ligand_rho']:+.4f} | "
              f"{100*v['frac_correct_sign']:.1f}% |")
    print()
    p = c["primary"]
    print(json.dumps(p, indent=2, default=str)[:2000])
    print()
    print("shipped-reference arm:", json.dumps(
        {k: (round(v["selected"], 5), round(v["gain"], 5)) if isinstance(v, dict) else v
         for k, v in ship["curve"].items()}), "full:", ship["full"])
    print()
    print("controls:", json.dumps({k: v.get("pass") for k, v in ctl.items()
                                   if isinstance(v, dict)}))
    print("refset filters:", json.dumps(rs["filters"]))
    # per-ligand table
    rows = []
    for lig in c["per_ligand"]:
        r = {"ligand": lig, **c["per_ligand"][lig]}
        for d in ("4", "5", "6", "7", "8"):
            r[f"sel_d{d}"] = c["per_ligand_mean"][d][lig]
        r["delta_8_4"] = r["sel_d8"] - r["sel_d4"]
        r["subset_pair_disagreement"] = c["subset_pair_disagreement"][lig]
        rows.append(r)
    df = pd.DataFrame(rows).sort_values("delta_8_4")
    df.to_csv(DATA / "refvalue_per_ligand.csv", index=False)
    print(f"\nper-ligand -> refvalue_per_ligand.csv ({len(df)} rows)")
    print(df[["ligand", "ref_depth", "oracle", "random", "sel_d4", "sel_d8",
              "delta_8_4", "subset_pair_disagreement"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
