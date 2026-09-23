"""FINDING 039 / playbook C5 - aromatase scores 0.863 in FINDING 022 and 0.631 in 024.

Two statistics on the same six pairs, not a disagreement: 022 reports the mean of
`lddt_sample0` (the engine's own rank-0 pose), 024 reports the mean over all five samples.
Checked on five targets and both aggregate rows.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
DP = REPO / "data" / "processed"

a = pd.read_csv(DP / "finetune" / "arm4_mix.csv")
sc = json.loads((DP / "finetune" / "scores_arm4_mix_base.json").read_text())
print("FILTER arm4_mix.csv: %d rows; split %s" % (len(a), a.split.value_counts().to_dict()))
te = a[a.split == "test"].copy()
miss = te[~te.pair.isin(sc)]
print("FILTER test split: %d rows -> %d scored; %d unscored: %s"
      % (len(te), te.pair.isin(sc).sum(), len(miss),
         miss[["pair", "uniprot"]].to_dict("records")))
te = te[te.pair.isin(sc)].copy()
te["mean5"] = [float(np.mean(sc[p]["lddt_pli"])) for p in te.pair]
te["s0"] = [float(sc[p]["lddt_sample0"]) for p in te.pair]

g = te.groupby("uniprot").agg(n=("pair", "size"), s0=("s0", "mean"), mean5=("mean5", "mean"))
out = {}
print("\n%-8s %3s  %-22s %-12s  %-22s %-12s" % ("uniprot", "n", "mean of sample0",
                                                "F022 quotes", "mean of all 5", "F024 quotes"))
for u, p22, p24 in [("P08684", 0.555, None), ("P11511", 0.863, 0.631),
                    ("P10614", 0.909, None), ("Q2IU02", 0.931, 0.924),
                    ("P20815", None, 0.506)]:
    if u not in g.index:
        continue
    r = g.loc[u]
    print("%-8s %3d  %-22.4f %-12s  %-22.4f %-12s" % (u, r.n, r.s0, p22, r.mean5, p24))
    out[u] = dict(n=int(r.n), mean_of_sample0=round(float(r.s0), 4),
                  mean_of_all_5=round(float(r.mean5), 4),
                  finding_022=p22, finding_024=p24,
                  gap=round(float(r.s0 - r.mean5), 4))

n3 = te[te.uniprot == "P08684"]
oth = te[te.uniprot != "P08684"]
print("\nCYP3A4  n=%d  sample0 %.4f  all-5 %.4f   (FINDING 022 prints 0.5550)"
      % (len(n3), n3.s0.mean(), n3.mean5.mean()))
print("others  n=%d  sample0 %.4f  all-5 %.4f   (FINDING 022 prints 0.8648)"
      % (len(oth), oth.s0.mean(), oth.mean5.mean()))
out["_cyp3a4"] = dict(n=len(n3), s0=round(float(n3.s0.mean()), 4),
                      mean5=round(float(n3.mean5.mean()), 4))
out["_others"] = dict(n=len(oth), s0=round(float(oth.s0.mean()), 4),
                      mean5=round(float(oth.mean5.mean()), 4))
out["_cavity_entry_counts_406set"] = {
    k: int(v) for k, v in a.uniprot.value_counts().reindex(
        ["P08684", "P11511", "P10614", "Q2IU02", "P20815"]).items()}
print("\n024's cavity-table denominators are crystals in the 406-entry set, not scored pairs:")
print("  ", out["_cavity_entry_counts_406set"])

(DP / "corrections_c5_aromatase.json").write_text(json.dumps(out, indent=2))
print("\nwrote data/processed/corrections_c5_aromatase.json")
