"""Within-ligand variability of every pocket term. FINDING 014's pre-check."""
import glob, sys
import pandas as pd, numpy as np
sys.path.insert(0,"/scratch/shenoy.am/cyp-finetune/src")
from cypstruct.qmscore.pocket import FEATURES, GROUPS
p=pd.concat([pd.read_csv(f) for f in sorted(glob.glob("/scratch/shenoy.am/zexp/pocket_terms_run5_s4*.csv"))])
g=p.groupby("name")
rows=[]
for f in FEATURES:
    v=p[f].astype(float)
    sd=g[f].apply(lambda s: np.nanstd(s.astype(float)))
    const=float((sd.fillna(0)<1e-9).mean())
    mu=g[f].apply(lambda s: np.nanmean(np.abs(s.astype(float))))
    cv=float(np.nanmedian((sd/(mu+1e-9)).values))
    grp=[k for k,c in GROUPS.items() if f in c][0]
    rows.append({"group":grp,"feature":f,"frac_const_within_ligand":round(const,3),
                 "median_within_CV":round(cv,3),"frac_missing":round(float(v.isna().mean()),3),
                 "mean":round(float(np.nanmean(v)),3)})
df=pd.DataFrame(rows)
pd.set_option("display.width",160)
print(df.to_string(index=False))
df.to_csv("/scratch/shenoy.am/zexp/pocket_within_cv.csv",index=False)
