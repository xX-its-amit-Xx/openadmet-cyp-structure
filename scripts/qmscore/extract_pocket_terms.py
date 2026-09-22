"""Compute `cypstruct.qmscore.pocket` terms for every pose in the run5 pool.

Reads only predictions. The crystal is never opened here, so nothing in the output can
be leaky by construction; truth is joined on afterwards by the evaluator.

    sbatch run_pocket.sh          # 1740 poses, ~87 ligands x 5 samples x 4 seeds
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path("/scratch/shenoy.am/cyp-finetune")
EXP = Path("/scratch/shenoy.am/zexp")
sys.path.insert(0, str(ROOT / "src"))


def read_canon(path: Path) -> str:
    seq, keep = [], False
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            keep = "cyp3a4" in line.lower()
            continue
        if keep:
            seq.append(line.strip())
    return "".join(seq)


def main() -> int:
    import pandas as pd

    from cypstruct.pose import load_structure
    from cypstruct.qmscore import pocket

    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", default="run5_s42,run5_s43,run5_s44,run5_s45")
    ap.add_argument("--fasta", default=str(EXP / "cyp_sequences.fasta"))
    ap.add_argument("--out", default=str(EXP / "pocket_terms.csv"))
    a = ap.parse_args()

    canon = read_canon(Path(a.fasta))
    if len(canon) != 503:
        raise SystemExit(f"canonical CYP3A4 should be 503 aa, got {len(canon)}")
    man = {m["name"]: m for m in json.loads((EXP / "manifest.json").read_text())}

    rows, fails = [], {}
    t0 = time.time()
    for tag in a.tags.split(","):
        preds = sorted((EXP / "out" / tag).glob("boltz_results_*/predictions"))[0]
        for d in sorted(preds.iterdir()):
            if not d.is_dir():
                continue
            smiles = man.get(d.name, {}).get("smiles")
            offset = None
            for cif in sorted(d.glob("*_model_*.cif")):
                rank = int(cif.stem.rsplit("_", 1)[1])
                try:
                    cx = load_structure(cif)
                    if offset is None:
                        offset, ident = pocket.resolve_offset(cx.prot_res, canon)
                        if ident < 0.80:
                            raise ValueError(f"numbering identity {ident:.2f} too low")
                    t = pocket.compute(
                        cx.lig_xyz, cx.lig_elem, cx.prot_xyz, cx.prot_key, cx.prot_res,
                        canon, fe=cx.fe, heme_xyz=cx.heme_xyz, heme_atom=cx.heme_atom,
                        heme_elem=cx.heme_elem, axial_sg=cx.axial_sg, smiles=smiles,
                        offset=offset)
                except Exception as exc:                      # noqa: BLE001
                    k = f"{type(exc).__name__}: {str(exc)[:60]}"
                    fails[k] = fails.get(k, 0) + 1
                    continue
                rec = {"name": d.name, "tag": tag, "rank": rank}
                rec.update(t.to_dict())
                rows.append(rec)
        print(f"{tag}: {len(rows)} poses, {time.time() - t0:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(a.out, index=False)
    summary = {
        "poses": len(df), "ligands": int(df["name"].nunique()),
        "median_offset": float(df["resnum_offset"].median()),
        "min_seq_identity": round(float(df["seq_identity"].min()), 3),
        "anchors_present_median": float(df["n_anchors_present"].median()),
        "phe_present_median": float(df["n_phe_present"].median()),
        "strain_computed": int(df["strain"].notna().sum()),
        "rings_found_frac": round(float((df["phe_n_rings"] > 0).mean()), 3),
        "failures": fails, "out": a.out, "seconds": round(time.time() - t0),
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
