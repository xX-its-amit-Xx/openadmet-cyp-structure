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

    from cypstruct.pose import (
        Complex,
        best_ligand_mapping,
        bisy_rmsd,
        lddt_pli,
        load_structure,
    )

    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="arm4_mix")
    ap.add_argument("--tag", default="base", help="base | ft")
    ap.add_argument("--min-identity", type=float, default=0.80,
                    help="flag a pair whose best offset still disagrees this much")
    a = ap.parse_args()

    def shift(cx, k):
        """Renumber a Complex's protein residues by k, coordinates untouched."""
        return Complex(
            name=cx.name, prot_xyz=cx.prot_xyz,
            prot_key=[(c, r + k, at) for (c, r, at) in cx.prot_key],
            prot_res={(c, r + k): n for (c, r), n in cx.prot_res.items()},
            lig_xyz=cx.lig_xyz, lig_elem=cx.lig_elem, lig_name=cx.lig_name,
            lig_chain=cx.lig_chain, fe=cx.fe, heme_xyz=cx.heme_xyz,
            heme_atom=cx.heme_atom, heme_elem=cx.heme_elem, axial_sg=cx.axial_sg)

    def best_offset(mdl, ref):
        """The renumbering that maximises residue-NAME agreement with the reference.

        THIS IS NOT COSMETIC. lddt_pli pairs protein atoms by residue number. Boltz
        numbers its output 1..N from the input sequence; a crystal uses auth numbering,
        which for CYP3A4 starts near 29. Where those disagree, every contact is compared
        against the WRONG residue and the score collapses to exactly 0.0 - silently, and
        only for the targets whose offset happens to be large.

        Measured on this set: 3NA0_2DC went 0.000 -> 0.986 at offset +43, 3DSJ_243
        0.000 -> 0.941 at +27, 3NXU_RIT 0.003 -> 0.755 at +22. Pairs that were already
        aligned come back at offset 0 and are unchanged, which is the control.

        Matching on residue names rather than coordinates keeps this superposition-free
        and cannot manufacture agreement: a wrong offset scores near-zero identity.
        """
        rres = {r: n for (c, r), n in ref.prot_res.items()}
        mres = {r: n for (c, r), n in mdl.prot_res.items()}
        best, bk = -1, 0
        for k in range(-80, 81):
            agree = sum(1 for r, n in mres.items() if rres.get(r + k) == n)
            if agree > best:
                best, bk = agree, k
        return bk, (best / max(len(mres), 1))

    # boltz names its results directory after the INPUT yaml directory's basename, which
    # is not the arm name once variants exist (arm4_mix_bonded -> boltz_results_arm4_mix_bonded).
    # Globbing it rather than reconstructing the name avoids a scorer that silently
    # reports "no predictions" for 420 files that are sitting right there.
    root = ROOT / "holdout_out" / f"{a.arm}_{a.tag}"
    cands = sorted(root.glob("boltz_results_*/predictions"))
    if not cands:
        msg = f"no predictions under {root} (looked for boltz_results_*/predictions)"
        raise SystemExit(msg)
    if len(cands) > 1:
        msg = f"ambiguous: {len(cands)} results dirs under {root}: {[str(c) for c in cands]}"
        raise SystemExit(msg)
    preds = cands[0]

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
        offset = ident = None
        for cif in sorted(d.glob("*_model_*.cif")):
            try:
                mdl = load_structure(cif)
                if offset is None:
                    # One offset per pair, from the first model: the numbering is a
                    # property of the input sequence, not of the diffusion sample.
                    offset, ident = best_offset(mdl, ref)
                mdl = shift(mdl, offset)
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
            # Recorded so a bad alignment is visible in the output rather than being
            # absorbed into a low score.
            "resnum_offset": offset, "seq_identity_at_offset": round(ident or 0.0, 3),
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
    low = [k for k, v in out.items() if v["seq_identity_at_offset"] < a.min_identity]
    vals0 = [v["lddt_sample0"] for v in out.values()]
    valsb = [v["lddt_best"] for v in out.values()]
    print(json.dumps({
        "arm": a.arm, "tag": a.tag, "pairs_scored": done, "pairs_failed": failed,
        "mean_lddt_sample0": round(st.mean(vals0), 4) if vals0 else None,
        "mean_lddt_best_of_n": round(st.mean(valsb), 4) if valsb else None,
        "pairs_below_min_identity": len(low),
        "low_identity_ids": sorted(low)[:10],
        "median_resnum_offset": st.median([v["resnum_offset"] for v in out.values()]) if out else None,
        "reasons": dict(sorted(reasons.items(), key=lambda kv: -kv[1])[:4]),
        "out": str(dst),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
