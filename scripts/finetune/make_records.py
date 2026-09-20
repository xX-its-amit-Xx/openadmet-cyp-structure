"""Convert our crystal structures into Boltz training records.

Runs ON EXPLORER, against the staged inputs. Uses `boltz.data.parse.mmcif.parse_mmcif` -
Boltz's own parser, not a reimplementation - because training data parsed differently from
how inference parses it degrades a fine-tune silently instead of failing. That class of
mismatch is the most expensive kind of bug available here.

Emits, per arm:
    <arm>/structures/<pdb>.npz    StructureV2, dumped by boltz's own writer
    <arm>/manifest.json           boltz Manifest over every parsed structure
    <arm>/manifest_train.json     the same, restricted to the train split
    <arm>/manifest_test.json      ditto for held-out

The manifests are real `boltz.data.types.Manifest` objects, not a convenience dict: the
datamodule calls `Manifest.load` and then reads `record.chains[i].msa_id` to find the MSA
and `record.id` to find the structure. A hand-rolled schema gets as far as loading and
then trains on no MSA at all.

SPLITS ARE PER-PDB, NOT PER-PAIR. build_arms.py splits ligand pairs, but one npz holds
every ligand in that entry, so a PDB with one train pair and one test pair would put the
held-out ligand's coordinates into training. Any PDB touching the test split goes to test
whole. The count of pairs demoted this way is reported.

The split itself comes from build_arms.py, which holds out whole Murcko-scaffold clusters
and checks held-out novelty against the challenge's measured 0.587 median (FINDING 017).
Nothing here re-splits; it only carries that decision forward.

    ./env/bin/python make_records.py --arm arm4_mix
"""
from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path

ROOT = Path("/scratch/shenoy.am/cyp-finetune")
CCD = ROOT / "boltz_cache" / "ccd.pkl"
# mols.tar unpacks WITH its own `mols/` prefix, so `tar -C mols` once produced
# mols/mols/*.pkl. That was flattened; the files are directly in mols/ now. The stale
# `mols/mols` path survived the flattening as an EMPTY DIRECTORY, so `.exists()` stayed
# True and this script kept handing the parser a moldir with nothing in it - 28 ligands
# absent from ccd.pkl died with "expected str, bytes or os.PathLike object, not
# NoneType". An existence check on a directory is not a check that it has contents.
MOLDIR = ROOT / "boltz_cache" / "mols"
MSA_NPZ = ROOT / "msa_npz"


def load_mols():
    """The CCD component dictionary the parser needs to recognise HEM and the ligands."""
    import pickle
    if CCD.exists():
        with CCD.open("rb") as fh:
            return pickle.load(fh)  # noqa: S301
    return None


def build_chain_infos(structure, msa_id):
    """One ChainInfo per chain, with the MSA attached to protein chains only.

    msa_id is a *string key* for proteins and the sentinel -1 everywhere else; the loader
    branches on `msa_id != -1 and msa_id != ""`, so a ligand or heme chain carrying a
    protein MSA id would try to load one and crash. `const.chain_type_ids["PROTEIN"]` is
    the discriminator rather than a hardcoded 0.
    """
    from boltz.data import const
    from boltz.data.types import ChainInfo

    protein = const.chain_type_ids["PROTEIN"]
    infos = []
    for i, ch in enumerate(structure.chains):
        is_prot = int(ch["mol_type"]) == protein
        infos.append(ChainInfo(
            chain_id=i,
            chain_name=str(ch["name"]),
            mol_type=int(ch["mol_type"]),
            cluster_id=msa_id if is_prot else "",
            msa_id=(msa_id if (is_prot and msa_id) else -1),
            num_residues=int(ch["res_num"]),
            valid=True,
            entity_id=int(ch["entity_id"]),
        ))
    return infos


def main() -> int:  # noqa: PLR0915
    import pandas as pd
    from boltz.data.parse.mmcif import parse_mmcif
    from boltz.data.types import Manifest, Record

    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="arm4_mix")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    arm_csv = ROOT / "finetune_arms" / f"{a.arm}.csv"
    df = pd.read_csv(arm_csv)
    if a.limit:
        df = df.head(a.limit)
    out = ROOT / "records" / a.arm
    (out / "structures").mkdir(parents=True, exist_ok=True)

    # Collapse pairs to PDBs, and let test win any disagreement (see module docstring).
    by_pdb = {}
    for r in df.itertuples():
        e = by_pdb.setdefault(r.pdb, {"splits": set(), "target": r.target_key, "ligands": []})
        e["splits"].add(r.split)
        e["ligands"].append(r.id)
    demoted = sum(1 for e in by_pdb.values() if e["splits"] == {"train", "test"})

    mols = load_mols()
    moldir = str(MOLDIR) if any(MOLDIR.glob("*.pkl")) else None
    print(
        "{}: {} pairs -> {} pdbs ({} mixed-split, demoted to test) mols={} moldir={} msa_npz={}".format(
            a.arm, len(df), len(by_pdb), demoted,
            "yes" if mols else "no", moldir, MSA_NPZ.exists()),
        flush=True)

    records, splits = [], {}
    ok = fail = no_msa = 0
    reasons = {}
    for pdb, e in sorted(by_pdb.items()):
        cif = ROOT / "rcsb" / f"{pdb}.cif"
        if not cif.exists():
            reasons["missing cif"] = reasons.get("missing cif", 0) + 1
            fail += 1
            continue
        try:
            # compute_interfaces=False is a WORKAROUND for an upstream bug, not a
            # preference. boltz/data/parse/mmcif.py defines compute_interfaces() as a
            # module function at line 319 AND takes a parameter of the same name; at
            # line 1206 the parameter shadows the function, so calling
            # compute_interfaces(atoms, chains) calls a bool and raises TypeError.
            # The default is True, so every caller using defaults hits it. False takes
            # the else-branch and returns an empty interface array.
            #
            # What that costs: interfaces is chain-chain contact metadata used to sample
            # interface-centred crops. Our training signal is the LIGAND pose, carried by
            # diffusion_loss_weight 4.0 on coordinates, and the cropper is handed the
            # ligand chain instead. Worth revisiting if interface-weighted loss is on.
            parsed = parse_mmcif(str(cif), mols=mols, moldir=moldir,
                                 use_assembly=False, compute_interfaces=False)
            parsed.data.dump(out / "structures" / "{}.npz".format(pdb))

            msa_id = str(e["target"])
            if not (MSA_NPZ / "{}.npz".format(msa_id)).exists():
                msa_id = ""
                no_msa += 1
            records.append(Record(
                id=pdb,
                structure=parsed.info,
                chains=build_chain_infos(parsed.data, msa_id),
                interfaces=[],
            ))
            splits[pdb] = "test" if "test" in e["splits"] else "train"
            ok += 1
        except Exception as exc:  # noqa: BLE001
            fail += 1
            k = "{}: {}".format(type(exc).__name__, str(exc)[:60])
            reasons[k] = reasons.get(k, 0) + 1
            if fail <= 3:
                traceback.print_exc(limit=2)
        if (ok + fail) % 50 == 0:
            print("  {}/{}  ok={} fail={}".format(ok + fail, len(by_pdb), ok, fail), flush=True)

    Manifest(records).dump(out / "manifest.json")
    train = [r for r in records if splits[r.id] == "train"]
    test = [r for r in records if splits[r.id] == "test"]
    Manifest(train).dump(out / "manifest_train.json")
    Manifest(test).dump(out / "manifest_test.json")

    print(json.dumps({
        "arm": a.arm, "pairs": len(df), "pdbs": len(by_pdb),
        "parsed": ok, "failed": fail,
        "mixed_split_demoted_to_test": demoted,
        "records_without_msa": no_msa,
        "train": len(train), "test": len(test),
        "failure_reasons": dict(sorted(reasons.items(), key=lambda kv: -kv[1])[:5]),
        "out": str(out),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
