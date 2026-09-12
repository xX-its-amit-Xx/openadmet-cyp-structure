"""Turn a scored pose pool into a submittable ZIP, and validate it before it is sent.

**Why this exists now, before the data has dropped.** The structure-track data is still
unreleased with 12 days to the interim deadline, and until this script existed there was no
tested path from "a pool of poses" to "a file the portal accepts". Every finding so far is
worth nothing without it, and the conversion has a trap in it that would surface at exactly
the wrong moment:

  **The validator requires a residue named literally `LIG`.** Boltz-2 emits the query ligand
  as `LIG1`, Chai-1 emits it as `LIG3`, and deposited structures use the PDB chemical
  component code. A submission built by copying co-folder output unchanged is rejected.

**Format — now the OFFICIAL spec**, read from the challenge Space's own `submission.py`
(`huggingface.co/spaces/openadmet/cyp-challenge`, fetched 2026-09-12), not inferred:

    "Submit a .zip archive containing exactly {STRUCTURE_DATASET_SIZE} .pdb files, one per
     compound, named after the compound identifier (e.g. x00011-1.pdb). Each file must be a
     full protein-ligand complex with the ligand residue named LIG."

Three things that settles:
  - one flat `.zip` of `<compound_id>.pdb` -- as built here;
  - the ligand residue must be named exactly `LIG` -- Boltz emits `LIG1`, Chai `LIG3`, so
    a zip of raw co-folder output is rejected. The converter renames and then asserts it;
  - **"full protein-ligand complex" resolves the heme question: KEEP it.** The default was
    already `--heme keep`; it is now the documented answer rather than a judgement call.

The Space also gates on FILE COUNT before anything else: a zip whose file count differs
from `STRUCTURE_DATASET_SIZE` is refused at upload with no further diagnosis. So
`--expect-n` checks that here, where the message can be useful.

⚠️ `STRUCTURE_TRACK_LIVE = False` and `STRUCTURE_DATASET_SIZE = 184  # TODO: Update when
final dataset is ready` as of 2026-09-12 -- 184 is the PXR count, a placeholder, and the
example id `x00011-1` is a PXR id too. Re-read both the moment the track goes live.

    python scripts/submit/build_submission.py build --tag val87b --out submissions/01_test.zip
    python scripts/submit/build_submission.py validate --zip submissions/01_test.zip
"""
from __future__ import annotations

import argparse
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED, SUBMISSIONS  # noqa: E402

LIG_RESNAME = "LIG"


def choose_poses(tag: str, arm: str = "unsteered") -> dict[str, str]:
    """Pick one pose per ligand with the FINDING 003 selector. Returns {ligand: sample}.

    score = 0.5 * z(pocket contacts) - z(mean RMSD to siblings), z-scored WITHIN ligand.
    Measured +0.0220 on held-out scaffold clusters, above the 99th percentile of the
    random-feature null (FINDING 007). Falls back to the medoid alone if the orientation
    features are missing, and says so rather than silently using a different rule.
    """
    scored = pd.read_csv(DATA_PROCESSED / f"poses_scored_{tag}.csv")
    scored = scored[scored.arm == arm]
    o = DATA_PROCESSED / f"orientation_features_{tag}_{arm}.csv"
    c = DATA_PROCESSED / f"consensus_features_{tag}_{arm}.csv"
    if not c.exists():
        raise SystemExit(f"missing {c}; run test_consensus_selector.py first")
    df = scored.merge(pd.read_csv(c), on=["ligand", "sample"], how="left")
    have_orient = o.exists()
    if have_orient:
        df = df.merge(pd.read_csv(o), on=["ligand", "sample"], how="left")
    df = df.dropna(subset=["mean_rmsd_to_others"])

    def zw(col):
        g = df.groupby("ligand")[col]
        return (df[col] - g.transform("mean")) / (g.transform("std") + 1e-9)

    if have_orient and "n_pocket_residues_touched" in df:
        df = df.dropna(subset=["n_pocket_residues_touched"])
        df["_s"] = 0.5 * zw("n_pocket_residues_touched") - zw("mean_rmsd_to_others")
        rule = "FINDING 003 (contacts + medoid)"
    else:
        df["_s"] = -zw("mean_rmsd_to_others")
        rule = "MEDOID ONLY - orientation features absent, weaker selector"
    print(f"selector: {rule}", flush=True)
    return {lig: g.loc[g._s.idxmax(), "sample"] for lig, g in df.groupby("ligand")}


def to_submission_pdb(cif_path: Path, out_pdb: Path, heme: str = "keep") -> dict:
    """Write a single-model PDB with the query ligand renamed to `LIG`.

    Returns a small report so the caller can assert on it rather than trust it.
    """
    import gemmi

    from cypstruct import pose as P

    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    while len(st) > 1:
        del st[1]                      # one model only

    n_lig = n_heme = 0
    for chain in st[0]:
        for res in list(chain):
            rname = res.name.strip().upper()
            info = gemmi.find_tabulated_residue(rname)
            if info and info.is_amino_acid():
                continue
            if P._looks_like_heme(res):
                n_heme += 1
                if heme == "drop":
                    chain.remove_residue(res)
                elif heme == "rename":
                    res.name = "HEM"
                continue
            # everything else that is not solvent is the query ligand
            if rname in ("HOH", "DOD"):
                chain.remove_residue(res)
                continue
            res.name = LIG_RESNAME
            n_lig += 1

    st.setup_entities()
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(out_pdb))
    return {"n_lig_residues": n_lig, "n_heme": n_heme, "out": str(out_pdb)}


def build(tag: str, arm: str, out_zip: Path, heme: str, pool_dir: Path | None,
          profile: str | None) -> None:
    import os
    import shutil
    import tempfile

    if profile:
        os.environ["MODAL_PROFILE"] = profile
    import modal

    picks = choose_poses(tag, arm)
    print(f"{len(picks)} ligands selected", flush=True)
    vol = modal.Volume.from_name("cyp-pool")
    tmp = Path(tempfile.mkdtemp(prefix="cypsub_", dir="D:/cyp_scratch"))
    reports = {}
    try:
        staged = tmp / "pdb"
        staged.mkdir(parents=True, exist_ok=True)
        for n, (lig, sample) in enumerate(sorted(picks.items()), 1):
            job = f"{lig}__{arm}__s1"
            src = tmp / f"{lig}.cif"
            try:
                src.write_bytes(b"".join(vol.read_file(f"/{tag}/{job}/{sample}.cif")))
            except Exception as exc:
                reports[lig] = {"error": f"{type(exc).__name__}"}
                continue
            try:
                reports[lig] = to_submission_pdb(src, staged / f"{lig}.pdb", heme)
            except Exception as exc:
                reports[lig] = {"error": f"convert: {type(exc).__name__}: {exc}"}
            src.unlink(missing_ok=True)
            if n % 20 == 0:
                print(f"  {n}/{len(picks)}", flush=True)

        pdbs = sorted(staged.glob("*.pdb"))
        out_zip.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for p in pdbs:
                zf.write(p, arcname=p.name)
        bad = {k: v for k, v in reports.items() if "error" in v}
        multi = {k: v for k, v in reports.items()
                 if "error" not in v and v["n_lig_residues"] != 1}
        print(f"\nwrote {len(pdbs)} PDBs -> {out_zip}")
        print(f"  errors: {len(bad)}   files without exactly one LIG residue: {len(multi)}")
        if bad:
            print("   first errors:", list(bad.items())[:3])
        if multi:
            print("   first multi/zero-LIG:", list(multi.items())[:3])
        (DATA_PROCESSED / f"submission_report_{tag}.json").write_text(
            json.dumps({"tag": tag, "heme": heme, "n": len(pdbs),
                        "reports": reports}, indent=1))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def validate(zip_path: Path, ligands_csv: Path, expect_n: int | None = None) -> int:
    """The checks the PXR validator ran, re-implemented so failures are found here.

    Deliberately strict about the ligand graph: a submission whose atoms do not match the
    requested SMILES is the failure mode that looks fine in a viewer and scores zero.
    """
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    df = pd.read_csv(ligands_csv)
    id_col = "structure" if "structure" in df.columns else "id"
    expected = dict(zip(df[id_col].astype(str), df["smiles"].astype(str)))

    errs: list[str] = []
    if not zip_path.exists():
        return _report([f"missing {zip_path}"])
    with zipfile.ZipFile(zip_path) as zf:
        names = [n for n in zf.namelist() if n.lower().endswith(".pdb")]
        got = {Path(n).stem for n in names}
        # The Space refuses on file count BEFORE any other check, with no diagnosis.
        if expect_n is not None and len(names) != expect_n:
            errs.append(f"zip has {len(names)} .pdb files, the portal expects {expect_n}")
        missing, extra = sorted(set(expected) - got), sorted(got - set(expected))
        if missing:
            errs.append(f"missing {len(missing)} structures, first: {missing[:5]}")
        if extra:
            errs.append(f"{len(extra)} unexpected structures, first: {extra[:5]}")

        import tempfile
        with tempfile.TemporaryDirectory(dir="D:/cyp_scratch") as td:
            for n in names:
                sid = Path(n).stem
                p = Path(td) / f"{sid}.pdb"
                p.write_bytes(zf.read(n))
                txt = p.read_text(errors="replace")
                lig_lines = [ln for ln in txt.splitlines()
                             if ln.startswith(("ATOM", "HETATM"))
                             and ln[17:20].strip().upper() == LIG_RESNAME]
                if not lig_lines:
                    errs.append(f"{sid}: no residue named {LIG_RESNAME}")
                    continue
                resids = {ln[22:27] for ln in lig_lines}
                if len(resids) != 1:
                    errs.append(f"{sid}: {len(resids)} {LIG_RESNAME} residues, expected 1")
                smi = expected.get(sid)
                if smi:
                    ref = Chem.MolFromSmiles(smi)
                    if ref is not None:
                        # Count heavy atoms by ATOMIC NUMBER, not via RemoveHs.
                        # Some challenge SMILES carry explicit hydrogens - metformin
                        # (MF8) is written [H]/N=C(/N)\N/C(=N\[H])/N(C)C - and RDKit
                        # deliberately KEEPS those, because they define the double-bond
                        # stereochemistry, so RemoveHs leaves the count at 11 for a
                        # 9-heavy-atom molecule. Counting Z > 1 is unambiguous and has
                        # no such exceptions. Without it the validator rejects a good
                        # pose, which is the worst false alarm available: it reads as a
                        # modelling failure and sends you to debug the co-folder.
                        n_heavy = sum(1 for at in ref.GetAtoms()
                                      if at.GetAtomicNum() > 1)
                        if len(lig_lines) != n_heavy:
                            errs.append(f"{sid}: ligand has {len(lig_lines)} atoms, "
                                        f"SMILES expects {n_heavy}")
    return _report(errs)


def _report(errs: list[str]) -> int:
    if not errs:
        print("VALID: all checks passed")
        return 0
    print(f"INVALID: {len(errs)} problem(s)")
    for e in errs[:20]:
        print("   -", e)
    return 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "validate"])
    ap.add_argument("--tag", default="val87b")
    ap.add_argument("--arm", default="unsteered")
    ap.add_argument("--out", default=str(SUBMISSIONS / "00_test.zip"))
    ap.add_argument("--zip", default=None)
    ap.add_argument("--ligands", default=str(DATA_PROCESSED / "validation_ligands.csv"))
    ap.add_argument("--heme", default="keep", choices=["keep", "drop", "rename"])
    ap.add_argument("--expect-n", type=int, default=None,
                    help="file count the portal expects (config.STRUCTURE_DATASET_SIZE)")
    ap.add_argument("--profile", default=None)
    a = ap.parse_args()
    if a.cmd == "build":
        build(a.tag, a.arm, Path(a.out), a.heme, None, a.profile)
    else:
        raise SystemExit(validate(Path(a.zip or a.out), Path(a.ligands), a.expect_n))
