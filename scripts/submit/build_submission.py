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

from cypstruct.paths import (  # noqa: E402
    DATA_PROCESSED, SCRATCH, SUBMISSIONS, guard_scratch)

LIG_RESNAME = "LIG"


def choose_poses(tag: str, arm: str = "unsteered") -> dict[str, str]:
    """Pick one pose per ligand. Returns {ligand: sample}.

    **Preferred: cross-engine agreement** (FINDING 011/012). Score each pose by its mean
    Chamfer distance, in the heme frame, to independent engines' poses of the SAME ligand;
    take the smallest. Measured **+0.0381** on CYP3A4 against the FINDING 003 selector's
    +0.0265 and Boltz confidence's -0.0417, and it **generalises**: +0.3006 across 17
    held-out P450 targets with no parameters fitted (FINDING 012).

    It needs a `xeng_{tag}.csv` of {ligand, sample, xeng} built from >= 4 GENUINELY
    independent reference poses - replicate jobs, deduplicated, because `diffusion_samples`
    does not diversify the ligand on OpenProtein (FINDING 009). Below that depth the
    feature is worthless: at one reference pose it measured -0.0055.

    Falls back to FINDING 003 (contacts + medoid, +0.0265) and then to the medoid alone,
    and always prints which rule ran rather than silently changing selectors.
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

    xf = DATA_PROCESSED / f"xeng_{tag}.csv"
    if xf.exists():
        xe = pd.read_csv(xf)
        df = df.merge(xe[["ligand", "sample", "xeng"]], on=["ligand", "sample"],
                      how="left")
        n_missing = int(df.xeng.isna().sum())
        if n_missing:
            print(f"  WARNING: {n_missing} poses have no cross-engine score", flush=True)
        df = df.dropna(subset=["xeng"])
        # unweighted and alone on purpose: adding the sibling-consensus term measured
        # WORSE at full depth (+0.0183 against +0.0336), and fitting weights collapses
        # it into the noise (+0.0149). There is nothing here to tune.
        df["_s"] = -zw("xeng")
        rule = "FINDING 011/012 cross-engine agreement (+0.0381)"
    elif have_orient and "n_pocket_residues_touched" in df:
        df = df.dropna(subset=["n_pocket_residues_touched"])
        df["_s"] = 0.5 * zw("n_pocket_residues_touched") - zw("mean_rmsd_to_others")
        rule = "FINDING 003 (contacts + medoid, +0.0265) - no xeng_*.csv found"
    else:
        df["_s"] = -zw("mean_rmsd_to_others")
        rule = "MEDOID ONLY - weakest selector; xeng and orientation features both absent"
    print(f"selector: {rule}", flush=True)
    return {lig: g.loc[g._s.idxmax(), "sample"] for lig, g in df.groupby("ligand")}


def choose_poses_blind(tag: str, ligands_csv: Path | None = None,
                       arm: str = "unsteered",
                       xeng_csv: Path | None = None,
                       fallback: bool = True) -> tuple[dict, dict]:
    """Pick one pose per ligand with NOTHING from a crystal. Returns (picks, report).

    `choose_poses` cannot run on a blind set for two independent reasons (playbook G2/G5):
    it reads `poses_scored_<tag>.csv`, which is produced by scoring against deposited
    structures, and it raises `SystemExit` without `consensus_features_<tag>_<arm>.csv`,
    which is produced by a Modal-only script. Neither exists on drop day. This is the same
    selection rule with the crystal-derived joins removed.

    The rule is unchanged and untuned: `cypstruct.xengine.select()`, i.e. `-z(xeng)` within
    the ligand, which is `argmin(xeng)` exactly. Board: selected 0.6164, oracle 0.6975,
    random 0.5769, gain **+0.0395** (FINDINGs 011/035/036/037).

    **Thin ligands are named, never silently degraded.** A ligand missing from
    `xeng_<tag>.csv` is one `build_xeng_feature --skip-thin` dropped for want of 4
    independent reference poses. Those fall back to FINDING 003 (contacts + medoid,
    +0.0265) *if* its two prediction-side feature files exist, and the report says exactly
    which ligands took which rule - because the difference is 0.013 LDDT-PLI per ligand and
    the old code path printed nothing at all.
    """
    import sys as _sys

    _sys.path.insert(0, str(REPO / "src"))
    from cypstruct import xengine as X

    xf = Path(xeng_csv) if xeng_csv else DATA_PROCESSED / f"xeng_{tag}.csv"
    if not xf.exists():
        raise SystemExit(
            f"missing {xf}. On a blind set the cross-engine feature is the selector; "
            "build it with `build_xeng_feature.py --tag "
            f"{tag} --pool <flat pose dir> --pool-flat --ref-npz <reference .npz>` "
            "after generating reference poses (`openprotein_cofold.py submit --sweep ...` "
            "then `refset`). Without it there is no blind selection rule at all.")
    xe = pd.read_csv(xf).dropna(subset=["xeng"])
    picks_df = X.select(xe)
    picks = dict(zip(picks_df.ligand.astype(str), picks_df["sample"].astype(str)))

    # the invariant the playbook asserts in three findings; cheap, so assert it here too
    argmin = xe.loc[xe.groupby("ligand").xeng.idxmin()]
    argmin_map = dict(zip(argmin.ligand.astype(str), argmin["sample"].astype(str)))
    if argmin_map != picks:
        differ = sorted(k for k in picks if argmin_map.get(k) != picks[k])
        raise SystemExit(f"argmin(xeng) != xengine.select() on {len(differ)} ligands "
                         f"({differ[:5]}) - the selector has changed shape, stop.")

    want = None
    if ligands_csv is not None:
        d = pd.read_csv(ligands_csv)
        idc = "structure" if "structure" in d.columns and "id" not in d.columns else "id"
        want = [str(v) for v in d[idc]]

    report = {"tag": tag, "xeng_csv": str(xf),
              "rule": "FINDING 011/012 cross-engine agreement (+0.0395) via "
                      "cypstruct.xengine.select()",
              "argmin_equals_select": True,
              "n_xeng_ligands": int(xe.ligand.nunique()),
              "n_xeng_poses": int(len(xe)),
              "selected_by_xeng": sorted(picks),
              "fell_back": [], "uncovered": []}

    if want is None:
        return picks, report

    thin = [lig for lig in want if lig not in picks]
    report["requested"] = len(want)
    if not thin:
        return picks, report

    # FINDING 003 fallback, both inputs prediction-side (no crystal):
    #   0.5 * z(n_pocket_residues_touched) - z(mean_rmsd_to_others)
    of = DATA_PROCESSED / f"orientation_features_{tag}_{arm}.csv"
    cf = DATA_PROCESSED / f"consensus_features_{tag}_{arm}.csv"
    if fallback and of.exists() and cf.exists():
        o = pd.read_csv(of)
        c = pd.read_csv(cf)
        f = o.merge(c, on=["ligand", "sample"], how="inner")
        f = f[f.ligand.astype(str).isin(thin)].dropna(
            subset=["n_pocket_residues_touched", "mean_rmsd_to_others"])
        if len(f):
            def zw(col):
                g = f.groupby("ligand")[col]
                return (f[col] - g.transform("mean")) / (g.transform("std") + 1e-9)
            f = f.assign(_s=0.5 * zw("n_pocket_residues_touched")
                         - zw("mean_rmsd_to_others"))
            got = f.loc[f.groupby("ligand")._s.idxmax()]
            for r in got.itertuples():
                picks[str(r.ligand)] = str(r.sample)
                report["fell_back"].append(str(r.ligand))
    report["fell_back"] = sorted(report["fell_back"])
    report["uncovered"] = sorted(lig for lig in want if lig not in picks)
    if report["fell_back"]:
        print(f"  {len(report['fell_back'])} ligand(s) fell back to FINDING 003 "
              f"(+0.0265, not +0.0395): {report['fell_back'][:10]}", flush=True)
    if report["uncovered"]:
        print(f"  !! {len(report['uncovered'])} ligand(s) have NO pose at all: "
              f"{report['uncovered'][:10]}", flush=True)
    return picks, report


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
          profile: str | None, blind: bool = False,
          ligands_csv: Path | None = None, pool_flat: bool = False) -> None:
    import os
    import shutil
    import tempfile

    blind_report = None
    if blind:
        picks, blind_report = choose_poses_blind(tag, ligands_csv, arm)
        print(f"selector: {blind_report['rule']}", flush=True)
        if pool_dir is None:
            raise SystemExit("--blind needs --pool-dir: there is no Modal volume for a "
                             "blind set and Modal is over its spend cap anyway")
    else:
        picks = choose_poses(tag, arm)
    print(f"{len(picks)} ligands selected", flush=True)

    # A LOCAL pool must work without Modal. `pool_dir` was a dead parameter - build()
    # always read from the Modal volume and was always called with None - so with Modal
    # over its spend cap the submission could not be built at all. That is a drop-day
    # blocker hiding in a code path nobody exercises until the day it matters.
    vol = None
    if pool_dir is None:
        if profile:
            os.environ["MODAL_PROFILE"] = profile
        import modal
        vol = modal.Volume.from_name("cyp-pool")
        print("pool source: Modal volume cyp-pool", flush=True)
    else:
        pool_dir = Path(pool_dir)
        print(f"pool source: local {pool_dir}", flush=True)
    # T9/G8: this used to be hardcoded to D:/cyp_scratch, on the drive with ~12 GB free
    # and the one the repo itself lives on. SCRATCH is C:/cyp_struct and guard_scratch
    # raises BEFORE the first write rather than wedging the box at 0 bytes.
    guard_scratch(2.0)
    tmp = Path(tempfile.mkdtemp(prefix="cypsub_", dir=str(SCRATCH)))
    reports = {}
    try:
        staged = tmp / "pdb"
        staged.mkdir(parents=True, exist_ok=True)
        for n, (lig, sample) in enumerate(sorted(picks.items()), 1):
            job = f"{lig}__{arm}__s1"
            src = tmp / f"{lig}.cif"
            try:
                if vol is not None:
                    src.write_bytes(b"".join(
                        vol.read_file(f"/{tag}/{job}/{sample}.cif")))
                else:
                    # accept both pool layouts: <pool>/<LIG>__<arm>__s1/<sample>.cif as
                    # the Modal volume is organised, and a flat <pool>/<sample>.cif
                    # Explorer's `collect` writes a FLAT directory whose file stems
                    # already carry the ligand (<lig>__s<seed>_m<k>.cif); the Modal
                    # volume nests by job. Both are accepted, and a flat pool whose
                    # stems are bare sample names still resolves by the second candidate.
                    cands = [pool_dir / job / f"{sample}.cif",
                             pool_dir / f"{sample}.cif",
                             pool_dir / lig / f"{sample}.cif",
                             pool_dir / f"{lig}__{sample}.cif"]
                    if pool_flat:
                        cands = [c for c in cands if c.parent == pool_dir]
                    hit = next((c for c in cands if c.exists()), None)
                    if hit is None:
                        raise FileNotFoundError(
                            f"{sample}.cif not found under {pool_dir}")
                    shutil.copyfile(hit, src)
            except Exception as exc:
                reports[lig] = {"error": f"{type(exc).__name__}: {exc}"}
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
                        "blind": bool(blind), "selection": blind_report,
                        "picks": picks, "reports": reports}, indent=1))
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
        guard_scratch(1.0)
        with tempfile.TemporaryDirectory(dir=str(SCRATCH)) as td:
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
    ap.add_argument("--pool-dir", default=None,
                    help="read poses from a local directory instead of the Modal volume")
    ap.add_argument("--blind", action="store_true",
                    help="select with cypstruct.xengine.select() over xeng_<tag>.csv "
                         "alone - no crystals, no Modal-derived feature files. This is "
                         "the drop-day path (playbook G5).")
    ap.add_argument("--pool-flat", action="store_true",
                    help="--pool-dir is a flat directory of <LIG>__<sample>.cif")
    a = ap.parse_args()
    if a.cmd == "build":
        build(a.tag, a.arm, Path(a.out), a.heme,
              Path(a.pool_dir) if a.pool_dir else None, a.profile,
              blind=a.blind, ligands_csv=Path(a.ligands) if a.blind else None,
              pool_flat=a.pool_flat)
    else:
        raise SystemExit(validate(Path(a.zip or a.out), Path(a.ligands), a.expect_n))
