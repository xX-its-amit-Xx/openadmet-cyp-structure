"""The ligand's INTERNAL CONFORMER as a scoring prior - torsions, not placement.

Pre-registered in `docs/PREREG_conformer_transfer.md`, committed at `f9d951e` before a
single number existed. Nothing here may be tuned after a score exists. Written up as
FINDING 032.

**The lever.** FINDING 030 refuted fragment pose transfer on ABSOLUTE PLACEMENT: the
superfamily's fragment centroids scatter 4.48 A about their own mean where the error to
fix is 2.5 A, so the prior is twice as wide as the thing it must resolve. FINDING 024
leaves **26% of CYP3A4's ligand error (1.60 A) in the ligand's own internal conformer**,
and a torsion is **frame-free** - it does not care where the pocket is or how the ligand
is turned - so 030's mechanism does not transfer to it. This is the largest untested term
in the decomposition.

Two arms, deliberately different in kind:

  A  donor torsion transfer - does the P450 crystallographic record say what torsion a
     shared rotatable bond should take? Circular statistics with an explicit symmetry
     period per torsion.
  B  prior-free internal plausibility - forget P450s: is this conformation chemically
     reasonable at all, against an ETKDG (CSD-derived torsion library) ensemble?

Both are **intramolecular and protein-free**. FINDING 029 measured clash and strain
against the co-folded pocket at -0.0632, anti-selective; nothing here touches a protein
atom, which is the only reason arm B is worth running at all.

Donor harvest, leakage filters, MCS and the heme frame are **imported from
`fragment_transfer.py`**, not rewritten, so the two findings share one implementation.

    python scripts/structure/conformer_transfer.py torsions --workers 10
    python scripts/structure/conformer_transfer.py armb
    python scripts/structure/conformer_transfer.py answer     # <- the gate
    python scripts/structure/conformer_transfer.py evaluate   # only if answer passes
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))

import fragment_transfer as FT  # noqa: E402  - donors, filters, MCS, frame

DATA = REPO / "data" / "processed"
SCRATCH = Path(os.environ.get("CT_SCRATCH", str(FT.SCRATCH.parent / "conformer_transfer")))

# ---- pre-registered constants; none of these may move ---------------------
MIN_MCS_ATOMS = FT.MIN_MCS_ATOMS          # 6, primary arm, same as 030
TANIMOTO_CUT = FT.TANIMOTO_CUT            # 0.90
FE_CUT = FT.FE_CUT                        # 10 A
MIN_DONORS = FT.MIN_DONORS                # 3
CYP3A_UNIPROT = FT.CYP3A_UNIPROT
SEED = FT.SEED                            # 20260922
TIE_DRAWS = FT.TIE_DRAWS                  # 64
NULL_DRAWS = FT.NULL_DRAWS                # 2000
BOOT_DRAWS = FT.BOOT_DRAWS                # 10000
SCRAMBLE_REDRAWS = FT.SCRAMBLE_REDRAWS    # 8

ETKDG_CONFS = 300
ETKDG_PRUNE = 0.5
AR_PCT_BAR = 0.65                         # answer-recognition pass rule
AR_P_BAR = 0.05

ROT_SMARTS = "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"


# ==========================================================================
# torsion machinery - definition, symmetry period, circular distance
# ==========================================================================

def _ranks(mol):
    """Symmetry classes: equal rank exactly for graph-automorphic atoms.

    Pre-registered as `CanonicalRankAtoms(breakTies=False)`. It is NOT resonance
    aware - see control C-SYM, which reports where that bites.
    """
    from rdkit import Chem
    return list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=False))


def _sym_order(mol, ranks, a, b, c, d):
    """lcm of the two ends' local symmetry orders, for the torsion (a,b,c,d)."""
    nb = mol.GetAtomWithIdx(b)
    nc = mol.GetAtomWithIdx(c)
    sb = sum(1 for x in nb.GetNeighbors()
             if x.GetIdx() != c and ranks[x.GetIdx()] == ranks[a])
    sc = sum(1 for x in nc.GetNeighbors()
             if x.GetIdx() != b and ranks[x.GetIdx()] == ranks[d])
    return math.lcm(max(sb, 1), max(sc, 1))


def _dihedral(P):
    """(..., 4, 3) -> degrees in (-180, 180]. Vectorised, right-handed IUPAC sign."""
    b0 = P[..., 0, :] - P[..., 1, :]
    b1 = P[..., 2, :] - P[..., 1, :]
    b2 = P[..., 3, :] - P[..., 2, :]
    b1n = b1 / np.maximum(np.linalg.norm(b1, axis=-1, keepdims=True), 1e-12)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    x = (v * w).sum(-1)
    y = (np.cross(b1n, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


def _fold(delta_deg, period_deg):
    """Folded circular deviation in [0, P/2]. `delta_deg` may be an array."""
    d = np.mod(np.mod(delta_deg, period_deg) + period_deg, period_deg)
    return np.minimum(d, period_deg - d)


def _circ_mean(angles_deg, period_deg):
    """Circular mean of angles that live on a period-P circle, plus resultant R."""
    k = 360.0 / period_deg
    th = np.radians(np.asarray(angles_deg, float) * k)
    s, c = np.sin(th).mean(), np.cos(th).mean()
    R = float(np.hypot(s, c))
    mu = math.degrees(math.atan2(s, c)) / k
    return mu, R


def _rot_bonds(mol):
    from rdkit import Chem
    patt = Chem.MolFromSmarts(ROT_SMARTS)
    out = set()
    for i, j in mol.GetSubstructMatches(patt, uniquify=True):
        b = mol.GetBondBetweenAtoms(i, j)
        if b is None or b.IsInRing():
            continue
        out.add((min(i, j), max(i, j)))
    return out


def _pattern_torsions(patt, qmol, qmatch, dmol, dmatch):
    """Shared torsions, in PATTERN index space, canonical a and d.

    A quadruple qualifies only if the b-c bond is an acyclic rotatable bond in BOTH
    molecules under the given matches, and a, d are the in-pattern neighbours of b and
    c with the smallest PATTERN index - a numbering both molecules share, so the same
    physical torsion is measured on both sides.
    """
    qrot, drot = _rot_bonds(qmol), _rot_bonds(dmol)
    adj = {i: [] for i in range(patt.GetNumAtoms())}
    for bnd in patt.GetBonds():
        u, v = bnd.GetBeginAtomIdx(), bnd.GetEndAtomIdx()
        adj[u].append(v)
        adj[v].append(u)
    tors = []
    for bnd in patt.GetBonds():
        b, c = bnd.GetBeginAtomIdx(), bnd.GetEndAtomIdx()
        for bb, cc in ((b, c), (c, b)):
            qb, qc = qmatch[bb], qmatch[cc]
            db, dc = dmatch[bb], dmatch[cc]
            if (min(qb, qc), max(qb, qc)) not in qrot:
                continue
            if (min(db, dc), max(db, dc)) not in drot:
                continue
            na = sorted(x for x in adj[bb] if x != cc)
            nd = sorted(x for x in adj[cc] if x != bb)
            if not na or not nd:
                continue
            tors.append((na[0], bb, cc, nd[0]))
            break
    # de-duplicate on the (b,c) pattern bond
    seen, keep = set(), []
    for t in tors:
        k = (min(t[1], t[2]), max(t[1], t[2]))
        if k in seen:
            continue
        seen.add(k)
        keep.append(t)
    return keep


# ==========================================================================
# stage 1 - arm A: query pose torsions and donor torsions
# ==========================================================================

def _query_mol_and_poses(lig, smiles, want_xtal, pdb):
    """Pose mol (coordinate order) + raw Cartesian coords of 20 poses (+ crystal)."""
    f = FT.SCRATCH / "frames" / f"{lig}.npz"
    if not f.exists():
        return None, None, None, "no-frames"
    d = np.load(f, allow_pickle=True)
    raw = np.asarray(d["raw"], float)
    elems = [str(e) for e in d["elems"]]
    names = [str(x) for x in d["names"]]
    mol, why = FT._mol_from_3d(raw[0], elems, smiles)
    if mol is None:
        return None, None, None, f"query-mol-{why}"
    xs = "not-requested"
    if want_xtal:
        cx, xs = _crystal_xyz(lig, pdb, mol)
        if cx is not None:
            raw = np.concatenate([raw, cx[None, :, :]], axis=0)
    return mol, raw, (names, xs), "ok"


def _crystal_xyz(lig, pdb, query_mol):
    """The query's OWN crystal ligand, raw Cartesian, in POSE atom order.

    DIAGNOSTIC ONLY - it is the answer. Atom order is transferred by substructure
    match, never index-for-index (the trap that once halved every LDDT-PLI here).
    """
    import pandas as pd
    from rdkit import Chem

    donors = pd.read_parquet(FT.SCRATCH / "donors.parquet")
    sub = donors[(donors["pdb"] == pdb) & (donors["lig"] == lig)]
    if len(sub) == 0:
        return None, "no-crystal-row"
    r = sub.iloc[0]
    cm = Chem.MolFromMolBlock(r["molblock"], sanitize=True)
    if cm is None:
        return None, "crystal-mol-unreadable"
    if cm.GetNumAtoms() != query_mol.GetNumAtoms():
        return None, "crystal-atom-count-mismatch"
    xyz = cm.GetConformer().GetPositions()
    m = cm.GetSubstructMatch(query_mol)
    if m and len(m) == query_mol.GetNumAtoms():
        return xyz[list(m)], "ok"
    m = query_mol.GetSubstructMatch(cm)
    if m and len(m) == query_mol.GetNumAtoms():
        inv = np.empty(len(m), int)
        inv[list(m)] = np.arange(len(m))
        return xyz[inv], "ok"
    return None, "no-atom-order-map"


def _run_query_torsions(args):
    lig, smiles, pdb, want_xtal = args
    import pandas as pd
    from rdkit import Chem, DataStructs, RDLogger
    RDLogger.DisableLog("rdApp.*")

    qmol, raw, extra, status = _query_mol_and_poses(lig, smiles, want_xtal, pdb)
    if qmol is None:
        return {"ligand": lig, "status": status}
    names, xtal_status = extra
    qranks = _ranks(qmol)
    q_n_rot = len(_rot_bonds(qmol))

    donors = pd.read_parquet(FT.SCRATCH / "donors.parquet")
    qfp = FT._fp(qmol)

    # ---- leakage filters L0-L5, identical to FINDING 030, each counted -----
    n0 = len(donors)
    counts = {"start": n0}
    keep = np.ones(n0, bool)

    l1 = donors["pdb"].values == pdb
    counts["L1_same_entry"] = int(l1.sum())
    keep &= ~l1

    l2 = donors["lig"].values == lig
    counts["L2_same_ligand_code"] = int((l2 & keep).sum())
    keep &= ~l2

    sim = {}
    for code, sub in donors.groupby("lig"):
        m = Chem.MolFromMolBlock(sub["molblock"].iloc[0], sanitize=True)
        sim[code] = (0.0 if m is None
                     else DataStructs.TanimotoSimilarity(qfp, FT._fp(m)))
    l3 = np.array([sim.get(c, 0.0) >= TANIMOTO_CUT for c in donors["lig"].values])
    counts["L3_tanimoto_ge_0.90"] = int((l3 & keep).sum())
    keep &= ~l3

    l5 = ~(donors["closest_fe"].values <= FE_CUT)
    counts["L5_not_active_site"] = int((l5 & keep).sum())
    keep &= ~l5

    unclass = donors["uniprot"].isna().values
    counts["L0_unclassified_target"] = int((unclass & keep).sum())
    keep &= ~unclass

    l4 = np.array([u in CYP3A_UNIPROT for u in donors["uniprot"].values])
    counts["L4_cyp3a_subfamily"] = int((l4 & keep).sum())

    out = {"ligand": lig, "status": "ok", "filters": counts, "names": names,
           "n_poses": int(len(raw)), "q_n_rot_bonds": q_n_rot,
           "crystal_row": (len(raw) - 1) if xtal_status == "ok" else None,
           "crystal_status": xtal_status,
           "n_heavy": int(qmol.GetNumAtoms())}

    mcs_cache: dict = {}
    for arm, mask in (("no3a", keep & ~l4), ("with3a", keep)):
        recs, meta = [], []
        tally = {"legal": 0, "mcs_ge6": 0, "mcs_matched": 0, "has_shared_torsion": 0}
        for i in np.nonzero(mask)[0]:
            tally["legal"] += 1
            row = donors.iloc[int(i)]
            code = row["lig"]
            if code not in mcs_cache:
                m = Chem.MolFromMolBlock(row["molblock"], sanitize=True)
                mcs_cache[code] = (FT._mcs_smarts(qmol, m, False)
                                   if m is not None else (None, 0))
            smarts, natoms = mcs_cache[code]
            if smarts is None or natoms < MIN_MCS_ATOMS:
                continue
            tally["mcs_ge6"] += 1
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue
            dmol = Chem.MolFromMolBlock(row["molblock"], sanitize=True)
            if dmol is None:
                continue
            qm = qmol.GetSubstructMatch(patt, useChirality=False)
            dm = dmol.GetSubstructMatch(patt, useChirality=False)
            if not qm or not dm:
                continue
            tally["mcs_matched"] += 1
            tors = _pattern_torsions(patt, qmol, qm, dmol, dm)
            if not tors:
                continue
            tally["has_shared_torsion"] += 1
            dranks = _ranks(dmol)
            dxyz = dmol.GetConformer().GetPositions()

            keys, periods, dangs = [], [], []
            for (a, b, c, d) in tors:
                qa, qb, qc, qd = qm[a], qm[b], qm[c], qm[d]
                da, db, dc, dd = dm[a], dm[b], dm[c], dm[d]
                s = math.lcm(_sym_order(qmol, qranks, qa, qb, qc, qd),
                             _sym_order(dmol, dranks, da, db, dc, dd))
                keys.append((qa, qb, qc, qd))
                periods.append(360.0 / s)
                dangs.append(float(_dihedral(dxyz[[da, db, dc, dd]][None, :, :])[0]))
            qidx = np.array([k for k in keys], int)                 # (T, 4)
            qang = _dihedral(raw[:, qidx, :])                       # (poses, T)
            dev = _fold(qang - np.asarray(dangs)[None, :],
                        np.asarray(periods)[None, :])
            recs.append(dev.mean(axis=1).astype(np.float32))
            meta.append({"pdb": row["pdb"], "lig": code, "chain": row["chain"],
                         "seqid": int(row["seqid"]), "uniprot": row["uniprot"],
                         "target_key": row["target_key"], "k": int(natoms),
                         "n_tors": len(tors), "tanimoto": float(sim.get(code, 0.0)),
                         "keys": [[int(x) for x in k] for k in keys],
                         "periods": [float(p) for p in periods],
                         "donor_angles": dangs})
        out[arm] = {"n_donors": len(recs), "meta": meta, "tally": tally,
                    "d": (np.stack(recs, 1).tolist() if recs else [])}
    return out


def stage_torsions(workers: int, xtal: bool) -> None:
    import pandas as pd

    val = pd.read_csv(DATA / "validation_ligands.csv")
    out = SCRATCH / ("tors_xtal" if xtal else "tors")
    out.mkdir(parents=True, exist_ok=True)
    jobs = [(r.id, r.smiles, r.pdb, xtal) for r in val.itertuples()
            if not (out / f"{r.id}.json").exists()]
    print(f"torsions[{'xtal' if xtal else 'primary'}]: {len(jobs)} to do", flush=True)
    if not jobs:
        return
    if workers <= 1:
        for j in jobs:
            r = _run_query_torsions(j)
            (out / f"{j[0]}.json").write_text(json.dumps(r))
            print(f"  {j[0]}: {r.get('status')} "
                  f"no3a={r.get('no3a', {}).get('n_donors')}", flush=True)
        return
    import concurrent.futures as cf
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(_run_query_torsions, jobs):
            (out / f"{r['ligand']}.json").write_text(json.dumps(r))
            print(f"  {r['ligand']}: {r.get('status')} "
                  f"no3a={r.get('no3a', {}).get('n_donors')}", flush=True)


# ==========================================================================
# stage 2 - arm B: prior-free internal plausibility. NO P450 data, NO protein.
# ==========================================================================

def _armb_one(args):
    lig, smiles, pdb = args
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, TorsionFingerprints
    from rdkit.Geometry import Point3D
    RDLogger.DisableLog("rdApp.*")

    qmol, raw, extra, status = _query_mol_and_poses(lig, smiles, True, pdb)
    if qmol is None:
        return {"ligand": lig, "status": status}
    names, xtal_status = extra

    ens = Chem.MolFromSmiles(smiles)
    if ens is None:
        return {"ligand": lig, "status": "ens-smiles-bad"}
    ens = Chem.RemoveHs(ens)
    # atom order by SUBSTRUCTURE MATCH, never index-for-index. The pose mol carries
    # explicit zero-H atoms (it was built from 3D), so it is the poorer query of the
    # two; match the ensemble INTO it and invert.
    match = ens.GetSubstructMatch(qmol, useChirality=False)
    if match and len(match) == qmol.GetNumAtoms():
        mapping = np.asarray(match, int)      # pose atom k -> ens atom mapping[k]
    else:
        m2 = qmol.GetSubstructMatch(ens, useChirality=False)
        if not m2 or len(m2) != ens.GetNumAtoms():
            return {"ligand": lig, "status": "ens-atom-map-failed"}
        mapping = np.empty(len(m2), int)
        mapping[list(m2)] = np.arange(len(m2))

    ensH = Chem.AddHs(ens)
    p = AllChem.ETKDGv3()
    p.randomSeed = SEED
    p.pruneRmsThresh = ETKDG_PRUNE
    p.useSmallRingTorsions = True
    cids = list(AllChem.EmbedMultipleConfs(ensH, ETKDG_CONFS, p))
    if not cids:
        return {"ligand": lig, "status": "embed-failed"}
    try:
        AllChem.MMFFOptimizeMoleculeConfs(ensH, mmffVariant="MMFF94s", maxIters=500)
    except Exception:
        pass
    ensE = Chem.RemoveHs(ensH)

    tl, tlr = TorsionFingerprints.CalculateTorsionLists(ensE)
    n_tors = len(tl) + len(tlr)
    if n_tors == 0:
        return {"ligand": lig, "status": "no-torsions", "names": names,
                "n_poses": int(len(raw)),
                "crystal_row": (len(raw) - 1) if xtal_status == "ok" else None}

    ens_ang = [TorsionFingerprints.CalculateTorsionAngles(ensE, tl, tlr, confId=c)
               for c in cids]

    # append each pose as a conformer of the ENSEMBLE molecule, atom order mapped
    pose_ids = []
    for pi in range(len(raw)):
        conf = Chem.Conformer(ensE.GetNumAtoms())
        for k in range(len(mapping)):
            x, y, z = raw[pi, k]
            conf.SetAtomPosition(int(mapping[k]), Point3D(float(x), float(y), float(z)))
        pose_ids.append(ensE.AddConformer(conf, assignId=True))

    tfd_min, tfd_mean, well = [], [], []
    for cid in pose_ids:
        pa = TorsionFingerprints.CalculateTorsionAngles(ensE, tl, tlr, confId=cid)
        vals = np.array([TorsionFingerprints.CalculateTFD(pa, ea) for ea in ens_ang])
        tfd_min.append(float(vals.min()))
        tfd_mean.append(float(vals.mean()))
        # etkdg_well: per acyclic torsion, folded distance to the NEAREST value that
        # torsion takes anywhere in the ensemble
        ds = []
        for ti, (_atoms, maxdev) in enumerate(tl):
            period = 360.0 if maxdev >= 179.0 else 180.0
            pv = pa[ti][0][0]
            ev = np.array([ea[ti][0][0] for ea in ens_ang])
            ds.append(float(_fold(pv - ev, period).min()))
        well.append(float(np.mean(ds)) if ds else float("nan"))

    strain = _mmff_strain(ens, mapping, raw)
    return {"ligand": lig, "status": "ok", "names": names,
            "n_poses": int(len(raw)), "n_confs": len(cids),
            "n_tfd_torsions": n_tors, "n_acyclic_torsions": len(tl),
            "crystal_row": (len(raw) - 1) if xtal_status == "ok" else None,
            "crystal_status": xtal_status,
            "etkdg_tfd_min": tfd_min, "etkdg_tfd_mean": tfd_mean,
            "etkdg_well": well, "mmff_strain": strain}


def _mmff_strain(ens, mapping, raw):
    """MMFF94s local strain: pose energy (H relaxed, heavy fixed) - free minimum.

    Declared in the pre-registration as a REPLICATION of FINDING 025's `strain` term
    (+0.0105 against a +0.0138 floor, negative within-ligand rho), not a candidate.

    Built on the SMILES-derived molecule, not on the 3D-perceived one, because the
    latter carries explicit zero-H atoms and MMFF would type it as a radical soup.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    out = []
    for pi in range(len(raw)):
        try:
            m = Chem.Mol(ens)
            m.RemoveAllConformers()
            conf = Chem.Conformer(m.GetNumAtoms())
            for k in range(len(mapping)):
                x, y, z = raw[pi, k]
                conf.SetAtomPosition(int(mapping[k]),
                                     Point3D(float(x), float(y), float(z)))
            m.AddConformer(conf, assignId=True)
            mh = Chem.AddHs(m, addCoords=True)
            props = AllChem.MMFFGetMoleculeProperties(mh, mmffVariant="MMFF94s")
            if props is None:
                out.append(float("nan"))
                continue
            ff = AllChem.MMFFGetMoleculeForceField(mh, props)
            for a in mh.GetAtoms():
                if a.GetAtomicNum() > 1:
                    ff.AddFixedPoint(a.GetIdx())
            ff.Minimize(maxIts=400)
            e_pose = ff.CalcEnergy()
            ff2 = AllChem.MMFFGetMoleculeForceField(
                mh, AllChem.MMFFGetMoleculeProperties(mh, mmffVariant="MMFF94s"))
            ff2.Minimize(maxIts=2000)
            out.append(float(e_pose - ff2.CalcEnergy()))
        except Exception:
            out.append(float("nan"))
    return out


def stage_armb(workers: int) -> None:
    import pandas as pd

    val = pd.read_csv(DATA / "validation_ligands.csv")
    out = SCRATCH / "armb"
    out.mkdir(parents=True, exist_ok=True)
    jobs = [(r.id, r.smiles, r.pdb) for r in val.itertuples()
            if not (out / f"{r.id}.json").exists()]
    print(f"armb: {len(jobs)} to do", flush=True)
    if not jobs:
        return
    if workers <= 1:
        for j in jobs:
            r = _armb_one(j)
            (out / f"{j[0]}.json").write_text(json.dumps(r))
            print(f"  {j[0]}: {r.get('status')}", flush=True)
        return
    # as_completed, NOT map: one macrocyclic ligand (erythromycin, 300 ETKDG
    # embeddings of a 14-membered ring) blocks an ordered map for half an hour and
    # every finished result behind it sits unwritten. Write on completion.
    import concurrent.futures as cf
    with cf.ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_armb_one, j): j[0] for j in jobs}
        for fut in cf.as_completed(futs):
            r = fut.result()
            (out / f"{r['ligand']}.json").write_text(json.dumps(r))
            print(f"  {r['ligand']}: {r.get('status')}", flush=True)


# ==========================================================================
# feature assembly
# ==========================================================================

def _arma_features(rec, arm, n_poses):
    """tor_cons / tor_best / tor_med / tor_cmean / tor_wgt for one query."""
    d = np.asarray(rec[arm]["d"], float)          # (poses, donors), degrees
    if d.size == 0:
        return None
    d = d[:n_poses]
    feats = {"tor_cons": d.mean(1), "tor_best": d.min(1),
             "tor_med": np.median(d, 1)}

    # per-query-torsion donor distributions, canonical match on both sides
    bykey: dict[tuple, list] = {}
    for m in rec[arm]["meta"]:
        for k, per, ang in zip(m["keys"], m["periods"], m["donor_angles"]):
            bykey.setdefault(tuple(k), []).append((per, ang))
    if not bykey:
        return None
    return feats, bykey


def _arma_dist_features(rec, arm, raw_angles):
    """tor_cmean and tor_wgt, which need the pose's own angle per torsion."""
    bykey: dict[tuple, list] = {}
    for m in rec[arm]["meta"]:
        for k, per, ang in zip(m["keys"], m["periods"], m["donor_angles"]):
            bykey.setdefault(tuple(k), []).append((per, ang))
    if not bykey:
        return None, None
    cm, wg, W = [], [], []
    for key, vals in bykey.items():
        period = min(v[0] for v in vals)          # coarsest period seen for this bond
        angs = [v[1] for v in vals]
        mu, R = _circ_mean(angs, period)
        pose_ang = raw_angles[key]                # (poses,)
        dev = _fold(pose_ang - mu, period)
        cm.append(dev)
        wg.append(dev * R)
        W.append(R)
    cm = np.stack(cm, 1)
    wg = np.stack(wg, 1)
    Wt = float(np.sum(W))
    return cm.mean(1), (wg.sum(1) / Wt if Wt > 1e-9 else cm.mean(1))


_RAW_CACHE: dict = {}


def _pose_angles(lig, keys, xtal):
    """Dihedrals of the requested atom quadruples for every pose of `lig`."""
    import pandas as pd
    ck = (lig, bool(xtal))
    if ck not in _RAW_CACHE:
        val = pd.read_csv(DATA / "validation_ligands.csv")
        row = val[val["id"] == lig].iloc[0]
        qmol, raw, _extra, _st = _query_mol_and_poses(lig, row["smiles"], xtal,
                                                      row["pdb"])
        _RAW_CACHE[ck] = None if qmol is None else raw
    raw = _RAW_CACHE[ck]
    if raw is None:
        return None
    return {k: _dihedral(raw[:, list(k), :]) for k in keys}


# ==========================================================================
# stage 3 - THE ANSWER-RECOGNITION TEST. It can stop the experiment.
# ==========================================================================

def _percentiles(vals, k):
    """Fraction of the 20 predicted poses whose feature value is WORSE (higher)."""
    c, p = vals[k], np.delete(vals, k)
    return float((p > c).mean()), bool(c < np.median(p)), float(c - p.mean())


def stage_answer() -> None:
    from scipy import stats

    out: dict = {"note": "answer-recognition: the crystal scored as a 21st pose. "
                         "DIAGNOSTIC ONLY - never enters a donor set or a selector.",
                 "pass_rule": {"mean_percentile_ge": AR_PCT_BAR,
                               "binomial_p_lt": AR_P_BAR}}

    # ---------- arm A ----------
    prim = {f.stem: json.loads(f.read_text())
            for f in (SCRATCH / "tors").glob("*.json")}
    rows = []
    for f in sorted((SCRATCH / "tors_xtal").glob("*.json")):
        r = json.loads(f.read_text())
        if r.get("status") != "ok" or r.get("crystal_row") is None:
            continue
        d = np.asarray(r["no3a"]["d"], float)
        if d.size == 0 or d.shape[1] < MIN_DONORS:
            continue
        rows.append((r, d))

    worst = 0.0
    for r, d in rows:
        p = prim.get(r["ligand"])
        if p is None or not p.get("no3a", {}).get("d"):
            continue
        dp = np.asarray(p["no3a"]["d"], float)
        if dp.shape[1] == d.shape[1] and dp.shape[0] <= d.shape[0]:
            worst = max(worst, float(np.abs(dp - d[:dp.shape[0]]).max()))
    out["C5_xtal_rerun_reproduces_primary"] = {
        "max_abs_delta_deg": worst, "pass": bool(worst < 1e-5), "ligands": len(rows)}

    a = {}
    for name in ("tor_cons", "tor_best", "tor_med", "tor_cmean", "tor_wgt"):
        pct, btm, dl = [], 0, []
        for r, d in rows:
            k = r["crystal_row"]
            if name in ("tor_cons", "tor_best", "tor_med"):
                v = {"tor_cons": d.mean(1), "tor_best": d.min(1),
                     "tor_med": np.median(d, 1)}[name]
            else:
                keys = sorted({tuple(x) for m in r["no3a"]["meta"] for x in m["keys"]})
                ang = _pose_angles(r["ligand"], keys, True)
                if ang is None:
                    continue
                cm, wg = _arma_dist_features(r, "no3a", ang)
                v = cm if name == "tor_cmean" else wg
            p, b, dd = _percentiles(v, k)
            pct.append(p)
            btm += int(b)
            dl.append(dd)
        a[name] = {"ligands": len(pct),
                   "mean_percentile_of_crystal": float(np.mean(pct)),
                   "median_percentile_of_crystal": float(np.median(pct)),
                   "crystal_better_than_median_pose": btm,
                   "mean_delta_deg": float(np.mean(dl)),
                   "sign_binom_p": float(stats.binomtest(btm, len(pct), 0.5).pvalue),
                   "PASS": bool(np.mean(pct) >= AR_PCT_BAR
                                and btm > len(pct) / 2
                                and stats.binomtest(btm, len(pct),
                                                    0.5).pvalue < AR_P_BAR)}
    out["armA"] = a

    # ---------- arm B ----------
    b = {}
    recs = [json.loads(f.read_text()) for f in sorted((SCRATCH / "armb").glob("*.json"))]
    for name in ("etkdg_tfd_min", "etkdg_tfd_mean", "etkdg_well", "mmff_strain"):
        pct, btm, dl, skipped = [], 0, [], 0
        for r in recs:
            if r.get("status") != "ok" or r.get("crystal_row") is None:
                skipped += 1
                continue
            v = np.asarray(r[name], float)
            if not np.isfinite(v).all() or v.std() < 1e-12:
                skipped += 1
                continue
            p, bb, dd = _percentiles(v, r["crystal_row"])
            pct.append(p)
            btm += int(bb)
            dl.append(dd)
        b[name] = {"ligands": len(pct), "skipped": skipped,
                   "mean_percentile_of_crystal": float(np.mean(pct)),
                   "median_percentile_of_crystal": float(np.median(pct)),
                   "crystal_better_than_median_pose": btm,
                   "mean_delta": float(np.mean(dl)),
                   "sign_binom_p": float(stats.binomtest(btm, len(pct), 0.5).pvalue),
                   "PASS": bool(np.mean(pct) >= AR_PCT_BAR
                                and btm > len(pct) / 2
                                and stats.binomtest(btm, len(pct),
                                                    0.5).pvalue < AR_P_BAR)}
    out["armB"] = b

    out["GATE"] = {
        "armA_primary_pass": a["tor_cons"]["PASS"],
        "armB_primary_pass": b["etkdg_tfd_min"]["PASS"],
        "any_pass": any(v["PASS"] for v in list(a.values()) + list(b.values())),
    }
    p = DATA / "conformer_answer_recognition.json"
    p.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print(f"\nwrote {p}", flush=True)


# ==========================================================================
# stage 4 - coverage, controls, and (only if the gate passes) the bars
# ==========================================================================

def stage_evaluate(force: bool = False) -> None:
    import pandas as pd
    from scipy import stats

    ar = json.loads((DATA / "conformer_answer_recognition.json").read_text())
    gate = ar["GATE"]["any_pass"]

    scored = pd.read_csv(DATA / "poses_scored_val87b.csv")
    scored = scored[scored["arm"] == "unsteered"]
    truth = {(r.ligand, r.sample): float(r.lddt_pli) for r in scored.itertuples()}
    xeng = pd.read_csv(DATA / "xeng_val87b.csv")
    xe = {(r.ligand, r.sample): float(r.xeng) for r in xeng.itertuples()}
    donors_all = pd.read_parquet(FT.SCRATCH / "donors.parquet")

    per, filt = {}, {}
    for f in sorted((SCRATCH / "tors").glob("*.json")):
        r = json.loads(f.read_text())
        per[r["ligand"]] = r
        for k, v in r.get("filters", {}).items():
            filt[k] = filt.get(k, 0) + int(v)
    armb = {}
    for f in sorted((SCRATCH / "armb").glob("*.json")):
        r = json.loads(f.read_text())
        armb[r["ligand"]] = r

    ligs = sorted(per)
    cov = {arm: np.array([per[l].get(arm, {}).get("n_donors", 0) for l in ligs])
           for arm in ("no3a", "with3a")}
    nrot = np.array([per[l].get("q_n_rot_bonds", 0) for l in ligs])

    coverage = {
        "queries": len(ligs),
        "query_mol_failures": sorted(l for l in ligs if per[l].get("status") != "ok"),
        "queries_with_zero_rotatable_bonds": int((nrot == 0).sum()),
        "rigid_query_ids": sorted(l for l, n in zip(ligs, nrot) if n == 0),
        "median_rotatable_bonds": float(np.median(nrot)),
        "filters_summed_over_queries": filt,
        "donor_pool": {
            "observations": int(len(donors_all)),
            "codes": int(donors_all["lig"].nunique()),
            "entries": int(donors_all["pdb"].nunique()),
            "targets": int(donors_all["target_key"].nunique())},
    }
    for arm in ("no3a", "with3a"):
        c = cov[arm]
        coverage[arm] = {
            "with_ge1_donor": int((c >= 1).sum()),
            "with_ge3_donors": int((c >= MIN_DONORS).sum()),
            "median_donors": float(np.median(c)), "mean_donors": float(c.mean()),
            "p90_donors": float(np.percentile(c, 90)), "max_donors": int(c.max()),
            "zero_donor_queries": int((c == 0).sum())}
    ntors = [m["n_tors"] for l in ligs for m in per[l].get("no3a", {}).get("meta", [])]
    coverage["shared_torsions_per_donor"] = {
        "mean": float(np.mean(ntors)) if ntors else 0.0,
        "median": float(np.median(ntors)) if ntors else 0.0,
        "pairs": len(ntors)}
    coverage["armB"] = {
        "ok": sum(1 for v in armb.values() if v.get("status") == "ok"),
        "no_torsions": sum(1 for v in armb.values()
                           if v.get("status") == "no-torsions"),
        "other_failures": sorted(v["ligand"] for v in armb.values()
                                 if v.get("status") not in ("ok", "no-torsions")),
        "mean_conformers": float(np.mean([v["n_confs"] for v in armb.values()
                                          if v.get("status") == "ok"]))}

    print("=" * 72)
    print("COVERAGE (reported before any selector number)")
    print(json.dumps(coverage, indent=2))
    print("=" * 72, flush=True)

    # C3 - the feature must not be constant within a ligand. Computed before the
    # gate check, because a control is owed whether or not a selector is built.
    c3 = {}
    for name in ("tor_cons", "etkdg_tfd_min", "etkdg_well", "mmff_strain"):
        sds = []
        for l in ligs:
            if name == "tor_cons":
                r = per[l]
                if r.get("no3a", {}).get("n_donors", 0) < MIN_DONORS:
                    continue
                v = np.asarray(r["no3a"]["d"], float)[:len(r["names"])].mean(1)
            else:
                r = armb.get(l)
                if r is None or r.get("status") != "ok":
                    continue
                v = np.asarray(r[name], float)[:len(per[l]["names"])]
                if not np.isfinite(v).all():
                    continue
            sds.append(float(np.std(v)))
        if sds:
            c3[name] = {"queries": len(sds),
                        "nonconstant_frac": float((np.array(sds) > 1e-9).mean()),
                        "min_within_ligand_sd": float(min(sds))}

    results = {"coverage": coverage,
               "answer_recognition_gate": ar["GATE"],
               "C3_nonconstant": c3,
               "C4_frame_identity": FT.check_frames(),
               "C_NUM": _numbering_control(scored),
               "C_SYM": _sym_control()}
    print("C3:", json.dumps(c3), flush=True)
    print("C4:", results["C4_frame_identity"])
    print("C_NUM:", results["C_NUM"])
    print("C_SYM:", json.dumps(results["C_SYM"]), flush=True)

    if not gate and not force:
        results["verdict"] = ("STOPPED at the answer-recognition gate, as "
                              "pre-registered. No selector was built.")
        p = DATA / "conformer_transfer_primary.json"
        p.write_text(json.dumps(results, indent=2))
        print("\nGATE FAILED - no selector built, as pre-registered.")
        print(f"wrote {p}", flush=True)
        return

    # ---- the bars (reached only if the gate passed, or --force for the record) ----
    results["bars"] = {}
    for arm in ("no3a", "with3a"):
        subset = [l for l, c in zip(ligs, cov[arm]) if c >= MIN_DONORS]
        if not subset:
            continue
        T, XE = [], []
        for l in subset:
            names = per[l]["names"]
            T.append(np.array([truth[(l, n)] for n in names]))
            XE.append(np.array([xe[(l, n)] for n in names]))
        featnames = ["tor_cons", "tor_best", "tor_med", "tor_cmean", "tor_wgt"]
        F = {k: [] for k in featnames}
        for l in subset:
            r = per[l]
            n = len(r["names"])
            d = np.asarray(r[arm]["d"], float)[:n]
            F["tor_cons"].append(d.mean(1))
            F["tor_best"].append(d.min(1))
            F["tor_med"].append(np.median(d, 1))
            keys = sorted({tuple(x) for m in r[arm]["meta"] for x in m["keys"]})
            ang = _pose_angles(l, keys, False)
            cm, wg = _arma_dist_features(r, arm, ang)
            F["tor_cmean"].append(cm[:n])
            F["tor_wgt"].append(wg[:n])
        if arm == "no3a":
            for k in ("etkdg_tfd_min", "etkdg_tfd_mean", "etkdg_well", "mmff_strain"):
                vals = []
                ok = True
                for l in subset:
                    r = armb.get(l)
                    if r is None or r.get("status") != "ok":
                        ok = False
                        break
                    vals.append(np.asarray(r[k], float)[:len(per[l]["names"])])
                if ok:
                    F[k] = vals
        results["bars"][arm] = _bars(subset, T, XE, F, arm, per, donors_all)

    p = DATA / "conformer_transfer_primary.json"
    p.write_text(json.dumps(results, indent=2, default=float))
    print(json.dumps(results.get("bars", {}), indent=2, default=float)[:14000])
    print(f"\nwrote {p}", flush=True)


def _bars(subset, T, XE, F, arm, per, donors_all):
    from scipy import stats

    n = len(subset)
    oracle = float(np.mean([t.max() for t in T]))
    baseline = float(np.mean([t.mean() for t in T]))
    block = {"subset_ligands": n, "oracle": oracle, "random_baseline": baseline}

    sel = {}
    for name, feats in F.items():
        if len(feats) != n:
            continue
        s = [-np.nan_to_num(f, nan=np.nanmax(f) if np.isfinite(f).any() else 0.0)
             for f in feats]
        sel[name] = FT._select(s, T, np.random.default_rng(SEED))
    inc_v, inc_p = FT._select([-x for x in XE], T, np.random.default_rng(SEED))
    sel["incumbent"] = (inc_v, inc_p)

    block["C3_nonconstant"] = {
        k: float((np.array([f.std() for f in F[k]]) > 1e-9).mean())
        for k in F if len(F[k]) == n}
    block["C3_min_within_sd"] = {
        k: float(min(f.std() for f in F[k])) for k in F if len(F[k]) == n}
    block["selection"] = {k: float(v.mean()) for k, (v, _) in sel.items()}
    block["gain_vs_random"] = {k: float(v.mean() - baseline) for k, (v, _) in sel.items()}

    rg = np.random.default_rng(SEED + 1)
    null = np.empty(NULL_DRAWS)
    for i in range(NULL_DRAWS):
        tot = 0.0
        for t in T:
            tot += t[int(np.argmax(rg.standard_normal(len(t))))]
        null[i] = tot / n - baseline
    block["null_random_feature"] = {
        "draws": NULL_DRAWS, "p95": float(np.percentile(null, 95)),
        "p99": float(np.percentile(null, 99)), "max": float(null.max())}
    block["bar1_pass"] = {k: bool(g > block["null_random_feature"]["p99"])
                          for k, g in block["gain_vs_random"].items()}
    block["null_exceedance"] = {k: float((null >= g).mean())
                                for k, g in block["gain_vs_random"].items()}

    wi = {}
    for k in F:
        if len(F[k]) != n:
            continue
        rhos = [stats.spearmanr(f, t).statistic for f, t in zip(F[k], T)
                if np.nanstd(f) > 1e-12]
        rhos = np.array([r for r in rhos if np.isfinite(r)])
        wi[k] = {"n": int(len(rhos)), "mean_rho": float(rhos.mean()),
                 "median_rho": float(np.median(rhos)),
                 "correct_sign_frac": float((rhos < 0).mean()),
                 "sign_binom_p": float(stats.binomtest(int((rhos < 0).sum()),
                                                       len(rhos), 0.5).pvalue)}
    block["within_ligand"] = wi

    paired = {}
    for name, (v, picks) in sel.items():
        if name == "incumbent":
            continue
        delta = v - inc_v
        tied = int(sum(int((picks[:, i] == inc_p[:, i]).all()) for i in range(n)))
        nz = delta[np.abs(delta) > 1e-12]
        bg = np.random.default_rng(SEED + 2)
        boot = np.array([delta[bg.integers(0, n, n)].mean() for _ in range(BOOT_DRAWS)])
        paired[name] = {
            "mean_delta": float(delta.mean()),
            "boot_ci95": [float(np.percentile(boot, 2.5)),
                          float(np.percentile(boot, 97.5))],
            "wilcoxon_p": float(stats.wilcoxon(nz).pvalue) if len(nz) else 1.0,
            "n_informative": int(len(nz)), "ligands_tied_same_pose": tied,
            "better": int((delta > 1e-12).sum()), "worse": int((delta < -1e-12).sum())}
    block["paired_vs_incumbent"] = paired

    head = np.array([t.max() for t in T]) - inc_v
    block["complementarity"] = {}
    for name in ("tor_cons", "etkdg_tfd_min"):
        if name not in sel:
            continue
        r = stats.pearsonr(sel[name][0] - inc_v, head)
        block["complementarity"][name] = {"pearson_r": float(r.statistic),
                                          "p": float(r.pvalue)}
    return block


# ==========================================================================
# controls
# ==========================================================================

def _numbering_control(scored):
    """C-NUM, FINDING 021: rescore 15 random poses from scratch and match the CSV."""
    import pandas as pd

    from cypstruct import pose as P

    rng = np.random.default_rng(SEED)
    exact_zero = int((scored["lddt_pli"] == 0.0).sum())
    val = pd.read_csv(DATA / "validation_ligands.csv")
    meta = {r.id: r for r in val.itertuples()}
    sub = scored.sample(15, random_state=SEED)
    worst, checked, fails = 0.0, 0, []
    zeros = scored[scored["lddt_pli"] == 0.0]
    ref_cache: dict = {}
    for r in sub.itertuples():
        m = meta.get(r.ligand)
        if m is None:
            continue
        try:
            if m.pdb not in ref_cache:
                from cypstruct.targets import fetch_cif
                ref_cache[m.pdb] = P.load_structure(fetch_cif(m.pdb), r.ligand)
            ref = ref_cache[m.pdb]
            cif = FT.POOL / f"{r.ligand}__unsteered__s1" / f"{r.sample}.cif"
            model = P.load_structure(cif)
            perm = P.best_ligand_mapping(m.smiles, model, ref)
            got = P.lddt_pli(model, ref, lig_perm=perm)
            worst = max(worst, abs(got - float(r.lddt_pli)))
            checked += 1
        except Exception as e:                      # noqa: BLE001
            fails.append(f"{r.ligand}:{type(e).__name__}")
    return {"poses_rescored": checked, "max_abs_delta": worst,
            "pass": bool(checked >= 10 and worst < 1e-6),
            "unsteered_poses_scoring_exactly_zero": exact_zero,
            "exact_zero_rows": [{"ligand": r.ligand, "sample": r.sample,
                                 "bisy_rmsd": float(r.bisy_rmsd)}
                                for r in zeros.itertuples()],
            "errors": fails[:5]}


def _sym_control():
    """C-SYM: does the pre-registered symmetry rule return the right periods?"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    cases = [("phenyl (CCc1ccccc1)", "CCc1ccccc1", 180.0),
             ("t-butyl (CCC(C)(C)C)", "CCC(C)(C)C", 120.0),
             ("carboxylate (CCC(=O)[O-])", "CCC(=O)[O-]", 180.0),
             ("carboxylic acid (CCC(=O)O)", "CCC(=O)O", 180.0),
             ("nitro (CC[N+](=O)[O-])", "CC[N+](=O)[O-]", 180.0),
             ("plain sp3-sp3 (CCCO)", "CCCO", 360.0)]
    out = {}
    for label, smi, want in cases:
        m = Chem.MolFromSmiles(smi)
        m = Chem.RemoveHs(m)
        rk = _ranks(m)
        got = None
        for (b, c) in sorted(_rot_bonds(m)):
            na = sorted(x.GetIdx() for x in m.GetAtomWithIdx(b).GetNeighbors()
                        if x.GetIdx() != c)
            nd = sorted(x.GetIdx() for x in m.GetAtomWithIdx(c).GetNeighbors()
                        if x.GetIdx() != b)
            if not na or not nd:
                continue
            got = 360.0 / _sym_order(m, rk, na[0], b, c, nd[0])
        out[label] = {"expected": want, "got": got, "pass": bool(got == want)}

    # relabelling control: a symmetry-equivalent substructure match must not move delta
    m = Chem.AddHs(Chem.MolFromSmiles("CCc1ccccc1"))
    AllChem.EmbedMolecule(m, randomSeed=SEED)
    m = Chem.RemoveHs(m)
    xyz = m.GetConformer().GetPositions()
    rk = _ranks(m)
    b, c = 1, 2
    orthos = [x.GetIdx() for x in m.GetAtomWithIdx(c).GetNeighbors()
              if x.GetIdx() != b]
    per = 360.0 / _sym_order(m, rk, 0, b, c, orthos[0])
    a0 = _dihedral(xyz[[0, b, c, orthos[0]]][None])[0]
    a1 = _dihedral(xyz[[0, b, c, orthos[1]]][None])[0]
    out["relabel_invariance"] = {
        "period": per, "angle_choice_1": float(a0), "angle_choice_2": float(a1),
        "folded_difference_deg": float(_fold(np.array([a0 - a1]), per)[0]),
        "pass": bool(_fold(np.array([a0 - a1]), per)[0] < 1e-6)}
    return out


def stage_scramble() -> None:
    """C1, applied to the ANSWER-RECOGNITION test rather than to a selector.

    The pre-registration attaches C1 to a selector. No selector is built, because
    the gate failed, so C1 is applied where the measurement actually lives: does the
    crystal's percentile move when the donors' torsion values are replaced by random
    values drawn from the same global pool? If scrambled reproduces matched, the
    result is a property of the poses and not of the crystallographic record.
    """
    from scipy import stats

    rng = np.random.default_rng(SEED + 7)
    per = {}
    for f in sorted((SCRATCH / "tors_xtal").glob("*.json")):
        r = json.loads(f.read_text())
        if r.get("status") != "ok" or r.get("crystal_row") is None:
            continue
        if np.asarray(r["no3a"]["d"], float).shape[-1] < MIN_DONORS:
            continue
        per[r["ligand"]] = r
    subset = sorted(per)

    pool = np.asarray([a for l in subset for m in per[l]["no3a"]["meta"]
                       for a in m["donor_angles"]], float)

    matched_pct, matched_btm = [], 0
    for l in subset:
        d = np.asarray(per[l]["no3a"]["d"], float)
        p, b, _ = _percentiles(d.mean(1), per[l]["crystal_row"])
        matched_pct.append(p)
        matched_btm += int(b)

    reps = []
    for rep in range(SCRAMBLE_REDRAWS):
        pct, btm = [], 0
        for l in subset:
            r = per[l]
            keys = sorted({tuple(x) for m in r["no3a"]["meta"] for x in m["keys"]})
            ang = _pose_angles(l, keys, True)
            acc = []
            for m in r["no3a"]["meta"]:
                dev = [_fold(ang[tuple(k)] - pool[rng.integers(0, len(pool))], p_)
                       for k, p_ in zip(m["keys"], m["periods"])]
                acc.append(np.stack(dev, 1).mean(1))
            v = np.stack(acc, 1).mean(1)
            p, b, _ = _percentiles(v, r["crystal_row"])
            pct.append(p)
            btm += int(b)
        reps.append((float(np.mean(pct)), btm))

    out = {"ligands": len(subset),
           "matched_mean_percentile": float(np.mean(matched_pct)),
           "matched_crystal_better_than_median": matched_btm,
           "matched_binom_p": float(stats.binomtest(matched_btm, len(subset),
                                                    0.5).pvalue),
           "scrambled_mean_percentile": float(np.mean([r[0] for r in reps])),
           "scrambled_crystal_better_than_median_mean": float(
               np.mean([r[1] for r in reps])),
           "sd_over_redraws": float(np.std([r[0] for r in reps])),
           "per_redraw": reps,
           "reading": ("if scrambled sits at 0.50 and matched does not, the "
                       "crystallographic record is doing the work - in whichever "
                       "direction it is doing it")}
    p = DATA / "conformer_scrambled_control.json"
    p.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print(f"\nwrote {p}", flush=True)


# ==========================================================================
# stage 5 - diagnostics. POST-HOC, labelled as such, never features.
# ==========================================================================

def _armb_pct(armb, lig):
    r = armb.get(lig)
    if r is None or r.get("status") != "ok" or r.get("crystal_row") is None:
        return float("nan")
    v = np.asarray(r["etkdg_tfd_min"], float)
    if not np.isfinite(v).all() or v.std() < 1e-12:
        return float("nan")
    return _percentiles(v, r["crystal_row"])[0]


def stage_diagnose() -> None:
    """Why an internal-conformer prior cannot rank these poses.

    | D1 | is there a torsional consensus at all? donor spread vs the error to fix |
    | D2 | the TORSION ORACLE - score by deviation from the CRYSTAL's own torsions |
    | D3 | what arm A actually ranks by - its agreement with the prior-free terms |
    | D4 | how much of the shared-torsion set the C-SYM resonance gap touches |

    None of this is a candidate feature; promoting any of it needs a fresh
    pre-registration. D2 reads the answer and is a diagnostic OF THE TERM, not a
    predictor.
    """
    import pandas as pd
    from rdkit import Chem
    from scipy import stats

    scored = pd.read_csv(DATA / "poses_scored_val87b.csv")
    scored = scored[scored["arm"] == "unsteered"]
    truth = {(r.ligand, r.sample): float(r.lddt_pli) for r in scored.itertuples()}
    xeng = pd.read_csv(DATA / "xeng_val87b.csv")
    xe = {(r.ligand, r.sample): float(r.xeng) for r in xeng.itertuples()}

    rows = []
    for f in sorted((SCRATCH / "tors_xtal").glob("*.json")):
        r = json.loads(f.read_text())
        if r.get("status") != "ok" or r.get("crystal_row") is None:
            continue
        d = np.asarray(r["no3a"]["d"], float)
        if d.size == 0 or d.shape[1] < MIN_DONORS:
            continue
        rows.append((r, d))

    armb = {}
    for f in sorted((SCRATCH / "armb").glob("*.json")):
        v = json.loads(f.read_text())
        armb[v["ligand"]] = v

    out: dict = {"note": "POST-HOC diagnostics; not pre-registered; not features"}

    # ---- D1 / D2 ---------------------------------------------------------
    donor_sd, pose_sd, xtal_err, per = [], [], [], []
    T, ORA, INC, FEA = [], [], [], []
    for r, d in rows:
        lig = r["ligand"]
        names = r["names"]
        n = len(names)
        keys = sorted({tuple(x) for m in r["no3a"]["meta"] for x in m["keys"]})
        ang = _pose_angles(lig, keys, True)
        if ang is None:
            continue
        bykey: dict[tuple, list] = {}
        for m in r["no3a"]["meta"]:
            for k, p_, a_ in zip(m["keys"], m["periods"], m["donor_angles"]):
                bykey.setdefault(tuple(k), []).append((p_, a_))
        dsd, psd, xer, ora = [], [], [], []
        for key, vals in bykey.items():
            period = min(v[0] for v in vals)
            angs = np.array([v[1] for v in vals])
            mu, R = _circ_mean(angs, period)
            dsd.append(float(np.degrees(math.sqrt(max(-2.0 * math.log(max(R, 1e-9)),
                                                      0.0))) / (360.0 / period)))
            pa = ang[key]
            mu_p, R_p = _circ_mean(pa[:n], period)
            psd.append(float(np.degrees(math.sqrt(max(-2.0 * math.log(max(R_p, 1e-9)),
                                                      0.0))) / (360.0 / period)))
            xer.append(_fold(pa[:n] - pa[r["crystal_row"]], period))
            ora.append(_fold(pa[:n] - pa[r["crystal_row"]], period))
        donor_sd.append(float(np.mean(dsd)))
        pose_sd.append(float(np.mean(psd)))
        xe_arr = np.stack(xer, 1).mean(1)
        xtal_err.append(float(xe_arr.mean()))
        t = np.array([truth[(lig, nm)] for nm in names])
        T.append(t)
        ORA.append(np.stack(ora, 1).mean(1))
        INC.append(np.array([xe[(lig, nm)] for nm in names]))
        FEA.append(d[:n].mean(1))
        per.append({"ligand": lig, "n_donors": int(d.shape[1]),
                    "n_shared_torsions": len(bykey),
                    "donor_circ_sd_deg": float(np.mean(dsd)),
                    "pose_circ_sd_deg": float(np.mean(psd)),
                    "mean_pose_to_crystal_torsion_dev_deg": float(xe_arr.mean()),
                    "crystal_percentile_tor_cons": _percentiles(
                        d.mean(1), r["crystal_row"])[0],
                    "crystal_percentile_etkdg_tfd_min": _armb_pct(armb, lig),
                    "oracle": float(t.max()), "random_mean": float(t.mean())})
    base = float(np.mean([t.mean() for t in T]))
    out["D1_consensus_width"] = {
        "ligands": len(donor_sd),
        "donor_torsion_circular_sd_deg": float(np.mean(donor_sd)),
        "median_deg": float(np.median(donor_sd)),
        "own_20_pose_torsion_circular_sd_deg": float(np.mean(pose_sd)),
        "mean_pose_to_crystal_torsion_deviation_deg": float(np.mean(xtal_err)),
        "reading": ("the analogue of FINDING 030's D4: if the donor spread is wide "
                    "against the pose-to-crystal deviation, there is nothing to "
                    "resolve with")}

    ov, _ = FT._select([-o for o in ORA], T, np.random.default_rng(SEED))
    iv, _ = FT._select([-x for x in INC], T, np.random.default_rng(SEED))
    rhos = np.array([stats.spearmanr(o, t).statistic for o, t in zip(ORA, T)
                     if np.std(o) > 1e-12])
    out["D2_torsion_oracle"] = {
        "ligands": len(T), "random_baseline": base,
        "pool_oracle": float(np.mean([t.max() for t in T])),
        "selected_by_true_torsions": float(ov.mean()),
        "gain_vs_random": float(ov.mean() - base),
        "incumbent": float(iv.mean()),
        "within_ligand_mean_rho": float(rhos.mean()),
        "correct_sign_frac": float((rhos < 0).mean()),
        "reading": ("LEAKY BY CONSTRUCTION - it scores each pose by how far its "
                    "torsions sit from the ANSWER's torsions. It is the ceiling of "
                    "every possible internal-conformer feature on this pool")}

    # ---- D3: what arm A agrees with --------------------------------------
    pairs = {"etkdg_tfd_min": [], "etkdg_well": [], "mmff_strain": []}
    for (r, d), f in zip(rows, FEA):
        b = armb.get(r["ligand"])
        if b is None or b.get("status") != "ok":
            continue
        n = len(f)
        for k in pairs:
            v = np.asarray(b[k], float)[:n]
            if np.isfinite(v).all() and v.std() > 1e-12 and f.std() > 1e-12:
                pairs[k].append(stats.spearmanr(f, v).statistic)
    out["D3_armA_vs_prior_free"] = {
        k: {"n": len(v), "mean_within_ligand_rho": float(np.mean(v)),
            "frac_positive": float((np.array(v) > 0).mean())}
        for k, v in pairs.items() if v}

    # ---- D4: resonance exposure of the C-SYM gap -------------------------
    res_patt = [Chem.MolFromSmarts(s) for s in
                ("[CX3](=[OX1])[OX2H1,OX1-]", "[NX3,NX4+](=[OX1])[OX1,OX1-]",
                 "[SX4](=[OX1])(=[OX1])", "[PX4](=[OX1])([OX2H1,OX1-])")]
    val = pd.read_csv(DATA / "validation_ligands.csv")
    smi = dict(zip(val["id"], val["smiles"]))
    touched = total = 0
    for r, _d in rows:
        m = Chem.MolFromSmiles(smi[r["ligand"]])
        if m is None:
            continue
        flagged = set()
        for p_ in res_patt:
            if p_ is None:
                continue
            for mt in m.GetSubstructMatches(p_):
                flagged.update(mt)
        keys = {tuple(x) for mm in r["no3a"]["meta"] for x in mm["keys"]}
        total += len(keys)
        touched += sum(1 for k in keys if set(k) & flagged)
    out["D4_resonance_exposure"] = {
        "shared_torsion_keys": total, "touching_a_resonance_group": touched,
        "fraction": (touched / total if total else 0.0),
        "reading": ("C-SYM shows CanonicalRankAtoms does not merge carboxyl, nitro "
                    "or sulfonyl oxygens, so those torsions are folded at 360 deg "
                    "where 180 deg is physical. This is how much of the set it "
                    "touches")}

    pd.DataFrame(per).to_csv(DATA / "conformer_per_ligand.csv", index=False)
    p = DATA / "conformer_diagnostics.json"
    p.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print(f"\nwrote {p}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=["torsions", "armb", "answer", "evaluate",
                                      "scramble", "controls", "diagnose"])
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--xtal", action="store_true",
                    help="append the query's own crystal - DIAGNOSTIC ONLY")
    ap.add_argument("--force", action="store_true",
                    help="run the bars even if the answer-recognition gate failed, "
                         "FOR THE RECORD ONLY")
    a = ap.parse_args()
    SCRATCH.mkdir(parents=True, exist_ok=True)
    if a.stage == "torsions":
        stage_torsions(a.workers, a.xtal)
    elif a.stage == "armb":
        stage_armb(a.workers)
    elif a.stage == "answer":
        stage_answer()
    elif a.stage == "scramble":
        stage_scramble()
    elif a.stage == "diagnose":
        stage_diagnose()
    elif a.stage == "controls":
        print(json.dumps({"C_SYM": _sym_control(),
                          "C4": FT.check_frames()}, indent=2))
    else:
        stage_evaluate(a.force)


if __name__ == "__main__":
    main()
