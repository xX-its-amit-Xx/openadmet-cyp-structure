"""How much induced fit does each pose DEMAND? — FINDING 029.

FINDING 028 measured that the co-folded CYP3A4 pocket is effectively rigid (chi1 sd 0.63
deg across 87 ligands against 20.63 deg in the crystals). Every scoring experiment here,
FINDING 025's 34-term physics scorer included, was therefore computed against a receptor
with no give. This script asks the inverse question: repack the pocket around EACH pose
independently and record what the repack costs. The cost is the feature, and it is
computable from a prediction alone - no crystal, no truth, nothing unavailable at
inference.

The degeneracy that forces the design is stated in `docs/PREREG_induced_fit_demand.md`:
a predicted pose sits inside the protein the same model built around it and so does not
clash with it (median 2.571 A, 3.79% below the 2.2 A cut, FINDING 028). The receptor is
therefore DE-MOULDED first - every rotameric chi1 snapped to the nearest staggered well,
{-60, +60, 180} - and the demand is the repack cost from that state.

    python scripts/structure/induced_fit_demand.py controls
    python scripts/structure/induced_fit_demand.py features --workers 12
    python scripts/structure/induced_fit_demand.py rotate   --workers 12
    python scripts/structure/induced_fit_demand.py evaluate

Everything pre-declared is in `docs/PREREG_induced_fit_demand.md`, committed before any
score was computed, and echoed into `PREREG` below. The acceptance criteria and grid
resolutions are FINDING 028's, taken verbatim by importing them.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

# the null blocks fit 87 tiny gradient-boosting models per repetition; BLAS/OpenMP
# threading costs more than it buys at 1,740 x 11, and the reps are parallelised instead
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

POOL = Path(os.environ.get("OE_POOL", "D:/cyp_scratch/val87b_unsteered"))
SCRATCH = Path(os.environ.get(
    "IFD_SCRATCH",
    r"C:\tb\tmp\1\claude\D--Users-ashenoy00000--windsurf-OpenADMET-cyp-structure"
    r"\bd288271-2aa5-4ed4-a343-ea31e5ad8c13\scratchpad\induced_fit"))
OUT = REPO / "data" / "processed"

# FINDING 028's chemistry tables and its stage-4 acceptance criteria, imported rather
# than re-typed so that a divergence is impossible.
from side_chain_diagnosis import (  # noqa: E402
    BACKBONE, CHI1_G, ROT, SC_ATOMS, PREREG as SCD_PREREG, dihedral)

PREREG = {
    "contact_radius_A": SCD_PREREG["contact_radius_A"],     # 6.0
    "clash_cut_A": SCD_PREREG["clash_cut_A"],               # 2.2
    "self_clash_cut_A": SCD_PREREG["self_clash_cut_A"],     # 2.6
    "chi_grid_deg": SCD_PREREG["chi_grid_deg"],             # 10
    "chi3_grid_deg": SCD_PREREG["chi3_grid_deg"],           # 20
    "staggered_wells_deg": [-60.0, 60.0, 180.0],
    "blockers": [212, 215, 304, 241, 213, 119],
    "rotate_control_deg": [15.0, 90.0],
    "rotate_control_axes": 5,
    "tie_draws": 64,
    "sign": {                       # +1 = larger is a BETTER pose. Chemistry, a priori.
        "n_demand": -1, "sum_pert": -1, "max_pert": -1, "n_uncleared": -1,
        "unsolved": -1, "clash_before": +1, "clash_after": +1, "d_clash": -1,
        "rot_strain": -1, "n_skel_block": -1, "skeleton_clear": +1,
    },
    "descriptive_only": ["rot_implaus_model", "n_set", "n_movable"],
}

WELLS = np.array(PREREG["staggered_wells_deg"], float)
FEATURES = list(PREREG["sign"])

# --------------------------------------------------------------------------
# AMENDMENT, declared and committed BEFORE any LDDT-PLI was looked at
# --------------------------------------------------------------------------
# The pre-registered de-moulding rule - chi1 snapped to the nearest staggered well - turns
# out to have thin dynamic range, because the model's own chi1 angles are ALREADY close to
# canonical (mean |chi1 - nearest well| ~9.5 deg per residue). On a 120-pose probe it left
# `n_demand` at 0 for 83% of poses. That is a property of the probe, measured without
# touching a single LDDT-PLI value, and it is reported in FINDING 029 as the
# pre-registered variant's result rather than hidden.
#
# The CROSS variant asks the same question of a receptor that was NOT built around the
# pose: place pose i in the protein of sibling sample j of the SAME ligand (superposed on
# CYP3A4_RIGID_CORE, the ligand-free frame declared in targets.py), de-mould, repack, and
# average the demand over all 19 siblings. Identical feature definitions, identical signs,
# identical bars, no new thresholds. Still computable at inference: it needs only other
# samples of the same prediction.
SETS = ("BLOCKERS", "POCKET", "XBLOCKERS", "XPOCKET")
SELF_SETS = ("BLOCKERS", "POCKET")


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------

def wrap180(a):
    return ((np.asarray(a, float) + 180.0) % 360.0) - 180.0


def nearest_well_delta(chi1: float) -> tuple[float, float]:
    """(delta to the nearest staggered well, that well). Universal sp3-sp3 rule."""
    d = wrap180(WELLS - chi1)
    k = int(np.argmin(np.abs(d)))
    return float(d[k]), float(WELLS[k])


def off_well(chi1_deviation: float) -> float:
    """|chi1 - nearest staggered well| for a state `chi1_deviation` deg off a well.

    Wells are 120 deg apart, so a 120 deg flip lands on another well and costs nothing
    while a 60 deg twist is maximally strained.
    """
    return float(abs(((chi1_deviation + 60.0) % 120.0) - 60.0))


def rot_many(pts: np.ndarray, p0: np.ndarray, p1: np.ndarray,
             degs: np.ndarray) -> np.ndarray:
    """Rotate (K,n,3) `pts` about p0->p1 by each angle in `degs` -> (M,K,n,3)."""
    k = p1 - p0
    k = k / (np.linalg.norm(k) + 1e-12)
    v = pts - p0
    th = np.radians(np.asarray(degs, float))[:, None, None, None]
    kv = np.cross(np.broadcast_to(k, v.shape), v)
    kkv = (v @ k)[..., None] * k
    return p0 + v * np.cos(th) + kv * np.sin(th) + kkv * (1.0 - np.cos(th))


def rotation_matrix(axis: np.ndarray, deg: float) -> np.ndarray:
    a = np.asarray(axis, float)
    a = a / (np.linalg.norm(a) + 1e-12)
    th = np.radians(deg)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


# --------------------------------------------------------------------------
# the demand measurement — one pose, one receptor set
# --------------------------------------------------------------------------

def demand(m, lig_xyz: np.ndarray, resset: list[int], byres, resname,
           hn: np.ndarray, skel_tree, counters: dict) -> dict:
    """Repack `resset` around `lig_xyz` starting from a DE-MOULDED pocket.

    Returns the pre-registered feature dict. Nothing here touches a crystal.
    """
    from scipy.spatial import cKDTree

    cut = PREREG["clash_cut_A"]
    selfcut = PREREG["self_clash_cut_A"]
    g1 = np.arange(-180.0, 180.0, PREREG["chi_grid_deg"])
    g3 = np.arange(-180.0, 180.0, PREREG["chi3_grid_deg"])

    xyz = m.prot_xyz.copy()
    movable, demould = [], {}
    for r in resset:
        rn = resname.get(r)
        if rn not in ROT:
            continue
        at = byres.get(r, {})
        if "CA" not in at or "CB" not in at:
            continue
        g = CHI1_G.get(rn)
        if g is None or g not in at or "N" not in at:
            continue
        c1 = dihedral(xyz[at["N"]], xyz[at["CA"]], xyz[at["CB"]], xyz[at[g]])
        dlt, _well = nearest_well_delta(c1)
        demould[r] = (c1, dlt)
        movable.append(r)
        if abs(dlt) > 1e-9:
            distal = [a for a in ROT[rn][0][2] if a in at]
            idx = [at[a] for a in distal]
            xyz[idx] = rot_many(xyz[idx][None], xyz[at["CA"]], xyz[at["CB"]],
                                np.array([dlt]))[0, 0]

    # atoms of the receptor set that we are allowed to move, plus their frame
    set_idx = [i for r in movable for a, i in byres[r].items()]
    if not set_idx:
        return None
    t_lig = cKDTree(np.asarray(lig_xyz, float))

    d_set = t_lig.query(xyz[set_idx], k=1)[0]
    clash_before = float(d_set.min())
    skeleton_clear = float(skel_tree.query(np.asarray(lig_xyz, float), k=1)[0].min())

    clashing, skel_block = [], 0
    for r in movable:
        idx = list(byres[r].values())
        dd = t_lig.query(xyz[idx], k=1)[0]
        if dd.min() >= cut:
            continue
        clashing.append(r)
        names = [a for a, _i in byres[r].items()]
        hard = [j for j, a in enumerate(names) if a in BACKBONE or a == "CB"]
        if hard and dd[hard].min() < cut:
            skel_block += 1

    # the environment that does NOT move: everything but the movable side chains
    mov_atom = {i for r in clashing for a, i in byres[r].items() if a not in BACKBONE}
    fixed = np.array([i for i in range(len(xyz)) if i not in mov_atom])
    env_fixed = np.vstack([xyz[fixed], hn]) if len(hn) else xyz[fixed]
    fixed_res = np.array([m.prot_key[i][1] for i in fixed] + [-1] * len(hn))

    per_res, new_pos = [], {}
    for r in clashing:
        rn = resname[r]
        at = byres[r]
        axes = [ax for ax in ROT[rn]
                if ax[0] in at and ax[1] in at and any(x in at for x in ax[2])]
        names = [a for a in SC_ATOMS[rn] if a in at]
        if not axes or not names:
            per_res.append(dict(resnum=int(r), resname=rn, cleared=False,
                                reason="no-axis"))
            counters["no_axis"] = counters.get("no_axis", 0) + 1
            continue
        keep = ~np.isin(fixed_res, [r - 1, r, r + 1])
        t_env = cKDTree(env_fixed[keep])
        base = xyz[[at[a] for a in names]]
        pos = {a: j for j, a in enumerate(names)}
        ax_pts = [(xyz[at[a0]], xyz[at[a1]],
                   [pos[a] for a in distal if a in pos]) for a0, a1, distal in axes]

        def scan(grids):
            """Torsions applied DISTAL FIRST so the composite is exact.

            FINDING 028's scan applies chi1 before chi2 while taking chi2's CB->CG axis
            from the pre-chi1 coordinates; because chi1 moves CG that distorts the CG-CD
            bond. Applying the distal torsion first removes the distortion - the axis
            atoms of every later torsion are proximal and have not moved.
            """
            nax = len(grids)
            cur = base[None]                       # (1, n, 3)
            deltas = np.zeros((1, nax))
            for k in range(nax - 1, -1, -1):       # chi3, then chi2, then chi1
                p0, p1, sel = ax_pts[k]
                degs = grids[k]
                nxt = np.repeat(cur[None], len(degs), axis=0)       # (M,K,n,3)
                if sel:
                    nxt[:, :, sel, :] = rot_many(cur[:, sel, :], p0, p1, degs)
                cur = nxt.reshape(-1, base.shape[0], 3)
                dl = np.repeat(deltas[None], len(degs), axis=0)
                dl[..., k] = np.asarray(degs, float)[:, None]
                deltas = dl.reshape(-1, nax)
            S = cur.shape[0]
            flat = cur.reshape(-1, 3)
            ok_lig = (t_lig.query(flat, k=1)[0].reshape(S, -1).min(1) >= cut)
            free = np.inf
            mags = np.abs(deltas).max(1)
            if ok_lig.any():
                free = float(mags[ok_lig].min())
            cand = np.flatnonzero(ok_lig)
            if not len(cand):
                return (np.inf, base, np.zeros(nax)), free
            order = cand[np.argsort(mags[cand], kind="stable")]
            # evaluate the self-clash test cheapest-first; the first survivor is the
            # minimum-max|delta| solution, which is FINDING 028's acceptance rule
            for start in range(0, len(order), 256):
                blk = order[start:start + 256]
                pts = cur[blk]
                dmin = t_env.query(pts.reshape(-1, 3), k=1)[0].reshape(len(blk),
                                                                      -1).min(1)
                good = np.flatnonzero(dmin >= selfcut)
                if len(good):
                    j = blk[good[0]]
                    return (float(mags[j]), cur[j].copy(), deltas[j].copy()), free
                counters["self_clash_bound"] = counters.get("self_clash_bound", 0) + 1
            return (np.inf, base, np.zeros(nax)), free

        grids = [g1 for _ in axes[:2]]
        res, free = scan(grids)
        if not np.isfinite(res[0]) and len(axes) >= 3:
            counters["chi3_escalations"] = counters.get("chi3_escalations", 0) + 1
            res2, free2 = scan([g1, g1, g3])
            free = min(free, free2)
            if np.isfinite(res2[0]):
                res = res2
        if np.isfinite(res[0]):
            new_pos[r] = ([at[a] for a in names], res[1])
            d1 = float(res[2][0]) if len(res[2]) else 0.0
            per_res.append(dict(resnum=int(r), resname=rn, cleared=True,
                                pert=float(res[0]), off_well=off_well(d1)))
            counters["cleared"] = counters.get("cleared", 0) + 1
        else:
            per_res.append(dict(
                resnum=int(r), resname=rn, cleared=False,
                reason=("self-clash-blocks" if np.isfinite(free)
                        else "no-clearing-rotamer")))
            counters["self-clash-blocks" if np.isfinite(free)
                     else "no-clearing-rotamer"] = counters.get(
                "self-clash-blocks" if np.isfinite(free)
                else "no-clearing-rotamer", 0) + 1

    moved = xyz.copy()
    for r, (idx, pts) in new_pos.items():
        moved[idx] = pts
    clash_after = float(t_lig.query(moved[set_idx], k=1)[0].min())

    cleared = [p for p in per_res if p["cleared"]]
    rim = 0.0
    for r in movable:
        rim += abs(demould[r][1])
    return dict(
        n_set=len(resset), n_movable=len(movable),
        n_demand=len(clashing),
        sum_pert=float(sum(p["pert"] for p in cleared)),
        max_pert=float(max([p["pert"] for p in cleared], default=0.0)),
        n_uncleared=int(sum(1 for p in per_res if not p["cleared"])),
        unsolved=int(clash_after < cut),
        clash_before=clash_before, clash_after=clash_after,
        d_clash=float(clash_after - clash_before),
        rot_strain=float(sum(p["off_well"] for p in cleared)),
        n_skel_block=int(skel_block),
        skeleton_clear=skeleton_clear,
        rot_implaus_model=float(rim))


def _prepare(m):
    """Per-structure indices reused by every receptor set and every rotation."""
    from scipy.spatial import cKDTree
    byres: dict[int, dict[str, int]] = {}
    for i, (_c, r, a) in enumerate(m.prot_key):
        byres.setdefault(r, {})[a] = i
    resname = {r: v for (_c, r), v in m.prot_res.items()}
    hm = np.asarray(m.heme_xyz, float)
    he = [str(e).upper() for e in m.heme_elem]
    hn = hm[[i for i, e in enumerate(he) if e != "FE"]] if len(hm) else hm
    skel = [i for i, (_c, _r, a) in enumerate(m.prot_key)
            if a in BACKBONE or a == "CB"]
    sk = np.vstack([m.prot_xyz[skel], hn]) if len(hn) else m.prot_xyz[skel]
    return byres, resname, hn, cKDTree(sk)


def pocket_set(m, lig_xyz, byres, resname) -> list[int]:
    """Rotameric residues within 6.0 A of THIS pose. Computed from the prediction."""
    from scipy.spatial import cKDTree
    d = cKDTree(np.asarray(lig_xyz, float)).query(m.prot_xyz, k=1)[0]
    near = {m.prot_key[i][1] for i in np.flatnonzero(d < PREREG["contact_radius_A"])}
    return sorted(r for r in near if resname.get(r) in ROT and r in byres)


# --------------------------------------------------------------------------
# stage: features
# --------------------------------------------------------------------------

def _pose_files(lig: str):
    job = POOL / f"{lig}__unsteered__s1"
    return sorted(job.glob("input_model_*.cif"),
                  key=lambda p: int(p.stem.rsplit("_", 1)[1]))


def features_ligand(lig: str) -> dict:
    from cypstruct import pose as P
    out = SCRATCH / f"feat_{lig}.json"
    if out.exists():
        return {"ligand": lig, "status": "cached"}
    rows, counters = [], {}
    for cif in _pose_files(lig):
        m = P.load_structure(cif)
        byres, resname, hn, skt = _prepare(m)
        L = np.asarray(m.lig_xyz, float)
        if not len(L):
            continue
        rec = {"ligand": lig, "sample": cif.stem}
        ps = pocket_set(m, L, byres, resname)
        for nm, rs in (("BLOCKERS", [r for r in PREREG["blockers"] if r in byres]),
                       ("POCKET", ps)):
            d = demand(m, L, rs, byres, resname, hn, skt, counters)
            if d is None:
                continue
            for k, v in d.items():
                rec[f"{nm}_{k}"] = v
        rows.append(rec)
    out.write_text(json.dumps({"rows": rows, "counters": counters}))
    return {"ligand": lig, "status": "ok", "poses": len(rows)}


def cross_ligand(lig: str) -> dict:
    """Demand of pose i inside the receptors of its 19 siblings (the AMENDMENT)."""
    from cypstruct import pose as P
    from cypstruct.targets import CYP3A4_RIGID_CORE

    core = set(CYP3A4_RIGID_CORE)
    out = SCRATCH / f"cross_{lig}.json"
    if out.exists():
        return {"ligand": lig, "status": "cached"}
    cifs = _pose_files(lig)
    ms = [P.load_structure(c) for c in cifs]
    keep = [k for k, m in enumerate(ms) if len(m.lig_xyz)]
    if len(keep) < 2:
        return {"ligand": lig, "status": "too-few-poses"}
    prep = [_prepare(m) for m in ms]
    cas = [{k: v for (_c, k), v in m.ca().items()} for m in ms]

    acc: dict[int, dict[str, list]] = {k: {} for k in keep}
    counters: dict = {}
    for i in keep:
        for j in keep:
            if i == j:
                continue
            sh = sorted(set(cas[i]) & set(cas[j]) & core)
            if len(sh) < 50:
                continue
            R, t, _f = P.kabsch(np.array([cas[i][k] for k in sh]),
                                np.array([cas[j][k] for k in sh]))
            L = np.asarray(ms[i].lig_xyz, float) @ R.T + t
            byres, resname, hn, skt = prep[j]
            for nm, rs in (("XBLOCKERS",
                            [r for r in PREREG["blockers"] if r in byres]),
                           ("XPOCKET", pocket_set(ms[j], L, byres, resname))):
                d = demand(ms[j], L, rs, byres, resname, hn, skt, counters)
                if d is None:
                    continue
                for k2, v in d.items():
                    acc[i].setdefault(f"{nm}_{k2}", []).append(float(v))
    rows = []
    for i in keep:
        if not acc[i]:
            continue
        rec = {"ligand": lig, "sample": cifs[i].stem,
               "n_sibling_receptors": len(acc[i].get("XPOCKET_n_demand", []))}
        for k2, v in acc[i].items():
            rec[k2] = float(np.mean(v))
        rows.append(rec)
    out.write_text(json.dumps({"rows": rows, "counters": counters}))
    return {"ligand": lig, "status": "ok", "poses": len(rows)}


def rotate_ligand(lig: str) -> dict:
    """The wrong-pose control: the SAME pose, rigidly rotated about its own centroid."""
    from cypstruct import pose as P
    out = SCRATCH / f"rot_{lig}.json"
    if out.exists():
        return {"ligand": lig, "status": "cached"}
    rng = np.random.default_rng(abs(hash(lig)) % (2 ** 31))
    cifs = _pose_files(lig)
    if not cifs:
        return {"ligand": lig, "status": "no-poses"}
    m = P.load_structure(cifs[0])
    byres, resname, hn, skt = _prepare(m)
    L = np.asarray(m.lig_xyz, float)
    if not len(L):
        return {"ligand": lig, "status": "no-ligand"}
    counters: dict = {}
    rows = []
    cen = L.mean(axis=0)

    # two receptors: the pose's OWN protein (the pre-registered control) and a SIBLING
    # sample's protein, which is the receptor the CROSS amendment actually scores in
    frames = [("self", m, byres, resname, hn, skt, L)]
    if len(cifs) > 1:
        from cypstruct.targets import CYP3A4_RIGID_CORE
        m2 = P.load_structure(cifs[1])
        b2, r2, h2, s2 = _prepare(m2)
        ca1 = {k: v for (_c, k), v in m.ca().items()}
        ca2 = {k: v for (_c, k), v in m2.ca().items()}
        sh = sorted(set(ca1) & set(ca2) & set(CYP3A4_RIGID_CORE))
        if len(sh) >= 50:
            R, t, _f = P.kabsch(np.array([ca1[k] for k in sh]),
                                np.array([ca2[k] for k in sh]))
            frames.append(("sibling", m2, b2, r2, h2, s2, L @ R.T + t))

    for fname, mm, br, rn2, hh, sk, L0 in frames:
        c0 = L0.mean(axis=0)
        d0 = demand(mm, L0, pocket_set(mm, L0, br, rn2), br, rn2, hh, sk, counters)
        if d0:
            rows.append(dict(ligand=lig, frame=fname, kind="original", deg=0.0,
                             rep=0, **d0))
        for deg in PREREG["rotate_control_deg"]:
            for rep in range(PREREG["rotate_control_axes"]):
                ax = rng.normal(size=3)
                Lr = (L0 - c0) @ rotation_matrix(ax, deg).T + c0
                d = demand(mm, Lr, pocket_set(mm, Lr, br, rn2), br, rn2, hh, sk,
                           counters)
                if d is None:
                    continue
                rows.append(dict(ligand=lig, frame=fname, kind=f"rot{int(deg)}",
                                 deg=float(deg), rep=rep, **d))
    out.write_text(json.dumps({"rows": rows, "counters": counters}))
    return {"ligand": lig, "status": "ok", "rows": len(rows)}


def _run(fn, ligs, workers):
    todo = list(ligs)
    print(f"{len(todo)} ligands, {workers} workers", flush=True)
    if workers <= 1:
        for l in todo:
            print(fn(l), flush=True)
        return
    import multiprocessing as mp
    with mp.Pool(workers) as pool:
        for r in pool.imap_unordered(fn, todo):
            print(r, flush=True)


def ligand_list():
    return sorted({p.name.split("__")[0] for p in POOL.iterdir() if p.is_dir()})


# --------------------------------------------------------------------------
# stage: controls
# --------------------------------------------------------------------------

def cmd_controls() -> None:
    import pandas as pd
    from cypstruct import pose as P
    from cypstruct import xengine as X
    from cypstruct.targets import fetch_sequences

    three = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
             "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
             "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V"}
    seq = fetch_sequences()["cyp3a4"]

    ref = X.load_reference(OUT / "reference_set_cyp3a4.npz")
    depth = X.reference_depth(ref)
    shipped = pd.read_csv(OUT / "xeng_val87b.csv")
    ship_map = {(r.ligand, r.sample): r.xeng for r in shipped.itertuples()}

    offsets, idents, blockers_ok, diffs, missing = {}, [], 0, [], 0
    n = 0
    for lig in ligand_list():
        for cif in _pose_files(lig):
            m = P.load_structure(cif)
            rn = {r: v for (_c, r), v in m.prot_res.items()}
            nums = sorted(rn)
            best, best_id = 0, -1.0
            for off in range(-30, 31):
                sh = [k for k in nums if 1 <= k + off <= len(seq)]
                if len(sh) < 30:
                    continue
                ident = sum(three.get(rn[k], "X") == seq[k + off - 1]
                            for k in sh) / len(sh)
                if ident > best_id:
                    best, best_id = off, ident
            offsets[best] = offsets.get(best, 0) + 1
            idents.append(best_id)
            blockers_ok += int(all(
                three.get(rn.get(r, ""), "X") == seq[r - 1]
                for r in PREREG["blockers"]))
            v = X.in_heme_frame(m)
            if v is None or lig not in ref:
                missing += 1
            else:
                mine = X.xeng_score(v, ref[lig])
                have = ship_map.get((lig, cif.stem))
                diffs.append(abs(mine - have) if have is not None else np.nan)
            n += 1

    scored = pd.read_csv(OUT / "poses_scored_val87b.csv")
    uns = scored[scored.arm == "unsteered"]
    diffs = np.array(diffs, float)
    rep = {
        "n_poses": n,
        "C1_numbering_offset_vs_canonical_sequence": {
            "offsets": {str(k): int(v) for k, v in offsets.items()},
            "min_residue_identity": float(np.min(idents)),
            "median_residue_identity": float(np.median(idents)),
            "all_six_blockers_match_sequence": int(blockers_ok),
            "fires": bool(len(offsets) >= 1)},
        "C2_exact_zero_lddt": {
            "n_unsteered_rows": int(len(uns)),
            "n_lddt_pli_exactly_zero": int((uns.lddt_pli == 0).sum()),
            "n_lddt_pli_null": int(uns.lddt_pli.isna().sum()),
            "n_unmapped": int((~uns.mapped.astype(bool)).sum())},
        "C3_incumbent_reproduction": {
            "reference_depth": {k: v for k, v in depth.items() if k != "note"},
            "n_compared": int(np.isfinite(diffs).sum()),
            "n_missing_reference": int(missing),
            "max_abs_diff": float(np.nanmax(diffs)) if len(diffs) else None,
            "median_abs_diff": float(np.nanmedian(diffs)) if len(diffs) else None,
            "reproduces_below_1e-6": bool(np.nanmax(diffs) < 1e-6)},
    }
    (OUT / "induced_fit_controls.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))


# --------------------------------------------------------------------------
# stage: evaluate
# --------------------------------------------------------------------------

def _load_features():
    import pandas as pd
    counters: dict = {}

    def grab(pattern):
        rows = []
        for f in sorted(SCRATCH.glob(pattern)):
            d = json.loads(f.read_text())
            rows.extend(d["rows"])
            for k, v in d["counters"].items():
                counters[k] = counters.get(k, 0) + v
        return pd.DataFrame(rows)

    self_df = grab("feat_*.json")
    cross_df = grab("cross_*.json")
    if len(cross_df):
        self_df = self_df.merge(cross_df, on=["ligand", "sample"], how="left")
    return self_df, counters


def cmd_rotate_report() -> None:
    import pandas as pd
    from scipy import stats
    rows, counters = [], {}
    for f in sorted(SCRATCH.glob("rot_*.json")):
        d = json.loads(f.read_text())
        rows.extend(d["rows"])
        for k, v in d["counters"].items():
            counters[k] = counters.get(k, 0) + v
    d = pd.DataFrame(rows)
    rep = {"n_ligands": int(d.ligand.nunique()), "counters": counters, "strata": {}}
    keys = ["n_demand", "clash_before", "n_uncleared", "skeleton_clear"]
    direction = {"n_demand": "greater", "clash_before": "less",
                 "n_uncleared": "greater", "skeleton_clear": "less"}
    for fr in sorted(d.frame.unique()):
        df = d[d.frame == fr]
        orig = df[df.kind == "original"].set_index("ligand")
        for kind in [k for k in df.kind.unique() if k != "original"]:
            sub = df[df.kind == kind]
            blk = {"n": int(len(sub))}
            n_pass = 0
            for k in keys:
                a = sub[k].values
                b = orig.loc[sub.ligand, k].values
                mw = stats.mannwhitneyu(a, b, alternative=direction[k])
                sign_ok = float((a > b).mean() if direction[k] == "greater"
                                else (a < b).mean())
                ok = bool(mw.pvalue < 0.01 and sign_ok >= 0.60)
                n_pass += int(ok)
                blk[k] = {"original_median": float(np.median(b)),
                          "rotated_median": float(np.median(a)),
                          "expected": direction[k],
                          "mannwhitney_p": float(mw.pvalue),
                          "frac_paired_sign_correct": round(sign_ok, 3),
                          "passes": ok}
            blk["n_criteria_passing"] = n_pass
            blk["CONTROL_FIRES"] = bool(n_pass >= 2)
            rep["strata"][f"{fr}:{kind}"] = blk
    rep["verdict"] = ("FIRES" if rep["strata"].get("self:rot90", {}).get("CONTROL_FIRES")
                      else "DOES NOT FIRE")
    (OUT / "induced_fit_rotation_control.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))


def _loo(X, y, names, seed: int = 0):
    """Leave-one-LIGAND-out out-of-fold prediction. FINDING 025's fit, unchanged."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    p = np.empty(len(y))
    for nm in np.unique(names):
        mk = names == nm
        g = HistGradientBoostingRegressor(max_iter=200,
                                          random_state=seed).fit(X[~mk], y[~mk])
        p[mk] = g.predict(X[mk])
    return p


def _null_pred(args):
    kind, rep, X, y, names, index = args
    rg = np.random.default_rng((3000 if kind == "shuffle" else 7000) + rep)
    if kind == "shuffle":
        Xs = X.copy()
        for nm in index:
            idx = index[nm]
            Xs[idx] = X[rg.permutation(idx)]
    else:
        Xs = rg.normal(size=X.shape)
    return kind, _loo(Xs, y, names)


def cmd_evaluate(n_null: int = 200, n_rand: int = 4000, boot: int = 10000) -> None:
    import pandas as pd
    from scipy import stats
    from sklearn.ensemble import HistGradientBoostingRegressor

    RNG = np.random.default_rng(29)
    feat, counters = _load_features()
    truth = pd.read_csv(OUT / "poses_scored_val87b.csv")
    truth = truth[truth.arm == "unsteered"][["ligand", "sample", "lddt_pli"]]
    xe = pd.read_csv(OUT / "xeng_val87b.csv")
    d = (feat.merge(truth, on=["ligand", "sample"], how="inner")
             .merge(xe, on=["ligand", "sample"], how="left")
             .sort_values(["ligand", "sample"]).reset_index(drop=True))
    d.to_csv(OUT / "induced_fit_features_val87b.csv", index=False)

    groups = list(d.groupby("ligand"))
    index = {nm: g.index.to_numpy() for nm, g in groups}
    y = d.lddt_pli.values
    names = d.ligand.values
    rnd_lig = np.array([g.lddt_pli.mean() for _, g in groups])
    orc_lig = np.array([g.lddt_pli.max() for _, g in groups])
    rnd, oracle = float(rnd_lig.mean()), float(orc_lig.mean())

    def per_lig(score, draws=PREREG["tie_draws"]):
        out = np.zeros(len(groups))
        for _ in range(draws):
            for i, (nm, g) in enumerate(groups):
                v = np.where(np.isfinite(score[index[nm]]), score[index[nm]], -np.inf)
                out[i] += g.lddt_pli.values[RNG.choice(np.flatnonzero(v == v.max()))]
        return out / draws

    def rho_block(score):
        rs = []
        for nm, g in groups:
            s = score[index[nm]]
            ok = np.isfinite(s)
            if ok.sum() < 3 or np.nanstd(s[ok]) == 0:
                continue
            r = stats.spearmanr(s[ok], g.lddt_pli.values[ok]).statistic
            if r == r:
                rs.append(r)
        rs = np.array(rs)
        return rs

    def summarise(score, label):
        per = per_lig(score)
        rs = rho_block(score)
        try:
            p = float(stats.wilcoxon(per, rnd_lig).pvalue)
        except ValueError:
            p = float("nan")
        return {"term": label, "selected": round(float(per.mean()), 4),
                "gain": round(float(per.mean() - rnd), 4),
                "rho": round(float(rs.mean()) if len(rs) else float("nan"), 3),
                "frac_rho_correct_sign": round(float((rs < 0).mean())
                                               if len(rs) else float("nan"), 3),
                "n_rho": int(len(rs)),
                "frac_ligands_beating_random": round(float((per > rnd_lig).mean()), 3),
                "p_wilcoxon": p}, per

    # ---- within-ligand constancy, BEFORE any term is tested -------------
    const = {}
    for s in SETS:
        for f in FEATURES + PREREG["descriptive_only"]:
            c = f"{s}_{f}"
            if c not in d.columns:
                continue
            frac = float(d.groupby("ligand")[c].nunique().le(1).mean())
            cv = d.groupby("ligand")[c].apply(
                lambda v: (v.std() / abs(v.mean())) if abs(v.mean()) > 1e-9 else 0.0)
            const[c] = {"frac_ligands_constant": round(frac, 3),
                        "median_within_ligand_cv": round(float(cv.median()), 4)}

    # ---- single terms, sign fixed a priori -------------------------------
    singles = []
    for s in SETS:
        for f, sg in PREREG["sign"].items():
            c = f"{s}_{f}"
            if c not in d.columns:
                continue
            if const[c]["frac_ligands_constant"] > 0.5:
                singles.append({"term": c, "sign": sg, "skipped":
                                "constant within >50% of ligands",
                                "frac_ligands_constant":
                                    const[c]["frac_ligands_constant"]})
                continue
            r, _ = summarise(sg * d[c].astype(float).values, c)
            r["sign"] = sg
            singles.append(r)
    # rho sign convention: `rho` above is computed on the SIGNED score, so a correct
    # term has a POSITIVE rho; report it that way rather than by the raw column.
    for r in singles:
        if "rho" in r:
            r["frac_rho_correct_sign"] = round(
                float((rho_block(r["sign"] * d[r["term"]].astype(float).values) > 0
                       ).mean()), 3)

    # ---- fitted ensembles, leave-one-LIGAND-out --------------------------
    ens, ens_per, ens_pred = {}, {}, {}
    for s in SETS:
        cols = [f"{s}_{f}" for f in FEATURES if f"{s}_{f}" in d.columns]
        if not cols:
            continue
        X = d[cols].astype(float).values
        pred = _loo(X, y, names)
        r, per = summarise(pred, f"ENSEMBLE {s} ({len(cols)} terms)")
        r["frac_rho_correct_sign"] = round(float((rho_block(pred) > 0).mean()), 3)
        ens[s], ens_per[s], ens_pred[s] = r, per, pred
    cols_all = [f"{s}_{f}" for s in SETS for f in FEATURES
                if f"{s}_{f}" in d.columns]
    Xall = d[cols_all].astype(float).values
    pred_all = _loo(Xall, y, names)
    r_all, per_all = summarise(pred_all, f"ENSEMBLE ALL ({len(cols_all)} terms)")
    r_all["frac_rho_correct_sign"] = round(float((rho_block(pred_all) > 0).mean()), 3)
    ens["ALL"], ens_per["ALL"], ens_pred["ALL"] = r_all, per_all, pred_all

    # ---- the incumbent, on identical poses -------------------------------
    inc_score = -d.xeng.fillna(np.inf).values
    inc_sum, inc_per = summarise(inc_score, "INCUMBENT xengine.select()")

    # ---- N1 random selection null ---------------------------------------
    n1 = np.array([per_lig(RNG.normal(size=len(d)), draws=1).mean() - rnd
                   for _ in range(n_rand)])
    # ---- N2 within-ligand shuffle / N3 matched-dimensionality Gaussian ---
    best_set = max([s for s in ens if s != "ALL"], key=lambda s: ens[s]["gain"])
    cols_b = [f"{best_set}_{f}" for f in FEATURES if f"{best_set}_{f}" in d.columns]
    Xb = d[cols_b].astype(float).values
    jobs = ([("shuffle", rep, Xb, y, names, index) for rep in range(n_null)] +
            [("gauss", rep, Xb, y, names, index) for rep in range(n_null)])
    import multiprocessing as mp
    with mp.Pool(12) as pool:
        preds = pool.map(_null_pred, jobs)
    n2 = np.array([per_lig(p, draws=8).mean() - rnd
                   for k, p in preds if k == "shuffle"])
    n3 = np.array([per_lig(p, draws=8).mean() - rnd
                   for k, p in preds if k == "gauss"])

    nulls = {
        "N1_random_selection": {"draws": n_rand, "p95": round(float(np.percentile(n1, 95)), 4),
                                "p99": round(float(np.percentile(n1, 99)), 4),
                                "max": round(float(n1.max()), 4)},
        "N2_within_ligand_shuffle": {"reps": n_null, "on": best_set,
                                     "mean": round(float(n2.mean()), 4),
                                     "sd": round(float(n2.std()), 4),
                                     "p95": round(float(np.percentile(n2, 95)), 4),
                                     "p99": round(float(np.percentile(n2, 99)), 4),
                                     "max": round(float(n2.max()), 4)},
        "N3_matched_dimension_gaussian": {"reps": n_null, "k": Xb.shape[1],
                                          "mean": round(float(n3.mean()), 4),
                                          "sd": round(float(n3.std()), 4),
                                          "p95": round(float(np.percentile(n3, 95)), 4),
                                          "p99": round(float(np.percentile(n3, 99)), 4),
                                          "max": round(float(n3.max()), 4)},
    }

    # ---- BAR 2, paired against the incumbent -----------------------------
    def paired(per, label):
        diff = per - inc_per
        bs = np.array([diff[RNG.integers(0, len(diff), len(diff))].mean()
                       for _ in range(boot)])
        nz = diff[np.abs(diff) > 1e-9]
        try:
            w = float(stats.wilcoxon(nz).pvalue) if len(nz) > 5 else float("nan")
        except ValueError:
            w = float("nan")
        # ties: ligands where both selectors' deterministic top-1 is the same pose
        return {"selector": label,
                "mean_difference": round(float(diff.mean()), 4),
                "bootstrap_ci95": [round(float(np.percentile(bs, 2.5)), 4),
                                   round(float(np.percentile(bs, 97.5)), 4)],
                "wilcoxon_p_untied": w,
                "n_ligands": int(len(diff)),
                "n_tied_value": int((np.abs(diff) <= 1e-9).sum()),
                "n_ligands_better": int((diff > 1e-9).sum()),
                "n_ligands_worse": int((diff < -1e-9).sum()),
                "passes": bool(diff.mean() > 0
                               and np.percentile(bs, 2.5) > 0
                               and w == w and w < 0.05)}

    # deterministic top-1 agreement (a value tie can hide a different pose)
    def top1(score):
        out = {}
        for nm, g in groups:
            v = np.where(np.isfinite(score[index[nm]]), score[index[nm]], -np.inf)
            out[nm] = g["sample"].values[int(np.argmax(v))]
        return out
    t_inc = top1(inc_score)
    agree = {s: int(sum(top1(ens_pred[s])[k] == t_inc[k] for k in t_inc))
             for s in ens_pred}

    bar2 = {s: paired(ens_per[s], f"ENSEMBLE {s}") for s in ens_per}
    for s in bar2:
        bar2[s]["n_same_pose_picked"] = agree[s]

    # ---- BAR 3, complementarity -----------------------------------------
    bar3 = {}
    for s in ens_per:
        gain = ens_per[s] - rnd_lig
        short = orc_lig - inc_per
        rp = stats.pearsonr(gain, short)
        rs = stats.spearmanr(gain, short)
        bar3[s] = {"pearson_r": round(float(rp.statistic), 3),
                   "pearson_p": float(rp.pvalue),
                   "spearman_r": round(float(rs.statistic), 3),
                   "spearman_p": float(rs.pvalue),
                   "corr_with_incumbent_outcome": round(
                       float(np.corrcoef(ens_per[s], inc_per)[0, 1]), 3),
                   "passes": bool(rp.statistic > 0 and rp.pvalue < 0.05)}

    best = max(ens, key=lambda s: ens[s]["gain"])
    bar1 = bool(ens[best]["gain"] > nulls["N2_within_ligand_shuffle"]["p95"]
                and ens[best]["gain"] > nulls["N3_matched_dimension_gaussian"]["p95"])
    best_single = max([r for r in singles if "gain" in r],
                      key=lambda r: r["gain"], default=None)
    bar1_single = bool(best_single
                       and best_single["gain"] > nulls["N1_random_selection"]["p95"])

    verdict = ("REFUTED" if not (bar1 or bar1_single)
               else "MEASURED" if not bar2[best]["passes"]
               else "SHIPS" if bar3[best]["passes"] else "MEASURED")

    rep = {
        "prereg": PREREG,
        "pool": {"n_poses": int(len(d)), "n_ligands": int(len(groups)),
                 "oracle": round(oracle, 4), "random": round(rnd, 4)},
        "criteria_fired": counters,
        "within_ligand_constancy": const,
        "single_terms": sorted([r for r in singles if "gain" in r],
                               key=lambda r: -r["gain"]) +
                        [r for r in singles if "gain" not in r],
        "ensembles": ens,
        "incumbent": inc_sum,
        "nulls": nulls,
        "BAR1_beats_random": {"best_single": best_single,
                              "single_clears_N1_p95": bar1_single,
                              "best_ensemble": best,
                              "ensemble_clears_N2_and_N3_p95": bar1},
        "BAR2_paired_vs_incumbent": bar2,
        "BAR3_complementarity": bar3,
        "VERDICT": verdict,
    }
    (OUT / "induced_fit_selector_val87b.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps({k: v for k, v in rep.items() if k != "prereg"}, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["controls", "features", "cross", "rotate",
                                      "rotate_report", "evaluate"])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--null-reps", type=int, default=200)
    a = ap.parse_args()
    SCRATCH.mkdir(parents=True, exist_ok=True)
    if a.stage == "controls":
        cmd_controls()
    elif a.stage == "features":
        ligs = ligand_list()[:a.limit] if a.limit else ligand_list()
        _run(features_ligand, [l for l in ligs
                               if not (SCRATCH / f"feat_{l}.json").exists()], a.workers)
    elif a.stage == "cross":
        ligs = ligand_list()[:a.limit] if a.limit else ligand_list()
        _run(cross_ligand, [l for l in ligs
                            if not (SCRATCH / f"cross_{l}.json").exists()], a.workers)
    elif a.stage == "rotate":
        ligs = ligand_list()[:a.limit] if a.limit else ligand_list()
        _run(rotate_ligand, [l for l in ligs
                             if not (SCRATCH / f"rot_{l}.json").exists()], a.workers)
    elif a.stage == "rotate_report":
        cmd_rotate_report()
    else:
        cmd_evaluate(n_null=a.null_reps)


if __name__ == "__main__":
    main()
