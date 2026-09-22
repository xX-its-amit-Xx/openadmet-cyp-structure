"""Expand the CYP3A4 pool by rotating each ligand INSIDE the model's own protein.

Pre-registered in `docs/PREREG_orientation_expansion.md`; nothing here may be tuned
after a scored number exists. Written up as FINDING 027.

**The one lever FINDING 024 leaves standing.** On CYP3A4 the co-folded pocket is right
(0.73 A), the ligand centroid is right (1.07 A, isotropic), and only the orientation is
wrong (30.1 deg against 6.1 deg for the rest of the P450 family). If that is true, the
poses the model never produces are reachable by a rigid rotation of the ligand about its
own centroid, in the protein the model already built, with **no new inference**.

Two stages, both CPU-only and resumable per ligand:

    python scripts/structure/orientation_expansion.py expand --workers 12
    python scripts/structure/orientation_expansion.py analyse

`expand` writes one `<LIG>.npz` per ligand into the scratch directory (~45 MB in total for
87 ligands x 20 poses x 1024 orientations) and skips ligands already done. `analyse` reads
them back and produces `data/processed/orientation_expansion.json`.

**The trap this script is built around.** An oracle rises mechanically with pool size, so
every table below reports the oracle *and* what the shipped selector actually reaches, and
the size-matched unfiltered null is computed from the same arrays so the filter cannot
take credit for the expansion's size. FINDING 013 is the precedent: a union pool added
+0.0375 of oracle that selection could not reach.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

POOL = Path(os.environ.get("OE_POOL", "D:/cyp_scratch/val87b_unsteered"))
SCRATCH = Path(os.environ.get(
    "OE_SCRATCH",
    r"C:\tb\tmp\1\claude\D--Users-ashenoy00000--windsurf-OpenADMET-cyp-structure"
    r"\bd288271-2aa5-4ed4-a343-ea31e5ad8c13\scratchpad\orientation_expansion"))

N_GRID = 1024                 # pre-registered
CLASH_CUT = 2.2               # pre-registered, A, heavy-atom, Fe excluded
CLASH_CUT_SOFT = 3.0          # pre-registered sensitivity only
NULL_SEED = 20260922          # pre-registered


# --------------------------------------------------------------------------
# the SO(3) grid
# --------------------------------------------------------------------------

def super_fibonacci(n: int = N_GRID) -> np.ndarray:
    """Deterministic near-uniform quaternions on SO(3) (Alexa, CVPR 2022).

    Deterministic on purpose: the identical grid is applied to every one of the 1,740
    poses, so no pose can be advantaged by a luckier draw, and the experiment has no RNG
    in it anywhere except the pre-registered size-matched null.
    """
    phi = np.sqrt(2.0)
    psi = 1.533751168755204288118041
    i = np.arange(n) + 0.5
    s = i / n
    r, R = np.sqrt(s), np.sqrt(1.0 - s)
    a = 2.0 * np.pi * i / phi
    b = 2.0 * np.pi * i / psi
    return np.column_stack([r * np.sin(a), r * np.cos(a), R * np.sin(b), R * np.cos(b)])


def quat_to_mat(q: np.ndarray) -> np.ndarray:
    """(N,4) unit quaternions -> (N,3,3) rotation matrices."""
    q = q / np.linalg.norm(q, axis=1, keepdims=True)
    w, x, y, z = q.T
    return np.stack([
        np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], 1),
        np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], 1),
        np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], 1),
    ], 1)


def grid_spacing(q: np.ndarray) -> dict:
    """Realised nearest-neighbour geodesic angle of the grid, in degrees."""
    c = np.abs(q @ q.T)
    np.fill_diagonal(c, 0.0)
    nn = 2.0 * np.degrees(np.arccos(np.clip(c.max(axis=1), -1, 1)))
    return {"n": int(len(q)), "nn_mean_deg": round(float(nn.mean()), 2),
            "nn_median_deg": round(float(np.median(nn)), 2),
            "nn_max_deg": round(float(nn.max()), 2)}


# --------------------------------------------------------------------------
# reference handling
# --------------------------------------------------------------------------

def load_reference(pdb_id: str, ligand_code: str):
    """Deposited structure, ONE chain, ligand pinned by CCD code."""
    import gemmi
    from cypstruct import pose as P
    from cypstruct.targets import fetch_cif
    cif = fetch_cif(pdb_id)
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    for chain in st[0]:
        if any(r.name.strip().upper() == ligand_code.upper() for r in chain):
            return P.load_structure(cif, ligand_code=ligand_code, assembly_chain=chain.name)
    return None


def renumber_to_reference(pred, ref):
    """FINDING 021. Shift the prediction's residue numbering onto the crystal's.

    Boltz numbers 1..N; a crystal keeps author numbering. `lddt_pli` pairs protein atoms
    by residue NUMBER, so a mismatch compares every contact against the wrong residue and
    the score reads **exactly 0.0** while looking like a target the model cannot fold.
    The offset maximises three-letter residue-NAME agreement, which cannot manufacture
    agreement where there is none.
    """
    pm = {num: name for (_c, num), name in pred.prot_res.items()}
    rm = {num: name for (_c, num), name in ref.prot_res.items()}
    if len(pm) < 30 or len(rm) < 30:
        return pred, 0, 0.0
    best, best_id = 0, -1.0
    for off in range(int(min(rm) - max(pm)), int(max(rm) - min(pm)) + 1):
        shared = [n for n in pm if n + off in rm]
        if len(shared) < 30:
            continue
        ident = sum(pm[n] == rm[n + off] for n in shared) / len(shared)
        if ident > best_id:
            best, best_id = off, ident
    if best_id < 0:
        return pred, 0, 0.0
    if best != 0:
        pred.prot_key = [(c, n + best, a) for (c, n, a) in pred.prot_key]
        pred.prot_res = {(c, n + best): v for (c, n), v in pred.prot_res.items()}
    return pred, best, best_id


# --------------------------------------------------------------------------
# one ligand
# --------------------------------------------------------------------------

def run_ligand(args) -> dict:
    lig_id, pdb, smiles = args
    import pandas as pd  # noqa: F401  (imported for the worker's side effects only)
    from scipy.spatial import cKDTree

    from cypstruct import pose as P
    from cypstruct import xengine as X
    from cypstruct.qmscore import geometry as G

    out = SCRATCH / f"{lig_id}.npz"
    if out.exists():
        return {"ligand": lig_id, "status": "cached"}

    ref = load_reference(pdb, lig_id)
    if ref is None or len(ref.lig_xyz) == 0:
        return {"ligand": lig_id, "status": "no-reference"}

    refset = X.load_reference(REPO / "data" / "processed" / "reference_set_cyp3a4.npz")
    xrefs = refset.get(lig_id, [])

    Q = super_fibonacci(N_GRID)
    Rm = quat_to_mat(Q)                                     # (G,3,3)

    job = POOL / f"{lig_id}__unsteered__s1"
    cifs = sorted(job.glob("input_model_*.cif"),
                  key=lambda p: int(p.stem.rsplit("_", 1)[1]))
    if not cifs:
        return {"ligand": lig_id, "status": "no-poses"}

    ref_pocket = P.pocket_residues_from_structure(ref, radius=8.0)
    n_lig_ref = len(ref.lig_xyz)

    samples, meta = [], []
    L_ld = np.zeros((len(cifs), N_GRID), np.float32)
    L_rm = np.zeros((len(cifs), N_GRID), np.float32)
    L_xe = np.zeros((len(cifs), N_GRID), np.float32)
    L_cl = np.zeros((len(cifs), N_GRID), np.float32)
    L_fd = np.zeros((len(cifs), N_GRID), np.float32)
    L_co = np.zeros((len(cifs), N_GRID), bool)

    for si, cif in enumerate(cifs):
        model = P.load_structure(cif)
        model, off, ident = renumber_to_reference(model, ref)
        perm = P.best_ligand_mapping(smiles, model, ref)
        mapped = perm is not None
        if perm is None:
            perm = np.arange(min(len(model.lig_xyz), n_lig_ref))

        L0 = np.asarray(model.lig_xyz, float)
        c = L0.mean(axis=0)
        # (G, n, 3) — every orientation of this ligand about its own centroid
        cand = np.einsum("gij,nj->gni", Rm, L0 - c) + c

        # ---- LDDT-PLI, masked pairs only ---------------------------------
        ref_idx = {(r, a): i for i, (_c, r, a) in enumerate(ref.prot_key)}
        mod_idx = {(r, a): i for i, (_c, r, a) in enumerate(model.prot_key)}
        shared = sorted(set(ref_idx) & set(mod_idx))
        ri = np.array([ref_idx[k] for k in shared])
        mi = np.array([mod_idx[k] for k in shared])
        lr = ref.lig_xyz[perm]
        dref = np.linalg.norm(lr[:, None, :] - ref.prot_xyz[ri][None, :, :], axis=2)
        mask = dref < 6.0
        ki, ji = np.nonzero(mask)
        dref_p = dref[ki, ji]
        Pm = model.prot_xyz[mi][ji]                          # (npairs, 3)

        def lddt(xyz_batch):                                 # (B, n, 3)
            d = np.linalg.norm(xyz_batch[:, ki, :] - Pm[None, :, :], axis=2)
            diff = np.abs(d - dref_p[None, :])
            return np.mean([(diff < t).mean(axis=1)
                            for t in (0.5, 1.0, 2.0, 4.0)], axis=0)

        L_ld[si] = lddt(cand[:, :len(perm), :])
        ld0 = float(lddt(L0[None, :len(perm), :])[0])

        # ---- BiSyRMSD: one site superposition per pose --------------------
        try:
            mca, rca = model.ca(), ref.ca()
            mnum = {r: v for (_c, r), v in mca.items()}
            rnum = {r: v for (_c, r), v in rca.items()}
            sh = sorted(set(mnum) & set(rnum) & set(ref_pocket))
            Rk, tk, _ = P.kabsch(np.array([mnum[r] for r in sh]),
                                 np.array([rnum[r] for r in sh]))
            lr_full = ref.lig_xyz[perm]
            al = cand[:, :len(perm), :] @ Rk.T + tk
            L_rm[si] = np.sqrt(((al - lr_full[None]) ** 2).sum(2).mean(1))
            rm0 = float(np.sqrt((((L0[:len(perm)] @ Rk.T + tk) - lr_full) ** 2)
                                .sum(1).mean()))
        except Exception:
            L_rm[si] = np.nan
            rm0 = float("nan")

        # ---- clash against the model's OWN protein + heme, Fe excluded ----
        heme = np.asarray(model.heme_xyz, float)
        he = [e.upper() for e in model.heme_elem]
        heme_nofe = heme[[i for i, e in enumerate(he) if e != "FE"]] if len(heme) else heme
        env = np.vstack([model.prot_xyz, heme_nofe]) if len(heme_nofe) else model.prot_xyz
        keep = np.linalg.norm(env - c, axis=1) < (np.linalg.norm(L0 - c, axis=1).max() + 8.0)
        tree = cKDTree(env[keep])
        dmin, _ = tree.query(cand.reshape(-1, 3), k=1)
        L_cl[si] = dmin.reshape(N_GRID, -1).min(axis=1)
        cl0 = float(tree.query(L0, k=1)[0].min())

        # ---- Fe coordination, vectorised over the grid --------------------
        els = [e.upper() for e in model.lig_elem]
        don = [i for i, e in enumerate(els) if e in ("N", "O", "S")]
        fe = model.fe
        if fe is not None and don:
            dd = np.linalg.norm(cand[:, don, :] - fe[None, None, :], axis=2)
            dj = dd.argmin(axis=1)
            L_fd[si] = dd[np.arange(N_GRID), dj]
            if model.axial_sg is not None:
                v1 = np.asarray(model.axial_sg, float) - fe
                v2 = cand[np.arange(N_GRID), np.array(don)[dj], :] - fe
                cth = (v2 @ v1) / (np.linalg.norm(v2, axis=1) * np.linalg.norm(v1) + 1e-12)
                ang = np.degrees(np.arccos(np.clip(cth, -1, 1)))
            else:
                ang = np.full(N_GRID, 180.0)
            L_co[si] = ((L_fd[si] >= G.COORD_LO) & (L_fd[si] <= G.COORD_HI)
                        & (ang >= G.TRANS_ANGLE_MIN))
        else:
            L_fd[si] = np.nan
        g0 = G.compute(L0, model.lig_elem, model.prot_xyz, model.heme_xyz,
                       model.heme_atom, model.axial_sg, heme_elem=model.heme_elem)

        # ---- xeng in the pose's OWN heme frame ----------------------------
        fr = X.heme_frame(model.heme_xyz, model.heme_elem, model.axial_sg)
        if fr is not None and xrefs:
            fez, nz, xz = fr
            yz = np.cross(nz, xz)
            A = np.column_stack([xz, yz, nz])                # lab -> heme frame
            V = (cand - fez[None, None, :]) @ A              # (G, n, 3)
            V0 = (L0 - fez) @ A
            acc = np.zeros(N_GRID)
            for r in xrefs:
                d = np.linalg.norm(V[:, :, None, :] - r[None, None, :, :], axis=3)
                acc += 0.5 * (d.min(axis=2).mean(axis=1) + d.min(axis=1).mean(axis=1))
            L_xe[si] = acc / len(xrefs)
            xe0 = float(np.mean([X.chamfer(V0, r) for r in xrefs]))
        else:
            L_xe[si] = np.nan
            xe0 = float("nan")

        samples.append(cif.stem)
        meta.append((ld0, rm0, xe0, cl0, float(g0.fe_donor_dist),
                     bool(g0.is_coordinated), int(off), float(ident), bool(mapped),
                     int(len(model.lig_xyz)), int(n_lig_ref)))

    np.savez_compressed(
        out, samples=np.array(samples), lddt=L_ld, bisy=L_rm, xeng=L_xe,
        clash=L_cl, fedon=L_fd, coord=L_co,
        meta=np.array(meta, dtype=object), pdb=np.array(pdb))
    return {"ligand": lig_id, "status": "ok", "poses": len(cifs)}


# --------------------------------------------------------------------------
# stage 1
# --------------------------------------------------------------------------

def cmd_expand(workers: int, limit: int | None) -> None:
    import pandas as pd
    SCRATCH.mkdir(parents=True, exist_ok=True)
    lig = pd.read_csv(REPO / "data" / "processed" / "validation_ligands.csv")
    have = {p.name.split("__")[0] for p in POOL.iterdir() if p.is_dir()}
    jobs = [(r.id, r.pdb, r.smiles) for r in lig.itertuples() if r.id in have]
    if limit:
        jobs = jobs[:limit]
    todo = [j for j in jobs if not (SCRATCH / f"{j[0]}.npz").exists()]
    print(f"{len(jobs)} ligands, {len(todo)} to do, {workers} workers", flush=True)
    print("grid:", grid_spacing(super_fibonacci(N_GRID)), flush=True)
    if not todo:
        return
    if workers <= 1:
        for j in todo:
            print(run_ligand(j), flush=True)
        return
    import multiprocessing as mp
    with mp.Pool(workers) as pool:
        for res in pool.imap_unordered(run_ligand, todo):
            print(res, flush=True)


# --------------------------------------------------------------------------
# stage 2
# --------------------------------------------------------------------------

def _select_xeng(xe: np.ndarray, truth: np.ndarray) -> float:
    """The shipped selector on one ligand: argmin of xeng (z-scoring within a ligand is
    monotone, so argmin of the raw distance is the identical choice)."""
    ok = np.isfinite(xe) & np.isfinite(truth)
    if not ok.any():
        return float("nan")
    return float(truth[ok][int(np.argmin(xe[ok]))])


def cmd_analyse() -> None:
    import pandas as pd
    from scipy import stats

    files = sorted(SCRATCH.glob("*.npz"))
    rng = np.random.default_rng(NULL_SEED)
    rows, ctrl = [], []

    for f in files:
        d = np.load(f, allow_pickle=True)
        lig = f.stem
        meta = d["meta"]
        ld, rm, xe = d["lddt"], d["bisy"], d["xeng"]
        cl, co = d["clash"], d["coord"]
        n_pose = ld.shape[0]

        o_ld = np.array([m[0] for m in meta], float)
        o_rm = np.array([m[1] for m in meta], float)
        o_xe = np.array([m[2] for m in meta], float)
        o_co = np.array([m[5] for m in meta], bool)
        for i in range(n_pose):
            ctrl.append(dict(ligand=lig, sample=str(d["samples"][i]), lddt=o_ld[i],
                             bisy=o_rm[i], xeng=o_xe[i], clash=meta[i][3],
                             offset=meta[i][6], identity=meta[i][7],
                             mapped=meta[i][8], n_lig=meta[i][9], n_ref=meta[i][10]))

        # pre-registered filter: clash >= 2.2 A, and keep coordination if the pose had it
        passf = cl >= CLASH_CUT
        need = o_co[:, None] & np.ones_like(passf)
        passf = passf & (~need | co)
        n_pass = int(passf.sum())

        # size-matched unfiltered null: same count per pose, no filter at all
        passn = np.zeros_like(passf)
        for i in range(n_pose):
            k = int(passf[i].sum())
            if k:
                passn[i, rng.choice(ld.shape[1], size=k, replace=False)] = True

        # sensitivity variants (reported, never used to decide)
        pass_soft = (cl >= CLASH_CUT_SOFT) & (~need | co)
        fd = d["fedon"]
        co_narrow = (fd >= 1.90) & (fd <= 2.45)
        pass_narrow = (cl >= CLASH_CUT) & (~need | co_narrow)

        def pack(m):
            return (np.concatenate([o_ld, ld[m]]), np.concatenate([o_rm, rm[m]]),
                    np.concatenate([o_xe, xe[m]]))

        f_ld, f_rm, f_xe = pack(passf)
        n_ld, n_rm, n_xe = pack(passn)
        s_ld, _s_rm, s_xe = pack(pass_soft)
        w_ld, _w_rm, w_xe = pack(pass_narrow)

        surv_ld, surv_xe = ld[passf], xe[passf]
        best_o_xe = np.nanmin(o_xe)
        rows.append(dict(
            surv_mean_lddt=float(np.nanmean(surv_ld)) if n_pass else np.nan,
            frac_surv_better_xeng=(float(np.mean(surv_xe < best_o_xe))
                                   if n_pass else np.nan),
            picks_rotated=bool(n_pass and np.nanmin(surv_xe) < best_o_xe),
            pass_per_pose=n_pass / n_pose,
            frac_pose_improved=float(np.mean(
                [(ld[i][passf[i]].max() > o_ld[i]) if passf[i].any() else False
                 for i in range(n_pose)])),
            grid_min_rmsd=float(np.nanmin(np.concatenate([o_rm, rm.ravel()]))),
            ligand=lig, n_pose=n_pose, n_grid=int(ld.shape[1]),
            n_cand=int(ld.size), n_pass=n_pass,
            pass_rate=n_pass / ld.size,
            pass_clash=float((cl >= CLASH_CUT).mean()),
            frac_pose_coord=float(o_co.mean()),
            orig_oracle=float(np.nanmax(o_ld)), orig_mean=float(np.nanmean(o_ld)),
            orig_sel=_select_xeng(o_xe, o_ld),
            exp_oracle=float(np.nanmax(f_ld)), exp_sel=_select_xeng(f_xe, f_ld),
            null_oracle=float(np.nanmax(n_ld)), null_sel=_select_xeng(n_xe, n_ld),
            soft_oracle=float(np.nanmax(s_ld)), soft_sel=_select_xeng(s_xe, s_ld),
            narrow_oracle=float(np.nanmax(w_ld)), narrow_sel=_select_xeng(w_xe, w_ld),
            unfilt_oracle=float(np.nanmax(np.concatenate([o_ld, ld.ravel()]))),
            orig_min_rmsd=float(np.nanmin(o_rm)), exp_min_rmsd=float(np.nanmin(f_rm)),
            null_min_rmsd=float(np.nanmin(n_rm)),
            orig_sub2=bool(np.nanmin(o_rm) < 2.0), exp_sub2=bool(np.nanmin(f_rm) < 2.0),
            null_sub2=bool(np.nanmin(n_rm) < 2.0),
            # within-ligand diagnostics
            wl_sd_lddt_orig=float(np.nanstd(o_ld)), wl_sd_lddt_exp=float(np.nanstd(f_ld)),
            wl_rho_exp=float(stats.spearmanr(f_xe, f_ld, nan_policy="omit").statistic),
            wl_rho_orig=float(stats.spearmanr(o_xe, o_ld, nan_policy="omit").statistic),
            exp_sel_rank=float(np.nanmean(f_ld <= _select_xeng(f_xe, f_ld))),
        ))

    df = pd.DataFrame(rows)
    cdf = pd.DataFrame(ctrl)
    df.to_csv(REPO / "data" / "processed" / "orientation_expansion_per_ligand.csv",
              index=False)

    # ---- controls ---------------------------------------------------------
    ship = pd.read_csv(REPO / "data" / "processed" / "poses_scored_val87b.csv")
    ship = ship[ship.arm == "unsteered"][["ligand", "sample", "lddt_pli", "bisy_rmsd"]]
    m = cdf.merge(ship, on=["ligand", "sample"], how="inner")
    xs = pd.read_csv(REPO / "data" / "processed" / "xeng_val87b.csv")
    m = m.merge(xs, on=["ligand", "sample"], how="left")
    c1 = np.abs(m.lddt - m.lddt_pli)
    c3 = np.abs(m.xeng_x - m.xeng_y) if "xeng_x" in m else np.abs(m.xeng - m.xeng)
    zeros = m[(m.lddt_pli == 0)]
    controls = {
        "C1_rows_compared": int(len(m)),
        "C1_max_abs_delta_lddt": float(c1.max()),
        "C1_frac_within_1e-9": float((c1 < 1e-9).mean()),
        "C2_exact_zero_rows": int((m.lddt == 0).sum()),
        "C2_exact_zeros_all_ejected": bool((zeros.bisy_rmsd > 10).all()) if len(zeros) else True,
        "C2_offsets": {str(k): int(v) for k, v in
                       cdf.groupby("offset").size().items()},
        "C2_min_identity": float(cdf.identity.min()),
        "C3_max_abs_delta_xeng": float(np.nanmax(c3)),
        "C4_atom_counts_agree": int((cdf.n_lig == cdf.n_ref).sum()),
        "C4_rows": int(len(cdf)),
        "C4_mapped": int(cdf.mapped.sum()),
    }

    # ---- headline ---------------------------------------------------------
    def paired(a, b, label):
        x = (df[a] - df[b]).values
        x = x[np.isfinite(x)]
        boot = np.array([np.mean(rng.choice(x, len(x), replace=True))
                         for _ in range(20000)])
        w = stats.wilcoxon(x) if np.any(x != 0) else None
        return {"label": label, "n": int(len(x)), "mean_diff": float(x.mean()),
                "ci95": [float(np.percentile(boot, 2.5)),
                         float(np.percentile(boot, 97.5))],
                "wilcoxon_p": float(w.pvalue) if w else 1.0,
                "better": int((x > 0).sum()), "worse": int((x < 0).sum()),
                "tied": int((x == 0).sum())}

    tests = [paired("exp_sel", "orig_sel", "filtered expansion vs original"),
             paired("null_sel", "orig_sel", "unfiltered size-matched null vs original"),
             paired("exp_sel", "null_sel", "filtered vs unfiltered null")]
    ps = np.array([t["wilcoxon_p"] for t in tests])
    order = np.argsort(ps)
    holm = np.empty(3)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (3 - rank) * ps[i])
        holm[i] = min(1.0, run)
    for t, h in zip(tests, holm):
        t["holm_p"] = float(h)

    rnd = float(np.mean([np.nanmean(np.load(f, allow_pickle=True)["meta"][:, 0]
                                    .astype(float)) for f in files]))

    summary = {
        "grid": grid_spacing(super_fibonacci(N_GRID)),
        "n_ligands": int(len(df)), "n_poses": int(df.n_pose.sum()),
        "n_candidates": int(df.n_cand.sum()), "n_surviving": int(df.n_pass.sum()),
        "pass_rate": float(df.n_pass.sum() / df.n_cand.sum()),
        "pass_rate_clash_only": float(df.pass_clash.mean()),
        "frac_poses_coordinating": float(df.frac_pose_coord.mean()),
        "controls": controls,
        "random_pose_mean": rnd,
        "oracle": {"original": float(df.orig_oracle.mean()),
                   "filtered_expanded": float(df.exp_oracle.mean()),
                   "unfiltered_size_matched": float(df.null_oracle.mean()),
                   "whole_grid_no_filter": float(df.unfilt_oracle.mean())},
        "selection": {"original": float(df.orig_sel.mean()),
                      "filtered_expanded": float(df.exp_sel.mean()),
                      "unfiltered_size_matched": float(df.null_sel.mean())},
        "sensitivity": {
            "clash_3.0A": {"oracle": float(df.soft_oracle.mean()),
                           "selection": float(df.soft_sel.mean())},
            "coord_1.90_2.45": {"oracle": float(df.narrow_oracle.mean()),
                                "selection": float(df.narrow_sel.mean())}},
        "sub2": {"original_reaching": int(df.orig_sub2.sum()),
                 "expanded_reaching": int(df.exp_sub2.sum()),
                 "null_reaching": int(df.null_sub2.sum()),
                 "rescued": int((df.exp_sub2 & ~df.orig_sub2).sum()),
                 "n": int(len(df))},
        "within_ligand": {
            "rho_xeng_lddt_original": float(df.wl_rho_orig.mean()),
            "rho_xeng_lddt_expanded": float(df.wl_rho_exp.mean()),
            "frac_rho_negative_original": float((df.wl_rho_orig < 0).mean()),
            "frac_rho_negative_expanded": float((df.wl_rho_exp < 0).mean()),
            "sd_lddt_original": float(df.wl_sd_lddt_orig.mean()),
            "sd_lddt_expanded": float(df.wl_sd_lddt_exp.mean()),
            "selected_percentile_in_expanded_pool": float(df.exp_sel_rank.mean())},
        "paired_tests": tests,
        "mechanism": {
            "survivors_per_pose_mean": float(df.pass_per_pose.mean()),
            "survivors_per_pose_median": float(df.pass_per_pose.median()),
            "ligands_with_zero_survivors": int((df.n_pass == 0).sum()),
            "surviving_candidate_mean_lddt": float(df.surv_mean_lddt.mean()),
            "original_pose_mean_lddt": float(df.orig_mean.mean()),
            "frac_survivors_with_better_xeng_than_best_original":
                float(df.frac_surv_better_xeng.mean()),
            "ligands_where_selector_picks_a_rotation": int(df.picks_rotated.sum()),
            "frac_poses_whose_best_rotation_beats_the_pose":
                float(df.frac_pose_improved.mean()),
            "oracle_gain_on_failed_ligands": float(
                (df.exp_oracle - df.orig_oracle)[df.orig_min_rmsd >= 2].mean()),
            "oracle_gain_on_solved_ligands": float(
                (df.exp_oracle - df.orig_oracle)[df.orig_min_rmsd < 2].mean()),
            "ligands_with_any_oracle_gain": int(
                ((df.exp_oracle - df.orig_oracle) > 1e-9).sum()),
            "failed_ligands": int((df.orig_min_rmsd >= 2).sum()),
            "failed_reaching_sub2_filtered": int(
                (df.exp_min_rmsd[df.orig_min_rmsd >= 2] < 2).sum()),
            "failed_reaching_sub2_whole_grid": int(
                (df.grid_min_rmsd[df.orig_min_rmsd >= 2] < 2).sum()),
            "failed_median_min_rmsd_original": float(
                df.orig_min_rmsd[df.orig_min_rmsd >= 2].median()),
            "failed_median_min_rmsd_filtered": float(
                df.exp_min_rmsd[df.orig_min_rmsd >= 2].median()),
            "failed_median_min_rmsd_whole_grid": float(
                df.grid_min_rmsd[df.orig_min_rmsd >= 2].median()),
        },
    }
    out = REPO / "data" / "processed" / "orientation_expansion.json"
    out.write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


# --------------------------------------------------------------------------
# stage 3 — the mechanism, including the measurement that decided the finding
# --------------------------------------------------------------------------

def cmd_diagnose() -> None:
    """Why the lever fails — and the control that validates the clash cutoff.

    The decisive measurement is the last one: take the CRYSTAL ligand, put it inside the
    MODEL's own protein, and ask whether it clashes. It uses the crystal, so it is a
    diagnostic and never a selection feature. If the true orientation does not fit in the
    pocket the model built, no search over orientations inside that pocket can find it.
    """
    import gemmi
    import pandas as pd
    from scipy.spatial import cKDTree

    from cypstruct import pose as P
    from cypstruct.targets import fetch_cif

    lig = pd.read_csv(REPO / "data" / "processed" / "validation_ligands.csv")
    rows, native = [], []
    for r in lig.itertuples():
        job = POOL / f"{r.id}__unsteered__s1"
        if not job.is_dir():
            continue
        ref = load_reference(r.pdb, r.id)
        if ref is None or len(ref.lig_xyz) == 0:
            continue

        # control: the crystal ligand against its OWN crystal protein. This is the
        # non-bonded floor the 2.2 A cutoff was pre-registered against.
        hm = np.asarray(ref.heme_xyz, float)
        he = [e.upper() for e in ref.heme_elem]
        hn = hm[[i for i, e in enumerate(he) if e != "FE"]] if len(hm) else hm
        envr = np.vstack([ref.prot_xyz, hn]) if len(hn) else ref.prot_xyz
        native.append(float(cKDTree(envr).query(ref.lig_xyz, k=1)[0].min()))

        pock = P.pocket_residues_from_structure(ref, radius=8.0)
        for cif in sorted(job.glob("input_model_*.cif"),
                          key=lambda p: int(p.stem.rsplit("_", 1)[1])):
            m = P.load_structure(cif)
            m, _o, _i = renumber_to_reference(m, ref)
            perm = P.best_ligand_mapping(r.smiles, m, ref)
            if perm is None:
                continue
            L0 = np.asarray(m.lig_xyz, float)
            c = L0.mean(0)
            mn = {k: v for (_c, k), v in m.ca().items()}
            rn = {k: v for (_c, k), v in ref.ca().items()}
            sh = sorted(set(mn) & set(rn) & set(pock))
            Rk, tk, _ = P.kabsch(np.array([rn[k] for k in sh]),
                                 np.array([mn[k] for k in sh]))    # crystal -> model
            Ltrue = ref.lig_xyz[perm] @ Rk.T + tk
            hm = np.asarray(m.heme_xyz, float)
            he = [e.upper() for e in m.heme_elem]
            hn = hm[[i for i, e in enumerate(he) if e != "FE"]] if len(hm) else hm
            env = np.vstack([m.prot_xyz, hn]) if len(hn) else m.prot_xyz
            keep = np.linalg.norm(env - c, axis=1) < (
                np.linalg.norm(L0 - c, axis=1).max() + 8.0)
            t = cKDTree(env[keep])
            rows.append(dict(
                ligand=r.id, sample=cif.stem,
                rg=float(np.sqrt(((L0 - c) ** 2).sum(1).mean())),
                clash_pred=float(t.query(L0, k=1)[0].min()),
                clash_true_in_place=float(t.query(Ltrue, k=1)[0].min()),
                clash_true_at_pred_centroid=float(
                    t.query(Ltrue - Ltrue.mean(0) + c, k=1)[0].min()),
                centroid_off=float(np.linalg.norm(Ltrue.mean(0) - c))))

    d = pd.DataFrame(rows)
    d.to_csv(SCRATCH / "truth_clash.csv", index=False)
    rg = float(d.rg.median())
    res = {
        "n_poses": int(len(d)),
        "ligand_radius_of_gyration_median": rg,
        "grid_half_spacing_deg": grid_spacing(super_fibonacci(N_GRID))["nn_mean_deg"] / 2,
        "rmsd_implied_by_grid_half_spacing": float(
            2 * rg * np.sin(np.radians(grid_spacing(super_fibonacci(N_GRID))["nn_mean_deg"] / 2) / 2)),
        "crystal_ligand_vs_own_crystal_protein": {
            "n": len(native), "min": float(np.min(native)),
            "p1": float(np.percentile(native, 1)),
            "median": float(np.median(native)),
            "frac_below_2.2A": float(np.mean(np.array(native) < CLASH_CUT))},
        "predicted_pose_in_own_protein": {
            "median": float(d.clash_pred.median()),
            "p5": float(d.clash_pred.quantile(0.05)),
            "frac_below_2.2A": float((d.clash_pred < CLASH_CUT).mean())},
        "crystal_ligand_in_MODEL_protein_in_place": {
            "median": float(d.clash_true_in_place.median()),
            "frac_below_2.2A": float((d.clash_true_in_place < CLASH_CUT).mean())},
        "crystal_ligand_in_MODEL_protein_at_predicted_centroid": {
            "median": float(d.clash_true_at_pred_centroid.median()),
            "frac_below_2.2A": float((d.clash_true_at_pred_centroid < CLASH_CUT).mean())},
        "centroid_offset_median": float(d.centroid_off.median()),
    }
    (REPO / "data" / "processed" / "orientation_expansion_diagnostics.json").write_text(
        json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["expand", "analyse", "diagnose"])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    if a.stage == "expand":
        cmd_expand(a.workers, a.limit)
    elif a.stage == "diagnose":
        cmd_diagnose()
    else:
        cmd_analyse()


if __name__ == "__main__":
    main()
