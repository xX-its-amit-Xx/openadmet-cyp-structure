"""Score the P450 pool and put the selector through leave-one-TARGET-out.

This is the experiment the whole harvest exists for. Every selector result so far has been
leave-one-scaffold-cluster-out inside a single pocket (CYP3A4, 87 ligands), where FINDING
007 puts the noise floor at +0.0138 and the working selector at +0.0220 - clear, but by
less than a factor of two, on one protein.

With 185 distinct targets the held-out unit can be the **protein**. A term that still pays
on a P450 whose pocket the fit has never seen is evidence of physics; a term that
evaporates was fitting CYP3A4's F/G loop. That is a qualitatively stronger test than
anything available at n=87, and it is the only way to find out which of the two we have.

Three stages, each cached, so the expensive one runs once:

    python scripts/structure/score_p450_pool.py score     # poses -> LDDT-PLI (slow)
    python scripts/structure/score_p450_pool.py features  # pocket contacts + consensus
    python scripts/structure/score_p450_pool.py validate  # leave-one-target-out

**The pocket is derived per target, never from CYP3A4.** `CYP3A4_POCKET` is a residue list
in P08684 numbering and is meaningless for a bacterial P450. Contacts are counted against
whatever residues line *that* structure's heme, computed from the prediction alone so the
feature stays inference-time legal (repo rule: no leaky features).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

UNI = DATA_PROCESSED / "p450_universe"
POSES = UNI / "poses"
SCORED = UNI / "p450_poses_scored.csv"
FEATS = UNI / "p450_features.csv"


def _ref_for(pdb: str, lig: str, chain: str):
    """Crystal reference for one pair, restricted to the chain the ligand sits in."""
    import gemmi

    from cypstruct import pose as P
    from cypstruct.targets import fetch_cif

    cif = fetch_cif(pdb)
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    use = chain
    if not any(c.name == use and any(r.name.strip().upper() == lig for r in c)
               for c in st[0]):
        use = next((c.name for c in st[0]
                    if any(r.name.strip().upper() == lig for r in c)), None)
    if use is None:
        raise ValueError(f"{lig} not found in {pdb}")
    return P.load_structure(cif, ligand_code=lig, assembly_chain=use)


def cmd_score(limit: int | None) -> dict:
    from cypstruct import pose as P

    df = pd.read_csv(UNI / "p450_cofold_set.csv")
    df["pair"] = df.pdb + "_" + df.id
    meta = {r.pair: r for r in df.itertuples()}

    done = set()
    if SCORED.exists():
        prev = pd.read_csv(SCORED)
        done = set(prev["pose"])
    else:
        prev = None

    rows, n_err = [], 0
    dirs = sorted(p for p in POSES.glob("*") if p.is_dir())
    if limit:
        dirs = dirs[:limit]
    print(f"{len(dirs)} pair directories", flush=True)

    for i, d in enumerate(dirs):
        m = meta.get(d.name)
        if m is None:
            continue
        files = [f for f in sorted(d.glob("*.cif")) if f.name not in done]
        if not files:
            continue
        try:
            ref = _ref_for(m.pdb, m.id, m.chain)
        except Exception as exc:
            n_err += 1
            print(f"  {d.name}: no reference ({type(exc).__name__})", flush=True)
            continue
        for f in files:
            try:
                mo = P.load_structure(f)
                if len(mo.lig_xyz) == 0:
                    continue
                perm = P.best_ligand_mapping(m.smiles, mo, ref)
                rows.append({
                    "pair": d.name, "pose": f.name, "ligand": m.id, "pdb": m.pdb,
                    "target_key": m.target_key, "uniprot": m.uniprot,
                    "coordinated": m.coordinated, "n_heavy": m.n_heavy,
                    "lddt_pli": P.lddt_pli(mo, ref, lig_perm=perm),
                    "bisy_rmsd": P.bisy_rmsd(mo, ref, lig_perm=perm),
                })
            except Exception:
                n_err += 1
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{len(dirs)}  scored={len(rows)} err={n_err}", flush=True)

    new = pd.DataFrame(rows)
    out = pd.concat([prev, new]) if prev is not None and len(new) else (
        new if prev is None else prev)
    out.to_csv(SCORED, index=False)
    if len(out) == 0:
        return {"scored": 0, "note": "no poses collected yet"}
    return {"poses": len(out), "pairs": out.pair.nunique(),
            "targets": out.target_key.nunique(), "errors": n_err,
            "pool_mean": round(float(out.lddt_pli.mean()), 4),
            "oracle": round(float(out.groupby("pair").lddt_pli.max().mean()), 4)}


def cmd_features() -> dict:
    """Pocket contacts and consensus RMSD, both computed from the prediction alone."""
    from cypstruct import pose as P

    sc = pd.read_csv(SCORED)
    rows = []
    for pair, grp in sc.groupby("pair"):
        d = POSES / pair
        loaded = {}
        for f in grp.pose:
            try:
                loaded[f] = P.load_structure(d / f)
            except Exception:
                pass
        if len(loaded) < 2:
            continue
        names = list(loaded)
        L = {n: np.asarray(loaded[n].lig_xyz) for n in names}
        for n in names:
            mo = loaded[n]
            prot = np.asarray(mo.prot_xyz)
            lig = L[n]
            # pocket contacts: ligand heavy atoms within 4.5 A of any protein atom.
            # Derived from THIS structure, so no CYP3A4 residue list leaks in.
            dm = np.linalg.norm(lig[:, None, :] - prot[None, :, :], axis=2)
            n_contacts = int((dm < 4.5).sum())
            # consensus: mean RMSD to the other samples of the same pair
            others = [L[o] for o in names if o != n]
            k = min(len(lig), min(len(o) for o in others))
            rms = [float(np.sqrt(((lig[:k] - o[:k]) ** 2).sum(1).mean())) for o in others]
            rows.append({"pair": pair, "pose": n, "n_pocket_contacts": n_contacts,
                         "mean_rmsd_to_others": float(np.mean(rms))})
    f = pd.DataFrame(rows)
    out = sc.merge(f, on=["pair", "pose"], how="inner")
    out.to_csv(FEATS, index=False)
    return {"rows": len(out), "pairs": out.pair.nunique()}


def cmd_validate(n_boot: int = 2000) -> dict:
    """Leave-one-TARGET-out, against a random-feature null measured on this same pool."""
    from cypstruct.select import orientation_consensus_score

    d = pd.read_csv(FEATS)
    d = d[d.groupby("pair").pose.transform("size") >= 3]
    if d.empty:
        return {"error": "not enough poses per pair yet"}
    rng = np.random.default_rng(0)

    def mean_pick(score_col: pd.Series) -> float:
        idx = d.assign(s=score_col.values).groupby("pair").s.idxmax()
        return float(d.loc[idx].lddt_pli.mean())

    rand = float(np.mean([
        d.groupby("pair").lddt_pli.apply(lambda s: s.sample(1, random_state=i).iloc[0])
        .mean() for i in range(100)]))
    oracle = float(d.groupby("pair").lddt_pli.max().mean())

    # leave-one-target-out: the weight is fit on every OTHER protein
    picks = []
    for tgt in d.target_key.unique():
        tr, te = d[d.target_key != tgt], d[d.target_key == tgt]
        if tr.empty or te.empty:
            continue
        best_w, best_v = 0.5, -1.0
        for w in (0.0, 0.25, 0.5, 0.75, 1.0):
            s = orientation_consensus_score(tr.n_pocket_contacts,
                                            tr.mean_rmsd_to_others, tr.pair,
                                            contact_weight=w)
            idx = tr.assign(s=np.asarray(s)).groupby("pair").s.idxmax()
            v = float(tr.loc[idx].lddt_pli.mean())
            if v > best_v:
                best_w, best_v = w, v
        s = orientation_consensus_score(te.n_pocket_contacts, te.mean_rmsd_to_others,
                                        te.pair, contact_weight=best_w)
        idx = te.assign(s=np.asarray(s)).groupby("pair").s.idxmax()
        picks.extend(te.loc[idx].lddt_pli.tolist())
    selected = float(np.mean(picks)) if picks else float("nan")

    # the null that FINDING 007 insists on: what does a RANDOM feature score here?
    null = []
    for _ in range(min(n_boot, 500)):
        s = rng.normal(size=len(d))
        null.append(mean_pick(pd.Series(s)) - rand)
    null = np.array(null)

    return {
        "pairs": int(d.pair.nunique()), "targets": int(d.target_key.nunique()),
        "poses": int(len(d)),
        "random": round(rand, 4), "oracle": round(oracle, 4),
        "selected_LOTO": round(selected, 4),
        "gain_vs_random": round(selected - rand, 4),
        "null_sd": round(float(null.std()), 4),
        "null_p95": round(float(np.percentile(null, 95)), 4),
        "null_p99": round(float(np.percentile(null, 99)), 4),
        "empirical_p": round(float((null >= (selected - rand)).mean()), 4),
        "verdict": ("BEATS the null on held-out TARGETS"
                    if selected - rand > np.percentile(null, 99)
                    else "inside the noise floor - do not call it promising"),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["score", "features", "validate"])
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.cmd == "score":
        print(json.dumps(cmd_score(a.limit or None), indent=2))
    elif a.cmd == "features":
        print(json.dumps(cmd_features(), indent=2))
    else:
        print(json.dumps(cmd_validate(), indent=2))
