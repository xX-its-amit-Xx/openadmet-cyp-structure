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


def renumber_to_reference(pred, ref):
    """Shift the prediction's residue numbering onto the crystal's.

    **Not cosmetic: without it most targets score exactly 0.000.** We fold the DEPOSITED
    CONSTRUCT sequence, which the engine numbers from 1, while the crystal keeps author
    numbering that usually starts elsewhere - 3TK3 runs 28..492 against a prediction
    numbered 1..476, an offset of 27. `align_by_residue` matches by residue NUMBER, so
    every contact then compares two different residues and LDDT-PLI collapses to zero.
    It looks exactly like a target the model cannot fold, which is why it survived a first
    pass: 20 targets, several at precisely 0.000, read as "co-folding is hard here".

    The offset maximises three-letter residue-name agreement over candidate shifts - no
    alignment library, and it cannot silently choose a frame where the sequences disagree.
    Returns (pred, offset, identity); identity below ~0.8 means the pairing is wrong and
    the pair should be dropped rather than scored.
    """
    pm = {num: name for (_ch, num), name in pred.prot_res.items()}
    rm = {num: name for (_ch, num), name in ref.prot_res.items()}
    if len(pm) < 30 or len(rm) < 30:
        return pred, 0, 0.0

    # Search EVERY shift that leaves the two ranges overlapping. The obvious narrow window
    # - anchoring the prediction's first residue to the crystal's first OBSERVED one - is
    # wrong whenever the crystal has a disordered N-terminus, which is most of the time:
    # 2JJO's first modelled residue is 19 and its true offset is 0, so a +-10 window around
    # 18 never tested the right answer and the pair was discarded as unmatchable. That
    # dropped 189 of 351 poses while looking like a data problem.
    best, best_id = 0, -1.0
    lo = min(rm) - max(pm)
    hi = max(rm) - min(pm)
    for off in range(int(lo), int(hi) + 1):
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

    # single-sequence poses are a separate, far worse population (oracle 0.0505 against
    # 0.8760, 75% catastrophic) and must not be averaged into the main pool - including
    # them dragged the reported P450 oracle from 0.7368 to 0.6254.
    skip_engines = {"protenix_v2_ss"}
    rows, n_err = [], 0
    dirs = sorted(p for p in POSES.glob("*") if p.is_dir())
    if limit:
        dirs = dirs[:limit]
    print(f"{len(dirs)} pair directories", flush=True)

    for i, d in enumerate(dirs):
        m = meta.get(d.name)
        if m is None:
            continue
        # poses live in per-engine subdirectories; recurse so a second engine's pool is
        # picked up rather than silently ignored
        files = [f for f in sorted(d.rglob("*.cif"))
                 if f.parent.name not in skip_engines
                 and f"{f.parent.name}/{f.name}" not in done]
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
                mo, off, ident = renumber_to_reference(mo, ref)
                if ident < 0.8:
                    n_err += 1
                    continue
                perm = P.best_ligand_mapping(m.smiles, mo, ref)
                rows.append({
                    "pair": d.name, "pose": f"{f.parent.name}/{f.name}",
                    "engine": f.parent.name, "ligand": m.id, "pdb": m.pdb,
                    "target_key": m.target_key, "uniprot": m.uniprot,
                    "coordinated": m.coordinated, "n_heavy": m.n_heavy,
                    "resnum_offset": off, "seq_identity": round(float(ident), 3),
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
                loaded[f] = P.load_structure(d / f)   # pose is "<engine>/<file>.cif"
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


def cmd_validate(n_null: int = 3000) -> dict:
    """The FINDING 012 generalisation test, as one command.

    Pool = protenix_v2 poses. Reference = esmfold2 poses of the SAME pair, one per
    replicate and deduplicated. Cross-engine by construction, on proteins the feature was
    never developed on, with **no parameters fitted** - so every protein is held out and
    there is no fold to leak across.

    This replaces an earlier `validate` that fitted three selector weights
    leave-one-target-out. That version measured +0.0149, inside the noise, because fitting
    weights on this much data overfits - the same result FINDING 002 got from a fitted
    ranker. The unweighted single term is what ships and what is tested here.
    """
    import re

    from scipy import stats

    from cypstruct import pose as P
    from cypstruct import xengine as X

    sc = pd.read_csv(SCORED)

    def refs(pair: str) -> list:
        out: dict[int, object] = {}
        d = POSES / pair / "esmfold2"
        for f in sorted(d.glob("*.cif")) if d.exists() else []:
            m = re.match(r"(.+)__r(\d+)s(\d+)\.cif$", f.name)
            rep = int(m.group(2)) if m else 0
            if rep in out:
                continue
            try:
                v = X.in_heme_frame(P.load_structure(f))
            except Exception:
                continue
            if v is not None:
                out[rep] = v
        return X._dedupe(list(out.values()))

    rows = []
    for pair in sc.pair.unique():
        r = refs(pair)
        if len(r) < 2:
            continue
        d = POSES / pair / "protenix_v2"
        for f in sorted(d.glob("*.cif")) if d.exists() else []:
            try:
                v = X.in_heme_frame(P.load_structure(f))
            except Exception:
                continue
            if v is None:
                continue
            rows.append({"pair": pair, "pose": f"protenix_v2/{f.name}",
                         "xeng": X.xeng_score(v, r)})
    if not rows:
        return {"error": "no pairs with >= 2 independent esmfold2 references yet"}

    m = sc.merge(pd.DataFrame(rows), on=["pair", "pose"])
    m = m[m.groupby("pair").pose.transform("size") >= 3]
    rng = np.random.default_rng(0)
    rand = float(np.mean([
        m.groupby("pair").lddt_pli.apply(lambda s: s.sample(1, random_state=i).iloc[0])
        .mean() for i in range(300)]))
    oracle = float(m.groupby("pair").lddt_pli.max().mean())
    null = np.array([
        float(m.assign(s=rng.normal(size=len(m)))
              .pipe(lambda t: t.loc[t.groupby("pair").s.idxmax()]).lddt_pli.mean()) - rand
        for _ in range(n_null)])
    sel = float(m.loc[m.groupby("pair").xeng.idxmin()].lddt_pli.mean())
    gain = sel - rand
    rs = np.array([stats.spearmanr(g.xeng, g.lddt_pli).statistic
                   for _p, g in m.groupby("pair") if g.xeng.nunique() >= 3])
    rs = rs[np.isfinite(rs)]

    per = []
    for t, g in m.groupby("uniprot"):
        if g.pair.nunique() < 4:
            continue
        r0 = float(np.mean([
            g.groupby("pair").lddt_pli.apply(lambda s: s.sample(1, random_state=i).iloc[0])
            .mean() for i in range(100)]))
        s0 = float(g.loc[g.groupby("pair").xeng.idxmin()].lddt_pli.mean())
        per.append((t, int(g.pair.nunique()), round(s0 - r0, 4)))

    return {
        "poses": len(m), "pairs": int(m.pair.nunique()),
        "proteins": int(m.uniprot.nunique()),
        "construct_sequences": int(m.target_key.nunique()),
        "catastrophic_frac": round(float((m.lddt_pli < 0.1).mean()), 4),
        "random": round(rand, 4), "oracle": round(oracle, 4),
        "selected": round(sel, 4), "gain": round(gain, 4),
        "null_p99": round(float(np.percentile(null, 99)), 4),
        "empirical_p": round(float((null >= gain).mean()), 4),
        "within_pair_rho": round(float(rs.mean()), 4),
        "proteins_positive": f"{sum(1 for x in per if x[2] > 0)}/{len(per)}",
        "per_protein": sorted(per, key=lambda x: -x[2]),
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
