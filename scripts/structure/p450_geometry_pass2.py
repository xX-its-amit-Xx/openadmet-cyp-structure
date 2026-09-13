"""Second pass over the P450 universe: the full active-site frame, cached so it is final.

Pass 1 (`harvest_p450_universe.py`) answered "which entries have a ligand in the pocket"
and threw the coordinates away. That was the right call for a survey and the wrong one
for everything after it: every geometric question since - the S-Fe-donor angle, the
out-of-plane displacement, the proximal-face fraction - has needed another 585 downloads.

So this pass extracts the whole active-site frame and **caches the atoms**: iron, the four
pyrrole nitrogens, the proximal cysteine sulfur, and every ligand atom with its element.
That is about fifty atoms per pair, a few MB for the entire superfamily, and it makes any
future geometric term a local computation instead of a network campaign.

The terms measured here are the ones `cypstruct.qmscore.geometry` is calibrated on. Those
windows currently come from 116 CYP3A4 entries. Recomputing them over ~500 pairs spanning
the whole P450 fold tests something the CYP3A4-only calibration cannot: whether they are
P450 chemistry or CYP3A4 idiosyncrasy.

    python scripts/structure/p450_geometry_pass2.py
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402
from cypstruct.targets import HEME_ALIASES  # noqa: E402

OUT = DATA_PROCESSED / "p450_universe"
COORD_MAX = 2.6


def frame(chain) -> dict | None:
    """Iron, porphyrin normal and proximal thiolate for one chain.

    The normal is taken from the four pyrrole nitrogens by SVD rather than from a cross
    product of two of them: with a ruffled or domed porphyrin any single pair gives a
    normal tilted by several degrees, and the sign convention then flips between entries.
    """
    fe = None
    ns: list[np.ndarray] = []
    for res in chain:
        if res.name.strip().upper() not in HEME_ALIASES:
            continue
        for at in res:
            e = at.element.name.upper()
            p = np.array([at.pos.x, at.pos.y, at.pos.z])
            if e == "FE":
                fe = p
            elif e == "N":
                ns.append(p)
        if fe is not None and len(ns) >= 4:
            break
    if fe is None or len(ns) < 4:
        return None

    # the four nitrogens closest to the iron are the pyrroles; extras are propionate
    # or solvent nitrogens that would tilt the plane if included
    ns = sorted(ns, key=lambda p: np.linalg.norm(p - fe))[:4]
    P = np.array(ns) - np.array(ns).mean(axis=0)
    normal = np.linalg.svd(P)[2][-1]
    normal /= np.linalg.norm(normal)

    sg, sg_d = None, None
    for res in chain:
        if res.name.strip().upper() != "CYS":
            continue
        for at in res:
            if at.name.strip().upper() == "SG":
                p = np.array([at.pos.x, at.pos.y, at.pos.z])
                d = float(np.linalg.norm(p - fe))
                if sg_d is None or d < sg_d:
                    sg, sg_d = p, d
    if sg is None:
        return None
    # Orient the normal to point DISTAL, i.e. away from the thiolate. Without this the
    # sign of every out-of-plane number is arbitrary and "proximal face" is meaningless.
    if np.dot(normal, sg - fe) > 0:
        normal = -normal
    return {"fe": fe, "normal": normal, "sg": sg, "sg_d": sg_d}


def measure(pid: str, cif: Path, want: set[tuple[str, str]]) -> tuple[list, list]:
    import gemmi

    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    st.remove_alternative_conformations()
    rows, atoms = [], []

    for chain in st[0]:
        fr = frame(chain)
        if fr is None:
            continue
        fe, normal = fr["fe"], fr["normal"]
        for res in chain:
            name = res.name.strip().upper()
            if (pid, name) not in want:
                continue
            xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res])
            els = [a.element.name.upper() for a in res]
            heavy = [i for i, e in enumerate(els) if e != "H"]
            if len(heavy) < 12:
                continue
            d_fe = np.linalg.norm(xyz - fe, axis=1)
            # signed height above the porphyrin plane; negative = proximal (buried) face
            h = (xyz - fe) @ normal

            cands = [i for i in heavy if els[i] in ("N", "O", "S")]
            don = min(cands, key=lambda i: d_fe[i]) if cands else None
            ang = None
            if don is not None:
                v1, v2 = fr["sg"] - fe, xyz[don] - fe
                ang = float(np.degrees(np.arccos(np.clip(
                    v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1))))
            rows.append({
                "pdb": pid, "chain": chain.name, "lig": name,
                "n_heavy": len(heavy),
                "closest_fe": float(d_fe[heavy].min()),
                "donor_elem": els[don] if don is not None else None,
                "donor_dist": float(d_fe[don]) if don is not None else None,
                "donor_height": float(h[don]) if don is not None else None,
                "s_fe_donor_deg": ang,
                "coordinated": bool(don is not None and d_fe[don] <= COORD_MAX),
                "cys_sg_fe": fr["sg_d"],
                "frac_proximal": float(np.mean(h[heavy] < 0)),
                "centroid_height": float(h[heavy].mean()),
                "max_height": float(h[heavy].max()),
            })
            atoms.append({
                "pdb": pid, "chain": chain.name, "lig": name,
                "elems": "".join(f"{els[i]}," for i in heavy),
                "xyz": np.asarray(xyz[heavy], dtype=np.float32).tobytes(),
                "fe": fe.astype(np.float32).tobytes(),
                "normal": normal.astype(np.float32).tobytes(),
                "sg": fr["sg"].astype(np.float32).tobytes(),
            })
    return rows, atoms


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default=str(OUT / "druglike_pairs.parquet"))
    a = ap.parse_args()

    pairs = pd.read_parquet(a.pairs)
    want: dict[str, set] = {}
    for r in pairs.itertuples():
        want.setdefault(r.pdb, set()).add((r.pdb, r.lig))
    print(f"{len(pairs)} pairs across {len(want)} entries", flush=True)

    done_path = OUT / "geometry_pass2.jsonl"
    seen = set()
    if done_path.exists():
        for ln in done_path.read_text().splitlines():
            if ln.strip():
                seen.add(json.loads(ln)["pdb"])

    tmp = Path(tempfile.gettempdir()) / "p450_cif"
    tmp.mkdir(exist_ok=True)
    atom_store: list[dict] = []
    todo = [p for p in want if p not in seen]
    print(f"to fetch: {len(todo)}", flush=True)

    with done_path.open("a") as fh:
        for i, pid in enumerate(todo):
            f = tmp / f"{pid}.cif"
            try:
                if not f.exists():
                    r = requests.get(
                        f"https://files.rcsb.org/download/{pid}.cif", timeout=120)
                    r.raise_for_status()
                    f.write_bytes(r.content)
                rows, atoms = measure(pid, f, want[pid])
                fh.write(json.dumps({"pdb": pid, "rows": rows}) + "\n")
                fh.flush()
                atom_store.extend(atoms)
            except Exception as exc:
                fh.write(json.dumps({"pdb": pid, "rows": [],
                                     "error": f"{type(exc).__name__}: {exc}"}) + "\n")
                fh.flush()
            finally:
                f.unlink(missing_ok=True)
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(todo)}", flush=True)

    recs = []
    for ln in done_path.read_text().splitlines():
        if ln.strip():
            recs.extend(json.loads(ln).get("rows", []))
    df = pd.DataFrame(recs)
    df.to_parquet(OUT / "p450_geometry.parquet")
    if atom_store:
        ap_path = OUT / "p450_atoms.parquet"
        new = pd.DataFrame(atom_store)
        if ap_path.exists():
            new = pd.concat([pd.read_parquet(ap_path), new]) \
                    .drop_duplicates(["pdb", "chain", "lig"])
        new.to_parquet(ap_path)
        print(f"cached atoms for {len(new)} pairs -> {ap_path.name}")

    print(f"\n{len(df)} measured pairs, {df.lig.nunique()} unique ligands")
    c = df[df.coordinated]
    for col in ("donor_dist", "s_fe_donor_deg", "donor_height"):
        v = c[col].dropna()
        if len(v):
            p5, p50, p95 = np.percentile(v, [5, 50, 95])
            print(f"  {col:16s} {p5:7.2f} / {p50:7.2f} / {p95:7.2f}   (n={len(v)})")
    print(f"  frac_proximal    mean {df.frac_proximal.mean():.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
