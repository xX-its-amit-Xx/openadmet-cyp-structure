"""Every deposited cytochrome P450 holo structure, not just the 87 CYP3A4 ligands.

**Why this and not more features.** Every negative result in this repo bottoms out in the
same place: n=87. FINDING 007 measured the noise floor at +0.0138 (95th pct) for a random
feature, which means a real effect of +0.015 is indistinguishable from luck at this sample
size. No amount of cleverness fixes that; only more ligands do. The P450 fold is one
superfamily with one cofactor and one axial cysteine, so a selection term that is real
CYP physics should transfer across it, and a term that only works on the 87 CYP3A4
ligands is probably fitting the 87.

The query is by **Pfam PF00067** rather than by a list of isoform accessions, because
naming the isoforms in advance is exactly the hand-typed list that `rcsb_holo_structures`
was written to avoid. It picks up bacterial P450s (BM3, cam, eryF) too, which is wanted:
they are the same chemistry with a different pocket, and they are the strongest available
test of whether a term generalises or memorises.

Two traps carried over from `build_reference_set.py`, which are not optional here:
multiple copies in the asymmetric unit (pair each ligand to the heme *in its own chain*),
and "the biggest HET group is the ligand" (it often is not).

CIFs stream through a temp directory and are deleted after parsing. A few thousand
entries at 1-3 MB each would be several GB, and the local disks have no room for that.

    python scripts/structure/harvest_p450_universe.py --limit 3000
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
from cypstruct.targets import HEME_ALIASES, IGNORE_HET  # noqa: E402

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
OUT = DATA_PROCESSED / "p450_universe"

COORD_MAX = 2.6      # Fe-donor: direct coordination
SITE_MAX = 6.5       # closest approach: in the active site at all


def p450_entries(limit: int = 3000, max_res: float = 3.2) -> list[str]:
    """All PF00067 entries with a resolution good enough to trust a 2 A distance."""
    q = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_annotation.annotation_id",
                "operator": "exact_match", "value": "PF00067"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less_or_equal", "value": max_res}},
        ]},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit},
                            "results_verbosity": "compact"},
    }
    r = requests.post(SEARCH, json=q, timeout=120)
    r.raise_for_status()
    return sorted(r.json().get("result_set", []))


def parse_entry(pdb_id: str, cif_path: Path) -> list[dict]:
    """One row per (chain, HET group), measured against the heme in that same chain."""
    import gemmi

    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    st.remove_alternative_conformations()
    rows: list[dict] = []

    for chain in st[0]:
        fe = None
        for h in (r for r in chain if r.name.strip().upper() in HEME_ALIASES):
            for at in h:
                if at.element.name.upper() == "FE":
                    fe = np.array([at.pos.x, at.pos.y, at.pos.z])
                    break
            if fe is not None:
                break
        if fe is None:
            continue

        # the proximal thiolate: nearest cysteine SG to the iron. Its presence at ~2.3 A
        # is the check that this really is a P450 active site and not a stray heme.
        sg_d = None
        for res in chain:
            if res.name.strip().upper() == "CYS":
                for at in res:
                    if at.name.strip().upper() == "SG":
                        d = float(np.linalg.norm(
                            np.array([at.pos.x, at.pos.y, at.pos.z]) - fe))
                        if sg_d is None or d < sg_d:
                            sg_d = d

        for res in chain:
            name = res.name.strip().upper()
            if name in HEME_ALIASES or name in IGNORE_HET:
                continue
            # gemmi.Residue has no is_amino_acid() in this build; the tabulated CCD
            # table does, and it also knows about modified residues and nucleotides.
            info = gemmi.find_tabulated_residue(name)
            if info is not None and (info.is_amino_acid() or info.is_water()
                                     or info.is_nucleic_acid()):
                continue
            xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res])
            if len(xyz) < 4:            # ions and single atoms are not ligands
                continue
            els = [a.element.name.upper() for a in res]
            d_fe = np.linalg.norm(xyz - fe, axis=1)
            closest = float(d_fe.min())
            if closest > SITE_MAX:
                continue                # a surface or crystallisation site, not the pocket

            # donor = closest N/O/S; the coordination chemistry the scorer calibrates on
            cands = [i for i, e in enumerate(els) if e in ("N", "O", "S")]
            don_i = min(cands, key=lambda i: d_fe[i]) if cands else None
            rows.append({
                "pdb": pdb_id, "chain": chain.name, "lig": name,
                "n_atoms": int(len(xyz)),
                "n_heavy": int(sum(1 for e in els if e != "H")),
                "closest_fe": closest,
                "donor_elem": els[don_i] if don_i is not None else None,
                "donor_dist": float(d_fe[don_i]) if don_i is not None else None,
                "coordinated": bool(don_i is not None and d_fe[don_i] <= COORD_MAX),
                "cys_sg_fe": sg_d,
            })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=3000)
    ap.add_argument("--max-res", type=float, default=3.2)
    a = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    ids = p450_entries(a.limit, a.max_res)
    print(f"PF00067 entries at <= {a.max_res} A: {len(ids)}", flush=True)
    (OUT / "entries.json").write_text(json.dumps(ids, indent=1))

    done_path = OUT / "rows.jsonl"
    seen = set()
    if done_path.exists():
        for ln in done_path.read_text().splitlines():
            if ln.strip():
                seen.add(json.loads(ln)["pdb"])
    print(f"already parsed: {len(seen)}", flush=True)

    tmp = Path(tempfile.gettempdir()) / "p450_cif"
    tmp.mkdir(exist_ok=True)
    n_new = n_err = 0
    todo = [p for p in ids if p not in seen]
    with done_path.open("a") as fh:
        for i, pid in enumerate(todo):
            f = tmp / f"{pid}.cif"
            try:
                if not f.exists():
                    r = requests.get(
                        f"https://files.rcsb.org/download/{pid}.cif", timeout=120)
                    r.raise_for_status()
                    f.write_bytes(r.content)
                rows = parse_entry(pid, f)
                # An entry with no ligand rows is still recorded, so the resume marker
                # advances and we do not re-download it on every run.
                fh.write(json.dumps({"pdb": pid, "rows": rows}) + "\n")
                fh.flush()
                n_new += 1
            except Exception as exc:
                fh.write(json.dumps({"pdb": pid, "rows": [],
                                     "error": f"{type(exc).__name__}: {exc}"}) + "\n")
                fh.flush()
                n_err += 1
            finally:
                f.unlink(missing_ok=True)     # never accumulate on the local disks
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(todo)}  new={n_new} err={n_err}", flush=True)

    recs = []
    for ln in done_path.read_text().splitlines():
        if ln.strip():
            recs.extend(json.loads(ln).get("rows", []))
    df = pd.DataFrame(recs)
    if len(df):
        df.to_parquet(OUT / "p450_ligands.parquet")
        print(f"\n{len(df)} (chain, ligand) observations "
              f"across {df.pdb.nunique()} entries, {df.lig.nunique()} unique ligands")
        print(f"coordinated: {int(df.coordinated.sum())} "
              f"({100*df.coordinated.mean():.1f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
