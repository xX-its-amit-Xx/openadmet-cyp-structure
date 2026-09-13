"""Turn the harvested P450 pairs into something that can actually be co-folded.

FINDING 008 produced 496 (entry, ligand) pairs with crystal ground truth, which is 5.7x
the 87 CYP3A4 ligands that every selection experiment has been stuck on. But a pair is
only useful as validation if we can *predict* it, and that needs two things the harvest
did not collect:

  * **the ligand as SMILES**, not just a three-letter CCD code, and
  * **the protein sequence of the entity the ligand is actually bound to** - not the
    canonical UniProt sequence. Crystallised P450s are truncated, mutated and tagged, and
    folding the wild-type sequence against a pose solved with an engineered construct
    quietly compares two different proteins.

Both come from the RCSB entry API, one request per entry, cached to disk.

The output is a csv shaped exactly like `validation_ligands.csv`, so everything already
written against that file works unchanged.

**Deliberately not done here:** MSAs. Each distinct sequence needs its own alignment, and
OpenProtein computes those server-side. That is a separate, slower step - this script only
has to produce the inputs for it.

    python scripts/structure/build_p450_cofold_set.py
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

OUT = DATA_PROCESSED / "p450_universe"
GRAPHQL = "https://data.rcsb.org/graphql"

# One query for everything about an entry: every polymer entity with its one-letter
# sequence, and every non-polymer entity with its CCD code and SMILES. Batched, because
# 500 entries at one REST call each is 500 round trips.
Q = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    polymer_entities{
      rcsb_id
      entity_poly{ pdbx_seq_one_letter_code_can }
      rcsb_polymer_entity_container_identifiers{ auth_asym_ids }
      rcsb_polymer_entity_align{ reference_database_accession }
    }
    nonpolymer_entities{
      rcsb_id
      rcsb_nonpolymer_entity_container_identifiers{ auth_asym_ids }
      nonpolymer_comp{
        chem_comp{ id name formula_weight }
        rcsb_chem_comp_descriptor{ SMILES_stereo }
      }
    }
  }
}
"""


def fetch(ids: list[str], chunk: int = 40) -> dict:
    """Entry metadata, cached. Re-running after a partial failure costs nothing."""
    cache_path = OUT / "entry_meta.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    todo = [i for i in ids if i not in cache]
    print(f"{len(todo)} entries to fetch ({len(cache)} cached)", flush=True)

    for i in range(0, len(todo), chunk):
        part = todo[i:i + chunk]
        try:
            r = requests.post(GRAPHQL, json={"query": Q, "variables": {"ids": part}},
                              timeout=120)
            r.raise_for_status()
            data = r.json()
            if data.get("errors"):
                print("  graphql errors:", str(data["errors"])[:300], flush=True)
            for e in (data.get("data", {}).get("entries") or []):
                if e:
                    cache[e["rcsb_id"]] = e
            cache_path.write_text(json.dumps(cache))
            print(f"  {min(i+chunk, len(todo))}/{len(todo)}", flush=True)
        except Exception as exc:
            print(f"  chunk {i} failed: {type(exc).__name__}: {exc}", flush=True)
        time.sleep(0.3)
    return cache


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-heavy", type=int, default=12)
    ap.add_argument("--max-len", type=int, default=600,
                    help="skip constructs too long to co-fold cheaply")
    a = ap.parse_args()

    geo = pd.read_parquet(OUT / "p450_geometry.parquet")
    geo = geo[geo.closest_fe <= 6.5]                 # in the pocket, not a surface site
    # one row per (entry, ligand): extra chains are copies of the same complex and would
    # weight a multi-chain entry more heavily than a monomeric one in every average
    geo = geo.sort_values("closest_fe").drop_duplicates(["pdb", "lig"])
    print(f"{len(geo)} in-pocket (entry, ligand) pairs", flush=True)

    meta = fetch(sorted(geo.pdb.unique()))

    rows, skipped = [], {}
    for r in geo.itertuples():
        e = meta.get(r.pdb)
        if not e:
            skipped["no metadata"] = skipped.get("no metadata", 0) + 1
            continue

        # the ligand: SMILES by CCD code
        smiles = None
        for ne in (e.get("nonpolymer_entities") or []):
            comp = (ne.get("nonpolymer_comp") or {})
            cc = (comp.get("chem_comp") or {})
            if (cc.get("id") or "").strip().upper() == r.lig:
                smiles = ((comp.get("rcsb_chem_comp_descriptor") or {})
                          .get("SMILES_stereo"))
                break
        if not smiles:
            skipped["no smiles"] = skipped.get("no smiles", 0) + 1
            continue

        # the protein: the LONGEST polymer entity, which for a P450 entry is the P450
        # itself rather than a redox partner, nanobody or crystallisation chaperone
        best_seq, best_acc = None, None
        for pe in (e.get("polymer_entities") or []):
            seq = ((pe.get("entity_poly") or {})
                   .get("pdbx_seq_one_letter_code_can") or "").replace("\n", "")
            if not seq or (best_seq and len(seq) <= len(best_seq)):
                continue
            best_seq = seq
            al = pe.get("rcsb_polymer_entity_align") or []
            best_acc = (al[0] or {}).get("reference_database_accession") if al else None
        if not best_seq:
            skipped["no sequence"] = skipped.get("no sequence", 0) + 1
            continue
        if len(best_seq) > a.max_len:
            skipped["too long"] = skipped.get("too long", 0) + 1
            continue
        if set(best_seq) - set("ACDEFGHIKLMNPQRSTVWY"):
            skipped["nonstandard residues"] = skipped.get("nonstandard residues", 0) + 1
            continue

        rows.append({
            "id": r.lig, "pdb": r.pdb, "chain": r.chain, "seqid": r.seqid,
            "smiles": smiles, "sequence": best_seq, "seq_len": len(best_seq),
            "uniprot": best_acc, "n_heavy": r.n_heavy,
            "coordinated": r.coordinated, "donor_dist": r.donor_dist,
            "closest_fe": r.closest_fe,
        })

    df = pd.DataFrame(rows)
    # md5, NOT the builtin hash(): Python salts string hashing per process, so hash()
    # would mint a different target_key every run and silently orphan the MSA jobs
    # already submitted against the previous keys.
    import hashlib
    df["target_key"] = df.sequence.map(
        lambda s: f"{len(s)}_{hashlib.md5(s.encode()).hexdigest()[:8]}")
    out = OUT / "p450_cofold_set.csv"
    df.to_csv(out, index=False)

    print(f"\n{len(df)} foldable pairs, {df.id.nunique()} unique ligands, "
          f"{df.target_key.nunique()} unique protein sequences")
    print(f"coordinated: {int(df.coordinated.sum())} ({100*df.coordinated.mean():.0f}%)")
    print(f"skipped: {skipped}")
    print(f"\n-> {out}")
    print("\nmost common targets:")
    print(df.groupby("uniprot").size().sort_values(ascending=False).head(8).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
