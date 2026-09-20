"""Convert the staged a3m MSAs into the .npz format Boltz's loader expects.

Runs ON EXPLORER. `load_input` does `np.load(msa_dir / f"{msa_id}.npz")` and feeds the
arrays straight into `MSA(**msa)` - it does NOT parse a3m. Staging 185 a3m files was
necessary but not sufficient; without this step every protein chain silently trains
single-sequence, which is exactly the failure mode that looks like a bad recipe.

File names are the `target_key` column from the arm CSVs (`<seq_len>_<hash>`), so the
record's msa_id can be that key verbatim.

    ./env/bin/python make_msa_npz.py
"""
from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path

ROOT = Path("/scratch/shenoy.am/cyp-finetune")
A3M = ROOT / "msa"
OUT = ROOT / "msa_npz"


def main() -> int:
    from boltz.data.parse.a3m import parse_a3m

    ap = argparse.ArgumentParser()
    ap.add_argument("--max-seqs", type=int, default=16384)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    srcs = sorted(A3M.glob("*.a3m"))
    ok = skip = fail = 0
    reasons: dict[str, int] = {}
    for src in srcs:
        dst = OUT / f"{src.stem}.npz"
        if dst.exists() and not a.force:
            skip += 1
            continue
        try:
            # taxonomy=None: pairing across species is for multimer complexes built from
            # different entities. Our targets are single-entity (homo-oligomers at most),
            # so there is nothing to pair and the taxonomy db would only cost a lookup.
            msa = parse_a3m(src, taxonomy=None, max_seqs=a.max_seqs)
            msa.dump(dst)
            ok += 1
        except Exception as exc:  # noqa: BLE001
            fail += 1
            k = f"{type(exc).__name__}: {str(exc)[:70]}"
            reasons[k] = reasons.get(k, 0) + 1
            if fail <= 2:
                traceback.print_exc(limit=2)
        if (ok + fail) % 25 == 0 and (ok + fail):
            print(f"  {ok + fail + skip}/{len(srcs)}", flush=True)

    print(json.dumps({
        "a3m_found": len(srcs), "converted": ok, "already_done": skip,
        "failed": fail, "reasons": dict(sorted(reasons.items(), key=lambda kv: -kv[1])[:4]),
        "out": str(OUT),
    }, indent=2))
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
