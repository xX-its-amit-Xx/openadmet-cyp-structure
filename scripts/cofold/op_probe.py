"""Answer, in one round trip, every unknown that gates the OpenProtein campaign.

Four questions, and each one changes what the campaign looks like:

  1. Does attaching the 6,979-sequence MSA break protenix-v2? The one successful fold did
     not record whether it carried an MSA, and if the MSA is what breaks it we would lose
     a whole tranche discovering that serially.
  2. Does `diffusion_samples=N` return N structures, or does the result list index over
     input complexes? FINDING 004 says the oracle is still climbing at 20 samples, so how
     samples are addressed decides whether this venue can supply a pool at all.
  3. Can one job carry several complexes? That is the difference between 18 jobs and 87.
  4. Which of the remaining nine engines accept protein + CCD heme + SMILES ligand?
     boltz2 and rosettafold-3 already fail; the others are untested.

Submitted together and left to run, because they are independent and each takes minutes.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

sys.path.insert(0, str(REPO / "scripts" / "cofold"))
from openprotein_cofold import build_complex, connect, get_msa  # noqa: E402

OUT = DATA_PROCESSED / "openprotein" / "probe.json"
# two real validation ligands, so a success is a usable pose and not just a smoke test
LIGS = [("1RD", "CC(C)(C)c1ccc(cc1)C(=O)Nc1ccc(cc1)n1ccnc1"),
        ("KLN", "Cn1cnc(c1)Cc1ccccc1")]


def main() -> int:
    import pandas as pd

    df = pd.read_csv(DATA_PROCESSED / "validation_ligands.csv")
    pairs = [(r.id, r.smiles) for r in df.head(2).itertuples()]

    s = connect()
    from cypstruct.targets import fetch_sequences
    seq = fetch_sequences()["cyp3a4"]
    msa = get_msa(s)

    def cx(smi, with_msa):
        return build_complex(seq, smi, msa if with_msa else None)

    # Every arm carries the MSA: the API REQUIRES one (or explicit single-sequence mode),
    # so the first round's "nomsa" arms were invalid requests rather than engine failures.
    # alphafold2 is dropped entirely - it warns that it ignores ligand chains, which makes
    # it useless for co-folding no matter whether it runs.
    probes = [
        ("protenix_v2__msa__s3", "protenix_v2", [cx(pairs[0][1], True)], 3),
        ("protenix_v2__batch2",  "protenix_v2", [cx(p[1], True) for p in pairs], 1),
        ("protenix__msa",        "protenix",    [cx(pairs[0][1], True)], 1),
        ("boltz_1x__msa",        "boltz_1x",    [cx(pairs[0][1], True)], 1),
        ("boltz_1__msa",         "boltz_1",     [cx(pairs[0][1], True)], 1),
        ("boltz_2__msa",         "boltz2",      [cx(pairs[0][1], True)], 1),
        ("esmfold2__msa",        "esmfold2",    [cx(pairs[0][1], True)], 1),
        ("rosettafold_3__msa",   "rosettafold_3", [cx(pairs[0][1], True)], 1),
    ]

    rec = json.loads(OUT.read_text()) if OUT.exists() else {}
    for name, engine, seqs, nsamp in probes:
        if name in rec:
            print(f"  {name:26s} already submitted {rec[name].get('job_id')}", flush=True)
            continue
        try:
            model = getattr(s.fold, engine)
            fut = model.fold(sequences=seqs, diffusion_samples=nsamp, num_recycles=3)
            rec[name] = {"job_id": str(fut.job_id), "engine": engine,
                         "n_seq": len(seqs), "samples": nsamp, "t": time.time(),
                         "ligands": [p[0] for p in pairs][:len(seqs)]}
            print(f"  {name:26s} -> {fut.job_id}", flush=True)
        except Exception as exc:
            rec[name] = {"error": f"{type(exc).__name__}: {exc}"}
            print(f"  {name:26s} SUBMIT-FAIL {type(exc).__name__}: {exc}", flush=True)
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(rec, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
