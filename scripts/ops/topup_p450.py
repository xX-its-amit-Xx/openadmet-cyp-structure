"""Top up the P450 campaign: one command the ops tick can run without thinking.

The generalisation set (FINDING 012) is the only experiment still producing new
information, and it grows in two ways - new proteins as MSAs finish, and deeper pools on
proteins already covered. Both are just `p450_campaign.py submit` with the right engine
and replicate count, and both are idempotent, but getting them right each cycle by hand is
how a step gets skipped at 3am.

What it does, in the order that matters:

1. **Protenix-v2 to depth** - this is the POOL being selected from. Pairs with fewer than
   3 poses are dropped from the generalisation test entirely, so thin pairs cost whole
   proteins: 48 thin pairs were gating 10 proteins when last measured.
2. **esmfold2 to depth** - this is the REFERENCE. Below 2 independent poses a pair cannot
   be scored at all, and below 4 the feature is unreliable (at one reference pose it
   measured -0.0055). The organometallics it cannot fold are skipped automatically.

It deliberately does NOT submit rosettafold_3 or protenix-v1: RF3 is half-deterministic in
single-sequence mode and places organometallics wrongly, and protenix-v1 is fully
deterministic, so replicates of either buy nothing (FINDING 011).

    python scripts/ops/topup_p450.py            # report the gaps, submit nothing
    python scripts/ops/topup_p450.py --submit
"""
from __future__ import annotations

import argparse
import collections
import glob
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

UNI = DATA_PROCESSED / "p450_universe"
POSES = UNI / "poses"


def replicate_counts(engine: str) -> dict[str, int]:
    out: dict[str, set] = collections.defaultdict(set)
    for p in glob.glob(str(POSES / "*" / engine / "*.cif")):
        pair = p.split(os.sep)[-3]
        m = re.match(r"(.+)__r(\d+)s(\d+)\.cif$", os.path.basename(p))
        if m:
            out[pair].add(int(m.group(2)))
    return {k: len(v) for k, v in out.items()}


def refresh_skip_lists() -> dict:
    """Learn which pairs each engine cannot fold, from the failures it already produced.

    Every engine tested so far has ligands it chokes on, and they are not the same set:
    esmfold2 refuses the Ir/Ru organometallics, and protenix_v2 turned out to fail on
    2FDW_D3G / 2FDY_D4G. Without this, each new replicate round re-submits the same
    doomed pairs - 8 protenix jobs and 15 esmfold2 jobs were burned that way before
    anyone looked at WHICH pairs were failing rather than how many.

    A pair is only listed once it has failed at least twice, so a single transient
    server-side error does not permanently exclude a pair that would otherwise fold.
    """
    import json

    state = UNI / "campaign.json"
    if not state.exists():
        return {}
    st = json.loads(state.read_text())
    fails: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for rec in st.get("folds", {}).values():
        if rec.get("done") != "failed":
            continue
        eng = rec.get("engine", "protenix_v2")
        for pair in rec["pairs"]:
            fails[eng][pair] += 1
    out = {}
    for eng, counter in fails.items():
        repeat = {p: c for p, c in counter.items() if c >= 2}
        if not repeat:
            continue
        path = UNI / f"{eng}_unsupported.json"
        prev = set(json.loads(path.read_text()).get("pairs", [])) if path.exists() else set()
        path.write_text(json.dumps({"pairs": sorted(repeat), "counts": repeat}, indent=1))
        out[eng] = {"pairs": len(repeat), "new": len(set(repeat) - prev)}
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--pool-replicates", type=int, default=12)
    ap.add_argument("--ref-replicates", type=int, default=6)
    ap.add_argument("--limit", type=int, default=150)
    a = ap.parse_args()

    skips = refresh_skip_lists()
    for eng, info in sorted(skips.items()):
        print(f"  {eng}: cannot fold {info['pairs']} pairs"
              + (f" (+{info['new']} newly learned)" if info["new"] else ""))

    scored = UNI / "p450_poses_scored.csv"
    if scored.exists():
        sc = pd.read_csv(scored)
        pool = sc[sc.engine == "protenix_v2"].groupby("pair").size()
        thin = set(pool[pool < 3].index)
        gated = sc[sc.pair.isin(thin)].uniprot.nunique()
        print(f"scored: {sc.pair.nunique()} pairs, {sc.uniprot.nunique()} proteins, "
              f"{sc.target_key.nunique()} construct sequences")
        print(f"  pairs with <3 protenix poses (dropped from the test): {len(thin)}"
              f"  -> gating {gated} proteins")
    pv = replicate_counts("protenix_v2")
    es = replicate_counts("esmfold2")
    print(f"  protenix_v2 replicates/pair: median "
          f"{int(pd.Series(pv).median()) if pv else 0}, "
          f"pairs below {a.pool_replicates}: {sum(1 for v in pv.values() if v < a.pool_replicates)}")
    print(f"  esmfold2 replicates/pair:    median "
          f"{int(pd.Series(es).median()) if es else 0}, "
          f"pairs below {a.ref_replicates}: {sum(1 for v in es.values() if v < a.ref_replicates)}")

    if not a.submit:
        print("\n(dry run - pass --submit to queue the work)")
        return 0

    runner = REPO / "scripts" / "cofold" / "p450_campaign.py"
    for engine, reps in (("protenix_v2", a.pool_replicates),
                         ("esmfold2", a.ref_replicates)):
        cmd = [sys.executable, str(runner), "submit", "--engine", engine,
               "--samples", "1", "--batch", "2", "--replicates", str(reps),
               "--limit", str(a.limit)]
        print(f"\n$ {' '.join(cmd[1:])}", flush=True)
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        tail = [ln for ln in r.stdout.splitlines() if "->" in ln or "skipping" in ln]
        print(f"  submitted {sum(1 for ln in tail if '->' in ln)} jobs"
              + (f"; {[ln for ln in tail if 'skipping' in ln][0]}"
                 if any('skipping' in ln for ln in tail) else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
