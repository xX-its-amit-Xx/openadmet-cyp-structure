"""Does ANY knob still produce pose diversity on OpenProtein?

FINDING 015: both protenix_v2 and esmfold2 now return a byte-identical pose for a given
input, so `--replicates` buys nothing at any count and pose depth is frozen at what is
already on disk. That kills the lever this campaign was built on. Before spending another
job on depth, find out whether anything still varies the output.

The knobs `protenix_v2.fold` exposes are `diffusion_samples` (FINDING 009: does not
diversify the ligand), `num_recycles`, `num_steps`, `templates`. This probe varies the two
untested ones on a handful of pairs and asks the only question that matters: **are the
returned structures distinct?**

Distinctness is judged on the LIGAND, in the heme frame, the same way the selector sees a
pose - not on the file bytes. Two runs can differ in a timestamp and be the same pose, and
they can share a header and place the ligand differently. Bytes were enough to prove
determinism (identical bytes cannot be different poses) but they are not enough to prove
diversity.

    python scripts/cofold/diversity_probe.py --pairs 4
"""
from __future__ import annotations

import argparse
import itertools
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "cofold"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

UNI = DATA_PROCESSED / "p450_universe"

# (num_recycles, num_steps) - the defaults are (3, 200) as this campaign submits them.
# (1, 200) was measured and DROPPED: -0.0596 against the default over 100 pairs. It is the
# only setting that degrades rather than diversifies (FINDING 016).
SETTINGS = [(3, 200), (10, 200), (3, 50), (3, 400)]


def main() -> int:
    from openprotein_cofold import _split_models, build_complex, connect

    from cypstruct import pose as P
    from cypstruct import xengine as X

    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=4)
    ap.add_argument("--settings", default=None,
                    help="comma-separated <recycles>x<steps>, e.g. '5x200,7x200'. "
                         "Defaults to the four validated settings. recycles must "
                         "be >= 2: num_recycles=1 measured -0.0596 (FINDING 016).")
    ap.add_argument("--engine", default="protenix_v2")
    ap.add_argument("--out", default=str(UNI / "diversity_probe"))
    a = ap.parse_args()

    s = connect()
    import json
    st = json.loads((UNI / "campaign.json").read_text())
    ready = {k for k, v in st["msa"].items() if v.get("status") == "SUCCESS"}
    cs = pd.read_csv(UNI / "p450_cofold_set.csv")
    cs["pair"] = cs.pdb + "_" + cs.id
    # pairs this engine has already proved it cannot fold - submitting them again buys
    # nothing and, worse, each doomed job costs the full retrieval budget while waiting
    unsup = UNI / f"{a.engine}_unsupported.json"
    skip = set(json.loads(unsup.read_text()).get("pairs", [])) if unsup.exists() else set()
    cs = cs[cs.target_key.isin(ready)].drop_duplicates("pair")
    if skip:
        before = len(cs)
        cs = cs[~cs.pair.isin(skip)]
        print(f"skipping {before - len(cs)} pairs {a.engine} cannot fold")
    cs = cs.head(a.pairs)
    settings = SETTINGS
    if a.settings:
        settings = []
        for tok in a.settings.split(","):
            rec, steps = tok.strip().lower().split("x")
            rec, steps = int(rec), int(steps)
            if rec < 2:
                print(f"refusing {tok!r}: num_recycles < 2 degrades by -0.0596 against "
                      "the default, it does not diversify (FINDING 016)")
                return 2
            settings.append((rec, steps))
    print(f"probing {len(cs)} pairs x {len(settings)} settings on {a.engine}\n")

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    futs = []
    for r in cs.itertuples():
        msa = s.load_job(st["msa"][r.target_key]["job_id"])
        for rec, steps in settings:
            f = getattr(s.fold, a.engine).fold(
                sequences=[build_complex(r.sequence, r.smiles, msa)],
                diffusion_samples=1, num_recycles=rec, num_steps=steps)
            futs.append((r.pair, rec, steps, f))
            print(f"  {r.pair} recycles={rec} steps={steps} -> {f.job_id}", flush=True)

    print("\nwaiting...")
    got: dict[str, dict] = {}
    failed = 0
    for pair, rec, steps, f in futs:
        # .get() on a job that is still running raises HTTPError, which is NOT a
        # determinism result. The first run of this probe lost 19 of 20 jobs that way
        # and then printed a confident verdict from the one that survived.
        res = None
        for _attempt in range(40):
            # Read the STATUS before retrying. A job that has already FAILED will never
            # become retrievable, and blindly retrying it costs the full 10-minute budget
            # per job - the full-set run stalled for 40 minutes on 2FDY_D4G, a pair
            # already known unfoldable, because this loop could not tell "not finished
            # yet" from "finished badly".
            try:
                status = str(f.job.status).upper()
            except Exception:
                status = ""
            if "FAIL" in status or "CANCEL" in status:
                break
            try:
                res = f.get()
                break
            except Exception:
                time.sleep(15)
        if res is None:
            failed += 1
            print(f"  {pair} r={rec} s={steps}: never retrievable", flush=True)
            continue
        d = out / pair
        d.mkdir(parents=True, exist_ok=True)
        names = _split_models(res[0].to_string(), pair, d, f"{rec}x{steps}")
        for n in names:
            try:
                v = X.in_heme_frame(P.load_structure(d / n))
            except Exception:
                continue
            if v is not None:
                got.setdefault(pair, {})[(rec, steps)] = v

    print(f"\n{'pair':<14}{'settings':>10}  distinct ligand placements")
    verdict = []
    for pair, by in sorted(got.items()):
        vs = list(by.values())
        keys = list(by.keys())
        distinct = X._dedupe(vs)
        # also report the largest pairwise displacement, so "distinct" is not just noise
        worst = 0.0
        for i, j in itertools.combinations(range(len(vs)), 2):
            if vs[i].shape == vs[j].shape:
                worst = max(worst, float(np.abs(vs[i] - vs[j]).max()))
        print(f"{pair:<14}{len(keys):>10}  {len(distinct)} distinct, "
              f"max atom displacement {worst:.3f} A")
        verdict.append(len(distinct))

    # A verdict needs enough settings to have actually varied something. The first run of
    # this probe concluded "no lever remains" from ONE pair with ONE setting, because 19
    # of 20 jobs raised HTTPError on retrieval. That cannot tell a deterministic engine
    # apart from a failed collection, and it is the more dangerous of the two errors
    # because the wrong answer is the one that sounds like a finding.
    usable = [p_ for p_, by in got.items() if len(by) >= 3]
    print(f"\npairs with >= 3 settings retrieved: {len(usable)} of {len(cs)}"
          f"   (jobs never retrievable: {failed})")
    if len(usable) < 2:
        print("REFUSING to conclude: too few settings came back to tell a deterministic")
        print("engine apart from a failed collection. Re-run when the queue is quiet.")
        return 2
    verdict = [len(X._dedupe(list(got[p_].values()))) for p_ in usable]
    med = float(np.median(verdict))
    print(f"median distinct placements across {len(settings)} settings: {med:.1f}")
    if med <= 1:
        print("VERDICT: num_recycles / num_steps do NOT diversify either. "
              "No replicate-style lever remains on this venue.")
    else:
        print(f"VERDICT: these knobs DO diversify - {med:.1f} distinct from "
              f"{len(settings)} settings. This is a usable depth lever.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
