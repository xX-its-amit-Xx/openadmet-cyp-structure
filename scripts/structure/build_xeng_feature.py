"""Produce the cross-engine feature the submission builder selects on.

This is the drop-day bridge between a pool of poses and `build_submission.py`. It writes
`data/processed/xeng_<tag>.csv` with {ligand, sample, xeng}, which `choose_poses` picks up
automatically and prefers over the older selectors.

**It refuses rather than degrades.** The feature is worth +0.0381 at >= 4 genuinely
independent reference poses and **-0.0055 at one** - so a silently thin reference set does
not give a weaker selector, it gives a harmful one. `reference_depth` is checked before
anything is written and the run aborts with the numbers if the pool is too shallow.

Three traps, all of them already paid for once:

* `diffusion_samples` does not diversify the ligand on OpenProtein - within one job every
  model shares a single ligand conformation, for Protenix and esmfold2 alike (FINDING 009).
* **Replicates buy NOTHING at all any more** (FINDING 015). Both OpenProtein engines now
  return a byte-identical pose per input: a 12->24 replicate doubling moved the oracle on
  0 of 489 pairs. The lever that still works is the SAMPLER - `num_recycles` / `num_steps`
  give 4 distinct placements from 4 settings (FINDING 016). `reference_poses` still
  dedupes, because historical pools are ~85% duplicates.
* More engines is NOT better. The best measured reference set is the two Protenix
  checkpoints; adding esmfold2 at matched depth drops +0.0380 to +0.0178 (FINDING 011).

Two reference sources, and the sweep is the drop-day one because it can be generated for
any target whereas the engines are frozen:

    # engine reference (historical pools)
    python scripts/structure/build_xeng_feature.py --tag val87b \\
        --pool D:/cyp_scratch/val87b_unsteered \\
        --refs protenix_v2 protenix

    # sampler-sweep reference (FINDING 016) - note the pool layout differs
    python scripts/structure/build_xeng_feature.py --tag p450 \\
        --pool data/processed/p450_universe/poses \\
        --pool-glob '*' --pattern 'protenix_v2/*.cif' \\
        --ref-sweep data/processed/p450_universe/diversity_probe --skip-thin
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

OP_ROOT = DATA_PROCESSED / "openprotein" / "op1"


def main() -> int:
    from cypstruct import pose as P
    from cypstruct import xengine as X

    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--pool", required=True,
                    help="directory of <LIG>__*/ dirs holding the poses to be SELECTED")
    ap.add_argument("--refs", nargs="+", default=["protenix_v2", "protenix"],
                    help="engine pool names under data/processed/openprotein/op1/")
    ap.add_argument("--pattern", default="*.cif")
    ap.add_argument("--pool-glob", default="*__*",
                    help="directory glob inside --pool. The default matches the val87b "
                         "layout (<LIG>__<arm>__<n>); use '*' together with --pattern "
                         "'<engine>/*.cif' for the P450 and sweep layouts.")
    ap.add_argument("--min-depth", type=int, default=4)
    ap.add_argument("--ref-sweep", default=None,
                    help="directory of <pair>/ dirs holding a SAMPLER SWEEP to use as the "
                         "reference instead of engine pools (FINDING 016). This is the "
                         "drop-day path: the sweep can be generated for any target, "
                         "whereas the engines are deterministic and cannot be deepened.")
    ap.add_argument("--frozen",
                    default=str(DATA_PROCESSED / "reference_set_cyp3a4.npz"),
                    help="frozen reference set, used when the raw pools are archived")
    ap.add_argument("--skip-thin", action="store_true",
                    help="omit ligands below --min-depth instead of refusing the whole "
                         "build; those ligands fall back to the older selector")
    ap.add_argument("--force", action="store_true",
                    help="write even if the reference set is too shallow (do not)")
    ap.add_argument("--pool-flat", action="store_true",
                    help="--pool is a FLAT directory of <LIG>__<rest>.cif rather than a "
                         "directory of per-ligand directories. This is what Explorer's "
                         "`boltz_depth.py collect` writes (<lig>__s<seed>_m<k>.cif), and "
                         "without it _emit() silently reads zero poses.")
    ap.add_argument("--ref-npz", default=None,
                    help="use this frozen reference set directly (what "
                         "`openprotein_cofold.py refset` writes for a blind ligand set). "
                         "Takes precedence over --refs and --ref-sweep.")
    a = ap.parse_args()

    if a.ref_npz:
        # The blind drop-day route: the reference was just generated for ligands that have
        # no historical pool at all, so there is nothing under openprotein/op1 to find.
        rp = Path(a.ref_npz)
        if not rp.exists():
            print(f"--ref-npz {rp} does not exist; generate it with "
                  "`openprotein_cofold.py refset`")
            return 1
        ref = X.load_reference(rp)
        print(f"reference: frozen {rp.name} ({len(ref)} ligands)")
        depth = X.reference_depth(ref)
        print(f"reference depth: {  {k: v for k, v in depth.items() if k != 'note'} }")
        thin = [k for k, v in ref.items() if len(v) < a.min_depth]
        if thin and a.skip_thin:
            for k in thin:
                ref.pop(k)
            print(f"skipping {len(thin)} ligands below depth {a.min_depth}: {thin[:10]}")
        elif thin and not a.force:
            print(f"REFUSING: {len(thin)} of {len(ref) + 0} ligands below depth "
                  f"{a.min_depth}: {thin[:10]}")
            print("At one reference pose this feature measured -0.0055 - it would make")
            print("selection WORSE, not weaker. Buy depth with sampler WAVES")
            print("(`openprotein_cofold.py submit --sweep ...`), or pass --skip-thin and")
            print("REPORT which ligands fell back to FINDING 003.")
            return 2
        return _emit(a, ref, P, X)

    # Prefer the raw pools; fall back to the frozen reference set if they have been
    # archived. A 2.2 GB pool deduplicates to ~240 KB of actual reference poses, so the
    # frozen copy is what makes the selector reproducible after the mmCIFs go to cold
    # storage - and stops the pool looking like dead weight that is safe to delete.
    # An ARCHIVED pool leaves an empty directory behind - the collectors recreate it on
    # every poll - and an empty directory still passes `exists()`. Testing existence would
    # therefore take the raw-pool branch, build an EMPTY reference set, and never reach the
    # frozen fallback. Check for actual poses, not for the container.
    if a.ref_sweep:
        # A sweep reference is keyed on the SETTING, not on a replicate index, and
        # num_recycles=1 is excluded because it measured -0.0596 against the default -
        # the one setting that degrades rather than diversifies (FINDING 016).
        root = Path(a.ref_sweep)
        frozen_sweep = DATA_PROCESSED / "reference_set_p450_sweep.npz"
        # An archived sweep leaves empty directories behind, and an empty directory still
        # passes exists() - the same trap that would have built an EMPTY reference set
        # from the archived engine pools. Check for actual poses, and fall back to the
        # frozen copy (1.1 MB, round-trip verified) rather than to nothing.
        if not any(root.glob("*/*.cif")) and frozen_sweep.exists():
            ref = X.load_reference(frozen_sweep)
            print(f"sweep poses absent; using frozen {frozen_sweep.name} "
                  f"({len(ref)} ligands)")
            depth = X.reference_depth(ref)
            print(f"reference depth: {  {k: v for k, v in depth.items() if k != 'note'} }")
            thin = [k for k, v in ref.items() if len(v) < a.min_depth]
            if thin and a.skip_thin:
                for k in thin:
                    ref.pop(k)
                print(f"skipping {len(thin)} ligands below depth {a.min_depth}")
            elif thin and not a.force:
                print(f"REFUSING: {len(thin)} ligands below depth {a.min_depth}; "
                      "pass --skip-thin to build for the rest")
                return 2
            return _emit(a, ref, P, X)
        ref = {}
        for dirp in sorted(root.glob("*")):
            if not dirp.is_dir():
                continue
            lig = dirp.name.split("__")[0]
            vs = []
            for f in sorted(dirp.glob("*.cif")):
                if "1x200" in f.name:
                    continue
                try:
                    v = X.in_heme_frame(P.load_structure(f))
                except Exception:
                    continue
                if v is not None:
                    vs.append(v)
            if vs:
                ref.setdefault(lig, []).extend(vs)
        ref = {k: X._dedupe(v) for k, v in ref.items()}
        if not ref:
            print(f"no sweep poses under {root}")
            return 1
        print(f"reference: sampler sweep from {root} ({len(ref)} ligands)")
        pools = []
    else:
        pools = [OP_ROOT / e for e in a.refs]
    populated = [p for p in pools if any(p.glob("*.cif"))] if pools else []
    if pools and len(populated) == len(pools):
        ref = X.reference_poses(pools, P.load_structure)
    elif not a.ref_sweep and Path(a.frozen).exists():
        ref = X.load_reference(a.frozen)
        print(f"raw pools absent; using frozen reference {Path(a.frozen).name}")
    elif not a.ref_sweep:
        print(f"neither the pools {[str(p) for p in pools]} nor {a.frozen} exist")
        return 1
    depth = X.reference_depth(ref)
    if not a.ref_sweep:
        print(f"reference engines: {a.refs}")
    print(f"reference depth: {  {k: v for k, v in depth.items() if k != 'note'} }")
    thin = [k for k, v in ref.items() if len(v) < a.min_depth]
    if thin and a.skip_thin:
        # Better than refusing the whole build when a handful of ligands are thin: emit
        # the feature only where it is trustworthy and leave the rest out, so
        # build_submission falls back to the older selector for those ligands rather
        # than losing the feature for every ligand. On the full sweep this is 11 of 490.
        for k in thin:
            ref.pop(k)
        print(f"\nskipping {len(thin)} ligands below depth {a.min_depth} "
              f"(e.g. {thin[:5]}); {len(ref)} keep the feature, the rest fall back")
    elif thin and not a.force:
        print(f"\nREFUSING: {len(thin)} of {len(ref)} ligands have fewer than "
              f"{a.min_depth} independent reference poses, e.g. {thin[:5]}.")
        print("At one reference pose this feature measured -0.0055 - it would make the")
        print("selection WORSE, not weaker.")
        # This message used to say "add replicate jobs (not samples)". FINDING 015 made
        # that advice dead: replicates are byte-identical on BOTH OpenProtein engines at
        # any count, so following it would burn jobs and change nothing.
        print("Add SAMPLER SETTINGS (num_recycles / num_steps, FINDING 016) - NOT")
        print("replicates, which are byte-identical on both engines (FINDING 015).")
        print(f"Or pass --skip-thin to build for the other {len(ref) - len(thin)} "
              "ligands and let these fall back.")
        return 2

    return _emit(a, ref, P, X)


def _pool_units(a):
    """(ligand, [pose files]) for either pool layout.

    The flat layout is not cosmetic: Explorer's `collect` writes one directory of
    `<lig>__s<seed>_m<k>.cif`, and `_emit`'s directory glob skips every one of them
    (`if not dirp.is_dir(): continue`), printing "no poses read from the pool" rather than
    anything about layout. Reshaping the directory by hand was the documented workaround
    and it is a copy of the whole pool on a disk with ~12 GB free.
    """
    root = Path(a.pool)
    if a.pool_flat:
        by: dict = {}
        for f in sorted(root.glob(a.pattern)):
            if f.is_file():
                by.setdefault(f.name.split("__")[0], []).append(f)
        return sorted(by.items())
    out = []
    for dirp in sorted(root.glob(a.pool_glob)):
        if dirp.is_dir():
            out.append((dirp.name.split("__")[0], sorted(dirp.glob(a.pattern))))
    return out


def _emit(a, ref, P, X) -> int:
    rows = []
    # The val87b pool names its directories <LIG>__<arm>__<n>; the P450 and sweep pools
    # name them <PAIR> with poses nested under an engine subdirectory. Hardcoding the
    # first layout made the builder silently find zero poses in the second - it printed
    # "no poses read from the pool" rather than anything about layout.
    n_units = n_no_ref = 0
    for lig, files in _pool_units(a):
        n_units += 1
        if lig not in ref or not ref[lig]:
            n_no_ref += 1
            continue
        for f in files:
            try:
                v = X.in_heme_frame(P.load_structure(f))
            except Exception:
                continue
            if v is None:
                continue
            rows.append({"ligand": lig, "sample": f.stem,
                         "xeng": X.xeng_score(v, ref[lig])})
    print(f"pool units: {n_units}   without a reference pose (will fall back): "
          f"{n_no_ref}")
    if not rows:
        print("no poses read from the pool - check --pool, --pattern and --pool-flat "
              "(a flat <LIG>__*.cif directory needs --pool-flat)")
        return 1

    df = pd.DataFrame(rows)
    out = DATA_PROCESSED / f"xeng_{a.tag}.csv"
    df.to_csv(out, index=False)
    print(f"\n{len(df)} poses over {df.ligand.nunique()} ligands -> {out.name}")
    print("build_submission.py will now select on this automatically.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
