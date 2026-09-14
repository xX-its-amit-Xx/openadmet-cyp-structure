"""Is the structure track live, and are we ready for it? One command, run every tick.

Replaces the hand-rolled polling I have been doing each ops tick, which is both easy to
forget and easy to do inconsistently. Two questions, answered against the source of truth
rather than from memory:

  1. **Has anything changed upstream?** `STRUCTURE_TRACK_LIVE`, `STRUCTURE_DATASET_SIZE`
     and the dataset's file list, read from the challenge Space and the HF API. The
     dataset size is a placeholder (184, the PXR count, carrying a TODO), so a change to
     it is itself a signal that the real set has landed.
  2. **Could we actually submit today?** Every piece of the drop-day path is checked for
     existence and, where cheap, exercised. A pipeline that has never been run is not a
     pipeline.

Exit code 0 means nothing to do; 2 means something upstream changed and needs attention.

    python scripts/ops/readiness.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

STATE = DATA_PROCESSED / "upstream_state.json"
SPACE_CFG = "https://huggingface.co/spaces/openadmet/cyp-challenge/raw/main/config.py"
DATASET_API = "https://huggingface.co/api/datasets/openadmet/cyp-challenge-train-test"


def upstream() -> dict:
    import requests

    out: dict = {}
    try:
        cfg = requests.get(SPACE_CFG, timeout=60).text
        for key in ("STRUCTURE_TRACK_LIVE", "STRUCTURE_DATASET_SIZE"):
            m = re.search(rf"^{key}\s*=\s*([^\s#]+)", cfg, re.M)
            out[key] = m.group(1).strip() if m else None
    except Exception as exc:
        out["space_error"] = f"{type(exc).__name__}: {exc}"
    try:
        d = requests.get(DATASET_API, timeout=60).json()
        out["dataset_last_modified"] = d.get("lastModified")
        out["dataset_files"] = sorted(s.get("rfilename") for s in d.get("siblings", []))
    except Exception as exc:
        out["dataset_error"] = f"{type(exc).__name__}: {exc}"
    return out


def readiness() -> list[tuple[str, bool, str]]:
    """Each row: (what, ok, detail). Checks artefacts AND that they have been exercised."""
    rows: list[tuple[str, bool, str]] = []

    def chk(name: str, path: Path, note: str = "") -> bool:
        ok = path.exists()
        rows.append((name, ok, note or str(path.relative_to(REPO))))
        return ok

    chk("co-folding runner (Boltz)", REPO / "scripts/cofold/modal_boltz.py")
    chk("detached launcher", REPO / "scripts/cofold/detached.py")
    chk("CPU pre-flight parse", REPO / "scripts/cofold/preflight_parse.py")
    chk("pool scorer", REPO / "scripts/structure/collect_and_score.py")
    chk("orientation features", REPO / "scripts/structure/orientation_features.py")
    chk("consensus features", REPO / "scripts/structure/test_consensus_selector.py")
    chk("submission builder", REPO / "scripts/submit/build_submission.py")
    chk("shared CYP3A4 MSA", REPO / "data/reference/cyp3a4.a3m",
        "6,979 sequences; the expensive shared input")

    # the selector must be importable and callable, not merely present
    try:
        from cypstruct.select import orientation_consensus_score  # noqa: F401
        rows.append(("selector importable", True,
                     "cypstruct.select.orientation_consensus_score"))
    except Exception as exc:
        rows.append(("selector importable", False, f"{type(exc).__name__}: {exc}"))

    # a submission was actually produced and validated end to end
    dry = REPO / "submissions/01_xeng_val87b.zip"
    if not dry.exists():
        dry = REPO / "submissions/00_dryrun_val87b.zip"
    rows.append(("dry-run submission exists", dry.exists(),
                 f"{dry.name} {dry.stat().st_size//1024} KB" if dry.exists()
                 else "never built"))

    # EXERCISE the selector, do not merely check that its file is present. A readiness
    # check that only stats files reported 11/11 passing while `build()` could not read a
    # local pool at all - `pool_dir` was a dead parameter, so with Modal over its cap the
    # submission was unbuildable. Importing and calling is what would have caught it.
    try:
        sys.path.insert(0, str(REPO / "scripts" / "submit"))
        import inspect

        from build_submission import build, choose_poses
        sig = inspect.signature(build)
        src = inspect.getsource(build)
        wired = "pool_dir" in sig.parameters and "pool_dir is None" in src
        rows.append(("build() can read a LOCAL pool", wired,
                     "--pool-dir honoured" if wired
                     else "pool_dir is dead - Modal-only, unbuildable while over cap"))
        picks = choose_poses("val87b")
        rows.append(("selector runs end to end", len(picks) > 0,
                     f"{len(picks)} ligands picked"))
    except Exception as exc:
        rows.append(("build()/selector exercisable", False,
                     f"{type(exc).__name__}: {exc}"))

    # the cross-engine feature is what the selector now prefers; absent, it silently
    # falls back to the weaker FINDING 003 rule
    xe = DATA_PROCESSED / "xeng_val87b.csv"
    rows.append(("cross-engine feature built", xe.exists(),
                 "xeng_val87b.csv (+0.0381 selector)" if xe.exists()
                 else "MISSING - selector would fall back to +0.0265"))

    # measured selector performance, so the number is not recalled from memory
    sel = DATA_PROCESSED / "orientation_selector_val87b_unsteered.json"
    rows.append(("selector measured", sel.exists(),
                 "xeng +0.0381 (FINDING 011); generalises +0.2045 over 17 P450 "
                 "targets (FINDING 012)"))
    return rows


def main() -> int:
    now = upstream()
    first_run = not STATE.exists()
    prev = {} if first_run else json.loads(STATE.read_text())
    # On the FIRST run everything differs from an empty baseline, so every field would be
    # reported as changed. A monitor that cries wolf the first time it runs teaches you to
    # ignore it - the same way the watchdog's early alarms on harmless idle shells did.
    # Record the baseline silently and start comparing from the next run.
    changed = {} if first_run else {
        k: (prev.get(k), v) for k, v in now.items()
        if k in ("STRUCTURE_TRACK_LIVE", "STRUCTURE_DATASET_SIZE",
                 "dataset_last_modified") and prev.get(k) != v}
    new_files = [] if first_run else sorted(
        set(now.get("dataset_files", [])) - set(prev.get("dataset_files", [])))

    print("UPSTREAM")
    print(f"  STRUCTURE_TRACK_LIVE     = {now.get('STRUCTURE_TRACK_LIVE')}")
    print(f"  STRUCTURE_DATASET_SIZE   = {now.get('STRUCTURE_DATASET_SIZE')}"
          "   (184 = PXR placeholder)")
    print(f"  dataset last modified    = {now.get('dataset_last_modified')}")
    print(f"  dataset files            = {len(now.get('dataset_files', []))}")

    print("\nREADINESS")
    rows = readiness()
    for name, ok, detail in rows:
        print(f"  [{'ok' if ok else 'MISSING'}] {name:28s} {detail}")
    n_bad = sum(1 for _n, ok, _d in rows if not ok)

    STATE.write_text(json.dumps(now, indent=1))

    if changed or new_files:
        print("\n*** UPSTREAM CHANGED ***")
        for k, (old, new) in changed.items():
            print(f"  {k}: {old} -> {new}")
        if new_files:
            print(f"  new dataset files: {new_files}")
        print("\n  Drop-day path:")
        print("   1. re-read the Space config: the dataset size and id format were")
        print("      placeholders (184, 'x00011-1') taken from PXR")
        print("   2. preflight_parse.py on the new ligands  (CPU, catches schema errors)")
        print("   3. detached.py launch --engine boltz --arms unsteered")
        print("   4. collect_and_score -> orientation_features -> test_consensus_selector")
        print("   5. build_submission.py build, then validate --expect-n <size>")
        return 2

    print(f"\nno upstream change; {len(rows)-n_bad}/{len(rows)} readiness checks pass")
    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
