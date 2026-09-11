"""Launch a co-folding batch DETACHED, so it survives the client that started it.

**The problem this solves, which has now bitten five times.** Every long-running local
process in this project dies when its turn ends: a Modal launcher, an OneDrive archive
twice, a Chai smoke test, and an 87-job Chai batch whose gate job succeeded and whose map
never ran because the client was killed 590 seconds in. `run_in_background` does not help
— the process is gone by the next turn either way.

Modal work started through `with app.run()` is **ephemeral**: it is tied to the client, so
when the client dies the app is torn down with it. That is why the 168-job Boltz batch
survived a full session teardown (it had already finished mapping) while everything
launched near a turn boundary did not.

**The fix is to deploy.** A deployed Modal app is server-side and permanent. `spawn_map`
hands it the whole input list and returns immediately; the work proceeds whether or not
anything local is alive. Progress is read back by counting `DONE.json` markers on the
output volume, so status is a property of the work rather than of a log file.

    python scripts/cofold/detached.py launch --engine chai --csv <csv> --tag chai87 --samples 10
    python scripts/cofold/detached.py status --engine chai --tag chai87
    python scripts/cofold/detached.py stop   --engine chai

Because every job is idempotent against its `DONE.json`, relaunching the same tag is safe
and resumes rather than repeats.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "cofold"))

import modal  # noqa: E402

ENGINES = {
    "chai": {"module": "modal_chai", "app": "cyp-cofold-chai", "fn": "cofold",
             "kind": "chai_cofold"},
    "boltz": {"module": "modal_boltz", "app": "cyp-cofold-boltz", "fn": "cofold",
              "kind": "boltz2_cofold"},
}


def _deploy(module: str) -> str:
    """`modal deploy` the module so its app becomes server-side and permanent."""
    path = REPO / "scripts" / "cofold" / f"{module}.py"
    # `modal deploy` prints a U+2713 check mark. On this Windows box the child inherits a
    # cp1252 stdout and dies encoding it, so the deploy "fails" for reasons that have
    # nothing to do with the deploy. Force UTF-8 on the child and decode leniently here.
    import os as _os

    env = {**_os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}
    cp = subprocess.run(["modal", "deploy", str(path)],
                        capture_output=True, timeout=1800, env=env)
    out = (cp.stdout or b"").decode("utf-8", "replace")
    err = (cp.stderr or b"").decode("utf-8", "replace")
    cp = subprocess.CompletedProcess(cp.args, cp.returncode, out, err)
    if cp.returncode != 0:
        raise RuntimeError(f"modal deploy failed:\n{cp.stdout[-1500:]}\n{cp.stderr[-2500:]}")
    return (cp.stdout or "")[-800:]


def _done_count(tag: str) -> tuple[int, int]:
    """(#jobs with DONE.json, #job dirs) on the output volume. Progress, not liveness."""
    vol = modal.Volume.from_name("cyp-pool")
    try:
        entries = list(vol.iterdir(f"/{tag}"))
    except Exception:
        return 0, 0
    jobs = sorted({e.path.rstrip("/").split("/")[-1] for e in entries})
    jobs = [j for j in jobs if "__" in j]
    done = 0
    for j in jobs:
        try:
            next(x for x in vol.iterdir(f"/{tag}/{j}")
                 if x.path.endswith("DONE.json"))
            done += 1
        except Exception:
            pass
    return done, len(jobs)


def launch(engine: str, csv: str, tag: str, samples: int, seeds: str,
           arms: str | None) -> None:
    from cypstruct import budget

    spec = ENGINES[engine]
    mod = __import__(spec["module"])

    seed_list = [int(s) for s in seeds.split(",")]
    kwargs = dict(csv_path=csv, tag=tag, samples=samples, seeds=seed_list)
    # BOTH engines take `arms`. Passing it only to Chai meant a Boltz launch silently
    # planned the STEERED arm as well, doubling a 9.6 GPU-h run to 18.5 - on an arm
    # FINDING 001 measured as a null. Boltz defaults to unsteered here for that reason;
    # ask for the steered arm explicitly if you ever want to re-measure it.
    if arms is None and engine == "boltz":
        arms = "unsteered"
    if arms is not None:
        kwargs["arms"] = tuple(a.strip() for a in arms.split(",") if a.strip())
    jobs = mod.plan(**kwargs)
    if not jobs:
        raise SystemExit(f"plan() produced no jobs for arms={arms!r}")

    # Skip what is already finished, so a relaunch resumes instead of repeating.
    done_before, _ = _done_count(tag)
    est = budget.estimate(spec["kind"], len(jobs), samples)
    ok, why, _ = budget.preflight_hours(spec["kind"], est)
    print(f"{len(jobs)} jobs planned ({done_before} already have DONE.json)  "
          f"est {est:.2f} GPU-h", flush=True)
    if not ok:
        raise SystemExit(f"PREFLIGHT REFUSED: {why}")

    print("deploying the app so the run outlives this process ...", flush=True)
    print("  ", _deploy(spec["module"]).strip().splitlines()[-1:], flush=True)

    fn = modal.Function.from_name(spec["app"], spec["fn"])
    run_id = f"{tag}-detached-{int(time.time())}"
    budget.record(run_id, "modal", spec["kind"], len(jobs) * samples, est,
                  note=f"detached spawn_map app={spec['app']} tag={tag}",
                  n_jobs=len(jobs), samples=samples, app=spec["app"], detached=True)

    # spawn_map returns immediately; the work is server-side from here on.
    if engine == "chai":
        # Chai's cofold takes a LIST of specs per call, so batch them.
        chunk = 4
        batches = [jobs[i:i + chunk] for i in range(0, len(jobs), chunk)]
        fn.spawn_map(batches)
        n_calls = len(batches)
    else:
        fn.spawn_map(jobs)
        n_calls = len(jobs)

    print(f"spawned {n_calls} detached call(s). This process can now exit safely.",
          flush=True)
    print(f"check with:  python scripts/cofold/detached.py status "
          f"--engine {engine} --tag {tag}", flush=True)


CAPACITY_MARKERS = ("waiting to be scheduled", "acquiring more capacity")


def waiting_for_capacity(app_id: str) -> str | None:
    """Is this app queued for a GPU rather than broken?

    Both batches once sat at zero running containers with their output counts frozen,
    which reads as a stall. The logs said otherwise: "waiting to be scheduled on a
    GPU_A100 worker". Queued is not stalled, and killing a queued run throws away its
    place in the queue along with the work. Always check before judging.
    """
    try:
        cp = subprocess.run(["modal", "app", "logs", app_id],
                            capture_output=True, text=True, timeout=90)
    except Exception:
        return None
    tail = ((cp.stdout or "") + (cp.stderr or ""))[-6000:]
    for line in reversed(tail.splitlines()):
        if any(m in line for m in CAPACITY_MARKERS):
            return line.strip()[:200]
    return None


def status(engine: str, tag: str) -> None:
    spec = ENGINES[engine]
    done, total = _done_count(tag)
    print(json.dumps({"tag": tag, "engine": engine, "app": spec["app"],
                      "jobs_with_output": done, "job_dirs_present": total}, indent=1))
    cp = subprocess.run(["modal", "app", "list", "--json"],
                        capture_output=True, text=True, timeout=180)
    try:
        apps = [a for a in json.loads(cp.stdout)
                if a.get("description") == spec["app"] and not a.get("stopped_at")]
        for a in apps:
            note = ""
            if str(a.get("tasks")) == "0":
                w = waiting_for_capacity(a["app_id"])
                note = f"  <- QUEUED FOR CAPACITY: {w}" if w else "  <- idle"
            print(f"  app {a['app_id']} state={a['state']} tasks={a['tasks']} "
                  f"created={a['created_at']}{note}")
        if not apps:
            print("  (no live app — finished, or never deployed)")
    except Exception:
        pass


def stop(engine: str) -> None:
    """Stop the deployed app. Deployed apps persist until told otherwise."""
    spec = ENGINES[engine]
    cp = subprocess.run(["modal", "app", "list", "--json"],
                        capture_output=True, text=True, timeout=180)
    for a in json.loads(cp.stdout):
        if a.get("description") == spec["app"] and not a.get("stopped_at"):
            r = subprocess.run(["modal", "app", "stop", "-y", a["app_id"]],
                               capture_output=True, text=True, timeout=180)
            print(f"  {a['app_id']}: {'stopped' if r.returncode == 0 else r.stderr[-200:]}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["launch", "status", "stop"])
    ap.add_argument("--engine", required=True, choices=sorted(ENGINES))
    ap.add_argument("--csv")
    ap.add_argument("--tag", default="")
    ap.add_argument("--samples", type=int, default=10)
    ap.add_argument("--seeds", default="1")
    ap.add_argument("--arms", default=None)
    a = ap.parse_args()
    if a.cmd == "launch":
        launch(a.engine, a.csv, a.tag, a.samples, a.seeds, a.arms)
    elif a.cmd == "status":
        status(a.engine, a.tag)
    else:
        stop(a.engine)
