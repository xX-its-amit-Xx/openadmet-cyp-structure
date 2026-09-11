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


def _deploy(module: str, profile: str | None = None) -> str:
    """`modal deploy` the module so its app becomes server-side and permanent."""
    path = REPO / "scripts" / "cofold" / f"{module}.py"
    # `modal deploy` prints a U+2713 check mark. On this Windows box the child inherits a
    # cp1252 stdout and dies encoding it, so the deploy "fails" for reasons that have
    # nothing to do with the deploy. Force UTF-8 on the child and decode leniently here.
    import os as _os

    env = {**_os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}
    if profile:
        # Set the profile per-command rather than flipping the global active one: other
        # Claude sessions share this machine and repointing their Modal client silently
        # would be an unpleasant surprise.
        env["MODAL_PROFILE"] = profile
    cp = subprocess.run(["modal", "deploy", str(path)],
                        capture_output=True, timeout=1800, env=env)
    out = (cp.stdout or b"").decode("utf-8", "replace")
    err = (cp.stderr or b"").decode("utf-8", "replace")
    cp = subprocess.CompletedProcess(cp.args, cp.returncode, out, err)
    if cp.returncode != 0:
        raise RuntimeError(f"modal deploy failed:\n{cp.stdout[-1500:]}\n{cp.stderr[-2500:]}")
    return (cp.stdout or "")[-800:]


def _done_count(tag: str, attempts: int = 3) -> tuple[int, int]:
    """(#jobs with DONE.json, #job dirs) from ONE recursive listing.

    The previous version issued one `iterdir` per job and counted any exception as
    "not done". Modal's volume reads intermittently return 500s, so the number wobbled
    badly - it read 28, then 37, then 25 within a few minutes on a run that was only ever
    moving forward. Every progress judgement built on it was therefore unreliable, and it
    twice made an advancing batch look stalled.

    One recursive listing, retried as a whole, is both far cheaper and far steadier: a
    transient failure now retries instead of silently subtracting from the count.
    """
    import time as _t

    vol = modal.Volume.from_name("cyp-pool")
    last = None
    for k in range(attempts):
        try:
            entries = list(vol.iterdir(f"/{tag}", recursive=True))
            jobs, done = set(), set()
            for e in entries:
                parts = e.path.strip("/").split("/")
                if len(parts) < 2 or "__" not in parts[1]:
                    continue
                jobs.add(parts[1])
                if parts[-1] == "DONE.json":
                    done.add(parts[1])
            return len(done), len(jobs)
        except Exception as exc:
            last = exc
            _t.sleep(2 * (k + 1))
    print(f"  [warn] could not list /{tag}: {type(last).__name__}", flush=True)
    return -1, -1          # sentinel: UNKNOWN, never mistake it for zero progress


def launch(engine: str, csv: str, tag: str, samples: int, seeds: str,
           arms: str | None, profile: str | None = None) -> None:
    import os as _os

    from cypstruct import budget

    if profile:
        # Credit AND volumes are per-workspace. Switching for credit means the new
        # workspace has no cyp-pool, no weight cache and no staged MSA - they are rebuilt
        # on first use. Archive before switching; see budget.WORKSPACES.
        _os.environ["MODAL_PROFILE"] = profile
        print(f"workspace: {profile} (volumes and credit are per-workspace)", flush=True)

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
    print(f"{len(jobs)} jobs planned ({done_before} already have DONE.json)  "
          f"est {est:.2f} GPU-h", flush=True)

    # Gate on REAL DOLLARS in the workspace we are about to spend them in. The
    # GPU-hour gate read "headroom available" at the exact moment Modal refused a launch
    # for exceeding its spend limit, because it was measuring a proxy I invented rather
    # than the quantity Modal bills. Hours are still printed, but dollars decide.
    spend = budget.modal_actual_spend()
    if spend.get("ok"):
        print(f"  workspace billed this month: ${spend['total_usd']:.2f}", flush=True)
    ok, why = budget.preflight_usd(est)
    if not ok:
        raise SystemExit(f"PREFLIGHT REFUSED (dollars): {why}")
    ok_h, why_h, _ = budget.preflight_hours(spec["kind"], est)
    if not ok_h:
        print(f"  note: local hour ledger also objects ({why_h})", flush=True)

    print("deploying the app so the run outlives this process ...", flush=True)
    print("  ", _deploy(spec["module"], profile).strip().splitlines()[-1:], flush=True)

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
    import os as _os

    env = {**_os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}
    try:
        cp = subprocess.run(["modal", "app", "logs", app_id],
                            capture_output=True, timeout=90, env=env)
    except Exception:
        return None
    out = (cp.stdout or b"").decode("utf-8", "replace")
    err = (cp.stderr or b"").decode("utf-8", "replace")
    tail = (out + err)[-6000:]
    for line in reversed(tail.splitlines()):
        if any(m in line for m in CAPACITY_MARKERS):
            return line.strip()[:200]
    return None


def status(engine: str, tag: str, profile: str | None = None) -> None:
    import os as _os

    if profile:
        _os.environ["MODAL_PROFILE"] = profile
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
                if w:
                    note = f"  <- QUEUED FOR CAPACITY: {w}"
                else:
                    # tasks=0 is a SNAPSHOT. Containers scale down between batches, so a
                    # single reading of zero says nothing - it has now looked like a stall
                    # three times on runs that were advancing fine. Re-check output
                    # progress over an interval before using the word idle.
                    d0, _ = _done_count(tag)
                    time.sleep(45)
                    d1, _ = _done_count(tag)
                    note = (f"  <- ADVANCING (+{d1 - d0} jobs in 45 s)" if d1 > d0
                            else "  <- idle (no container, no progress in 45 s)")
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
    ap.add_argument("--profile", default=None,
                    help="Modal workspace to run in (see cypstruct.budget.WORKSPACES). "
                         "Volumes are per-workspace, so a switch starts from empty caches.")
    a = ap.parse_args()
    if a.cmd == "launch":
        launch(a.engine, a.csv, a.tag, a.samples, a.seeds, a.arms, a.profile)
    elif a.cmd == "status":
        status(a.engine, a.tag, a.profile)
    else:
        stop(a.engine)
