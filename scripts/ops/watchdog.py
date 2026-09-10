"""Kill runaway remote jobs. Run me on a schedule; I am the guard of last resort.
The failure this exists for: a Modal app that keeps a GPU reserved on a wedged
subprocess, or an array that retries a poison input, quietly burning a month of credits
while nobody is watching. The PXR campaign's Modal runners had *no* timeout on their
inner `subprocess.run` and *no* container cap, so a hung child held an A100 for the full
8-hour function reservation with no signal that anything was wrong.
Checks, in order of how much they can save:
  1. **Modal apps running longer than their declared budget** -> stop them.
  2. **Ledger entries stuck in `launched`** past their expected finish -> flag (and stop
     the matching app if one is live).
  3. **Monthly spend approaching the cap** -> report loudly; refuse-to-launch is enforced
     in `budget.preflight_hours`, this is the early warning.
  4. **Local disk** -> the box has repeatedly hit 0 bytes free on C:; a full disk makes
     every downstream failure look like something else.
Default is dry-run. Pass --kill to actually stop apps; the guard should not be able to
terminate real work merely because a schedule fired.
    python scripts/ops/watchdog.py                 # report
    python scripts/ops/watchdog.py --kill          # report and stop runaways
    python scripts/ops/watchdog.py --max-age 4     # stricter age limit, in hours
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
from cypstruct import budget  # noqa: E402
from cypstruct.storage import health  # noqa: E402
# Apps this project owns. The watchdog must never touch an app it does not recognise —
# the user runs other Modal work from this account.
OWNED_APP_PREFIXES = ("cyp-cofold-", "cyp-qm-", "cyp-finetune-")
def modal_apps() -> list[dict]:
    """Live Modal apps, via the CLI's JSON output."""
    try:
        cp = subprocess.run(["modal", "app", "list", "--json"],
                            capture_output=True, text=True, timeout=180)
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        return [{"_error": f"could not list modal apps: {exc}"}]
    if cp.returncode != 0:
        return [{"_error": cp.stderr[-500:]}]
    try:
        return json.loads(cp.stdout)
    except json.JSONDecodeError:
        return [{"_error": "modal app list did not return JSON"}]
def _app_name(app: dict) -> str:
    """App name from `modal app list --json`, failing loudly on an unexpected schema."""
    if "description" in app:
        return str(app.get("description") or "")
    if "name" in app:
        return str(app.get("name") or "")
    raise KeyError(
        f"modal app record has neither 'description' nor 'name': {sorted(app)}. "
        "The CLI JSON schema changed; fix _app_name before trusting this watchdog.")
def _age_hours(created: str | float | None) -> float:
    if created is None:
        return 0.0
    if isinstance(created, (int, float)):
        return (time.time() - float(created)) / 3600
    from datetime import datetime, timezone
    txt = str(created).strip().replace("Z", "+00:00")
    try:
        dt = datetime.fromisoformat(txt)
    except ValueError:
        try:
            dt = datetime.strptime(txt[:19], "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
        except ValueError:
            # Unparseable timestamps must NOT read as age 0 - that silently exempts an
            # app from every age check. Report it as ancient so a human looks.
            return float("inf")
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return (datetime.now(timezone.utc) - dt).total_seconds() / 3600
def stop_app(app_id: str) -> str:
    # `modal app stop` prompts for confirmation and aborts without a TTY, so -y is not
    # optional here: without it the watchdog reports "stop failed" forever while the app
    # keeps running.
    cp = subprocess.run(["modal", "app", "stop", "-y", app_id],
                        capture_output=True, text=True, timeout=180)
    return "stopped" if cp.returncode == 0 else f"stop failed: {cp.stderr[-300:].strip()}"


CRASHLOOP_MARKERS = ("crash-looping", "containers are repeatedly failing to start",
                     "Runner failed with exception")


def crash_looping(app_id: str, tail_chars: int = 6000) -> str | None:
    """Detect a container that dies on import and is being restarted forever.

    This is the most expensive silent failure mode and the one an app-state check misses
    completely: a crash-looping app reports **zero tasks**, so it looks idle and harmless
    while Modal restarts it indefinitely. It is also easy to cause — a module-level path
    or import that only resolves on the developer's machine kills every container at
    import time. That exact bug hit this project on its first launch.

    Returns the offending log line, or None.
    """
    try:
        cp = subprocess.run(["modal", "app", "logs", app_id],
                            capture_output=True, text=True, timeout=90)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    tail = (cp.stdout or "")[-tail_chars:] + (cp.stderr or "")[-tail_chars:]
    for line in tail.splitlines():
        if any(m in line for m in CRASHLOOP_MARKERS):
            return line.strip()[:300]
    return None
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kill", action="store_true",
                    help="actually stop runaway apps (default is report only)")
    ap.add_argument("--max-age", type=float, default=8.0,
                    help="hours after which a running owned app is considered a runaway")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    report: dict = {"checked": time.strftime("%Y-%m-%d %H:%M:%S"),
                    "dry_run": not a.kill, "max_age_hours": a.max_age}
    # 1 + 2: live apps and stale ledger rows
    apps = modal_apps()
    err = next((x["_error"] for x in apps if "_error" in x), None)
    if err:
        report["modal_error"] = err
        apps = []
    # `modal app list --json` keys are: app_id, description, state, tasks, created_at,
    # stopped_at. An earlier version of this file read "Name"/"App ID"/"Created at" — the
    # column headings from the human-readable table — so every lookup returned None and the
    # watchdog reported zero runaways no matter what was running. A guard that cannot fail
    # loudly is worse than no guard, so `_app_name` raises on an unrecognised schema.
    owned = [x for x in apps if _app_name(x).startswith(OWNED_APP_PREFIXES)]
    runaways, idle_shells = [], []
    for x in owned:
        state = str(x.get("state", "")).lower()
        if x.get("stopped_at"):
            continue
        try:
            tasks = int(x.get("tasks") or 0)
        except (TypeError, ValueError):
            tasks = 0
        age = _age_hours(x.get("created_at"))
        rec = {"app_id": x.get("app_id"), "name": _app_name(x),
               "age_hours": round(age, 2), "state": state, "tasks": tasks}
        # An EPHEMERAL app with zero tasks is an empty client-side shell, not a runaway:
        # no containers means no GPU and no spend, and `modal app stop` refuses it anyway
        # (it returns "Aborted!"). Alarming on those produces a check that fails forever
        # and trains everyone to ignore the watchdog. Cost tracks TASKS, so act on tasks.
        if tasks == 0 and state == "ephemeral":
            # Zero tasks usually means an empty client shell. But a CRASH-LOOPING app also
            # reports zero tasks, so check the logs before writing it off as harmless.
            loop = crash_looping(rec["app_id"]) if rec["app_id"] else None
            if loop:
                rec["crash_loop"] = loop
                rec["action"] = (stop_app(rec["app_id"]) if a.kill
                                 else "would stop (dry run)")
                runaways.append(rec)
                continue
            if age > a.max_age:
                idle_shells.append(rec)
            continue
        if state in ("deployed", "running", "ephemeral") and age > a.max_age:
            if a.kill and rec["app_id"]:
                rec["action"] = stop_app(rec["app_id"])
                budget.close(rec["app_id"], "killed",
                             note=f"watchdog: {rec['age_hours']}h > {a.max_age}h")
            else:
                rec["action"] = "would stop (dry run)"
            runaways.append(rec)
    report["owned_apps_running"] = len(owned)
    report["runaways"] = runaways
    report["idle_ephemeral_shells"] = idle_shells   # informational; these cost nothing
    stale = budget.stale(max_age_hours=a.max_age)
    report["stale_ledger_rows"] = [
        {"run_id": r["run_id"], "kind": r["kind"], "venue": r["venue"],
         "age_hours": round((time.time() - r["started"]) / 3600, 2),
         "est_gpu_hours": r["est_gpu_hours"]}
        for r in stale]
    # 3: spend
    spend = {v: budget.spent(v) for v in ("modal", "boltz", "openprotein")}
    cap = budget.CAPS["modal_gpu_hours"]
    used = spend["modal"]["gpu_hours"]
    report["spend"] = spend
    report["modal_cap_fraction"] = round(used / cap, 3) if cap else None
    if cap and used > 0.8 * cap:
        report["ALERT_spend"] = (f"modal GPU hours at {used:.1f}/{cap:.0f} "
                                 f"({used/cap:.0%} of the monthly cap)")
    # 4: disk
    h = health()
    report["disk"] = h
    if h["free_C_gb"] < 1.0 or h["free_D_gb"] < 1.0:
        report["ALERT_disk"] = (f"low disk: C: {h['free_C_gb']} GB, D: {h['free_D_gb']} GB. "
                                "Bulk writes will fail; push artifacts to OneDrive.")
    out = REPO / "data" / "processed" / "watchdog_last.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    if a.json:
        print(json.dumps(report, indent=2))
    else:
        print(f"[watchdog {report['checked']}] "
              f"owned apps running: {report['owned_apps_running']}, "
              f"runaways: {len(runaways)}, stale ledger rows: {len(report['stale_ledger_rows'])}")
        print(f"  modal GPU-h this month: {used:.2f} / {cap:.0f}")
        print(f"  disk  C: {h['free_C_gb']} GB   D: {h['free_D_gb']} GB   "
              f"rclone cache {h['rclone_vfs_cache_gb']} GB")
        for k in ("ALERT_spend", "ALERT_disk", "modal_error"):
            if report.get(k):
                print(f"  !! {report[k]}")
        for r in runaways:
            print(f"  !! runaway {r['name']} ({r['age_hours']}h): {r['action']}")
    return 1 if (runaways or report.get("ALERT_spend") or report.get("ALERT_disk")) else 0
if __name__ == "__main__":
    raise SystemExit(main())

