"""Spend and runaway guards for remote compute.

This exists because the specific failure the user named is real and cheap to hit: a
Modal app that retries a failing input forever, or a poll loop with no exit, burns a
month of credits overnight with nothing to show for it.

Three independent guards, because any one of them can be defeated:

1. **Per-function limits, declared at the decorator.** `timeout`, `retries`,
   `max_containers` are set explicitly on every remote function in this repo.
   Modal's defaults (no container cap) are the wrong defaults for a credit budget.
2. **A local ledger.** Every launch records estimated GPU-seconds before it starts.
   `preflight()` refuses a launch that would push the running total past the cap.
   Refusing before spending is the only check that helps.
3. **A watchdog.** `scripts/ops/watchdog.py` lists live Modal apps and kills anything
   older than its declared budget, so a wedged run dies without needing me awake.

The ledger is deliberately a plain JSON file, not a service: it has to be readable and
editable by a human at 3am when something has gone wrong.
"""
from __future__ import annotations

import json
import os
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .paths import DATA_PROCESSED

LEDGER = DATA_PROCESSED / "compute_ledger.json"

# Monthly ceilings. Modal credits refresh monthly; these are set well under the grant
# so that a bug costs a fraction of the month rather than all of it.
CAPS = {
    "modal_gpu_hours": 60.0,      # across all GPU functions, per calendar month
    "modal_cpu_hours": 200.0,
    "boltz_api_jobs": 4000,       # hosted Boltz predictions per month
    "openprotein_jobs": 2000,
}

# Cost model. A co-folding *job* runs the trunk once and then diffuses N samples, so cost
# is NOT linear in samples — it is (fixed trunk) + (small per-sample increment). Charging
# per sample overestimates a 20-sample job roughly fourfold, which matters because the
# preflight then refuses runs that would actually fit. Estimates are seeded from the PXR
# campaign's observed runtimes and are overwritten by `calibrate_from_ledger()` once this
# project has measured its own.
COST_HINTS_JOB = {           # GPU-hours of fixed cost per job (trunk + MSA load + I/O)
    "boltz2_cofold": 0.030,
    "af3_cofold": 0.060,
    "protenix_cofold": 0.045,
    "chai_cofold": 0.035,
}
COST_HINTS_SAMPLE = {        # additional GPU-hours per diffusion sample
    "boltz2_cofold": 0.004,
    "af3_cofold": 0.008,
    "protenix_cofold": 0.006,
    "chai_cofold": 0.005,
}
# Non-cofold work, still charged per unit.
COST_HINTS = {
    "xtb_pocket_singlepoint": 0.002,   # CPU-hour
    "finetune_step": 0.004,
}


def estimate(kind: str, n_jobs: int, samples_per_job: int = 1) -> float:
    """GPU-hours for `n_jobs` jobs of `kind`, each producing `samples_per_job` samples."""
    if kind in COST_HINTS_JOB:
        return n_jobs * (COST_HINTS_JOB[kind]
                         + samples_per_job * COST_HINTS_SAMPLE.get(kind, 0.004))
    return n_jobs * samples_per_job * COST_HINTS.get(kind, 0.02)


def calibrate_from_ledger(kind: str) -> dict:
    """Replace the cost hints for `kind` with what this project actually measured.

    Guessed cost models drift; a preflight built on a stale guess either blocks work
    that fits or waves through work that does not.
    """
    rows = [r for r in _load()
            if r.get("kind") == kind and r.get("actual_gpu_hours") and r.get("meta", {}).get("n_jobs")]
    if len(rows) < 2:
        return {"kind": kind, "status": "not enough measured runs", "n": len(rows)}
    import statistics
    per_job = [r["actual_gpu_hours"] / r["meta"]["n_jobs"] for r in rows]
    med = statistics.median(per_job)
    COST_HINTS_JOB[kind] = round(med * 0.7, 4)
    COST_HINTS_SAMPLE[kind] = round(med * 0.3 / max(1, rows[-1]["meta"].get("samples", 5)), 5)
    return {"kind": kind, "status": "calibrated", "n": len(rows),
            "median_gpu_h_per_job": round(med, 4)}


# Minimum wall-clock time before a run may be judged "not advancing", per job kind.
# A co-folding job that diffuses 20 samples takes minutes; concluding it has stalled
# because no output exists after two of them is how a healthy run gets killed. That
# happened once: the ops tick stopped a 168-job batch 2.5 minutes after launch, read
# "no output directories yet" as a stall, and destroyed a run that was fine.
MIN_AGE_BEFORE_STALL_CALL = {
    "boltz2_cofold": 900,      # seconds; ~5 min per job plus image pull and queueing
    "af3_cofold": 1800,
    "protenix_cofold": 1200,
    "chai_cofold": 1200,
    "_default": 900,
}


def may_judge_stalled(kind: str, started: float) -> tuple[bool, str]:
    """Is this run old enough that "no progress" actually means something?

    Returns (may_judge, reason_if_not). Call this BEFORE killing anything for lack of
    output. Liveness checks that fire faster than one unit of work can complete do not
    detect stalls, they cause them.
    """
    need = MIN_AGE_BEFORE_STALL_CALL.get(kind, MIN_AGE_BEFORE_STALL_CALL["_default"])
    age = time.time() - started
    if age < need:
        return False, (f"only {age/60:.1f} min old; a {kind} job needs about "
                       f"{need/60:.0f} min before absence of output means anything")
    return True, ""


@dataclass
class Entry:
    run_id: str
    venue: str                 # modal | boltz | openprotein | explorer
    kind: str                  # a key of COST_HINTS, or free text
    units: int
    est_gpu_hours: float
    started: float
    month: str
    status: str = "launched"   # launched | done | killed | failed
    actual_gpu_hours: float | None = None
    note: str = ""
    meta: dict = field(default_factory=dict)


def _load() -> list[dict]:
    if LEDGER.exists():
        try:
            return json.loads(LEDGER.read_text())
        except json.JSONDecodeError:
            bad = LEDGER.with_suffix(".corrupt.json")
            LEDGER.rename(bad)
            return []
    return []


def _save(rows: list[dict]) -> None:
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    tmp = LEDGER.with_suffix(".tmp")
    tmp.write_text(json.dumps(rows, indent=1))
    tmp.replace(LEDGER)


def month_key(ts: float | None = None) -> str:
    return time.strftime("%Y-%m", time.localtime(ts or time.time()))


def spent(venue: str = "modal", month: str | None = None) -> dict:
    """GPU/CPU hours booked this month, counting launched-but-unfinished runs too.

    Counting in-flight work is the point: a check that only sums *completed* runs
    happily green-lights ten more launches while ten are already burning.
    """
    m = month or month_key()
    rows = [r for r in _load() if r.get("month") == m and r.get("venue") == venue]
    est = sum(r.get("actual_gpu_hours") or r.get("est_gpu_hours", 0.0)
              for r in rows if r.get("status") in ("launched", "done"))
    return {"month": m, "venue": venue, "n_runs": len(rows), "gpu_hours": round(est, 3),
            "in_flight": sum(1 for r in rows if r.get("status") == "launched")}


def preflight_hours(kind: str, est: float, venue: str = "modal",
                    cap_key: str = "modal_gpu_hours") -> tuple[bool, str, float]:
    """Same gate as `preflight`, but taking an already-computed GPU-hour estimate."""
    used = spent(venue)["gpu_hours"]
    cap = CAPS.get(cap_key, 1e9)
    if used + est > cap:
        return False, (f"would exceed {cap_key}: {used:.2f} used + {est:.2f} estimated "
                       f"> {cap:.2f} cap."), est
    if est > cap * 0.5:
        return False, (f"single launch estimated at {est:.2f} GPU-h is over half the "
                       f"{cap:.2f} monthly cap - split it into batches."), est
    return True, "", est


def preflight(kind: str, units: int, venue: str = "modal",
              cap_key: str = "modal_gpu_hours") -> tuple[bool, str, float]:
    """Would this launch fit inside the monthly cap? Returns (ok, reason, est_hours)."""
    per = COST_HINTS.get(kind, 0.02)
    est = per * units
    used = spent(venue)["gpu_hours"]
    cap = CAPS.get(cap_key, 1e9)
    if used + est > cap:
        return False, (f"would exceed {cap_key}: {used:.2f} used + {est:.2f} estimated "
                       f"> {cap:.2f} cap. Reduce units or raise the cap deliberately."), est
    if est > cap * 0.5:
        return False, (f"single launch estimated at {est:.2f} GPU-h is over half the "
                       f"{cap:.2f} monthly cap — split it into batches."), est
    return True, "", est


def record(run_id: str, venue: str, kind: str, units: int, est: float,
           note: str = "", **meta) -> Entry:
    e = Entry(run_id=run_id, venue=venue, kind=kind, units=units, est_gpu_hours=est,
              started=time.time(), month=month_key(), note=note, meta=meta)
    rows = _load()
    rows.append(asdict(e))
    _save(rows)
    return e


def close(run_id: str, status: str = "done", actual_gpu_hours: float | None = None,
          note: str = "") -> None:
    rows = _load()
    for r in rows:
        if r.get("run_id") == run_id:
            r["status"] = status
            if actual_gpu_hours is not None:
                r["actual_gpu_hours"] = actual_gpu_hours
            if note:
                r["note"] = (r.get("note", "") + " | " + note).strip(" |")
            r["ended"] = time.time()
    _save(rows)


def stale(max_age_hours: float = 6.0) -> list[dict]:
    """Runs still marked 'launched' well past when they should have finished.

    This is what the watchdog acts on. A run that has been 'launched' for eight hours
    is either wedged or looping; either way nobody is reading its output.
    """
    now = time.time()
    return [r for r in _load()
            if r.get("status") == "launched"
            and (now - r.get("started", now)) > max_age_hours * 3600]


def summary() -> dict:
    out = {"caps": CAPS, "by_venue": {}}
    for v in ("modal", "boltz", "openprotein", "explorer"):
        out["by_venue"][v] = spent(v)
    out["stale_runs"] = [{"run_id": r["run_id"], "kind": r["kind"],
                          "age_hours": round((time.time() - r["started"]) / 3600, 2)}
                         for r in stale()]
    return out


if __name__ == "__main__":
    print(json.dumps(summary(), indent=2))

# --------------------------------------------------------------------------
# ACTUAL spend, in dollars, from Modal itself
#
# The guard in this module failed on 2026-09-11 and it is worth being precise about how.
# It tracked GPU-HOURS against a cap I invented (60/month), computed from MY OWN estimates
# of what each run would cost. Modal bills in DOLLARS against a workspace spend limit I
# never queried. So the ledger read "44 of 60 GPU-h, headroom available" at the exact
# moment Modal refused a launch with "Workspace has exceeded its spend limit".
#
# Worse, the ledger counted runs that were planned and then killed, so it over-reported
# the very number it was gating on while under-reporting the one that mattered.
#
# A budget guard that measures a proxy instead of the billed quantity is not a guard.
# --------------------------------------------------------------------------

MONTHLY_USD_CAP = 25.0     # act well below the real limit; raise deliberately, not by drift


def modal_actual_spend(period: str = "this month") -> dict:
    """Real dollars billed by Modal, from `modal billing report`.

    This is the authoritative number. The local ledger is a planning aid; this is what
    actually stops the work.
    """
    import subprocess

    try:
        cp = subprocess.run(["modal", "billing", "report", "--for", period, "--json"],
                            capture_output=True, timeout=420,
                            env={**os.environ, "PYTHONIOENCODING": "utf-8",
                                 "PYTHONUTF8": "1"})
        raw = (cp.stdout or b"").decode("utf-8", "replace")
        start = raw.find("[")
        rows = json.loads(raw[start:]) if start >= 0 else []
    except Exception as exc:
        return {"ok": False, "error": f"{type(exc).__name__}: {exc}"}

    by_app: dict[str, float] = {}
    total = 0.0
    for r in rows:
        try:
            c = float(r.get("cost", 0.0))
        except (TypeError, ValueError):
            continue
        total += c
        by_app[r.get("description", "?")] = by_app.get(r.get("description", "?"), 0.0) + c
    return {"ok": True, "period": period, "total_usd": round(total, 2),
            "by_app": {k: round(v, 2) for k, v in sorted(by_app.items(),
                                                         key=lambda kv: -kv[1])},
            "cap_usd": MONTHLY_USD_CAP,
            "over_cap": total > MONTHLY_USD_CAP}


def preflight_usd(est_gpu_hours: float, usd_per_gpu_hour: float = 2.2) -> tuple[bool, str]:
    """Gate a launch on REAL dollars remaining, not on estimated hours.

    `usd_per_gpu_hour` is calibrated from the val87b batch: 168 jobs x 20 samples billed
    **$26.24**, against an estimate of 18.5 GPU-h, i.e. about $1.4 per estimated GPU-hour
    on A100-40GB. The default is deliberately above that so the gate errs toward refusing.
    """
    spend = modal_actual_spend()
    if not spend.get("ok"):
        return False, f"cannot read Modal billing ({spend.get('error')}); refusing to launch blind"
    projected = spend["total_usd"] + est_gpu_hours * usd_per_gpu_hour
    if projected > MONTHLY_USD_CAP:
        return False, (f"${spend['total_usd']:.2f} already billed this month; this launch "
                       f"adds about ${est_gpu_hours * usd_per_gpu_hour:.2f}, projecting "
                       f"${projected:.2f} against a ${MONTHLY_USD_CAP:.2f} cap")
    return True, ""

# --------------------------------------------------------------------------
# Workspaces
#
# Modal credit lives per WORKSPACE, and so do Volumes. Switching workspaces to get more
# credit therefore also loses read access to `cyp-pool`, the weight caches and the staged
# MSA - they must be re-created in the new workspace. Archive first, switch second.
#
# Known workspaces (2026-09-11, from Amit):
#   ashenoydiscovery  - primary; OVER its spend limit at $37.58 billed. Holds every
#                       existing volume: cyp-pool, cyp-boltz-cache, cyp-chai-cache.
#   xx-its-amit-xx    - the GitHub account, ~$30 credit, $0.00 billed. Clean slate.
#   re-agent-hackathon - ~$140 credit. NOT configured on this box; needs its token added
#                       to ~/.modal.toml. Amit's guidance: reserve it for FINE-TUNING /
#                       pretraining, not for inference batches.
# --------------------------------------------------------------------------

WORKSPACES = {
    "ashenoydiscovery": {"credit_usd": 0.0, "note": "primary; over spend limit; holds all volumes"},
    "xx-its-amit-xx": {"credit_usd": 30.0, "note": "GitHub account; use for inference batches"},
    "re-agent-hackathon": {"credit_usd": 140.0, "note": "NOT CONFIGURED; reserve for fine-tuning"},
}


def active_profile() -> str:
    """The Modal profile a subprocess would use right now."""
    env = os.environ.get("MODAL_PROFILE")
    if env:
        return env
    try:
        import tomllib
        cfg = tomllib.loads(Path.home().joinpath(".modal.toml").read_text())
        for name, block in cfg.items():
            if isinstance(block, dict) and block.get("active"):
                return name
    except Exception:
        pass
    return "?"


def workspace_env(profile: str) -> dict:
    """Environment for running a Modal command against a specific workspace.

    Set per-command rather than flipping the global active profile: other Claude sessions
    share this machine, and silently repointing their Modal client would be rude at best.
    """
    return {**os.environ, "MODAL_PROFILE": profile,
            "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}

