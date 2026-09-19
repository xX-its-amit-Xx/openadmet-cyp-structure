"""Pull every MSA off OpenProtein and stage it for Explorer.

**This is the step that killed the PXR fine-tune.** Explorer's GPU compute nodes have no
direct internet (proxy `http://10.99.0.130:3128`), which is structurally incompatible with
any `--use_msa_server` flow: three successive debug jobs died on weight/MSA download before
anyone traced it to the network. The plan's own note is "budget MSA staging as its own
task", so that is what this is.

185 MSAs live as finished OpenProtein jobs. They are the single most expensive input we
own - each is a search against the full sequence database - and they exist only as remote
job results. Pulled to disk they become a portable asset that works offline, survives the
OpenProtein account, and can be rsynced to any cluster.

Resumable by construction: a target whose `.a3m` already exists on disk is skipped, so an
interrupted run costs nothing. That matters because this is 185 sequential network calls
and something will interrupt it.

    python scripts/finetune/stage_msas.py pull          # to data/reference/msa/
    python scripts/finetune/stage_msas.py push          # rsync to Explorer /scratch
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "cofold"))

from cypstruct.paths import DATA_PROCESSED, REFERENCE  # noqa: E402

UNI = DATA_PROCESSED / "p450_universe"
MSA_DIR = REFERENCE / "msa"
REMOTE = "explorer:/scratch/shenoy.am/cyp-finetune"


def cmd_pull(limit: int | None) -> dict:
    from openprotein_cofold import connect
    s = connect()
    st = json.loads((UNI / "campaign.json").read_text())
    MSA_DIR.mkdir(parents=True, exist_ok=True)
    ready = {k: v for k, v in st["msa"].items() if v.get("status") == "SUCCESS"}
    print(f"{len(ready)} MSA jobs to stage")

    done = fail = skip = 0
    for i, (target, rec) in enumerate(sorted(ready.items())):
        out = MSA_DIR / f"{target}.a3m"
        if out.exists() and out.stat().st_size > 0:
            skip += 1
            continue
        if limit and done >= limit:
            break
        try:
            job = s.load_job(rec["job_id"])
            rows = list(job.get())
            # an MSA comes back as rows of (id, sequence); write plain a3m so any tool
            # can read it, rather than a pickled object tied to this SDK version
            lines = []
            for r in rows:
                if isinstance(r, (list, tuple)) and len(r) >= 2:
                    lines.append(f">{r[0]}")
                    lines.append(str(r[1]))
                else:
                    lines.append(str(r))
            out.write_text("\n".join(lines))
            done += 1
        except Exception as exc:
            print(f"  {target}: {type(exc).__name__}: {str(exc)[:80]}", flush=True)
            fail += 1
        if (done + fail) and (done + fail) % 20 == 0:
            print(f"  pulled {done}, failed {fail}, skipped {skip}", flush=True)
    total = len(list(MSA_DIR.glob("*.a3m")))
    mb = sum(f.stat().st_size for f in MSA_DIR.glob("*.a3m")) / 1024**2
    return {"pulled": done, "failed": fail, "already_had": skip,
            "on_disk": total, "size_mb": round(mb, 1), "dir": str(MSA_DIR)}


def cmd_push() -> dict:
    """rsync to Explorer. Login node has internet; compute nodes do not - which is the
    whole reason this exists."""
    dest = REMOTE + "/msa/"
    subprocess.run(["ssh", "explorer", "mkdir -p /scratch/shenoy.am/cyp-finetune/msa"],
                   check=False, timeout=120)
    # rsync parses "D:\..." as host:path and refuses with "source and destination cannot
    # both be remote". Under Git Bash the drive letter has to become a POSIX root or the
    # colon is read as a hostname separator.
    src = str(MSA_DIR)
    if len(src) > 2 and src[1] == ":":
        src = "/" + src[0].lower() + src[2:].replace("\\", "/")
    cp = subprocess.run(
        ["rsync", "-az", "--info=stats2", f"{src}/", dest],
        capture_output=True, text=True, timeout=5400)
    ok = cp.returncode == 0
    return {"ok": ok, "dest": dest,
            "stderr": cp.stderr[-300:] if not ok else "",
            "summary": [l for l in cp.stdout.splitlines() if "files transferred" in l
                        or "Total file size" in l][:3]}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pull", "push"])
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    out = cmd_pull(a.limit or None) if a.cmd == "pull" else cmd_push()
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
