"""Bulk artifact storage that does not destroy this machine.

**The problem, measured 2026-09-09.** The `O:` drive is an rclone mount started with
`--vfs-cache-mode full --cache-dir C:\\Temp\\rclone-cache` and **no `--vfs-cache-max-size`**.
Every byte read *or written* through `O:\\` is therefore copied into a cache on `C:` and
kept indefinitely. Simply listing and reading a few OneDrive folders during orientation
grew that cache to **8 GB and drove C: to 0 bytes free**, which broke unrelated tools on
the box until it was cleaned up.

So the inherited advice "write cofold artifacts straight to O:" is actively harmful here:
a full pose pool written through the mount would fill C: several times over.

**What this module does instead.** Bulk artifacts are staged in a small directory on the
local disk and pushed to OneDrive with `rclone copy` against the **remote** (`onedrive:`),
which talks to the Graph API directly and never touches the VFS cache. Staging is deleted
after a verified transfer, so peak local usage is one batch, not one campaign.

Reads work the same way in reverse: `pull()` fetches just what is needed.

`O:\\` remains fine for *browsing* and for small files. It is not a bulk data path.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
import warnings
from pathlib import Path

from .paths import free_gb

RCLONE = shutil.which("rclone") or r"C:\Temp\tools\rclone.exe"
RCLONE_CONF = os.environ.get("RCLONE_CONFIG", r"C:\Temp\tools\rclone.conf")
REMOTE = "onedrive:rclone-offload/cyp-structure"

# Staging lives on the drive with the most room among the local disks. Small by design.
_STAGE_CANDIDATES = [Path("D:/cyp_stage"), Path("C:/cyp_stage"), Path("/tmp/cyp_stage")]
MAX_STAGE_GB = 2.0


def stage_dir() -> Path:
    """Pick a staging directory on whichever local volume currently has room."""
    best, best_free = None, -1.0
    for c in _STAGE_CANDIDATES:
        try:
            c.mkdir(parents=True, exist_ok=True)
        except OSError:
            continue
        f = free_gb(c)
        if f > best_free:
            best, best_free = c, f
    if best is None or best_free < 0.5:
        raise RuntimeError(
            f"no local volume has 0.5 GB free for staging (best {best_free:.2f} GB). "
            "Free space or run this step on Modal/Explorer."
        )
    return best


def _run(args: list[str], timeout: int = 3600) -> subprocess.CompletedProcess:
    return subprocess.run([RCLONE, *args, "--config", RCLONE_CONF],
                          capture_output=True, text=True, timeout=timeout)


def push(local: str | Path, remote_subpath: str, *, move: bool = True,
         transfers: int = 8, timeout: int = 3600) -> dict:
    """Copy (or move) a local directory to OneDrive, bypassing the O: VFS cache.

    Returns a small receipt dict. Raises on transfer failure — a silent partial upload
    followed by a local delete is the one outcome that actually loses data, so the
    delete only happens through `rclone move`, which verifies before removing.
    """
    local = Path(local)
    if not local.exists():
        raise FileNotFoundError(local)
    # REMOTE already ends in the project name, so passing "cyp-structure/pool/x" nests it
    # twice and the archive lands somewhere nobody looking beside pool/op1 would ever find
    # it. That happened to a 1.7 GiB archive and had to be moved server-side afterwards.
    sub = remote_subpath.strip("/")
    leaf = REMOTE.rstrip("/").rsplit("/", 1)[-1]
    if sub == leaf or sub.startswith(leaf + "/"):
        sub = sub[len(leaf):].lstrip("/")
        warnings.warn(
            f"remote_subpath started with {leaf!r}, which REMOTE already ends in; "
            f"using {sub!r} so the archive lands beside the others, not nested under "
            f"a second {leaf!r}", stacklevel=2)
    dest = f"{REMOTE}/{sub}"
    t0 = time.time()
    verb = "move" if move else "copy"
    cp = _run([verb, str(local), dest, "--transfers", str(transfers),
               "--stats-one-line", "--stats", "30s", "-v"], timeout=timeout)
    if cp.returncode != 0:
        raise RuntimeError(f"rclone {verb} failed ({cp.returncode}): {cp.stderr[-2000:]}")
    return {"local": str(local), "remote": dest, "verb": verb,
            "seconds": round(time.time() - t0, 1)}


def pull(remote_subpath: str, local: str | Path, *, timeout: int = 3600) -> Path:
    """Fetch a remote subtree to a local path (again bypassing the mount)."""
    local = Path(local)
    local.mkdir(parents=True, exist_ok=True)
    src = f"{REMOTE}/{remote_subpath.strip('/')}"
    cp = _run(["copy", src, str(local), "--transfers", "8"], timeout=timeout)
    if cp.returncode != 0:
        raise RuntimeError(f"rclone copy failed ({cp.returncode}): {cp.stderr[-2000:]}")
    return local


def ls(remote_subpath: str = "") -> list[str]:
    """List a remote subtree. Cheap way to implement skip-done without downloading."""
    src = f"{REMOTE}/{remote_subpath.strip('/')}" if remote_subpath else REMOTE
    cp = _run(["lsf", "-R", src], timeout=600)
    if cp.returncode != 0:
        return []
    return [ln for ln in cp.stdout.splitlines() if ln.strip()]


class Batch:
    """Stage a batch of artifacts locally, then push and clear.

    Usage:
        with Batch("pool/boltz2/cyp3a4") as b:
            (b.path / "lig001.cif").write_text(cif)
            ...
        # on exit: rclone move to OneDrive, staging removed

    The context manager exists so the "delete the staging" step cannot be forgotten,
    which is how a staging directory quietly becomes a second copy of the campaign.
    """

    def __init__(self, remote_subpath: str, max_gb: float = MAX_STAGE_GB):
        self.remote_subpath = remote_subpath
        self.max_gb = max_gb
        self.path = stage_dir() / remote_subpath.replace("/", "_")
        self.receipt: dict | None = None

    def __enter__(self) -> "Batch":
        self.path.mkdir(parents=True, exist_ok=True)
        return self

    def size_gb(self) -> float:
        return sum(f.stat().st_size for f in self.path.rglob("*") if f.is_file()) / 1024**3

    def check(self) -> None:
        """Call periodically inside a long loop; raises before the disk fills."""
        if self.size_gb() > self.max_gb:
            raise RuntimeError(
                f"staging {self.path} exceeded {self.max_gb} GB. Flush more often: "
                "push a batch and start a new one."
            )
        if free_gb(self.path) < 0.5:
            raise RuntimeError(f"only {free_gb(self.path):.2f} GB free at {self.path}")

    def __exit__(self, exc_type, exc, tb) -> bool:
        if exc_type is not None:
            # leave staging in place so the partial work is recoverable
            return False
        if any(self.path.rglob("*")):
            self.receipt = push(self.path, self.remote_subpath, move=True)
            (self.path.parent / f"{self.path.name}.receipt.json").write_text(
                json.dumps(self.receipt, indent=1))
        shutil.rmtree(self.path, ignore_errors=True)
        return False


def health() -> dict:
    """Disk + cachestatus. Print this before any bulk step."""
    cache = Path("C:/Temp/rclone-cache/vfs")
    cache_gb = 0.0
    if cache.exists():
        try:
            cache_gb = sum(f.stat().st_size for f in cache.rglob("*") if f.is_file()) / 1024**3
        except OSError:
            cache_gb = float("nan")
    return {
        "free_C_gb": round(free_gb("C:/"), 2),
        "free_D_gb": round(free_gb("D:/"), 2),
        "rclone_vfs_cache_gb": round(cache_gb, 2),
        "warning": ("rclone VFS cache is unbounded (no --vfs-cache-max-size); "
                    "never write bulk data through O:\\ — use storage.push()"),
    }


if __name__ == "__main__":
    print(json.dumps(health(), indent=2))
