"""Paths + the storage policy, which on this box is a hard constraint, not advice.

Measured 2026-09-09 (`Get-PSDrive`):  **C: 2.1 GB free · D: 1.5 GB free · O: 4.94 TB free.**
The inherited CYP-challenge CLAUDE.md still says "C: ~16 GB" — that is stale by an
order of magnitude. A single 4-engine × 20-sample pose pool over a few hundred
ligands is tens of GB, so *nothing bulky may touch C: or D:*.

Rules enforced here:
  - POOL / QM / FINETUNE artifacts -> O: (cold, write-once, read-occasionally).
  - Repo (code + small JSON/parquet summaries) -> D:, and only that.
  - Genuinely latency-sensitive scratch -> C:/cyp_struct, but call `guard_scratch()`
    first; it raises rather than let a run wedge the box at 0 bytes free.
  - Anything heavy and hot belongs on a remote venue (Modal / Explorer), not here.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = _HERE.parents[2]

# --- repo-local (small files only) -----------------------------------------
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_EXTERNAL = PROJECT_ROOT / "data" / "external"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
SUBMISSIONS = PROJECT_ROOT / "submissions"
DOCS = PROJECT_ROOT / "docs"

# --- cold storage (OneDrive via rclone) ------------------------------------
_O = Path("O:/rclone-offload/cyp-structure")
OFFLOAD = _O if os.path.exists("O:/") else DATA_EXTERNAL / "_offload"

# NOTE: these are the *conceptual* remote locations. Do NOT write bulk data to them
# through the O: mount — its VFS cache is unbounded and lives on the near-full C:.
# Use `cypstruct.storage.push()/Batch()`, which goes to the rclone remote directly.
POOL = OFFLOAD / "pool"            # cofolded poses, per engine/target/ligand
QM_OUT = OFFLOAD / "qm"            # QM/NNP outputs (wavefunction-derived features)
FINETUNE = OFFLOAD / "finetune"    # training shards + checkpoints

# Experimental reference structures stay LOCAL: ~120 CYP3A4 mmCIFs is a few hundred MB,
# and every scorer-calibration pass re-reads them. Round-tripping that through a network
# mount on every run would be slow and would refill the C: cache each time.
REFERENCE = PROJECT_ROOT / "data" / "reference"

# --- hot scratch (tiny, guarded) -------------------------------------------
SCRATCH = Path("C:/cyp_struct") if os.name == "nt" else Path("/tmp/cyp_struct")
CLAIMS = SCRATCH / "claims"        # mkdir-based locks so parallel agents don't collide

# Refuse to start heavy local work below this much free space.
MIN_FREE_GB = 1.0


def free_gb(path: str | Path) -> float:
    """Free space in GB on the volume holding `path`."""
    try:
        return shutil.disk_usage(str(path)).free / 1024**3
    except OSError:
        return 0.0


def guard_scratch(need_gb: float = 1.0, path: Path | None = None) -> None:
    """Raise before writing if the target volume cannot spare `need_gb`.

    Called at the top of anything that writes more than a few MB locally. The
    failure mode this prevents is real: at 1.5 GB free, one careless pose dump
    fills D: and every other tool on the box starts erroring in confusing ways.
    """
    p = Path(path or SCRATCH)
    p.mkdir(parents=True, exist_ok=True)
    have = free_gb(p)
    if have < max(need_gb, MIN_FREE_GB):
        raise RuntimeError(
            f"refusing to write to {p}: {have:.2f} GB free, need {need_gb:.2f} GB. "
            f"Send this output to OFFLOAD ({OFFLOAD}) or run it on Modal/Explorer instead."
        )


def claim(slug: str) -> bool:
    """Atomically claim a work item. True if we got it, False if someone else has it.

    mkdir is atomic on NTFS, so this is a real lock across concurrent agents.
    """
    CLAIMS.mkdir(parents=True, exist_ok=True)
    try:
        (CLAIMS / slug).mkdir()
        return True
    except FileExistsError:
        return False


def ensure_dirs() -> None:
    for p in (DATA_RAW, DATA_EXTERNAL, DATA_PROCESSED, SUBMISSIONS,
              POOL, QM_OUT, FINETUNE, REFERENCE):
        p.mkdir(parents=True, exist_ok=True)


ensure_dirs()
