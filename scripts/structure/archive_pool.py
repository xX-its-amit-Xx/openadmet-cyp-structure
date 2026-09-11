r"""Archive a co-folding pose pool from the Modal volume to OneDrive, durably.

**Why this exists.** Right now the only copy of a pose pool lives on a Modal volume.
That is convenient but not durable: it is tied to one account's storage, it is not
browsable, and it is not where the rest of this user's finished work lives. The pool is
also expensive — the val87b batch cost real GPU hours — so it should survive.

**Why it does NOT write through `O:\`.** The `O:` drive is an rclone mount configured
with `--vfs-cache-mode full` and **no `--vfs-cache-max-size`**, with its cache on the
near-full `C:`. Every byte written through the drive letter is copied to C: and kept.
That mount is what drove C: to zero bytes three separate times in this project, once
silently truncating 40% of a scoring run. So this script streams to the SAME OneDrive
storage through the `onedrive:` rclone remote, which talks to the Graph API directly and
never touches the VFS cache.

Shape on OneDrive:
    onedrive:rclone-offload/cyp-structure/pool/<tag>/<ligand>__<arm>__s<seed>/*.cif

Resumable: a job already present remotely is skipped, so an interrupted archive resumes
rather than restarts.

    python scripts/structure/archive_pool.py --tag val87b
    python scripts/structure/archive_pool.py --tag val87b --include-npz   # also the arrays
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

import modal  # noqa: E402

from cypstruct import storage  # noqa: E402
from cypstruct.paths import DATA_PROCESSED, free_gb, safe_workers  # noqa: E402


def scratch_root() -> Path:
    """Staging directory on a local volume with room. Never the session temp dir."""
    d = Path("D:/cyp_scratch")
    try:
        d.mkdir(parents=True, exist_ok=True)
        if free_gb(d) >= 5.0:
            return d
    except OSError:
        pass
    c = Path("C:/cyp_scratch")
    c.mkdir(parents=True, exist_ok=True)
    if free_gb(c) < 2.0:
        raise RuntimeError("no local volume has room to stage an archive batch")
    return c


def read_with_retry(vol, path: str, attempts: int = 4) -> bytes | None:
    """Read one volume file, tolerating transient server errors.

    Modal's volume reads intermittently return `500 Internal Server Error: error while
    getting bucket object`. The first archive attempt aborted the entire run on the first
    such error and transferred nothing. A transient remote fault should cost one file and
    a retry, never a whole batch, so failures here are isolated and reported rather than
    raised.
    """
    import time as _t

    for k in range(attempts):
        try:
            return b"".join(vol.read_file(path))
        except Exception as exc:
            if k == attempts - 1:
                print(f"    [skip] {path}: {type(exc).__name__}", flush=True)
                return None
            _t.sleep(1.5 * (k + 1))
    return None



def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--include-npz", action="store_true",
                    help="also archive the per-sample plddt/pae/pde arrays (much larger)")
    ap.add_argument("--batch", type=int, default=24,
                    help="jobs to stage before each push; keeps peak local usage small")
    a = ap.parse_args()

    vol = modal.Volume.from_name("cyp-pool")
    jobs = sorted({e.path.rstrip("/").split("/")[-1] for e in vol.iterdir(f"/{a.tag}")})
    jobs = [j for j in jobs if "__" in j]
    print(f"{len(jobs)} jobs in /{a.tag}", flush=True)

    already = {ln.split("/")[0] for ln in storage.ls(f"pool/{a.tag}") if "/" in ln}
    todo = [j for j in jobs if j not in already]
    print(f"{len(already)} already archived, {len(todo)} to go", flush=True)
    if not todo:
        print("nothing to do")
        return

    root = scratch_root()
    pushed = n_files = 0
    failed: list[str] = []
    t0 = time.time()

    for i in range(0, len(todo), a.batch):
        chunk = todo[i: i + a.batch]
        stage = Path(tempfile.mkdtemp(prefix=f"cyparch_{a.tag}_", dir=str(root)))
        try:
            # Volume reads are network-bound, so fetch jobs concurrently. Sequential
            # archiving ran at ~12 jobs per 260 s, which is far too slow to finish inside
            # one turn - and this data is the single most expensive artifact in the
            # project, so getting a durable copy quickly matters.
            def fetch(job: str) -> tuple[str, int, list[str]]:
                dest = stage / job
                dest.mkdir(parents=True, exist_ok=True)
                got, bad = 0, []
                try:
                    entries = list(vol.iterdir(f"/{a.tag}/{job}"))
                except Exception as exc:
                    return job, 0, [f"{job} <iterdir {type(exc).__name__}>"]
                for e in entries:
                    fn = e.path.split("/")[-1]
                    keep = (fn.endswith(".cif") or fn.endswith(".json")
                            or (a.include_npz and fn.endswith(".npz")))
                    if not keep:
                        continue
                    data = read_with_retry(vol, e.path)
                    if data is None:
                        bad.append(f"{job}/{fn}")
                        continue
                    (dest / fn).write_bytes(data)
                    got += 1
                return job, got, bad

            nw = safe_workers(8, ram_per_worker_gb=0.5)
            with ThreadPoolExecutor(max_workers=nw) as ex:
                for _job, got, bad in ex.map(fetch, chunk):
                    n_files += got
                    failed.extend(bad)
            if free_gb(stage) < 2.0:
                raise RuntimeError(f"only {free_gb(stage):.2f} GB free while staging")
            # rclone MOVE: it verifies the transfer before deleting the local copy, so a
            # partial upload cannot silently destroy the staged data.
            storage.push(stage, f"pool/{a.tag}", move=True)
            pushed += len(chunk)
            print(f"  archived {pushed}/{len(todo)} jobs, {n_files} files, "
                  f"{time.time()-t0:.0f}s", flush=True)
        finally:
            shutil.rmtree(stage, ignore_errors=True)

    receipt = {"tag": a.tag, "jobs_archived": pushed, "files": n_files,
               "n_failed_reads": len(failed), "failed_sample": failed[:15],
               "remote": f"{storage.REMOTE}/pool/{a.tag}",
               "include_npz": a.include_npz,
               "seconds": round(time.time() - t0, 1)}
    (DATA_PROCESSED / f"archive_{a.tag}.json").write_text(json.dumps(receipt, indent=1))
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
