"""AlphaFold 3 co-folding on Modal, weights fetched from Kaggle.

Ported from the PXR campaign's runner with four defects fixed. Those defects are worth
naming because each one is a way to spend GPU hours on nothing:

1. **The inner `subprocess.run` had no `timeout=`.** A hung AF3 child held an A100 for the
   full 8-hour function reservation, silently. Fixed: an explicit inner timeout, always
   strictly less than the Modal function timeout so the container reports rather than dies.
2. **No `retries=`, no container cap.** Fixed: `retries=1`, `max_containers` bounded.
3. **Exceptions were swallowed per ligand**, so 184 consecutive failures still ran to
   completion and paid full wall-clock. Fixed: a consecutive-failure circuit breaker.
4. **Two entry points referenced an undefined variable** and raised `NameError` at call
   time. Not ported.

**Licence.** AF3 model parameters are provided by Google DeepMind under terms that permit
non-commercial use and **prohibit redistribution**. They are pulled at runtime from the
user's own Kaggle account into a private Modal volume and never leave it. Fine-tuning the
weights is outside those terms — see docs/FINETUNE_PLAN.md.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import modal

# Modal serialises this module and imports it INSIDE the container, where the file lands
# at /root/modal_boltz.py — so `parents[2]` raises IndexError and every container dies on
# import. Modal then restarts them, which is precisely the crash-loop that burns credits
# while looking like "still starting". Resolve the repo path only when it exists.
_here = Path(__file__).resolve()
REPO = _here.parents[2] if len(_here.parents) > 2 else None
if REPO is not None and (REPO / "src" / "cypstruct").is_dir():
    sys.path.insert(0, str(REPO / "src"))

APP_NAME = "cyp-cofold-af3"
GPU = os.environ.get("CYP_AF3_GPU", "A100-40GB")
MAX_CONTAINERS = int(os.environ.get("CYP_AF3_MAX_CONTAINERS", "4"))
FN_TIMEOUT = int(os.environ.get("CYP_AF3_TIMEOUT", "10800"))
INNER_TIMEOUT = FN_TIMEOUT - 600          # always leave the container room to report
MAX_CONSECUTIVE_FAILURES = 5

image = (
    modal.Image.from_registry("nvidia/cuda:12.6.3-devel-ubuntu24.04", add_python="3.12")
    .apt_install("git", "wget", "clang", "gcc", "g++", "make", "zlib1g-dev", "zstd",
                 "hmmer", "python3.12-dev")
    .run_commands(
        "git clone --depth 1 https://github.com/google-deepmind/alphafold3.git /opt/af3",
        "cd /opt/af3 && pip install --upgrade pip && pip install .",
        "cd /opt/af3 && build_data",
    )
    .pip_install("kaggle", "rdkit", "gemmi")
)

app = modal.App(APP_NAME, image=image)
weights_vol = modal.Volume.from_name("cyp-af3-weights", create_if_missing=True)
out_vol = modal.Volume.from_name("cyp-pool", create_if_missing=True)


@app.function(volumes={"/weights": weights_vol},
              secrets=[modal.Secret.from_name("kaggle")],
              timeout=3600, retries=1, max_containers=1)
def fetch_weights(kaggle_ref: str = "google-deepmind/alphafold-3/jax/default") -> dict:
    """Pull AF3 parameters from Kaggle into the private weights volume, once.

    `kaggle_ref` is overridable because model slugs move; a hardcoded slug that 404s is a
    confusing failure to debug from inside a container.
    """
    import subprocess

    dest = Path("/weights/af3")
    if any(dest.glob("*.bin")) or any(dest.glob("*.npz")):
        weights_vol.commit()
        return {"ok": True, "cached": True,
                "files": sorted(p.name for p in dest.iterdir())[:20]}
    dest.mkdir(parents=True, exist_ok=True)
    cp = subprocess.run(
        ["kaggle", "models", "instances", "versions", "download", kaggle_ref,
         "-p", str(dest), "--untar"],
        capture_output=True, text=True, timeout=3000)
    files = sorted(p.name for p in dest.rglob("*") if p.is_file())
    weights_vol.commit()
    return {"ok": bool(files), "cached": False, "n_files": len(files),
            "files": files[:20], "stderr": "" if files else cp.stderr[-2000:]}


@app.function(gpu=GPU, volumes={"/weights": weights_vol, "/out": out_vol},
              timeout=FN_TIMEOUT, retries=1, max_containers=MAX_CONTAINERS)
def cofold_batch(batch: list[dict], tag: str, msa: str, seeds: list[int]) -> list[dict]:
    """Run a chunk of ligands. Commits after each one so a preemption loses at most one.

    AF3 is invoked with `--norun_data_pipeline` and the MSA supplied inline: the target is
    a single sequence, so running the data pipeline per ligand would recompute the same
    alignment hundreds of times.
    """
    import subprocess

    model_dir = "/weights/af3"
    results, consecutive_failures = [], 0

    for lig in batch:
        sid = lig["id"]
        dest = Path("/out") / tag / sid
        if (dest / "DONE.json").exists():
            results.append({"id": sid, "status": "cached"})
            continue

        work = Path("/tmp") / sid
        (work / "in").mkdir(parents=True, exist_ok=True)
        (work / "out").mkdir(parents=True, exist_ok=True)
        job = {
            "name": sid, "modelSeeds": list(seeds), "dialect": "alphafold3", "version": 2,
            "sequences": [
                {"protein": {"id": "A", "sequence": lig["sequence"],
                             "unpairedMsa": msa, "pairedMsa": "", "templates": []}},
                {"ligand": {"id": "H", "ccdCodes": ["HEM"]}},
                {"ligand": {"id": "L", "smiles": lig["smiles"]}},
            ],
        }
        (work / "in" / f"{sid}.json").write_text(json.dumps(job))

        t0 = time.time()
        try:
            cp = subprocess.run(
                ["python", "/opt/af3/run_alphafold.py",
                 "--input_dir", str(work / "in"), "--output_dir", str(work / "out"),
                 "--model_dir", model_dir, "--norun_data_pipeline"],
                capture_output=True, text=True, timeout=INNER_TIMEOUT)
            rc, err = cp.returncode, cp.stderr[-3000:]
        except subprocess.TimeoutExpired:
            rc, err = -9, f"run_alphafold.py exceeded {INNER_TIMEOUT}s"
        except Exception as exc:
            rc, err = -1, f"{type(exc).__name__}: {exc}"

        dest.mkdir(parents=True, exist_ok=True)
        n_cif = 0
        for p in (work / "out").rglob("*"):
            if p.suffix in (".cif", ".json"):
                (dest / f"{p.parent.name}__{p.name}").write_bytes(p.read_bytes())
                n_cif += p.suffix == ".cif"

        ok = rc == 0 and n_cif > 0
        rec = {"id": sid, "returncode": rc, "n_structures": n_cif,
               "seconds": round(time.time() - t0, 1),
               "status": "ok" if ok else "failed",
               "stderr_tail": "" if ok else err}
        (dest / "DONE.json").write_text(json.dumps(rec, indent=1))
        out_vol.commit()
        results.append(rec)

        consecutive_failures = 0 if ok else consecutive_failures + 1
        if consecutive_failures >= MAX_CONSECUTIVE_FAILURES:
            # Every input failing the same way means the configuration is wrong, not the
            # ligands. Continuing would pay full wall-clock to produce nothing.
            results.append({"status": "aborted",
                            "reason": f"{consecutive_failures} consecutive failures; "
                                      "circuit breaker tripped"})
            break

    return results


def submit(csv_path: str, tag: str, seeds: tuple[int, ...] = (1, 2, 3, 4),
           chunk: int = 6, dry_run: bool = False) -> None:
    import csv as _csv

    from cypstruct import budget
    from cypstruct.chem import standardize
    from cypstruct.targets import fetch_sequences

    seq = fetch_sequences()["cyp3a4"]
    ligs = []
    with open(csv_path) as fh:
        for row in _csv.DictReader(fh):
            smi = standardize(row.get("smiles") or row.get("SMILES") or "")
            lid = row.get("id") or row.get("Molecule_Name")
            if smi and lid:
                ligs.append({"id": lid, "smiles": smi, "sequence": seq})

    est = budget.estimate("af3_cofold", len(ligs), len(seeds))
    ok, why, _ = budget.preflight_hours("af3_cofold", est)
    print(f"ligands: {len(ligs)}  seeds: {len(seeds)}  est {est:.2f} GPU-h")
    if not ok:
        raise SystemExit(f"PREFLIGHT REFUSED: {why}")
    if dry_run:
        print(json.dumps(ligs[:2], indent=1)[:800])
        return

    msa_path = REPO / "data" / "reference" / "cyp3a4.a3m"
    if not msa_path.exists():
        raise SystemExit(
            f"missing {msa_path}. AF3 runs with --norun_data_pipeline and needs the MSA "
            "inline; export it from the Boltz cache volume first "
            "(scripts/cofold/modal_boltz.py warms it).")
    msa = msa_path.read_text()

    run_id = f"{tag}-af3-{int(time.time())}"
    budget.record(run_id, "modal", "af3_cofold", len(ligs) * len(seeds), est,
                  note=f"app={APP_NAME}", n_jobs=len(ligs), samples=len(seeds), app=APP_NAME)
    try:
        with app.run():
            w = fetch_weights.remote()
            print("weights:", json.dumps(w)[:400], flush=True)
            if not w.get("ok"):
                raise SystemExit("AF3 weights not available; not launching the array")
            chunks = [ligs[i:i + chunk] for i in range(0, len(ligs), chunk)]
            n_ok = n_fail = 0
            for res in cofold_batch.starmap(
                    [(c, tag, msa, list(seeds)) for c in chunks],
                    return_exceptions=True):
                if isinstance(res, Exception):
                    n_fail += 1
                    continue
                n_ok += sum(1 for r in res if r.get("status") in ("ok", "cached"))
                n_fail += sum(1 for r in res if r.get("status") == "failed")
                if any(r.get("status") == "aborted" for r in res):
                    print("  !! circuit breaker tripped in a chunk", flush=True)
        budget.close(run_id, "done", note=f"{n_ok} ok, {n_fail} failed")
        print(f"DONE: {n_ok} ok, {n_fail} failed")
    except BaseException as exc:
        budget.close(run_id, "failed", note=str(exc)[:300])
        raise


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["submit", "plan", "weights"])
    ap.add_argument("--csv")
    ap.add_argument("--tag", default="af3")
    ap.add_argument("--seeds", default="1,2,3,4")
    a = ap.parse_args()
    if a.cmd == "weights":
        with app.run():
            print(json.dumps(fetch_weights.remote(), indent=2))
    else:
        submit(a.csv, a.tag, tuple(int(s) for s in a.seeds.split(",")),
               dry_run=(a.cmd == "plan"))
