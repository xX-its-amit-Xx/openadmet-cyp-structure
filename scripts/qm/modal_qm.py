"""Tier-1 ligand QM (GFN2-xTB) on Modal CPU, one job per LIGAND.

Why this is a CPU app and why it is cheap. Tier 1 computes properties of the MOLECULE,
not of the pose (`docs/QM_SCORER_DESIGN.md` §3). The 87-ligand x 20-sample x 2-arm pool
is ~3,500 poses but **87 jobs here**, and each one is a handful of GFN2-xTB
optimizations on a drug-sized molecule — seconds to a couple of minutes on two cores.
There is no GPU anywhere in this file and there must never be one.

Toolchain, decided by measurement rather than by assumption (`--probe` reruns it):

  * **conda-forge `xtb` 6.7.1 via micromamba — CHOSEN.** Installs clean, builds in ~22 s,
    and its CLI gives geometry optimization, charge/multiplicity flags and Mulliken
    charges to a `charges` file. Every tier-1 quantity falls out of it with no SCF driven
    by hand.
  * `pip install tblite` — works (GFN2 single point on water returned -5.0705 Eh), but the
    Python API exposes single points only, so protonation geometries would need a
    separate optimizer. Kept as a documented fallback, not used.
  * `pip install xtb` (xtb-python) — needs the compiled library that conda-forge ships
    anyway, so it offers nothing over the CLI in a micromamba image.
  * `pip install qmdesc` — **rejected.** Two independent reasons: it imports
    `pkg_resources`, removed in setuptools >= 81, so it fails at import on any current
    image (`ModuleNotFoundError: No module named 'pkg_resources'`); and it pulls torch
    plus ~3 GB of CUDA wheels to predict charges for a CPU job. Its descriptors are also
    charges and bond orders, not proton affinities, so it could not supply the
    sigma-donor term that is the point of §2B.

Runaway guards, same posture as `scripts/cofold/modal_boltz.py`:
  - explicit `timeout`, `retries=1`, `max_containers` on the remote function;
  - `.map()` over a finite job list, never a polling loop;
  - per-ligand idempotency against the output volume, so a resumed run repeats no work;
  - `cypstruct.budget` preflight under kind `xtb_pocket_singlepoint`, charged against
    the **CPU** cap, which refuses the launch before it starts;
  - the run is in the local ledger before it launches, so `scripts/ops/watchdog.py` can
    kill it even if this process dies.

Usage:
    python scripts/qm/modal_qm.py probe                        # toolchain + PA sanity
    python scripts/qm/modal_qm.py plan   --csv ligands.csv --tag t1
    python scripts/qm/modal_qm.py submit --csv ligands.csv --tag t1
    python scripts/qm/modal_qm.py collect --tag t1             # -> data/processed/qm/
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import modal

# Modal imports this module INSIDE the container, where it lands at /root/modal_qm.py.
# `parents[2]` then raises IndexError and every container dies on import, which Modal
# answers by restarting them — the crash-loop that looks like "still starting" while it
# burns credits. Resolve the repo path only when it actually exists.
_here = Path(__file__).resolve()
REPO = _here.parents[2] if len(_here.parents) > 2 else None
if REPO is not None and (REPO / "src" / "cypstruct").is_dir():
    sys.path.insert(0, str(REPO / "src"))

APP_NAME = "cyp-qm-tier1"

# --- resource envelope ------------------------------------------------------
MAX_CONTAINERS = int(os.environ.get("CYP_QM_MAX_CONTAINERS", "8"))
FN_TIMEOUT = int(os.environ.get("CYP_QM_TIMEOUT", "1800"))     # seconds per ligand
CPU_PER_JOB = float(os.environ.get("CYP_QM_CPU", "2.0"))

# xtb is pinned. GFN2 total energies shift slightly between xtb releases, and the proton
# affinities here are differences of total energies cached to JSON and compared across
# runs — an unpinned image would silently make old and new numbers incomparable.
image = (
    modal.Image.micromamba(python_version="3.11")
    # `requests` is not used by this app, but `cypstruct/__init__.py` imports `targets`,
    # which imports it at module scope. Importing anything from the package therefore
    # pulls it in, and without it every container dies at import with
    # ModuleNotFoundError before a single xtb call is made.
    .micromamba_install("xtb=6.7.1", "rdkit=2024.09.5", "numpy", "requests",
                        channels=["conda-forge"])
    .env({"OMP_NUM_THREADS": "2", "MKL_NUM_THREADS": "2",
          "OMP_STACKSIZE": "4G", "XTBPATH": "/tmp"})
    .add_local_python_source("cypstruct")
)

app = modal.App(APP_NAME, image=image)

# Tier-1 output: one small JSON per ligand. Separate from `cyp-pool` (which holds bulky
# predicted structures) because these are re-read constantly by the scorer and are tiny.
qm_vol = modal.Volume.from_name("cyp-qm", create_if_missing=True)


# ==========================================================================
# remote
# ==========================================================================


@app.function(volumes={"/qm": qm_vol}, timeout=FN_TIMEOUT, retries=1,
              max_containers=MAX_CONTAINERS, cpu=CPU_PER_JOB)
def ligand_qm_job(spec: dict) -> dict:
    """Tier 1 for ONE ligand. Idempotent against the output volume.

    The whole point of the cost architecture is that this is called once per molecule.
    If you ever find yourself mapping this over poses, stop: the pose-dependent half is
    `ligand_qm.read_against_pose`, which is a dict lookup and costs nothing.
    """
    import traceback

    from cypstruct.qmscore.ligand_qm import compute_ligand_qm

    lig = spec["ligand_id"]
    dest = Path("/qm") / spec["tag"]
    dest.mkdir(parents=True, exist_ok=True)
    out = dest / f"{lig}.json"
    if out.exists() and not spec.get("force"):
        try:
            rec = json.loads(out.read_text())
            if rec.get("status") in ("ok", "partial"):
                return {"ligand_id": lig, "status": "cached",
                        "n_donors": len(rec.get("donors", [])),
                        "seconds": rec.get("seconds")}
        except json.JSONDecodeError:
            pass       # truncated cache entry: fall through and recompute

    t0 = time.time()
    try:
        q = compute_ligand_qm(
            spec["smiles"], lig,
            n_conformers=spec.get("n_conformers", 4),
            solvent=spec.get("solvent") or None,
            also_solvated=spec.get("also_solvated", False),
            n_bde_sites=spec.get("n_bde_sites", 0),
            workdir=Path("/tmp") / lig,          # container-local, never the volume
        )
        q.save(out)
        qm_vol.commit()
        best = None
        live = [d for d in q.donors if not d.rejected]
        if live:
            best = max((d.proton_affinity_kcal, d.atom_idx) for d in live
                       if d.proton_affinity_kcal == d.proton_affinity_kcal)[1] \
                if any(d.proton_affinity_kcal == d.proton_affinity_kcal for d in live) else None
        return {"ligand_id": lig, "status": q.status, "n_donors": len(q.donors),
                "n_rejected": sum(1 for d in q.donors if d.rejected),
                "best_donor_idx": best, "n_xtb_calls": q.n_xtb_calls,
                "seconds": round(time.time() - t0, 1), "error": q.error[:200]}
    except Exception:
        return {"ligand_id": lig, "status": "failed",
                "seconds": round(time.time() - t0, 1),
                "error": traceback.format_exc()[-800:]}


@app.function(timeout=1200, retries=1, max_containers=1, cpu=2.0)
def probe() -> dict:
    """Verify the toolchain and check that proton affinities order as chemistry predicts.

    The four prototypes are the ones §2B of the design names as the discrimination the
    term has to reproduce: imidazole > pyridine > thiazole > oxazole. If that ordering
    ever breaks, the sigma-donor term is not measuring what it claims to and nothing
    downstream should be trusted.
    """
    import subprocess

    from cypstruct.qmscore.ligand_qm import compute_ligand_qm, xtb_available

    v = subprocess.run(["xtb", "--version"], capture_output=True, text=True)
    ver = next((ln.strip() for ln in (v.stdout + v.stderr).splitlines()
                if "version" in ln.lower()), "?")

    probes = {
        "imidazole": "c1c[nH]cn1",
        "pyridine": "c1ccncc1",
        "thiazole": "c1cscn1",
        "oxazole": "c1cocn1",
        "2,6-lutidine": "Cc1cccc(C)n1",
        "ketoconazole": "CC(=O)N1CCN(CC1)c1ccc(OC[C@@H]2CO[C@](Cn3ccnc3)(O2)c2ccc(Cl)cc2Cl)cc1",
        "clotrimazole": "Clc1ccccc1C(n1ccnc1)(c1ccccc1)c1ccccc1",
        "testosterone": "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O",
    }
    out = {"xtb": ver, "xtb_available": xtb_available(), "ligands": {}}
    for name, smi in probes.items():
        q = compute_ligand_qm(smi, name, n_conformers=3,
                              workdir=Path("/tmp") / name.replace(",", "_"))
        out["ligands"][name] = {
            "status": q.status, "n_heavy": q.n_heavy, "seconds": q.seconds,
            "n_xtb_calls": q.n_xtb_calls, "error": q.error[:200],
            "donors": [{"idx": d.atom_idx, "pattern": d.smarts_pattern,
                        "rejected": d.rejected,
                        "pa_kcal": d.proton_affinity_kcal,
                        "fukui_minus": round(d.fukui_minus, 4)
                        if d.fukui_minus == d.fukui_minus else None,
                        "vbur_metal": round(d.vbur_at_metal, 1)
                        if d.vbur_at_metal == d.vbur_at_metal else None}
                       for d in q.donors],
            "easiest_ch": sorted(
                [{"idx": c.atom_idx, "env": c.environment, "bde": c.bde_proxy_kcal}
                 for c in q.carbons], key=lambda c: c["bde"])[:3],
        }
    return out


@app.function(timeout=900, retries=1, max_containers=1, cpu=2.0)
def calibrate_h_atom() -> dict:
    """GFN2 total energy of an isolated H atom, for the optional xTB BDE route.

    A constant of the method, not a fitted parameter. Recompute it whenever the pinned
    xtb version changes and update `ligand_qm._H_ATOM_ENERGY`.
    """
    from cypstruct.qmscore.ligand_qm import _xtb

    r = _xtb(Path("/tmp/hatom"), "1\n\nH 0.0 0.0 0.0\n", charge=0, uhf=1)
    return {"ok": r["ok"], "energy_hartree": r.get("energy"), "error": r.get("error", "")[:200]}


# ==========================================================================
# local driver
# ==========================================================================


def plan(csv_path: str, tag: str, **kw) -> list[dict]:
    """Build the finite, deduplicated job list.

    Deduplication is on the **standardised SMILES**, not on the row id: the same
    compound appears under several ids across the pool and the reference set, and tier-1
    descriptors depend only on the molecule. On the 87-ligand pool this is the difference
    between paying once and paying per appearance.
    """
    import csv as _csv

    from cypstruct.chem import standardize

    jobs: dict[str, dict] = {}
    with open(csv_path, newline="") as fh:
        for row in _csv.DictReader(fh):
            lid = (row.get("id") or row.get("idx") or row.get("ligand")
                   or row.get("Molecule_Name") or "").strip()
            raw = (row.get("smiles") or row.get("SMILES") or "").strip()
            smi = standardize(raw) if raw else None
            if not smi or not lid:
                continue
            if smi in jobs:
                jobs[smi].setdefault("aliases", []).append(lid)
                continue
            jobs[smi] = dict(ligand_id=lid, smiles=smi, tag=tag, **kw)
    return list(jobs.values())


def submit(csv_path: str, tag: str, dry_run: bool = False, **kw) -> None:
    from cypstruct import budget

    jobs = plan(csv_path, tag, **kw)
    if not jobs:
        raise SystemExit("no parseable ligands in the CSV; nothing to do")

    # Charged against the CPU cap, not the GPU cap. `xtb_pocket_singlepoint` is the
    # existing kind for xTB work in `budget.COST_HINTS` (0.002 CPU-hour per unit); a
    # tier-1 ligand is a few optimizations rather than one single point, so it is
    # charged as several units.
    units = len(jobs) * UNITS_PER_LIGAND
    est = budget.estimate("xtb_pocket_singlepoint", units)
    ok, why, _ = budget.preflight_hours("xtb_pocket_singlepoint", est, venue="modal",
                                        cap_key="modal_cpu_hours")
    print(f"planned jobs: {len(jobs)} unique ligands  ({units} xtb units)  "
          f"est {est:.3f} CPU-h")
    if not ok:
        raise SystemExit(f"PREFLIGHT REFUSED: {why}")
    if dry_run:
        print(json.dumps(jobs[:2], indent=1)[:1200])
        return

    run_id = f"{tag}-qm-{int(time.time())}"
    budget.record(run_id, "modal", "xtb_pocket_singlepoint", units, est,
                  note=f"app={APP_NAME} tag={tag}", n_jobs=len(jobs), app=APP_NAME)
    done = fail = 0
    errs: list[str] = []
    try:
        with app.run():
            t0 = time.time()
            for rec in ligand_qm_job.map(jobs, order_outputs=False,
                                         return_exceptions=True):
                if isinstance(rec, Exception):
                    fail += 1
                    if len(errs) < 5:
                        errs.append(f"{type(rec).__name__}: {rec}"[:250])
                        print(f"  [exception] {errs[-1]}", flush=True)
                    continue
                if rec.get("status") in ("ok", "partial", "cached"):
                    done += 1
                else:
                    fail += 1
                    if len(errs) < 5:
                        errs.append(f"{rec.get('ligand_id')}: {rec.get('error', '')[:200]}")
                        print(f"  [failed] {errs[-1]}", flush=True)
                if (done + fail) % 10 == 0:
                    print(f"  {done} ok / {fail} failed of {len(jobs)}", flush=True)
            wall = time.time() - t0
        # Book what it ACTUALLY cost, so `budget.calibrate_from_ledger` can replace the
        # guessed cost hint with a measured one for the next run.
        actual = wall * CPU_PER_JOB / 3600.0
        budget.close(run_id, "done", actual_gpu_hours=round(actual, 4),
                     note=f"{done} ok, {fail} failed, {wall/60:.1f} min wall")
        print(f"DONE: {done} ok, {fail} failed  ({wall/60:.1f} min wall, "
              f"~{actual:.3f} CPU-h billed)")
    except BaseException as exc:
        budget.close(run_id, "failed", note=str(exc)[:300])
        raise


# A tier-1 ligand costs roughly: 1 neutral optimization + 3 Fukui single points +
# one optimization per candidate donor (~1.4 on the reference set). Charged as 6 units
# of `xtb_pocket_singlepoint` to stay conservative against the cap.
UNITS_PER_LIGAND = 6


def collect(tag: str, dest: str | Path | None = None) -> None:
    """Pull the per-ligand JSONs down into the repo.

    These stay on D: rather than going to OneDrive: the whole tier-1 output for a few
    hundred ligands is a few MB of JSON that the scorer re-reads on every calibration
    pass. `storage.push()` is for bulk artifacts; this is not one. Nothing here is
    written through the `O:\\` drive letter, whose rclone VFS cache is unbounded and
    lives on the near-full C:.
    """
    from cypstruct.paths import DATA_PROCESSED, free_gb

    out = Path(dest) if dest else DATA_PROCESSED / "qm" / tag
    out.mkdir(parents=True, exist_ok=True)
    if free_gb(out) < 0.5:
        raise SystemExit(f"refusing to write to {out}: {free_gb(out):.2f} GB free")

    vol = modal.Volume.from_name("cyp-qm")
    n = 0
    for f in vol.iterdir(f"/{tag}"):
        data = b"".join(vol.read_file(f.path))
        (out / Path(f.path).name).write_bytes(data)
        n += 1
    print(f"collected {n} ligand QM records -> {out}")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=["probe", "plan", "submit", "collect", "calibrate"])
    ap.add_argument("--csv")
    ap.add_argument("--tag", default="tier1")
    ap.add_argument("--solvent", default="", help="ALPB solvent, e.g. water. Blank = gas.")
    ap.add_argument("--also-solvated", action="store_true",
                    help="additionally compute every PA in ALPB water (doubles cost)")
    ap.add_argument("--n-conformers", type=int, default=4)
    ap.add_argument("--n-bde-sites", type=int, default=0,
                    help="replace the graph BDE proxy with real GFN2 homolysis on the "
                         "K most abstractable carbons (K extra optimizations per ligand)")
    ap.add_argument("--force", action="store_true", help="ignore cached outputs")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    if a.cmd == "probe":
        with app.run():
            print(json.dumps(probe.remote(), indent=1))
    elif a.cmd == "calibrate":
        with app.run():
            print(json.dumps(calibrate_h_atom.remote(), indent=1))
    elif a.cmd == "collect":
        collect(a.tag)
    else:
        if not a.csv:
            raise SystemExit("--csv is required for plan/submit")
        submit(a.csv, a.tag, dry_run=(a.dry_run or a.cmd == "plan"),
               solvent=a.solvent, also_solvated=a.also_solvated,
               n_conformers=a.n_conformers, n_bde_sites=a.n_bde_sites, force=a.force)
