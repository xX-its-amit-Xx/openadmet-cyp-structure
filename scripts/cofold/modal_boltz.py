"""Boltz-2 co-folding on Modal, with steering constraints and hard runaway guards.

Why self-hosted rather than the hosted Boltz API: the hosted endpoint accepts templates
but does not expose `constraints`. Constraints are the entire point here. OpenADMET's own
co-folding analysis finds that on CYP3A4 the heme is placed correctly and the backbone is
accurate — the failure is **ligand orientation**, including near-180-degree flips, with
Boltz-2 reaching only ~58% of poses under 2 A BiSyRMSD. Our parse of 116 deposited CYP3A4
entries shows **83 of 101 unique ligands coordinate the iron, through nitrogen in almost
every case**, and `cypstruct.chem.coordinating_atoms` identifies the right donor for
**82 of those 83** (98.8% recall).

So for the dominant binding mode we can name the atom that must face the iron. That turns
an orientation failure into a constraint. This module submits both a STEERED and an
UNSTEERED arm for every ligand so the effect is measured rather than assumed.

Runaway guards, since an unattended Modal app is the specific thing that burns a month of
credits overnight:
  - explicit `timeout`, `retries=1`, and `max_containers` on the remote function;
  - `.map()` over a finite list, never a polling loop;
  - per-input idempotency against the output volume, so a resumed run does no repeat work;
  - a preflight against `cypstruct.budget`, which refuses launches over the monthly cap;
  - the run is registered in the local ledger BEFORE it starts, so
    `scripts/ops/watchdog.py` can kill it even if this process dies.

Usage:
    modal run scripts/cofold/modal_boltz.py --help-plan          # print what would run
    python scripts/cofold/modal_boltz.py submit --csv ligands.csv --tag pilot1 --samples 20
    python scripts/cofold/modal_boltz.py collect --tag pilot1
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

APP_NAME = "cyp-cofold-boltz"

# --- resource envelope -----------------------------------------------------
GPU = os.environ.get("CYP_BOLTZ_GPU", "A100-40GB")
MAX_CONTAINERS = int(os.environ.get("CYP_BOLTZ_MAX_CONTAINERS", "8"))
FN_TIMEOUT = int(os.environ.get("CYP_BOLTZ_TIMEOUT", "2400"))   # seconds per ligand job

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "wget", "build-essential")
    .pip_install(
        "torch==2.5.1",
        extra_options="--index-url https://download.pytorch.org/whl/cu124",
    )
    .pip_install("boltz[cuda]==2.2.1", "numpy<2", "rdkit", "gemmi", "pyyaml")
    .env({"HF_HOME": "/cache/hf", "BOLTZ_CACHE": "/cache/boltz"})
)

app = modal.App(APP_NAME, image=image)

# Weights + CCD + the (single) target MSA. Persisted so containers do not re-download.
cache_vol = modal.Volume.from_name("cyp-boltz-cache", create_if_missing=True)
# Predicted structures. Read back by `collect`, then pushed to OneDrive.
out_vol = modal.Volume.from_name("cyp-pool", create_if_missing=True)


def _clean_a3m(raw: bytes) -> bytes:
    """Strip NUL bytes and normalise the tail of an a3m alignment.

    The MSA server's output arrived with a single trailing NUL byte. Boltz's a3m
    parser maps every character through `const.prot_letter_to_token`, so one stray
    byte raises a KeyError after the schema has parsed and the GPU is already
    allocated. A 4 MB alignment, 6,979 sequences, every one of them valid, killed
    by its last byte.

    Applied on write AND on read, so a cache written before this existed is healed
    rather than needing to be recomputed.
    """
    return raw.replace(b"\x00", b"").rstrip() + b"\n"


# ==========================================================================
# remote
# ==========================================================================


@app.function(volumes={"/cache": cache_vol}, timeout=3600, retries=1, max_containers=1)
def warm_cache(sequence: str) -> dict:
    """One-time: download Boltz weights + CCD, and compute the target MSA ONCE.

    The MSA is the expensive shared input and the target is a single sequence, so
    computing it per ligand would repeat the same ColabFold query thousands of times
    and get us rate-limited for no benefit. It is computed here and reused by every job.
    """
    import subprocess

    Path("/cache/boltz").mkdir(parents=True, exist_ok=True)
    msa_dir = Path("/cache/msa")
    msa_dir.mkdir(parents=True, exist_ok=True)
    tgt = msa_dir / "cyp3a4.a3m"

    if not tgt.exists():
        work = Path("/tmp/msa_probe")
        work.mkdir(parents=True, exist_ok=True)
        y = work / "probe.yaml"
        y.write_text(
            "version: 1\nsequences:\n  - protein:\n      id: A\n"
            f"      sequence: {sequence}\n"
        )
        # `boltz predict` with the MSA server writes the alignment into its output dir;
        # a trivial protein-only prediction is the cheapest way to materialise it.
        cp = subprocess.run(
            ["boltz", "predict", str(y), "--use_msa_server", "--out_dir", str(work),
             "--diffusion_samples", "1", "--recycling_steps", "1",
             "--cache", "/cache/boltz"],
            capture_output=True, text=True, timeout=3000)
        found = list(work.rglob("*.a3m"))
        if found:
            tgt.write_bytes(_clean_a3m(found[0].read_bytes()))
        else:
            return {"ok": False, "stderr": cp.stderr[-3000:], "stdout": cp.stdout[-2000:]}

    else:
        # heal a cache written before the sanitiser existed
        cleaned = _clean_a3m(tgt.read_bytes())
        if cleaned != tgt.read_bytes():
            tgt.write_bytes(cleaned)

    cache_vol.commit()
    return {"ok": True, "msa_bytes": tgt.stat().st_size,
            "weights": sorted(p.name for p in Path("/cache/boltz").glob("*"))[:20]}


@app.function(
    gpu=GPU,
    volumes={"/cache": cache_vol, "/out": out_vol},
    timeout=FN_TIMEOUT,
    retries=1,                 # one retry, not Modal's silent-forever default posture
    max_containers=MAX_CONTAINERS,
)
def cofold(spec: dict) -> dict:
    """Run one (ligand x arm) job. Idempotent: returns immediately if output exists."""
    import subprocess
    import traceback

    job_id = spec["job_id"]
    dest = Path("/out") / spec["tag"] / job_id
    done_marker = dest / "DONE.json"
    if done_marker.exists():
        return {"job_id": job_id, "status": "cached", **json.loads(done_marker.read_text())}

    t0 = time.time()
    work = Path("/tmp") / job_id
    work.mkdir(parents=True, exist_ok=True)
    # Boltz resolves a relative `msa:` path against the PROCESS working directory, not
    # against the yaml's directory. Writing the alignment next to the yaml and naming it
    # relatively therefore fails with "MSA file cyp3a4.a3m not found" even though the file
    # is right there. Substitute an absolute path at run time.
    msa_src = Path("/cache/msa/cyp3a4.a3m")
    local_msa = work / "cyp3a4.a3m"
    if msa_src.exists():
        # sanitise on READ as well as on write: the cached copy predates the fix
        local_msa.write_bytes(_clean_a3m(msa_src.read_bytes()))
    y = work / "input.yaml"
    y.write_text(spec["yaml"].replace("MSA_PATH", str(local_msa)))

    cmd = ["boltz", "predict", str(y), "--out_dir", str(work),
           "--diffusion_samples", str(spec.get("samples", 5)),
           "--recycling_steps", str(spec.get("recycling", 3)),
           "--sampling_steps", str(spec.get("sampling_steps", 200)),
           "--output_format", "mmcif", "--cache", "/cache/boltz",
           "--override"]
    if spec.get("seed") is not None:
        cmd += ["--seed", str(spec["seed"])]
    if not msa_src.exists():
        cmd += ["--use_msa_server"]

    try:
        cp = subprocess.run(cmd, capture_output=True, text=True,
                            timeout=max(300, FN_TIMEOUT - 180))
        rc, err = cp.returncode, cp.stderr[-4000:]
    except subprocess.TimeoutExpired:
        rc, err = -9, "boltz predict exceeded the in-container timeout"
    except Exception:
        rc, err = -1, traceback.format_exc()[-4000:]

    dest.mkdir(parents=True, exist_ok=True)
    n_struct = 0
    for p in work.rglob("*"):
        if p.suffix in (".cif", ".json", ".npz") and "predictions" in str(p):
            (dest / p.name).write_bytes(p.read_bytes())
            n_struct += p.suffix == ".cif"

    rec = {"job_id": job_id, "returncode": rc, "n_structures": n_struct,
           "seconds": round(time.time() - t0, 1), "arm": spec.get("arm"),
           "ligand": spec.get("ligand"), "seed": spec.get("seed"),
           "status": "ok" if (rc == 0 and n_struct) else "failed",
           "stderr_tail": "" if (rc == 0 and n_struct) else err}
    (dest / "DONE.json").write_text(json.dumps(rec, indent=1))
    out_vol.commit()
    return rec


@app.function(volumes={"/cache": cache_vol}, timeout=1800, retries=1, max_containers=1,
              cpu=2.0)
def probe_atom_names(smiles_list: list[str]) -> dict:
    """Ask Boltz what it will call each heavy atom of these SMILES ligands.

    Replicating the naming rule locally is not reliable: the name is
    `SYMBOL + CanonicalRankAtoms(mol)[i] + 1`, and the rank depends on whether hydrogens
    were added before ranking and on the RDKit version doing the ranking. Guessing wrong
    costs a whole GPU batch and shows up only as a KeyError deep in the schema parser.

    So we ask the authority. This runs on CPU, takes seconds, and gives names for the
    entire ligand set at once.
    """
    import traceback

    out: dict = {}
    try:
        from boltz.data.parse.schema import parse_boltz_schema
        from boltz.data.types import MSA  # noqa: F401  (import probe)
    except Exception:
        return {"_error": "could not import boltz parser: " + traceback.format_exc()[-1500:]}

    from rdkit import Chem
    from rdkit.Chem import AllChem

    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                out[smi] = {"_error": "unparseable"}
                continue
            # Mirror Boltz's own construction exactly, then read the names it assigns.
            work = Chem.AddHs(mol)
            AllChem.EmbedMolecule(work, randomSeed=1)
            canon_noh = list(AllChem.CanonicalRankAtoms(mol))
            canon_h = list(AllChem.CanonicalRankAtoms(work))
            out[smi] = {
                "n_heavy": mol.GetNumAtoms(),
                "names_no_h": {i: mol.GetAtomWithIdx(i).GetSymbol().upper()
                               + str(canon_noh[i] + 1) for i in range(mol.GetNumAtoms())},
                "names_with_h": {i: work.GetAtomWithIdx(i).GetSymbol().upper()
                                 + str(canon_h[i] + 1) for i in range(mol.GetNumAtoms())},
            }
        except Exception:
            out[smi] = {"_error": traceback.format_exc()[-600:]}
    return out


@app.function(volumes={"/cache": cache_vol}, timeout=1800, retries=1, max_containers=1,
              cpu=2.0)
def probe_yaml(variants: dict) -> dict:
    """Run Boltz's REAL schema parser over candidate yamls, on CPU, and report which parse.

    This is the cheap adjudicator. A constraint that names a nonexistent ligand atom
    raises deep inside `token_spec_to_ids`, and on the GPU path that costs an entire
    batch to discover. Here it costs seconds, and it settles questions like "does Boltz
    rank atoms before or after adding hydrogens" by asking rather than by reasoning.
    """
    import pickle
    import traceback
    from pathlib import Path as _P

    import yaml as _y

    out: dict = {}
    try:
        from boltz.data.parse.schema import parse_boltz_schema
    except Exception:
        return {"_error": "import failed: " + traceback.format_exc()[-1200:]}

    mol_dir = _P("/cache/boltz/mols")
    ccd = {}
    for cand in (_P("/cache/boltz/ccd.pkl"), _P("/cache/boltz/ccd.json")):
        if cand.exists():
            try:
                ccd = pickle.loads(cand.read_bytes())
            except Exception:
                ccd = {}
            break

    # The parser's signature has changed across Boltz releases (`boltz2=` exists in some
    # versions and not others). Introspect rather than pin, so this probe keeps working
    # when the image is rebuilt on a newer Boltz.
    import inspect

    params = list(inspect.signature(parse_boltz_schema).parameters)
    out["_signature"] = params

    def _call(label, data):
        kwargs = {}
        if "boltz_2" in params:
            kwargs["boltz_2"] = True
        elif "boltz2" in params:
            kwargs["boltz2"] = True
        args = [label, data, ccd]
        if "mol_dir" in params:
            args.append(mol_dir)
        return parse_boltz_schema(*args, **kwargs)

    for label, text in variants.items():
        try:
            data = _y.safe_load(text.replace("MSA_PATH", "/cache/msa/cyp3a4.a3m"))
            target = _call(label, data)
            n = None
            try:
                n = int(len(target.structure.atoms))
            except Exception:
                pass
            out[label] = {"ok": True, "n_atoms": n}
        except Exception as exc:
            out[label] = {"ok": False, "error": f"{type(exc).__name__}: {exc}"[:300]}
    return out


# ==========================================================================
# local driver
# ==========================================================================


def build_yaml(sequence: str, smiles: str, *, steer: bool, donor_atom_name: str | None,
               template_cif: str | None, axial_cys: int = 442,
               pocket_residues: list[int] | None = None,
               msa_path: str = "MSA_PATH") -> str:
    """Compose a Boltz-2 input.

    The heme is always present and always bonded to the axial cysteine — that bond is
    not a steering choice, it is the resting state of the enzyme, and leaving it to be
    inferred invites a heme that floats.

    The STEERED arm adds, on top of that:
      - a template from the nearest holo CYP3A4 (fixes the plastic F/G region to a real
        observed conformation rather than an average),
      - a pocket restraint listing the I-helix and Phe-cluster residues,
      - a forced contact from the heme iron to the predicted donor atom at 2.4 A.
    The last line is the one that fixes orientation for the majority binding mode.
    """
    L = ["version: 1", "sequences:",
         "  - protein:", "      id: A", f"      sequence: {sequence}",
         f"      msa: {msa_path}",
         "  - ligand:", "      id: H", "      ccd: HEM",
         "  - ligand:", "      id: L", f"      smiles: '{smiles}'"]

    cons = ["constraints:",
            "  - bond:",
            f"      atom1: [A, {axial_cys}, SG]",
            "      atom2: [H, 1, FE]"]
    if steer:
        if pocket_residues:
            contacts = ", ".join(f"[A, {r}]" for r in pocket_residues)
            cons += ["  - pocket:", "      binder: L",
                     f"      contacts: [{contacts}]",
                     "      max_distance: 8.0", "      force: false"]
        if donor_atom_name:
            # `contact` tokens take TWO elements, [CHAIN_ID, RES_IDX-or-ATOM_NAME] —
            # unlike `bond`, whose atoms take three. Passing three here raises
            # "too many values to unpack (expected 2)" inside the Boltz schema parser,
            # which `boltz predict` catches per input and skips, so the process exits 0
            # with no structures. That is why the first run reported rc=0, n=0 in 25s.
            cons += ["  - contact:",
                     "      token1: [H, FE]",
                     f"      token2: [L, {donor_atom_name}]",
                     "      max_distance: 2.4", "      force: true"]
    L += cons
    if steer and template_cif:
        L += ["templates:", f"  - cif: {template_cif}", "    force: true", "    threshold: 2.0"]
    return "\n".join(L) + "\n"


def plan(csv_path: str, tag: str, samples: int, seeds: list[int],
         arms: tuple[str, ...] = ("steered", "unsteered")) -> list[dict]:
    """Build the finite job list. No job is created without a resolvable job_id."""
    import csv as _csv

    from cypstruct.chem import boltz_atom_name, coordinating_atoms, standardize
    from cypstruct.targets import CYP3A4_POCKET, fetch_sequences

    seq = fetch_sequences()["cyp3a4"]
    pocket = sorted(set(CYP3A4_POCKET["i_helix"] + CYP3A4_POCKET["phe_cluster"]
                        + CYP3A4_POCKET["polar"]))
    jobs = []
    with open(csv_path) as fh:
        for row in _csv.DictReader(fh):
            lid = row.get("id") or row.get("idx") or row.get("Molecule_Name")
            smi = standardize(row.get("smiles") or row.get("SMILES") or "")
            if not smi or not lid:
                continue
            donors = coordinating_atoms(smi, top_k=1)
            # Boltz names ligand atoms from the CCD/SMILES parse; for a SMILES ligand the
            # atom name is the element plus a 1-based ordinal over that element.
            donor_name = boltz_atom_name(smi, donors[0].atom_idx) if donors else None
            for arm in arms:
                if arm == "steered" and donor_name is None:
                    continue     # nothing to steer with; the unsteered arm covers it
                for seed in seeds:
                    jobs.append(dict(
                        job_id=f"{lid}__{arm}__s{seed}", tag=tag, ligand=lid, arm=arm,
                        seed=seed, samples=samples,
                        yaml=build_yaml(seq, smi, steer=(arm == "steered"),
                                        donor_atom_name=donor_name,
                                        template_cif=None, pocket_residues=pocket),
                    ))
    return jobs


def submit(csv_path: str, tag: str, samples: int = 10, seeds: tuple[int, ...] = (1, 2, 3),
           dry_run: bool = False) -> None:
    from cypstruct import budget

    jobs = plan(csv_path, tag, samples, list(seeds))
    units = len(jobs) * samples
    est = budget.estimate("boltz2_cofold", len(jobs), samples)
    ok, why, _ = budget.preflight_hours("boltz2_cofold", est)
    print(f"planned jobs: {len(jobs)}  ({units} ligand-samples)  est {est:.2f} GPU-h")
    if not ok:
        raise SystemExit(f"PREFLIGHT REFUSED: {why}")
    if dry_run:
        print(json.dumps(jobs[0], indent=1)[:1500])
        return

    run_id = f"{tag}-{int(time.time())}"
    budget.record(run_id, "modal", "boltz2_cofold", units, est,
                  note=f"app={APP_NAME} tag={tag}", n_jobs=len(jobs),
                  samples=samples, app=APP_NAME)
    try:
        with app.run():
            from cypstruct.targets import fetch_sequences
            print("warming cache / computing the shared target MSA ...", flush=True)
            w = warm_cache.remote(fetch_sequences()["cyp3a4"])
            print("  ", json.dumps(w)[:400], flush=True)
            if not w.get("ok"):
                raise SystemExit("cache warm failed; not launching the array")

            done = fail = 0
            for rec in cofold.map(jobs, order_outputs=False, return_exceptions=True):
                if isinstance(rec, Exception):
                    fail += 1
                    continue
                done += rec.get("status") in ("ok", "cached")
                fail += rec.get("status") == "failed"
                if (done + fail) % 20 == 0:
                    print(f"  {done} ok / {fail} failed of {len(jobs)}", flush=True)
        budget.close(run_id, "done", note=f"{done} ok, {fail} failed")
        print(f"DONE: {done} ok, {fail} failed")
    except BaseException as exc:
        budget.close(run_id, "failed", note=str(exc)[:300])
        raise


def collect(tag: str) -> None:
    """Pull results out of the Modal volume and push them to OneDrive."""
    from cypstruct.storage import Batch

    with Batch(f"pool/boltz2/{tag}") as b:
        vol = modal.Volume.from_name("cyp-pool")
        n = 0
        for entry in vol.iterdir(f"/{tag}"):
            for f in vol.iterdir(f"/{tag}/{entry.path.split('/')[-1]}"):
                data = b"".join(vol.read_file(f.path))
                p = b.path / f.path.lstrip("/")
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_bytes(data)
                n += 1
            b.check()
        print(f"staged {n} files -> pushing to OneDrive")
    print("collected.")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["submit", "collect", "plan"])
    ap.add_argument("--csv")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--samples", type=int, default=10)
    ap.add_argument("--seeds", default="1,2,3")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    seeds = tuple(int(s) for s in a.seeds.split(","))
    if a.cmd in ("submit", "plan"):
        submit(a.csv, a.tag, a.samples, seeds, dry_run=(a.dry_run or a.cmd == "plan"))
    else:
        collect(a.tag)
