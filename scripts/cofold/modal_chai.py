"""Chai-1 co-folding on Modal — a THIRD ENGINE FOR THE TAIL, not a second steered arm.

Read the steering paragraph before using this. It is the reason this file exists in the
shape it does.

**Chai-1, not Chai-2.** Chai-2 is not open: the `chaidiscovery` GitHub org publishes only
`chai-lab` ("Chai-1, SOTA model for biomolecular structure prediction"), the HuggingFace
org publishes exactly one model (`chaidiscovery/chai-1`), and there is no `chai-2` project
on PyPI. Chai-2 is referenced in chai-lab's README only as a citation to an antibody-design
preprint, and was announced as early access "to select partners". It is also an antibody
*design* model, not a co-folder, so it would not be a drop-in here even if it were open.
Verified 2026-09-10 — see the URL block below.

**STEERING: Chai-1 CANNOT express the constraint this project's thesis depends on.**
Our steering signal is an atom-level distance from the heme iron to a named ligand donor
nitrogen (measured window 1.90-2.45 A). Chai-1's restraint file has an `@ATOM` syntax and
it parses, but `load_manual_restraints_for_chai1` reads only the residue position and
residue name for `contact` and `pocket` restraints and never references `atom_nameA` /
`atom_nameB`; only `covalent` restraints use them. Worse, ligand tokenization gives every
ligand atom the same `residue_index`, and `add_distance_restraint` asserts exactly one
matching residue token. So:

  - the default arm here is **`unsteered`**, and that is deliberate;
  - an opt-in **`pocket`** arm emits residue-level pocket restraints (ligand chain must sit
    near the CYP3A4 pocket residues). Honest expectation: near zero. OpenADMET's own
    analysis says co-folding already puts the ligand in the right *place* on CYP3A4 and
    gets its *orientation* wrong, so a restraint that only constrains place is aimed at the
    axis that is not failing. It is available so that claim is measured, not assumed;
  - an opt-in **`covsteer`** arm writes a `covalent` restraint from the heme FE to the
    predicted donor atom. This is the only atom-level lever Chai-1 exposes. It is
    **UNVERIFIED end to end**: both restraint examples Chai ships are protein->ligand, and
    Chai's own docs warn "Chai-1 was not trained on intra-chain bonds ... we have not
    evaluated whether specifying such bond information yields expected behaviors."

**Silent-drop hazard.** `TokenDistanceRestraint._generate` wraps restraint construction in
`except Exception: logger.error(f"Error {e} generating distance constraints: ...")` and then
falls through to an all-ignore matrix. A malformed or unsatisfiable restraint does not fail
the run — it produces an unsteered prediction that looks like a successful steered one. So
`cofold` greps stderr for that message and marks the job `status="restraints_dropped"`.
Without that check a steered arm quietly becomes a second unsteered arm and the A/B reads
as "steering does nothing".

**The heme is a SMILES here, not a CCD code.** Chai-1's `>ligand|` entity routes straight
to `get_lig_residues(smiles=...)`; the string is handed to RDKit and never looked up in the
CCD. CCD codes are accepted only for modified polymer residues `AAA(SEP)AAA` and for
`>glycan|` entities. So HEM goes in as the PDB chemical-component dictionary's OpenEye
canonical SMILES (43 heavy atoms, C34 N4 O4 Fe — checked to parse in RDKit locally). This
is a real quality difference from the Boltz and Protenix runners, which both get a proper
CCD heme with its ideal geometry, and it is the first thing to suspect if Chai's heme
placement is worse than theirs.

Runaway guards, same posture as the other runners in this directory:
  - explicit `timeout`, `retries=1`, `max_containers`;
  - an inner `subprocess.run(timeout=...)` STRICTLY LESS than the Modal function timeout,
    because the PXR AF3 runner omitted it and a hung child held an A100 for 8 hours;
  - a consecutive-failure circuit breaker per chunk, because 184 consecutive failures once
    ran to completion and paid full wall-clock;
  - `DONE.json` idempotency against the output volume, so a resume does no repeat work;
  - a `budget.preflight_hours` gate that REFUSES before launching, and a ledger row written
    BEFORE the launch so `scripts/ops/watchdog.py` can clean up if this process dies;
  - the first job is run ALONE and inspected before the rest are mapped. That both warms
    the weight cache without a cold-start stampede and turns a configuration error into one
    wasted job instead of a batch.

Sources (fetched 2026-09-10/11, quote-checked, not recalled):
  https://github.com/chaidiscovery/chai-lab
  https://raw.githubusercontent.com/chaidiscovery/chai-lab/main/README.md
  https://raw.githubusercontent.com/chaidiscovery/chai-lab/main/chai_lab/chai1.py
  https://raw.githubusercontent.com/chaidiscovery/chai-lab/main/chai_lab/data/dataset/inference_dataset.py
  https://raw.githubusercontent.com/chaidiscovery/chai-lab/main/chai_lab/data/parsing/msas/aligned_pqt.py
  https://raw.githubusercontent.com/chaidiscovery/chai-lab/main/chai_lab/data/parsing/restraints.py
  https://github.com/chaidiscovery/chai-lab/blob/main/examples/msas/README.md
  https://github.com/chaidiscovery/chai-lab/blob/main/examples/restraints/README.md
  https://github.com/chaidiscovery/chai-lab/blob/main/examples/covalent_bonds/README.md
  https://raw.githubusercontent.com/chaidiscovery/chai-lab/v0.6.1/requirements.in
  https://pypi.org/pypi/chai_lab/json                 (0.6.1, 2025-03-18)
  https://huggingface.co/api/models?author=chaidiscovery   (only chai-1; ungated)
  https://api.github.com/orgs/chaidiscovery/repos      (no chai-2 repo)
  https://files.rcsb.org/ligands/view/HEM.cif          (HEM SMILES + the FE atom name)

Licence: Apache 2.0 for BOTH code and model weights, commercial use permitted. Weights come
from `https://chaiassets.com/chai1-inference-depencencies/...` over plain HTTPS with no
token and no gate, so unlike AF3 nothing here depends on a personal account.

Usage:
    python scripts/cofold/modal_chai.py plan   --csv data/processed/validation_ligands.csv --tag chai87
    python scripts/cofold/modal_chai.py probe  --csv data/processed/validation_ligands.csv --tag chai87
    python scripts/cofold/modal_chai.py submit --csv data/processed/validation_ligands.csv --tag chai87 \
        --samples 5 --seeds 1 --arms unsteered

Scoring is `scripts/structure/collect_and_score.py --tag chai87`; this module writes the
`/<tag>/<ligand>__<arm>__s<seed>/` layout and the `confidence_<stem>.json` sidecars that
scorer reads. (Note it labels every pose `engine="boltz2"` internally — harmless for a
single-engine tag, but do not merge tags across engines without fixing that.)
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import modal

# Modal serialises this module and imports it INSIDE the container, where the file lands
# at /root/modal_chai.py — so `parents[2]` raises IndexError and every container dies on
# import. Modal then restarts them, which is precisely the crash-loop that burns credits
# while looking like "still starting". Resolve the repo path only when it exists.
_here = Path(__file__).resolve()
REPO = _here.parents[2] if len(_here.parents) > 2 else None
if REPO is not None and (REPO / "src" / "cypstruct").is_dir():
    sys.path.insert(0, str(REPO / "src"))

APP_NAME = "cyp-cofold-chai"

# --- resource envelope -----------------------------------------------------
# chai-lab's README names hardware classes, not a GB-per-token figure ("A100 80GB or
# H100 80GB or L40S 48GB ... A10 and A30 will work for smaller complexes"). A CYP3A4
# MONOMER plus heme plus one drug-sized ligand is ~520 tokens, comfortably a "smaller
# complex", so A100-40GB is the default and `low_memory=True` (chai's own default) is
# left on. UNVERIFIED: no published memory number for this token count.
GPU = os.environ.get("CYP_CHAI_GPU", "A100-40GB")
MAX_CONTAINERS = int(os.environ.get("CYP_CHAI_MAX_CONTAINERS", "6"))
FN_TIMEOUT = int(os.environ.get("CYP_CHAI_TIMEOUT", "3600"))      # seconds per chunk
INNER_TIMEOUT = max(600, FN_TIMEOUT - 600)   # always leave the container room to REPORT
MAX_CONSECUTIVE_FAILURES = 5
CHUNK = int(os.environ.get("CYP_CHAI_CHUNK", "4"))

# PDB chemical-component dictionary, HEM, OpenEye SMILES_CANONICAL. The CACTVS descriptor
# in the same file uses `|` dative-bond notation, which RDKit refuses outright (checked);
# this one parses to 43 heavy atoms = C34 N4 O4 Fe, which is heme b.
# https://files.rcsb.org/ligands/view/HEM.cif
HEM_SMILES = (
    "Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C)C=C)"
    "C(=C(C7=C2)C)C=C)C)CCC(=O)O"
)

# Chai lettered chains follow FASTA order. We always write protein, heme, ligand.
CHAIN_PROTEIN, CHAIN_HEME, CHAIN_LIGAND = "A", "B", "C"

image = (
    modal.Image.debian_slim(python_version="3.11")
    # chai_lab 0.6.1 declares `torch>=2.3.1,<2.7` (requirements.in at the v0.6.1 tag).
    # Install torch from the CUDA index FIRST so the chai install does not silently pull
    # a CPU wheel or drag the pin somewhere else.
    .pip_install(
        "torch==2.6.0",
        extra_options="--index-url https://download.pytorch.org/whl/cu124",
    )
    # Pinned on chai-lab's own advice: "API is quite stable, but we recommend pinning the
    # version in your requirements, i.e.: chai_lab==0.6.1". `main` has drifted past the
    # release (notably rdkit 2024.9 -> 2025.09) while still reporting __version__ 0.6.1,
    # and ligand atom NAMES are rdkit-assigned, so an unpinned install would silently move
    # the names any atom-level restraint refers to.
    .pip_install("chai_lab==0.6.1", "gemmi", "pandas", "pyarrow")
    .env({"CHAI_DOWNLOADS_DIR": "/cache/chai", "HF_HOME": "/cache/hf"})
)

app = modal.App(APP_NAME, image=image)

# Weights (6 component .pt files + the conformer pickle + traced ESM2-3B) and the converted
# MSA. Persisted so containers do not re-download ~4 GB each.
cache_vol = modal.Volume.from_name("cyp-chai-cache", create_if_missing=True)
# Predicted structures. SHARED with the Boltz and AF3 runners — the scorer reads
# /<tag>/<ligand>__<arm>__s<seed>/ out of this one volume.
out_vol = modal.Volume.from_name("cyp-pool", create_if_missing=True)


# ==========================================================================
# input construction (pure, importable, testable on CPU)
# ==========================================================================


def build_fasta(sequence: str, smiles: str) -> str:
    """Chai-1 entity-prefixed FASTA: protein (A), heme (B), ligand (C).

    The parser (`read_inputs`, inference_dataset.py) accepts exactly five prefixes —
    protein, ligand, rna, dna, glycan — matched case-insensitively on the text before the
    first `|`, and allows exactly one `|`-part after it. Anything else raises.

    The heme is NOT optional. It is the catalytic cofactor and the CYP3A4 pocket is
    literally the space above its distal face; folding without it asks the model to invent
    a cavity that has nothing in it.
    """
    return (
        f">protein|name=cyp3a4\n{sequence}\n"
        f">ligand|name=hem\n{HEM_SMILES}\n"
        f">ligand|name=lig\n{smiles}\n"
    )


RESTRAINT_HEADER = ("chainA,res_idxA,chainB,res_idxB,connection_type,confidence,"
                    "min_distance_angstrom,max_distance_angstrom,comment,restraint_id")


def build_restraints(sequence: str, arm: str, *, pocket_residues: list[int],
                     donor_atom_name: str | None, max_distance: float = 8.0) -> str | None:
    """Restraint CSV for a steered arm, or None when the arm needs no file.

    `res_idx` is the one-letter residue code followed by its 1-based index, e.g. `I301`,
    and Chai validates it against the sequence — so a wrong letter fails loudly, which is
    the one mercy in a file format whose other failure mode is silent.

    `pocket` rows: `chainA` is the BINDER chain with an EMPTY `res_idxA` (Chai's own
    example is `B,,A,C387,pocket,...`), one row per pocket residue.

    `covsteer` rows: a `covalent` restraint, the only connection type whose `@ATOM` suffix
    survives parsing. `@FE` on the heme chain is the CCD iron name; the ligand side uses
    the name RDKit assigns, which is version-sensitive — see `chai_probe` before trusting it.
    Chai's docs explicitly say intra/inter-ligand bonds are unevaluated. Treat any result
    from this arm as a spike, not evidence.
    """
    rows: list[str] = []
    if arm == "pocket":
        for i, r in enumerate(pocket_residues):
            if not (1 <= r <= len(sequence)):
                continue
            aa = sequence[r - 1]
            rows.append(f"{CHAIN_LIGAND},,{CHAIN_PROTEIN},{aa}{r},pocket,1.0,0.0,"
                        f"{max_distance},cyp3a4-distal-pocket,pocket_{i}")
    elif arm == "covsteer":
        if not donor_atom_name:
            return None
        rows.append(f"{CHAIN_HEME},@FE,{CHAIN_LIGAND},@{donor_atom_name},covalent,"
                    f"1.0,0.0,0.0,fe-coordination-UNVERIFIED,fecoord_0")
    if not rows:
        return None
    return RESTRAINT_HEADER + "\n" + "\n".join(rows) + "\n"


def chai_atom_name(smiles: str, atom_idx: int) -> str | None:
    """Best available guess at the name chai-lab's RDKit pass gives a heavy atom.

    DELIBERATELY A GUESS, and labelled one. Chai's covalent-bond docs tell you to go and
    look: "you will need check how the specific version of rdkit that we use in chai-lab
    ... assigns atom names". We cannot replicate that reliably from outside — the Boltz
    equivalent in `cypstruct.chem.boltz_atom_name` needed a CPU probe against the engine's
    own parser to settle whether ranking happens before or after AddHs, and got a different
    answer than either naive rule.

    So this returns the element-plus-occurrence form and `chai_probe` exists to check it
    against the installed chai-lab. Only the `covsteer` arm consumes this, and that arm is
    opt-in precisely because of this uncertainty.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or atom_idx >= mol.GetNumAtoms():
        return None
    sym = mol.GetAtomWithIdx(int(atom_idx)).GetSymbol().upper()
    n = 1 + sum(1 for i in range(int(atom_idx))
                if mol.GetAtomWithIdx(i).GetSymbol().upper() == sym)
    return f"{sym}{n}"


# The script the container actually runs. Kept as a string and written to disk so the call
# happens in a CHILD process with its own timeout: an in-process `run_inference` cannot be
# interrupted, and an uninterruptible child is how an A100 gets held for hours in silence.
_RUNNER_SRC = '''
import json, sys
from pathlib import Path
from chai_lab.chai1 import run_inference

cfg = json.loads(Path(sys.argv[1]).read_text())
msa_dir = Path(cfg["msa_dir"]) if cfg.get("msa_dir") else None
constraint = Path(cfg["constraint_path"]) if cfg.get("constraint_path") else None

# Every argument after `fasta_file` is keyword-only in chai_lab.chai1.run_inference.
# `use_msa_server` and `msa_directory` are mutually exclusive - run_inference asserts it.
run_inference(
    fasta_file=Path(cfg["fasta"]),
    output_dir=Path(cfg["out"]),
    msa_directory=msa_dir,
    use_msa_server=False,
    constraint_path=constraint,
    use_esm_embeddings=True,
    num_trunk_recycles=int(cfg.get("recycles", 3)),
    num_diffn_timesteps=int(cfg.get("timesteps", 200)),
    num_diffn_samples=int(cfg.get("samples", 5)),
    num_trunk_samples=int(cfg.get("trunk_samples", 1)),
    seed=int(cfg["seed"]) if cfg.get("seed") is not None else None,
    device="cuda:0",
    low_memory=True,
)
print("CHAI_RUN_OK")
'''

# chai-lab logs this and CONTINUES when a restraint cannot be built. Grepping for it is the
# only way to tell a steered run from a steered-looking run.
_DROP_MARKER = "generating distance constraints"


# ==========================================================================
# remote
# ==========================================================================


@app.function(volumes={"/cache": cache_vol}, timeout=3600, retries=1, max_containers=1,
              cpu=4.0)
def stage_msa(a3m_text: str) -> dict:
    """Convert the shared CYP3A4 a3m to the `.aligned.pqt` Chai actually reads. Once.

    The target is ONE sequence. Calling the ColabFold MSA server per ligand would fire the
    same query 168 times, get us rate-limited, and produce 168 identical alignments — so
    the alignment computed for the Boltz campaign (`data/reference/cyp3a4.a3m`, 6,979
    sequences) is reused verbatim here.

    Chai does not read a3m at inference time. It looks for `<sha256(seq.upper())>.aligned.pqt`
    in `msa_directory`, one file per unique chain sequence, with columns
    `sequence, source_database, pairing_key, comment` and a FIRST ROW whose
    `source_database` is literally `query` (asserted in `parse_aligned_pqt_to_msa_context`).
    `merge_a3m_in_directory` builds exactly that, including the hashed filename, so we let
    chai's own converter do it rather than hand-rolling the schema.

    The a3m is named `uniref90_hits.a3m` because the converter infers the source database
    from the filename stem and silently defaults unrecognised names to UNIREF90 anyway.
    Our alignment is actually a ColabFold uniref+envdb mixture, so the label is an
    approximation — `source_database` is a featurized column, not bookkeeping, and this is
    the honest place where this runner is less than exact. It is also why we do not bother
    with pairing keys: one protein chain has nothing to pair against.
    """
    import traceback

    out_dir = Path("/cache/msa_pqt")
    out_dir.mkdir(parents=True, exist_ok=True)
    existing = sorted(p.name for p in out_dir.glob("*.aligned.pqt"))
    if existing:
        cache_vol.commit()
        return {"ok": True, "cached": True, "files": existing}

    work = Path("/tmp/a3m_in")
    work.mkdir(parents=True, exist_ok=True)
    # Strip NUL bytes defensively. A single trailing NUL in an MSA server export once
    # killed a 4 MB, 6,979-sequence alignment inside Boltz's a3m parser AFTER the GPU had
    # been allocated. Cheap insurance; see `_clean_a3m` in modal_boltz.py.
    cleaned = a3m_text.replace("\x00", "").rstrip() + "\n"
    (work / "uniref90_hits.a3m").write_text(cleaned)

    try:
        from chai_lab.data.parsing.msas.aligned_pqt import merge_a3m_in_directory

        merge_a3m_in_directory(work, output_directory=out_dir)
    except Exception:
        return {"ok": False, "error": traceback.format_exc()[-2500:]}

    files = sorted(p.name for p in out_dir.glob("*.aligned.pqt"))
    cache_vol.commit()
    return {"ok": bool(files), "cached": False, "files": files,
            "a3m_sequences": cleaned.count(">")}


@app.function(volumes={"/cache": cache_vol}, timeout=5400, retries=2, max_containers=1,
              cpu=4.0)
def prefetch_weights() -> dict:
    """Download Chai-1's weights into the cache volume, on CPU, with retry. Once.

    The smoke test's gate job spent **960 seconds of A100 time** downloading weights and
    then died on a truncated transfer:
    `IncompleteRead(3564830720 bytes read, 2114109746 more expected)` - roughly 3.5 GB of
    a ~5.6 GB payload. Two things were wrong with that. A multi-gigabyte download does not
    belong on a GPU container, where it bills at GPU rates to do no compute. And a partial
    file left in the cache poisons every later run, because the downloader may see a file
    present and not re-fetch it.

    So: fetch on CPU, with retries, and clear anything implausibly small before trying, so
    a new attempt starts clean rather than inheriting a corrupt cache.
    """
    import shutil
    import traceback

    cache = Path("/cache/chai")
    cache.mkdir(parents=True, exist_ok=True)

    def inventory() -> dict:
        return {str(f.relative_to(cache)): f.stat().st_size
                for f in cache.rglob("*") if f.is_file()}

    before = inventory()
    removed = []
    for f in list(cache.rglob("*")):
        if (f.is_file() and f.suffix in (".pt", ".pt2", ".safetensors")
                and f.stat().st_size < 1_000_000):
            removed.append(str(f.relative_to(cache)))
            f.unlink(missing_ok=True)

    err = ""
    try:
        import chai_lab.utils.paths as _paths

        if hasattr(_paths, "download_models_if_needed"):
            _paths.download_models_if_needed()
        else:
            for comp in ("feature_embedding.pt", "token_embedder.pt", "trunk.pt",
                         "diffusion_module.pt", "confidence_head.pt"):
                getter = getattr(_paths, "chai1_component", None)
                if getter is None:
                    break
                try:
                    getter(comp)
                except Exception:
                    pass
    except Exception:
        err = traceback.format_exc()[-2500:]

    after = inventory()
    total = sum(after.values())
    cache_vol.commit()
    return {"ok": total > 1_000_000_000 and not err,
            "total_bytes": total, "n_files": len(after),
            "new_files": sorted(set(after) - set(before))[:20],
            "removed_truncated": removed,
            "cache_free_bytes": shutil.disk_usage("/cache").free,
            "error": err}



@app.function(volumes={"/cache": cache_vol}, timeout=1800, retries=1, max_containers=1,
              cpu=2.0)
def chai_probe(specs: list[dict]) -> dict:
    """CPU dry-run: parse every FASTA and restraint file with chai-lab's OWN parsers.

    This is the same gate `scripts/cofold/preflight_parse.py` applies to Boltz, and it
    exists for the same reason: a schema error discovered on the GPU path costs a batch,
    and here it may not even announce itself — Chai logs a restraint failure and keeps
    going. Seconds of CPU against hours of A100.

    It answers three questions that cannot be answered by reasoning:
      1. does the entity-prefixed FASTA parse, and does the heme SMILES survive RDKit;
      2. does the restraint CSV parse and validate its residue codes against the sequence;
      3. what does the INSTALLED rdkit actually name the ligand heavy atoms, so the
         `covsteer` arm's `@ATOM` reference can be checked rather than hoped for.

    Import paths are probed rather than assumed, and reported, so this keeps working (or
    fails loudly with a name) when chai-lab is upgraded.
    """
    import traceback

    out: dict = {"_imports": {}}

    read_inputs = None
    for mod, attr in (("chai_lab.data.dataset.inference_dataset", "read_inputs"),):
        try:
            read_inputs = getattr(__import__(mod, fromlist=[attr]), attr)
            out["_imports"][f"{mod}.{attr}"] = "ok"
        except Exception as exc:
            out["_imports"][f"{mod}.{attr}"] = f"{type(exc).__name__}: {exc}"[:200]

    parse_restraints = None
    for mod, attr in (("chai_lab.data.parsing.restraints", "parse_pairwise_table"),
                      ("chai_lab.data.parsing.restraints", "parse_restraints"),
                      ("chai_lab.data.parsing.restraints", "_parse_res_idx")):
        try:
            parse_restraints = getattr(__import__(mod, fromlist=[attr]), attr)
            out["_imports"][f"{mod}.{attr}"] = "ok"
            break
        except Exception as exc:
            out["_imports"][f"{mod}.{attr}"] = f"{type(exc).__name__}: {exc}"[:200]

    try:
        from rdkit import Chem

        out["_rdkit"] = Chem.rdBase.rdkitVersion
        out["_hem_heavy_atoms"] = (
            Chem.MolFromSmiles(HEM_SMILES).GetNumAtoms()
            if Chem.MolFromSmiles(HEM_SMILES) else None)
    except Exception:
        out["_rdkit"] = traceback.format_exc()[-300:]

    work = Path("/tmp/chai_probe")
    work.mkdir(parents=True, exist_ok=True)
    for spec in specs:
        jid = spec["job_id"]
        rec: dict = {}
        f = work / f"{jid}.fasta"
        f.write_text(spec["fasta"])
        if read_inputs is not None:
            try:
                inputs = read_inputs(f)
                rec["n_entities"] = len(inputs)
                rec["entity_types"] = [str(getattr(i, "entity_type", "?")) for i in inputs]
                rec["fasta_ok"] = True
            except Exception as exc:
                rec["fasta_ok"] = False
                rec["fasta_error"] = f"{type(exc).__name__}: {exc}"[:300]
        if spec.get("restraints"):
            c = work / f"{jid}.restraints"
            c.write_text(spec["restraints"])
            rec["n_restraint_rows"] = spec["restraints"].count("\n") - 1
            if parse_restraints is not None:
                try:
                    parse_restraints(c)
                    rec["restraints_ok"] = True
                except Exception as exc:
                    rec["restraints_ok"] = False
                    rec["restraints_error"] = f"{type(exc).__name__}: {exc}"[:300]
        out[jid] = rec
    return out


@app.function(
    gpu=GPU,
    volumes={"/cache": cache_vol, "/out": out_vol},
    timeout=FN_TIMEOUT,
    retries=1,                 # one retry, not Modal's silent-forever default posture
    max_containers=MAX_CONTAINERS,
)
def cofold(batch: list[dict]) -> list[dict]:
    """Run a chunk of (ligand x arm x seed) jobs. Commits after each so preemption is cheap.

    Idempotent per job: an existing `DONE.json` short-circuits before any GPU work.
    """
    import subprocess
    import traceback

    results: list[dict] = []
    consecutive_failures = 0
    msa_dir = Path("/cache/msa_pqt")

    for spec in batch:
        job_id = spec["job_id"]
        dest = Path("/out") / spec["tag"] / job_id
        done_marker = dest / "DONE.json"
        if done_marker.exists():
            try:
                rec = json.loads(done_marker.read_text())
            except json.JSONDecodeError:
                rec = {}
            results.append({**rec, "job_id": job_id, "status": "cached"})
            continue

        t0 = time.time()
        work = Path("/tmp") / job_id
        # chai-lab's example notes "Inference expects an empty directory". A retry into a
        # populated directory is a confusing failure, so start clean every time.
        if work.exists():
            import shutil

            shutil.rmtree(work, ignore_errors=True)
        (work / "out").mkdir(parents=True, exist_ok=True)

        fasta = work / "input.fasta"
        fasta.write_text(spec["fasta"])
        cpath = None
        if spec.get("restraints"):
            cpath = work / "input.restraints"
            cpath.write_text(spec["restraints"])

        cfg = {
            "fasta": str(fasta),
            "out": str(work / "out"),
            "msa_dir": str(msa_dir) if any(msa_dir.glob("*.aligned.pqt")) else None,
            "constraint_path": str(cpath) if cpath else None,
            "samples": spec.get("samples", 5),
            "recycles": spec.get("recycles", 3),
            "timesteps": spec.get("timesteps", 200),
            "trunk_samples": spec.get("trunk_samples", 1),
            "seed": spec.get("seed"),
        }
        (work / "cfg.json").write_text(json.dumps(cfg))
        (work / "run_chai.py").write_text(_RUNNER_SRC)

        try:
            cp = subprocess.run(
                [sys.executable, str(work / "run_chai.py"), str(work / "cfg.json")],
                capture_output=True, text=True, timeout=INNER_TIMEOUT)
            rc = cp.returncode
            err = cp.stderr[-4000:]
            combined = (cp.stdout or "") + (cp.stderr or "")
        except subprocess.TimeoutExpired:
            rc, err, combined = -9, f"chai run_inference exceeded {INNER_TIMEOUT}s", ""
        except Exception:
            rc, err, combined = -1, traceback.format_exc()[-4000:], ""

        # Harvest. Chai writes pred.model_idx_N.cif and scores.model_idx_N.npz.
        dest.mkdir(parents=True, exist_ok=True)
        n_struct = 0
        for p in sorted((work / "out").rglob("*")):
            if not p.is_file():
                continue
            if p.suffix == ".cif":
                (dest / p.name).write_bytes(p.read_bytes())
                n_struct += 1
                side = _confidence_sidecar(p)
                if side is not None:
                    (dest / f"confidence_{p.stem}.json").write_text(
                        json.dumps(side, indent=1))
            elif p.suffix in (".npz", ".json"):
                (dest / p.name).write_bytes(p.read_bytes())

        # A dropped restraint is NOT a crash. Detect it or the steered arm lies.
        dropped = (spec.get("arm") not in (None, "unsteered")
                   and _DROP_MARKER in combined)
        ok = rc == 0 and n_struct > 0
        status = "ok" if ok else "failed"
        if ok and dropped:
            status = "restraints_dropped"

        rec = {"job_id": job_id, "returncode": rc, "n_structures": n_struct,
               "seconds": round(time.time() - t0, 1), "engine": "chai1",
               "arm": spec.get("arm"), "ligand": spec.get("ligand"),
               "seed": spec.get("seed"), "status": status,
               "restraints_dropped": bool(dropped),
               "stderr_tail": "" if ok else err}
        (dest / "DONE.json").write_text(json.dumps(rec, indent=1))
        out_vol.commit()
        results.append(rec)

        consecutive_failures = 0 if ok else consecutive_failures + 1
        if consecutive_failures >= MAX_CONSECUTIVE_FAILURES:
            # Every input failing the same way means the configuration is wrong, not the
            # ligands. Continuing would pay full wall-clock to produce nothing.
            results.append({"status": "aborted", "job_id": job_id,
                            "reason": f"{consecutive_failures} consecutive failures; "
                                      "circuit breaker tripped"})
            break

    return results


def _confidence_sidecar(cif: Path) -> dict | None:
    """Turn Chai's `scores.model_idx_N.npz` into the JSON the pool scorer reads.

    `collect_and_score.py` looks for `confidence_<cif stem>.json` and reads, in order,
    `complex_ipde` (negated, because it is an ERROR), then `confidence_score`, `ptm`,
    `iptm`, `complex_plddt`. Writing that file here is what lets a Chai tag be scored by
    the same code path as a Boltz tag instead of needing a parallel scorer.

    Keys are read defensively rather than assumed: the npz is whatever `SampleRanking`
    serialises in the installed version, so we take every 0-d/1-element numeric array we
    find and then pick a `confidence_score` from the first key that exists.
    """
    idx = cif.stem.split("model_idx_")[-1]
    npz = cif.parent / f"scores.model_idx_{idx}.npz"
    if not npz.exists():
        return None
    try:
        import numpy as np

        d = np.load(npz, allow_pickle=True)
        flat: dict = {}
        for k in d.files:
            v = d[k]
            try:
                if v.size == 1:
                    flat[k] = float(v.reshape(-1)[0])
            except Exception:
                continue
        if not flat:
            return None
        for k in ("aggregate_score", "confidence_score", "ptm", "iptm", "complex_plddt"):
            if k in flat:
                flat["confidence_score"] = flat[k]
                break
        return flat
    except Exception:
        return None


# ==========================================================================
# local driver
# ==========================================================================


def plan(csv_path: str, tag: str, samples: int, seeds: list[int],
         arms: tuple[str, ...] = ("unsteered",)) -> list[dict]:
    """Build the finite job list. No job is created without a resolvable job_id.

    `job_id` is `<ligand>__<arm>__s<seed>` because `collect_and_score.py` splits on `__`
    to recover the ligand and the arm. Do not put a double underscore in an arm name.
    """
    import csv as _csv

    from cypstruct.chem import coordinating_atoms, standardize
    from cypstruct.targets import CYP3A4_POCKET, fetch_sequences

    seq = fetch_sequences()["cyp3a4"]
    pocket = sorted(set(CYP3A4_POCKET["i_helix"] + CYP3A4_POCKET["phe_cluster"]
                        + CYP3A4_POCKET["polar"]))
    jobs: list[dict] = []
    with open(csv_path) as fh:
        for row in _csv.DictReader(fh):
            lid = row.get("id") or row.get("idx") or row.get("Molecule_Name")
            smi = standardize(row.get("smiles") or row.get("SMILES") or "")
            if not smi or not lid:
                continue
            donors = coordinating_atoms(smi, top_k=1)
            donor_name = chai_atom_name(smi, donors[0].atom_idx) if donors else None
            fasta = build_fasta(seq, smi)
            for arm in arms:
                restraints = build_restraints(seq, arm, pocket_residues=pocket,
                                              donor_atom_name=donor_name)
                if arm != "unsteered" and restraints is None:
                    continue      # nothing to steer with; the unsteered arm covers it
                for seed in seeds:
                    jobs.append(dict(
                        job_id=f"{lid}__{arm}__s{seed}", tag=tag, ligand=lid, arm=arm,
                        seed=seed, samples=samples, fasta=fasta, restraints=restraints,
                    ))
    return jobs


def _load_msa() -> str:
    msa_path = (REPO / "data" / "reference" / "cyp3a4.a3m") if REPO else None
    if msa_path is None or not msa_path.exists():
        raise SystemExit(
            f"missing {msa_path}. Chai reads a precomputed MSA (converted to .aligned.pqt) "
            "rather than calling the MSA server per ligand; the alignment is the one the "
            "Boltz campaign already computed. Warm it with "
            "`python scripts/cofold/modal_boltz.py` first, or copy it into data/reference/.")
    return msa_path.read_text()


def submit(csv_path: str, tag: str, samples: int = 5, seeds: tuple[int, ...] = (1,),
           arms: tuple[str, ...] = ("unsteered",), dry_run: bool = False,
           probe_only: bool = False) -> None:
    from cypstruct import budget

    jobs = plan(csv_path, tag, samples, list(seeds), arms=arms)
    if not jobs:
        raise SystemExit("no jobs planned - check the csv columns (id, smiles)")
    est = budget.estimate("chai_cofold", len(jobs), samples)
    ok, why, _ = budget.preflight_hours("chai_cofold", est)
    print(f"planned jobs: {len(jobs)}  ({len(jobs) * samples} ligand-samples)  "
          f"arms={arms}  est {est:.2f} GPU-h")
    if not ok:
        raise SystemExit(f"PREFLIGHT REFUSED: {why}")
    if dry_run:
        j = dict(jobs[0])
        print(json.dumps({k: (v[:600] if isinstance(v, str) else v)
                          for k, v in j.items()}, indent=1))
        return

    if probe_only:
        with app.run():
            res = chai_probe.remote(jobs[:40])
        imports = res.pop("_imports", {})
        rdk = res.pop("_rdkit", None)
        hem = res.pop("_hem_heavy_atoms", None)
        bad = {k: v for k, v in res.items()
               if v.get("fasta_ok") is False or v.get("restraints_ok") is False}
        print(f"probe: {len(res) - len(bad)}/{len(res)} clean   rdkit={rdk}  "
              f"hem_heavy_atoms={hem}")
        print("imports:", json.dumps(imports, indent=1))
        for k, v in sorted(bad.items())[:20]:
            print(f"  FAIL {k}: {json.dumps(v)[:200]}")
        raise SystemExit(1 if bad else 0)

    run_id = f"{tag}-chai-{int(time.time())}"
    budget.record(run_id, "modal", "chai_cofold", len(jobs) * samples, est,
                  note=f"app={APP_NAME} tag={tag} arms={','.join(arms)}",
                  n_jobs=len(jobs), samples=samples, app=APP_NAME)
    n_ok = n_fail = n_dropped = 0
    first_errors: list[str] = []
    try:
        with app.run():
            print("staging the shared MSA as .aligned.pqt ...", flush=True)
            w = stage_msa.remote(_load_msa())
            print("  ", json.dumps(w)[:500], flush=True)
            if not w.get("ok"):
                raise SystemExit("MSA staging failed; not launching the array")

            # Pull the weights on CPU BEFORE any GPU is allocated. The first smoke test
            # burned 960 s of A100 time downloading and then died on a truncated transfer
            # (3.5 GB of a ~5.6 GB payload). A multi-gigabyte download billed at GPU rates
            # for doing no compute is pure waste, and the partial file it leaves behind
            # poisons the cache for every later run.
            print("prefetching Chai weights on CPU ...", flush=True)
            pw = prefetch_weights.remote()
            print("  ", json.dumps({k: v for k, v in pw.items()
                                    if k != "error"})[:400], flush=True)
            if not pw.get("ok"):
                print("   error tail:", (pw.get("error") or "")[-600:], flush=True)
                raise SystemExit("weight prefetch failed; not launching the array")

            # Then run ONE job alone. With the weights already cached this is a real
            # schema/environment check rather than a download, and it turns an error into
            # one wasted job instead of a whole batch.
            print("gating on a single job before mapping the rest ...", flush=True)
            gate = cofold.remote([jobs[0]])
            print("  ", json.dumps(gate)[:900], flush=True)
            if not any(r.get("status") in ("ok", "cached", "restraints_dropped")
                       for r in gate):
                raise SystemExit("gate job produced no structures; not launching the array")
            n_ok += 1

            rest = jobs[1:]
            chunks = [rest[i:i + CHUNK] for i in range(0, len(rest), CHUNK)]
            for res in cofold.map(chunks, order_outputs=False, return_exceptions=True):
                if isinstance(res, Exception):
                    n_fail += 1
                    # Print the first few. Counting exceptions and printing nothing is how
                    # a run that failed every input reports "0 ok, N failed" with no way
                    # to tell why. A failure count without a reason is not a diagnostic.
                    if len(first_errors) < 5:
                        first_errors.append(f"{type(res).__name__}: {res}"[:300])
                        print(f"  [exception] {first_errors[-1]}", flush=True)
                    continue
                for r in res:
                    st = r.get("status")
                    if st in ("ok", "cached"):
                        n_ok += 1
                    elif st == "restraints_dropped":
                        n_ok += 1
                        n_dropped += 1
                    elif st == "failed":
                        n_fail += 1
                        if len(first_errors) < 5:
                            first_errors.append(r.get("stderr_tail", "")[-300:])
                            print(f"  [failed] {r.get('job_id')}: {first_errors[-1]}",
                                  flush=True)
                    elif st == "aborted":
                        print(f"  !! circuit breaker: {r.get('reason')}", flush=True)
                print(f"  {n_ok} ok / {n_fail} failed of {len(jobs)}", flush=True)
        note = f"{n_ok} ok, {n_fail} failed"
        if n_dropped:
            note += f", {n_dropped} with SILENTLY DROPPED restraints"
        budget.close(run_id, "done", note=note)
        print("DONE:", note)
        if n_dropped:
            print("!! Those jobs produced structures but chai discarded the restraints. "
                  "They are UNSTEERED poses wearing a steered label - do not A/B them.")
    except BaseException as exc:
        budget.close(run_id, "failed", note=str(exc)[:300])
        raise


def collect(tag: str) -> None:
    """Pull results out of the Modal volume and push them to OneDrive.

    Usually you do NOT want this: `scripts/structure/collect_and_score.py --tag <tag>`
    streams one job at a time, scores it and discards it, with a peak footprint of ~10 MB.
    This box has repeatedly sat near zero free bytes on both C: and D:, so a bulk download
    is the exception, for when the poses themselves need to be archived.
    """
    from cypstruct.storage import Batch

    with Batch(f"pool/chai1/{tag}") as b:
        vol = modal.Volume.from_name("cyp-pool")
        n = 0
        for entry in vol.iterdir(f"/{tag}"):
            job = entry.path.rstrip("/").split("/")[-1]
            for f in vol.iterdir(f"/{tag}/{job}"):
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

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=["submit", "plan", "probe", "collect"])
    ap.add_argument("--csv")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--samples", type=int, default=5)
    ap.add_argument("--seeds", default="1")
    ap.add_argument("--arms", default="unsteered",
                    help="comma-separated: unsteered | pocket | covsteer. "
                         "pocket and covsteer are opt-in; read the module docstring for "
                         "why the default is unsteered only.")
    a = ap.parse_args()

    if a.cmd == "collect":
        collect(a.tag)
    else:
        seeds = tuple(int(s) for s in a.seeds.split(","))
        arms = tuple(s.strip() for s in a.arms.split(",") if s.strip())
        unknown = set(arms) - {"unsteered", "pocket", "covsteer"}
        if unknown:
            raise SystemExit(f"unknown arm(s): {sorted(unknown)}")
        submit(a.csv, a.tag, a.samples, seeds, arms=arms,
               dry_run=(a.cmd == "plan"), probe_only=(a.cmd == "probe"))
