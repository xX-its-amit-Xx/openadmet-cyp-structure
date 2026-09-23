"""Boltz-2 co-folding on Explorer (SLURM) at conditioning MATCHED to the Modal pool.

This is the generation path that exists because the other two are gone: Modal is over its
spend cap, and `FINDING_034` measured that OpenProtein's `boltz2` **rejects an uploaded
MSA** (`JobStatus.FAILURE, "internal server error"`), so the only depth it can sell is
single-sequence -- poses whose ceiling matches the existing pool and whose average is
0.057 LDDT-PLI lower. Mixed in, they cost 0.0065 of selected score. Depth bought from a
degraded configuration is not depth.

Explorer can sell matched depth, and almost all of the work was already done by the
fine-tuning campaign, which is why this module is thin:

  * `/scratch/shenoy.am/cyp-finetune/env/bin/boltz` -- **boltz 2.2.1**, the exact version
    Modal pinned (`boltz[cuda]==2.2.1`), on torch 2.5.1.
  * `boltz_cache/boltz2_conf.ckpt` -- pulled by `boltz_cache/get.sh` from
    `huggingface.co/boltz-community/boltz-2`, i.e. the same artefact `boltz predict`
    downloads for itself. Weights, CCD and `mols.tar` are all local.
  * `submit_holdout.sh` -- a working `sbatch` template, already exercised: five COMPLETED
    inference jobs at ~43 minutes each.

**The structural hazard, and why nothing here touches the network.** Explorer's GPU nodes
have no direct internet (proxy `http://10.99.0.130:3128`), which is incompatible with
`--use_msa_server`; the PXR campaign lost ~12 days to exactly this. Every input this
module needs -- the alignment, the weights, the CCD -- is staged to `/scratch` from a
login node first, and `--cache` points the job at the local copy so no code path reaches
for a download and hangs.

**The alignment is the one that matters.** `data/reference/cyp3a4.a3m` is the *same file*
the Modal pool was folded against -- 6,979 sequences, md5 `6de0ee13...`. It is staged
verbatim and its md5 is re-verified on the far side **inside a Slurm job**, never on the
login node, whose memory cgroup SIGKILLs long reads and lets a truncated read report
success.

**Conditioning is not re-implemented, it is imported.** `build_yaml` is lifted out of
`scripts/cofold/modal_boltz.py` by source extraction rather than copied, so the YAML this
writes is produced by the very function that produced the existing pool. Only two things
differ by construction, and both are the point of the exercise: the `msa:` path, and the
seed.

Usage (each stage is resumable and prints counts):

    python scripts/explorer/boltz_depth.py plan                 # YAMLs + manifest, local
    python scripts/explorer/boltz_depth.py stage                # rsync to /scratch
    python scripts/explorer/boltz_depth.py verify               # md5 of the staged MSA, IN A JOB
    python scripts/explorer/boltz_depth.py submit --seed 101    # sbatch, returns job ids
    python scripts/explorer/boltz_depth.py poll                 # ADVANCING progress, not liveness
    python scripts/explorer/boltz_depth.py collect              # rsync mmCIFs back
    python scripts/explorer/boltz_depth.py score                # LDDT-PLI / BiSyRMSD / xeng

**A BLIND set: `--csv` and `--tag`.** Without them `plan` reads `validation_ligands.csv`
and filters to the 14-ligand predicted-Type-I stratum, which is the FINDING 035
pre-registration and not a property of the generator -- on a test set it would silently
write inputs for the wrong ligands. `--csv` takes every id in the file; `--tag` names the
subdirectory every stage keys off, on both sides, and defaults to the historical
`stratum`/`smoke` so an existing campaign is untouched. `plan` writes
`plans/<tag>.json`, which is where `collect` gets its expected ligand list and `score`
learns the set has no crystals and refuses:

    python scripts/explorer/boltz_depth.py plan --csv data/processed/test_ligands.csv \\
        --tag drop
    python scripts/explorer/boltz_depth.py stage   --tag drop
    python scripts/explorer/boltz_depth.py submit  --tag drop --seed 101 --samples 20
    python scripts/explorer/boltz_depth.py collect --tag drop      # then STOP: no `score`
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "cofold"))

from cypstruct.paths import DATA_PROCESSED, REFERENCE, SCRATCH  # noqa: E402

# --- fixed by the pre-registration; do not tune -----------------------------
COORD_MAX = 2.6                 # repo constant (build_reference_set.py)
ARM = "explorer_boltz2_msa"     # the new arm's label, distinct from the pool's "unsteered"
SAMPLES = 20                    # diffusion_samples, matching the existing pool exactly
RECYCLING = 3
SAMPLING_STEPS = 200
AXIAL_CYS = 442

REMOTE_ROOT = "/scratch/shenoy.am/cyp-depth"
FT_ROOT = "/scratch/shenoy.am/cyp-finetune"      # env, weights, CCD -- reused, not rebuilt
SSH = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=30", "explorer"]

LOCAL = SCRATCH / "matched_depth"                # C:\cyp_struct\... -- never D:
YAML_DIR = LOCAL / "yaml"
POSE_DIR = LOCAL / "poses"
MSA_SRC = REFERENCE / "cyp3a4.a3m"
JOBS_JSON = DATA_PROCESSED / "matched_depth_jobs.json"
POSES_CSV = DATA_PROCESSED / "matched_depth_poses.csv"


def _run(cmd: list[str], timeout: int = 600) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def _ssh(script: str, timeout: int = 600) -> subprocess.CompletedProcess:
    return _run(SSH + [script], timeout=timeout)


def _win_to_posix(p: Path) -> str:
    """rsync/tar read `D:\\x` as host:path. Under Git Bash the drive must become a POSIX
    root or the colon is parsed as a hostname separator (stage_msas.py hit this first)."""
    s = str(p)
    if len(s) > 2 and s[1] == ":":
        s = "/" + s[0].lower() + s[2:].replace("\\", "/")
    return s.replace("\\", "/")


def _git_bash() -> str:
    """Absolute path to GIT Bash.

    Never the bare name `bash`: on this box `bash` on PATH resolves to
    `C:\\Windows\\System32\\bash.exe`, the WSL launcher, which answers every command with
    "Windows Subsystem for Linux has no installed distributions" on **stdout** and exits
    1. That is the worst possible failure shape -- a real exit code with plausible-looking
    output and an empty stderr -- and it silently ate the first three staging attempts
    here, which reported `ok: false` with no error at all.
    """
    for c in (r"C:\Program Files\Git\bin\bash.exe",
              r"C:\Program Files\Git\usr\bin\bash.exe",
              r"C:\Program Files (x86)\Git\bin\bash.exe"):
        if Path(c).exists():
            return c
    raise RuntimeError("Git Bash not found; refusing to fall back to WSL's bash stub")


def _bash(script: str, timeout: int = 3600) -> subprocess.CompletedProcess:
    """Run a POSIX snippet through Git Bash.

    Needed because `tar | ssh` is a pipeline and Windows `cmd` cannot express one.
    """
    return subprocess.run([_git_bash(), "-c", script],
                          capture_output=True, text=True, timeout=timeout)


def _push_dir(local: Path, remote: str, timeout: int = 3600) -> dict:
    """Copy a directory's contents to Explorer with tar-over-ssh.

    **rsync does not work from this box at all.** Every invocation -- through subprocess
    and typed straight into the shell alike -- dies with `dup() in/out/err failed` before
    a single byte moves, an MSYS2 rsync/OpenSSH stdio incompatibility rather than anything
    about the cluster. `scp` and `tar | ssh` both work, so the transport is tar: it
    preserves the directory, moves many small files in one connection, and needs no
    remote rsync at all.

    The trap this walks past is in the return code, not the transport: the first thing
    tried here printed `rsync rc=0` from an `echo $?` that had captured `tail`'s status,
    not rsync's -- a failed copy reporting success. Hence PIPESTATUS below, and hence the
    far-side count that follows every push.
    """
    lp = _win_to_posix(local)
    cmd = (f"tar czf - -C '{lp}' . | ssh -o BatchMode=yes explorer "
           f"'mkdir -p {remote} && tar xzf - -C {remote}'; "
           f"echo TAR_RC=${{PIPESTATUS[0]}} SSH_RC=${{PIPESTATUS[1]}}")
    r = _bash(cmd, timeout=timeout)
    rc = {k: v for k, v in
          (t.split("=") for t in r.stdout.split() if "_RC=" in t)}
    return {"tar_rc": rc.get("TAR_RC"), "ssh_rc": rc.get("SSH_RC"),
            "ok": rc.get("TAR_RC") == "0" and rc.get("SSH_RC") == "0",
            "stderr": r.stderr[-300:]}


def _push_file(local: Path, remote: str, timeout: int = 1800) -> dict:
    r = _bash(f"scp -o BatchMode=yes '{_win_to_posix(local)}' explorer:{remote}; "
              f"echo SCP_RC=$?", timeout=timeout)
    rc = [t for t in r.stdout.split() if t.startswith("SCP_RC=")]
    ok = bool(rc) and rc[0] == "SCP_RC=0"
    return {"ok": ok, "rc": rc[0] if rc else None, "stderr": r.stderr[-300:]}


# ---------------------------------------------------------------------------
# conditioning: imported from the runner that built the pool, never re-written
# ---------------------------------------------------------------------------

def _import_build_yaml():
    """Lift `build_yaml` out of modal_boltz.py without importing modal.

    `modal` is not installed locally and modal_boltz builds an App at module scope, so a
    plain import fails. Extracting the function's AST node and exec'ing that gives the
    byte-identical function body -- which is the whole point: the new arm's YAML must come
    from the same code that produced the 20 poses it will be compared against, not from a
    copy that has drifted.
    """
    src = (REPO / "scripts" / "cofold" / "modal_boltz.py").read_text()
    tree = ast.parse(src)
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "build_yaml":
            ns: dict = {}
            exec(compile(ast.Module(body=[node], type_ignores=[]),
                         "modal_boltz.build_yaml", "exec"), ns)
            return ns["build_yaml"], ast.get_source_segment(src, node)
    raise RuntimeError("build_yaml not found in modal_boltz.py")


def type_i_ligands() -> pd.DataFrame:
    """The predicted-Type-I stratum. Prediction-side rule, nothing from a crystal."""
    L = pd.read_csv(DATA_PROCESSED / "binding_mode_labels_cyp3a4.csv")
    sel = L[L.pred_fe_donor_median > COORD_MAX].copy()
    return sel.sort_values("id").reset_index(drop=True)


def smoke_ligand() -> str:
    """The smoke-test ligand, chosen a priori and NOT from the stratum.

    Fixed rule: the first predicted-Type-II ligand alphabetically. Exercising the path on
    a stratum ligand would generate poses before the pre-registration was committed.
    """
    L = pd.read_csv(DATA_PROCESSED / "binding_mode_labels_cyp3a4.csv")
    return sorted(L[L.pred_fe_donor_median <= COORD_MAX].id.tolist())[0]


def counts() -> dict:
    L = pd.read_csv(DATA_PROCESSED / "binding_mode_labels_cyp3a4.csv")
    sel = type_i_ligands()
    return {
        "validation_ligands": int(len(L)),
        "filter_pred_fe_donor_median_gt_2.6": int(len(sel)),
        "excluded_predicted_type_II": int((L.pred_fe_donor_median <= COORD_MAX).sum()),
        "stratum": sorted(sel.id.tolist()),
        "smoke_ligand_not_in_stratum": smoke_ligand(),
    }


# ---------------------------------------------------------------------------
# plan -- write the YAMLs locally
# ---------------------------------------------------------------------------

PLAN_DIR = LOCAL / "plans"


def _tag_for(smoke: bool, tag: str | None) -> str:
    """The subdirectory every stage keys off.

    Defaults to the two historical names so an existing campaign keeps running against
    exactly the directories it already wrote: `stratum` and `smoke`. `--tag` is additive.
    """
    if tag:
        if not tag.replace("_", "").replace("-", "").isalnum():
            raise SystemExit(f"--tag {tag!r} must be alphanumeric (plus - and _): "
                             "it becomes a directory name on both sides")
        return tag
    return "smoke" if smoke else "stratum"


def _save_plan(tag: str, rec: dict) -> None:
    PLAN_DIR.mkdir(parents=True, exist_ok=True)
    (PLAN_DIR / f"{tag}.json").write_text(json.dumps(rec, indent=1))


def _load_plan(tag: str) -> dict | None:
    f = PLAN_DIR / f"{tag}.json"
    return json.loads(f.read_text()) if f.exists() else None


def _ligand_table(csv: str | None) -> pd.DataFrame:
    """The ligand table for a tag, blind-set or validation-set alike.

    A blind CSV has `id,smiles` and no `pdb`. That is the whole difference, and it is why
    `score` cannot run on one (G2) while `plan`/`collect` can.
    """
    df = pd.read_csv(csv or (DATA_PROCESSED / "validation_ligands.csv"))
    idc = "structure" if "structure" in df.columns and "id" not in df.columns else "id"
    if idc not in df.columns:
        raise SystemExit(f"{csv}: needs an `id` (or `structure`) column, got "
                         f"{list(df.columns)}")
    if "smiles" not in df.columns:
        raise SystemExit(f"{csv}: needs a `smiles` column, got {list(df.columns)}")
    df = df.rename(columns={idc: "id"}).astype({"id": str})
    return df.drop_duplicates("id").reset_index(drop=True)


def plan(smoke: bool = False, csv: str | None = None,
         tag: str | None = None) -> dict:
    from cypstruct.chem import standardize
    from cypstruct.targets import fetch_sequences

    build_yaml, _ = _import_build_yaml()
    seq = fetch_sequences()["cyp3a4"]
    tbl = _ligand_table(csv)
    vl = tbl.set_index("id")

    # `--csv` takes EVERY id in the file. The Type-I filter is a property of the
    # FINDING 035 pre-registration, not of the generator, and applying it to a blind test
    # set would silently drop most of the submission. The RUNBOOK claimed "point plan at
    # the test-set ligand CSV; nothing else changes" long before this flag existed.
    if csv:
        if smoke:
            raise SystemExit("--smoke selects a validation ligand by rule; it is "
                             "meaningless with --csv. Use --tag to name a small run.")
        ids = sorted(tbl.id.tolist())
        filt = {"source_csv": csv, "rule": "all ids in the csv", "n": len(ids)}
    else:
        ids = [smoke_ligand()] if smoke else sorted(type_i_ligands().id.tolist())
        filt = {"source_csv": "validation_ligands.csv",
                "rule": "smoke ligand" if smoke else
                        "predicted Type-I stratum (pred_fe_donor_median > 2.6)",
                "n": len(ids)}
    sub = _tag_for(smoke, tag)
    out = YAML_DIR / sub
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob("*.yaml"):
        old.unlink()

    remote_msa = f"{REMOTE_ROOT}/msa/cyp3a4.a3m"
    written, skipped = [], []
    for lid in ids:
        if lid not in vl.index:
            skipped.append(lid)
            continue
        smi = standardize(str(vl.loc[lid, "smiles"]))
        if not smi:
            skipped.append(lid)
            continue
        # steer=False: the existing pool's `unsteered` arm. The heme bond is NOT steering
        # -- build_yaml emits it unconditionally, as the resting state of the enzyme.
        y = build_yaml(seq, smi, steer=False, donor_atom_name=None,
                       template_cif=None, axial_cys=AXIAL_CYS,
                       pocket_residues=None, msa_path=remote_msa)
        (out / f"{lid}.yaml").write_text(y)
        written.append(lid)

    rec = {"mode": sub, "tag": sub, "csv": csv, "filter": filt,
           "yaml_written": len(written), "ligands": written,
           "skipped": skipped, "dir": str(out),
           "msa_path_in_yaml": remote_msa,
           "blind": csv is not None and "pdb" not in tbl.columns,
           "counts": (counts() if not smoke and not csv else None)}
    _save_plan(sub, rec)
    return rec


def conditioning() -> dict:
    """Everything that has to match, with how it was checked. Printed into the PREREG."""
    build_yaml, src = _import_build_yaml()
    msa_bytes = MSA_SRC.read_bytes()
    seqs = msa_bytes.count(b">")
    local_md5 = hashlib.md5(msa_bytes).hexdigest()

    # the pool's own conditioning, read back out of the runner rather than remembered
    pool = pd.read_csv(DATA_PROCESSED / "poses_scored_val87b.csv")
    un = pool[pool.arm == "unsteered"]
    return {
        "msa": {"file": str(MSA_SRC), "md5": local_md5, "sequences": seqs,
                "nul_bytes": msa_bytes.count(b"\x00"),
                "bytes": len(msa_bytes)},
        "build_yaml_sha1": hashlib.sha1(src.encode()).hexdigest(),
        "axial_cys": AXIAL_CYS,
        "diffusion_samples": SAMPLES,
        "recycling_steps": RECYCLING,
        "sampling_steps": SAMPLING_STEPS,
        "existing_pool": {
            "rows": int(len(un)), "ligands": int(un.ligand.nunique()),
            "samples_per_ligand": sorted(set(un.groupby("ligand").size().tolist())),
        },
    }


# ---------------------------------------------------------------------------
# stage -- push inputs to /scratch from a LOGIN node (which has internet)
# ---------------------------------------------------------------------------

def stage(smoke: bool = False, tag: str | None = None) -> dict:
    sub = _tag_for(smoke, tag)
    mk = _ssh(f"mkdir -p {REMOTE_ROOT}/msa {REMOTE_ROOT}/yaml/{sub} "
              f"{REMOTE_ROOT}/out {REMOTE_ROOT}/logs && echo OK", timeout=120)
    if "OK" not in mk.stdout:
        return {"ok": False, "stage": "mkdir", "err": mk.stderr[-400:]}

    res = {"ok": True, "mkdir": True}
    # the alignment: 4 MB, pushed once and reused by every job
    res["msa"] = _push_file(MSA_SRC, f"{REMOTE_ROOT}/msa/")
    res["yaml"] = _push_dir(YAML_DIR / sub, f"{REMOTE_ROOT}/yaml/{sub}")
    res["ok"] = res["msa"]["ok"] and res["yaml"]["ok"]

    # never trust the transport's own word: count and hash on the far side
    chk = _ssh(f"echo YAML=$(ls {REMOTE_ROOT}/yaml/{sub}/*.yaml 2>/dev/null | wc -l); "
               f"md5sum {REMOTE_ROOT}/msa/cyp3a4.a3m; "
               f"echo SEQS=$(grep -c '^>' {REMOTE_ROOT}/msa/cyp3a4.a3m)", timeout=300)
    res["remote_check"] = chk.stdout.strip().splitlines()
    res["msa_md5_matches_local"] = (
        hashlib.md5(MSA_SRC.read_bytes()).hexdigest() in chk.stdout)
    return res


def push_sbatch() -> dict:
    """Copy the sbatch templates up. Kept as repo files so they are reviewable/versioned."""
    out = {}
    for name in ("boltz_depth.sbatch", "verify_inputs.sbatch"):
        out[name] = _push_file(REPO / "scripts" / "explorer" / name, f"{REMOTE_ROOT}/")
    # CRLF would make the shebang unparseable on Linux; strip it rather than trusting
    # git's autocrlf to have done the right thing on a Windows checkout.
    fix = _ssh(f"cd {REMOTE_ROOT} && sed -i 's/\\r$//' *.sbatch && chmod +x *.sbatch "
               f"&& head -1 boltz_depth.sbatch | cat -A | head -1", timeout=180)
    out["dos2unix_and_chmod"] = fix.stdout.strip()
    return out


# ---------------------------------------------------------------------------
# verify -- big-file checks run INSIDE a job, never on the login node
# ---------------------------------------------------------------------------

def verify() -> dict:
    """Submit the md5/import check as a batch job and wait for it.

    The login node's memory cgroup SIGKILLs long or memory-hungry reads, and a killed
    process can still exit 0 -- so a login-node md5 of a 2.3 GB checkpoint is exactly the
    check that reports success while having read nothing. This runs it with an explicit
    --mem on a compute node.
    """
    sub = _ssh(f"cd {REMOTE_ROOT} && sbatch --parsable verify_inputs.sbatch", timeout=180)
    jid = sub.stdout.strip().splitlines()[-1] if sub.stdout.strip() else ""
    if not jid.isdigit():
        return {"ok": False, "submit_stdout": sub.stdout[-400:],
                "submit_stderr": sub.stderr[-400:]}
    st = _wait(jid, poll=20, limit=40)
    log = _ssh(f"cat {REMOTE_ROOT}/logs/verify_{jid}.log", timeout=180)
    return {"ok": st.get("state", "").startswith("COMPLETED"),
            "job": jid, "state": st.get("state"), "log": log.stdout[-3000:]}


def _wait(jid: str, poll: int = 60, limit: int = 120) -> dict:
    """Poll for ADVANCING progress, never for liveness."""
    last = None
    for i in range(limit):
        q = _ssh(f"squeue -j {jid} -h -o '%T %M' 2>/dev/null; "
                 f"sacct -j {jid} -n -X -o State%20,Elapsed 2>/dev/null | head -1",
                 timeout=120)
        txt = q.stdout.strip()
        state = txt.split()[0] if txt else "UNKNOWN"
        if state != last:
            print(f"  [{i}] {jid}: {txt}", flush=True)
            last = state
        if state.startswith(("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
                             "OUT_OF_ME", "NODE_FAIL")):
            return {"state": state, "raw": txt}
        time.sleep(poll)
    return {"state": "TIMEOUT_WAITING", "raw": last}


# ---------------------------------------------------------------------------
# submit
# ---------------------------------------------------------------------------

def submit(seed: int, smoke: bool = False, partition: str = "gyorilab,gpu",
           samples: int | None = None, tag: str | None = None) -> dict:
    sub = _tag_for(smoke, tag)
    n = samples if samples is not None else SAMPLES
    cmd = (f"cd {REMOTE_ROOT} && sbatch --parsable --partition={partition} "
           f"boltz_depth.sbatch {sub} {seed} {n}")
    r = _ssh(cmd, timeout=300)
    jid = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else ""
    rec = {"job": jid, "sub": sub, "seed": seed, "samples": n,
           "partition": partition, "submitted_at": time.time(),
           "stdout": r.stdout[-300:], "stderr": r.stderr[-300:]}
    if jid.isdigit():
        d = json.loads(JOBS_JSON.read_text()) if JOBS_JSON.exists() else {"jobs": []}
        d["jobs"].append(rec)
        JOBS_JSON.write_text(json.dumps(d, indent=1))
    return rec


def poll(jid: str | None = None) -> dict:
    """ADVANCING progress: how many mmCIFs exist, not whether the job is alive."""
    d = json.loads(JOBS_JSON.read_text()) if JOBS_JSON.exists() else {"jobs": []}
    jids = [jid] if jid else [j["job"] for j in d["jobs"] if j["job"].isdigit()]
    q = _ssh("squeue -u shenoy.am -o '%.10i %.9P %.14j %.2t %.10M %R'; echo '---'; "
             f"find {REMOTE_ROOT}/out -name '*_model_*.cif' 2>/dev/null | wc -l; "
             f"echo '---'; ls -la {REMOTE_ROOT}/logs/ 2>/dev/null | tail -8",
             timeout=180)
    states = {}
    for j in jids:
        s = _ssh(f"sacct -j {j} -n -X -o State%20,Elapsed,NodeList%15", timeout=120)
        states[j] = s.stdout.strip()
    return {"queue_and_progress": q.stdout, "states": states}


def tail(jid: str, n: int = 40) -> str:
    r = _ssh(f"tail -{n} {REMOTE_ROOT}/logs/*_{jid}.log 2>/dev/null", timeout=180)
    return r.stdout + r.stderr


# ---------------------------------------------------------------------------
# collect
# ---------------------------------------------------------------------------

def collect(smoke: bool = False, tag: str | None = None) -> dict:
    """Pull the predicted mmCIFs back and flatten them to <ligand>__s<seed>_<model>.cif."""
    sub = _tag_for(smoke, tag)
    dest = POSE_DIR / sub
    dest.mkdir(parents=True, exist_ok=True)
    # boltz writes out/<sub>_s<seed>/boltz_results_<sub>/predictions/<lig>/<lig>_model_k.cif
    # Pull with tar-over-ssh (rsync is unusable from this box) and filter to the mmCIFs so
    # the per-pose npz confidence blobs stay on /scratch -- D: has ~13 GB and is not a
    # place to land a pool.
    # scope the find to THIS sub's output dirs. Pulling all of out/ would sweep the smoke
    # run's poses into the stratum set (and vice versa) -- the flattening step keys on the
    # ligand directory name and would not notice.
    r = _bash(
        f"ssh -o BatchMode=yes explorer "
        f"\"cd {REMOTE_ROOT}/out && find . -path './{sub}_s*' -name '*_model_*.cif' "
        f"-print0 | tar czf - --null -T -\" | tar xzf - -C '{_win_to_posix(dest)}'; "
        f"echo SSH_RC=${{PIPESTATUS[0]}} TAR_RC=${{PIPESTATUS[1]}}", timeout=7200)
    if "SSH_RC=0 TAR_RC=0" not in r.stdout:
        return {"ok": False, "rc": r.stdout[-200:], "err": r.stderr[-500:]}

    flat = POSE_DIR / (sub + "_flat")
    flat.mkdir(parents=True, exist_ok=True)
    n = 0
    seen: dict[str, int] = {}
    for f in sorted(dest.rglob("*_model_*.cif")):
        # .../out/<sub>_s<seed>/boltz_results_*/predictions/<lig>/<lig>_model_k.cif
        seed = "0"
        for part in f.parts:
            if part.startswith(sub + "_s"):
                seed = part.split("_s")[-1]
        lig = f.parent.name
        k = f.stem.rsplit("_model_", 1)[-1]
        tgt = flat / f"{lig}__s{seed}_m{k}.cif"
        tgt.write_bytes(f.read_bytes())
        seen[lig] = seen.get(lig, 0) + 1
        n += 1
    pl = _load_plan(sub) or {}
    planned = pl.get("ligands") or []
    res = {"ok": True, "tag": sub, "collected": n, "flat_dir": str(flat),
           "per_ligand": seen, "ligands": len(seen)}
    if planned:
        # Check 4 of Step 2: every planned ligand present, at the same depth. A count of
        # files says nothing about which ligand is missing 20 of them.
        res["planned"] = len(planned)
        res["missing_ligands"] = sorted(set(planned) - set(seen))
        depths = sorted(set(seen.values()))
        res["poses_per_ligand"] = depths
        res["uniform_depth"] = len(depths) == 1 and not res["missing_ligands"]
    return res


# ---------------------------------------------------------------------------
# score -- identical machinery to FINDING 034's scorer, different input dir
# ---------------------------------------------------------------------------

def score(smoke: bool = False, tag: str | None = None) -> dict:
    """Score against crystals. **Cannot run on a blind set** -- see G2.

    Every ligand needs a deposited structure; without one `load_crystal` returns None and
    the ligand is counted into `skipped["no_crystal"]`, so a blind set returns a CSV with
    zero rows and a cheerful `ok`. That is the shape of failure this repo keeps paying
    for, so it is now refused up front rather than reported as an empty success.
    """
    from cypstruct import pose as P
    from cypstruct import xengine as X
    from cypstruct.qmscore import geometry as G
    from two_ligand_cofold import load_crystal, renumber_to_reference

    sub = _tag_for(smoke, tag)
    pl = _load_plan(sub) or {}
    tbl = _ligand_table(pl.get("csv"))
    if "pdb" not in tbl.columns:
        raise SystemExit(
            f"tag {sub!r} was planned from {pl.get('csv')!r}, which has no `pdb` column: "
            "there are no crystals to score against and every ligand would be counted as "
            "`no_crystal`, returning zero rows. Run `collect` and build the selection "
            "feature from the flat pose directory instead (playbook G2).")
    refset = X.load_reference(DATA_PROCESSED / "reference_set_cyp3a4.npz")
    vl = tbl.set_index("id")
    ids = pl.get("ligands") or (
        [smoke_ligand()] if smoke else sorted(type_i_ligands().id.tolist()))
    flat = POSE_DIR / (sub + "_flat")

    rows = []
    skipped = {"no_crystal": 0, "no_ligand": 0, "load_failed": 0, "mapping_failed": 0,
               "no_heme_frame": 0, "low_renumber_identity": 0, "no_files": 0}
    for lid in ids:
        pdb = vl.loc[lid, "pdb"]
        smiles = vl.loc[lid, "smiles"]
        ref = load_crystal(pdb, lid)
        if ref is None or len(ref.lig_xyz) == 0:
            skipped["no_crystal"] += 1
            continue
        ref_pocket = P.pocket_residues_from_structure(ref, radius=8.0)
        files = sorted(flat.glob(f"{lid}__*.cif"))
        if not files:
            skipped["no_files"] += 1
            continue
        for f in files:
            try:
                cx = P.load_structure(f)
            except Exception:
                skipped["load_failed"] += 1
                continue
            if len(cx.lig_xyz) == 0:
                skipped["no_ligand"] += 1
                continue
            cx, off, ident = renumber_to_reference(cx, ref)
            if ident < 0.95:
                skipped["low_renumber_identity"] += 1
            perm = P.best_ligand_mapping(smiles, cx, ref)
            mapped = perm is not None
            if perm is None:
                skipped["mapping_failed"] += 1
                perm = np.arange(min(len(cx.lig_xyz), len(ref.lig_xyz)))
            ld = P.lddt_pli(cx, ref, lig_perm=perm)
            bs = P.bisy_rmsd(cx, ref, align_resnums=ref_pocket, lig_perm=perm)
            geo = G.compute(cx.lig_xyz, cx.lig_elem, cx.prot_xyz,
                            cx.heme_xyz, cx.heme_atom, cx.axial_sg)
            hf = X.in_heme_frame(cx)
            if hf is None:
                skipped["no_heme_frame"] += 1
            xe = (X.xeng_score(hf, refset.get(lid, [])) if hf is not None
                  else float("nan"))
            rows.append({
                "ligand": lid, "pdb": pdb, "arm": ARM, "sample": f.stem.split("__")[1],
                "lddt_pli": ld, "bisy_rmsd": bs, "xeng": xe,
                "fe_donor_dist": geo.fe_donor_dist,
                "s_fe_donor_angle": geo.s_fe_donor_angle,
                "is_coordinated": geo.is_coordinated,
                "frac_proximal": geo.frac_atoms_proximal,
                "renumber_offset": off, "renumber_identity": round(ident, 4),
                "mapped": mapped,
                "md5": hashlib.md5(f.read_bytes()).hexdigest(),
                "n_ref_poses": len(refset.get(lid, [])),
            })
    res = pd.DataFrame(rows)
    if sub == "stratum":
        out_csv = POSES_CSV
    elif sub == "smoke":
        out_csv = DATA_PROCESSED / "matched_depth_smoke_poses.csv"
    else:
        out_csv = DATA_PROCESSED / f"matched_depth_{sub}_poses.csv"
    res.to_csv(out_csv, index=False)
    return {"arm": ARM, "rows": int(len(res)),
            "ligands": int(res.ligand.nunique()) if len(res) else 0,
            "distinct_md5": int(res.md5.nunique()) if len(res) else 0,
            "skipped": skipped,
            "renumber_offsets": (res.renumber_offset.value_counts().to_dict()
                                 if len(res) else {}),
            "min_renumber_identity": (float(res.renumber_identity.min())
                                      if len(res) else None),
            "n_exact_zero_lddt": int((res.lddt_pli == 0).sum()) if len(res) else 0,
            "mean_lddt": float(res.lddt_pli.mean()) if len(res) else None,
            "out": str(out_csv)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["plan", "conditioning", "stage", "push", "verify",
                                    "submit", "poll", "tail", "collect", "score"])
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--csv", default=None,
                    help="ligand csv (id,smiles) to plan over. Without it, plan keeps "
                         "its historical behaviour: the 14-ligand predicted-Type-I "
                         "stratum of validation_ligands.csv.")
    ap.add_argument("--tag", default=None,
                    help="subdirectory name used by plan/stage/submit/collect/score on "
                         "BOTH sides. Defaults to 'stratum' ('smoke' with --smoke).")
    ap.add_argument("--seed", type=int, default=101)
    ap.add_argument("--samples", type=int, default=None)
    ap.add_argument("--job", default=None)
    ap.add_argument("--partition", default="gyorilab,gpu")
    a = ap.parse_args()

    if a.cmd == "plan":
        out = plan(a.smoke, a.csv, a.tag)
    elif a.cmd == "conditioning":
        out = conditioning()
    elif a.cmd == "stage":
        out = stage(a.smoke, a.tag)
    elif a.cmd == "push":
        out = push_sbatch()
    elif a.cmd == "verify":
        out = verify()
    elif a.cmd == "submit":
        out = submit(a.seed, a.smoke, a.partition, a.samples, a.tag)
    elif a.cmd == "poll":
        out = poll(a.job)
    elif a.cmd == "tail":
        print(tail(a.job or ""))
        return 0
    elif a.cmd == "collect":
        out = collect(a.smoke, a.tag)
    else:
        out = score(a.smoke, a.tag)
    print(json.dumps(out, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
