"""IDEA B3 — differential fuzzing: perturb the INPUT, not the engine.

Cross-engine agreement is the only selector that has beaten random here (+0.0381), and
FINDING 016 showed the same mechanism works when the perturbation is a sampler setting
rather than a different architecture. This asks whether it also works when the
perturbation is *chemical*: fold tautomers of the same ligand and treat their poses as
independent opinions about where the molecule sits.

Why it is not merely another sweep. Sampler settings perturb the *search*; tautomers
perturb the *question*. If a pose is robust to which tautomer you ask about, that is
evidence about the pocket rather than about the sampler - a different and possibly
less-correlated axis. FINDING 013 says added diversity usually enlarges the oracle without
becoming findable, so the prior is not good; the point is that this axis is free and has
never been tried.

**The property that makes it valid at all:** tautomer and protonation changes preserve the
HEAVY-atom composition (verified 38/38 on our ligand set), so a heme-frame Chamfer
distance between a tautomer's pose and the original pool is comparing like with like. A
chirality flip would also preserve heavy atoms but changes the molecule, so it belongs
here only as a negative control - if flipped-chirality poses score as good a reference as
tautomers, the "reference" is measuring nothing chemical.

    python scripts/cofold/fuzz_reference.py submit --pairs 60
    python scripts/cofold/fuzz_reference.py collect
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "cofold"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402

UNI = DATA_PROCESSED / "p450_universe"
OUT = UNI / "fuzz"
STATE = OUT / "state.json"
ENGINE = "protenix_v2"
N_TAUT = 3          # perturbed variants per ligand, beyond the original


def _state() -> dict:
    return json.loads(STATE.read_text()) if STATE.exists() else {"folds": {}}


def _save(d: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    try:
        base = json.loads(STATE.read_text()) if STATE.exists() else {}
    except (OSError, json.JSONDecodeError):
        base = {}
    merged = dict(base.get("folds", {}))
    merged.update(d.get("folds", {}))
    d["folds"] = merged
    if STATE.exists():
        try:
            shutil.copyfile(STATE, STATE.with_suffix(".json.prev"))
        except OSError:
            pass
    tmp = STATE.with_suffix(f".json.tmp{os.getpid()}")
    tmp.write_text(json.dumps(d, indent=1))
    os.replace(tmp, STATE)


def variants(smiles: str, n: int) -> list[tuple[str, str]]:
    """Distinct tautomers, plus a chirality-flipped NEGATIVE CONTROL where possible."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem.MolStandardize import rdMolStandardize
    RDLogger.DisableLog("rdApp.*")

    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return []
    heavy = m.GetNumHeavyAtoms()
    canon = Chem.MolToSmiles(m)
    out, seen = [], {canon}
    for i, t in enumerate(rdMolStandardize.TautomerEnumerator().Enumerate(m)):
        s = Chem.MolToSmiles(t)
        # a variant that changed the heavy-atom cloud is not a perturbation of the same
        # question, and its poses could not be Chamfer-compared to the original pool
        if s in seen or t.GetNumHeavyAtoms() != heavy:
            continue
        seen.add(s)
        out.append((f"taut{i}", s))
        if len(out) >= n:
            break
    centres = Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
    if centres:
        mm = Chem.MolFromSmiles(canon)
        idx, _ = centres[0]
        a = mm.GetAtomWithIdx(idx)
        a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW
                       if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW
                       else Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        s = Chem.MolToSmiles(mm)
        if s != canon:
            out.append(("chiralflip", s))      # negative control, not a reference member
    return out


def cmd_submit(n_pairs: int, limit: int | None) -> dict:
    from openprotein_cofold import build_complex, connect
    s = connect()
    st = _state()
    camp = json.loads((UNI / "campaign.json").read_text())
    ready = {k for k, v in camp["msa"].items() if v.get("status") == "SUCCESS"}
    cs = pd.read_csv(UNI / "p450_cofold_set.csv")
    cs["pair"] = cs.pdb + "_" + cs.id
    unsup = UNI / f"{ENGINE}_unsupported.json"
    skip = set(json.loads(unsup.read_text()).get("pairs", [])) if unsup.exists() else set()
    cs = cs[cs.target_key.isin(ready) & ~cs.pair.isin(skip)].drop_duplicates("pair")
    cs = cs.head(n_pairs)

    claimed = {v["key"] for v in st["folds"].values()}
    n = 0
    for r in cs.itertuples():
        vs = variants(r.smiles, N_TAUT)
        if not vs:
            continue
        msa = s.load_job(camp["msa"][r.target_key]["job_id"])
        for kind, smi in vs:
            key = f"{r.pair}:{kind}"
            if key in claimed:
                continue
            if limit and n >= limit:
                break
            try:
                fut = getattr(s.fold, ENGINE).fold(
                    sequences=[build_complex(r.sequence, smi, msa)],
                    diffusion_samples=1, num_recycles=3)
            except Exception as exc:
                st["folds"][f"unsent:{key}"] = {
                    "key": key, "pair": r.pair, "kind": kind, "smiles": smi,
                    "done": "failed", "error": f"{type(exc).__name__}: {exc}"[:150]}
                _save(st); n += 1
                continue
            st["folds"][str(fut.job_id)] = {
                "key": key, "pair": r.pair, "kind": kind, "smiles": smi,
                "submitted": time.time()}
            _save(st); n += 1
            if n % 25 == 0:
                print(f"  submitted {n}", flush=True)
    return {"submitted": n, "pairs": int(cs.shape[0]), "records": len(st["folds"])}


def cmd_collect() -> dict:
    from openprotein_cofold import _split_models, connect
    s = connect()
    st = _state()
    ok = pend = fail = 0
    for job_id, rec in st["folds"].items():
        if rec.get("done") is True:
            ok += 1; continue
        if rec.get("done") == "failed":
            fail += 1; continue
        try:
            f = s.load_job(job_id)
            status = str(f.job.status).upper()
            if "FAIL" in status or "CANCEL" in status:
                rec["done"] = "failed"; rec["error"] = status; fail += 1; continue
            if "SUCCESS" not in status:
                pend += 1; continue
            res = f.get()
        except Exception:
            pend += 1; continue
        d = OUT / "poses" / rec["pair"]
        d.mkdir(parents=True, exist_ok=True)
        try:
            rec["files"] = _split_models(res[0].to_string(), rec["pair"], d, rec["kind"])
            rec["done"] = True; ok += 1
        except Exception as exc:
            rec["done"] = "failed"; rec["error"] = f"parse: {exc}"[:150]; fail += 1
    _save(st)
    return {"collected": ok, "pending": pend, "failed": fail}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["submit", "collect"])
    ap.add_argument("--pairs", type=int, default=60)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    fn = cmd_submit(a.pairs, a.limit or None) if a.cmd == "submit" else cmd_collect()
    print(json.dumps(fn, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
