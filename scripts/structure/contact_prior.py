"""Score a pose by how closely its pocket contacts match what real CYP3A4 ligands touch.

Two orientation features have now been tested and rejected (contact-consensus, which was
redundant with the contact count; and azimuth about the iron axis, which was orthogonal
but carried no signal). Both asked whether a pose agrees with *the other samples*. This
asks something different: does it agree with *the crystallography*.

The prior comes from the deposited structures themselves. For each of ~116 CYP3A4 entries
we already have locally, take the set of residues within 4.5 A of its ligand. Across the
set that gives a per-residue frequency: how often does a real bound ligand touch Ser119,
Phe304, Ala305, and so on. A predicted pose that engages the residues real ligands engage
is more plausible than one that makes the same *number* of contacts elsewhere - which is
all `n_pocket_residues_touched` can currently tell.

**Leakage control.** The prior for scoring ligand L is built from every deposited structure
EXCEPT the entry L itself, so a ligand never votes on its own pose. Without that this
feature would trivially recover the answer and report a spectacular, meaningless gain.
That is the exact failure mode CLAUDE.md warns about under "no leaky features".

    python scripts/structure/contact_prior.py build
    python scripts/structure/contact_prior.py test --tag val87b
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from cypstruct import pose as P  # noqa: E402
from cypstruct.paths import DATA_PROCESSED, REFERENCE  # noqa: E402

PRIOR = DATA_PROCESSED / "cyp3a4_contact_prior.json"
CUTOFF = 4.5


def build() -> None:
    """Per-residue ligand-contact frequency across every deposited CYP3A4 holo structure."""
    import gemmi

    ref_geo = pd.read_parquet(DATA_PROCESSED / "cyp3a4_reference_geometry.parquet")
    holo = ref_geo[(ref_geo.ligand != "")
                   & (ref_geo.binding_class != "peripheral")]
    per_entry: dict[str, dict] = {}
    cifs = sorted((REFERENCE / "rcsb").glob("*.cif"))
    print(f"{len(cifs)} cached mmCIFs, {holo.pdb_id.nunique()} with a non-peripheral ligand",
          flush=True)

    for n, cif in enumerate(cifs, 1):
        pdb = cif.stem.upper()
        rows = holo[holo.pdb_id == pdb]
        if rows.empty:
            continue
        for lig_code in sorted(rows.ligand.unique()):
            try:
                st = gemmi.read_structure(str(cif))
                st.setup_entities()
                chain = next((c.name for c in st[0]
                              if any(r.name.strip().upper() == lig_code for r in c)), None)
                if chain is None:
                    continue
                cx = P.load_structure(cif, ligand_code=lig_code, assembly_chain=chain)
            except Exception:
                continue
            if len(cx.lig_xyz) == 0 or len(cx.prot_xyz) == 0:
                continue
            d = np.linalg.norm(cx.prot_xyz[:, None, :] - cx.lig_xyz[None, :, :],
                               axis=2).min(1)
            res = sorted({cx.prot_key[i][1] for i in np.where(d <= CUTOFF)[0]})
            per_entry[f"{pdb}:{lig_code}"] = {"pdb": pdb, "ligand": lig_code,
                                              "residues": res}
        if n % 25 == 0:
            print(f"  {n}/{len(cifs)} parsed, {len(per_entry)} contact sets", flush=True)

    counts = Counter()
    for v in per_entry.values():
        counts.update(v["residues"])
    n_entries = max(1, len(per_entry))
    freq = {int(r): c / n_entries for r, c in counts.items()}
    PRIOR.write_text(json.dumps({"n_entries": n_entries, "cutoff": CUTOFF,
                                 "freq": freq, "per_entry": per_entry}, indent=1))
    top = sorted(freq.items(), key=lambda kv: -kv[1])[:15]
    print(f"\n{len(per_entry)} contact sets over {len({v['pdb'] for v in per_entry.values()})} entries")
    print("most frequently contacted residues:")
    for r, f in top:
        print(f"   {r:4d}  {f:.2f}")
    print(f"-> {PRIOR}")


def contacts(tag: str, arm: str, workers: int, profile: str | None) -> None:
    """Per-pose residue contact sets, pulled once and cached.

    The orientation pass computed this fingerprint and then threw away everything except
    its cardinality, which is why `n_pocket_residues_touched` could only say HOW MANY
    residues a pose touches and never WHICH. The prior needs the identities.
    """
    import os
    import shutil
    import tempfile
    import threading
    from concurrent.futures import ThreadPoolExecutor, as_completed

    if profile:
        os.environ["MODAL_PROFILE"] = profile
    import modal

    from cypstruct.paths import free_gb, safe_workers

    scored = pd.read_csv(DATA_PROCESSED / f"poses_scored_{tag}.csv")
    scored = scored[scored.arm == arm]
    jobs = sorted({f"{r.ligand}__{arm}__s1" for r in scored.itertuples()})
    vol = modal.Volume.from_name("cyp-pool")
    root = Path("D:/cyp_scratch")
    root.mkdir(parents=True, exist_ok=True)
    if free_gb(root) < 5.0:
        root = Path("C:/cyp_scratch")
        root.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=f"cprior_{tag}_", dir=str(root)))
    nw = safe_workers(workers, ram_per_worker_gb=1.0)
    print(f"{len(jobs)} jobs, {nw} worker(s)", flush=True)
    out: dict[str, list[int]] = {}
    lock = threading.Lock()

    def do(job: str) -> None:
        lig = job.split("__")[0]
        w = tmp / job
        w.mkdir(parents=True, exist_ok=True)
        try:
            for e in vol.iterdir(f"/{tag}/{job}"):
                fn_ = e.path.split("/")[-1]
                if fn_.endswith(".cif"):
                    (w / fn_).write_bytes(b"".join(vol.read_file(e.path)))
            for cif in sorted(w.glob("*.cif")):
                try:
                    m = P.load_structure(cif)
                except Exception:
                    continue
                if len(m.lig_xyz) == 0 or len(m.prot_xyz) == 0:
                    continue
                d = np.linalg.norm(m.prot_xyz[:, None, :] - m.lig_xyz[None, :, :],
                                   axis=2).min(1)
                res = sorted({m.prot_key[i][1] for i in np.where(d <= CUTOFF)[0]})
                with lock:
                    out[f"{lig}|{cif.stem}"] = res
        except Exception:
            pass
        finally:
            shutil.rmtree(w, ignore_errors=True)

    try:
        with ThreadPoolExecutor(max_workers=nw) as ex:
            futs = [ex.submit(do, j) for j in jobs]
            for k, f in enumerate(as_completed(futs), 1):
                f.result()
                if k % 20 == 0:
                    print(f"  {k}/{len(jobs)}", flush=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    dest = DATA_PROCESSED / f"pose_contacts_{tag}_{arm}.json"
    dest.write_text(json.dumps(out))
    print(f"{len(out)} poses -> {dest}")

def test(tag: str, arm: str = "unsteered") -> None:
    from scipy import stats

    pri = json.loads(PRIOR.read_text())
    per_entry = pri["per_entry"]

    scored = pd.read_csv(DATA_PROCESSED / f"poses_scored_{tag}.csv")
    scored = scored[scored.arm == arm]
    orient = DATA_PROCESSED / f"orientation_features_{tag}_{arm}.csv"
    cons = DATA_PROCESSED / f"consensus_features_{tag}_{arm}.csv"
    feats = pd.read_csv(orient).merge(pd.read_csv(cons), on=["ligand", "sample"], how="left")

    # per-pose contact sets, recomputed the same way as the prior
    cfile = DATA_PROCESSED / f"pose_contacts_{tag}_{arm}.json"
    if not cfile.exists():
        raise SystemExit(f"missing {cfile}; run `contacts` first")
    pose_contacts = json.loads(cfile.read_text())

    rows = []
    for key, res in pose_contacts.items():
        lig, sample = key.split("|", 1)
        # LEAVE-ONE-LIGAND-OUT prior: drop every entry whose ligand is this one
        freq = Counter()
        n = 0
        for ek, ev in per_entry.items():
            if ev["ligand"] == lig:
                continue
            n += 1
            freq.update(ev["residues"])
        if n == 0:
            continue
        s = set(res)
        # mean prior frequency of the residues this pose touches
        score = float(np.mean([freq.get(r, 0) / n for r in s])) if s else 0.0
        # and the recall of the commonly-contacted residues
        common = {r for r, c in freq.items() if c / n >= 0.5}
        recall = len(s & common) / max(1, len(common))
        rows.append({"ligand": lig, "sample": sample,
                     "prior_mean_freq": score, "prior_recall": recall,
                     "n_contacts": len(s)})
    pf = pd.DataFrame(rows)
    df = scored.merge(feats, on=["ligand", "sample"]).merge(pf, on=["ligand", "sample"])
    df = df.dropna(subset=["lddt_pli", "mean_rmsd_to_others", "n_pocket_residues_touched"])
    print(f"{len(df)} poses / {df.ligand.nunique()} ligands\n")

    zw = lambda c: ((df[c] - df.groupby("ligand")[c].transform("mean"))
                    / (df.groupby("ligand")[c].transform("std") + 1e-9))

    def ev(score, label):
        d = df.assign(_s=score)
        pick, rnd = [], []
        for _l, g in d.groupby("ligand"):
            if len(g) < 2:
                continue
            pick.append(g.loc[g._s.idxmax(), "lddt_pli"])
            rnd.append(g.lddt_pli.mean())
        pick, rnd = np.array(pick), np.array(rnd)
        w = stats.wilcoxon(pick, rnd)
        print(f"  {label:46s} delta={(pick-rnd).mean():+.4f} "
              f"beat={100*(pick>rnd).mean():4.1f}% p={w.pvalue:.4f}")
        return (pick - rnd).mean()

    ev(zw("prior_mean_freq"), "crystallographic contact prior (mean freq)")
    ev(zw("prior_recall"), "recall of commonly-contacted residues")
    base = 0.5 * zw("n_pocket_residues_touched") - zw("mean_rmsd_to_others")
    ev(base, "FINDING 003 baseline")
    for w in (0.25, 0.5, 1.0):
        ev(base + w * zw("prior_recall"), f"  baseline + {w} * prior_recall")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "contacts", "test"])
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--profile", default=None)
    ap.add_argument("--tag", default="val87b")
    ap.add_argument("--arm", default="unsteered")
    a = ap.parse_args()
    if a.cmd == "build":
        build()
    elif a.cmd == "contacts":
        contacts(a.tag, a.arm, a.workers, a.profile)
    else:
        test(a.tag, a.arm)
