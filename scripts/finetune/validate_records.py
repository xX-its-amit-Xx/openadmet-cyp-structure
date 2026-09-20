"""Featurize every record once, and drop the ones that raise.

Runs ON EXPLORER, before training. The 350-step arm4_mix run died three minutes in on
`KeyError: 'C8'` inside `process_atom_features` - a ligand whose mmCIF atom names do not
all appear in its reference conformer. One record out of 388, four hours of GPU wasted.

Two ways to survive that, and only one of them is honest:

- Catch inside `__getitem__` and fall back to another index. Boltz does this. At N=303
  one bad record silently becomes a large share of an epoch, and nothing says so.
- Do the pass up front, on CPU, and write manifests containing only records that actually
  featurize. Training then has no surprises and the drop count is a number we can report.

This is the second. `manifest_train.json` / `manifest_test.json` are REWRITTEN in place;
the originals are kept as `*_full.json` on first run, so re-running is idempotent and the
unfiltered set is never lost.

    ./env/bin/python validate_records.py --arm arm4_mix --workers 8
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path("/scratch/shenoy.am/cyp-finetune")
sys.path.insert(0, str(ROOT))


def main() -> int:
    import torch
    from boltz.data.types import Manifest

    from boltz2_data import Boltz2FinetuneDataset

    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="arm4_mix")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-tokens", type=int, default=384)
    ap.add_argument("--max-atoms", type=int, default=3072)
    ap.add_argument("--max-drop-frac", type=float, default=0.15,
                    help="abort rather than train on a set this badly degraded")
    a = ap.parse_args()

    arm_dir = ROOT / "records" / a.arm
    report: dict[str, dict] = {}

    class Probe(torch.utils.data.Dataset):
        """Returns the error text instead of raising, so one bad item cannot kill the pass."""

        def __init__(self, ds):
            self.ds = ds

        def __len__(self):
            return len(self.ds)

        def __getitem__(self, i):
            rid = self.ds.records[i].id
            try:
                self.ds[i]
            except Exception as exc:  # noqa: BLE001
                return (rid, f"{type(exc).__name__}: {str(exc)[:90]}")
            return (rid, "")

    for split in ("train", "test"):
        full = arm_dir / f"manifest_{split}_full.json"
        live = arm_dir / f"manifest_{split}.json"
        if not full.exists():
            full.write_text(live.read_text())
        manifest = Manifest.load(full)

        ds = Boltz2FinetuneDataset(
            manifest, arm_dir, ROOT / "msa_npz", ROOT / "boltz_cache" / "mols",
            max_tokens=a.max_tokens, max_atoms=a.max_atoms,
            training=False,   # deterministic: a record that fails here fails every epoch
        )
        loader = torch.utils.data.DataLoader(
            Probe(ds), batch_size=1, num_workers=a.workers,
            collate_fn=lambda b: b[0], shuffle=False)

        ok, bad = [], {}
        for n, (rid, err) in enumerate(loader, 1):
            if err:
                bad[rid] = err
            else:
                ok.append(rid)
            if n % 50 == 0:
                print(f"  {a.arm}/{split} {n}/{len(ds)}  bad={len(bad)}", flush=True)

        keep = [r for r in manifest.records if r.id in set(ok)]
        frac = 1 - len(keep) / max(len(manifest.records), 1)
        Manifest(keep).dump(live)

        reasons: dict[str, int] = {}
        for e in bad.values():
            reasons[e.split(":")[0] + ": " + e.split(": ", 1)[-1][:40]] = (
                reasons.get(e.split(":")[0] + ": " + e.split(": ", 1)[-1][:40], 0) + 1)
        report[split] = {
            "total": len(manifest.records), "kept": len(keep), "dropped": len(bad),
            "dropped_frac": round(frac, 4),
            "dropped_ids": sorted(bad),
            "reasons": dict(sorted(reasons.items(), key=lambda kv: -kv[1])[:6]),
        }
        if frac > a.max_drop_frac:
            print(json.dumps(report, indent=2))
            msg = (f"{a.arm}/{split}: dropped {frac:.1%} of records, over the "
                   f"{a.max_drop_frac:.0%} limit. Fix the parser, do not train on this.")
            raise SystemExit(msg)

    (arm_dir / "validation_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "dropped_ids"}
                      for k, v in report.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
