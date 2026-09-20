"""Rewrite the held-out YAMLs with an explicit Cys-SG -> heme-FE bond.

Runs ON EXPLORER. This is the experiment the fine-tuning control pointed at: the base
checkpoint puts 44/85 held-out poses under 2 A and drops 39/85 into catastrophe, and it
does that with NO heme bond. This repo has already measured that an explicit bond is what
moves Boltz-2 to crystallographic Fe-donor geometry - median 2.23 A with 84% inside the
window, where an unbonded pilot put zero poses in range. If the bond converts
catastrophes into successes it is worth many times the -0.0125 fine-tuning returned, and
it costs no training at all.

**Finding the ligating cysteine without looking at the answer.** P450s carry an
absolutely conserved heme-binding motif near the C-terminus, FxxGxxxCxG, whose cysteine
is the proximal thiolate. That is a SEQUENCE fact, so it is computable for a blind target
from its sequence alone - unlike reading the residue out of the crystal, which would be
leaky and would not transfer to the challenge set. The search is deliberately restricted
to the C-terminal 40% and takes the last match, because scattered `G...C.G` can occur by
chance earlier in a 490-residue protein.

Everything else about these YAMLs is byte-identical to build_holdout_yaml.py, so the
comparison against the existing base control isolates the bond.

    ./env/bin/python build_holdout_yaml_bonded.py --arm arm4_mix
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

ROOT = Path("/scratch/shenoy.am/cyp-finetune")

# FxxGxxxCxG and its common relaxations. Tried in order; the first that matches wins, so
# the strictest evidence is preferred and the permissive pattern is only a fallback.
MOTIFS = (
    r"F..G...C.G",
    r"[FWY]..G...C.G",
    r"G...C.G",
)


def find_axial_cys(seq: str) -> tuple[int | None, str]:
    """1-based index of the heme-ligating cysteine, and which motif found it."""
    start = int(len(seq) * 0.55)          # C-terminal ~45%; the motif sits ~50 res from the end
    tail = seq[start:]
    for pat in MOTIFS:
        hits = list(re.finditer(pat, tail))
        if hits:
            m = hits[-1]                   # last match: the motif is the C-terminal-most one
            cys_in_match = m.group().index("C")
            return start + m.start() + cys_in_match + 1, pat
    return None, ""


def main() -> int:
    import pandas as pd
    from boltz.data.types import Manifest

    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="arm4_mix")
    a = ap.parse_args()

    arm_dir = ROOT / "records" / a.arm
    test_ids = {r.id for r in Manifest.load(arm_dir / "manifest_test.json").records}
    df = pd.read_csv(ROOT / "finetune_arms" / f"{a.arm}.csv")
    df = df[(df["split"] == "test") & (df["pdb"].isin(test_ids))]

    out = ROOT / "holdout_yaml" / f"{a.arm}_bonded"
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob("*.yaml"):
        old.unlink()

    written = no_cys = no_msa = 0
    by_motif: dict[str, int] = {}
    positions = []
    for r in df.itertuples():
        a3m = ROOT / "msa" / f"{r.target_key}.a3m"
        if not a3m.exists():
            no_msa += 1
            continue
        seq = str(r.sequence)
        cys, motif = find_axial_cys(seq)
        if cys is None:
            # No bond rather than a guessed one. A wrong bond is far worse than none: it
            # would pin the heme to an arbitrary residue and wreck a pose that was fine.
            no_cys += 1
            continue
        by_motif[motif] = by_motif.get(motif, 0) + 1
        positions.append(round(cys / len(seq), 3))

        smiles = str(r.smiles).replace("'", "")
        yaml = (
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: {seq}\n"
            f"      msa: {a3m}\n"
            "  - ligand:\n"
            "      id: B\n"
            "      ccd: HEM\n"
            "  - ligand:\n"
            "      id: C\n"
            f"      smiles: '{smiles}'\n"
            "constraints:\n"
            "  - bond:\n"
            f"      atom1: [A, {cys}, SG]\n"
            "      atom2: [B, 1, FE]\n"
        )
        (out / f"{r.pdb}_{r.id}.yaml").write_text(yaml)
        written += 1

    print(json.dumps({
        "arm": a.arm, "pairs_written": written,
        "skipped_no_cys": no_cys, "skipped_no_msa": no_msa,
        "motif_used": by_motif,
        "cys_relative_position": {
            "min": min(positions) if positions else None,
            "median": sorted(positions)[len(positions) // 2] if positions else None,
            "max": max(positions) if positions else None,
        },
        "out": str(out),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
