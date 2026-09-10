"""CYP target definitions — sequences, the heme, and the pocket that actually matters.

Primary target for the OpenADMET structure track is **CYP3A4** (UniProt P08684).
The other three challenge isoforms are kept here because the activity track needs
them and because cross-isoform folds are a useful control.

Three CYP-specific facts drive every downstream design decision in this repo:

1. **The iron already has five ligands.** Four pyrrole nitrogens of protoporphyrin IX
   plus the *proximal* thiolate of **Cys442**. Only the sixth, distal coordination
   site is available to a ligand. Any pose that puts ligand density on the proximal
   face is physically wrong no matter what the co-folder's confidence says, and a
   scorer that does not know this cannot tell the two faces apart.

2. **Type II inhibitors coordinate the iron.** An sp2 nitrogen lone pair (imidazole,
   triazole, pyridine) donates into Fe(III) d_z2 at **~2.0-2.3 A**. This is a dative
   bond with a hard directional requirement, not a soft contact. A prior pilot in the
   sibling activity repo cofolded 24 CYP ligands and got **zero** poses inside that
   window (closest 2.80 A) — the co-folders miss the coordination geometry entirely.
   That is the single largest, most CYP-specific scoring signal available to us.

3. **The pocket is plastic, and it is plastic in one place.** OpenADMET's own cryoEM
   announcement (2026-09-08) states the remodelling that reshapes the site is
   concentrated in the **F/G loop**, which is exactly the region left unmodelled in
   many X-ray structures — and names this as the reason co-folding does poorly here.
   So F/G-loop residues are excluded from rigid reference-frame alignment (they move
   legitimately) and are handled as an explicit conformational-ensemble axis instead.

Residue numbering below follows full-length UniProt P08684. Crystal/cryoEM constructs
are N-terminally truncated (the transmembrane anchor is removed) but deposit with
UniProt numbering, so these indices transfer; `pocket_residues()` re-validates against
whatever structure is loaded rather than trusting the index blindly.
"""
from __future__ import annotations

import json
from pathlib import Path

import requests

from .paths import DATA_EXTERNAL, REFERENCE

SEQ_CACHE = DATA_EXTERNAL / "cyp_sequences.fasta"

CYP_TARGETS: dict[str, dict] = {
    "cyp3a4": {"uniprot": "P08684", "axial_cys": 442, "length": 503},
    "cyp2c9": {"uniprot": "P11712", "axial_cys": 435, "length": 490},
    "cyp2d6": {"uniprot": "P10635", "axial_cys": 443, "length": 497},
    "cyp1a2": {"uniprot": "P05177", "axial_cys": 458, "length": 516},
}

PRIMARY = "cyp3a4"

# --- heme ------------------------------------------------------------------
# What we ASK a co-folder to build in (CCD code for heme b / iron-protoporphyrin IX):
HEME_CCD = "HEM"
# What we ACCEPT when parsing deposited structures. CYPs do not all deposit as HEM —
# e.g. CYP2C9 1OG5/1OG2 use HEC — and matching only "HEM" silently reports those as
# apo, which quietly shrinks the validation set while looking like a clean result.
HEME_ALIASES = frozenset({"HEM", "HEC", "HEA", "HEB", "HDD", "DHE", "SRM", "VER", "1CP", "MH0"})

# Solvent/cryo junk that shows up as "ligand" in deposited CIFs and must never be
# treated as the ligand of interest.
IGNORE_HET = frozenset({
    "HOH", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "1PE", "MPD", "TRS",
    "CL", "NA", "K", "MG", "CA", "ZN", "ACT", "DMS", "IMD", "FMT", "NO3", "IOD",
    "BME", "CIT", "MES", "EPE", "TLA", "NH4", "CO3", "AZI", "ACY", "P6G", "P33",
})

# --- CYP3A4 active site (UniProt P08684 numbering) --------------------------
# Sources: the standard CYP3A4 crystallographic literature (1TQN/1W0E/2V0M/3NXU
# families). These are the residues that line the distal cavity and repeatedly
# appear in ligand-contact analyses. Treated as a PRIOR, not a fact: every script
# that uses them re-derives actual contacts from the structure at hand.
CYP3A4_POCKET = {
    # I-helix, lying directly across the distal face of the heme
    "i_helix": [301, 304, 305, 309, 312],           # Ile301, Phe304, Ala305, Thr309
    # the phenylalanine cluster forming the roof / peripheral site
    "phe_cluster": [57, 108, 213, 215, 219, 220, 241, 304],
    # polar anchors most often seen H-bonding substrates
    "polar": [106, 119, 212, 214, 373, 374, 375],   # Arg106, Ser119, Arg212, Asp214
    # heme propionate salt bridges (structural, not ligand-contacting)
    "heme_anchor": [372, 375, 440, 441, 442],
    # substrate access channel mouth
    "access": [50, 74, 120, 210, 211, 482, 483, 484],
}

# The plastic region. Excluded from rigid alignment; modelled as an ensemble axis.
# F-helix -> F' -> G' -> G-helix. This is the span OpenADMET flags as the remodelling
# hotspot and as the reason co-folding underperforms on CYP3A4.
CYP3A4_FG_LOOP = (202, 260)

# Residues used as the ALIGNMENT FRAME for pose comparison: the rigid structural
# core. Deliberately the heme-proximal half, which does not move between holo forms.
CYP3A4_RIGID_CORE = [
    *range(280, 340),   # I-helix and flanks
    *range(350, 400),   # K-helix / meander
    *range(425, 470),   # heme-binding loop incl. Cys442, the most rigid element
]


def pocket_atoms(isoform: str = PRIMARY) -> list[int]:
    """Flat, deduped, sorted list of pocket residue numbers for an isoform."""
    if isoform != "cyp3a4":
        raise NotImplementedError(
            f"pocket prior is curated for cyp3a4 only; got {isoform}. "
            "Derive contacts from a holo structure instead of guessing."
        )
    out: set[int] = set()
    for v in CYP3A4_POCKET.values():
        out.update(v)
    return sorted(out)


def in_fg_loop(resnum: int, isoform: str = PRIMARY) -> bool:
    """True if a residue sits in the plastic F/G region (excluded from rigid fits)."""
    if isoform != "cyp3a4":
        return False
    lo, hi = CYP3A4_FG_LOOP
    return lo <= resnum <= hi


# --- sequences -------------------------------------------------------------
def fetch_sequences(force: bool = False) -> dict[str, str]:
    """Canonical UniProt sequences for the four isoforms, cached to a FASTA."""
    if SEQ_CACHE.exists() and not force:
        return _read_fasta(SEQ_CACHE)
    seqs: dict[str, str] = {}
    for iso, meta in CYP_TARGETS.items():
        acc = meta["uniprot"]
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=30)
        r.raise_for_status()
        seqs[iso] = "".join(r.text.splitlines()[1:])
    SEQ_CACHE.parent.mkdir(parents=True, exist_ok=True)
    with open(SEQ_CACHE, "w") as fh:
        for iso, seq in seqs.items():
            fh.write(f">{iso}|{CYP_TARGETS[iso]['uniprot']}\n{seq}\n")
    return seqs


def _read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split("|")[0]
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line.strip()
    return seqs


def verify_axial_cys(isoform: str = PRIMARY) -> bool:
    """Sanity-check that the recorded axial cysteine index really is a Cys.

    Cheap, and it catches an off-by-one in the numbering convention before that
    error propagates into every geometric term in the scorer.
    """
    seq = fetch_sequences()[isoform]
    idx = CYP_TARGETS[isoform]["axial_cys"]
    return len(seq) >= idx and seq[idx - 1] == "C"


# --- experimental structures ----------------------------------------------
_RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"


def rcsb_holo_structures(isoform: str = PRIMARY, limit: int = 200) -> list[str]:
    """Query RCSB for every deposited structure of this isoform.

    Deliberately a live query rather than a hardcoded PDB list: the CYP3A4 holo
    set is large and still growing, a stale list silently shrinks the validation
    set, and a hand-typed list is exactly where a wrong code hides for weeks.
    """
    acc = CYP_TARGETS[isoform]["uniprot"]
    query = {
        "query": {
            "type": "group", "logical_operator": "and", "nodes": [
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers"
                                 ".reference_sequence_identifiers.database_accession",
                    "operator": "exact_match", "value": acc}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers"
                                 ".reference_sequence_identifiers.database_name",
                    "operator": "exact_match", "value": "UniProt"}},
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit},
                            "results_verbosity": "compact"},
    }
    r = requests.post(_RCSB_SEARCH, json=query, timeout=60)
    r.raise_for_status()
    return sorted(r.json().get("result_set", []))


def fetch_cif(pdb_id: str, dest: Path | None = None) -> Path:
    """Download one mmCIF into cold storage (never onto the near-full local disks)."""
    dest = Path(dest or REFERENCE / "rcsb")
    dest.mkdir(parents=True, exist_ok=True)
    out = dest / f"{pdb_id.upper()}.cif"
    if out.exists() and out.stat().st_size > 0:
        return out
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", timeout=120)
    r.raise_for_status()
    out.write_bytes(r.content)
    return out


if __name__ == "__main__":
    seqs = fetch_sequences()
    print("sequences:")
    for iso, s in seqs.items():
        cys_ok = verify_axial_cys(iso)
        print(f"  {iso:8s} {CYP_TARGETS[iso]['uniprot']}  {len(s):4d} aa   "
              f"axial Cys{CYP_TARGETS[iso]['axial_cys']} verified={cys_ok}")
    print(f"\ncyp3a4 pocket prior: {len(pocket_atoms())} residues")
    print(f"F/G plastic span: {CYP3A4_FG_LOOP}  (excluded from rigid alignment)")
    ids = rcsb_holo_structures()
    print(f"\nRCSB entries for CYP3A4 (P08684): {len(ids)}")
    print("  " + " ".join(ids[:40]) + (" ..." if len(ids) > 40 else ""))
    (REFERENCE / "rcsb_cyp3a4_entries.json").write_text(json.dumps(ids, indent=1))
