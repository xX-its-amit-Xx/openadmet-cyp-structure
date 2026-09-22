"""Does an APO CYP3A4 crystal predict the holo F/G lesion (210-216)? The caveat FINDING
028's first addendum left open.

That addendum refuted HOLO templates over 210-216: CYP3A4's own deposited crystals
disagree with each other by 4.25 A median / 7.20 A p90 (holo-holo, different ligand,
1,424 pairs) where the model predicts the same span to 1.03 A median / 4.98 A p90. But
the 107-entry harvest contained ZERO apo entries, so FINDING 024's separate claim -- that
the only always-legal blind template is an APO structure -- was recorded as untested
rather than refuted. This script closes it.

Two questions:

  Q1  How well does an apo entry predict a holo crystal's 210-216 backbone? Superpose
      each apo entry onto each of the 87 validation crystals on CYP3A4_RIGID_CORE (a
      ligand-free frame declared in targets.py long before any of this) and measure the
      CA deviation over the span alone.

  Q2  Would an apo 210-216 backbone + CB actually ADMIT the query's crystal ligand?
      Graft the superposed apo span into the query frame and apply the same a-priori
      2.2 A clearance test the first addendum used.

DECISION RULE, fixed before any number was read (echoed into the output JSON):

  An apo template is LICENSED only if its CA deviation over 210-216 against the 87
  validation crystals is BELOW THE MODEL -- median < 1.03 A AND p90 < 4.98 A. The bar is
  the model, not the holo templates. Being better than the 4.25 A holo-holo spread is
  NOT sufficient and must not be reported as a win.

Superposition, span, clearance cutoff and the rigid-core frame are IMPORTED from
side_chain_diagnosis.py rather than rewritten, so every number here is commensurable with
the first addendum by construction rather than by assertion.

    python scripts/structure/apo_template_ceiling.py            # everything
    python scripts/structure/apo_template_ceiling.py --classify # apo discovery only
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

# Imported, not re-implemented: the span, the rigid-core frame, the per-entry chain
# reader (one chain per entry, backbone+CB of the span), the crystal loader and
# FINDING 021's renumbering. Same code path as the first addendum.
from side_chain_diagnosis import (  # noqa: E402
    CORE, LESION, PREREG, _entry_chain, load_reference, renumber_to_reference,
    residue_atoms,
)

RCSB = REPO / "data" / "reference" / "rcsb"
POOL = Path("D:/cyp_scratch/val87b_unsteered")

# The span's identity in UniProt P08684 numbering. Checked residue by residue in EVERY
# entry before any geometry is trusted -- FINDING 021's trap, and a numbering slip here
# would fake the entire result.
EXPECTED_SPAN = {210: "LEU", 211: "LEU", 212: "ARG", 213: "PHE",
                 214: "ASP", 215: "PHE", 216: "LEU"}

# The control span used by the first addendum: the I-helix, inside the rigid core. In a
# correct superposition it must agree to 0.4-1.3 A even on pairs whose lesion is 3-12 A
# apart. If it does not, the superposition is wrong, not the biology.
CONTROL_SPAN = list(range(300, 321))
CONTROL_OK = (0.4, 1.3)

PREREG_APO = {
    "span": [LESION[0], LESION[-1]],
    "frame": "CYP3A4_RIGID_CORE (ligand-free, declared in targets.py)",
    "min_core_CA_shared": 50,
    "min_span_residues": 4,
    "clash_cut_A": PREREG["clash_cut_A"],          # 2.2 A, inherited, a priori
    "bar_is_the_model": {"median_A": 1.03, "p90_A": 4.98},
    "holo_holo_reference_not_the_bar": {"median_A": 4.25, "p90_A": 7.20},
    "decision_rule": (
        "LICENSED only if the apo->holo spread over 210-216 is below the MODEL on BOTH "
        "median (<1.03 A) and p90 (<4.98 A). Beating the 4.25 A holo-holo spread is not "
        "sufficient and is not a win."),
    "apo_definition": (
        "UniProt P08684 entry whose only non-water heteroatom components are HEM and "
        "members of cypstruct.targets.IGNORE_HET (buffer/cryoprotectant). Reported "
        "alongside the closest non-HEM het atom to the heme Fe, so a cryoprotectant "
        "sitting in the pocket cannot hide."),
    "control_I_helix": {"span": [CONTROL_SPAN[0], CONTROL_SPAN[-1]],
                        "required_range_A": list(CONTROL_OK)},
}


# --------------------------------------------------------------------------
# apo discovery
# --------------------------------------------------------------------------

def classify_entries() -> dict:
    """Every deposited P08684 entry, split apo / holo by its heteroatom content."""
    import gemmi
    from cypstruct.targets import IGNORE_HET, rcsb_holo_structures, fetch_cif

    ids = rcsb_holo_structures(limit=500)
    rows, downloaded = {}, []
    for pid in ids:
        cif = RCSB / f"{pid}.cif"
        if not cif.exists():
            cif = fetch_cif(pid)          # small; lands beside the rest, deletes nothing
            downloaded.append(pid)
        st = gemmi.read_structure(str(cif))
        st.setup_entities()
        het, fe, hetxyz = {}, [], []
        for ch in st[0]:
            for r in ch:
                nm = r.name.strip().upper()
                info = gemmi.find_tabulated_residue(nm)
                if info and (info.is_amino_acid() or info.is_water()):
                    continue
                het[nm] = het.get(nm, 0) + 1
                for at in r:
                    if nm == "HEM":
                        if at.name.strip() == "FE":
                            fe.append([at.pos.x, at.pos.y, at.pos.z])
                    else:
                        hetxyz.append((nm, [at.pos.x, at.pos.y, at.pos.z]))
        nondrug = sorted(k for k in het if k not in IGNORE_HET and k != "HEM")
        near = None
        if fe and hetxyz:
            F = np.array(fe)
            d = [(float(np.linalg.norm(np.array(x) - F, axis=1).min()), n)
                 for n, x in hetxyz]
            dd, nn = min(d)
            near = {"code": nn, "dist_to_nearest_Fe_A": round(dd, 2)}
        try:
            res = float(st.resolution)
        except Exception:
            res = float("nan")
        rows[pid] = {"het": sorted(het), "non_ignored_het": nondrug,
                     "apo": not nondrug, "resolution": res,
                     "spacegroup": str(st.spacegroup_hm),
                     "closest_non_HEM_het_to_Fe": near}
    apo = sorted(p for p, v in rows.items() if v["apo"])
    return {"n_entries_queried": len(ids), "n_downloaded_now": len(downloaded),
            "downloaded": downloaded, "apo_ids": apo, "n_apo": len(apo),
            "n_holo": len(ids) - len(apo), "entries": rows}


# --------------------------------------------------------------------------
# geometry helpers (superposition reused from the first addendum's frame)
# --------------------------------------------------------------------------

def span_dev(donor_ca: dict, query_ca: dict, span: list[int]):
    """Superpose donor onto query on CYP3A4_RIGID_CORE CAs; return per-residue deviation
    over `span`, the core fit, and the I-helix control, all in the SAME superposition."""
    from cypstruct import pose as P
    sh = sorted(set(donor_ca) & set(query_ca) & CORE)
    if len(sh) < PREREG_APO["min_core_CA_shared"]:
        return None
    R, t, fit = P.kabsch(np.array([donor_ca[k] for k in sh]),
                         np.array([query_ca[k] for k in sh]))
    use = [r for r in span if r in donor_ca and r in query_ca]
    d = {r: float(np.linalg.norm(donor_ca[r] @ R.T + t - query_ca[r])) for r in use}
    cs = [r for r in CONTROL_SPAN if r in donor_ca and r in query_ca]
    ctrl = (float(np.sqrt(np.mean([np.linalg.norm(donor_ca[r] @ R.T + t - query_ca[r])
                                   ** 2 for r in cs]))) if len(cs) >= 10 else float("nan"))
    return {"R": R, "t": t, "core_fit": float(fit), "n_core": len(sh),
            "dev": d, "control_rms": ctrl, "n_control": len(cs)}


def rms(vals) -> float:
    v = [x for x in vals if np.isfinite(x)]
    return float(np.sqrt(np.mean(np.square(v)))) if v else float("nan")


def stat_block(x) -> dict | None:
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    if not len(x):
        return None
    return {"n": int(len(x)), "median_A": float(np.median(x)),
            "p25_A": float(np.percentile(x, 25)), "p75_A": float(np.percentile(x, 75)),
            "p90_A": float(np.percentile(x, 90)), "max_A": float(x.max()),
            "frac_over_1A": float((x > 1.0).mean()), "frac_over_2A": float((x > 2.0).mean())}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--classify", action="store_true",
                    help="stop after apo discovery and the modelling audit")
    a = ap.parse_args()

    import pandas as pd
    from cypstruct import pose as P

    out_json = REPO / "data" / "processed" / "apo_template_ceiling.json"
    cls_path = REPO / "data" / "processed" / "apo_entry_classification.json"

    # ---- 1. which entries are apo ---------------------------------------
    cls = classify_entries()
    cls_path.write_text(json.dumps(cls, indent=1))
    apo_ids = cls["apo_ids"]
    print(f"[classify] {cls['n_entries_queried']} P08684 entries, "
          f"{cls['n_apo']} apo: {apo_ids} "
          f"({cls['n_downloaded_now']} downloaded now)", flush=True)

    # ---- 2. FIRST: do the apo entries model the span at all? -------------
    apo = {}
    audit, name_fail = [], []
    for pid in apo_ids:
        e = _entry_chain(RCSB / f"{pid}.cif", None)
        if e is None:
            audit.append({"pdb": pid, "status": "no usable chain"})
            continue
        # the numbering trap: check the span really is LLRFDFL before trusting anything
        import gemmi
        st = gemmi.read_structure(str(RCSB / f"{pid}.cif"))
        st.setup_entities()
        names = {}
        for ch in st[0]:
            if ch.name != e["chain"]:
                continue
            for r in ch:
                if r.seqid.num in EXPECTED_SPAN:
                    names[r.seqid.num] = r.name.strip().upper()
        bad = {k: v for k, v in names.items() if v != EXPECTED_SPAN[k]}
        if bad:
            name_fail.append({"pdb": pid, "mismatch": bad})
        mod = [r for r in LESION if r in e["ca"]]
        audit.append({
            "pdb": pid, "chain": e["chain"], "resolution": round(e["resolution"], 2),
            "spacegroup": e["spacegroup"],
            "span_residues_modelled": len(mod), "modelled": mod,
            "missing": [r for r in LESION if r not in e["ca"]],
            "CA_B_factors": {str(r): round(e["b"][r], 1) for r in mod},
            "median_CA_B": round(float(np.median([e["b"][r] for r in mod])), 1) if mod else None,
            "min_CA_occupancy": round(float(min([e["occ"][r] for r in mod])), 2) if mod else None,
            "residue_names": {str(k): v for k, v in sorted(names.items())},
            "residue_names_match_L210_L211_R212_F213_D214_F215_L216": not bad,
            "closest_non_HEM_het_to_Fe": cls["entries"][pid]["closest_non_HEM_het_to_Fe"],
        })
        apo[pid] = e
    print(json.dumps({"apo_span_audit": audit,
                      "C_numbering_mismatches": name_fail}, indent=1), flush=True)

    usable = [p for p in apo if sum(1 for r in LESION if r in apo[p]["ca"])
              >= PREREG_APO["min_span_residues"]]
    result = {
        "what": "can an APO CYP3A4 crystal supply the holo 210-216 backbone the model "
                "gets wrong? The caveat FINDING 028's first addendum left open.",
        "prereg": PREREG_APO,
        "apo_discovery": {k: v for k, v in cls.items() if k != "entries"},
        "apo_span_audit": audit,
        "C_numbering_mismatches_apo": name_fail,
        "n_apo_modelling_the_span": len(usable),
        "apo_modelling_the_span": usable,
    }
    if not usable:
        result["verdict"] = ("NO APO ENTRY MODELS 210-216. A template cannot carry what "
                             "the crystal does not contain; apo templates are refuted "
                             "without any spread being computed.")
        out_json.write_text(json.dumps(result, indent=1))
        print(json.dumps(result["verdict"], indent=1))
        return
    if a.classify:
        out_json.write_text(json.dumps(result, indent=1))
        return

    # ---- 3. the 87 validation crystals ----------------------------------
    lig = pd.read_csv(REPO / "data" / "processed" / "validation_ligands.csv")
    have = {p.name.split("__")[0] for p in POOL.iterdir() if p.is_dir()} \
        if POOL.exists() else set(lig.id)
    queries, qbad = {}, []
    for r in lig.itertuples():
        if r.id not in have:
            continue
        ref, _cif, _ch = load_reference(r.pdb, r.id)
        if ref is None or not len(ref.lig_xyz):
            continue
        ca = {k: v for (_c, k), v in ref.ca().items()}
        nm = {n: v for (_c, n), v in ref.prot_res.items() if n in EXPECTED_SPAN}
        bad = {k: v for k, v in nm.items() if v != EXPECTED_SPAN[k]}
        if bad:
            qbad.append({"ligand": r.id, "pdb": r.pdb, "mismatch": bad})
        own = {}
        for rn in LESION:
            for nmm, xyz in residue_atoms(ref, rn).items():
                if nmm in ("N", "CA", "C", "O", "CB"):
                    own[(rn, nmm)] = xyz
        queries[r.id] = {"pdb": str(r.pdb), "ca": ca, "lig": np.asarray(ref.lig_xyz, float),
                         "span_modelled": [x for x in LESION if x in ca], "own_bbcb": own}
    print(f"[queries] {len(queries)} validation crystals loaded, "
          f"{len(qbad)} numbering mismatches", flush=True)

    # holo donors, for an apples-to-apples contrast in the SAME pairing (donor -> the 87
    # query crystals), so apo is not compared against a differently-constructed number
    holo_ids = [p for p in cls["entries"] if not cls["entries"][p]["apo"]]
    holo = {}
    for pid in holo_ids:
        e = _entry_chain(RCSB / f"{pid}.cif", None)
        if e is not None and sum(1 for r in LESION if r in e["ca"]) == len(LESION):
            holo[pid] = e
    print(f"[donors] {len(usable)} apo, {len(holo)} holo with the span fully modelled",
          flush=True)

    # ---- 4. Q1: the spread, and the I-helix control ----------------------
    rows = []
    for kind, donors in (("apo", {p: apo[p] for p in usable}), ("holo", holo)):
        for pid, e in donors.items():
            for lgid, q in queries.items():
                if pid.upper() == q["pdb"].upper():
                    continue                      # never grade a donor on itself
                sp = span_dev(e["ca"], q["ca"], LESION)
                if sp is None:
                    rows.append(dict(kind=kind, donor=pid, ligand=lgid, pdb=q["pdb"],
                                     n_span=0, rms=np.nan, dropped="core CA < 50"))
                    continue
                if len(sp["dev"]) < PREREG_APO["min_span_residues"]:
                    rows.append(dict(kind=kind, donor=pid, ligand=lgid, pdb=q["pdb"],
                                     n_span=len(sp["dev"]), rms=np.nan,
                                     dropped="span < 4 residues in common"))
                    continue
                rows.append(dict(
                    kind=kind, donor=pid, ligand=lgid, pdb=q["pdb"],
                    n_span=len(sp["dev"]), n_core=sp["n_core"],
                    core_fit=sp["core_fit"], rms=rms(sp["dev"].values()),
                    max_dev=float(max(sp["dev"].values())),
                    control_ihelix_rms=sp["control_rms"], dropped=""))
    pairs = pd.DataFrame(rows)
    pairs.to_csv(REPO / "data" / "processed" / "apo_template_pairs.csv", index=False)

    ok = pairs[pairs.rms.notna()]
    apo_ok, holo_ok = ok[ok.kind == "apo"], ok[ok.kind == "holo"]
    q1 = {
        "apo_to_holo_query": stat_block(apo_ok.rms),
        "holo_to_holo_query_SAME_pairing": stat_block(holo_ok.rms),
        "per_apo_donor": {p: stat_block(apo_ok[apo_ok.donor == p].rms)
                          for p in sorted(apo_ok.donor.unique())},
        "best_single_apo_donor_by_median": (
            min(((p, float(apo_ok[apo_ok.donor == p].rms.median()))
                 for p in apo_ok.donor.unique()), key=lambda z: z[1])
            if len(apo_ok) else None),
        "oracle_best_apo_donor_per_query": stat_block(
            apo_ok.groupby("ligand").rms.min()) if len(apo_ok) else None,
        "note_oracle": "the per-query minimum is NOT available blind; reported only to "
                       "show the ceiling of apo templating even with an oracle donor",
        "median_core_fit_A": float(ok.core_fit.median()),
    }

    # Every apo donor returned the SAME frac_over_2A. That is either a bug or a fact,
    # so it is checked rather than reported: do the six apo entries occupy ONE loop
    # conformation? Two ways of asking -- their agreement with each other directly, and
    # their agreement with each other per query.
    aa = []
    for i in range(len(usable)):
        for j in range(i + 1, len(usable)):
            sp = span_dev(apo[usable[i]]["ca"], apo[usable[j]]["ca"], LESION)
            if sp and len(sp["dev"]) == len(LESION):
                aa.append(dict(a=usable[i], b=usable[j],
                               rms=rms(sp["dev"].values()),
                               control_ihelix_rms=sp["control_rms"]))
    aad = pd.DataFrame(aa)
    piv = apo_ok.pivot_table(index="ligand", columns="donor", values="rms")
    rng = piv.max(axis=1) - piv.min(axis=1)
    q1["apo_are_one_conformation"] = {
        "apo_apo_pairs": int(len(aad)),
        "apo_apo_span_spread": stat_block(aad.rms) if len(aad) else None,
        "apo_apo_I_helix_control": stat_block(aad.control_ihelix_rms) if len(aad) else None,
        "per_query_range_across_the_6_apo_donors": stat_block(rng),
        "queries_where_all_6_agree_within_0.5A": int((rng < 0.5).sum()),
        "n_queries": int(len(piv)),
        "queries_where_EVERY_apo_donor_is_under_1A": int((piv.max(axis=1) < 1.0).sum()),
        "queries_where_EVERY_apo_donor_is_over_2A": int((piv.min(axis=1) > 2.0).sum()),
        "queries_intermediate": int(len(piv) - (piv.max(axis=1) < 1.0).sum()
                                    - (piv.min(axis=1) > 2.0).sum()),
        "note": "identical frac_over_2A across donors is explained, not a filter that "
                "never fires: the apo entries are one conformer, so they succeed and "
                "fail on the same queries",
    }
    q1["per_query_apo_median"] = {
        str(k): round(float(v), 2) for k, v in piv.median(axis=1).sort_values().items()}

    # every filter, with counts
    filters = {
        "pairs_attempted": int(len(pairs)),
        "dropped_core_CA_under_50": int((pairs.dropped == "core CA < 50").sum()),
        "dropped_span_under_4_residues": int(
            (pairs.dropped == "span < 4 residues in common").sum()),
        "pairs_used": int(len(ok)),
        "apo_pairs_used": int(len(apo_ok)), "holo_pairs_used": int(len(holo_ok)),
        "queries_with_span_fully_modelled": int(sum(
            len(q["span_modelled"]) == len(LESION) for q in queries.values())),
        "queries_with_span_under_4_residues": int(sum(
            len(q["span_modelled"]) < PREREG_APO["min_span_residues"]
            for q in queries.values())),
        "C_numbering_mismatches_validation_crystals": qbad,
    }

    # the I-helix control, on exactly the pairs the first addendum would have flagged
    big = ok[ok.rms > 3.0]
    ctrl = {
        "all_pairs_I_helix_rms": stat_block(ok.control_ihelix_rms),
        "pairs_with_lesion_over_3A": int(len(big)),
        "I_helix_on_those": stat_block(big.control_ihelix_rms),
        "frac_of_those_inside_0.4_1.3A": float(
            ((big.control_ihelix_rms >= CONTROL_OK[0]) &
             (big.control_ihelix_rms <= CONTROL_OK[1])).mean()) if len(big) else None,
        "worst_examples": [
            {"donor": r.donor, "query": f"{r.ligand}/{r.pdb}", "kind": r.kind,
             "lesion_rms_A": round(float(r.rms), 2),
             "I_helix_rms_A": round(float(r.control_ihelix_rms), 2),
             "core_fit_A": round(float(r.core_fit), 2)}
            for _i, r in big.sort_values("rms", ascending=False).head(6).iterrows()],
        "apo_pairs_with_lesion_over_3A": int((apo_ok.rms > 3.0).sum()),
    }

    # ---- 4b. the head-to-head the decision rule actually asks for --------
    # The 1.03 A bar is the model's own error over the SAME span in the SAME frame. It
    # is recomputed here rather than quoted, on the same 58 scorable queries, so the
    # comparison is per query and not between two differently-pooled distributions.
    mcsv = REPO / "data" / "processed" / "apo_model_span_error.csv"
    if mcsv.exists():
        md = pd.read_csv(mcsv)
    else:
        mrows = []
        for lgid, q in queries.items():
            job = POOL / f"{lgid}__unsteered__s1"
            span = q["span_modelled"]
            if len(span) < PREREG_APO["min_span_residues"] or not job.exists():
                continue
            ref, _cif, _ch = load_reference(q["pdb"], lgid)
            for cif in sorted(job.glob("input_model_*.cif")):
                m = P.load_structure(cif)
                m, _o, _i = renumber_to_reference(m, ref)
                mn = {k: v for (_c, k), v in m.ca().items()}
                sp = span_dev(mn, q["ca"], LESION)
                if sp is None or len(sp["dev"]) < PREREG_APO["min_span_residues"]:
                    continue
                mrows.append(dict(ligand=lgid, pdb=q["pdb"], sample=cif.stem,
                                  n_span=len(sp["dev"]), rms=rms(sp["dev"].values()),
                                  control_ihelix_rms=sp["control_rms"]))
        md = pd.DataFrame(mrows)
        md.to_csv(mcsv, index=False)

    model_per_lig = md.groupby("ligand").rms.median()
    apo_per_lig = apo_ok.groupby("ligand").rms.median()
    hh = pd.concat([model_per_lig.rename("model"), apo_per_lig.rename("apo")],
                   axis=1).dropna()
    q1["model_recomputed_here"] = stat_block(md.rms)
    q1["head_to_head_per_query"] = {
        "n_queries": int(len(hh)),
        "queries_where_apo_beats_the_model": int((hh.apo < hh.model).sum()),
        "queries_where_the_model_beats_apo": int((hh.model < hh.apo).sum()),
        "median_apo_minus_model_A": float((hh.apo - hh.model).median()),
        "apo_median_of_per_query_medians_A": float(hh.apo.median()),
        "model_median_of_per_query_medians_A": float(hh.model.median()),
        "binomial_p_apo_better": float(
            __import__("scipy.stats", fromlist=["x"]).binomtest(
                int((hh.apo < hh.model).sum()), int(len(hh)), 0.5).pvalue),
        "wilcoxon_p": float(
            __import__("scipy.stats", fromlist=["x"]).wilcoxon(hh.apo, hh.model).pvalue),
        "note": "a per-query paired comparison, not two pooled distributions",
    }

    # ---- 5. Q2: would the grafted apo span admit the crystal ligand? -----
    cut = PREREG_APO["clash_cut_A"]
    hard7 = ["5AW", "ERY", "A1A4T", "1RD", "MWY", "X7P", "QEP"]
    clear_rows = []
    for lgid, q in queries.items():
        L = q["lig"]
        own = np.array(list(q["own_bbcb"].values()))
        own_min = (float(np.linalg.norm(own[:, None, :] - L[None, :, :], axis=2).min())
                   if len(own) else float("nan"))
        for pid in usable:
            e = apo[pid]
            if pid.upper() == q["pdb"].upper():
                continue
            sp = span_dev(e["ca"], q["ca"], LESION)
            if sp is None:
                continue
            pts = np.array([v @ sp["R"].T + sp["t"] for v in e["bbcb"].values()])
            if not len(pts):
                continue
            dmin = float(np.linalg.norm(pts[:, None, :] - L[None, :, :], axis=2).min())
            clear_rows.append(dict(ligand=lgid, pdb=q["pdb"], donor=pid,
                                   n_bbcb=len(pts), min_contact=dmin,
                                   clears=bool(dmin >= cut),
                                   own_crystal_min_contact=own_min,
                                   hard7=lgid in hard7))
    cr = pd.DataFrame(clear_rows)
    cr.to_csv(REPO / "data" / "processed" / "apo_clearance.csv", index=False)

    per_lig = cr.groupby("ligand").agg(n_donors=("clears", "size"),
                                       n_clearing=("clears", "sum"),
                                       best=("min_contact", "max")).reset_index()
    own_ctrl = cr.groupby("ligand").own_crystal_min_contact.first()
    q2 = {
        "cutoff_A": cut,
        "n_query_ligands": int(len(per_lig)),
        "n_apo_donors": len(usable),
        "grafts_tested": int(len(cr)),
        "grafts_clearing": int(cr.clears.sum()),
        "frac_grafts_clearing": float(cr.clears.mean()),
        "ligands_with_at_least_one_clearing_apo_donor": int((per_lig.n_clearing > 0).sum()),
        "ligands_with_NO_clearing_apo_donor": int((per_lig.n_clearing == 0).sum()),
        "ligands_all_apo_donors_clear": int(
            (per_lig.n_clearing == per_lig.n_donors).sum()),
        "median_min_contact_A": float(cr.min_contact.median()),
        "CONTROL_own_crystal_span_vs_own_ligand": {
            "n": int(own_ctrl.notna().sum()),
            "median_A": float(own_ctrl.median()),
            "min_A": float(own_ctrl.min()),
            "frac_below_cut": float((own_ctrl < cut).mean()),
            "note": "the query's OWN deposited 210-216 backbone+CB against its own "
                    "ligand; this is the number a correct span must beat. It is only "
                    "meaningful where the query models the span, so it is stratified.",
        },
        "CONTROL_own_crystal_span_FULLY_MODELLED_ONLY": (lambda s: {
            "n": int(len(s)), "median_A": float(s.median()) if len(s) else None,
            "min_A": float(s.min()) if len(s) else None,
            "frac_below_cut": float((s < cut).mean()) if len(s) else None})(
            own_ctrl[[lg for lg in own_ctrl.index
                      if len(queries[lg]["span_modelled"]) == len(LESION)]].dropna()),
        "by_loop_state_of_the_query": (lambda: {
            st: {"n_ligands": int(len(v)),
                 "frac_grafts_clearing": float(
                     cr[cr.ligand.isin(v)].clears.mean()) if len(v) else None,
                 "ligands_with_no_clearing_donor": int(
                     (per_lig[per_lig.ligand.isin(v)].n_clearing == 0).sum())}
            for st, v in {
                "query in the apo state (apo span dev < 1 A)": [
                    lg for lg in piv.index if piv.loc[lg].max() < 1.0],
                "query in a different state (apo span dev > 2 A)": [
                    lg for lg in piv.index if piv.loc[lg].min() > 2.0],
                "query span not scorable (< 4 residues modelled)": [
                    lg for lg in queries
                    if len(queries[lg]["span_modelled"]) < PREREG_APO["min_span_residues"]],
            }.items()})(),
        "seven_unrescuable": [
            {"ligand": lg, "pdb": str(cr[cr.ligand == lg].pdb.iloc[0]),
             "own_span_residues_modelled": len(queries[lg]["span_modelled"]),
             "donors_tested": int((cr.ligand == lg).sum()),
             "donors_clearing": int(cr[cr.ligand == lg].clears.sum()),
             "best_donor": str(cr[cr.ligand == lg].sort_values(
                 "min_contact", ascending=False).donor.iloc[0]),
             "best_min_contact_A": round(float(cr[cr.ligand == lg].min_contact.max()), 2),
             "own_crystal_min_contact_A": (
                 round(float(own_ctrl[lg]), 2) if np.isfinite(own_ctrl.get(lg, np.nan))
                 else None)}
            for lg in hard7 if lg in set(cr.ligand)],
    }

    # ---- 6. verdict against the rule fixed in advance ---------------------
    bar = PREREG_APO["bar_is_the_model"]
    med = q1["apo_to_holo_query"]["median_A"]
    p90 = q1["apo_to_holo_query"]["p90_A"]
    licensed = bool(med < bar["median_A"] and p90 < bar["p90_A"])
    result.update({
        "filters": filters,
        "Q1_spread_over_210_216": q1,
        "C_I_helix_control": ctrl,
        "Q2_clearance": q2,
        "verdict": {
            "apo_median_A": med, "apo_p90_A": p90,
            "model_median_A": bar["median_A"], "model_p90_A": bar["p90_A"],
            "licensed": licensed,
            "statement": ("LICENSED" if licensed else
                          "REFUTED -- an apo template does not predict 210-216 better "
                          "than the model, which is the bar that was fixed in advance"),
        },
    })
    out_json.write_text(json.dumps(result, indent=1))
    print(json.dumps({k: result[k] for k in
                      ("filters", "Q1_spread_over_210_216", "C_I_helix_control",
                       "Q2_clearance", "verdict")}, indent=1))


if __name__ == "__main__":
    main()
