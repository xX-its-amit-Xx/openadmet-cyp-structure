"""FINDING 031 analysis — did a second ligand copy change the first copy's pose?

Reads `data/processed/two_ligand_poses_<engine>.csv` (written by
`scripts/cofold/two_ligand_cofold.py score`) and answers, in the order the
pre-registration fixed them:

  1. the paired primary endpoint, arm B - arm A and arm C - arm A, on the copy
     nearest the iron, with the FINDING 024 decomposition alongside the raw score
  2. where the second copy went, against the measured three-site geometry
  3. whether the F/G span 210-216 moved between arms
  4. pool oracle, and only then pool selection

Every threshold used here is quoted from `docs/PREREG_two_ligand_cofold.md`. None was
chosen after a scored result was seen.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts" / "cofold"))

from cypstruct.paths import DATA_PROCESSED  # noqa: E402
import two_ligand_cofold as T  # noqa: E402

RNG = np.random.default_rng(20260922)

# pre-registered acceptance thresholds, quoted
SHIPS_DLDDT = 0.020
SHIPS_P = 0.05
SHIPS_SPECIFICITY = 0.010
SHIPS_DROT = -5.0
FG_MOVED_A = 0.50
PERIPHERAL_MIN_N = 5


def boot_ci(x, n=10000):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if len(x) < 3:
        return (np.nan, np.nan)
    d = RNG.choice(x, (n, len(x)), replace=True).mean(axis=1)
    return (float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)))


def paired(df, col, arm_a="single", arm_b="double"):
    """Per-ligand mean within each arm, then the paired difference."""
    p = df.pivot_table(index="ligand", columns="arm", values=col, aggfunc="mean")
    if arm_a not in p or arm_b not in p:
        return None
    d = (p[arm_b] - p[arm_a]).dropna()
    if len(d) < 3:
        return None
    try:
        w = stats.wilcoxon(d, zero_method="wilcox")
        pval = float(w.pvalue)
    except ValueError:
        pval = 1.0
    lo, hi = boot_ci(d.values)
    return {"metric": col, "arm": arm_b, "n": int(len(d)),
            "mean_a": round(float(p[arm_a].mean()), 4),
            "mean_b": round(float(p[arm_b].mean()), 4),
            "delta": round(float(d.mean()), 4),
            "median_delta": round(float(d.median()), 4),
            "ci95": [round(lo, 4), round(hi, 4)],
            "wilcoxon_p": round(pval, 5),
            # sign counts, NOT "better/worse": for the RMSD and angle columns a
            # positive difference is worse, and naming them better/worse inverted
            # the reading of half the table.
            "n_pos": int((d > 0).sum()), "n_neg": int((d < 0).sum()),
            "tied": int((d == 0).sum())}


# ---------------------------------------------------------------------------
# F/G span: compare the two PREDICTIONS to each other
# ---------------------------------------------------------------------------

def fg_between_arms(engine: str, ligands: pd.DataFrame) -> pd.DataFrame:
    from cypstruct import pose as P

    out = T.OUT_ROOT / engine
    rows = []
    for r in ligands.itertuples():
        ref = T.load_crystal(r.pdb, r.id)
        if ref is None:
            continue
        pocket = set(P.pocket_residues_from_structure(ref, radius=8.0))
        loaded = {}
        for arm in T.ARMS:
            f = sorted(out.glob(f"{r.id}__{arm}__r0*.cif"))
            if not f:
                continue
            cx = P.load_structure(f[0])
            cx, off, ident = T.renumber_to_reference(cx, ref)
            loaded[arm] = cx
        if "single" not in loaded:
            continue
        base = loaded["single"]
        row = {"ligand": r.id}
        # each arm's own distance to the crystal over 210-216
        for arm, cx in loaded.items():
            try:
                al, _fit, _n = P.align_by_residue(cx, ref, resnums=sorted(pocket))
                v, n_at = T.backbone_rmsd(al, ref, set(T.FG_SPAN))
            except ValueError:
                v, n_at = np.nan, 0
            row[f"fg_to_crystal_{arm}"] = v
            row["n_fg_atoms"] = n_at
        # arm-to-arm, superposed on the shared pocket CA set (the pre-registered rule),
        # and on the pocket set with 210-216 REMOVED (a declared sensitivity, because
        # including the measured span in the fit absorbs part of the motion)
        for arm in ("double", "decoy"):
            if arm not in loaded:
                continue
            for tag, resset in (("", sorted(pocket)),
                                ("_excl", sorted(pocket - set(T.FG_SPAN)))):
                if len(resset) < 3:
                    row[f"fg_single_vs_{arm}{tag}"] = np.nan
                    continue
                try:
                    al, _f, _n = P.align_by_residue(loaded[arm], base, resnums=resset)
                    v, _na = T.backbone_rmsd(al, base, set(T.FG_SPAN))
                except ValueError:
                    v = np.nan
                row[f"fg_single_vs_{arm}{tag}"] = v
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# pool effects: ORACLE first, then selection
# ---------------------------------------------------------------------------

def pool_effects(df: pd.DataFrame, n_draws: int = 4000) -> dict:
    def block(sub, name):
        sub = sub.dropna(subset=["lddt_pli"])
        g = sub.groupby("ligand")
        oracle = float(g.lddt_pli.max().mean())
        rnd = float(sub.lddt_pli.mean())
        sel = np.nan
        x = sub.dropna(subset=["xeng"])
        if len(x):
            # `xengine.select` is argmax of -z(xeng) within the ligand. z-scoring within
            # a group is a monotone decreasing map of xeng, so the choice is argmin(xeng)
            # — identical, and defined for a ligand that has only one pose, where the
            # z-score is 0/0. The equivalence is asserted below rather than assumed.
            picked = x.loc[x.groupby("ligand").xeng.idxmin()]
            from cypstruct import xengine as X
            multi = x.groupby("ligand").filter(lambda g: len(g) > 1)
            if len(multi):
                ship = X.select(multi, ligand_col="ligand", xeng_col="xeng")
                mine = picked[picked.ligand.isin(ship.ligand)]
                assert set(ship.file) == set(mine.file), "select() disagreement"
            sel = float(picked.lddt_pli.mean())
        # random-draw null over the same pool
        draws = []
        gs = {k: v.lddt_pli.values for k, v in g}
        for _ in range(n_draws):
            draws.append(np.mean([RNG.choice(v) for v in gs.values()]))
        draws = np.array(draws)
        return {"pool": name, "n_ligands": int(sub.ligand.nunique()),
                "poses_per_ligand": round(float(g.size().mean()), 2),
                "oracle": round(oracle, 4), "random_mean": round(rnd, 4),
                "selected": round(sel, 4) if sel == sel else None,
                "gain_vs_random": (round(sel - rnd, 4) if sel == sel else None),
                "null_p99": round(float(np.percentile(draws, 99) - rnd), 4)}

    a = df[df.arm == "single"]
    ab = df[df.arm.isin(["single", "double"])]
    return {"arm_A_only": block(a, "arm A (single) only"),
            "arm_A_plus_B": block(ab, "arm A + arm B")}


def main(engine: str = "protenix_v2") -> dict:
    src = DATA_PROCESSED / f"two_ligand_poses_{engine}.csv"
    df = pd.read_csv(src)
    ligands, counts = T.pilot_ligands()

    # ---- controls -------------------------------------------------------
    c1 = {"min_renumber_identity": round(float(df.renumber_identity.min()), 4),
          "renumber_offsets": sorted(df.renumber_offset.unique().tolist()),
          "n_exact_zero_lddt": int((df.lddt_pli == 0.0).sum()),
          "exact_zero_with_bisy_gt_10": int(
              ((df.lddt_pli == 0.0) & (df.bisy_rmsd > 10)).sum())}
    # C2 determinism: distinct md5 per (ligand, arm) across replicates
    g = df.groupby(["ligand", "arm"]).md5.nunique()
    c2 = {"groups": int(len(g)), "distinct_files_median": float(g.median()),
          "groups_with_1_distinct": int((g == 1).sum()),
          "groups_with_2_distinct": int((g >= 2).sum()),
          "lddt_sd_within_group": round(float(
              df.groupby(["ligand", "arm"]).lddt_pli.std().mean()), 6)}
    two = df[df.arm.isin(["double", "decoy"])]
    c3 = {"rows": int(len(two)),
          "with_two_entities": int((two.n_entities == 2).sum()),
          "double_with_two_query_copies": int(
              (df[df.arm == "double"].n_query_copies == 2).sum()),
          "double_rows": int((df.arm == "double").sum()),
          "inter_copy_min_p5": round(float(np.nanpercentile(
              df[df.arm == "double"].inter_copy_min_dist, 5)), 2),
          "inter_copy_min_median": round(float(np.nanmedian(
              df[df.arm == "double"].inter_copy_min_dist)), 2),
          "inter_copy_below_2.5A": int(
              (df[df.arm == "double"].inter_copy_min_dist < 2.5).sum()),
          "crystal_floor_A": 2.90}
    c5 = {"mapped": int(df.mapped.sum()), "rows": int(len(df)),
          "mapping_rate": round(float(df.mapped.mean()), 4)}
    # the two pre-registered copy-choice rules agree or they do not; if they agree
    # everywhere the "ceiling" row below is identical to the primary, which is a
    # result about the prediction and not a copy-paste error (too-clean-numbers rule)
    dd = df[df.arm == "double"]
    c6 = {"double_rows": int(len(dd)),
          "nearest_fe_copy_is_also_best_match": int(dd.same_copy_fe_and_best.sum()),
          "max_lddt_gap_between_the_two_rules": round(float(
              (dd.lddt_pli_bestmatch - dd.lddt_pli).max()), 6)}

    # ---- primary --------------------------------------------------------
    metrics = ["lddt_pli", "bisy_rmsd", "rot_deg", "trans_a", "conformer_a",
               "pocket_ca_rmsd", "raw_rmsd"]
    primary = {"double": [paired(df, m, "single", "double") for m in metrics],
               "decoy": [paired(df, m, "single", "decoy") for m in metrics]}
    primary = {k: [r for r in v if r] for k, v in primary.items()}

    # README trap 1: count DISTINCT poses, never jobs. Four nominal replicates are not
    # four opinions here, so every paired number is repeated on the deduplicated pool.
    ded = df.drop_duplicates(["ligand", "arm", "md5"])
    depth = df.groupby(["ligand", "arm"]).md5.nunique()
    dedup_note = {"nominal_reps_per_cell": 4,
                  "distinct_poses_per_cell_mean": round(float(depth.mean()), 2),
                  "distinct_poses_per_cell_min": int(depth.min()),
                  "distinct_poses_per_cell_max": int(depth.max()),
                  "rows_before": int(len(df)), "rows_after": int(len(ded))}
    dedup_primary = {
        "double": [paired(ded, m, "single", "double") for m in metrics],
        "decoy": [paired(ded, m, "single", "decoy") for m in metrics]}
    dedup_primary = {k: [r for r in v if r] for k, v in dedup_primary.items()}

    # secondary: the 2-way oracle over copies (ceiling, never a result)
    ceiling = paired(df, "lddt_pli_bestmatch", "single", "double")

    # specificity: arm B vs arm C directly, paired
    spec = {m: paired(df, m, "decoy", "double") for m in ("lddt_pli", "rot_deg")}

    # ---- localisation ---------------------------------------------------
    d2 = df[df.arm == "double"].groupby("ligand").agg(
        second_class=("second_class", lambda s: s.mode().iat[0]),
        second_fe=("second_fe_dist", "mean"),
        per_contacts=("second_peripheral_contacts", "mean"))
    loc = {"counts": d2.second_class.value_counts().to_dict(),
           "n_ligands": int(len(d2)),
           "median_second_fe_dist": round(float(d2.second_fe.median()), 2),
           "peripheral_called": bool(
               (d2.second_class == "peripheral_groove").sum() >= PERIPHERAL_MIN_N),
           "peripheral_bar": PERIPHERAL_MIN_N}
    dc = df[df.arm == "decoy"].groupby("ligand").second_class.agg(
        lambda s: s.mode().iat[0])
    loc["decoy_counts"] = dc.value_counts().to_dict()

    # pre-registered follow-up: where the second copy is in the ACTIVE SITE, does the
    # scored copy's orientation get better or worse?
    piv_l = df.pivot_table(index="ligand", columns="arm", values="lddt_pli",
                           aggfunc="mean")
    piv_r = df.pivot_table(index="ligand", columns="arm", values="rot_deg",
                           aggfunc="mean")
    sub = {}
    for cls in d2.second_class.unique():
        ids = d2.index[d2.second_class == cls]
        sub[cls] = {
            "n": int(len(ids)),
            "d_lddt": round(float((piv_l.loc[ids, "double"]
                                   - piv_l.loc[ids, "single"]).mean()), 4),
            "d_rot_deg": round(float((piv_r.loc[ids, "double"]
                                      - piv_r.loc[ids, "single"]).mean()), 2)}

    # per-ligand table for the write-up, plus the WITHIN-ligand spread that decides
    # whether any of these paired differences can be read at all
    per = df.pivot_table(index="ligand", columns="arm",
                         values=["lddt_pli", "rot_deg", "bisy_rmsd"],
                         aggfunc="mean").round(4)
    spread = df.groupby(["ligand", "arm"]).agg(
        n=("lddt_pli", "size"), sd=("lddt_pli", "std"),
        rot_sd=("rot_deg", "std")).reset_index()
    per.columns = ["_".join(c) for c in per.columns]
    per = per.join(d2)
    per.to_csv(DATA_PROCESSED / f"two_ligand_per_ligand_{engine}.csv")
    spread.to_csv(DATA_PROCESSED / f"two_ligand_within_ligand_{engine}.csv", index=False)
    within = {"mean_within_cell_sd_lddt": round(float(spread.sd.mean()), 4),
              "mean_within_cell_sd_rot_deg": round(float(spread.rot_sd.mean()), 2),
              "replicates_per_cell_median": float(spread.n.median()),
              "note": ("the paired delta must be read against this, not against the "
                       "between-ligand spread")}

    # ---- F/G ------------------------------------------------------------
    fg = fg_between_arms(engine, ligands)
    fg.to_csv(DATA_PROCESSED / f"two_ligand_fg_{engine}.csv", index=False)
    fgs = {}
    for c in fg.columns:
        if c == "ligand":
            continue
        v = pd.to_numeric(fg[c], errors="coerce").dropna()
        if len(v):
            fgs[c] = {"n": int(len(v)), "median": round(float(v.median()), 3),
                      "mean": round(float(v.mean()), 3), "max": round(float(v.max()), 3)}
    fg_moved = None
    if "fg_single_vs_double" in fg and "fg_single_vs_decoy" in fg:
        md = float(pd.to_numeric(fg.fg_single_vs_double, errors="coerce").median())
        mc = float(pd.to_numeric(fg.fg_single_vs_decoy, errors="coerce").median())
        fg_moved = bool(md >= FG_MOVED_A and md > mc)
        fgs["verdict_inputs"] = {"median_single_vs_double": round(md, 3),
                                 "median_single_vs_decoy": round(mc, 3),
                                 "bar": FG_MOVED_A, "MOVED": fg_moved}

    # ---- pools ----------------------------------------------------------
    # on the DEDUPLICATED pool, because an oracle over duplicates is an oracle over
    # one pose written twice (FINDING 015)
    pools = pool_effects(ded)
    pools["_raw_with_duplicates"] = pool_effects(df)

    # ---- verdict --------------------------------------------------------
    pl = {r["metric"]: r for r in primary["double"]}
    pc = {r["metric"]: r for r in primary["decoy"]}
    dB, dC = pl["lddt_pli"]["delta"], pc["lddt_pli"]["delta"]
    ships = (dB >= SHIPS_DLDDT and pl["lddt_pli"]["wilcoxon_p"] < SHIPS_P
             and (dB - dC) >= SHIPS_SPECIFICITY
             and pl["rot_deg"]["delta"] <= SHIPS_DROT)
    refuted = ((dB <= 0 and pl["lddt_pli"]["ci95"][1] < SHIPS_DLDDT) or (dC >= dB))
    verdict = "SHIPS" if ships else ("REFUTED" if refuted else "MEASURED")

    out = {"engine": engine, "pilot_counts": counts,
           "controls": {"C1_numbering": c1, "C2_determinism": c2,
                        "C3_two_copies": c3, "C5_symmetry_mapping": c5,
                        "C6_copy_choice_rules_agree": c6},
           "within_ligand_spread": within,
           "pose_depth_after_dedup": dedup_note,
           "primary_paired": primary,
           "primary_paired_deduplicated": dedup_primary,
           "ceiling_bestmatch_copy": ceiling,
           "specificity_double_vs_decoy": spec,
           "second_copy_localisation": loc, "by_localisation_class": sub,
           "fg_span_210_216": fgs, "pools_oracle_then_selection": pools,
           "acceptance": {"delta_B_lddt": dB, "delta_C_lddt": dC,
                          "wilcoxon_p": pl["lddt_pli"]["wilcoxon_p"],
                          "delta_B_rot_deg": pl["rot_deg"]["delta"],
                          "SHIPS": ships, "REFUTED": refuted,
                          "FG_MOVED": fg_moved,
                          "PERIPHERAL_OBSERVED": loc["peripheral_called"]},
           "verdict": verdict}
    (DATA_PROCESSED / f"two_ligand_analysis_{engine}.json").write_text(
        json.dumps(out, indent=1, default=str))
    return out


if __name__ == "__main__":
    print(json.dumps(main(sys.argv[1] if len(sys.argv) > 1 else "protenix_v2"),
                     indent=1, default=str))
