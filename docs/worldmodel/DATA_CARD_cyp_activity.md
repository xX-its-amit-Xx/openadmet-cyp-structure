# DATA CARD — `cyp_activity.parquet`

A stratified, assay-traceable activity table for CYP3A4 and the four other major
drug-metabolising human P450s, built 2026-09-20 from the verified raw corpus at
`/scratch/shenoy.am/worldmodel/raw/`.

| | |
|---|---|
| **Path** | `/scratch/shenoy.am/worldmodel/derived/cyp_activity.parquet` (zstd, 6.5 MB) |
| **Rows** | **261,164** |
| **Columns** | 45 |
| **Unique compounds** (standard InChIKey) | 51,051 |
| **Build scripts** | `/scratch/shenoy.am/worldmodel/scripts/{cypcommon,chembl_extract,bdb_extract,build_cyp_activity,qc_cyp_activity}.py` |
| **Machine-readable stats** | `/scratch/shenoy.am/worldmodel/derived/cyp_activity_stats.json` |
| **Range-violation dump** | `/scratch/shenoy.am/worldmodel/derived/cyp_activity_range_violations.tsv` |
| **Licence** | mixed, per row — see `source_licence`. CC BY-SA 3.0 (ChEMBL) is viral and dominates 70% of rows. |

The table is **not** aggregated. One row is one reported measurement, and every row keeps
`assay_id`, `assay_description`, `confidence_score`, `standard_type`, `standard_relation`
and the source document. Nothing was dropped for being ugly; it was flagged.

---

## 1. Counts

### By source

| `source` | Rows | Licence |
|---|---:|---|
| `chembl` (ChEMBL 37) | 183,577 | CC BY-SA 3.0 |
| `openadmet_challenge` | 46,110 | Apache-2.0 |
| `bindingdb` (202609) | 30,393 | CC BY 4.0 |
| `openadmet_octant` | 1,084 | CC BY 4.0 |

### By target

| Target | ChEMBL id | UniProt | Rows | Unique compounds | Unique assays |
|---|---|---|---:|---:|---:|
| **CYP3A4** | CHEMBL340 | P08684 | **83,785** | 45,571 | 7,501 |
| CYP2D6 | CHEMBL289 | P10635 | 53,253 | 34,986 | 5,105 |
| CYP2C9 | CHEMBL3397 | P11712 | 48,702 | 33,220 | 4,893 |
| CYP1A2 | CHEMBL3356 | P05177 | 41,933 | 29,979 | 3,794 |
| CYP2C19 | CHEMBL3622 | P33261 | 33,491 | 24,833 | 3,723 |

Source × target (ChEMBL / BindingDB / challenge / Octant):
CYP3A4 56,367 / 10,122 / 16,212 / 1,084 · CYP2D6 36,234 / 6,663 / 10,356 / — ·
CYP2C9 33,399 / 5,787 / 9,516 / — · CYP1A2 27,652 / 4,255 / 10,026 / — ·
CYP2C19 29,925 / 3,566 / — / — (the OpenADMET challenge does not cover CYP2C19).

### By measurement class (`measurement_class`)

`potency` 146,541 · `other` 61,548 · `percent_effect` 20,858 · `single_concentration`
17,504 · `effect_size` 14,303 · `kinetic_or_other` 410.

Top `standard_type`: AC50 83,967 · IC50 70,476 · Potency 23,977 · log2fc_single_conc
17,504 · Inhibition 16,048 · pIC50 15,392 · Emax 14,303 · Ki 5,688.

**`standard_type` is not a modality.** Of the 83,967 AC50 rows, effectively all are the
PubChem qHTS cytochrome panel re-deposited into ChEMBL (`PUBCHEM_BIOASSAY: Cytochrome
panel assay …`, ~16.8 k rows per CYP). ChEMBL files the same qHTS measurement under both
`AC50` and `Potency`. Treating AC50 and IC50 as one column silently merges a 2009 qHTS
screen with three decades of hand-curated DMPK literature.

---

## 2. ChEMBL `confidence_score` — which subset is safe

| Score | Rows | Meaning |
|---|---:|---|
| **9** | **77,824** (51,936 of them `potency`) | direct single protein — **this is the safe subset** |
| 8 | 105,753 | homologous single protein |

Nothing below 8 survives, because the five target IDs are themselves single-protein
targets. **Filter `source == 'chembl' and confidence_score == 9` for anything that claims
to be a measurement on one enzyme.** For CYP3A4 specifically that leaves 9,894 ChEMBL
IC50 rows. Boltz-2's affinity head paid a 75% loss to this filter; here it is 58%
(77,824 of 183,577 ChEMBL rows kept).

`confidence_score` is null for BindingDB and OpenADMET rows — those sources do not
express the concept. Do not fill it.

---

## 3. Censored values — flagged, never silently converted

`censored` (bool) and `censor_direction` (`right` / `left`) are derived from
`standard_relation`, which is preserved verbatim.

| | Rows censored | Fraction |
|---|---:|---:|
| **All rows** | 36,922 | **14.1%** |
| Potency rows only | — | **22.9%** |
| ChEMBL | 21,611 | 11.8% (18.2% of ChEMBL potency) |
| **BindingDB** | 15,311 | **50.4%** |
| OpenADMET challenge / Octant | 0 | 0% |

Direction: 33,819 right-censored (`>`, `>=`), 3,103 left-censored.

**BindingDB is half bounds.** Half of its CYP rows are `>10000 nM` style non-results.
A loader that strips the `>` turns 15,311 "did not inhibit" statements into a phantom
mode at pIC50 = 5.0 — the Davis pathology, at scale, inside our own corpus.

**Zero censored rows in OpenADMET is a false comfort — see §6.**

---

## 4. Probe substrate — how far it could be recovered

`probe`, `probe_class` (clinical_drug / fluorogenic / luminogenic), `probe_all` (every
probe named in the text) and `probe_status`.

| `probe_status` | Rows |
|---|---:|
| `canonical` — a known probe for *this* enzyme | **29,412** |
| `sole_noncanon` — one probe named, not canonical here | 659 |
| `ambiguous_panel` — several probes named, none canonical here | 3,973 |
| `none` — no probe named | 227,120 |

**30,071 rows (11.5%) carry a probe assignment. 231,093 (88.5%) do not.**

By source: **ChEMBL 30,071 of 183,577 (16.4%). BindingDB 0 of 30,393. OpenADMET 0 of
47,194.** By target: CYP3A4 12,520 (14.9%), CYP1A2 5,299 (12.6%), CYP2C19 4,202 (12.5%),
CYP2D6 4,362 (8.2%), CYP2C9 3,688 (7.6%).

Top CYP3A4 probes: BFC 3,650 · midazolam 3,650 · testosterone 2,634 · luciferin 789 ·
BQ 652 · fluorescein 404 · DBF 245 · BzRes 84 · nifedipine 78 · erythromycin 77.
(BFC and midazolam landing on the same integer is coincidence, checked: 156 vs 879
distinct assays. 1,742 of the BFC rows are one DRUGMATRIX panel.)

Two things the assignment had to handle, and one it cannot:

- **Panel descriptions name several probes.** 13,843 rows name more than one, because
  one description covers a CYP panel ("CYP3A4 … BFC … CYP2D6 … MAMC"). The same text is
  attached to rows on different targets. Assignment is therefore **target-aware**: the
  probe chosen is the one canonical for that row's enzyme. A first-match parser assigns
  midazolam to CYP2D6 rows — an earlier pass of this build did exactly that.
- **Generic tokens must lose.** "coumarin" and "fluorescein" are matched last, so
  "7-benzyloxy-4-trifluoromethyl coumarin" resolves to BFC.
- **BindingDB carries no protocol at all.** Its only description field is the target name
  ("Cytochrome P450 3A4"). Probe recovery there is not hard, it is impossible — the
  information is not in the file. This is the strongest argument in the corpus for not
  pooling BindingDB IC50 with anything.

---

## 5. The trap, quantified: within-compound cross-probe pIC50 spread

CYP3A4, `standard_type == IC50`, uncensored, in range, probe canonical for CYP3A4.
Per (compound, probe) median, then max − min across probes within a compound.

| | |
|---|---:|
| Rows used | 3,274 |
| Compounds with ≥ 2 distinct probes | **304** |
| Spread p25 / **median** / p75 | 0.10 / **0.31** / 0.61 |
| Spread p90 / max | 0.98 / 2.96 |
| Fraction > 0.5 log | 33.2% |
| Fraction > 1.0 log | **8.9%** |
| Fraction > 2.0 log | 0.3% |
| Fraction exactly 0.00 | 7.2% |

Head-to-head, midazolam vs testosterone (n = 146): median signed Δ **−0.003**, median
|Δ| **0.243**, 4.1% differ by more than a log, 8.2% are byte-identical.

Across all five CYPs (332 compound-target pairs): median 0.35, p90 1.05, 12.3% > 1 log.

**This revises the received claim.** `DATA_AFFINITY.md` §7 says CYP3A4 IC50 against
midazolam, testosterone and a fluorogenic probe differ "by more than an order of
magnitude on the same compound". On the measurable subset that is the **tail, not the
centre**: the typical cross-probe disagreement is **0.31 log units**, and midazolam vs
testosterone is centred on **zero with a 0.24-log typical deviation**. Order-of-magnitude
disagreement is real (the max is 2.96 logs) but it is ~9% of compounds, not the rule.

Three reasons not to read that as reassurance:

1. **n = 304 of 45,571 CYP3A4 compounds.** The measurable subset is 0.7% of the target,
   and it is the subset someone chose to run twice — enriched for well-behaved,
   literature-notable compounds.
2. **7.2% of "cross-probe" spreads are exactly zero**, which is not a measurement result.
   It means the same number was filed under two differently-worded assays.
3. The centred midazolam–testosterone comparison says the two *clinical* probes agree on
   average. It says nothing about a fluorogenic probe against a clinical one, where the
   pairs are too thin (BFC vs BzRes n = 25, 24% > 1 log; BFC vs BQ n = 19) to carry weight.

**Use as the label-noise floor: ~0.3 log units typical, ~1.0 at p90.** A model whose
held-out RMSE on pooled CYP3A4 IC50 is below ~0.5 is fitting deposition artefacts.

---

## 6. Replicates, and why the naive number is wrong

Deduplication was on (`source`, `compound_id`, `target_name`, `assay_id`,
`standard_type`, `standard_relation`, `standard_value`, `standard_units`, `pactivity`).
**3,764 exact duplicates removed** (ChEMBL 2,599, BindingDB 1,165). Genuine replicates
were kept.

| Grouping | Pairs with ≥2 rows | Median range | p90 | > 1 log |
|---|---:|---:|---:|---:|
| (InChIKey, target, standard_type) | 23,471 | **0.00** | 0.49 | 2.9% |
| CYP3A4 IC50 only | 2,966 | **0.00** | 0.65 | 5.5% |
| ChEMBL only | 8,923 | **0.00** | 0.50 | 3.6% |
| **CYP3A4 IC50, same probe, different assay** | 235 | **0.28** | 1.30 | **13.2%** |
| **ChEMBL CYP3A4 IC50, different documents** | 235 | 0.00 | 1.42 | 13.6% |

A median range of exactly 0.00 is not agreement between labs. **57% of ChEMBL CYP3A4
IC50 compound–target pairs that appear in two or more different documents have identical
values**, and 59.8% of all replicated ChEMBL pairs come from a single document. The
corpus's "replicates" are dominated by one measurement re-deposited — patent plus paper,
BindingDB into ChEMBL, qHTS into both.

**Report `replicate_spread_CYP3A4_IC50_same_probe` (median 0.28, p90 1.30, 13.2% > 1 log)
as the replicate noise. Do not report the 0.00.**

---

## 7. Physical range — inspected, not dropped

Assertion: pActivity ∈ [2, 13]. `range_flag` ∈ {`ok`, `below_2`, `above_13`,
`no_pactivity`}. **Nothing was deleted.** All violations are in
`cyp_activity_range_violations.tsv`.

| | Rows |
|---|---:|
| `ok` | 146,137 |
| `below_2` | **399** |
| `above_13` | 0 |
| `no_pactivity` (non-concentration type, or unconvertible units) | 114,628 |

Breakdown of the 399: ChEMBL AC50 105, ChEMBL IC50 79, ChEMBL Ki 16, ChEMBL Km 1,
BindingDB IC50 33, BindingDB Ki 8, **OpenADMET challenge pIC50 131, Octant pIC50 26**.

What inspection found:

- **200 of the 399 already carry a ChEMBL `data_validity_comment`** (199 "Outside typical
  range", 1 "Potential transcription error"). ChEMBL flagged them; most pipelines never
  read that column.
- **The unit-confusion signature is present.** 142 violations have an in-range sibling
  measurement of the same (compound, target, type); **44 of those sit 3.0 ± 0.7 log units
  away** — the mM/nM 1000× error, sign-preserving, invisible to a wide outlier filter.
- The worst ChEMBL Ki rows are absurd on their face (5×10¹⁰ to 7×10¹⁴ nM, i.e. pKi −1.7
  to −5.8) and are almost certainly wrong-unit or wrong-column deposits.
- **The 157 OpenADMET violations are not errors — they are the assay floor**, see below.

---

## 8. OpenADMET: the closest thing to our evaluation distribution

`source == 'openadmet_challenge'` (Apache-2.0) and `source == 'openadmet_octant'`
(CC BY 4.0), both clearly separated by the `source` column.

- 46,110 challenge rows: direct-inhibition pIC50 (4 CYPs), TDI-condition pIC50,
  17,504 single-concentration log2 fold-changes, 14,303 Emax values.
- 1,084 Octant CYP3A4 pIC50 with SE.
- The 750-compound blinded test set is **deliberately not in this table** — it has no
  labels and including it invites leakage.

Three warnings:

1. **Neither release states its probe substrate.** The Octant README states a
   fluorescence readout; the challenge release states neither probe nor readout. Both
   have `probe = null`; `probe_class = 'fluorogenic'` is set only for Octant, on the
   README's authority. Do not assume the challenge assay is midazolam-based.
2. **The Octant IC50 is not a reversible IC50.** It is measured after a 30-minute
   active-enzyme pre-incubation and therefore folds in time-dependent inhibition
   (`preincubation = True`). It is not comparable to a literature reversible IC50, and
   OpenADMET says so.
3. **The pIC50 floor is soft censoring with no relation flag.** Challenge direct-inhibition
   pIC50 has min 1.91, p05 2.25, and **24% of values below 4** — a fitted tail below the
   lowest tested concentration, filed as `standard_relation = '='`. These are bounds
   wearing the costume of measurements. `censored` is `False` for all 46,110 OpenADMET
   rows because the release carries no relation column; that zero in §3 means "not
   recorded", not "not censored".

**Overlap with the literature is almost nil, and where it exists it disagrees.** Only
**287 of the challenge's 6,141 InChIKey skeletons appear anywhere in ChEMBL** (78 in
BindingDB; Octant 50). On the 75 overlapping (compound, target) pairs:
Pearson r **0.57**, Spearman ρ **0.50**, **RMSE 0.78**, median signed Δ **−0.37**
(OpenADMET reads ~0.37 log less potent), 17.3% differ by more than a log.

That is the honest transfer estimate from public literature data to the challenge's own
assay, and it is barely better than the cross-probe spread in §5. Octant and the
challenge share 1,075 of 1,084 skeletons — they are one platform, not two datasets.

---

## 9. Schema

Identity — `source`, `source_licence`, `record_id`, `target_name`, `target_chembl_id`,
`uniprot`, `compound_id`, `smiles`, `inchikey`.

Measurement — `standard_type`, `measurement_class`, `standard_relation`, `censored`,
`censor_direction`, `standard_value`, `standard_units`, `concentration_M`, `pactivity`,
`pactivity_units` (`-log10(M)`), `value_error`, `pchembl_value`, `range_flag`.

Assay provenance — `assay_id`, `assay_chembl_id`, `assay_description`, `assay_type`,
`assay_test_type`, `confidence_score`, `assay_organism`, `bao_format`, `probe`,
`probe_class`, `probe_all`, `probe_status`, `probe_assigned`, `assay_system`
(HLM / hepatocytes / recombinant / cell-based), `preincubation`.

Document — `doc_id`, `doi`, `pubmed_id`, `year`, `journal`.

QC — `data_validity_comment`, `activity_comment`, `potential_duplicate`.

---

## 10. Recipes

```python
import pandas as pd
df = pd.read_parquet('/scratch/shenoy.am/worldmodel/derived/cyp_activity.parquet')

# strictest literature subset: single protein, uncensored, physical, probe known
strict = df[(df.source == 'chembl') & (df.confidence_score == 9) &
            (df.standard_type == 'IC50') & (~df.censored) &
            (df.range_flag == 'ok') & (df.probe_status == 'canonical')]

# never pool probes; model probe as a fixed effect or split on it
strict.groupby(['target_name', 'probe']).size()

# censored rows are bounds — keep them, with a censored likelihood
bounds = df[df.censored & (df.range_flag == 'ok')]

# evaluation-distribution subset
oa = df[df.source.str.startswith('openadmet')]
```

**Do not** merge ChEMBL with BindingDB by summing: 10,875 of BindingDB's 11,714
skeletons are already in ChEMBL. **Do not** average `pactivity` across `assay_id`.
**Do not** use `standard_value` without `standard_units` and `standard_relation`.

---

## 11. Known gaps

- **Papyrus, BRENDA, BioLiP2, HiQBind, CatPred-DB, SKEMPI are not in this table.** Papyrus
  is pre-aggregated ChEMBL 34 (the aggregation destroys exactly the spread §5 measures);
  BRENDA/CatPred carry kcat/Km for P450s but overwhelmingly on surrogate dye substrates;
  HiQBind/BioLiP2 contribute a low-tens number of CYP complexes with a matched affinity.
  Each is a separate build if wanted.
- **PubChem AID 1851 is in here only via ChEMBL's re-deposit**, not from PubChem directly.
- `assay_system` is parsed from free text and is null wherever the text is silent.
- No compound standardisation beyond source-supplied SMILES; InChIKeys are source-supplied
  for ChEMBL/BindingDB and RDKit-computed for OpenADMET. Salts and tautomers are not
  reconciled — join on the 14-character skeleton when that matters.

## 12. Build notes

Everything ran CPU-only on Explorer under Slurm (`short` partition, 24–64 GB), never on
the login node: `unzip`/`tar` of the 9 GB BindingDB TSV and 30 GB ChEMBL SQLite were
batch jobs, per the §Traps warning that the login-node cgroup SIGKILLs long reads while
reporting success. Venv at `/scratch/shenoy.am/worldmodel/venv` (needs
`module load python/3.13.5` for `libpython3.13.so`); `/scratch/shenoy.am/cyp-finetune/env`
was not touched. No bulk data was written to the Windows box.

One bug worth recording: `Series.astype(str)` on a pandas nullable `string`/`boolean`
dtype **keeps `pd.NA`**, so a concatenated dedup key silently became NA for every row with
one missing field — and `duplicated()` treats those as equal. The first build "removed
157,755 duplicates" and deleted every OpenADMET row. The tell was the number being too
large, not an error.
