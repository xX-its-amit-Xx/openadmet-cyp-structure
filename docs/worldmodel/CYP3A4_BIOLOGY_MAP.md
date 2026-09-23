# CYP3A4 — the biology map, gene to pathway to cell type

**Written 2026-09-22 for the OpenADMET CYP3A4 structure track.** This is the companion to
[`CYP3A4_EVOLUTION.md`](CYP3A4_EVOLUTION.md), which asks *why CYP3A4 is the way it is* from
phylogeny, selection and endogenous chemistry. This document asks a different question:
**what is the object, mechanistically, from DNA to protein complex to cell type**, and which
parts of that map are load-bearing for predicting a ligand pose.

It **extends** the evolution document and does not repeat it. Where the two overlap — the CYP3A
cluster's gene order, the induction machinery, zonation, multiple occupancy, the heme — the
evolution doc's section is cited and only *new* material is added here. The machine-readable form of
everything below is
[`data/processed/cyp3a4_biology_graph.json`](../../data/processed/cyp3a4_biology_graph.json).

**How to read it.** Every non-obvious claim carries a URL, DOI, PDB ID or database accession.
Four status marks are used throughout and they are meant literally:

| mark | meaning |
|---|---|
| **ESTABLISHED** | multiple independent sources, or one primary measurement with a method I can name |
| **CONTESTED** | the literature genuinely disagrees. Both sides are given; no side is picked |
| **[INFERENCE]** | my reasoning, not a reported result |
| **NO EVIDENCE FOUND** | I looked and found nothing. Where the search itself was the limit, it says *(search-limited)* — that means unsearched, **not** negative |

These four marks also exist as a machine-readable `status` field on **every node and every edge** of
the companion JSON, so the grading is queryable at the point of use and not only readable here (§10).

Several of the things this document was asked to cover turn out to have little or no literature —
microproteins, lncRNAs at the locus, most PTMs, chaperones. Those are reported as gaps. A gap
correctly reported is worth more here than a confident guess.

Each section ends with a short **"what this licenses for the structure track"** line: concretely,
which partners are worth co-folding with, which compounds represent the liver substrate space, and
which paralogs are worth docking against.

**Two operational notes.** (1) The session's WebSearch budget was exhausted early, so all literature
work was done through direct Europe PMC, PubMed E-utilities, UniProt, Ensembl, RCSB/PDBe, Reactome,
KEGG, STRING and PubChem REST endpoints. That is good for reading known literature and weaker at
exhaustive discovery of older, non-open-access work; affected items are flagged *(search-limited)*.
(2) Several numbers below are **primary measurements made for this document** over the full PDB or
over local repo data, not literature claims; those are marked **MEASURED**.

---

---

## 1. Gene to protein

`CYP3A4_EVOLUTION.md` §0 gives the four-gene locus, the ancestral gene order and the paralog
identities. This section adds the identifiers, the coordinate arithmetic, the transcript inventory
and the *cis*-element map with positions.

**Assembly convention:** GRCh38/hg38, chr7, 1-based inclusive. CYP3A4 is on the **minus strand**, so
"upstream" (negative TSS-relative positions) means **increasing** chr7 coordinates.

**The anchor everything hangs on.** The canonical TSS is **chr7:99,784,184** (the 5′ end of
ENST00000651514 / NM_017460.6). Independent check: CYP3A4\*1B (rs2740574), universally called "−290"
in the literature, sits at chr7:99,784,473 = **289 bp** upstream. So the literature's +1 and the
RefSeq/Ensembl 5′ end agree to within 1 nt, and the conversion for a TSS-relative position −n is
`chr7 = 99,784,184 + n`. **ESTABLISHED.**

### 1.1 Identifiers and coordinates

| field | value |
|---|---|
| HGNC | **HGNC:2637**; previous symbol **CYP3A3** ([rest.genenames.org](https://rest.genenames.org/fetch/symbol/CYP3A4)) |
| NCBI Gene | **1576** ([esummary](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=1576&retmode=json)) |
| Ensembl | **ENSG00000160868** |
| UniProt | **P08684**, 503 aa (= ENSP00000498939) |
| CCDS / OMIM / UCSC | CCDS5674 / 124010 / uc064fwx.1 |
| **MANE Select** | **NM_017460.6 → NP_059488.2** (= ENST00000651514.1) |
| second RefSeq | NM_001202855.3 → NP_001189784.1; legacy NM_000776 merged into NM_017460 |
| aliases | CP33, CP34, CYP3A, CYPIIIA3, CYPIIIA4, HLP, NF-25, P450C3, P450PCN1, **VDDR3** |

**Three annotations give three slightly different spans — CONTESTED only in the trivial sense that
they differ in annotation scope:** NCBI 99,756,967–99,784,184 (**27,218 bp**, the canonical
footprint), Ensembl 99,756,954–99,784,327 (27,374 bp, including non-canonical 5′/3′ extensions),
GENCODE v26 as used by GTEx 99,756,960–99,784,265 (27,306 bp).

### 1.2 The CYP3A cluster, with distances

| gene | Ensembl / NCBI | biotype | start | end | strand |
|---|---|---|---|---|---|
| **CYP3A5** | ENSG00000106258 | coding | 99,647,442 | 99,680,043 | − |
| **CYP3A51P** (= CYP3AP1, CYP3A5P1) | ENSG00000282277 / 1578 | unprocessed pseudogene | 99,685,145 | 99,700,034 | − |
| CYP3A7-CYP3A51P readthrough | ENSG00000282301 | coding | 99,684,957 | 99,735,102 | − |
| **CYP3A7** | ENSG00000160870 | coding | 99,704,105 | 99,735,196 | − |
| **CYP3AP2** (= CYP3A5P2) | NCBI 79424, **not in Ensembl** | pseudogene | 99,748,249 | 99,752,490 | − |
| **CYP3A4** | ENSG00000160868 | coding | 99,756,954 | 99,784,327 | − |
| CYP3A137P | ENSG00000261511 | unprocessed pseudogene | 99,820,018 | 99,820,086 | + |
| **CYP3A43** | ENSG00000021461 | coding | 99,827,922 | 99,867,801 | **+** |
| CYP3A52P | ENSG00000260524 | unprocessed pseudogene | 99,872,168 | 99,872,258 | + |

Sources: [Ensembl region overlap](https://rest.ensembl.org/overlap/region/human/7:99600000-99900000?feature=gene;content-type=application/json),
[NCBI esummary 1578, 79424](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=1578,79424&retmode=json).
Nomenclature note: NCBI has retired `CYP3AP1`/`CYP3A5P1`; the current symbol is **CYP3A51P**.

Intergenic gaps (edge-to-edge, **[INFERENCE]** — arithmetic on the coordinates above): CYP3A5 →
CYP3A51P 5.1 kb; CYP3A51P → CYP3A7 4.1 kb; CYP3A7 → CYP3AP2 13.1 kb; CYP3AP2 → CYP3A4 4.5 kb;
**CYP3A4 → CYP3A43 43.6 kb**. Array span CYP3A5 start → CYP3A43 end = **220,360 bp**.
[Gellner *et al.* 2001](https://doi.org/10.1097/00008571-200103000-00002) (PMID 11266076) sequenced
the locus and reported "231 kb" with three pseudogenes; the difference is annotation-edge choice and
pseudogene naming, so quote the 231 kb figure with its 2001 attribution rather than as a current
number. **CONTESTED (trivially).**

**The orientation fact that matters.** CYP3A5, CYP3A51P, CYP3A7, CYP3AP2 and CYP3A4 are all minus
strand, tandem head-to-tail; **CYP3A43 is plus strand**. **[INFERENCE]:** that puts CYP3A4's TSS
(99,784,184) and CYP3A43's TSS (99,827,922) at the *facing* ends of the two genes, i.e. **CYP3A4 and
CYP3A43 are divergent, head-to-head, with their promoters ~43.6 kb apart, and every known CYP3A4
upstream element (XREM, CLEM4, R2, R4) lies inside that shared intergenic block.** That is the
structural explanation for the otherwise strange 4C result that deleting R4 *raises* CYP3A4 and
*lowers* CYP3A43 (§1.5).

### 1.3 Paralog co-expression — GTEx v8 median TPM

| tissue | CYP3A4 | CYP3A5 | CYP3A7 | CYP3A43 |
|---|---|---|---|---|
| **liver** | **335.3** | **155.9** | **24.5** | 4.29 |
| small intestine (terminal ileum) | 48.2 | 61.5 | 0.36 | 0.11 |
| prostate | 0.14 | **42.5** | 0.45 | 0.66 |
| stomach | 0.14 | 45.3 | — | — |
| oesophagus mucosa | 0.31 | 33.7 | 0.0 | 0.03 |
| pancreas | 2.41 | 28.4 | — | 1.47 |
| kidney cortex | 0.21 | 8.15 | 0.69 | 0.07 |

([GTEx API](https://gtexportal.org/api/v2/expression/medianGeneExpression?datasetId=gtex_v8&gencodeId=ENSG00000160868.14&format=json))

- **Liver co-expresses CYP3A4 and CYP3A5** at 335 vs 156 TPM — ESTABLISHED, and directly relevant:
  any liver microsomal measurement of "CYP3A" is a mixture of two enzymes with **measurably different
  pocket plasticity** (`CYP3A4_EVOLUTION.md` §4.3).
- **CYP3A7 is present in adult GTEx liver at 24.5 TPM.** Not negligible, and consistent with reports
  that CYP3A7 persists in some adults ([Li & Lampe 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6739124/)).
- **CYP3A43's prostate dominance is CONTESTED.** The original report is unambiguous — "the highest
  expression level of CYP3A43 mRNA is observed in the prostate"
  ([Gellner 2001](https://doi.org/10.1097/00008571-200103000-00002)) — but **GTEx v8 does not
  reproduce it**: CYP3A43 ranks liver 4.29 > testis 1.77 > pancreas 1.47 > prostate 0.66, and in
  prostate the dominant CYP3A is **CYP3A5 at 42.5 TPM**, ~64× CYP3A43. Both should be reported.

### 1.4 Exon/intron structure and the transcript inventory

Canonical ENST00000651514 = NM_017460.6:

| metric | value |
|---|---|
| exons | **13**, all coding (CDS opens in exon 1, closes in exon 13) |
| mRNA | **2,781 bp**; CDS 104–1,615 (**1,512 nt** incl. stop); protein **503 aa** |
| **5′UTR** | **103 nt** |
| **3′UTR** | **1,166 nt** |
| introns | **24,437 bp** over 12 |
| ATG / stop | chr7:99,784,081 / chr7:99,758,133 |

Exon 1 is 174 bp (chr7:99,784,011–99,784,184) and exon 13 is 1,262 bp (chr7:99,756,967–99,758,228,
of which 1,166 nt is 3′UTR). Intron 6 is 1,265 bp — the one that matters for CYP3A4\*22.

**Ensembl annotates 35 transcripts** ([REST lookup](https://rest.ensembl.org/lookup/id/ENSG00000160868?expand=1;content-type=application/json)):
**27 protein-coding**, **7 NMD**, **1 retained-intron** (ENST00000480043, 651 nt), **1 with CDS not
defined**. ⚠️ Most of the coding set is a recent long-read-derived expansion (ENST000008592xx /
ENST000009843xx); treat individual members as annotation, not as validated isoforms. **No validated
alternative first exon or alternative promoter for CYP3A4 was found — NO EVIDENCE FOUND.**

UniProt P08684 additionally records **chimeric trans-spliced transcripts joining CYP3A43 exon 1 to
CYP3A4 exons**, expressed at very low hepatic level — which is a second, independent hint that the
two divergent genes share a regulatory neighbourhood.

**CYP3A4\*22 (rs35599367) is the one splice variant with real mechanism.** ESTABLISHED.
chr7:99,768,693 (GRCh38), an intron-6 variant **191 bp** upstream of exon 7 (Wang & Sadee say 192 —
an off-by-one). The T allele roughly doubles an alternatively spliced transcript with **partial
intron 6 retention** in human liver (mean 19% vs 5% of transcripts, *P* = 0.006) but **not in small
intestine**, reproduced in a HepG2 minigene (18% vs 12%, *P* = 0.003) and absent in intestinal
LS-174T ([Wang & Sadee 2016, *Pharmacogenet Genomics* 26:40](https://doi.org/10.1097/FPC.0000000000000183),
PMC4674354). The retained segment introduces an in-frame stop **and an AATAAA inside intron 6**, and
the product lacks the heme-binding signature. It is **ENST00000480043**, whose terminal exon runs
233 nt into intron 6 — the same event, with a 233 vs 255 nt annotation-boundary difference.

**3′UTR length variants are real and functionally characterised.** ESTABLISHED. NM_017460.6 itself
annotates **two** polyA signals and sites: an AATAAA at mRNA 2,049–2,054 with a site at 2,072
(**short 3′UTR ≈ 457 nt**, annotated "major") and one at 2,746–2,751 with a site at 2,781
(**long 3′UTR 1,166 nt**). [Li, Gaedigk, Hart, Leeder & Zhong 2012, *Mol Pharmacol*
81:222](https://doi.org/10.1124/mol.111.074393) (PMC3250109) showed the **short-3′UTR transcript is
preferentially expressed in developed liver and differentiated hepatocytes, preferentially induced by
rifampicin and phenobarbital, more stable (3.7 h vs 2.2 h) and gives 2.4-fold more protein**. This is
a genuine developmental and inducible APA switch. *(No PolyA_DB/APA-atlas record retrieved —
search-limited.)*

### 1.5 The promoter and the *cis*-element map, with positions

| element | vs TSS | chr7 (GRCh38) **[INFERENCE from the anchor]** | factors | status |
|---|---|---|---|---|
| **proximal ER6** (pER6 / proximal PXRE) | **−172 / −149** | 99,784,333–99,784,356 | PXR:RXRα, **VDR:RXR**, CAR (weakly) | ESTABLISHED |
| C/EBP proximal sites | −121 / −130 | ~99,784,305–99,784,314 | C/EBPα, C/EBPβ-LAP | CONTESTED (review-sourced) |
| DR1 | −237 / −211 | 99,784,395–99,784,421 | HNF4α, antagonised by COUP-TFII | CONTESTED |
| **NFSE**, containing rs2740574 (\*1B) | −290 (TSS) / −392 (ATG) | **99,784,473** | "unspecific" nuclear proteins | CONTESTED |
| **XREM** | **−7,836 / −7,607** | **99,791,791–99,792,020** | PXR:RXRα, CAR:RXR, HNF4α | ESTABLISHED |
| — footprint **FP3** (≈ dNR1, DR3) | −7,738 / −7,715 | 99,791,899–99,791,922 | PXR motif | ESTABLISHED |
| — footprint **FP4** (≈ dNR2, ER6) | −7,698 / −7,682 | 99,791,866–99,791,882 | PXR motif | ESTABLISHED |
| **CLEM4** | **−11,400 / −10,500** | **99,794,684–99,795,584** | **HNF-1α, HNF-4α, USF1, AP-1** | ESTABLISHED |
| — CLEM4 −11,129_−11,128insTGT | | 99,795,312–99,795,313 | disrupts **USF1**, −36% enhancer activity, MAF 3.1% | ESTABLISHED |
| PPARα regions PBR-I/II/III | −2,915/−2,903; −3,062/−3,050; −7,784/−7,764; −8,816/−8,804 | see report | PPARα:RXRα at DR1-like motifs | CONTESTED (positions review-sourced) |
| distal C/EBP enhancer | −5,950 / −5,663 | 99,789,847–99,790,134 | C/EBPα, C/EBPβ-LAP/LIP | CONTESTED |
| FXR element | "345-bp element in the 5′-flanking region", position not published | — | FXR:RXR | element ESTABLISHED, position NO EVIDENCE FOUND |

**Goodwin, Hodgson & Liddle 1999** ([*Mol Pharmacol* 56:1329](https://doi.org/10.1124/mol.56.6.1329),
PMID 10570062) is the primary source for XREM, verbatim: rifampicin induction "was observed only with
the longest construct, which encompassed bases −13000 to +53… Further deletion mutants localized the
induction to bases −7836 to −7607… four protected sites (FP1, FP2, FP3, and FP4). Two of these sites,
FP3 (bases −7738 to −7715) and FP4 (bases −7698 to −7682), overlapped binding motifs for … hPXR."
⚠️ The now-standard dNR1/dNR2 labels map onto FP3/FP4 in the downstream literature, but that mapping
could not be verified against the original figure — **CONTESTED**. XREM boundaries appear in print as
−7836/−7607, −7.2/−7.8 kb, −7836/−7208 and −7800/−7600; all describe the same module.

**Matsumura 2004** ([*Mol Pharmacol* 65:326](https://doi.org/10.1124/mol.65.2.326), PMID 14742674) is
the primary source for CLEM4 and names **HNF-1α, HNF-4α, USF1 and AP-1** by gel shift, with
"essentially all sites required for maximal enhancer activity". **It does not name CAR. NO EVIDENCE
FOUND for CAR at CLEM4.**

**Sequence-level detail is thin and should not be invented.** The one sequence statement that could be
sourced is that the proximal ER6, the distal DR3 and the far-module ER6 **all share the 5′ half-site
`TGAACT`**, with the proximal module an imperfect ER6, the distal an imperfect DR3 and the far module
a **perfect ER6** ([Liu 2008, *Biochem J*](https://pmc.ncbi.nlm.nih.gov/articles/PMC4114763/)). The
full 24-/18-nt element sequences are **NO EVIDENCE FOUND** in open sources — Goodwin 1999 is
paywalled. They must be read off the paper or extracted from the genome at the coordinates above.

**The 4C map — Collins & Wang 2020** ([*Pharmacogenet Genomics* 30:107](https://doi.org/10.1097/FPC.0000000000000402),
PMC8457025), 4C + 3C + CRISPR deletion + reporter in primary human hepatocytes: four regions contact
the CYP3A4 promoter. **R2** (~10 kb upstream, overlapping CLEM4) and **R4** (~5 kb upstream of the
*CYP3A43* promoter) have demonstrated regulatory roles; R1 and R3 do not. **Deleting R4 increased
CYP3A4 and decreased CYP3A43**, "likely reflecting competitive domain interactions within the CYP3A
cluster". **rs62471956** in R4 raises transcriptional activity, associates with higher CYP3A43 and
lower CYP3A4 across 136 livers, and is in **complete LD with CYP3A4\*22**.
⚠️ **Caveat worth propagating:** the paper did **not** ChIP PXR, CAR, HNF4A or PPARA — it did p300
ChIP-qPCR only. Anyone citing it for "HNF4A binds R2" is wrong. It also gives **no genomic coordinates
and no genome build** for R1–R4.

### 1.6 Which factor acts under which stimulus

| stimulus | factor | element(s) | source |
|---|---|---|---|
| **rifampicin** | **PXR:RXRα**, cooperatively at both modules; recruits HNF4α, SRC-1, p300, NCOA6 | XREM + proximal ER6 | [Goodwin 1999](https://doi.org/10.1124/mol.56.6.1329); [Li & Chiang 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1524881/) |
| **phenobarbital** | **CAR:RXR**, overlapping PXR; CAR binds the proximal ER6 only weakly | XREM | [Robertson 2003](https://doi.org/10.1124/mol.64.1.42); [Sueyoshi & Negishi 2001](https://doi.org/10.1146/annurev.pharmtox.41.1.123) |
| **dexamethasone** | **GR — indirectly**, by inducing PXR and RXRα mRNA (max at 100 nM, 6–12 h, RU486-sensitive) | none mapped | [Pascussi 2000](https://doi.org/10.1124/mol.58.2.361). **A direct GRE in CYP3A4: NO EVIDENCE FOUND** |
| **1,25(OH)₂D₃** | **VDR:RXR**, with PXR explicitly excluded | **proximal ER6** | [Thummel 2001](https://doi.org/10.1124/mol.60.6.1399) |
| **lithocholic acid** | **PXR and VDR** (VDR is the higher-affinity LCA sensor) | ER6 / DR3 / DR4 | [Staudinger 2001](https://doi.org/10.1073/pnas.051551698); [Makishima 2002](https://doi.org/10.1126/science.1070477) |
| **chenodeoxycholic acid** | **FXR** | two sites in a 345-bp 5′ element | [Gnerre 2004](https://doi.org/10.1097/00008571-200410000-00001) |
| **St John's wort / hyperforin** | **PXR** — hyperforin is a direct ligand, **Kᵢ = 27 nM** | ER6 / XREM | [Moore 2000](https://doi.org/10.1073/pnas.130155097) |
| **carbamazepine** | **PXR**, by CBZ *and* its 10,11-epoxide | ER6 / XREM | [Oscarson 2006](https://doi.org/10.1016/j.clpt.2006.08.013) |
| **nifedipine and other dihydropyridines** | **PXR**, requiring dexamethasone-induced PXR; **both** PXREs needed; explicitly **not CAR** | proximal ER6 + XREM | [Drocourt 2001](https://doi.org/10.1124/dmd.29.10.1325) |
| **paclitaxel** | **PXR/SXR** — ⚠️ the canonical Synold 2001 citation could not be verified through the endpoints available. **Verify before citing.** | | |
| **IL-6 / inflammation** | **C/EBPβ-LIP** and **NF-κB p65**, see §3 | | |
| **TCDD / 3-methylcholanthrene** | **AhR — represses** CYP3A4 and PXR | no XRE mapped | [Rasmussen 2017](https://doi.org/10.1016/j.toxlet.2017.05.029) |

**Two corrections to the standard story, both worth carrying.** (i) **AhR is a repressor of CYP3A4,
not an inducer** — AhR knockdown *increases* both basal and rifampicin-induced CYP3A4. (ii) **HNF4α's
mode of action is genuinely CONTESTED**: Tirona 2003 shows HNF4α is *required* for PXR/CAR induction
([*Nat Med* 9:220](https://doi.org/10.1038/nm815)), but Li & Chiang 2006 found that mutating the
putative HNF4α site in XREM **did not affect** basal promoter activity and instead showed
rifampicin-activated PXR *recruiting* HNF4α to both modules; and Liu 2008 found HNF4α *increases*
distal-module activity while *decreasing* CLEM4 activity. Required, yes; fixed-site, unresolved.
**NO EVIDENCE FOUND** for SREBP or ESR1 acting on CYP3A4, or for a mapped LXRE.

### 1.7 ENCODE ChIP-seq at the locus — independent, orthogonal confirmation

277 TF clusters over 102 factors in chr7:99,755,000–99,800,000
([UCSC API, `encRegTfbsClustered`](https://api.genome.ucsc.edu/getData/track?genome=hg38;track=encRegTfbsClustered;chrom=chr7;start=99755000;end=99800000)).
⚠️ **The track pools all ENCODE cell types**; ENCODE's own region search over the regulatory block
returns 895 datasets of which 191 are HepG2 and 46 are liver, i.e. ~26% hepatic. The erythroid/lymphoid
factors in the list (GATA1, TAL1, SPI1, PAX5, EBF1, BATF, BCL11A), almost all clustered inside
intron 6/7, are K562/GM12878 and **must not** be read as liver regulation. With that caveat:

- **Proximal promoter:** POLR2A (score 1000), TAF1, SP1, **RXRA (472)** and **HNF4A (123)** peaks span
  the proximal ER6.
- **XREM is the densest, highest-scoring cluster in the whole region**: RXRA, HNF4A, SP1 and JUND all
  at score **1000**, plus ATF3 995, FOXA2 711, FOXA1 651, HNF4G 577, **NR2F2/COUP-TFII 483**.
- **CLEM4 reproduces Matsumura 2004 almost exactly**: **USF1 at 1000**, USF2, **HNF1A 590**, **HNF4A
  652**, **JUND and ATF3 both 1000** (= AP-1) — every factor the 2004 gel shifts named — plus FOXA1,
  FOXA2, SP1 and RXRA at 1000.

That a 1999–2004 map built from DNase footprinting and gel shifts is reproduced peak-for-peak by
modern ChIP-seq is the strongest validation the regulatory section has. Two unexploited leads:
**NR3C1 has two peaks at the locus** (~−9.3 kb and inside CLEM4) despite no GRE ever being mapped,
and **BHLHE40**, a clock repressor, has a maximal-score peak inside CLEM4.

**What this licenses for the structure track.** Almost nothing directly, and that is worth saying
plainly: no amount of promoter mapping changes a ligand pose. Two indirect items do matter.
(i) **CYP3A4\*22 produces a transcript with no heme-binding signature**, so the \*22 protein is not a
variant enzyme to model — it is *less* enzyme. Any attempt to explain inter-individual variation
structurally should not look for a \*22 conformer. (ii) **Liver expresses CYP3A4 and CYP3A5 together
at 335 and 156 TPM**, so CYP3A5 is not an exotic paralog but a co-resident with measurably different
lid plasticity — which makes the six CYP3A5 structures (and the clotrimazole and azamulin pairs,
8SPD/8SG5 and 6OOA/7SV2) the single best leak-free negative control set available for a lid-sensitive
scorer.

---

## 2. What acts on the transcript

Much of what is usually asserted about this layer does not survive checking, so this section is
organised around what has a **3′UTR reporter with a seed mutant** behind it and what does not.

### 2.1 microRNAs — direct 3′UTR targeting IS established, and three commonly cited names are errors

**ESTABLISHED, and it corrects a claim that appears in reviews and in this project's own first pass.**
The founding paper is titled "MicroRNAs regulate CYP3A4 expression via **direct and indirect**
targeting" ([Pan, Gao & Yu 2009, *DMD* 37:2112](https://doi.org/10.1124/dmd.109.027680),
PMC2769037): **miR-27b** represses a CYP3A4 3′UTR–luciferase reporter *and* a VDR 3′UTR reporter,
with CYP3A4 protein down >30%. Both arms operate. "CYP3A4 miRNA regulation is only indirect, through
PXR/VDR/RXRα" is **wrong**.

**Direct — a validated CYP3A4 3′UTR site, ranked by strength of evidence:**

| miRNA | evidence | status |
|---|---|---|
| **miR-27b-3p** | reporter + seed mutant + western + activity (2009); mimic/inhibitor in HepG2 **and** primary human hepatocytes, and **20.0% of CYP3A activity variance across 55 livers** ([Liu 2016, *Sci Rep* 6:26544](https://doi.org/10.1038/srep26544), PMC4876377); bioengineered miR-27b-3p lowers CYP3A4 protein and midazolam 1′-hydroxylase ([*Acta Pharm Sin B* 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6543075/)); circulating miR-27b tracks 4β-OHC/cholesterol (*P* = 0.04, n = 28) with **no association to mRNA, i.e. translational action** ([Ekström 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4777245/)) | **ESTABLISHED, three labs** |
| **miR-206** | reporter with MRE mutant; **5.8%** of CYP3A activity variance (Liu 2016) | ESTABLISHED, one lab |
| **miR-627** | GFP-3′UTR reporter + mutant + western, two independent papers ([Wei 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC3942699/); [Sun 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5011008/)) | ESTABLISHED |
| **miR-200a-3p, miR-150-5p** | mutants + mimics + inhibitors + endogenous CYP3A4 ([Huang 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6546834/)) | supported |
| miR-27a; miR-628-3p, miR-641; miR-122-5p | single-paper reporter evidence | weaker |
| miR-577, miR-1, miR-532-3p | **all three from ONE paper**, HEK293T luciferase + mutant only, **no endogenous CYP3A4 readout anywhere** | weak |
| **svRNAb** (a vault-RNA-derived small RNA, not a miRNA) | targets the CYP3A4 3′UTR; [Persson 2009 *Nat Cell Biol*](https://pubmed.ncbi.nlm.nih.gov/19749744/), replicated in HepG2 and 19 human livers | ESTABLISHED, and usually forgotten |

**Indirect — the upstream regulator, not CYP3A4:** **miR-148a-3p → PXR 3′UTR**
([Takagi 2008, *JBC* 283:9674](https://pubmed.ncbi.nlm.nih.gov/18268015/)) — ⚠️ CONTESTED, a later
study found all the correlations null across 24 livers; **miR-30c-1-3p, miR-18a-5p, miR-140-3p → PXR**;
**miR-34a → RXRA** (a coding-region MRE); **miR-142** (authors state indirect); **miR-107 and
miR-1260**, whose associations **vanish once ESR1 and other TFs enter the model**.

⚠️ **Three names in circulation are errors and should not be repeated.**
**miR-1260b was the NEGATIVE CONTROL** in Liu 2016 ("not expected to bind CYP3A4 or CYP3A5 mRNA").
**miR-21 was tested and failed** — "the luciferase activities of neither the… wild-type nor the mutant
plasmid were affected by miR-21 and miR-130a". **miR-298 in Pan 2009 is mouse mmu-miR-298.**
Separately, **rs4646437 is an intron variant** (chr7:99,767,460, dbSNP `intron_variant`), not a 3′UTR
variant, and **no miRNA mechanism has ever been proposed for it** in any of its 139 papers.

**The quantitative punchline, and it is a strange one.** In the one study that partitions liver
variance three ways (n = 55 livers, Liu 2016): **miR-27b 20.0%, CYP3A4 mRNA 9.5%, miR-206 5.8%.**
**One miRNA explains more than twice as much CYP3A activity variance as CYP3A4 mRNA itself, and
about three times what CYP3A4\*22 explains.**

**[INFERENCE], connecting to §1.4:** the short-3′UTR isoform loses ~709 nt of 3′UTR, is preferentially
expressed in differentiated adult hepatocytes, is preferentially induced by rifampicin, is more
stable (3.7 h vs 2.2 h) and gives 2.4-fold more protein. If the miR-27b/miR-206 sites lie in the
distal 709 nt, **APA would be a miRNA-escape mechanism** and would explain all of that at once. The
site positions within the 3′UTR were not recoverable from open sources, so this is a hypothesis with
an obvious experiment attached, not a result.

### 2.2 RNA-binding proteins — NO EVIDENCE FOUND

No RBP has demonstrated binding to the CYP3A4 transcript in the literature or in the CLIP-derived
resources queried. Reported as **NO EVIDENCE FOUND (search-limited)**: ENCODE's HepG2 eCLIP panel
covers ~150 RBPs, so the data to answer this may already exist and simply have not been analysed for
this gene.

### 2.3 uORFs and microproteins — a uAUG exists, and it has zero ribosome support

**The 5′UTR is 103 nt and contains exactly one uAUG, at nt 27, opening an 11-codon uORF
(MHIAQQRATQS) that stops at nt 60–62 — weak Kozak, terminating 41 nt before an optimal-Kozak main
AUG.** That was computed independently twice from NM_017460.6.

**But it has no experimental support: NO EVIDENCE FOUND.** It is absent from the GENCODE Ribo-seq ORF
phase I (7,264 ORFs) and phase II (28,359) catalogues and from sORFs.org — all three checked with
working positive controls (BUD31, 350 kb away, has 7 catalogued uORFs) — and there is no literature
across six query formulations. **This is the right shape of answer for this section: a sequence
feature that is real, an annotation that is absent, and no claim made in between.**

### 2.4 lncRNAs — exactly one, and it is contested

**CYP3A4-AS1 / AC069294.1** (ENSG00000273407, HGNC:58107, approved 2024) is a **350 nt single-exon
antisense transcript inside the CYP3A4 gene body**, chr7:99,766,543–99,766,892 (+). One functional
paper exists: [Collins & Wang 2022, *Pharmacogenet Genomics* 32:16](https://pmc.ncbi.nlm.nih.gov/articles/PMC8578198/)
— knockdown in Huh7 raises CYP3A4 mRNA ~3×, overexpression lowers it 89%, and **CYP3A4\*1G
(rs2242480) associates with 1.26× more lncRNA and −31% CYP3A4**.

⚠️ **Three reasons to hold it loosely.** Its liver abundance is **0.53 TPM against CYP3A4's 335 TPM —
630-fold less than its target**; there are **zero replications**; and an independent WGCNA of the same
GTEx v8 liver data, hunting exactly this class of transcript, named LINC02499, HNF4A-AS1, DBH-AS1 and
AC027682.6 and **did not find it**
([Huang 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9775998/)). Everything else in the cluster
(lnc-CYP3A7-\*, lnc-CYP3A43-\*, transcripts over CYP3AP2) is **annotation-only with no functional
literature**, and **CYP3AP1/CYP3AP2 ceRNA or sponge activity: NO EVIDENCE FOUND.**

### 2.5 mRNA stability and translational control

- **APA changes stability, measurably**: short-3′UTR **3.7 h** vs long-3′UTR **2.2 h**, with 2.4-fold
  more protein ([Li 2012](https://doi.org/10.1124/mol.111.074393)).
- **miR-27b acts translationally, not on transcript level** — hepatic miR-27b tracks CYP3A activity
  with **no association to mRNA** (Ekström 2015).
- **IL-6 repression is transcriptional, not a stability effect** (§3.5). **NO EVIDENCE FOUND** for
  IL-6- or TNF-driven destabilisation of CYP3A4 mRNA.
- **m6A — CONTESTED and weak, not established.** The METTL3 → CYP3A4 story rests on two papers **from
  the same laboratory**, so it is not independent convergence, and the Mettl3/Hnf4a axis in them is
  explicitly developmental. Against it,
  [Nakano 2020, *Biochem Pharmacol* 171:113697](https://pubmed.ncbi.nlm.nih.gov/31706844/) found that
  3-deazaadenosine raised CYP1A2, CYP2B6 and CYP2C8 in HepaRG while **CYP3A4 was conspicuously not a
  responder**, and m6A *decreased* CYP2C8 — the opposite sign. **Mark this contested; do not repeat it
  as a mechanism.**
- **A-to-I editing — a measured null worth recording.** REDIportal gives **141 sites, 86 of them
  exonic, all inside a single AluSx in the 3′UTR, median editing 0.42%, zero CDS sites and zero
  recoding**, and no inverted-repeat duplex can form. So CYP3A4 mRNA *is* edited, trivially, in a way
  that cannot change the protein.
- **IRES and codon-usage effects: NO EVIDENCE FOUND** for CYP3A4 specifically *(search-limited)*.
- **The largest post-transcriptional control is at the protein level.** CYP3A4 protein and activity
  oscillate circadianly **with no mRNA rhythm**, even from a constitutive CMV promoter, because gp78
  oscillates antiphase (`CYP3A4_EVOLUTION.md` §5.4).

### 2.6 The enhancers act by pre-formed looping, and rifampicin does not restructure it

**ESTABLISHED, and it replaces an "untested" claim.** The 4C study was run **± rifampicin**: "10
million primary culture human hepatocytes… were treated with DMSO or rifampicin (25 μM) for 18 hrs",
and "treatment with the CYP3A4 inducer Rif yielded **similar interaction patterns**, with decreased or
increased signals occurring at some regions"
([Collins & Wang 2020](https://doi.org/10.1097/FPC.0000000000000402), PMC8457025). **The contacts are
largely pre-formed; PXR activation works inside an existing loop rather than building one.**

From the same paper, a clean division of labour: "**although deletion of XREM had no effect on basal
CYP3A4 expression, it reduced induction of CYP3A4 by Rif**" — and deleting XREM cut rifampicin
induction from **4.6 ± 1.5-fold to 1.3 ± 0.4-fold**. **XREM is induction-specific; CLEM4/R2 carries
the constitutive load.**

**A second 4C paper adds the shared enhancer.** [Collins, Nworu, Mohammad *et al.* 2022, *Clin Transl
Sci* 15:2720](https://pubmed.ncbi.nlm.nih.gov/36045613/) show a **distal regulatory region ~90 kb from
the CYP3A4 promoter that loops to CYP3A4, CYP3A5 and CYP3A7, with CRISPR deletion dropping all
three** — one shared enhancer for three genes — and the African-specific **rs115025140 giving 1.8×
mRNA across 246 livers**. ⚠️ **Reproducibility flag:** the paper gives **no chr7 coordinates and names
no genome assembly** for that region, so any downstream coordinate for the DRR is an inference.

**eRNA at XREM — CONTESTED / WEAK, not absent.** FANTOM5 CAGE gives **liver adult pool1 0.143 TPM** at
0.94× median depth and **hepatocyte donor3 0.568 TPM**, with **HepG2 at 0.000**; the top-expressing
libraries are all non-hepatic. That is trace-level. ⚠️ And **FANTOM5 is entirely unstimulated**, so an
*inducible* eRNA at XREM is **untested by any dataset**.

### 2.7 Three things nobody has measured, and two search traps

**Not measured by anyone — keep this list explicit, because each is a one-experiment gap:**

1. **No transcription-shutoff half-life for CYP3A4 mRNA in primary human hepatocytes.** The existing
   estimates disagree by **4–8×** and have never been reconciled.
2. **No test of whether the \*22 intron-6-retention transcript is an NMD substrate.** Ensembl biotypes
   it `retained_intron`, not NMD, and it carries its own internal polyA signal — but nobody has run
   the UPF1 knockdown.
3. **No nuclear run-on and no Pol II ChIP at CYP3A4.** So "**rifampicin increases the transcription
   rate**" — repeated everywhere — **has never been measured as a rate.** Everything is steady-state
   mRNA.

**Two search traps that will bite anyone re-running this work.** (i) **Europe PMC phrase queries
URL-encoded with `%20` silently return `hitCount: 0`** rather than an error — use `+`. I hit this
independently five times in this session and it is indistinguishable from a real negative.
(ii) **The entire CYP3A4\*22 splicing literature is indexed under "intron 6", not "alternative
splicing", and the 4C papers are only findable via "chromatin conformation capture".** A keyword
search on the obvious term returns nothing for both.

**What this licenses for the structure track.** One thing, and it is worth stating because it closes
off a whole class of worry: **no post-transcriptional mechanism produces a structurally different
CYP3A4 protein.** \*22 makes a transcript with no heme-binding signature (i.e. less enzyme, not
another enzyme); the miRNAs and APA change abundance and translation rate; the A-to-I editing is
confined to a 3′UTR Alu at 0.42% and cannot recode anything. **Every route into inter-individual
variation is a route into *how much* enzyme there is, not *which* enzyme.** So the single protein we
model is the right one, and variability in CYP3A4 activity is not a reason to expect conformational
heterogeneity in the cryoEM references.

---

## 3. Epigenetics

### 3.1 The CYP3A4 promoter is not a CpG island, and the whole cluster sits in a desert

**ESTABLISHED, by annotation and by sequence, independently reproduced twice for this document.**

- **UCSC `cpgIslandExt` (hg38) returns zero CpG islands across the entire CYP3A cluster.** The gap runs
  **chr7:99,617,262 → 99,918,979 — a 301,717 bp island desert** containing CYP3A5, CYP3A7, CYP3A4,
  CYP3A4-AS1 and CYP3A43 in their entirety. Nearest island to the CYP3A4 TSS: **134,795 bp**. Sanity
  check: widening to a 1.5 Mb window returns 39 islands, so the track and the API work.
  *(I re-ran the narrower query independently: 0 islands in chr7:99,750,000–99,800,000.)*
- **Direct sequence computation** over chr7:99,770,000–99,800,000 (30 kb spanning the TSS, XREM and
  CLEM4): **GC 41.4%, CpG obs/exp 0.14, 181 CpGs (6.0/kb)**. **0 of 597** 200-bp windows pass
  Gardiner-Garden & Frommer and **0 of 590** 500-bp windows pass Takai & Jones. The first 500 bp
  upstream of the TSS contains **2 CpGs**; the next 500 bp contains **none**.
- **NO EVIDENCE FOUND for a quotable published sentence** saying so. The claim rests on annotation
  plus computation, which is stronger than a citation.

**Two consequences.** Every "CYP3A4 promoter methylation" result is about a handful of **isolated,
non-island CpGs** — which is why positions are always quoted individually. And **only 10 Illumina EPIC
probes exist in the whole CYP3A4 gene region**, because arrays cluster on islands. Methylation arrays
are close to blind here.

### 3.2 The complete published per-CpG inventory

| position (vs TSS) | finding | source |
|---|---|---|
| **−383** | hypermethylated in neonates vs adolescents, **P = 0.00001**, 48 pediatric + 34 prenatal livers | [Vyhlidal 2016, *DMD* 44:1020](https://pubmed.ncbi.nlm.nih.gov/26772622/) |
| **−1547, −10,762** | the two Kacevska CpGs reported to associate with CYP3A4 expression; −10,762 is **inside CLEM4** | Kacevska 2012, as quoted in [Habano 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4587720/) |
| **−1521, −1569, −10,813, −10,851, −10,895** | hypermethylated in gastric cancer vs adjacent healthy (*p* = 0.003 to ~0) — **independently replicates the same two hotspots** | [Golestanian 2022](https://pubmed.ncbi.nlm.nih.gov/35331105/) |
| **−5998, −5731, −5725** | C/EBP-site CpGs hypermethylated by diazinon *while transcription rose 27-fold* — **paradoxical sign, CONTESTED** | Golestanian 2022 |
| −36, −82/−86, −258, −296, −367/−372/−374 | associated with risperidone response, n = 288. **CONTESTED** — peripheral blood, not liver, no replication | [Shi 2017](https://pubmed.ncbi.nlm.nih.gov/28696411/) |
| cg19046783 | Spearman *r* = 0.52, explains **18%** of dose-normalised tacrolimus AUC. **CONTESTED** — n = 23, 1 of 10 probes, no multiple-testing correction | [Koudijs 2025](https://pubmed.ncbi.nlm.nih.gov/40632895/) |

**The structural coincidence is the strongest argument the sparse signal is real:** the hotspots land
on the known regulatory modules — the −1.5 kb proximal region, the −5.95 kb C/EBPβ element, and
CLEM4 at −10.8 kb.

**Kacevska 2012** ([*Biochimie* 94:2338](https://pubmed.ncbi.nlm.nih.gov/22906825/)) is the canonical
reference, and what it says is an *association*: bisulfite sequencing of a ~12 kb region in **72 adult
and 7 fetal livers** found "highly variable CpG methylation sites… which correspond to important
CYP3A4 transcription factor binding sites including the **proximal promoter, XREM and CLEM4** as well
as in separate **C/EBP and HNF4α binding regions**", with **fetal hypermethylation** relative to
adult. ⚠️ Its full text is unobtainable — paywalled, no PMC deposit, the only green OA copy is a
deleted record. **A percent-methylation-by-CpG-by-age table for CYP3A4 does not exist in accessible
form, and nobody downstream should write one.**

⚠️ **Two attributions that circulate and are wrong.** There is **no Habano paper on CYP3A4 promoter
methylation in intestine vs liver** — Habano 2011 shows the methylated gene in colon is **PXR**, with
CYP3A4's own 5′-distal methylation "**equal across all six cell lines**", i.e. it cannot explain the
expression differences; and **Habano 2015 is a well-powered NEGATIVE** for CYP3A4 across 20 livers:
"although considerable inter-individual differences were observed in CYP3A4 expression (CV 22.9%), we
**did not detect a significant relationship between DNA methylation and mRNA expression**". And
**Okino/Tokizane are about CYP1A1 and CYP1B1, not CYP3A4** — both of which have CpG-island-like
promoters, so the mechanism does not transfer.

**Does demethylation reactivate CYP3A4?** Globally yes, causally unproven. 5-aza-dC changes CYP3A4/5/7
in HepG2 and raises CYP3A4 mRNA 28–116-fold in colon lines — but **via PXR demethylation**. The
cleanest causal result is that **in vitro methylation of the CYP3A4 enhancer blocks PXR binding and
abolishes rifampicin induction** ([Wang 2021](https://pubmed.ncbi.nlm.nih.gov/33048369/)). Against it,
Habano 2015 found CYP3A4 **not** among the DAC-responsive DME genes. **Verdict: demethylation
reactivates the CYP3A4 *circuit*; there is no clean demonstration that demethylating CYP3A4's own
promoter is sufficient.**

### 3.3 Histone marks and accessibility — the signal is at the enhancers, not the promoter

**ESTABLISHED across four independent liver datasets.** In Roadmap **E066 (adult liver)**, a single
contiguous **7.4 kb H3K27ac domain spans the entire XREM–CLEM4 interval at fold-enrichment 33.5
(q = 184.9)**, against 9.56 at the promoter. H3K4me1 is enhancer-skewed (CLEM4 8.59, XREM 5.22,
promoter 3.16 — the weakest). The **E066 chromHMM 18-state model calls `1_TssA` at the promoter and an
unbroken `9_EnhA1`/`10_EnhA2` active-enhancer block from −4.2 kb to −12.6 kb**, super-enhancer-like,
with XREM and CLEM4 both inside it. **H3K27me3 and H3K9me3: zero peaks in 30 kb** (files verified
non-empty; nearest peaks 198 kb and 68 kb away).

Three ENCODE liver donors reproduce the ratio: distal/promoter H3K27ac fold-enrichment **6.2×, 3.3×
and 2.2×**. **ATAC-seq** (ENCODE's only human liver experiment, ENCSR124NNL) puts **CLEM4 as the most
accessible element in 56 kb — signal 8.26, q = 329.7 — about 4× the promoter's q (83.0)**, with XREM
also called. ENCODE liver DNase calls promoter, XREM and CLEM4 as reproducible DHS in **4 of 6
donors**, while **embryonic liver has only 2 peaks and neither is the promoter or CLEM4** — consistent
with postnatal activation.

⚠️ **HepG2 is epigenetically dead at the CYP3A4 promoter**, on two platforms: **exactly one H3K27ac
peak in 56 kb** (at CLEM4), **no H3K4me3 peaks anywhere**, and no promoter DHS. Roadmap E118 agrees.
**Do not use HepG2 as the CYP3A4 chromatin model** — which is awkward, because much of the reporter
literature is HepG2.

**PXR-ligand-induced changes are real and are the MLL/SET answer, weaker than usually stated.**
Rifampicin in LS174T (⚠️ **colon cells, not hepatocytes**) gives, at 96 h, **H3K4me3 ×1.8 distal /
×3.5 proximal, H3K27me3 −45% / −53%, H3 acetylation up**, with **NCOA6 recruited 4.9×/6.4× and p300
3.7×/3.4×**; siPXR raises **H3K27me3 8.9×/6.6×** and siNCOA6 raises it 5.3×/14.0× while abolishing
induction entirely ([Yan 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5508193/)). NCOA6 is the
**ASCOM** adaptor that delivers **MLL3/MLL4 (KMT2C/KMT2D)** and the **UTX/KDM6A** demethylase — but
⚠️ **no paper ChIPs KMT2C or KMT2D at CYP3A4; the MLL3/MLL4 link is INFERRED from NCOA6 alone**, and
**SET7/9 at CYP3A4: NO EVIDENCE FOUND.**

**New, and it cuts against the Roadmap null:** EZH2 inhibition (GSK126) and **CBX4 knockdown both
raise CYP3A4**, CBX4 binds PXR, and ChIP places **both CBX4 and PXR at the promoter and the enhancer**
([Kamiyoshihara 2026, *DMD*](https://pubmed.ncbi.nlm.nih.gov/42685388/)). **[INFERENCE]:** H3K27me3
sits below peak-calling threshold in bulk adult liver *because* PXR/NCOA6 keeps it low — consistent
with siPXR raising it 8.9-fold — rather than being absent.

**CTCF: an established negative in liver.** Zero CTCF peaks in chr7:99,750,000–99,806,000 across
**four independent ENCODE liver CTCF ChIP-seq files**. HepG2 does have a strong site in CYP3A4 intron
1. Over 260 kb the cluster carries **71 cCREs — 51 distal enhancer-like, and only 3 CTCF-only**, none
within 40 kb of the CYP3A4 TSS. **[INFERENCE]:** enhancer-dense and CTCF-sparse, consistent with one
shared distal enhancer serving three promoters rather than each gene in its own insulated sub-TAD.
**Hi-C or TAD data at chr7q22.1: NO EVIDENCE FOUND — a genuine, exploitable gap.**

### 3.4 The CYP3A7 → CYP3A4 neonatal switch

**The ontogeny is ESTABLISHED and quantitative.** CYP3A7 is **75% of total hepatic CYP protein at
0–12 days**, falls to **2% by 1–2 years**, and has an **age₅₀ of 26 days**; CYP3A4 is **not detected**
in neonates and reaches **49% of total CYP at >1–2 y** and ~16% in adults
([Subash & Prasad 2024, *DMD*](https://pmc.ncbi.nlm.nih.gov/articles/PMC11585312/), 50 pediatric + 8
adult donors). In hepatoblasts the single-cell split is **CYP3A7 69.9% of cells vs CYP3A4 0.09%
(10 of 11,673)**. Critically,
[Lacroix 1997](https://pubmed.ncbi.nlm.nih.gov/9266706/) showed **total CYP3A protein stays nearly
constant** while CYP3A4 reaches 30–40% of adult by one month: **it is a swap, not a net induction.**

**Is it methylation-driven? Genuinely CONTESTED — four datasets, four partly incompatible answers.**

1. **For:** a cytosine **−383 bp from CYP3A4 is hypermethylated in neonates (P = 0.00001)** and
   CYP3A7 proximal-promoter cytosines are hypomethylated in neonates, both tracking the mRNA
   trajectories (Vyhlidal 2016); Kacevska 2012 finds fetal hypermethylation at promoter/XREM/CLEM4.
2. **Against, in adults:** across 70 livers, **DNA methylation contributes nothing to CYP3A4 once TFs
   are in the model — TFs 69.71%, CYP3A4\*22 2.56%, total R² 72.27%**
   ([Collins & Wang 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11596782/)). For CYP3A7,
   methylation does add 4.09%.
3. **Against, in mouse:** "**DNA was not hypermethylated in the Cyp3a locus at any age**"
   ([Li 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2672803/)); the switch tracks histone marks.
4. **Against, mechanistically:** published back-to-back with Vyhlidal in the same issue,
   [Giebel/Hines 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4931893/) finds the CYP3A4 locus
   **bivalent and poised in BOTH fetal and postnatal liver — the gene is never closed** — and states
   outright that the CYP3A7 mark pattern "**is inconsistent with known CYP3A7 developmental
   expression**" and that their human H3K27me3 pattern "**contrasts with** both of the patterns
   previously observed for mouse".

**Synthesis, and it is TF-first [INFERENCE, but well supported].** All four groups agree the locus is
**poised, not closed** — nothing needs opening. What changes in *trans* is the **abundance of PXR and
CAR**: prenatally CAR is 757 ± 480 and PXR 271 ± 190 molecules/ng and both rise postnatally, while
**RXRα and HNF4α are flat across the window** (Vyhlidal 2006), and postnatal CYP3A4 correlates with
PXR at *r*² = 0.610 and CAR at 0.723. **So "HNF4A arrives at birth and turns on CYP3A4" cannot be
right as stated: HNF4A is necessary but is not the timer.** And histone marks at this locus are
demonstrably *downstream* of TF binding — PXR knockdown abolishes all rifampicin-driven mark changes.

**CYP3A7\*1C clinches the cis side.** The allele is literally **a CYP3A4 promoter segment pasted into
CYP3A7**, giving CYP3A7 the **CYP3A4 proximal ER6** — and **PXR and CAR bind CYP3A4-ER6 with higher
affinity and transactivate only ER6-containing constructs**
([Burk 2002, *JBC* 277:24280](https://pubmed.ncbi.nlm.nih.gov/11940601/)). Carriers have **~50% lower
DHEAS** and, by GWAS, **−49.2% urinary oestrone-3-glucuronide (P = 3.1 × 10⁻¹⁸)**. **The switch is a
cis-encoded difference in promoter tuning to trans-factors that rise after birth**, not a locus-level
on/off device — and the whole thing is cis-encoded within 700 kb, because a human artificial
chromosome carrying the CYP3A locus **recapitulates stage-specific human CYP3A expression in mice**
([Kazuki 2013](https://pubmed.ncbi.nlm.nih.gov/23125282/)).

⚠️ **The glucocorticoid-surge story is in vitro only.** Dexamethasone induces CYP3A4 in fetal human
hepatocytes in a **GR-mediated, PXR-independent** way (RU486-sensitive; rifampicin does *not* work,
because **PXR mRNA is not expressed in human fetal liver cells**). But **NO EVIDENCE FOUND** links the
measured neonatal cortisol surge to measured CYP3A4 induction in humans, and awkwardly, **GR is the
best correlate of CYP3A7 — the gene being switched off.**

⚠️ **The shared enhancer points away from promoter competition.** The **DRR** (R1, ~90 kb from the
CYP3A4 promoter) is a real shared enhancer: p300 ChIP fold-enrichment 68–130 across three donors, and
**CRISPR deletion reduces CYP3A4, CYP3A5 and CYP3A7 together**
([Collins 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9652438/)). Under a winner-take-all
competition model, deleting the shared enhancer should *relieve* one gene; it hurts all three.
**NO EVIDENCE FOUND for actual promoter competition — nobody ran the experiment**, and all of this
work is adult liver, none developmental. **The biggest unexplained hole is the mechanism of CYP3A7
shut-off**; the only serious candidate is **NFI**, whose isoforms activate the CYP3A4 promoter and
repress CYP3A7's with developmentally distinct complexes
([Riffel/Vyhlidal 2009](https://pubmed.ncbi.nlm.nih.gov/19706729/)) — **unfollowed-up for 17 years.**

### 3.5 What actually drives the 10–100× inter-individual variability

**"10–100 fold" is a review composite, not one measurement.** The cleanest primary numbers are **31×**
in testosterone 6β-hydroxylase activity across 46 livers (Westlind 1999), **>30×** in 20 intestines
(Paine 1997), **7× basal and 11× induced** *in vivo* across 367 twins (Rahmioglu 2011), and
**p95/p5 = 218× in GTEx v8 liver TPM across 226 donors** (Q3/Q1 only 8.2×). Structurally,
**in-vitro expression varies far more than in-vivo clearance, and mRNA > protein > activity**.
CYP3A4 phenotype is **unimodal** — there is no poor-metaboliser mode.

**The heritability figure everyone quotes is an artefact of the method. CONTESTED, and this is the
most important correction in this section.**

| study | design | estimate |
|---|---|---|
| Ozdemir 2000 | ⚠️ **not a twin study** — a meta-analysis of repeat-dosing studies (n = 161 across 16 studies) using Kalow's rGC | **rGC 0.96** for midazolam clearance |
| Rahmioglu 2011 | TwinsUK, n = 367, quinine probe after 14 d St John's wort | **h² = 66% for INDUCED activity**; ⚠️ uninduced was not recorded |
| **Matthaei 2020** | **the real classical twin study** — 43 MZ + 14 DZ pairs, midazolam | **additive genetic 15% of midazolam AUC variance; COMMON ENVIRONMENT 48%** |

rGC counts as "genetic" **any** source of variance stable *within* a person — diet, chronic
comedication, BMI, smoking, chronic low-grade inflammation, stable epigenetic state. Kalow's own
method paper concedes this; Ozdemir's finding of *higher apparent genetic control at night than by
day* is direct evidence the statistic absorbs non-genetic structure; and Vesell's group showed in 1984
that intratwin correlation was **nearly twice as high for twins living together as apart**. **And
common environment is exactly the term rGC cannot see.** Converging nulls: a 310-twin GWAS found **no
genome-wide hit**; a 466-liver eQTL study found **only one trans-SNP** touching testosterone
hydroxylation; and in 134 cancer patients with a 14-fold activity range **no CYP3A4/CYP3A5 variant was
a significant predictor (all P > 0.29)**.

**Defensible line: the 66–96% heritability figures are an upper bound from a method that cannot
separate genes from stable environment, not an established fact with missing heritability waiting to
be found.**

**Known variants explain ~7–25% in total.** **CYP3A4\*22** (NFE 4.9%, AFR 0.9%, EAS 0%) was discovered
*by* allelic-expression-imbalance screening and explains **7% of mRNA and 12% of activity variance**;
**CYP3A4\*1B is a non-functional LD marker for CYP3A5\*1** (it is 64% in AFR where CYP3A5\*1 is 69.5%
and ~0 in EAS); **CYP3A4\*1G** is CONTESTED and its best mechanism is the CYP3A4-AS1 lncRNA of §2.4,
also in high LD with CYP3A5\*1; and **CYP3A5\*3 is a first-order confound** because midazolam,
nifedipine, quinine and 4β-OHC all report CYP3A4 + CYP3A5 combined. The largest single trans effects
found — **POR\*28** and **PPARA rs4253728 (8–9% of activity)** — **are not in CYP3A4 at all.**

⚠️ **Do not report any GTEx eQTL null for CYP3A4.** An apparent "zero liver eQTLs" result was traced
to an **API coverage artefact**: the same endpoint returns zero Liver rows for CYP3A5 and never
returns rs776746 in any tissue, and the `independentSqtl` endpoint returns zero rows for everything.
The GTEx *expression* data are unaffected.

**Sex: ESTABLISHED in tissue (~2×), CONTESTED in vivo.** CYP3A4 protein and mRNA are **2× higher in
women** across 94 surgical livers (P < 0.0001) with 50% higher verapamil N-dealkylation, and the
difference is *larger* in the drug-naive subgroup ([Wolbold 2003](https://pubmed.ncbi.nlm.nih.gov/14512885/))
— but a 52-subject stable-isotope midazolam study found "**no significant difference**" in systemic or
oral clearance. **[INFERENCE]:** the effect is real on hepatic protein but washes out in systemic
clearance of a high-extraction, partly flow-limited probe. That mismatch is itself a warning against
equating expression variance with clearance variance.

**Imprinting: NO EVIDENCE FOUND.** No CYP3A gene appears at any confidence level in the Geneimprint
human catalogue (all 282 rows parsed), and no study has tested parent-of-origin expression here.
⚠️ The one search hit that looks relevant — rat androgen "imprinting" of Cyp3a2 — is developmental
hormone programming, a keyword false positive. **Allele-specific expression, by contrast, is
ESTABLISHED and cis-genetic**: strong allelic imbalance goes with *low* total CYP3A4 mRNA, the
opposite of an imprinting model, and AEI screening is how \*22 was found.

**What this licenses for the structure track.** Three things, all deflationary and all useful.
(i) **No epigenetic mechanism produces a different protein** — see §2's closing line; the whole
epigenetic layer modulates abundance. (ii) **The regulatory action is at enhancers 8–11 kb away, not
at the promoter**, which is why HepG2 — the workhorse of the reporter literature — has no CYP3A4
promoter chromatin at all; anyone reasoning from HepG2 data about CYP3A4 should know that.
(iii) **The heritability correction matters for how we talk about the challenge.** If ~48% of
midazolam AUC variance is common environment rather than genetics, then inter-individual differences
in CYP3A4 handling are mostly not encoded in the protein — which is consistent with there being
exactly one structure to predict, and with `FINDING_022`'s framing that CYP3A4's difficulty is
geometric, not genetic.

---

## 4. Pathways, primary and secondary

`CYP3A4_EVOLUTION.md` §3 inventories the endogenous reactions from UniProt and argues from them to
pocket design. This section adds the **database membership with stable IDs**, the **Rhea reaction
IDs**, the downstream phase-II/transporter coupling, and one headline negative about what the
pathway databases actually contain.

### 4.1 Reactome — complete, and surprisingly sparse

Queried live via `reactome.org/ContentService/data/mapping/UniProt/P08684/pathways?species=9606`.
(The `/pathways/low/entity/P08684/allForms/9606` form **404s** — it wants a Reactome stId, not a
UniProt accession. CYP3A4's ReferenceGeneProduct is **R-HSA-52639**.)

| stId | pathway | hierarchy |
|---|---|---|
| **R-HSA-211981** | **Xenobiotics** | Metabolism → Biological oxidations (R-HSA-211859) → Phase I (R-HSA-211945) → Cytochrome P450 by substrate type (R-HSA-211897) → Xenobiotics |
| **R-HSA-5423646** | **Aflatoxin activation and detoxification** | a *direct* child of Biological oxidations |
| **R-HSA-211945** | Phase I – Functionalization of compounds | direct, via the two "binds inhibitor" reactions |
| **R-HSA-9027307** | Biosynthesis of maresin-like SPMs | Metabolism of lipids → SPMs → DHA-derived → maresins |
| **R-HSA-9749641 / 9754706 / 9757110** | **Aspirin / Atorvastatin / Prednisone ADME** | Drug ADME (R-HSA-9748784), a top-level pathway with 8 children |

Sixteen reactions in total, including the aflatoxin trio (R-HSA-156526, 5423664, 5423672), the
atorvastatin 2-OH/4-OH quartet (R-HSA-9756138/9756162/9756169/9756180), salicylate → 2,3- and
2,5-DHBA, prednisone/prednisolone oxidation, loperamide N-demethylation and two inhibitor-binding
reactions.

⚠️ **A parsing trap worth propagating.** R-HSA-211882 is *named* "CYP3A7 can 6beta-hydroxylate
testosterone", but its catalyst is a **CYP3A set containing CYP3A43, CYP3A4, CYP3A5 and CYP3A7**.
Filtering Reactome by display name loses CYP3A4's only steroid reaction in the database.

**The headline negative — MEASURED by enumerating participants.** Under "Cytochrome P450 arranged by
substrate type" there are six substrate-class children. **CYP3A4 is in exactly one.**

| sibling pathway | CYP3A4 present? |
|---|---|
| R-HSA-211976 Endogenous sterols | **NO** (11A1, 11B1/2, 19A1, 21A2, 27A1, 39A1, 46A1, 51A1, 7A1, 7B1, 8B1, 1B1, 4V2) |
| R-HSA-211935 Fatty acids | **NO** |
| R-HSA-211979 Eicosanoids | **NO** |
| R-HSA-211916 Vitamins | **NO** (24A1, 26A1/B1/C1, 27B1, 2R1) |
| R-HSA-211958 Miscellaneous | **NO** — but **CYP3A43 is** |
| R-HSA-211981 Xenobiotics | **YES** |

**Reactome models CYP3A4 as a xenobiotic enzyme and nothing else.** Its 4β-cholesterol hydroxylase,
testosterone/cortisol 6β-hydroxylase, 14,15-epoxygenase and vitamin D activities — all of which
UniProt documents with Rhea IDs and experimental evidence — **are absent**. Build a CYP3A4 map from
Reactome alone and you get the drug half and none of the endobiotic half. (Likewise, **Reactome does
not model PXR → CYP3A4 induction at all**: NR1I2 maps only to the generic
R-HSA-383280 Nuclear Receptor transcription pathway, and CYP3A4 is not a participant in it.)

### 4.2 KEGG disagrees with Reactome, and is broader

Verbatim from `https://rest.kegg.jp/get/hsa:1576`:

```
PATHWAY  hsa00140 Steroid hormone biosynthesis     hsa00591 Linoleic acid metabolism
         hsa00830 Retinol metabolism               hsa00980 Metabolism of xenobiotics by CYP450
         hsa00982 Drug metabolism - cytochrome P450 hsa00983 Drug metabolism - other enzymes
         hsa01100 Metabolic pathways               hsa04976 Bile secretion
         hsa05204 Chemical carcinogenesis - DNA adducts
         hsa05207 Chemical carcinogenesis - receptor activation
```

⚠️ **hsa00590 arachidonic acid metabolism does NOT list CYP3A4**, despite four arachidonate Rhea
reactions in UniProt. Also verbatim: `ORTHOLOGY K17689 [EC:1.14.13.32 1.14.14.55 1.14.14.56
1.14.14.57 1.14.14.73 1.14.14.-]` — note **EC 1.14.14.57 taurochenodeoxycholate 6α-hydroxylase** and
**EC 1.14.14.56 1,8-cineole 2-exo-monooxygenase**; `DISEASE H01143 Vitamin D-dependent rickets`;
`DRUG_TARGET Cobicistat`.

**The two databases genuinely disagree about what CYP3A4 does.** KEGG places it in steroid, bile,
linoleate and retinol metabolism; Reactome places it only in xenobiotics. Use both or neither.

### 4.3 Endogenous reactions with Rhea IDs

From the 43 CATALYTIC ACTIVITY blocks in [P08684](https://rest.uniprot.org/uniprotkb/P08684.txt),
all ESTABLISHED with ECO:0000269 experimental evidence:

- **Steroids.** testosterone → 6β **RHEA:46296**, 1β **RHEA:53260**, 2β **RHEA:53264**, 11β
  **RHEA:53444**; DHT → 18-OH **RHEA:53212**, 19-OH **RHEA:53200**; androstenedione → 6β
  **RHEA:47256**; **progesterone → 6β RHEA:47252 and 16α RHEA:47260**; **cortisol → 6β-OH
  RHEA:55828**; cortisone → 6β **RHEA:55832**.
- **Estrogens, all four positions.** E2 → 2-OH **RHEA:47212**, 4-OH **RHEA:47280**, 16α
  **RHEA:47332**, 16β **RHEA:47336**; estrone → 2-OH **RHEA:47208**, 4-OH **RHEA:47292**, 16α
  **RHEA:47204**.
- **Cholesterol, all five positions, one paper (PMID 21576599).** 4β **RHEA:46128**, 22R
  **RHEA:46140**, 24R **RHEA:46144**, 25 **RHEA:50256**, 26 **RHEA:50264**. CYP3A4 carries the
  UniProt AltName **"Cholesterol 25-hydroxylase"**.
- **Fatty acids.** arachidonate bisallylic 7-/10-/13-HETE **RHEA:52288 / 52296 / 52292**; linoleate
  11-OH **RHEA:52284**. **Epoxygenation is ω-6-selective:** arachidonate → (14R,15S)-EET
  **RHEA:49860** and (14S,15R)-EET **RHEA:49856**; EPA → (17R,18S)-epoxy-ETE **RHEA:39779**; DHA →
  (19R,20S)- **RHEA:52120** and (19S,20R)- **RHEA:52124**.
  ⚠️ **"CYP3A4 makes EETs" is only correct for 14,15-EET.** UniProt documents no 5,6-, 8,9- or
  11,12-EET from arachidonate, and GO lists only **GO:0008404 arachidonate 14,15-epoxygenase**.
  UniProt's own wording: "epoxidation of double bonds of PUFA with a **preference for the last double
  bond**."
- **Anandamide, all three.** 8,9- **RHEA:53140**, 11,12- **RHEA:53144**, 14,15- **RHEA:53148**.
- **Retinoids.** retinol → retinal **RHEA:42092**; **all-trans-retinoate → 4-OH-retinoate
  RHEA:51984** (GO:0008401, IDA).
- **Xenobiotics with Rhea.** quinine 3-monooxygenase **RHEA:20149**; albendazole **RHEA:55924**;
  fenbendazole **RHEA:55928**; **1,8-cineole RHEA:32895**; 1,4-cineole **RHEA:49160**.

**⚠️ Vitamin D — the usual premise needs correcting.** There are **no Rhea IDs** for any CYP3A4
vitamin D reaction; the evidence is GO/IDA and literature only (GO:0062181 1α,25-(OH)₂D₃
**23-hydroxylase**, GO:0070576 **24-hydroxylase**, GO:0030343 **25-hydroxylase**). And the position
depends on the substrate:

- On **calcitriol**, CYP3A4 does **23- and 24-hydroxylation**, giving 1,23R,25-(OH)₃D₃,
  1,24S,25-(OH)₃D₃ and 1,23S,25-(OH)₃D₃, with **opposite stereochemical preference to CYP24A1**; it
  is the dominant source of that activity in human intestine and liver because CYP24A1 is essentially
  absent there ([Xu *et al.* 2006, *Mol Pharmacol* 70:1774](https://doi.org/10.1124/mol.105.017392),
  PMID 16207822).
- On **25(OH)D₃**, the dominant CYP3A4 route is **4β-hydroxylation**, not 23/24
  ([Wang *et al.* 2012, *Mol Pharmacol* 81:498](https://doi.org/10.1124/mol.111.074393), PMC3310418),
  confirmed in vivo — rifampin raised 4β,25-(OH)₂D₃ by **60% in healthy volunteers** with a 10% fall
  in 1α,25-(OH)₂D₃ ([*JBMR* 28:1101](https://doi.org/10.1002/jbmr.1839), PMID 23212742).
- Mutagenesis on vitamin D turnover localised the effect to **residues 119, 120, 301, 305 and 479**
  ([Gupta *et al.* 2005, *JCEM* 90:1210](https://doi.org/10.1210/jc.2004-0966), PMID 15546903) — a
  residue set that overlaps both this repo's pocket definition and the positively selected codon 479.
- **VDDR3**: UniProt VAR_084337 **I301T**, "10-fold increased activity towards calcitriol… decreased
  activity for non-vitamin D substrates", altering **SRS-4**
  ([Roizen *et al.* 2018, *JCI* 128:1913](https://doi.org/10.1172/JCI98680), PMID 29461981).

**Bile acids are in KEGG's EC list but absent from UniProt's reaction block.** The primary evidence is
[Araya & Wikvall 1999, *BBA* 1438:47](https://doi.org/10.1016/s1388-1981(99)00031-1) (PMID 10216279):
"recombinant expressed **CYP3A4 was the only enzyme** that was active towards these bile acids and the
enzyme catalyzed an efficient **6α-hydroxylation of both taurochenodeoxycholic acid and lithocholic
acid**", TCDCA V_max 18.2 nmol/nmol P450/min, K_m 90 µM, **cytochrome b₅ required for maximal
activity**, and a K_m in human liver microsomes of 716 µM "which might give an explanation for the
limited formation of 6α-hydroxylated bile acids in healthy humans."

### 4.4 Downstream coupling — phase II and transporters

| CYP3A4 product | conjugating enzyme | status |
|---|---|---|
| 1′-hydroxymidazolam | **UGT2B4 + UGT2B7** (O-gluc), **UGT1A4** (N-gluc) | ESTABLISHED ([PMID 20713656](https://pubmed.ncbi.nlm.nih.gov/20713656/)) |
| 4-hydroxymidazolam | UGT1A4 | ESTABLISHED |
| lithocholic acid | **SULT2A1**, and *not* SULT2B1b or SULT1E1 | ESTABLISHED ([PMID 30918069](https://pubmed.ncbi.nlm.nih.gov/30918069/)) |
| **AFB1 → exo-8,9-epoxide** | **GSTM1, GSTP1, GSTT1** | ESTABLISHED — the classic CYP3A4→GST coupling, and the one Reactome models |
| 1,25-(OH)₂D₃ | **UGT1A4 ≫ UGT2B4, UGT2B7**; 25-O-glucuronide | ESTABLISHED ([PMID 18177842](https://pubmed.ncbi.nlm.nih.gov/18177842/)) |
| **6β-hydroxytestosterone** | — | **NO EVIDENCE FOUND.** It appears only as a CYP3A4 activity probe, never as a characterised phase-II substrate. Do not assert one. |
| **4β-hydroxycholesterol** | — | **NO EVIDENCE FOUND.** 25-OH-cholesterol → SULT2B1b and 27-OH → SULT2A1 are established; extending to 4β-OHC is **[INFERENCE]**, not evidence. |

**The regulon is real and it is PXR's.** [Sonoda, Xie, Evans *et al.* 2002, *PNAS*
99:13801](https://doi.org/10.1073/pnas.192448899): "PXR protects the body from hepatotoxicity of
secondary bile acids such as LCA by inducing… CYP3A… activation of PXR also increases… DHEA
sulfotransferase (SULT2A1)… and PAPSS2… **PXR serves as a master regulator of the phase I and II
responses**." Rifampin in primary human hepatocytes induces **UGT1A4 7.9-fold and UGT1A1 4.8-fold**
([PMC5869567](https://pmc.ncbi.nlm.nih.gov/articles/PMC5869567/)). **ABCC2/MRP2** is induced by FXR,
PXR *and* CAR through an **ER-8 element 440 bp upstream**
([Kast *et al.* 2002, *JBC* 277:2908](https://pubmed.ncbi.nlm.nih.gov/11706036/)). Note the
architecture: **CYP3A4 = ER6 + XREM at −7.8 kb; ABCB1 = DR4 at −8 kb; ABCC2 = ER-8 at −440 bp.**
Three element geometries, one receptor.

**The CYP3A4 / P-gp "shared substrate space" is CONTESTED, and this is worth stating plainly.** The
claim originates with [Wacher, Wu & Benet 1995](https://pubmed.ncbi.nlm.nih.gov/7619215/) and is
repeated as a "concerted barrier". Against it: [Lin & Yamazaki
2003](https://pubmed.ncbi.nlm.nih.gov/12489979/) argue intestinal P-gp is "unlikely to be
quantitatively important unless a very small oral dose is given" because it saturates;
[Darwich *et al.* 2010](https://pubmed.ncbi.nlm.nih.gov/21189140/) found the F_G effect "limited to
certain conditions… a very limited area of the parameter space matching very few therapeutic drugs";
[Estudante *et al.* 2013](https://pubmed.ncbi.nlm.nih.gov/24044638/) decoupled transport from
metabolism for verapamil. **And the two gradients run in opposite directions along the gut**
(**[INFERENCE]** from two primary papers): CYP3A falls proximal→distal (31 → 17 pmol/mg,
[Paine *et al.* 1997, *JPET* 283:1552](https://pubmed.ncbi.nlm.nih.gov/9400028/)) while P-gp rises
([Mouly & Paine 2003, *Pharm Res* 20:1595](https://pubmed.ncbi.nlm.nih.gov/14620512/)). Co-localised
at cell resolution, anti-correlated at organ resolution. **NO EVIDENCE FOUND** for any published
*quantitative* cheminformatic test of the overlap against a chance baseline — "broad overlap" has
never been formally measured.

### 4.5 Upstream regulators as a pathway

The *cis*-elements are in §1.5. The receptor-level summary, each ESTABLISHED unless marked:

- **PXR/NR1I2** — dominant xenosensor, binds both proximal ER6 and the XREM DR3/ER6; the deletion
  series goes 3-fold → **~50-fold** on hPXR co-transfection.
- **CAR/NR1I3** — two high-affinity motifs at ≈−7,720 and ≈−150, requiring cooperativity between
  promoter and XREM. ⚠️ **"the human CAR response elements also mediate trans-activation of CYP3A4 by
  hPXR"** ([Goodwin 2002](https://pubmed.ncbi.nlm.nih.gov/12130689/)) — so element mutagenesis
  **cannot separate CAR from PXR**.
- **VDR** — the intestine-specific arm, at the proximal ER6.
- **GR/NR3C1** — dual and biphasic: nanomolar dexamethasone gives a 3–4× component by inducing PXR
  and CAR, supramicromolar gives 15–30× by direct PXR ligation
  ([Pascussi 2001](https://pubmed.ncbi.nlm.nih.gov/11737189/)).
- **HNF4A** — the liver-specificity factor; required, mechanism CONTESTED (§1.6).
- **RXRA** — obligate partner at every element.
- **FXR/NR1H4** — **indirect and repressive, via SHP**: GW4064 cuts CYP3A4 mRNA 75% while raising SHP
  ~3× ([Zhang, Pan & Jeong 2015, *DMD* 43:743](https://doi.org/10.1124/dmd.114.062836), PMC4407707).
  **FXR splits the regulon: it represses the CYP and induces the efflux pump.**
- **SHP/NR0B2** — represses PXR, CAR *and* GR transactivation, and "reduced PXR recruitment of HNF4α
  and SRC-1 to the CYP3A4 chromatin"; reciprocally, activated PXR inhibits the SHP promoter — a
  feed-forward loop in which PXR relieves its own repressor.
- **Coactivators:** **NCOA1/SRC-1**, **PGC-1α**, and **NCOA6 + p300** with measured ChIP enrichment
  after rifampicin (NCOA6 +3.1× distal / +3.2× proximal; p300 +2.2× / +1.6×), alongside a chromatin
  signature of **H3K4me3 +2.6× proximal / +1.8× distal, H3K27me3 −30% / −60%, H3ac up at both**
  ([PMC5508193](https://pmc.ncbi.nlm.nih.gov/articles/PMC5508193/)). **NCOA3/RAC3** replaces SMRT on
  rifampicin.
- **Corepressor: NCOR2/SMRT** — ESTABLISHED, binds PXR's LBD via ID2; overexpression inhibits and
  silencing enhances PXR transactivation of the CYP3A4 promoter
  ([Johnson 2006](https://pubmed.ncbi.nlm.nih.gov/16219912/)).
- **MED1 and NCOR1 at CYP3A4: NO EVIDENCE FOUND.** Both appear only in Reactome's *generic*
  nuclear-receptor complexes. Do not assert them at this locus.

### 4.6 Where CYP3A4 sits in the DDI network

FDA's **clinical index** list (content current 2023-06-05): CYP3A sensitive index substrates
**midazolam, triazolam**; strong index inhibitors **clarithromycin, itraconazole**; moderate
**erythromycin, fluconazole, verapamil**; strong index inducers **carbamazepine, phenytoin,
rifampin**. The fuller example lists moved to the
[FDA examples page](https://www.fda.gov/drugs/drug-interactions-labeling/healthcare-professionals-fdas-examples-drugs-interact-cyp-enzymes-and-transporter-systems)
(updated 2026-05-29) and run to ~30 sensitive substrates, ~19 strong inhibitors, ~9 strong inducers.
Thresholds, verbatim: strong inhibitor **≥5-fold AUC increase**, moderate ≥2 to <5, weak ≥1.25 to <2;
strong inducer **≥80% AUC decrease**, moderate 50–80%, weak 20–50%.
⚠️ **The table is not a DAG** — ivacaftor, ticagrelor and lomitapide are substrates *and* weak
inhibitors; ritonavir is a strong inhibitor (in combination) *and* a weak inducer.

**The "~30–50% of drugs" number is three different quantities that got merged.** The reconciliation:

| source | number | what it actually counts |
|---|---|---|
| [Shimada 1994, *JPET* 270:414](https://pubmed.ncbi.nlm.nih.gov/8035341/) | CYP3A ≈ **30% of total hepatic P450** | immunoquantified **protein abundance**, not drugs — most circulating "30%" figures are this, misquoted |
| [Rendic & Guengerich 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4303333/) | **27% for drugs**, 13% for general chemicals | fraction of catalogued **reactions** |
| [Saravanakumar 2019, *Clin Pharmacokinet* 58:1281](https://pmc.ncbi.nlm.nih.gov/articles/PMC6773482/) | **40% of the top-200 prescribed; 64% of 2005–2016 approvals** (CYP3A4+3A5: 43% → 74%) | fraction of **drugs** with CYP3A4 as the major pathway |
| Zanger & Schwab 2013 | "70–80% of all drugs" | for **all ~12 CYPs collectively** — a very common misattribution |

**If one defensible modern number is needed: Saravanakumar 2019 — 40% of the top-200 prescribed and
64% of recent approvals — and note it says the CYP3A share is *rising* in new approvals.**

**What this licenses for the structure track.** Three things. (i) **The compound panel must not be
built from Reactome or from the PDB alone.** Reactome models CYP3A4 as purely xenobiotic; the local
CYP3A4 proxy set is **72 of 87 Type II heme-coordinating ligands and only 15 Type I** (MEASURED from
`data/processed/poses_scored_val87b.csv`). Both samples are biased away from the endogenous,
non-coordinating, lipophilic substrate space — steroids, oxysterols, bile acids, retinoids, PUFA
amides — which is what the enzyme was built for and which `FINDING_022` implies is where the model is
least calibrated. The `ranked_compounds` list in the JSON is built to span exactly that gap.
(ii) **CYP3A5 is a genuine catalytic sibling on the vitamin D and steroid reactions**, and the
mutagenesis that localises vitamin-D regiochemistry names **119, 120, 301, 305, 479** — three of which
are in this repo's pocket definition and one of which is a positively selected codon. That is a
ready-made hypothesis for where a regiochemistry-sensitive scorer should look. (iii) **Cytochrome b₅
is required for maximal bile-acid 6α-hydroxylation** — the one place the literature makes b₅
substrate-specific rather than generic. It still does not make b₅ worth co-folding (§6.9), but it is
the reaction to look at if anyone revisits that.

---

## 5. Cell types and states

`CYP3A4_EVOLUTION.md` §5.2 states the zonation and the intestinal gradient qualitatively. This
section puts numbers on both, corrects two things that are usually said, and supplies the ranked
co-expression list.

### 5.1 Zonation — pericentral, ESTABLISHED, and **steeper** than CYP2E1 or GLUL

The best quantitative source is
[Yakubovsky, Bahar Halpern & Itzkovitz 2026, *Nature*, "A spatial atlas of the healthy human liver
from live donors"](https://pmc.ncbi.nlm.nih.gov/articles/PMC13216088/), Visium HD, 3 live healthy
donors, 8 lobular layers (Supplementary Tables 4 and 6):

| gene | zone 3 : zone 1 |
|---|---|
| CYP1A2 | 58× |
| **CYP3A4** | **21.1×** |
| CYP2E1 | 11.7× |
| GLUL | 4.7× |
| **CYP3A5** | **1.8× the other way — PERIPORTAL** |

CYP3A4 is listed **first** in their pericentral landmark gene set, with a Visium portal-bias score of
−2.500 (5.66× pericentral) against CYP2E1's −2.081 and GLUL's −1.726.

⚠️ **Correction to a commonly stated premise: CYP3A4's gradient is steeper than CYP2E1's and GLUL's,
not shallower.** And **CYP3A4 and CYP3A5 are zonated in opposite directions**, which matters because
every non-isoform-specific "CYP3A" measurement averages them.

**Protein agrees.** Whole-slide automated quantification across 8 portality layers in n = 6 humans
gives a "moderate to strong pericentral signal confined to zone 3 and extending into zone 2", with
**CYP3A4 staining ~50% of the lobular surface** against CYP2E1's 72.2 ± 9.4% and glutamine
synthetase's 5–13% ([Albadry 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11137285/)); the same
paper notes "the strongest periportal to perivenous gradient was observed in humans". The classical
human IHC is [Ratanasavanh 1991, *Hepatology* 13:1142](https://pubmed.ncbi.nlm.nih.gov/1904834/):
"P-450 IIIA was restricted to **centrilobular and midzonal** hepatocytes in normal adult liver" — and,
importantly, **pan-lobular in fetal liver**, i.e. **zonation is acquired postnatally**.

⚠️ **The disagreement is real but methodological.** Every *in situ* method says pericentral. The two
dissociation-based scRNA studies with only a computational axis do not: MacParland 2018 puts CYP3A4
in an **interzonal** cluster (Cluster 15), and Aizarani 2019 assigns it to "periportal module 1"
(*P*adj 4.2 × 10⁻²⁰) — but **that same table also places ASS1, a canonical periportal gene, on the
central side**, so the axis fails its own control. *(A further caution: full-text retrieval found
Aizarani 2019 does not mention CYP3A4 at all, so the module assignment should be checked before being
cited.)*

Single-cell fraction-expressing across pooled CELLxGENE human liver datasets is directionally
consistent but compressed, as fraction-based metrics always are: **centrilobular 66.98%, unzonated
54.00%, periportal 51.51%, midzonal 46.28%** — a ~1.7× gradient in *fraction*, against a 21× gradient
in *level*.

### 5.2 Which cells, with numbers

**GTEx v8 median TPM:** **liver 335.3**, **small intestine terminal ileum 48.3**, pancreas 2.41, skin
~1.9, testis 1.03, adrenal 1.01, everything else < 0.8; all 11 non-cerebellar brain regions
0.027–0.095. **Liver : ileum = 6.95 : 1; liver : adrenal = 333 : 1.**
**HPA v24 consensus nTPM:** liver 3367, small intestine 754, duodenum 713, pancreas 25.3, colon 2.9,
**lung 0.0, placenta 0.0**; **Tau 0.93**, "Tissue enriched (liver)", single-cell **"Group enriched
(Hepatocytes, Enterocytes)"**, "Not detected in immune cells", cluster **"Hepatocytes — Bile
production & excretion"**. **Duodenum : colon ≈ 250×.**

⚠️ **Read the single-cell numbers with the ambient-RNA caveat.** Soup from lysed hepatocytes inflates
CYP3A4 in every other liver cell type: Andrews snRNA-seq reports Kupffer cells **33%** CYP3A4-positive
while Guilliams CITE-seq of the same organ reports macrophages **2.5%** and T cells **0.03%**. **Where
they disagree, the low number is the biological one.** HPA single-cell nCPM: hepatocytes **3988.6**,
enterocytes **2039.8**, then stellate 32.3 (123× down), cholangiocytes 28.5, Kupffer 16.4, endothelium
1.0, T/B cells ≤0.8. **ESTABLISHED: CYP3A4 is absent from cholangiocytes, stellate cells, Kupffer
cells, LSECs and lymphocytes at any pharmacologically meaningful level.**

**Enterocytes: the gradient is along the crypt–villus axis as much as along the gut.**
Single-cell across the maturation axis gives **crypt stem 1.26% → transit-amplifying 0.80% → mature
enterocyte 51.28%** — a **~40–60× induction at terminal differentiation** — and **colonocytes 0.79%**
([Burclaff 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9043569/)). Along the gut,
[Paine 1997, *JPET* 283:1552](https://pubmed.ncbi.nlm.nih.gov/9400033/) remains the reference:
median microsomal CYP3A **duodenum 31 → ileum 17 pmol/mg**, midazolam 1′-OH V_max **liver 850 >
duodenum 644 > jejunum 426 > ileum 68** pmol/min/mg (ileum **9.5× below duodenum**), with
**K_m ≈ 4 µM everywhere — so the gradient is enzyme abundance, not affinity**, and duodenal mucosa is
per-mg comparable to liver (CL_int 157 vs 200 µL/min/mg). CYP3A is **~80% of all immunoquantified
intestinal P450** ([Paine 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC2222892/)). ⚠️ A
single-cell atlas ordering (ileum ≈ duodenum > jejunum) **contradicts** this, but rests on 405 jejunal
cells from 1–2 donors and on a saturating fraction metric; **cite Paine 1997 for the gradient.**

**Development: the ratio inverts by three orders of magnitude.** In fetal hepatoblasts (n = 11,673),
**CYP3A7 69.88% of cells vs CYP3A4 0.09% (10 cells)** — and two further fetal atlases return **zero**
CYP3A4-positive cells. In adult hepatocytes on the same platform, **CYP3A4 72.90% vs CYP3A7 1.52%**.
⚠️ Popescu 2019, Segal 2019 and Wesley 2022 contain **zero CYP mentions** — excellent for hepatoblast
trajectories, useless for CYP3A4.

**Extrahepatic claims that are artefacts:** the "adrenal hepatocyte 45.7% CYP3A4+" signal traces to a
neuroblastoma atlas with liver contamination; kidney, prostate and placental "CYP3A" is **CYP3A5 or
CYP3A7**, not CYP3A4 (kidney proximal tubule CYP3A5 603.7 nCPM vs CYP3A4 1.2; prostate CYP3A5 229.7 vs
CYP3A4 0.3); lung alveolar type 1 and 2 are both **0.0**.

### 5.3 Ranked co-expression — ARCHS4, Pearson, verbatim

**Source:** ARCHS4 (Ma'ayan lab), Pearson correlation across the human RNA-seq sample matrix, queried
by `POST https://maayanlab.cloud/matrixapi/coltop` with `{"id":"CYP3A4","count":101}` — the same call
the ARCHS4 gene page's own JavaScript makes, which is why a plain page fetch returns nothing.
**I re-ran it independently and obtained byte-identical values, so it is deterministic and
reproducible.** Top 40, in rank order:

| # | gene | r | # | gene | r |
|---|---|---|---|---|---|
| 1 | **OTC** | 0.7276 | 21 | A1CF | 0.6135 |
| 2 | **ALDOB** | 0.7158 | 22 | TM4SF5 | 0.6120 |
| 3 | **ABCG8** | 0.6938 | 23 | F11 | 0.6095 |
| 4 | **APOC3** | 0.6853 | 24 | **NR1H4 (FXR)** | 0.6081 |
| 5 | **SLC2A2** | 0.6778 | 25 | C8A | 0.6058 |
| 6 | **CYP2C9** | 0.6524 | 26 | **UGT1A4** | 0.6041 |
| 7 | FLJ22763 | 0.6490 | 27 | RP11-400G3.5 | 0.5987 |
| 8 | **ABCG5** | 0.6418 | 28 | AFM | 0.5961 |
| 9 | SLC17A4 | 0.6405 | 29 | **UGT2A3** | 0.5949 |
| 10 | APOA4 | 0.6401 | 30 | F13B | 0.5948 |
| 11 | MTND4P20 | 0.6400 | 31 | PRAP1 | 0.5947 |
| 12 | PLA2G12B | 0.6378 | 32 | SMLR1 | 0.5938 |
| 13 | **FABP1** | 0.6341 | 33 | **CYP4F2** | 0.5934 |
| 14 | **CYP3A7** | 0.6337 | 34 | APCS | 0.5928 |
| 15 | F9 | 0.6290 | 35 | ENPP7 | 0.5927 |
| 16 | G6PC | 0.6245 | 36 | RP4-608O15.3 | 0.5922 |
| 17 | **CYP2C19** | 0.6204 | 37 | GBA3 | 0.5887 |
| 18 | PCK1 | 0.6172 | 38 | AGXT2 | 0.5878 |
| 19 | CFHR2 | 0.6162 | 39 | HAO1 | 0.5873 |
| 20 | MOGAT2 | 0.6138 | 40 | CFHR5 | 0.5843 |

Further down: **UGT1A1 #41, CYP3A5 #53, NR1I2/PXR #64, HNF4A #80.**

**Two things are worth noticing.** First, **the top neighbours are the hepatocyte/enterocyte identity
programme — urea cycle, gluconeogenesis, apolipoproteins, coagulation factors, sterol transporters —
not other CYPs.** CYP3A4's expression is a readout of *being a differentiated hepatocyte*, which is
also what the crypt-to-villus 40–60× jump says. Second, **FABP1 is #13 (r = 0.634)** — independent
support for the 2026 directed-substrate-transfer result in §6.6 — and **NR1H4/FXR is #24**, its own
upstream repressor.

**A second named source, and two dead ones.** STRING v12's **co-expression channel only** (`ascore`,
not the combined score) ranks: **CYP3A7 0.624, ALDOB 0.450, CYP2C9 0.446, ADH4 0.410, OTC 0.356,
APOC3 0.352, CYP2C19 0.326, SLC2A2 0.302, SULT2A1 0.291, APOB 0.290** — a different compendium and a
different metric, agreeing on OTC, ALDOB, APOC3, SLC2A2 and CYP2C9. ⚠️ **COXPRESdb is dead** (every
API pattern 404s), and **HumanBase's liver network returns CYP3A4 edges all at the 0.1065 mincut with
olfactory receptors and POTE paralogs in the top 40 — noise; do not use it.**

### 5.4 Induction and repression states, with magnitudes

**Thresholds first**, verbatim from ICH M12 (FDA final guidance, August 2024, pp. 32–33): a **strong**
inducer decreases a sensitive index substrate's AUC by **≥80%**, **moderate** ≥50 to <80%, **weak**
≥20 to <50%.

| state | magnitude | source |
|---|---|---|
| **rifampicin, PHH** | CYP3A4 **mRNA 23×, protein 12×, activity 13×** at 10 µM; E_max 12.3, EC₅₀ 0.847 µM | [Dixit 2007](https://pubmed.ncbi.nlm.nih.gov/17639026/) |
| **rifampicin, clinical** | **oral midazolam AUC −96%**; in the same subjects **oral −91% vs IV −52%**; oral CL 1.56 → 34.4 L·h/kg (~22×) with **bioavailability −88%** | [Backman 1996](https://pubmed.ncbi.nlm.nih.gov/8549036/), [Kirby 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3100903/), [Gorski 2003](https://pubmed.ncbi.nlm.nih.gov/12966371/) |
| St John's wort / hyperforin | oral MDZ CL **+109%**, oral AUC >−50% vs IV −20%; cyclosporine AUC −46% in transplant patients | [Wang 2001](https://pubmed.ncbi.nlm.nih.gov/11673747/), [Bauer 2003](https://pmc.ncbi.nlm.nih.gov/articles/PMC1894728/) |
| carbamazepine | PHH mRNA **8.3× (range 3.5–14.5)**; CAR-preferring, a weak PXR activator | [Sugiyama 2016](https://pubmed.ncbi.nlm.nih.gov/26711482/), [Faucette 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC4091905/) |
| dexamethasone | **biphasic**: 3–4× at nanomolar (GR → ↑PXR/RXRα) then 15–30× supramicromolar (direct PXR) | [Pascussi 2001](https://pubmed.ncbi.nlm.nih.gov/11737189/) |
| efavirenz | moderate; maraviroc AUC −50% vs rifampicin's −70%; simvastatin acid −58% | ICH M12 index moderate inducer |
| **pregnancy** | **third-trimester midazolam CL/F(unbound) +72%** (paired, n = 13); canonical figure ~+100%; ⚠️ dextromethorphan gives only **+35–38% and flat across all three trimesters** | [Hebert 2008](https://pubmed.ncbi.nlm.nih.gov/18288078/), [Tracy 2005](https://pubmed.ncbi.nlm.nih.gov/15696014/) |
| **IL-6 / inflammation** | PHH CYP3A4 mRNA to **<5% of control** at 10 ng/mL, n = 9 donors; **IL-1β and LPS equally potent** | [Aitken & Morgan 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC2171046/) |
| IL-6 blockade, clinical | tocilizumab cut **simvastatin AUC_last to 43%** of baseline at 1 week, n = 12 | [Schmitt 2011](https://pubmed.ncbi.nlm.nih.gov/21430660/) |
| surgery / acute inflammation | **CYP3A metabolic ratio −61% within the same person by day 3** after elective hip surgery (n = 30) — while CYP2B6 rose 120% | [Lenoir 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8247903/) |
| cancer | ERBT CYP3A function vs CRP **r = −0.64, P < 1e-5**; ~50% reduction with abnormal liver chemistry | [Rivory 2002](https://pmc.ncbi.nlm.nih.gov/articles/PMC2364233/), [Baker 2004](https://pubmed.ncbi.nlm.nih.gov/15623611/) |
| cirrhosis | midazolam total clearance **−41%**, oral bioavailability **38% → 76%** | [Pentikäinen 1989](https://pubmed.ncbi.nlm.nih.gov/2723115/) |
| **fasting** | **36 h fasting INCREASES midazolam clearance by 12% (p < 0.01)**; a 3-day high-fat diet does nothing | [Lammers 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6244726/) |

**Four mechanistic and methodological points worth carrying.**

1. **IL-6 dismantles the induction machinery, not just the gene** — it "rapidly and markedly decreases
   the expression of **PXR** and **CAR** mRNAs" while leaving AhR and GR untouched, and blunts both
   rifampicin and phenobarbital induction ([Pascussi 2000](https://pubmed.ncbi.nlm.nih.gov/10924340/)).
   The repression is **gp130-dependent but NOT JAK/STAT** (§1.6), plus NF-κB p65 sequestering RXRα.
2. **Within an individual, hepatic and intestinal induction are inversely related** — Spearman
   **r = −0.68, p < 0.0001, n = 52** (Gorski 2003). A person is induced strongly at the liver *or* the
   gut, not both. This is the best candidate explanation for the spread in CYP3A DDI magnitude, and it
   means **any fold-induction number is meaningless without the route**.
3. ⚠️ **A PXR EC₅₀ does not predict the clinical direction.** Ritonavir is a *more potent* PXR agonist
   than rifampicin (0.44 vs 0.87 µM) and induces CYP3A4 mRNA 19-fold in PHH **while suppressing
   activity**; clinically, inhibition wins and midazolam AUC *rises*. Any model reading PXR activation
   as a proxy for clinical CYP3A effect gets ritonavir backwards.
4. **CONTESTED — MASLD/MASH.** Two well-executed studies disagree: one reports **2.4-fold higher
   plasma midazolam in NASH, 4β-hydroxycholesterol 37–51% lower and hepatic CYP3A4 mRNA 69% lower**
   ([Woolsey 2015](https://pubmed.ncbi.nlm.nih.gov/26231377/)); the other, on banked microsomes from
   52 graded livers, finds **mRNA essentially unchanged and protein and activity not significant
   (p = 0.112, p = 0.180)** ([Fisher 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2769034/)).
   In vivo activity vs banked microsomes; both stand.
   ⚠️ And **fasting has the opposite sign to the FGF21 mechanism**: FGF21 *represses* CYP3A4 by
   lowering nuclear PXR, yet measured 36 h fasting *raises* activity 12%. **Nobody has reconciled
   these. NO EVIDENCE FOUND** for a large fasting effect in either direction.
   **Microbiome: mouse evidence only** — germ-free mice differentially express 112 hepatic genes
   dominated by xenobiotic metabolism, and indole-3-propionic acid is a PXR ligand in vivo, but **no
   human study perturbs the microbiome and measures CYP3A4 activity. INFERRED and plausible,
   unproven.**

**What this licenses for the structure track.** Three concrete things.
(i) **The compound panel should be drawn from the liver hepatocyte programme, not from the PDB.** The
top co-expression neighbours are urea-cycle, gluconeogenic, apolipoprotein and sterol-transport genes
(ABCG5/ABCG8 at ranks 3 and 8), which is the chemistry of a pericentral hepatocyte — sterols, bile
acids, lipids. Against that, the local CYP3A4 structural proxy set is **72 of 87 Type II
heme-coordinating ligands and only 15 Type I** (MEASURED from `data/processed/poses_scored_val87b.csv`).
The `ranked_compounds` list is built to cover that gap.
(ii) **CYP3A5 is the paralog worth docking against, and the reason is now sharper than sequence
identity.** It is co-expressed in the same liver at 156 TPM, **zonated in the opposite direction**,
and has six structures including two ligands (clotrimazole, azamulin) crystallised in both paralogs.
CYP3A7 is the second choice — 24.5 TPM in adult liver, two structures, and a pocket with nearly the
same residues that **cannot deform**, which makes the CYP3A4/CYP3A7 pair the cleanest available test
of whether a scorer is sensitive to lid *plasticity* rather than to sequence.
(iii) **Nothing here argues for modelling more than one conformational state per ligand.** Every
state in this section — induction, repression, zonation, development — changes **how much** enzyme is
present, not its fold.

---

## 6. The protein's partners

This is the section the structure track cares about most, so it is organised by **what kind of
evidence exists**, not by how famous the partner is. The short version: CYP3A4 has exactly two
partners with strong, quantitative, reproducible biochemistry (POR and CYB5A), one with a crystal
structure of *itself* and a disputed functional sign (PGRMC1), a genuine but under-measured class of
P450–P450 mixed oligomers, a cytosolic substrate-delivery partner discovered in 2026 (FABP1), and a
long tail of database-only associations that should not be treated as interactions.

**The single most important structural fact: no experimental structure of CYP3A4 in complex with any
protein partner exists.** All 131 PDB entries cross-referenced from UniProt
[P08684](https://rest.uniprot.org/uniprotkb/P08684.txt) are the soluble Δ-N-terminal catalytic
domain alone. There is no CYP3A4–POR, CYP3A4–b₅, CYP3A4–PGRMC1 or CYP3A4–CYP3A4 complex to fold
against or to validate against.

### 6.1 NADPH–cytochrome P450 reductase (POR) — obligatory, structurally uncharacterised in complex

| | |
|---|---|
| UniProt | [P16435](https://rest.uniprot.org/uniprotkb/P16435.txt), 677 aa, TM helix **22–42**, ER membrane, cytoplasmic side |
| Cofactors | FMN 86–91, 138–141, 173–182, 208; FAD 424, 454–457, 472–474, 478, 488–491, 676 |
| Structures **of POR alone** | 3QE2 (1.75 Å), 3QFC (1.80 Å), 3QFS/3QFT (1.40 Å), 3QFR (2.40 Å), 5EMN (2.20 Å), 5FA6 (2.30 Å), 3FJO (2.50 Å), 1B1C (1.93 Å), and the 2024–25 FMN-domain series 9EAX–9EB1 at 1.10–1.50 Å ([PDBe SIFTS](https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/P16435)) |
| Structure **with CYP3A4** | **none** |
| Evidence type | reconstitution kinetics, Nanodisc 1:1 complexes, fluorescence titration, mutagenesis of POR |

**ESTABLISHED.** A functionally homogeneous **1:1 CYP3A4:POR complex** in Nanodiscs is fully
catalytically competent ([Denisov *et al.* 2007, *JBC* 282:7066,
PMID 17213193](https://doi.org/10.1074/jbc.m609589200)) — oligomerisation is not required for
electron transfer. The **membrane lipid matters**: liver microsomal lipid enhances both the activity
and the redox coupling of co-localised POR–CYP3A4 in Nanodiscs relative to plain
phosphatidylcholine ([Liu *et al.* 2017, *FEBS J*
284:2302](https://doi.org/10.1111/febs.14129), PMID 28618157). POR's **hinge and linker are
functional elements, not spacers**: a minimum 8-residue linker between the membrane segment and the
FMN domain is required for a catalytically competent P450 complex
([doi:10.1016/j.jinorgbio.2024.112667](https://doi.org/10.1016/j.jinorgbio.2024.112667),
PMID 39032346, CYP2B4), and hinge mutations change the interaction
([Campelo *et al.* 2018, *IJMS* 19:3914](https://doi.org/10.3390/ijms19123914), PMID 30563285).

**CONTESTED — is POR only an electron donor?** MS footprinting of CYP2A6 ± the POR FMN domain found
*increased* surface exposure rather than interface protection, which the authors read as **long-range
allosteric modulation** of the P450 by its reductase
([*DMD* 54:100210](https://doi.org/10.1016/j.dmd.2025.100210), PMID 41494244). POR mutations also
alter **CYP1A2 substrate regioselectivity**
([*DMD* 2026](https://doi.org/10.1016/j.dmd.2026.100300), PMID 42105715). Neither experiment was done
on CYP3A4, so for CYP3A4 specifically this is **INFERRED by analogy**.

**Where the interface is, and why it matters here.** This repo has already measured it: the 12
published proximal-face interface positions (**K96, K115, K121, K127, K130, K141, K421, R422, K424,
K428, K440, R446**) sit **9.7 Å from the iron, 8.7 Å from the nearest ligand heavy atom and 10.7 Å
from the nearest modelled F/G residue**, across the porphyrin, and ligands are 100% distal
(`FINDING_024`, §"Cytochrome b5 as a co-folded partner").

### 6.2 Cytochrome b₅ (CYB5A) — the tightest-binding partner, and a 30-year mechanism argument

| | |
|---|---|
| UniProt | [P00167](https://rest.uniprot.org/uniprotkb/P00167.txt), 134 aa, TM **109–131**, ER, cytoplasmic side; soluble isoform 2 |
| Structures | 2I96 (NMR, soluble domain 1–108). **No complex with any P450.** |
| Affinity for CYP3A4 | in the **2.5–61 nM** band measured for nine P450s |

**ESTABLISHED, and the number is new.** Using a fluorescent site-directed mutant (Alexa 488-T70C-b₅),
[Kim, Kim, Tateishi & Guengerich 2021, *DMD* 49:850](https://doi.org/10.1124/dmd.121.000475)
(PMID 34330716) measured binding to P450s 1A2, 2B6, 2C8, 2C9, 2E1, 2S1, 4A11, **3A4** and 17A1 with
**K_d 2.5–61 nM**; only weak binding to 2D6 and none to 2A6. Two further results from that paper are
directly actionable here:

1. **CYP3A4's affinity for b₅ is *decreased* by a bound substrate or inhibitor** — the ligand and the
   redox partner are not independent.
2. Titrating POR into a CYP3A4•b₅ complex **partially restored** b₅ fluorescence with
   K_d,apparent **89 nM**, and gel filtration independently supported a **ternary
   CYP3A4–b₅–POR complex**. The authors conclude ternary complexes "are relevant in P450 3A4
   reactions as opposed to a shuttle mechanism."

**CONTESTED — ternary or mutually exclusive?** Directly against that,
[Urban, Perret & Pompon 2025, *Sci Rep* 15](https://doi.org/10.1038/s41598-025-09060-5)
(PMID 40628863) argue from rapid kinetics with redox-modified b₅ analogues for **mutually exclusive
CPR–P450 and b₅–P450 binary complexes**, with b₅ acting as both acceptor and donor while remaining
"partially confined" with specific CYP3A4 partners — a redox rather than conformational mechanism.
They entertain a ternary complex but do not require it. **Report both; do not pick.**

**CONTESTED — electron transfer or allosteric effect?** The classic result is that b₅ stimulates
CYP3A4 *without* donating an electron
([Yamazaki *et al.* 1996, *JBC* 271:27438](https://doi.org/10.1074/jbc.271.44.27438), PMID 8910324,
title: "Lack of electron transfer from cytochrome b5 in stimulation of catalytic activities of
cytochrome P450 3A4"), extended by
[Yamazaki *et al.* 2001, *JBC* 276:30885](https://doi.org/10.1074/jbc.M105011200) (PMID 11413149,
"Stimulation of cytochrome P450 reactions by apo-cytochrome b5"). Apo-b₅ reproducing holo-b₅'s effect
has been confirmed for CYP2C9
([Locuson *et al.* 2007, *DMD* 35:1174](https://doi.org/10.1124/dmd.107.014910), PMID 17446262). The
modern reading is that both happen and the balance is substrate-specific. A concrete recent example:
perillyl-alcohol activation of midazolam 1′-hydroxylation (2.7-fold on V_max/K_m) occurs **only in
recombinant CYP3A4 reconstituted with b₅**, not without
([*Pharmaceutics* 16:1581](https://doi.org/10.3390/pharmaceutics16121581), PMID 39771560).

**A small but useful detail:** injecting isatin onto a preformed CYB5A/CYP3A4 complex on an SPR chip
raised its dissociation rate by 30%, and 270 mM isatin raised K_d threefold
([Ershov *et al.* 2017, *Biomed Khim* 63:170](https://doi.org/10.18097/pbmc20176302170),
PMID 28414290) — i.e. a small molecule can modulate the protein–protein interface itself.

### 6.3 PGRMC1 — the only protein UniProt annotates as a CYP3A4 interactor, and its sign is disputed

| | |
|---|---|
| UniProt | [O00264](https://rest.uniprot.org/uniprotkb/O00264.txt), 195 aa, TM **25–43**, heme axial ligand **Tyr113** |
| Structure | **4X8Y**, 1.95 Å, cytosolic domain 72–195 with heme — of PGRMC1 alone, not of a complex |
| UniProt P08684 SUBUNIT | "Interacts with PGRMC1; the interaction requires PGRMC1 homodimerization." |

**ESTABLISHED that PGRMC1 dimerises through heme–heme stacking and that dimerisation is required for
its P450 interactions.** [Kabe *et al.* 2016, *Nat Commun* 7:11030](https://doi.org/10.1038/ncomms11030)
(PMID 26988023) solved the cytosolic domain at 1.95 Å, showed five-coordinate Fe–Tyr113, showed the
open heme face mediates dimerisation, and showed CO blocks dimerisation by occupying the sixth
coordination site. Haem-mediated dimerisation "is required for interactions with EGFR and cytochromes
P450." The abstract names no individual P450; UniProt's PGRMC1 entry lists **CYP1A1 and CYP3A4**.

**CONTESTED — direction of effect, and whether the interaction is direct.**
[Szczesna-Skorupa & Kemper 2011, *Mol Pharmacol* 79:340](https://doi.org/10.1124/mol.110.068478)
(PMID 21081644) titled their paper "Progesterone receptor membrane component 1 **inhibits** the
activity of drug-metabolizing cytochromes P450 **and binds to cytochrome P450 reductase**". They
report concentration-dependent inhibition of CYP2C2, CYP2C8 and **CYP3A4**, *partially reversed by
over-expressing POR*, while CYP51 (a sterol P450) activity *decreased* on PGRMC1 knockdown. Two
readings survive: PGRMC1 is a direct P450 modulator with opposite signs for sterol vs drug P450s, or
PGRMC1 competes for POR and the "CYP3A4 interaction" is largely indirect. The POR-rescue result
favours the second. **No structure of a PGRMC1–P450 complex exists.**

**Related MAPR proteins: NO EVIDENCE FOUND.** Searching PGRMC2, NENF/neudesin and CYB5D2 against
cytochrome P450 activity returned only lamprey progesterone-receptor phylogeny and a membrane-steroid-
receptor cancer review. There is **no evidence** that any MAPR protein other than PGRMC1 modulates
CYP3A4.

### 6.4 CYP3A4 oligomers and P450–P450 hetero-oligomers — real, and functionally consequential

**ESTABLISHED that CYP3A4 self-associates in a concentration-dependent, membrane-density-dependent
way, and that this changes measurable kinetics.** Davydov and colleagues, by luminescence resonance
energy transfer plus crosslinking
([*JBC* 290:3850, PMID 25533469](https://doi.org/10.1074/jbc.m114.615443)), report a
monomer ⇌ oligomer equilibrium and **mixed oligomers for CYP3A4/CYP3A5, CYP3A4/CYP2E1 and
CYP3A5/CYP2E1**, with the crosslinking/LRET distances supporting an interface "similar to that
observed in the crystallographic dimers". Downstream consequences, all measured:

| observation | oligomeric CYP3A4 | monomeric (Nanodisc) | source |
|---|---|---|---|
| dithionite reduction kinetics | three-exponential | near-monoexponential and fast | [PMID 16229479](https://doi.org/10.1021/bi0509346) |
| fraction reducible by the BMR flavin domain | **50%** (proteoliposomes 55 ± 6%) | almost complete | [PMID 20026040](https://doi.org/10.1016/j.bbabio.2009.12.008) |
| αNF activation of 7-BFC | present, requires L/P ≤ 140 | undetectable at L/P ≥ 750 | [*Biochem J* 453:219, PMID 23651100](https://doi.org/10.1042/bj20130398) |

**The functionally most interesting hetero-oligomer is CYP3A4–CYP2E1.** Alcohol-induced CYP2E1 in
human liver microsomes *increases* CYP3A4's metabolic rate and *attenuates* its homotropic
cooperativity ([*ABB* 698:108677, PMID 33197431](https://doi.org/10.1016/j.abb.2020.108677)) — a
physiological perturbation changing a CYP3A4 kinetic signature. Conversely,
CYP3A4 **did not** affect CYP2C9 metabolism at saturating POR in an immobilised-P450 SPR system
([*DMD* 44:1364, PMID 26961240](https://doi.org/10.1124/dmd.115.067637)) — a clean negative worth
keeping.

Two conformational-plasticity numbers from the same programme, useful as priors: pressure-induced
transition at **P½ = 1.45 ± 0.33 kbar, K_eq° = 0.13 ± 0.06**
([PMID 17555301](https://doi.org/10.1021/bi602400y)) and **ΔV° = −36.8 ± 5.0 mL/mol** with only
**15%** of the enzyme in the pressure-promoted conformation at ambient pressure
([PMID 27074675](https://doi.org/10.1016/j.bpj.2016.02.026)).

**NO EVIDENCE FOUND** for a CYP3A4–CYP1A2 or CYP3A4–HMOX1 complex specifically. The BRET work on
heme oxygenase-1 heteromers covers CYP1A1, CYP1A2 and CYP2D6
([PMID 33148696](https://doi.org/10.1074/jbc.ra120.015911),
[PMID 33394027](https://doi.org/10.1042/bcj20200768)) and explicitly notes "not all P450s" form them;
CYP3A4 was not among those tested.

### 6.5 The membrane — the partner that is always present and never modelled

**ESTABLISHED.** CYP3A4 is a single-pass ER membrane protein (UniProt P08684) whose F′–G′ face is
bilayer-embedded; the DHEA-S structure 8GK3 resolves a steroid sulfate sitting **on** that face. In
practice:

- **Anionic lipid is used and matters.** Recent mechanistic work embeds CYP3A4 in nanodiscs of
  **1:1 POPC:POPG** ([*Mol Pharmacol* 104:97, PMID 37536953](https://doi.org/10.1124/molpharm.123.000698)).
- **Lipid identity changes measured ligand binding.** POPC nanodiscs give "much higher binding of
  fragments" than DMPC or DPhPC
  ([*Chem Biol Drug Des* 105:e70080, PMID 40087816](https://doi.org/10.1111/cbdd.70080)).
- **Native liver lipid beats synthetic PC** for POR–CYP3A4 coupling
  ([PMID 28618157](https://doi.org/10.1111/febs.14129)).
- The **allosteric midazolam site is located in the F′/G′ region** in lipid nanodiscs
  ([*Biochemistry* 59:1074, PMID 31961139](https://doi.org/10.1021/acs.biochem.9b01001)).

**NO EVIDENCE FOUND, in the searches run,** for a measured membrane-insertion depth or heme-tilt
angle for CYP3A4 specifically, or for a specific lipid headgroup requirement beyond "anionic helps,
native lipid helps".

### 6.6 Cytosolic partners: FABP1, and a hint about GSTs

**NEW and ESTABLISHED (2026).**
[McCarty, Proffett & Guengerich, *J Med Chem* 69:7313–7328](https://doi.org/10.1021/acs.jmedchem.5c03754)
(PMID 41778753, PMC13036779) show that liver fatty-acid-binding protein **FABP1 binds the CYP3A4
substrates diazepam and sulfinpyrazone with K_d < 2.5 µM**, weakly attenuates CYP3A4 oxidation of
both, and — the important part — **kinetic modelling favours a model of directed substrate transfer
to CYP3A4**, not simple sequestration. Reconstituting human liver microsomes with cytosol attenuated
diazepam metabolism but *stimulated* sulfinpyrazone oxidation, "possibly due to GST A1-1".

This is the first well-controlled evidence that a **soluble protein hands substrate to CYP3A4**, and
it is a different mechanism class from everything above. A broader review of intracellular lipid
binding proteins is careful to say the mechanistic contribution is "largely undefined"
([*DMD* 51:700, PMID 37012074](https://doi.org/10.1124/dmd.122.001010)).

### 6.7 Chaperones, ERAD machinery, and what is *not* there

- **ERAD machinery — ESTABLISHED** (developed in §8): gp78/AMFR and CHIP/STUB1 E3 ligases with
  UBC7/UBE2G2 and UbcH5a ([PMID 19103148](https://doi.org/10.1016/j.abb.2008.12.001);
  [PMID 21270532](https://doi.org/10.4161/cbt.11.6.14834)).
- **Folding chaperones — NO EVIDENCE FOUND.** Searches for Hsp70/HSPA, Hsp90/GRP94, calnexin, BiP or
  DNAJ binding CYP3A4 returned nothing CYP3A4-specific; the closest hit is a general ER
  quantity-control review that names CYP3A4 only as an ERAD substrate
  ([PMID 27743014](https://doi.org/10.1007/s00232-016-9931-0)). This is a real gap in the literature,
  not a gap in the search alone — but treat it as unproven either way.

### 6.8 Database-only associations — read these as noise, not as partners

- **IntAct/UniProt** list only **PGRMC1, FANCG and LIAT1** as CYP3A4 protein interactors. FANCG and
  LIAT1 come from high-throughput screens with no follow-up and no mechanism; I found nothing that
  makes either biologically interpretable. **NO EVIDENCE FOUND** for function.
- [Ershov *et al.* 2025, *Mol Biol Rep* 52:207](https://doi.org/10.1007/s11033-025-10341-5)
  (PMID 39937310) aggregated 287 PPIs and 246 unique interactors across 33 microsomal P450s and found
  **UBC, PGRMC1 and FANCG** to be the most frequent common interactors — i.e. exactly the three that
  a degradation tag, a heme chaperone and a screen artefact would produce. 20 "moonlighting" proteins
  appear in the CYP3A4 subinteractome.
- [**STRING**](https://string-db.org/api/json/interaction_partners?identifiers=CYP3A4&species=9606)'s
  top 50 partners are dominated by UGTs, other CYPs and GSTs at combined scores >0.95 — but the
  *experiments* channel for almost all of them is 0.045, and the score is carried by "database"
  (pathway co-membership) and "textmining". The only STRING partners with a non-trivial experimental
  channel are **UGT2B7 (0.626), CYP3A5 (0.378), CYP2C9 (0.317), CYB5A (0.314) and POR (0.294)**.
  **This is a worked example of why a STRING score is not an interaction.**

### 6.9 Ranked partners, and what the repo's own measurements say about co-folding them

| rank | partner | why it matters | complex structure? | co-fold verdict |
|---|---|---|---|---|
| 1 | **a second copy of the query ligand** | the only "partner" with CYP3A4 structures showing it (2V0M, 4K9T, 8SO1, 8SO2, 7SV2, 8SG5, 8GK3) and the only one that touches the F/G roof | yes, many | ⚠️ **SPENT AND NEGATIVE — see the box below** |
| 2 | **membrane lipid (POPC / POPG)** | the F′–G′ face is bilayer-embedded; lipid identity changes measured binding | no CYP3A4 complex; 8GK3 puts DHEA-S on that face | **test next, on a LOWERED prior** |
| 3 | **CYB5A** | K_d 2.5–61 nM, ternary with POR, and substrate binding changes the affinity | **no** | do not spend — proximal face, 9.7 Å from Fe across the heme (`FINDING_024`) |
| 4 | **POR / POR FMN domain** | obligatory electron donor, possible allosteric modulator | **no** | do not spend — same geometry argument |
| 5 | **CYP3A4 itself (homodimer)** | the oligomer interface *is* the F′–G′ surface that carries the peripheral site | no | speculative; interface unvalidated |
| 6 | **CYP2E1** | mixed oligomer that measurably raises CYP3A4 rate and kills its cooperativity | no | speculative |
| 7 | **PGRMC1** | the only UniProt-annotated protein interactor; heme-dimer structure exists (4X8Y) | no complex | no — sign disputed, may act via POR |
| 8 | **FABP1** | directed substrate transfer, K_d < 2.5 µM for two CYP3A4 substrates | no | no — cytosolic, does not shape the pocket |

**What this licenses for the structure track.** The repo has already measured, in `FINDING_024`, that
the b₅/POR proximal-face interface sits **9.7 Å from the iron and 10.7 Å from the nearest modelled
F/G residue, across the porphyrin**, that protein accuracy at the pocket (0.73 Å) is uncorrelated
with ligand error (ρ = +0.03), and that no CYP3A4–b₅ complex exists to validate against. Nothing in
this section overturns that — if anything §6.1–6.2 strengthen it, because the partners with the best
biochemistry are precisely the ones bound on the wrong face. **The partners worth co-folding are the
ones that touch the F/G roof, and there are only two: a second ligand copy, and lipid.** Both are
cheap: both are just extra ligand entities in a Boltz YAML, neither needs a second chain, and both
are aimed squarely at residues 211–216 and Phe304, which `FINDING_028` measured as doing 83% of the
pose exclusion.

> ### ⚠️ Update, 2026-09: rank 1 was run and is REFUTED (`FINDING_031`)
>
> The second-ligand co-fold has been executed off the back of this list — **48 scored OpenProtein jobs, 49 including the feasibility probe** —
> and it is a **negative**. The second copy landed in the **active site on 12 of 15 ligands** rather
> than in the peripheral groove, and **the first copy got worse by −0.0664 LDDT-PLI, with rotation
> error moving +4.45° in the wrong direction**.
>
> The geometric reasoning in this section held up unusually well, which is why the result is worth
> reading carefully rather than just discarding: the model placed the second copy at a **median
> 8.81 Å from the iron**, **reproducing 2V0M's second ketoconazole unprompted**, and **45 of 60
> second copies did contact the groove residue set identified in §7.2**. The prediction about *where*
> a second ligand would go was right. The prediction that it would *help* was wrong — because it
> engaged those residues **from inside the cavity**, displacing the first copy, rather than from the
> surface.
>
> **Consequence for rank 2.** Membrane lipid is now the only untested partner on this list, and its
> prior is **LOWERED**, not merely unchanged. The one thing that has now been measured is that
> **adding mass on the distal side makes the pose worse**. The surviving hypothesis is narrower than
> the one written above: that lipid restrains residues 211–216 **from outside** the cavity, where a
> second ligand demonstrably acts from inside. It should be run as a falsification test, not as an
> expected gain.
>
> **Consequence for how this document should be used.** This is the first experiment this map has
> caused, and it went the other way. The `status` field now attached to every node and edge in the
> JSON exists for exactly this reason: so that a CONTESTED or INFERRED claim is visible at the point
> where someone is deciding whether to spend compute on it, not only here in the prose.

---

## 7. The second site and cooperativity

`CYP3A4_EVOLUTION.md` §4.4–4.5 already argues that CYP3A4 has a *gradient* of sites rather than one
pocket, and gives the Denisov/Sligar 1-/2-/3-testosterone decomposition. This section does not repeat
that. It adds **primary measurements over the whole PDB** made for this document, and it separates
what is established from what is routinely mis-cited.

**Method for everything marked MEASURED below:** the complete PDB lists for CYP3A4 (P08684, 122
entries), CYP3A5 (P20815, 6) and CYP3A7 (P24462, 2) were pulled from UniProt; ligand stoichiometry
came from the RCSB API for all 130; every mmCIF was downloaded and, per ligand copy, the distance to
the heme Fe, residues within 4.5 Å, occupancy, B-factor, `label_alt_id` and inter-copy minimum atom
distance were computed.

### 7.1 Every CYP3A structure with more than one copy of a ligand — MEASURED

Criterion: drug-like non-polymer ligand (MW ≥ 110, heme and curated cryoprotectants/buffers
excluded), copies **per protein chain** > 1.

| PDB | protein | ligand | copies/chain | Å | year | where the copies are (Fe → closest ligand atom) |
|---|---|---|---|---|---|---|
| **2V0M** | 3A4 | ketoconazole | **2** | 2.80 | 2007 | #1 active site, Fe–N **2.1 Å**; #2 at **9.3–9.8 Å**, upper cavity / channel mouth. **Not** the 1W0F peripheral site |
| **4K9T** | 3A4 | desoxyritonavir analogue | **2** | 2.50 | 2013 | Fe **2.2 Å** and **5.5 Å**, both inside; ⟨B⟩ 57 / 68 |
| **4K9U** | 3A4 | desoxyritonavir analogue | **2** | 2.85 | 2013 | Fe **2.4 Å** and **5.9 Å**; ⚠️ ⟨B⟩ **154 and 200 Å²** — treat as unreliable |
| **8SO1** | 3A4 | caffeine | **3** | 2.05 | 2023 | active site 3.7 Å (occ 0.67) · channel 13.5 Å (occ 0.91, ⟨B⟩ 117) · peripheral 19.7 Å (occ 0.87) |
| **8SO2** | 3A4 | caffeine | **6** | 2.15 | 2023 | active-site **stack of three** at 3.4 / 7.0 / 10.5 Å · channel 13.1 Å · peripheral ×2 at 19.9 and 21.1 Å |
| *7LXL* | *3A4* | *tethered testosterone dimer* | *1 molecule = 2 steroid units* | *2.75* | *2021* | *3.7 Å; ⟨B⟩ 172* |
| **8SG5** | **3A5** | clotrimazole | **3** | 2.80 | 2023 | Fe–N 1.9 Å · 9.4 Å · 12.9 Å, reproduced in **all 8 chains** to ±0.2 Å |
| **7SV2** | **3A5** | azamulin | **2** | 2.46 | 2022 | 4.2 Å and 9.8 Å, all 4 chains |
| **8GK3** | **3A7** | DHEA-S | **2–5** | 2.60 | 2023 | active site 4.0 · channel 8.5 · surface 24 / 27 / 35 Å, most spanning two chains |

**All of these are genuine simultaneous occupancy, not deposited alternates.** Minimum inter-copy
atom distances: 2.90 Å (4K9T), 3.12 Å (4K9U, 8SO2), 3.22/3.44 Å (8SO2), 3.30 Å (7SV2), 3.50/3.81 Å
(8SG5), 3.58 Å (2V0M). **Zero atom pairs below 2.5 Å in any entry** — these are van der Waals and
π-stacking separations.

**Two measured negatives worth as much as the table.**

1. **No CYP3A4 structure has two *different* drug-like ligands bound simultaneously. Zero out of
   122.** The only apparent multi-ligand hits are 6DAJ (inhibitor + a buffer **glutamine** 23 Å away
   on the surface) and 4D6Z / 9YK4 (inhibitor + **imidazole** from buffer ligating the Fe). The whole
   structural case for *heterotropic* two-ligand occupancy rests on **homotropic** structures plus
   one covalently tethered dimer.
2. **8SPD vs 8SG5 is a controlled contrast.** Same ligand (clotrimazole), same laboratory, adjacent
   years: **1 copy in CYP3A4, 3 copies in CYP3A5**. Multi-occupancy is not a generic CYP3A property.

### 7.2 The peripheral site: reproducible, and its residue list needs correcting

**ESTABLISHED that the groove is a reproducible ligand-binding site** — five entries, two chemotypes,
two laboratories, 19 years apart (MEASURED):

| PDB | ligand | Fe → closest atom | occ | ⟨B⟩ |
|---|---|---|---|---|
| 1W0F | progesterone | **17.0 Å** | 1.00 | 45.1 |
| 5A1P | progesterone | 17.1 Å | 1.00 | 80.0 |
| 5A1R | progesterone | 17.2 Å | 1.00 | 75.6 |
| 8SO1 | caffeine | 19.7 Å | 0.87 | 70.5 |
| 8SO2 | caffeine ×2 | 19.9 / 21.1 Å | 0.97 / 0.70 | 61.8 / 69.6 |

**Residues within 4.5 Å (MEASURED):** progesterone in 1W0F contacts **Arg212, Phe213, Asp214,
Asp217, Phe219, Phe220, Ile238, Val240** (5A1P/5A1R add Cys239); caffeine in 8SO1/8SO2 contacts
**Asp217, Phe219, Phe220, Ile238, Cys239, Val240**.

Against the residue list this site is usually given:

- **Phe219, Phe213, Asp214 — confirmed.** Phe219 is the stacking platform (caffeine stacks parallel
  at ~3.8 Å).
- **Asp217 — confirmed and usually omitted.** Its carboxyl is 3.5 Å from the caffeine carbonyl in
  both caffeine structures; it is arguably the most important polar residue at the site.
- **Leu211 — REFUTED as a peripheral-site residue.** It never contacts progesterone or caffeine at
  the peripheral site in any of the five entries. It is an **active-site** residue (it contacts
  ketoconazole #1 in 2V0M).
- **Glu218 — REFUTED.** Not a contact in any entry.
- **Phe215, Phe241, Phe304 — REFUTED as peripheral-site residues**, confirmed as the Phe cluster the
  site sits above and as channel/cavity residues.
- **Add Phe220, Ile238, Cys239, Val240.**

**There is a third site, not two.** The access channel is itself a resolved locus: fluorol in 8DYC
sits **12.3 Å** from Fe (occ 0.88) contacting Y53, F57, D76, R106, P107, F108, F215, F220, I223,
**T224**, P227, I230, A370, M371, R372, E374, and caffeine occupies the same channel in 8SO1/8SO2.
So: **active site 0–11 Å, channel ~12–14 Å, peripheral groove ~17–21 Å.**

### 7.3 A citation error that propagates

**Harlow & Halpert 1998, *PNAS* 95:6636** ([PMID 9618464](https://doi.org/10.1073/pnas.95.12.6636))
is the classic mutagenesis cited in support of the peripheral site. The L211F/D214E ("CYP3A5-like")
double mutant increased testosterone and progesterone 6β-hydroxylation at low [S], **decreased
α-naphthoflavone stimulation**, and **abolished homotropic cooperativity with both steroids**. That
result is ESTABLISHED.

**But the paper argues for the opposite mechanism.** Its own abstract states the hypothesis that "the
most likely location of effector binding is **in the active site along with the substrate**" and that
the substitutions were designed "to mimic the action of the effector by **reducing the size of the
active site**". L211 and D214 came from a homology model of the *active site*, six years before
1W0F — and the measurement in §7.2 shows **L211 is not at the peripheral site at all**. Harlow &
Halpert 1998 is evidence for mechanism (a), simultaneous occupancy of one cavity. It is routinely
mis-cited for mechanism (b) because D214 happens to appear in both residue lists. **REFUTED as
peripheral-site evidence.**

**Genuine peripheral-site mutagenesis does exist, and it is recent:**

- **F219A** ([*JBC* 299:105117, PMID 37524132](https://doi.org/10.1016/j.jbc.2023.105117)):
  progesterone S₅₀ **103 → 158 µM**, n_H **1.65 → 1.25**, high-spin **70% → 45%**; for caffeine it
  converts one-site to two-site binding (K_d 0.18 mM and 15 mM). The strongest solution evidence that
  this site modulates cooperativity.
- **F213A/F213S/F213Y in POPC Nanodiscs**
  ([*Biochemistry* 58:1411, PMID 30785734](https://doi.org/10.1021/acs.biochem.8b01268)): reduced
  progesterone activation of carbamazepine epoxidation, worst in F213A — **in monomeric enzyme**.

### 7.4 The counter-evidence, given straight

1. **The site lies on a crystal-packing interface, conceded by its proponents.** The 2023 caffeine
   paper states that "the peripheral site lies at the interface" of the crystallographic dimer and
   that caffeine association there "might be assisted by crystal packing or dimerization of CYP3A4",
   and that whether the crystallographic dimer resembles a natural one "remains to be established".
   The nuance: the *proximal* peripheral caffeine is 6.4 Å from the symmetry-related ligand (too far
   for packing to explain) while the *distal* one is 4.4 Å with a symmetry Glu244 at 3.5 Å — matching
   the measured occupancies, **0.97 proximal vs 0.70 distal**.
2. **Williams *et al.* hedged in the original**: "may be involved in the initial recognition of
   substrates or allosteric effectors" ([*Science* 305:683](https://doi.org/10.1126/science.1099736)).
3. **Progesterone has never been crystallised in the active site.** 5A1P/5A1R co-crystallised under
   two conditions and both "contained only one PRG molecule bound to the previously identified
   peripheral site" ([*Biochemistry* 54:4083](https://doi.org/10.1021/acs.biochem.5b00510)). Read one
   way that is replication; read the other, three attempts to crystallise the complex of a canonical
   cooperative substrate have never put it where it reacts.
4. **B-factors are never low**: progesterone ⟨B⟩ 45 / 80 / 76.
5. **5TE8 is the sharpest mechanistic objection.** Midazolam binding causes "a dramatic conformational
   switch in the F-G fragment… resulting in a collapse of the active site cavity"
   ([*PNAS* 114:486](https://doi.org/10.1073/pnas.1616198114)); MEASURED, in all three chains
   midazolam is in the **active site** (3.0–3.5 Å from Fe) and **the peripheral site is empty**, and
   its contacts include L216, D217, P218 — the same F-helix segment that forms the peripheral groove.
   **The peripheral site and the productive active site compete for the same F/F′ segment.**
6. **Solution localisation of a ligand to the peripheral site: NO EVIDENCE FOUND.** No ITC, NMR,
   crosslinking or HDX experiment places a ligand there. What exists is indirect and **INFERRED**:
   the F219A perturbations, and stopped-flow on bromoergocryptine supporting "a three-step BEC binding
   model according to which the drug binds first at a peripheral site without perturbing the heme
   spectrum and then translocates into the active site"
   ([*JBC* 287:3510, PMID 22157006](https://doi.org/10.1074/jbc.m111.317081)).

**Verdict: CONTESTED, and it should stay contested.** The groove is real and reproducible; whether
occupancy of it is functional in a membrane rather than in a lattice is not settled by any structure,
and the best functional evidence is mutational — which cannot distinguish "this residue lines an
allosteric site" from "this residue controls F′-helix dynamics".

### 7.5 Homotropic cooperativity

**ESTABLISHED that it occurs.** In purified recombinant CYP3A4, positive cooperativity for
**testosterone, 17β-estradiol, amitriptyline and aflatoxin B1**, with the same patterns in a
CYP3A4–reductase fusion protein
([Ueng *et al.* 1997, *Biochemistry* 36:370](https://doi.org/10.1021/bi962359z), PMID 9003190). The
best-documented single Hill coefficient is progesterone, **n_H = 1.65 ± 0.07, S₅₀ = 103 ± 3 µM**
(spectral titration, Δ3-22 CYP3A4, *JBC* 2023). 7-Benzyloxyquinoline is cooperative in human liver
microsomes and **attenuated by CYP2E1 enrichment**
([PMID 33197431](https://doi.org/10.1016/j.abb.2020.108677)). Carbamazepine and phenanthrene are
documented ([PMID 21992114](https://doi.org/10.1021/bi200924t);
[PMID 8204577](https://doi.org/10.1021/bi00187a009)). Midazolam is homotropic **on regioselectivity,
not rate** (§7.7).

⚠️ **A methodological warning that matters for a negatives-publishing project:** Hill coefficients
from *steady-state turnover* and from *spectral binding titration* are not the same measurement, and
§7.7 shows the two can move in opposite directions for the same ligand pair. Hill coefficients for
diazepam, pyrene, nifedipine and felodipine were **NOT FOUND within the search budget** — unsearched,
not negative.

### 7.6 Heterotropic activation — and where it fails to replicate

| pair | direction / magnitude | system | source |
|---|---|---|---|
| αNF → phenanthrene | ↑V_max, **K_m unchanged**; reciprocal, neither competitively inhibits the other | purified | [*Biochemistry* 33:6450](https://doi.org/10.1021/bi00187a009) — the founding experiment for simultaneous occupancy |
| αNF → aflatoxin B1 | **↑ 8,9-epoxidation AND ↓ 3α-hydroxylation, same enzyme, same effector** | purified + fusion | [PMID 9003190](https://doi.org/10.1021/bi962359z) |
| **αNF → 7-BFC** | **>2× in HLM; NO ACTIVATION in Supersomes or Baculosomes**; restored to **225%** at L/P = 140; **undetectable at L/P ≥ 750** | density-controlled | [*Biochem J* 453:219, PMID 23651100](https://doi.org/10.1042/bj20130398) |
| αNF → 7-BFC, second failure | **eliminated** by CYP2E1 enrichment of HLM | HLM ± CYP2E1 | [PMID 33197431](https://doi.org/10.1016/j.abb.2020.108677) |
| progesterone → carbamazepine epoxidation | activation, **reduced by F213A/S/Y** | POPC Nanodiscs, **monomeric** | [PMID 30785734](https://doi.org/10.1021/acs.biochem.8b01268) |
| progesterone → midazolam | shifts 1′-OH : 4-OH ratio, **no rate increase** | Nanodiscs | [PMID 34015213](https://doi.org/10.1021/acs.biochem.1c00161) |
| carbamazepine ↔ midazolam | NMR T₁ + docking put the two **stacked** in the active site | purified / NMR | [PMID 21992114](https://doi.org/10.1021/bi200924t) |
| atorvastatin lactone → dronedarone | ARVL is the dominant effector; MD implicates **F213 and F219** | Nanodiscs, monomeric | [PMID 29200287](https://doi.org/10.1021/acs.biochem.7b01012) |
| VU0448187 → midazolam | **>100% activation in vitro**, "diminished effect" toward in vivo | microsomes → hepatocytes → in vivo | [PMID 24003250](https://doi.org/10.1124/dmd.113.052662) |
| testosterone → nifedipine | **NOT FOUND within the search budget** | | |
| quinidine → diclofenac | **NOT FOUND.** Searches surfaced only the analogous dapsone → CYP2C9 activation. **Do not cite this pair as CYP3A4 without verifying it.** | | |
| flavonoids (quercetin, galangin) | **NOT FOUND within the search budget** | | |
| clean clinical *in vivo* activation | **NO EVIDENCE FOUND.** Treat "heterotropic activation is largely an in-vitro phenomenon" as CONTESTED but currently unrefuted. | | |

**The αNF/7-BFC row is the most consequential entry in this section.** The prototypical heterotropic
activation of the most-studied human drug-metabolising enzyme **does not occur** in the commercial
recombinant systems most laboratories use, and its magnitude is a monotonic function of enzyme
surface density in the membrane. That is a reproducibility failure with a mechanism attached, and it
is the reason "CYP3A4 cooperativity" should never be reported without the lipid-to-protein ratio and
the enzyme source.

### 7.7 The mechanism debate, laid out fairly

**(a) Multiple simultaneous occupancy of one cavity — ESTABLISHED, and now structurally proven.**
Shou 1994's non-competitive reciprocal kinetics and Ueng 1997's opposite-direction effects on two
aflatoxin sites, settled structurally by §7.1: 8SO2 stacks three caffeines parallel to the heme.
⚠️ But the authors of that structure judge the six-caffeine complex "likely to be nonfunctional
because, by clogging the active site, the caffeine trimer could preclude product dissociation".
**A multi-occupancy crystal structure is not by itself evidence of a catalytically relevant state.**

**(b) A distinct peripheral allosteric site — CONTESTED.** See §7.2–7.4.

**(c) Conformational heterogeneity — ESTABLISHED as a contributor**, and the 2020s resolution favours
**induced fit over conformational selection**: stopped-flow on four irreversible inhibitors in
monodisperse anionic Nanodiscs gave single-site binding with multi-step kinetics, globally fit by IF
([*Mol Pharmacol* 104:97, PMID 37536953](https://doi.org/10.1124/molpharm.123.000698)). Atkins'
review remains the canonical statement that "the functional effect of any specific ligand as an
activator or inhibitor can be substrate dependent"
([*Annu Rev Pharmacol Toxicol* 45:291](https://doi.org/10.1146/annurev.pharmtox.45.120403.100004)).

**(d) Oligomerisation — CONTESTED, and this is the live fault line.** Davydov's evidence is that αNF
activation of 7-BFC *requires* oligomer-competent surface density, that LRET shows a sigmoidal
density dependence which αNF abolishes, and that a physiological perturbation (alcohol-induced
CYP2E1) changes CYP3A4 cooperativity in human liver microsomes. Sligar's evidence is that
**homotropic and heterotropic allostery both persist in strictly monomeric CYP3A4** in Nanodiscs —
the 2018–2022 series is titled, in part, "…Mediated by **Monomeric** CYP3A4".

*How to state it fairly:* the two programmes measure different observables. Davydov shows that
**αNF activation of turnover** is oligomer-dependent and vanishes in dilute membranes; Sligar shows
that **midazolam/progesterone/carbamazepine regioselectivity switching and testosterone spin-state
cooperativity** persist in a monomer. The synthesis the evidence supports: **CYP3A4 has an intrinsic
monomeric allosteric mechanism, and oligomerisation amplifies it, sometimes by more than the
intrinsic effect, in a way that depends on membrane surface density.**

**And (b) and (d) are not independent.** The crosslinking and LRET distances for the CYP3A4 oligomer
support an interface "similar to that observed in the crystallographic dimers"
([PMID 25533469](https://doi.org/10.1074/jbc.m114.615443)) — which is **the same F′–G′ surface the
peripheral site sits on**. The peripheral site, the oligomer interface and the F/G loop OpenADMET
names as the reason co-folding fails are **one surface, not three problems**.

**2020s verdict: no single mechanism is favoured.** The honest 2026 answer is still the title of the
2008 review: *multiple binding sites, multiple conformers — or both?*
([PMID 19040328](https://doi.org/10.1517/17425250802500028)).

### 7.8 Binding and catalytic cooperativity are decoupled — the conceptual key

Three independent methods, three decades, one conclusion:

- 1994, kinetics: αNF raises V_max **without changing K_m** ([PMID 8204577](https://doi.org/10.1021/bi00187a009)).
- 2011, Nanodisc binding: the αNF/testosterone pair shows **negligible binding cooperativity**
  ([*JBC* 286:5540](https://doi.org/10.1074/jbc.M110.182055)).
- 2021, Nanodisc regiochemistry: the second midazolam binds and **"no increase in catalytic rate is
  observed"** — it changes only which C–H is oxidised, shifting 1′-OH : 4-OH to 1:2 (35% of the minor
  product) ([PMID 34015213](https://doi.org/10.1021/acs.biochem.1c00161)).

**For CYP3A4 the second ligand changes the chemistry, not the affinity.** A model or a selector
scored on affinity-like signals is aimed at the wrong observable.

⚠️ **And the location of the second midazolam is CONTESTED:** the Nanodisc work puts MDZ2 "at an
allosteric site at the membrane surface"; NMR T₁ relaxation plus docking puts two midazolams
**stacked in the active site** ([PMID 21992114](https://doi.org/10.1021/bi200924t)). The only
midazolam crystal structure (5TE8) has one MDZ per chain in the active site with the peripheral site
empty (MEASURED) and supports neither.

### 7.9 Alternate ligand conformations in CYP3A crystals — the number the challenge needs

Parsing `label_alt_id` for every non-buffer ligand copy across all 130 CYP3A4/3A5/3A7 entries:

> **Exactly one entry in the entire PDB models an alternate conformation of a CYP3A ligand: 6BDI,
> ligand DEJ, altlocs A/B at occupancy 0.41 / 0.59. 1 / 130 = 0.8%.**

Ligands refined at partial occupancy number 37 copies, almost all confined to four entries (8GK3 ×22
at 0.87–0.91, 8SO2 ×6 at 0.70–0.97, 8SO1 ×3 at 0.67–0.91, plus 9BBB 0.92, 9PLK 0.91, 8DYC 0.88).
Every other CYP3A4 ligand in the PDB is modelled at full occupancy. This is consistent with this
repo's own count over the 87-ligand validation set (`FINDING_W002`: 2 of 87 with ligand altlocs, one
of them a PEG cryoprotectant).

**Two readings, both defensible, and we should state both.**

1. *The cryoEM observation is genuinely new information.* Crystallisation selects one
   lattice-compatible conformer and refinement convention discourages depositing two; cryoEM of a
   symmetric trimer with no lattice does not. The multi-conformer density is then the first honest
   look at an ensemble that was always there.
2. *The X-ray record is evidence of absence and the cryoEM interpretation carries the burden of
   proof.* Refinement at 1.7–2.1 Å (5VCC, 6MA8, 9YK4, 8EXB, 8SO1) would have resolved altlocs if they
   were there; instead CYP3A4 ligand disorder has been absorbed into **B-factors and occupancy**
   rather than discrete conformers — and B-factors of 100–200 Å² are routine here (4K9U, 7LXL, 9PLJ,
   6DAJ).

**What this licenses for the structure track.** Three things, all concrete. (i) The partner most
worth co-folding is **a second copy of the ligand**: multi-occupancy is measured in 6 of 122 CYP3A4
entries and in the two best CYP3A5 entries, it is clash-free in every case, and the extra copy sits
against exactly the Phe cluster (F108/R212/F213/F215/F241/F304) that `FINDING_028` measured as doing
83% of the pose exclusion. (ii) **Do not build a heterotropic two-ligand experiment on structural
precedent** — zero CYP3A4 entries have two different drug-like ligands, so the reference would be
invented. (iii) The 1-in-130 alternate-conformation rate means the challenge's multi-conformer claim
has **no X-ray precedent to calibrate against**; if the reference is multi-modal while our scoring is
uni-modal, that is a plausible and testable contributor to the pure-rotation error
`FINDING_024` measured (30.1° with the protein right to 0.73 Å), and `FINDING_W002` has already shown
the submission format forecloses expressing an ensemble anyway.

---

## 8. Post-translational modifications

`CYP3A4_EVOLUTION.md` §5.4 covers the systems view — ERAD as a control point, circadian proteolysis
through gp78, and the negative finding that CYP3A4 has no credible moonlighting function. This
section adds the **residue-level map and its structural placement**, which is what a pose model can
actually use.

**Start from the baseline: UniProt annotates essentially nothing.** In
[P08684](https://rest.uniprot.org/uniprotkb/P08684.txt), `MOD_RES`, `CARBOHYD`, `LIPID` and
`CROSSLNK` features are all **zero**. The only PTM comment is "Polyubiquitinated in the presence of
AMFR and UBE2G1 and also STUB1/CHIP and UBE2D1 (in vitro)". ⚠️ That line names **UBE2G1** where the
entire primary literature uses **UBC7/UBE2G2** — an annotation slip; cite the primary papers.
Everything below therefore comes from the literature, not from the database.

### 8.1 Phosphorylation — a degron **triad**, not a pair

**ESTABLISHED.** [Wang, Liao, Hoe, Acharya, Deng, Krutchinsky & Correia 2009, *JBC*
284:5671](https://doi.org/10.1074/jbc.M806104200) (PMID 19095658): PKC targets **Thr264 and Ser420**;
liver cytosolic kinases additionally target **Ser478**; single S478A is significantly stabilised and
the **S478A/T264A/S420A triple mutant gives 4-fold stabilisation**. The commonly quoted pair
"Thr264 and Ser478" **drops Ser420**, which is a full member.

The complete in-vitro site map is [Wang *et al.* 2012, *Mol Cell Proteomics*
11:M111.010132](https://doi.org/10.1074/mcp.M111.010132) (PMID 22101235): **PKA sites** S116, S119,
S134, S259, **S478**; **PKC sites** T92, S100, T103, S116, S119, S131, S134, T136, S139, S259,
**T264**, T284, S398, **S420**. Verbatim: "whereas **PKA phosphorylates Ser-478**, **PKC phosphorylates
both Thr-264 and Ser-420**."

**The stoichiometry is the surprise.** Percent of peptide phosphorylated (native / CuOOH-inactivated):
**T284 40.1 / 52.9**, **T264 24.2 / 23.2**, S420 1.26 / 14.1, S131 0.98 / 10.5, S259 – / 10.1, and
**S478 only 0.32 / 1.29**. The paper's own conclusion: Ser478 is "the most important of these residues
for UBC7/gp78-mediated CYP3A4 ubiquitination, **despite the low extent of its PKA-mediated
phosphorylation**."

**Crucially for us: phosphorylation is a degradation signal, not an activity switch.** Verbatim:
mutation of T264, S478 and/or S420 to Ala, or S478 to Asp, "**has no appreciable effect on either its
structure or function**" — CO-reduced spectra, testosterone 6β-hydroxylase and BFC O-debenzylase all
normal. The older Oesch-Bartlomowicz thesis that P450 phosphorylation is a functional switch has
worked examples in **CYP2B1** and **CYP2E1 Ser129**, and **NO EVIDENCE FOUND** for any CYP3A4 residue.
So "phosphorylation changes CYP3A4 catalytic activity" is **unsupported at the residue level and
actively contradicted** by the mutant spectra.

**One observation worth flagging, made in the primary paper rather than by us.** **CYP3A5 is a natural
"S478D"** — it carries Asp at the CYP3A4-478 position and Lys at the CYP3A4-264 position — and is
constitutively more ubiquitinated than CYP3A4 in HepG2 **yet is not degraded faster**. Position 478 is
therefore simultaneously (i) one of the three codons under positive selection across the catarrhine
CYP3A phylogeny, (ii) the paralog-distinguishing regioselectivity residue whose S478D swap cuts
aflatoxin B1 turnover 80–90% and changes which metabolite is made (`CYP3A4_EVOLUTION.md` §1.3), and
(iii) the pivotal phosphodegron whose phosphomimetic form *is* the CYP3A5 amino acid.
**[INFERENCE]:** I have found no paper connecting all three, and the connection may be coincidence at
one convenient serine — but it is a concrete, testable hypothesis that phospho-Ser478 is a transient
chemical mimic of the CYP3A5 cavity floor.

**PhosphoSitePlus could not be retrieved** *(tooling, not absence)* — the legacy endpoints 403 and
the modern one is a JS shell. ⚠️ The correct PSP protein id for CYP3A4 is **1290776**; the commonly
quoted 6402 is **NHSL3**, a different protein. From iPTMnet instead: **Y347, T409, Y432 and S437** are
high-throughput-discovery calls with **no curated reference**, and the two PRIDE phospho-sites
(**T138** p=0.81, **S139** p=0.77, PXD012174) sit on a peptide the API itself flags as **not unique to
CYP3A4**. §8.5 argues on geometry that most of the HTP set is implausible.

### 8.2 Ubiquitination — eight mapped lysines

**ESTABLISHED in vitro.** Wang 2012 identified "eight CYP3A4 Lys residues (**Lys-115, Lys-127,
Lys-168, Lys-173, Lys-282, Lys-466, Lys-487, and Lys-492**) specifically ubiquitinated in vitro by
UBC7/gp78 and/or UbcH5a/CHIP":

- **UBC7/gp78:** K115, K168, K282, K492
- **UbcH5a/CHIP:** K127, K168, K173, K466, K487, K492 — **K168 and K492 shared**
- **Chain linkage:** gp78 largely **K48**; CHIP also **K11** and **K63** (these are ubiquitin's own
  lysines, not CYP3A4's).

⚠️ **A database artefact not to propagate:** iPTMnet and PSP list "**S116 — ubiquitination**", which
is chemically implausible for an isopeptide bond and sits immediately beside the genuine K115. It is
an off-by-one.

⚠️ **All eight are from reconstituted in-vitro E2/E3 systems.** The only in-cell corroboration is a
footnote that rat liver CYP3A Lys-168 is ubiquitinated in cultured hepatocytes. **No in-cell or
in-liver ubiquitination site map for human CYP3A4 exists — NO EVIDENCE FOUND.**

**The phosphodegron mechanism** ([Wang, Kim, Trnka, Liu, Burlingame & Correia 2015, *JBC*
290:3308](https://doi.org/10.1074/jbc.M114.611525), PMID 25451919): the phospho- and ubiquitinated
residues "reside within the cytosol-accessible surface loop and/or conformationally assembled acidic
Asp/Glu clusters… such Ser/Thr phosphorylation would generate **P450 phosphodegrons** for molecular
recognition by the E2-E3 complexes." The named clusters are **K282 flanked by Ser281 and
Glu283-Thr284-Glu285-Ser286**, near **Glu258-Ser259-Glu262-Asp263-Thr264**, and **K115** in a surface
cluster of acidic/phosphorylatable residues.

**The machinery, all ESTABLISHED:** **gp78/AMFR + UBC7/UBE2G2**; **CHIP/STUB1 + UbcH5a/UBE2D1 +
Hsc70/Hsp40**; **p97/VCP** (with Npl4, Ufd1) for retrotranslocation; the 26S proteasome. A clean
negative worth keeping: pulse-chase shows a 52% CYP3A loss over 6 h that **both proteasome inhibitors
block** while **neither NH₄Cl nor 3-MA does** ([Faouzi *et al.* 2007, *Biochemistry*
46:7793](https://doi.org/10.1021/bi700510d), PMID 17550236) — **CYP3A goes via ERAD, not via the
autophagic-lysosomal route**; that route is assigned to CYP2B1, CYP2C11 and CYP2E1. In vivo,
liver-conditional **gp78 knockout stabilises Cyp3a in a functionally active form**
([*Mol Pharmacol* 96:641](https://doi.org/10.1124/mol.119.117069), PMID 31492698).

### 8.3 Nitration — five tyrosines, and the lesion is at the redox interface

**ESTABLISHED in vitro.** [Lin, Kenaan, Zhang & Hollenberg 2012, *Chem Res Toxicol*
25:2642](https://doi.org/10.1021/tx3002753) (PMID 23016756): peroxynitrite nitrates
**Tyr99, Tyr307, Tyr347, Tyr430 and Tyr432**. "**Tyr99, Tyr347, and Tyr430 are on the proximal side**
of CYP3A4 and are in close contact with three acidic residues in the FMN domain of CPR."
The functional lesion is **electron transfer from POR, not the heme**: the heme absolute spectrum is
unchanged, activity loss is much greater with POR than with t-BuOOH, and **both Y430F and Y430V lose
CPR-dependent activity** ("both the aromatic and the hydroxyl groups of Tyr are required"). Contrast
with CYP2B6/2E1, where peroxynitrite *does* modify the heme — CYP3A4 is the exception. ⚠️ Purified
enzyme and chemical peroxynitrite; **not demonstrated in cells or liver**.

**NO-dependent degradation: a clean species split.** Human CYP3A4 does **not** show it — NOS
inhibitors fail to block cytokine-driven CYP3A4 downregulation in human hepatocytes
([PMID 18206661](https://pubmed.ncbi.nlm.nih.gov/18206661/)) — while rat CYP3A1 does. The
residue-resolved nitration→degradation work is on **CYP2B6 (Y190, Y317, Y380)** and does not transfer.

**A better-evidenced NO mechanism, and one with a direct bearing on co-folding.**
[Das Sinha, Islam, Biswas & Stuehr 2025, *JBC* 301:110772](https://doi.org/10.1016/j.jbc.2025.110772)
(PMID 41022321): CYP3A4 and CYP2D6 exist as a **~60:40 heme-free : heme-bound mixture** in CHO and
HepG2 cells; **1–10 nM NO drives heme into the apo-pool** giving 2–3× more heme-replete enzyme (via a
GAPDH–heme complex and Hsp90), while 25–100 nM NO causes heme loss. **Heme occupancy is not 1.0** —
worth remembering wherever a pipeline hard-codes HEM into every input.

**S-nitrosylation of Cys442: NO EVIDENCE FOUND.** UniProt has no MOD_RES at 442, only the heme-Fe
axial annotation, and the established NO–P450 chemistry is Fe–NO, not thiolate S-nitrosylation.
⚠️ **Do not confuse with PXR**, which *is* S-nitrosylated — that is the transcription factor.

### 8.4 Covalent adducts — three named residues, and they are not where you would guess

| compound | mechanism | residue | source |
|---|---|---|---|
| **raloxifene** | apoprotein alkylation by the diquinone methide, +471 Da | **Cys239** | [PMID 17497897](https://pubmed.ncbi.nlm.nih.gov/17497897/) |
| **bergamottin** | apoprotein adduct by oxygenated 6′,7′-dihydroxybergamottin, +388 Da, on peptide 272–282 | **Gln273** | [*DMD* 40:998, PMID 22344702](https://doi.org/10.1124/dmd.111.043802) |
| **ritonavir** | covalent bond to apoprotein, 0.93 ± 0.04 mol/mol, on peptide 255–268 | **Lys257** (the only Lys in that span) | [*Mol Pharmacol* 86:665, PMID 25274602](https://doi.org/10.1124/mol.114.094862) |
| 17α-ethynylestradiol | ~50% heme destruction **and** apoprotein binding | **none identified** | PMID 11907170 |
| mibefradil | **heme destruction only**, "lack of stable heme and/or apoprotein adducts" | **none** | PMID 21447734 |
| delavirdine | irreversible binding to a ~50 kDa protein | **none identified** | PMID 9765359 |

**Two corrections worth carrying.** (i) **Lys257 is the ritonavir adduct; Cys239 is the raloxifene
adduct** — two separate findings from the same group, routinely conflated. (ii) **Cys239 is a
recurrent electrophile sink, not a raloxifene quirk**: it is also hit by chelerythrine
([PMID 38330258](https://pubmed.ncbi.nlm.nih.gov/38330258/)) and by pyrene-iodoacetamide, and an MD
study finds the diquinone methide reaches it **during unbinding**
([*JCIM* 62:6172, PMID 36457253](https://pubmed.ncbi.nlm.nih.gov/36457253/)).

⚠️ **CONTESTED — the ritonavir mechanism.** [*DMD* 41:1813, PMID
23886699](https://doi.org/10.1124/dmd.113.052001) reports ~50% heme loss and a heme–protein adduct;
the 2014 paper finds radioactive ritonavir co-eluting with the polypeptide "**and not with heme**".
Sevrioukova's review is titled "The Mechanism-Based Inactivation of CYP3A4 by Ritonavir: What
Mechanism?" and concludes it is unsettled ([PMID 36077262](https://doi.org/10.3390/ijms23179927)).

⚠️ **Two qualifiers that should temper any adduct-based modelling.** Adduction is not inactivation —
pyrene-iodoacetamide "alkylated residue **Cys239** exclusively" yet **did not inhibit** CYP3A4
([PMID 32193357](https://pubmed.ncbi.nlm.nih.gov/32193357/)). And adducts are promiscuous in real
tissue — raloxifene adducts were found on **78 proteins** in human liver microsomes, including CYPs
not subject to time-dependent inhibition, so "adducts can be **benign**"
([PMID 39442081](https://pubmed.ncbi.nlm.nih.gov/39442081/)).

### 8.5 Acetylation, methylation, SUMOylation, glycosylation

- **Acetylation: HTP-MS only** — K70, K127, K413, no PMID, no enzyme. ⚠️ **K127 is also a
  CHIP-ubiquitination site**, so acetylation/ubiquitination competition there is plausible but
  **untested [INFERENCE]**.
- **Methylation: HTP-MS only** — R158, R161, no PMID.
- **SUMOylation: NO EVIDENCE FOUND, explicitly.** Every "CYP3A4 SUMOylation" hit is about **PXR**
  being SUMOylated. Do not cite PXR SUMOylation as CYP3A4 SUMOylation.
- **Glycosylation: ESTABLISHED ABSENT.** No CARBOHYD feature, and mechanistically impossible — the
  N-terminal helix (TRANSMEM 2–22) is the sole anchor and the ~480-residue catalytic domain faces the
  **cytosol**, where OST is not.

### 8.6 Where every PTM site sits — computed from 1TQN, not guessed

Distances and burial were computed for this document from `data/reference/rcsb/1TQN.cif` (ligand-free
CYP3A4, chain A, 468 observed residues). Burial is a **percentile of the 10 Å neighbour-count
distribution over all 468 residues** (min 37, p25 113, median 158, p75 191, max 237), because raw
counts mean nothing on their own. **MEASURED.**

| residue | modification | d(Fe) Å | d(heme) Å | burial %ile | region |
|---|---|---|---|---|---|
| **K115** | Ub (gp78) | 16.4 | 10.4 | **13%** | B′–C loop (SRS1) |
| **S119** | p (PKA/PKC) | **8.1** | **3.6** | 63% | **active site**, in this repo's `CYP3A4_POCKET["polar"]` |
| **K127** | Ub (CHIP) + acetyl | 16.3 | 7.7 | 32% | C-helix, **proximal face / POR interface** |
| **K168** | Ub (gp78 + CHIP) | 31.2 | 25.2 | **3% — joint most exposed** | D–E loop |
| **C239** | raloxifene / chelerythrine adduct | 19.5 | 16.4 | 53% | **G′/G helix, SRS3 — inside the F/G span** |
| **K257** | ritonavir adduct | 23.1 | 18.2 | 68% | **H-helix — inside the F/G span** |
| **S259** | p (PKA/PKC) | 29.7 | 23.7 | 16% | H–I loop, in the E258-S259-E262-D263-T264 acidic cluster |
| **T264** | **p (PKC) DEGRON** | 31.6 | 25.0 | **3% — joint most exposed** | H–I loop, cytosolic surface |
| **Q273** | bergamottin adduct | 19.3 | 13.0 | 67% | I-helix N-terminus |
| **K282 / T284** | Ub (gp78) / **p (PKC), 40–53% stoichiometry** | — | — | **UNMODELLED** | **the only internal gap in 1TQN, residues 282–285** |
| **Y307** | nitration | **8.8** | **5.7** | 77% | I-helix, SRS4 — active site |
| **Y99 / Y347 / Y430** | nitration | 16.5 / 16.9 / 16.3 | 9.1 / 13.8 / 9.2 | 32 / 37 / 44% | **proximal face, POR contacts** |
| **S420** | **p (PKC) DEGRON** | 23.0 | 19.0 | 30% | **meander, immediately after the PERF motif — proximal face** |
| **S478** | **p (PKA) MAJOR DEGRON** | 20.0 | 16.6 | 27% | **β4 hairpin, adjacent to SRS6 (476–486)** |
| K466 / K487 / K492 | Ub (CHIP / both) | 28.6 / 22.8 / 28.9 | 21.9 / 18.1 / 22.3 | 22 / 7 / 27% | β3/β4 and C-terminal tail |
| *C442* | *heme axial thiolate (reference)* | *2.3* | — | — | *heme signature* |

**Five things the geometry says.**

1. **The PTM landscape and the ligand-binding landscape are almost disjoint.** Of 38 modified
   residues, **exactly one (Ser119) is in this repo's pocket definition**, and only three (C239,
   K257, S259) fall inside the F/G span 202–260. Every validated degron and **all eight ubiquitinated
   lysines** sit at burial percentiles **3–32%**, i.e. on the solvent-exposed surface, exactly as the
   E2/E3 phosphodegron model requires. **PTM state should not perturb pose prediction in the distal
   cavity.**
2. **T264 and K168 tie as the two most solvent-exposed residues in the whole protein (3rd
   percentile)** — and they are the two best-supported modification sites. A satisfying internal
   consistency check.
3. ⚡ **The highest-stoichiometry PTM site in CYP3A4 lies in the only disordered segment of 1TQN.**
   The crystal's single internal gap is **residues 282–285**, precisely the K282 (gp78-ubiquitinated)
   / **T284 (40–53% phosphorylated)** cluster that Correia's group independently flagged as "a
   'disordered' region not featured in the crystal structure". Two labs, two methods, the same four
   residues. **Any co-folding model is guessing at 282–285, and that is where the dominant covalent
   modification lives.**
4. **A coherent proximal-face cluster.** Y99, Y347 and Y430 (all three named as POR-contacting),
   Y432, phospho-S420 immediately after the **PERF motif** (P416-E417-R418-F419), and the C-helix
   cluster K127/S131/S134/T138/S139. **[INFERENCE]:** phosphorylating S420 drops a negative charge
   directly beside the conserved **E-R-R triad** — mechanistically coherent with Y430F abolishing
   POR-dependent turnover, though not a published claim.
5. **The adduct residues line the egress path, not the catalytic cavity.** C239 (16.4 Å from heme),
   K257 (18.2 Å) and Q273 (13.0 Å) are all too far from the heme to be oxidised in situ but buried
   enough to be unreachable from bulk solvent. This matches the MD finding that raloxifene's
   metabolite hits Cys239 **during unbinding**. **A useful structural prior: the substrate egress
   channel, not the active site, is the electrophile sink.**

**And a sanity check on the HTP set.** The unreferenced HTP phospho-calls are precisely the
implausible ones: **S119 (3.6 Å from heme), S437 (3.4 Å, inside the heme signature FGSGPRNCIG),
Y307 (5.7 Å, I-helix SRS4)** and **S398 (6.4 Å)** all sit at 50–80th-percentile burial **inside the
catalytic cavity**, where no kinase can reach them in the folded holo-enzyme. The validated degron
triad — T264, S420, S478 — is untouched by that criticism; all three are surface.

**What this licenses for the structure track.** (i) **Nothing about the PTMs argues for modelling
them** — they are surface, sub-stoichiometric for the functional ones, and demonstrably do not change
structure or activity. That is a useful negative: do not spend on modified-residue inputs.
(ii) **Residues 282–285 are unmodelled in the reference crystal and are the site of the dominant
modification.** Any per-residue attribution of pose error (as in `FINDING_028`) should exclude that
span rather than silently treat it as agreement. (iii) **Cys239 and Lys257 are both inside the F/G
span and both are documented electrophile sinks reached during ligand egress** — which is independent,
chemistry-based evidence that the F/G region is the part of the protein a ligand actually touches on
its way in and out, and therefore the part a scorer must get right. (iv) **Heme occupancy in cells is
~40%, not 100%** — a caution for any pipeline that assumes the holo form.

---

## 9. What this whole map licenses, and what it forecloses

**The three partners worth co-folding, in order — with rank 1 now spent.**

1. ⚠️ **A second copy of the query ligand — RUN AND REFUTED (2026-09).** It was the only "partner"
   with CYP3A4 structural precedent (6 of 122 entries, clash-free in every case, minimum inter-copy
   distance 2.90 Å), the only one that touches the F/G roof, and one extra ligand entity in a YAML.
   **48 scored OpenProtein jobs later it is a negative:** the second copy landed in the **active site on 12
   of 15 ligands** instead of the peripheral groove, and the first copy got **worse by −0.0664
   LDDT-PLI with rotation error moving +4.45° the wrong way**. What did hold was the geometry — median
   **8.81 Å from the iron**, **2V0M's second ketoconazole reproduced unprompted**, and **45 of 60
   second copies contacting the §7.2 groove residues** — but from **inside** the cavity, displacing the
   first ligand. **Right about where, wrong about whether it helps.** See the box in §6.9 and `FINDING_031`.
2. **Membrane lipid (POPC, or 1:1 POPC:POPG) — now the only untested partner, on a LOWERED prior.**
   The F′–G′ face is bilayer-embedded, 8GK3 resolves a steroid sulfate sitting on it, lipid identity
   changes measured fragment binding, and the allosteric midazolam site is located in that region in
   nanodiscs. Also just ligand entities. ⚠️ But the one thing now measured is that **adding mass on
   the distal side makes the pose worse**, so the surviving hypothesis is the narrower one: that lipid
   restrains residues 211–216 **from outside** the cavity, where a second ligand demonstrably acts
   from inside. Run it as a falsification test, not as an expected gain.
3. **Nothing else.** CYB5A has the best biochemistry of any protein partner (K_d 2.5–61 nM, a ternary
   complex with POR, substrate-dependent affinity) and POR is obligatory — and both bind the
   **proximal face, 9.7 Å from the iron and 10.7 Å from the nearest modelled F/G residue, across the
   porphyrin**, on a target where protein accuracy is uncorrelated with ligand error (ρ = +0.03).
   `FINDING_024` already spent that argument; §6 only strengthens it.

**The paralogs worth docking against.** **CYP3A5** first — co-expressed in the same liver at 156 TPM,
**zonated in the opposite direction**, six structures, and two ligands (clotrimazole in 8SPD/8SG5,
azamulin in 6OOA/7SV2) crystallised in **both** paralogs, which is a ready-made leak-free pair.
**CYP3A7** second — two structures, ~93% identity, a pocket with nearly the same residues that
**cannot deform**, which makes the 3A4/3A7 pair the cleanest test of whether a scorer is sensitive to
lid *plasticity* rather than to sequence identity. **CYP3A43** has no structure and no known
substrate, so it is a target for a different question entirely.

**The compound space.** The `ranked_compounds` list in the JSON spans what CYP3A4 actually handles in
the liver rather than what happens to have been crystallised. The reason is measured: the local
CYP3A4 structural proxy set is **72 of 87 Type II heme-coordinating ligands and only 15 Type I**
(from `data/processed/poses_scored_val87b.csv`), and Reactome models CYP3A4 as a purely xenobiotic
enzyme while UniProt documents cholesterol at five positions, estradiol at four, anandamide
epoxidation and vitamin D 23/24-hydroxylation with Rhea IDs and experimental evidence. **Both of our
usual sampling frames are biased away from the endogenous, non-coordinating, lipophilic substrate
space — which is exactly the space the evolution document argues the pocket was built for.**

**What the map forecloses, and this is the most useful negative in it.** Sections 2, 3 and 5 between
them cover every known route into CYP3A4's enormous inter-individual variability — splice variants,
alternative polyadenylation, six validated miRNAs, a contested lncRNA, promoter methylation, histone
marks, enhancer looping, zonation, development, induction, inflammation, pregnancy, fasting and
disease. **Not one of them produces a structurally different CYP3A4 protein.** \*22 makes a
transcript with no heme-binding signature; the miRNAs and APA change abundance and translation rate;
A-to-I editing is confined to a 3′UTR Alu at 0.42% and cannot recode; the PTMs are surface,
sub-stoichiometric where functional, and demonstrably do not change structure or activity. **Every
route into variation is a route into *how much* enzyme there is, not *which* enzyme.** So there is one
protein to predict, and the challenge's multi-conformer observation — for which the X-ray record gives
**1 alternate-conformation entry in 130** — is about ligand placement, not about protein variants.

---

## 10. Reproducibility notes and open gaps

### Search traps that will bite anyone re-running this

- **Europe PMC phrase queries URL-encoded with `%20` silently return `hitCount: 0`** rather than an
  error. Use `+`. This was hit independently five times in this session and is indistinguishable from
  a real negative.
- **The entire CYP3A4\*22 splicing literature is indexed under "intron 6", not "alternative
  splicing"**, and the **4C papers are only findable via "chromatin conformation capture"**. The
  obvious keyword returns nothing in both cases.
- **ARCHS4's co-expression is a POST endpoint** (`POST https://maayanlab.cloud/matrixapi/coltop`,
  `{"id":"CYP3A4","count":101}`). Fetching the gene page returns a JavaScript shell with no data.
- **PhosphoSitePlus's protein id for CYP3A4 is 1290776.** The commonly quoted 6402 is **NHSL3**, a
  different protein.
- **The Reactome endpoint is `/data/mapping/UniProt/P08684/pathways`.** The `/pathways/low/entity/…`
  form 404s — it wants a Reactome stId, not a UniProt accession.
- **GTEx's eQTL association endpoints have incomplete coverage.** They return zero Liver rows for
  CYP3A5 and never return rs776746 in any tissue, and `independentSqtl` returns zero rows for
  everything. **Do not report a GTEx eQTL null for CYP3A4.** The expression endpoints are fine.
- **COXPRESdb is dead** (every API pattern 404s), and **HumanBase's liver network** returns CYP3A4
  edges all at the 0.1065 mincut with olfactory receptors and POTE paralogs in the top 40 — noise.
- **HepG2 is not a CYP3A4 chromatin model**: no promoter H3K27ac, no H3K4me3 anywhere, no promoter
  DHS. Much of the reporter literature is HepG2.

### Things nobody has measured

- **No transcription-shutoff half-life for CYP3A4 mRNA in primary human hepatocytes.** Existing
  estimates disagree by 4–8× and have never been reconciled.
- **No test of whether the \*22 intron-6-retention transcript is an NMD substrate.**
- **No nuclear run-on and no Pol II ChIP at CYP3A4**, so "rifampicin increases the transcription rate"
  has never been measured as a rate — everything is steady-state mRNA.
- **No Hi-C, TAD or insulator analysis of the CYP3A locus.** Exploitable.
- **No test of parent-of-origin expression at CYP3A**, and no CYP3A gene appears in the Geneimprint
  catalogue.
- **No snATAC-seq or zonation-resolved CYP3A4 accessibility.**
- **No classical-twin heritability estimate for *basal* (as opposed to induced) CYP3A4 activity.**
- **No demonstration that CYP3A7 and CYP3A4 compete for the shared DRR**, and nobody has deleted it at
  the relevant developmental stage. **The mechanism of CYP3A7 shut-off is the field's biggest hole**;
  the only candidate is NFI, unfollowed-up since 2009.
- **No human study perturbing the microbiome and measuring CYP3A4 activity.**
- **No MLL3/MLL4 ChIP at CYP3A4** — the ASCOM link is inferred from NCOA6 alone.
- **No in-cell or in-liver ubiquitination site map** for human CYP3A4; all eight lysines are from
  reconstituted in-vitro E2/E3 systems.
- **No named UGT or SULT for 6β-hydroxytestosterone or for 4β-hydroxycholesterol.**
- **No quantitative test of the CYP3A4/P-gp substrate-overlap claim against a chance baseline.** It
  has been asserted since 1995 and never measured.
- **No experimental structure of CYP3A4 with any protein partner.** All 131 PDB entries are the
  soluble Δ-N catalytic domain alone.

### The machine-readable graph

[`data/processed/cyp3a4_biology_graph.json`](../../data/processed/cyp3a4_biology_graph.json) holds
the same content as a node/edge graph over eight layers (gene, regulation, transcript, protein,
complex, pathway, celltype, compound), with `meta.ranked_partners` and `meta.ranked_compounds` as the
two lists meant to drive the next round of co-folding.

**Every node and every edge carries an explicit `status` field** — `ESTABLISHED`, `CONTESTED`,
`INFERRED` or `NO_EVIDENCE_FOUND` — so the grading in this document is queryable rather than only
readable. Current distribution over **265 graded claims**: nodes 108 established / 11 contested / 1
no-evidence-found; edges 129 established / 11 contested / 5 inferred. `meta.status_vocabulary` and
`meta.status_note` carry the definitions. `meta.ranked_partners` entries carry a `verdict` pill from
`{TEST_NEXT, SPENT_NEGATIVE, DO_NOT_SPEND, SPECULATIVE}` and a `prior` pill from
`{RAISED, UNCHANGED, LOWERED, REFUTED, NOT_APPLICABLE}`, with the prose in `rationale`.

The reason the `status` field exists is in the §6.9 box: this map has now caused an experiment, and
that experiment was a negative. A CONTESTED claim needs to be visible where someone is deciding
whether to spend compute on it, not only here in the prose.

Node `detail` strings are tooltip-length (≤400 characters, asserted). The build asserts that every
edge endpoint resolves to a node id, that **no node is left without any edge**, that every node has at
least one citation, that every status and verdict is in its vocabulary, and that the file re-parses
and re-validates identically after writing.
