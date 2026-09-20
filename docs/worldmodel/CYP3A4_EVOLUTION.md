# CYP3A4 — a causal decomposition

**Written 2026-09-20 as part of ask 2 of the shared-embedding world model (`docs/worldmodel/README.md`).**

This is not a summary of what CYP3A4 does. It is an argument for *why it is the way it is*, built
from phylogeny, selection statistics, endogenous chemistry, structure, and regulatory context, and
ending in a prediction about what its pocket should accept that no one has assayed.

**How to read it.** Everything stated as fact carries a DOI, PMID, PDB ID or database accession.
Everything I am reasoning to rather than reporting is tagged **[INFERENCE]**. Where the literature
genuinely disagrees I say so and give both sides rather than picking one — the disagreements are
usually the informative part. Where a popular story turns out to be thin, I say that too; §2 is
mostly that.

---

## 0. The object, precisely identified

Four genes, one ~220 kb tandem array on **chromosome 7q22.1** (GRCh38, Ensembl):

| gene | Ensembl | coords (chr7) | strand | UniProt | length |
|---|---|---|---|---|---|
| **CYP3A5** | ENSG00000106258 | 99,647,442–99,680,043 | − | **P20815** | 502 aa |
| **CYP3A7** | ENSG00000160870 | 99,704,105–99,735,196 | − | **P24462** | 503 aa |
| **CYP3A4** | ENSG00000160868 | 99,756,954–99,784,327 | − | **P08684** | 503 aa |
| **CYP3A43** | ENSG00000021461 | 99,827,922–99,867,801 | **+** | **Q9HB55** | 503 aa |

Physical order along the chromosome is **5 – 7 – 4 – 43**, which is exactly the order in the
reconstructed ancestral catarrhine locus (Qiu et al. 2008, Fig. 4). The array is interleaved with
"detritus exons" — fragments of dead CYP3A copies — and with the pseudogenes CYP3AP1/CYP3AP2. Three
of the four genes are on the minus strand; **CYP3A43 is inverted**, the signature of the event that
accompanied its birth.

Two structural facts to hold onto from the start:

- CYP3A4 and CYP3A7 are **~93% identical** in protein sequence, and CYP3A4/CYP3A5 ~84%. The
  functional differences below are therefore produced by a few dozen residues, most of them in and
  above the pocket.
- The proximal heme thiolate is **Cys442** in all four. That is the one thing none of them can change.

---

## 1. Phylogeny and selection: what happened, and what it was for

### 1.1 Deep time — CYP3A is a *relocated* gene, and eutherians threw one away

The definitive comparative-genomic work is **Qiu, Taudien, Herlyn, Schmitz, Zhou, Chen, Roberto,
Rocchi, Platzer & Wojnowski (2008), *Pharmacogenetics and Genomics* 18:53–66,
doi:10.1097/fpc.0b013e3282f313f8 (PMID 18216722)** — 16 vertebrate genomes, ~450 My, 107 CYP3
protein sequences from 31 species, Bayesian + ML phylogeny cross-checked with 10 retroposon
insertions and 44 indels as rare genomic changes. Its results are the backbone of this section.

- The CYP3 **family** (CYP3A/3B/3C/3D) diverged from other P450s ≥800 Mya. All four subfamilies
  coexist only in **ray-finned fish**; the CYP3B/3C split traces to the fish-specific whole-genome
  duplication (Yan & Cai 2010, *PLoS ONE* 5:e14276, doi:10.1371/journal.pone.0014276).
- **Amniote ancestors carried two CYP3A genes**, ancestral to CYP3A80 and CYP3A37. Both survive in
  birds, squamates, monotremes and marsupials.
- **At the origin of Eutheria, the CYP3A80 lineage was lost outright**, and the surviving CYP3A37
  ortholog was **translocated from one genomic neighbourhood (CYPHR1) to another (CYPHR2)**. Every
  CYP3A gene in every placental mammal — mouse Cyp3a11/13/16/25/41/44/57/59, dog CYP3A12/3A26, pig
  CYP3A29, and all four human genes — descends from that single relocated copy.

**[INFERENCE]** A translocation plus the loss of its paralog is not a neutral event: it means the
surviving copy landed under a new set of *cis* elements. The eutherian CYP3A story therefore starts
with a **regulatory** change, not a catalytic one. This matters for §5: the XREM/CLEM enhancer
architecture that makes CYP3A4 the most inducible drug-metabolizing enzyme in the human genome sits
in a neighbourhood the gene *moved into*. The promiscuity we care about may be downstream of where
the gene ended up as much as of what it binds.

### 1.2 The primate expansion — birth, death, conversion, all in ~65 My

Qiu et al. reconstructed the locus in human, chimp, orangutan, rhesus, olive baboon, marmoset and
galago (BAC libraries for the last two; GenBank EF600755–EF600758, EF589790–EF589808):

- **Galago (Strepsirrhini):** ≥3 tandem CYP3A genes, one **inverted**. Basal, and divergent enough
  that none of them is a one-to-one ortholog of anything human.
- **Marmoset (New World monkey):** three genes, **CYP3A5 / CYP3A90 / CYP3A21**. It expanded the
  *CYP3A5* branch, not the CYP3A4 branch. **No human CYP3A gene has a one-to-one marmoset ortholog**,
  and no CYP3A43 ortholog exists there at all. Qiu et al. use this to argue explicitly that
  **marmoset is a poor model of human CYP3A pharmacology**, despite its popularity.
- **Rhesus / baboon (Old World monkeys):** four genes — CYP3A5, 3A7, 3A4, 3A43 — in human order,
  plus a ~40 kb low-copy-repeat insertion between CYP3A7 and CYP3A4 that humans lack.
- **Hominidae:** a **fifth** gene, **CYP3A67**, produced by duplication of CYP3A7 early in the great
  apes. Chimpanzee and orangutan have it; **humans lost it by unequal crossover** with breakpoints
  downstream of CYP3A7 and CYP3A67 (not, as earlier work claimed, by recombination with CYP3AP1).
  Qiu et al. PCR-screened 99 Central European and 38 Bantu samples and found no polymorphic
  retention — the loss is fixed.
- **CYP3A43 sits at the base of the anthropoid CYP3A phylogeny** — it is the sister of all other
  anthropoid CYP3A genes, so it predates the catarrhine/platyrrhine split (>40 Mya) and was
  *lost* in marmoset. It is now pseudogenising **independently in three lineages**: human and
  chimpanzee produce mostly aberrant transcripts, and a frameshifting exon-9 deletion segregates in
  rhesus.
- **Gene conversion is rampant.** Nine converted gene pairs across rat, dog, pig, medaka and
  primates. The most consequential: in the common hominid ancestor, **exon 6 and part of intron 6 of
  CYP3A7 were overwritten by the corresponding piece of CYP3A4** — traceable by a shared 15 bp
  intron-6 indel present in human/chimp/orangutan CYP3A4 and CYP3A7 but absent from rhesus, baboon,
  CYP3A67 and all CYP3A5. **The paralogs are not independent lineages; they trade sequence.**

**[INFERENCE]** The gene conversion is the single most under-appreciated fact in this section. It
means that (a) phylogenies built from CYP3A coding sequence are locally wrong wherever conversion
has acted, and (b) CYP3A4 and CYP3A7 are more similar in the converted block than their divergence
time implies. Anyone using CYP3A7 as an "independent" check on a CYP3A4 model is partly checking
CYP3A4 against itself.

### 1.3 Where positive selection actually falls — the residues

PAML site and branch models on 18 catarrhine CYP3A coding sequences (Qiu et al. 2008, Table 1):

| test | result |
|---|---|
| M1a vs **M2a** | 2Δℓ = 11.44, *P* < 0.01. ω₂ = **2.63** on ~7% of sites |
| M7 vs **M8** | 2Δℓ = 15.78, *P* < 0.001. ω = **2.29** on ~11% of sites |
| M0 vs free-ratio | 2Δℓ = 53.06 vs cv 46.19, *P* < 0.05 → ω is lineage-specific |

**Sites under selection across the whole phylogeny (BEB, P(ω>1) ≥ 0.95 under M8): codons 437, 478
and 479.** Under M2a the same three sites at P = 0.87–0.92. Note the corollary in the same table:
**89% of sites are under slight-to-strong purifying selection**. This is a conserved protein with
three hot residues, not a fast-evolving one.

Two branches carry ω = **∞** (nonsynonymous changes with *zero* synonymous changes):

- **Branch Ho7 — CYP3A7 in the hominoid stem, 15.5 nonsynonymous / 0.0 synonymous.** Eighteen
  candidate sites: **H28R, S29T, V50A, F57Y, F74I, G77C, R78Q, V81M, L108F, S116N, S215P, I220V,
  S286T, V296M, E333K, A337T, N437S, R484L.** Eight of the eighteen lie in the first 100 residues —
  the N-terminal membrane anchor and redox-partner interface, not the pocket.
- **Branch Hs4 — human CYP3A4 after the chimpanzee split, 6.0 nonsynonymous / 0.0 synonymous.** Six
  sites: **R54H, R78Q, I129L, I224T, R478S, T489V.**

**R78Q arose independently on both branches** — a parallel substitution whose function is still
unknown. Position 108 flips L→F on the CYP3A7 branch, restoring the ancestral phenylalanine.

Why those residues matter functionally, per Qiu et al.'s own review of prior mutagenesis:

- **F108** is one of the phenylalanines forming the ordered hydrophobic "roof" above the CYP3A4
  active site. Its identity changes activity toward midazolam, aflatoxin B1, testosterone,
  progesterone and lapachol.
- **S478D** (the CYP3A5 residue) leaves testosterone 6β-hydroxylation untouched but cuts aflatoxin
  B1 turnover by 80–90% **and changes which aflatoxin metabolite is made**.
- **L479T** (CYP3A5-like) does the same to aflatoxin and halves testosterone hydroxylase activity;
  **L479F** (CYP3A7-like) changes the 7-hexoxycoumarin product profile.
- **P215** sits in the F–F′ connector and helps orient progesterone, but its mutagenesis changed
  neither activity nor cooperativity — the one site of the four where the functional story fails.

**So three of four selected residues that have been mutagenised alter regioselectivity, not
affinity.** Selection on primate CYP3A has been acting on **where the substrate gets oxidised**,
not on whether it binds.

**[INFERENCE — and this is the load-bearing inference of §1].** Positions 437, 478 and 479 are in
the **β4 / meander / C-terminal loop region**, directly beneath and behind the heme, and 479 is two
residues from the 480–482 stretch that the CYP3A7•DHEA-S structure shows contacting the steroid
(§4.4). These are **product-orientation residues**, not gatekeeping residues. Combine that with §4's
finding that the F/G roof is where the big conformational motion lives, and the picture is: the
**entrance evolved by changing its dynamics; the floor evolved by changing its chemistry.** Those
are two separately tunable knobs, and primate evolution turned the second one.

### 1.4 Does CYP3A stand out among P450s? Partly

**Ryan, Bradley and colleagues (2020), *Xenobiotica* 50:1406–1412,
doi:10.1080/00498254.2020.1785580 (PMID 32558606)** classified the 15 CYP1–3 subfamilies as
xenobiotic-metabolizing (XM) or endogenous-metabolizing (EM) across primates and tested two
predictions.

- Prediction 1 — XM subfamilies should vary more in gene copy number — **failed**. EM and XM
  subfamilies varied about equally.
- Prediction 2 — XM subfamilies should show diversifying selection, concentrated on substrate
  recognition sites (SRSs) — **partly held**. Exactly four subfamilies showed diversifying selection:
  **CYP2C, CYP2D, CYP2E and CYP3A**; no EM subfamily did. Three of those four (2C, 2D, **3A**)
  showed a significant link between SRSs and the selected codons.

Independently, **Roca-Umbert et al. (2019), *BMC Evolutionary Biology* 19:39,
doi:10.1186/s12862-019-1366-7 (PMID 30704392)** applied Hierarchical Boosting to 91 taste-receptor
and P450 genes in CEU/CHB/YRI. Their headline negative is worth as much as the positive: the famous
population differences at **TAS2R38 and TAS2R16 are consistent with drift, not selection**, and
CYP1/CYP2 genes showed no sweep signature. Only three P450s did: **CYP3A4 and CYP3A43 in Europeans**,
and CYP27A1 in East Asians.

So the evidence is: CYP3A is genuinely one of the few xenobiotic P450 subfamilies under detectable
diversifying selection in primates, **and** the CYP3A locus is one of the few P450 loci carrying a
recent human sweep. That is a real signal. What it was selected *for* is §2 and §5.

### 1.5 The fetal/adult switch — why CYP3A7 exists

The facts:

- CYP3A7 is the dominant CYP3A in **human fetal liver**; UniProt P24462's only tissue-specificity
  annotation is "expressed in fetal liver (at protein level)". Activity peaks around birth, falls
  sharply in the first week, and persists at low level through the first year, while CYP3A4 rises.
- CYP3A7's distinctive endogenous reaction is **16α-hydroxylation of dehydroepiandrosterone
  3-sulfate (DHEA-S)** — UniProt P24462 annotates it, and CYP3A4 does not have it. The product
  16α-OH-DHEA-S returns to the placenta, is desulfated, and is aromatised by CYP19A1 to **estriol**,
  the defining estrogen of human pregnancy.
- CYP3A7 also 4-hydroxylates and 18-hydroxylates **all-trans retinoic acid** (P24462), the ligand of
  the embryonic retinoic acid receptors — i.e. it is positioned to terminate a morphogen gradient.
- **Qiu et al.'s branch Ho7 (ω = ∞) is exactly the branch on which hepatic CYP3A7 expression became
  restricted to the fetal period.** Rhesus, olive baboon and hamadryas baboon express CYP3A7 in adult
  liver; humans, chimpanzees and orangutans do not. The 15 nonsynonymous / 0 synonymous burst and
  the expression restriction happened on the same branch, alongside the exon-6 conversion from
  CYP3A4.

**[INFERENCE] The switch is not "a fetal enzyme and an adult enzyme."** It is *one* CYP3A gene that
was expressed throughout life in the Old-World-monkey ancestor, then split its job in two at the
hominoid stem: CYP3A7 kept the **steroid-sulfate and retinoid** work and was confined to the window
where estriol and RA gradients matter, while CYP3A4 took the adult compartment. The likely driver is
that the two jobs became **incompatible**: an enzyme tuned to clear retinoic acid efficiently is a
liability in an adult that needs retinoid signalling, and an enzyme tuned for broad adult xenobiotic
clearance is a teratogenic hazard in a fetus. Temporal separation is the cheapest way out of a
pleiotropy trap, and the 18 substitutions are the cost of specialising the fetal copy afterwards.

The **ancestral CYP3A function is therefore visible in CYP3A7, and it is steroid chemistry** —
specifically the oxidation of a **sulfated, anionic steroid at an unactivated ring carbon** (C16 of
the D-ring). Note what that requires: a pocket that can hold a rigid four-ring hydrophobe *and*
place an anionic sulfate somewhere tolerable, and that can bring an sp³ C–H 5–6 Å from the iron.
Hold that thought for §6.

### 1.6 Non-primate CYP3A orthologs — a warning about model species

Because every eutherian CYP3A descends from one translocated gene and most extant CYP3A genes are
products of **recent, lineage-specific** duplication (Qiu et al. 2008), CYP3A gene names across
species are near-meaningless as functional labels. Mouse has ~8 Cyp3a genes; dog has CYP3A12 and
CYP3A26; marmoset's three are none of ours. **There is no mouse CYP3A4.** Humanised-liver mice and
CYP3A7-transgenic mice exist precisely because of this.

**[INFERENCE]** For the OpenADMET structure challenge this has one practical consequence: the
P450-superfamily reference set in `data/processed/p450_universe/` contains many CYP3A entries whose
relationship to CYP3A4 is *paralogy within a recent radiation*, not orthology. Sequence identity
will overstate functional transferability, and structural transfer should beat sequence transfer —
consistent with FINDING 022's observation that co-folding has essentially solved the rest of the
P450 family while CYP3A4 remains the exception.

---

## 2. The diet / xenobiotic hypothesis: how load-bearing is it?

### 2.1 The claim

**Gonzalez & Nebert (1990), *Trends in Genetics* 6:182–186, doi:10.1016/0168-9525(90)90174-5
(PMID 2196721)** — "Evolution of the P450 gene superfamily: animal–plant 'warfare', molecular drive
and human genetic differences in drug oxidation." The claim: P450 diversity is the product of
continuous coevolution between plants making phytoalexins and animals answering with new
detoxifying enzymes, and human pharmacogenetic polymorphism is the residue of that arms race. The
frame goes back to **Ehrlich & Raven's (1964)** escape-and-radiate model of butterfly/plant
coevolution and is restated at review level roughly annually (e.g. *Nat Prod Bioprospect* 3:1–7,
doi:10.1007/s13659-013-0004-0, which adds the nice twist that plants evolved **pro-toxins** that
only become toxic *after* P450 activation).

The 1990 paper is a hypothesis paper. It contains no comparative-genomic test, because in 1990 none
was possible. **Everything used to support it since was collected later, and it does not all point
the same way.**

### 2.2 Evidence for

- **Koala.** Johnson et al. (2018), *Nature Genetics* 50:1102–1111, doi:10.1038/s41588-018-0153-5
  (PMID 29967444). A eucalyptus obligate with a large **CYP2C** expansion, plus expanded vomeronasal
  and bitter-taste receptors. This is the best single case: a hyper-toxic monotypic diet with a
  matching, family-specific P450 bloom. Note it is **CYP2C, not CYP3A**.
- **Carnivora.** *Animals* 12:2821, doi:10.3390/ani12202821 (PMID 36290207). Comprehensive CYP1–3
  genomics across carnivorans finds **CYP2C and CYP3A expansion specifically in the omnivores** —
  brown bear, black bear, dog, badger — and not in the obligate carnivores. This is the closest
  thing to a direct diet↔CYP3A copy-number correlation that exists, and it is at the *level of
  feeding guild within one order*, which is the right unit of comparison.
- **Woodrats.** *Comp Biochem Physiol C* 280:109870, doi:10.1016/j.cbpc.2024.109870 (PMID 38428625).
  *Neotoma lepida* populations tolerant of creosote resin show dose-dependent induction of five
  biotransformation families, and the induced genes correspond to **species-specific duplication
  events**. Duplication and functional relevance line up in the same species.
- The SRS result of §1.4: selection in primate CYP3A is concentrated on substrate-recognition sites.

### 2.3 Evidence against, or at least against its being load-bearing

- **Copy number does not track the xenobiotic/endogenous distinction in primates.** The 2020
  *Xenobiotica* study's first prediction failed outright: EM and XM subfamilies vary equally in gene
  copy number (doi:10.1080/00498254.2020.1785580).
- **Specialised herbivores can *lose* CYP3A, not gain it.** *Mol Phylogenet Evol* 217:108550,
  doi:10.1016/j.ympev.2026.108550 (PMID 41580058). Bamboo lemurs (*Prolemur*, *Hapalemur*) eat
  cyanogenic bamboo at rates unmatched in Primates — and show **elevated gene loss in CYP2B, CYP2C,
  CYP2D, CYP2J and CYP3A** plus relaxed selection in CYP2F/CYP2J. Extreme dietary toxin load,
  *contracted* P450 repertoire. The authors' reading is that specialists adapt to *one* toxin
  (cyanide, handled by non-P450 chemistry) and shed the generalist toolkit. Diet predicts P450
  repertoire — but not monotonically, and not in the direction the warfare story implies.
- **Ecology fails to predict P450 repertoire in a clean test system.** *Genome Biol Evol* 13:evab261,
  doi:10.1093/gbe/evab261 (PMID 34850870) surveyed P450s across ten bee species chosen to contrast
  chemical exposure (orchid bees collect volatile compounds for perfume bouquets; others do not).
  **No relationship between a species' ecology and its P450 repertoire.** What they did find is a
  *structural* regularity: P450 clades divide into "stable" and "unstable", xenobiotic-metabolizing
  genes preferentially sit in unstable clades, and unstable clades show both dynamic turnover and
  adaptive signal. That is a **birth-and-death** result, not a warfare result. The same
  stable/labile dichotomy recurs in mosquitoes (*Insects* 16:184, doi:10.3390/insects16020184), where
  the authors explicitly note the demarcation between detox and essential P450s is **fuzzy**.
- **Most claimed diet-driven selection signals in humans are drift.** Roca-Umbert et al. 2019
  (doi:10.1186/s12862-019-1366-7) specifically dissolved the TAS2R38/TAS2R16 story and found no
  sweep at CYP1 or CYP2. Two of the field's favourite diet-adaptation exemplars did not survive a
  properly powered selection scan.
- **The detox framing is incoherent for CYP3A4's most famous reaction.** CYP3A4 **bioactivates**
  aflatoxin B1 to the 8,9-*exo*-epoxide that forms the guanine adduct causing hepatocellular
  carcinoma. And the positively selected residues 478/479 are precisely the ones that change
  **which aflatoxin metabolite is produced** (§1.3). An enzyme that makes a potent hepatocarcinogen
  out of a common dietary mould toxin is not obviously the product of selection for detoxification.
  Either the selection was for something else and aflatoxin activation is a by-product, or it was
  for aflatoxin handling and the epoxide is an unavoidable cost of the chemistry.

### 2.4 Verdict

**The hypothesis is real but it is not load-bearing for CYP3A4 specifically, and it is stated far
more confidently than the data support.**

What survives: (i) xenobiotic-metabolizing P450 subfamilies genuinely sit in phylogenetically
unstable clades with high duplication/loss turnover and elevated adaptive signal; (ii) CYP3A is one
of exactly four primate CYP1–3 subfamilies with detectable diversifying selection on SRSs; (iii)
within Carnivora, CYP3A copy number tracks omnivory.

What does not survive: (i) the general claim that repertoire size tracks dietary toxin load —
bamboo lemurs and orchid bees falsify it in opposite directions; (ii) the claim that a large
xenobiotic repertoire is *evidence of* dietary selection rather than of neutral birth-and-death; and
(iii) the specific claim that **CYP3A4's promiscuity** is a diet adaptation. Nobody has shown that.
The two hard selection events in primate CYP3A (Ho7, Hs4) are lineage-specific rather than
geographically structured, which Qiu et al. themselves read as pointing at **homeostasis or diet
rather than environment** — and they explicitly floated increased animal food and cooking in the
last 2 My as a driver for Hs4. That is a plausible story with no direct evidence behind it, and they
label it speculation. I will too.

**[INFERENCE]** The story that actually fits the data is **"endogenous first, promiscuity as
by-product."** CYP3A's deep function is steroid and lipid oxidation (§3). A pocket built to hold a
rigid steroid with an anionic head, with a flexible lipid-facing lid, is *incidentally* a pocket
that holds most lipophilic drugs. Selection then acted on regioselectivity (437/478/479) rather than
on breadth. On this reading the drug promiscuity was never selected for at all — it is the shadow
cast by a large hydrophobic steroid-and-lipid site with a membrane-embedded entrance. The competing
reading, that breadth was selected for, requires that breadth be both heritable and
fitness-relevant, and I can find no measurement of either.

---

## 3. Endogenous substrates — the strongest evidence about the pocket's design

This is the section that constrains §6, because it is the only evidence about what the pocket held
*before* there were drugs. The UniProt **P08684** catalytic-activity block is the best curated
inventory; every reaction below is annotated there with primary references.

### 3.1 Steroids — the core competence

- **Testosterone → 1β-, 2β-, 6β-hydroxy.** 6β-hydroxytestosterone is the canonical CYP3A4 marker
  reaction. Note there are **three** distinct positions on the same rigid scaffold — 6β on the
  B-ring β-face, 2β on the A-ring, 1β adjacent. That is a pocket holding the steroid loosely enough
  to sample several faces.
- **Progesterone → 6β- and 16α-hydroxy.** Two positions at *opposite ends* of the molecule.
  Progesterone is co-crystallised in **1W0F** (2004) and **5A1P/5A1R** (2015,
  doi:10.1021/acs.biochem.5b00510).
- **Androstenedione → 6β-hydroxy; dihydrotestosterone → 18- and 19-hydroxy.**
- **Cortisol → 6β-hydroxycortisol; cortisone → 6β-hydroxycortisone.** Urinary 6β-OH-cortisol/cortisol
  is a classical *in vivo* CYP3A phenotyping ratio. CYP3A4 is therefore a constitutive participant in
  glucocorticoid turnover, which closes a loop with §5 (glucocorticoids drive PXR).
- **Estradiol → 2-, 4-, 16α- and 16β-hydroxy; estrone → 2-, 4- and 16α-hydroxy.** Four positions on
  estradiol. 4-hydroxyestradiol is the genotoxic catechol estrogen.

**[INFERENCE]** Take the steroid set as a whole: **CYP3A4 hydroxylates the steroid nucleus at C1,
C2, C4, C6, C16, C18 and C19.** No steroidogenic P450 does that — CYP17A1, CYP19A1 and CYP11B1 each
do one or two positions with tight regiocontrol. A site that reaches seven positions on one scaffold
is not a site that **orients**; it is a site that **contains**. The geometry is: hold the hydrophobe
somewhere above the heme, and let thermal motion present whichever C–H is nearest. This is the
single most important structural inference in this document, and §4 confirms it geometrically.

### 3.2 Cholesterol and oxysterols — five positions, and a signalling product

P08684 annotates cholesterol → **4β-, 22R-, 24R-, 25- and 26-hydroxycholesterol**. Five distinct
oxidations of one substrate, spanning the B-ring (4β) and the entire isooctyl side chain (22R, 24R,
25, 26).

**4β-hydroxycholesterol is the important one.** It is essentially CYP3A4-specific, has a plasma
half-life of days, and rises 2–4-fold within a week of a strong CYP3A4 inducer, which makes it the
standard **endogenous biomarker of CYP3A4 activity** (review: *Steroids* 233–234:109819,
doi:10.1016/j.steroids.2026.109819, PMID 42251959). It is also an **LXR agonist** — i.e. CYP3A4's
product is itself a nuclear-receptor ligand feeding back into cholesterol, fatty-acid and glucose
homeostasis. An inverse association between plasma 4β-OHC and blood pressure has been reported in a
645-person cohort (*JAHA* 15:e043913, doi:10.1161/jaha.125.043913), and the 4β-OHC/cholesterol ratio
is ~25% higher in CYP3A5 expressers (*Eur J Clin Pharmacol* 82:180, doi:10.1007/s00228-026-04114-7).

**[INFERENCE]** This is not detoxification in any sense. CYP3A4 takes the most abundant sterol in
the body and converts it to a signalling molecule that regulates lipid metabolism. If you want an
"original function" candidate that is not steroid-hormone catabolism, this is it: **CYP3A4 as a
constitutive oxysterol generator sitting between cholesterol load and LXR tone.** That framing also
explains why the enzyme is expressed at ~30–40% of hepatic P450 content in an organ that is not
usually being poisoned.

### 3.3 Fatty acids, eicosanoids and endocannabinoids

Also annotated on P08684:

- Linoleate → 11-hydroxy.
- **Arachidonate → 7-, 10-, 13-hydroxy and (14R,15S)/(14S,15R)-EET.** Both epoxide enantiomers.
- **EPA → (17R,18S)-epoxide; DHA → (19R,20S)- and (19S,20R)-epoxide.**
- **Anandamide (N-arachidonoylethanolamine) → 8,9-, 11,12- and 14,15-epoxides.** Three distinct
  epoxyeicosatrienoyl ethanolamides, themselves bioactive at CB receptors.

**[INFERENCE]** Epoxidation of a *cis* double bond deep inside a 20-carbon flexible chain is a
completely different geometric demand from hydroxylating C6β of a rigid steroid. To do both, the
pocket must be able to present either a **rigid plate** or a **coiled chain** to the ferryl oxygen.
The only way to do that with one site is to make the site large, largely apolar, and flexible at the
lid — exactly what §4 describes. The anandamide epoxidations in particular say the pocket tolerates
a **polar ethanolamide head** while oxidising a distant apolar midpoint; that is a "substrate hangs
out of the entrance" geometry, and the channel structures (§4.5) show exactly that.

### 3.4 Retinoids and vitamin D — the ones with fitness consequences

- **All-*trans*-retinol → retinal; all-*trans*-retinoate → 4-hydroxyretinoate** (P08684). CYP3A7
  additionally makes 18-hydroxyretinoate (P24462). CYP3A enzymes are therefore in the retinoic-acid
  termination pathway alongside the dedicated CYP26 family.
- **Vitamin D.** The decisive evidence is human genetics, not biochemistry. The **CYP3A4 p.I301T**
  gain-of-function substitution causes **vitamin D–dependent rickets type 3 (VDDR3)** — low serum
  calcium, low 25(OH)D₃ and low 1α,25(OH)₂D₃ (case report: *Cureus* 15:e49976,
  doi:10.7759/cureus.49976). The enzymology was resolved in 2026: I301T does not merely accelerate
  4-hydroxylation; it produces **new** metabolites, 11α,25-(OH)₂D₃ and 1α,11α,25-(OH)₃D₃, with far
  lower VDR affinity. Docking implicates a new H-bond from the 3β-OH of 25(OH)D₃ to **Thr301** and
  contacts of the 25-OH with **Arg372/Glu374**, redirecting oxidation to C11α (*FEBS J* 293:660–676,
  doi:10.1111/febs.70277, PMID 41046353).

**This is the strongest single argument in the whole document that CYP3A4 is not "a drug enzyme."**
A one-residue change in the pocket of CYP3A4 produces a **Mendelian disorder of calcium
homeostasis**. Any enzyme whose active-site geometry is under that kind of physiological constraint
has been under purifying selection for an endogenous substrate. It also retro-fits the selection
scans: Qiu et al. and Thompson et al. disagreed about whether the human CYP3A sweep was driven by
**salt/water homeostasis** (CYP3A5, kidney) or **rickets/vitamin D** (CYP3A4). VDDR3 shows the
vitamin D arm is at least mechanically possible. It does not show it happened.

Note the direction, because it is counter-intuitive: **more CYP3A4 activity is bad for vitamin D
status.** The enzyme inactivates the vitamin. So any selection for vitamin D retention at high
latitude would be selection to *reduce* CYP3A4 activity, not raise it.

### 3.5 Bile acids and PXR

CYP3A4 6α-hydroxylates the cytotoxic secondary bile acid **lithocholic acid**, and LCA and
3-keto-LCA are themselves PXR ligands that induce CYP3A4 — a closed feed-forward detoxification
loop through the gut–liver axis (*Front Immunol* 16:1692684, doi:10.3389/fimmu.2025.1692684;
*J Endocr Soc* 9:bvaf119, doi:10.1210/jendso/bvaf119). This is a genuine endogenous-toxin
detoxification role that owes nothing to plants.

### 3.6 And one plant metabolite, sitting right there in the reference set

P08684 also annotates **1,8-cineole → 2-*exo*-hydroxy-1,8-cineole** and the same for 1,4-cineole.
1,8-Cineole (eucalyptol) is a bicyclic monoterpene ether — the principal volatile of *Eucalyptus*.
So the canonical, manually curated reaction list for human CYP3A4 includes a **eucalyptus terpene**
alongside cortisol and anandamide. It is a fair emblem of §2's ambiguity: the substrate range does
look like something built for plant chemistry, and it also looks like something built for steroids,
because those two requirements specify almost the same pocket.

---

## 4. Structure and plasticity

### 4.1 The structural corpus, as of today

Queried against RCSB by UniProt accession (2026-09-20):

| accession | entries | notes |
|---|---|---|
| **P08684** (CYP3A4) | **122** | all X-ray; 1.7–3.4 Å |
| **P20815** (CYP3A5) | **6** | 5VEU, 6MJM, 7LAD, 7SV2, 8SG5, 9MS2 |
| **P24462** (CYP3A7) | **2** | 7MK8, 8GK3 |
| **Q9HB55** (CYP3A43) | **0** | no structure exists |

**Every deposited CYP3A structure is a crystal structure. There is no cryo-EM CYP3A4 in the PDB.**
That is the gap the OpenADMET release fills, and it is why this challenge is not a re-scoring of
existing data.

Landmarks:

| PDB | Å | year | contents | DOI |
|---|---|---|---|---|
| **1TQN** | 2.05 | 2004 | ligand-free (Yano) | 10.1074/jbc.C400293200 |
| **1W0E / 1W0F / 1W0G** | 2.80 / 2.65 / 2.73 | 2004 | apo / **progesterone** / metyrapone (Williams) | 10.1126/science.1099736 |
| **2J0D** | 2.75 | 2006 | **erythromycin** | 10.1073/pnas.0603236103 |
| **2V0M** | 2.80 | 2007 | **ketoconazole** (2 copies) | 10.1073/pnas.0603236103 |
| **3NXU** | 2.00 | 2010 | **ritonavir** | 10.1073/pnas.1010693107 |
| **3TJS** | 2.25 | 2012 | desthiazolylmethyloxycarbonyl ritonavir | 10.1016/j.abb.2012.02.018 |
| **3UA1** | 2.15 | 2011 | **bromoergocryptine** | 10.1074/jbc.M111.317081 |
| **4I3Q** | 2.60 | 2013 | water-coordinated (near-apo) | 10.1021/jm400288z |
| **4I4G / 4I4H / 4K9T–4K9X** | 2.4–2.9 | 2013 | desoxyritonavir analogues | 10.1021/jm400288z, 10.1021/bi4005396 |
| **5A1P / 5A1R** | 2.5 / 2.45 | 2015 | **progesterone** (± citrate) | 10.1021/acs.biochem.5b00510 |
| **5TE8** | 2.70 | 2016 | **midazolam** | 10.1073/pnas.1616198114 |
| **5G5J** | 2.60 | 2017 | **metformin** | 10.1016/j.chembiol.2017.08.009 |
| **5VCC / 5VCD** | **1.70 / 1.95** | 2017 | glycerol-only; highest-resolution CYP3A4 | 10.1021/acs.biochem.7b00334 |
| **5VC0 / 5VCE / 5VCG** | 2.7 / 2.2 / 2.2 | 2017 | ritonavir, bromoergocryptine (Cys-depleted) | 10.1021/acs.biochem.7b00334 |
| **6MA6 / 6MA7 / 6MA8** | 2.18 / 2.09 / 1.83 | 2019 | metyrapone / **fluconazole** / PMSF | 10.1021/acs.biochem.8b01221 |
| **6OO9 / 6OOA / 6OOB** | 2.25 / 2.52 / 2.20 | 2019 | **mibefradil / azamulin / bergamottin** (suicide substrate) | 10.3390/ijms20174245 |
| **7LXL** | 2.75 | 2021 | **testosterone dimer** | 10.1016/j.ejmech.2021.113496 |
| **8DYC** | 2.40 | 2022 | **fluorol**, bound in the *channel* | 10.3390/ijms232012591 |
| **8SO1 / 8SO2** | 2.05 / 2.15 | 2023 | **caffeine ×3 and ×6** | 10.1016/j.jbc.2023.105117 |
| **8SPD** | 2.90 | 2023 | clotrimazole | 10.1124/dmd.123.001464 |
| **9BBB** | 2.50 | 2024 | **cobicistat** | 10.1016/j.abb.2024.110071 |
| **9GK1** | 2.95 | 2025 | serial synchrotron, **room temperature** | 10.1016/j.abb.2025.110419 |
| **9PLJ / 9PLK** | 2.72 / 2.25 | 2025 | **Δ⁹-THC** / darifenacin | 10.1016/j.jbc.2025.110709 |
| **9MS1 / 9BV5–9BVC** | 2.7–3.2 | 2025 | imidazole-series inhibitors, paired with CYP3A5 | 10.1038/s41467-025-58749-8 |
| — CYP3A5 — | | | | |
| **6MJM** | 2.20 | 2019 | substrate-free CYP3A5 | 10.1074/jbc.RA119.007928 |
| **5VEU** | 2.91 | 2017 | CYP3A5 + ritonavir | 10.1124/mol.117.109744 |
| **7LAD / 7SV2 / 8SG5 / 9MS2** | 2.2–2.8 | 2021–25 | clobetasol / azamulin / clotrimazole / SJYHJ-111 | 10.1021/jacs.1c07066; 10.1016/j.jbc.2022.101909; 10.1124/dmd.123.001464; 10.1038/s41467-025-58749-8 |
| — CYP3A7 — | | | | |
| **7MK8** | 2.15 | 2021 | first CYP3A7 (hexamutant, DTT in site) | 10.3390/ijms22115831 |
| **8GK3** | 2.60 | 2023 | **CYP3A7 + four DHEA-S** | 10.1016/j.jbc.2023.104993 |

### 4.2 Apo vs holo — and a genuine, unresolved disagreement

The two 2004 papers reached opposite conclusions **in the same year on the same protein**:

- **Williams et al., *Science* 305:683–686** (1W0E/1W0F/1W0G): "a surprisingly **small** active site,
  with **little conformational change** associated with the binding of either compound."
- **Yano et al., *JBC* 279:38091–38094** (1TQN): "a **relatively large** substrate-binding cavity…
  consistent with its capacity to oxidize bulky substrates such as cyclosporin, statins, taxanes and
  macrolide antibiotics," and larger near the heme iron than CYP2C8's.

Then **Ekroos & Sjögren (2006), *PNAS* 103:13682–13687, doi:10.1073/pnas.0603236103** broke the tie
by co-crystallising big ligands: "**In contrast to previous reports, the protein undergoes dramatic
conformational changes upon ligand binding with an increase in the active site volume by >80%.**"
Ketoconazole binds as **two copies**; erythromycin shows multiple binding modes; and the two
complexes are **two distinct open conformations**, i.e. the expansion is ligand-shaped, not a single
"open state." Their closing warning has aged extremely well: the flexibility "challenges any attempt
to apply computational design tools without the support of relevant experimental data."

Then it reversed again. **Sevrioukova & Poulos (2017), *PNAS* 114:486–491,
doi:10.1073/pnas.1616198114 (5TE8)** found that midazolam **contracts** the site by ~20%: Leu216
sweeps **12 Å** toward the drug's fluorophenyl ring, the F′ and G helices partly unwind, the F helix
gains a turn, G′ rotates 20°, the B–C loop moves 5 Å out, the C-terminal loop moves 3 Å in, and the
D/E/H helices and the I-helix N-terminus shift 0.6–1.6 Å. Midazolam's imidazole nitrogen H-bonds
**Ser119**, and its C1 — the major site of metabolism — sits **4.4 Å** from the iron.

**The honest summary is: there is no single CYP3A4 active-site volume.** Reported behaviour spans a
~20% contraction (midazolam) to a >80% expansion (ketoconazole/erythromycin) relative to ligand-free,
with two independent open conformations at the top end. Absolute Å³ values quoted in the secondary
literature (~950–1400 Å³ apo, ~2000 Å³ expanded) come from different cavity-detection algorithms with
different probe radii and are **not comparable across papers**; the 2017 midazolam paper reports its
result only as a percentage and gives no absolute volume at all.

**[INFERENCE]** The 2004 disagreement was not a dispute about the protein; it was a dispute about
which crystal form and which ligand. A protein whose reported cavity volume depends that strongly on
what is inside it is a protein whose **apo structure carries little information about its holo
structure**. For the challenge, that is the mechanism behind claim A in CLAUDE.md: the backbone is
right and the ligand orientation is wrong, because the pocket that would orient the ligand *is
created by the ligand*.

### 4.3 The F/G region and the phenylalanine roof

The moving part is the **F–F′–G′–G** fragment (~residues 202–260), which in the membrane-bound
enzyme is the part that **inserts into the bilayer**. Above the heme sits an ordered hydrophobic
cluster — the **phenylalanine cluster** — described in Williams et al. 2004 and repeatedly since:
**Phe57, Phe108, Phe213, Phe215, Phe219, Phe220, Phe241, Phe304**, with Phe189/Phe203 and
Leu210/Leu211/Leu249/Val253/Leu295/Met296 completing the patch in the CYP3A7 structures.

**Phe304 is the gatekeeper.** It adopts different rotamers in CYP3A4 and CYP3A5 depending on ligand
size, and controls access of large ligands to the active site proper. In **all three** CYP3A7
structures it is locked in the "up" rotamer packed into the F′/G′ roof — and this is the structural
explanation offered for why CYP3A7 is the poor catalyst of the family: **"decreased structural
plasticity rather than the active site microenvironment defines the ligand binding ability of
CYP3A7"** (*IJMS* 22:5831, doi:10.3390/ijms22115831). Elongated β3–β4 strands, an H-bond across the
substrate channel and steric constraints in the C-terminal loop add further rigidity, changing the
heme environment and impeding the low→high spin transition essential for reactivity.

Two residues from that cluster, **F108 and F215, are on Qiu's positively selected list for branch
Ho7** (§1.3). The hominoid burst on CYP3A7 hit the roof of the pocket.

**[INFERENCE] The functional axis separating CYP3A4 from CYP3A5 and CYP3A7 is plasticity, not
chemistry.** CYP3A5's site is "**taller and narrower**" than CYP3A4's, forcing ritonavir into a
distinctly different conformation, driven by F–G-region differences and β-sheet ionic interactions
(*Mol Pharmacol* 93:14–24, doi:10.1124/mol.117.109744). CYP3A7's site has almost the same residues
as CYP3A4's but cannot move. Consistently with that, the UniProt reaction inventory for CYP3A5
(P20815) is close to a strict **subset** of CYP3A4's: estradiol 2-/4-OH, estrone 2-/4-OH,
testosterone 6β, androstenedione 6β, progesterone 6β, retinol, retinoate. It has **none** of
CYP3A4's cholesterol, arachidonate, EPA/DHA, anandamide or linoleate chemistry. **The paralog that
cannot deform is the paralog that lost the flexible-chain substrates.**

### 4.4 Multiple occupancy — the defining kinetic and structural property

This is where CYP3A4 is genuinely unusual among P450s, and it is now directly visible.

**Structures with more than one copy of the same ligand:**

- **2V0M** — two ketoconazole (2006).
- **8SO1 / 8SO2** — **three and six caffeine molecules** (2023, doi:10.1016/j.jbc.2023.105117). In
  the ternary complex one caffeine is positioned for C8-hydroxylation. In the senary complex **three
  caffeines stack parallel to the heme**, with the proximal one poised for 3-N-demethylation — but
  the stack's own hydrophobic interactions are extensive enough that product release may be blocked.
  Caffeine is simultaneously found **in the substrate channel and on the outer peripheral surface**.
  Aromatic stacking dominates at all three sites; mutagenesis confirmed active-site **Arg212**,
  intrachannel **Thr224** and peripheral **Phe219**.
- **7LXL** — a covalently tethered testosterone dimer (2021).
- **8GK3** — CYP3A7 with **four DHEA-S molecules** (2023, doi:10.1016/j.jbc.2023.104993). Twelve
  chains in the asymmetric unit, six ligand-free and six carrying four DHEA-S each at full occupancy:
  - **DHEA-S1** in the active site proper, **4.1 Å above the heme iron**, C17 ketone toward the iron,
    C3 **sulfate** H-bonded (2.5–3.1 Å) to the **backbone amide of Gly481** in β4. The site of
    metabolism, **C16, is 5.8 Å from the iron**. Contacts: Y57, R105, S119, N214, P215, L216, A305,
    T309, V369, A370, G480, G481, L482 — **almost entirely hydrophobic plus one backbone H-bond**.
  - **DHEA-S2** stacked above it, C17→Fe **8.7 Å**, steroid long axis parallel to the F helix,
    sandwiched between **Pro215** and **Phe108**, sulfate H-bonded to **Lys224** (and Arg106 in one
    chain).
  - **DHEA-S3** in the access channel between the A′ and F′ helices, sulfate reaching the surface.
  - **DHEA-S4** on the **F′–G′ hydrophobic face that is normally embedded in the ER membrane**.

  Cα RMSD between liganded and unliganded chains is only 0.66–0.85 Å — **four steroids bind with
  almost no backbone rearrangement**, the differences confined to A/A′/B/F′/G′/G and β1–2.

**The kinetics say the same thing, quantitatively.** Denisov, Sligar and colleagues, using monomeric
CYP3A4 in Nanodiscs with a defined 1:1 reductase complex (*JBC* 282:7066–7076,
doi:10.1074/jbc.m609589200), resolved the contributions of the **one-, two- and three-testosterone**
species:

- **1 TST bound:** accelerates NADPH consumption but produces **essentially no product** — completely
  uncoupled.
- **2 TST bound:** product formation reaches its **maximum rate**.
- **3 TST bound:** no further rate increase, but **coupling efficiency improves again**.

Two follow-ups matter for interpretation. Careful Nanodisc work on the classic
α-naphthoflavone/testosterone heterotropic pair found **negligible binding cooperativity**: the
apparent activation is the additive effect of a second molecule occupying the cavity, not a specific
allosteric interaction (*JBC* 286:5540–5545, doi:10.1074/jbc.M110.182055; *Arch Biochem Biophys*
488:146–152, doi:10.1016/j.abb.2009.06.013). And Denisov & Sligar's synthesis (*Arch Biochem
Biophys* 519:91–102, doi:10.1016/j.abb.2011.12.017) frames the whole phenomenon as **"functional
cooperativity in monomeric proteins"** that "does not require substantial binding cooperativity" —
what it requires is one or more binding sites of higher affinity than the single catalytically
productive site near the iron.

**[INFERENCE] This is the mechanistic heart of CYP3A4 and it should change how we score poses.** One
substrate in a cavity that big rattles; the ferryl oxidant is consumed making water instead of
product. The *second* molecule is not an allosteric effector in the classical sense — it is a
**steric chaperone**, a wall that converts a loose cavity into a tight one. Yano's 2004 abstract said
this before any structure showed it: "The lower constraints on the motions of small substrates near
the site of oxygen activation may diminish the efficiency of substrate oxidation, which may, in
turn, be improved by space restrictions imposed by the presence of a second substrate molecule."
Twenty years later, 8SO1/8SO2/8GK3 show the second molecule.

Consequence for the challenge: **a single-ligand co-folded pose of a small CYP3A4 substrate may be
modelling a catalytically irrelevant species.** If the crystallographic or cryo-EM ligand density
corresponds to a multiply-occupied or partially-occupied site, and the prediction places one ligand
in the middle of the cavity, LDDT-PLI will be wrong for a reason that has nothing to do with the
model's physics. This is a concrete, testable explanation for part of the oracle gap in FINDING 011,
and it fits OpenADMET's own statement that several of their ligands show **multiple mutually
exclusive conformations** rather than one pose.

### 4.5 The peripheral site and the channel — and the argument about them

**The peripheral site.** Williams et al. 2004 found progesterone in **1W0F** bound ~17 Å from the
heme, **above the phenylalanine cluster**, on the outer surface, and proposed it as a site for
"initial recognition of substrates or allosteric effectors." It has been argued about ever since —
crystallisation artefact vs. real effector site. Evidence that it is real: caffeine occupies the same
peripheral region in 8SO1/8SO2, and mutation of **Phe219** there measurably affects caffeine
association (doi:10.1016/j.jbc.2023.105117); steroid bioconjugation to an engineered allosteric site
changes substrate binding and coupling (PMID 29958895). Evidence that it is not *required*: the 5TE8
paper reports that the progesterone site **disintegrates on midazolam binding**. My reading is that
it is a **real but low-specificity lipid-facing surface pocket**, not a dedicated regulatory site —
see below.

**The channel is itself a binding site.** The fluorol structure **8DYC**
(doi:10.3390/ijms232012591) puts the ligand **not in the active site but in the substrate channel**,
stabilised by H-bonds to polar channel residues including **Thr224** and **Arg372**, and the authors
argue for an explicit sequential model: hydrophobic ligands **dock first on the nearby peripheral
surface, migrate into the channel, and only then enter the active site**. The DHEA-S structure shows
all three stations occupied simultaneously (S4 surface, S3 channel, S1/S2 active site) — a
**freeze-frame of the whole pathway**. 8GK3 also maps four exit channels from the active site:
channel 1 (the main one, holding DHEA-S3), channel 2 between the B–B′ and B′–C loops, channel 3
between the F and I helices, and channel 4 between the B–B′ loop and β1–1.

**Dynamics agree.** Long unbiased MD (*Commun Chem* 9:17, doi:10.1038/s42004-025-01815-5, PMID
41513816) captured spontaneous binding from bulk to the catalytic centre: the ligand sits on the
surface for **>5 µs** with the **F–F′ loop down**, the loop then lifts to admit the ligand and
*restricts its exit*, the ligand rearranges internally, and the loop closes again. The final pose
matched the crystallographic one and the Markov-state-model K_d matched experiment. **The F–F′ loop
is a one-way gate, not a passive lid.**

**The most counter-intuitive dynamics result.** HDX-MS on nanodisc-embedded CYP3A4 with a ligand
panel, cross-checked by PCA of MD trajectories (*J Inorg Biochem* 244:112211,
doi:10.1016/j.jinorgbio.2023.112211, PMID 37080138): effects concentrate in the **F and G helices**
(most ligands **increase** F-helix and connecting-loop flexibility while **decreasing** C-terminal
G-helix flexibility) and propagate to E–F–G, C–D and H–I, with ligand-specific differences in the
A″–A′ loop, B–C region, E helix, K–β1 region, proximal loop and C-terminal loop. Type II
(heme-ligating) ligands give correlated responses in the C–D region and G-helix C-terminus. And,
against the textbook: **"a wide range of ligands modestly increase CYP3A4 dynamics throughout the
protein, including effects remote from the active site"** — the opposite of the usual
binding-rigidifies-the-site paradigm.

**[INFERENCE] CYP3A4 does not have one binding site; it has a gradient.** Surface → channel →
cavity, with affinity and specificity both increasing along it, and a gate in the middle whose
opening is stochastic on a multi-microsecond timescale. "The" binding site is a convenient fiction.
Any scoring function that assumes a single well-defined pocket with single occupancy is modelling
something that does not exist. This also reframes the peripheral-site debate: it is not an allosteric
site *or* an artefact — it is **the first rung of the ladder**, and whether a given ligand is
resolved there depends on whether that rung happens to be its deepest well.

### 4.6 The heme, and what coordination does and does not tell you

Cys442 thiolate, proximal. The resting state is six-coordinate low-spin with an axial water;
substrate binding displaces it, shifting to high spin (Type I difference spectrum) and raising the
reduction potential enough for the first electron transfer. Type II ligands — azoles, pyridines,
imidazoles, thiazoles — **coordinate the iron directly through nitrogen**, giving the red-shifted
Soret and, in general, potent inhibition.

The CYP3A4 corpus is rich in these: ketoconazole (2V0M), fluconazole (6MA7), clotrimazole (8SPD),
metyrapone (1W0G, 6MA6), ritonavir (3NXU, thiazole N), cobicistat (9BBB, thiazole N), azamulin
(6OOA), and the whole bipyridyl/Ir(III) photocaged series (7UAY/7UAZ/8EW*/8EXB).

The cobicistat structure is the best worked example of how *little* coordination alone buys you
(*Arch Biochem Biophys* 758:110071, doi:10.1016/j.abb.2024.110071): cobicistat ligates the heme via
its thiazole nitrogen exactly as ritonavir does, but binds with **2-fold lower affinity and rate**
(K_s 0.030 µM, 0.72 s⁻¹) because its bulky morpholine clashes with the **F–F′ connector**, which
becomes disordered, and it cannot H-bond **Ser119**. Yet its IC₅₀ (0.24 µM) equals ritonavir's
(0.22 µM), rescued by phenyl-group hydrophobic and aromatic contacts. **Coordination, the Ser119
H-bond and the F–F′ contacts are three separable contributions, and they trade off against each
other.**

**[INFERENCE]** This is direct structural corroboration of this repo's FINDING 008/019 and claim I:
the iron anchor is a **saturated** constraint. Nearly every Type II CYP3A4 ligand reaches
coordination geometry; the discrimination lives in the F–F′ contacts and the Ser119 H-bond, i.e.
**above and to the side of the iron, in the part that moves**. We measured that empirically (84% of
Boltz poses inside the coordination window; the angle term inert to four decimals) before reading it
here. The structural literature says the same thing for a physical reason.

### 4.7 Oligomerisation and redox partners

CYP3A4 functions as a monomer in Nanodiscs and the 1:1 CYP3A4:CPR complex is fully competent
(doi:10.1074/jbc.m609589200), so oligomerisation is not required for catalysis, though it occurs in
microsomes. Cytochrome b₅ stimulates most CYP3A4 reactions; CYP3A7 is reported to be less
b₅-dependent, in a substrate-specific way whose molecular basis is unknown, and Qiu et al. note that
**eight of the eighteen Ho7-selected residues lie in the first 100 amino acids** — the region that
governs ER targeting and redox-partner interaction. There is also mounting evidence that the
reductase FMN domain is an **allosteric modulator** and not merely an electron donor: MS footprinting
of CYP2A6 ± FMN domain found **increased** surface exposure rather than interface protection
(*DMD* 54:100210, doi:10.1016/j.dmd.2025.100210).

**[INFERENCE]** If the hominoid CYP3A7 burst hit the N-terminus and redox interface as hard as it hit
the pocket roof, then part of what distinguishes the paralogs is **how they are wired**, not what
they bind. That is invisible to any structure of the soluble catalytic domain — and invisible to
co-folding, which models a truncated construct with no membrane and no CPR.

---

## 5. Regulation and context

### 5.1 The induction machinery

CYP3A4 is the most inducible drug-metabolizing enzyme in the human genome, and the architecture is
layered:

- **Proximal promoter ER6** (~−170 bp) — the original PXR response element.
- **XREM** (xenobiotic-responsive enhancer module, ~−7.2 to −7.8 kb) — DR3 + ER6.
- **CLEM4** (~−11 kb).
- **PXR (NR1I2)** is the dominant sensor; **CAR (NR1I3)** contributes; **VDR** mediates
  1,25-(OH)₂D₃ induction, prominently in intestine; **FXR**, HNF4α, C/EBP and glucocorticoid
  signalling (which raises PXR itself) modulate.

Two chromatin-conformation studies have since rewritten the map:

- **4C/3C + CRISPR deletion in primary human hepatocytes** (*Pharmacogenet Genomics* 30:107–116,
  doi:10.1097/fpc.0000000000000402) found **four** regions contacting the CYP3A4 promoter. R2 is
  XREM/CLEM4; **R4 is novel and does something remarkable — deleting it *increases* CYP3A4 and
  *decreases* CYP3A43**, apparently by competitive domain–domain interaction within the cluster.
  rs62471956 in R4 is in **complete LD with CYP3A4\*22**.
- A shared **distal regulatory region (DRR)** contacts the CYP3A4 promoter and, when CRISPR-deleted,
  **reduces CYP3A4, CYP3A5 and CYP3A7 together** (*Clin Transl Sci* 15:2720–2731,
  doi:10.1111/cts.13398). rs115025140 (African-specific) raises CYP3A4 mRNA 1.8× and protein 1.6×;
  rs776744/rs776742 raises CYP3A5 mRNA 1.39×.

**[INFERENCE] The four CYP3A genes are not four independent units — they are one regulatory domain
with shared enhancers and internal competition.** That reframes §1: "gene duplication followed by
subfunctionalisation" understates it. The paralogs share *cis* elements, and they have also traded
*coding* sequence by gene conversion. They are better modelled as a single locus with four
alternative catalytic outputs than as four genes. And it supplies the missing mechanism for the
CYP3A7 fetal switch: shifting which promoter in a shared-enhancer cluster wins is a far cheaper
evolutionary move than building a new regulatory program.

Downregulation matters too: **IL-6 and inflammation suppress CYP3A4**, which is why acute illness
changes drug clearance.

### 5.2 Zonation and tissue context

CYP3A4 is a **pericentral (zone 3)** hepatocyte enzyme, in the same Wnt/β-catenin-controlled
compartment as CYP2E1 and CYP1A2 — the low-oxygen, high-xenobiotic-metabolism end of the sinusoid.
The modern spatial-omics literature (e.g. *Nature* 653:1148–1157, doi:10.1038/s41586-026-10377-y;
*Nat Metab* 8:741–756, doi:10.1038/s42255-026-01459-2; review *Annu Rev Pathol* 21:185–212,
doi:10.1146/annurev-pathmechdis-042624-091820) confirms strong zonation of ~half the measurable
hepatocyte proteome, reports that **key functions are pericentrally shifted in humans relative to
mice**, and shows that **pericentral proteins are the most vulnerable to loss of zonation** in
disease. Practical corollary: "CYP3A4 abundance" in a homogenised liver is an average over a steep
spatial gradient, and it collapses in MASLD/MASH.

**Intestine vs liver.** CYP3A4 is the dominant enterocyte P450, concentrated in the proximal small
intestine (duodenum/jejunum), falling along the gut and effectively absent from colon. It is the
main reason oral CYP3A substrates have low bioavailability, and the target of the grapefruit-juice
interaction — furanocoumarins such as **bergamottin** are mechanism-based inactivators, and
bergamottin is co-crystallised as a suicide substrate in **6OOB**. Intestinal and hepatic CYP3A4 are
not regulated identically (VDR is relatively more important in gut), and intestinal ontogeny is
slower and less certain than hepatic (*CPT PSP* 13:1570–1581, doi:10.1002/psp4.13192). UniProt
P08684 additionally records protein-level expression in prostate, bile duct, nasal mucosa, kidney,
adrenal cortex, gallbladder, pancreatic intercalated ducts, parathyroid chief cells and ovarian
corpus luteum — a distribution that is hard to explain as a dedicated first-pass barrier.

CYP3A5 is relatively more abundant in **kidney** and intestine than in liver — which is the whole
basis of the salt-homeostasis hypothesis below.

### 5.3 Population genetics — hard numbers

All frequencies from **gnomAD v4 genomes**, queried 2026-09-20.

**CYP3A5\*3 — rs776746 (chr7:99,672,916 T>C, GRCh38).** Intron-3 variant creating a cryptic splice
site → exon 3B inclusion → premature stop → no functional protein. Global AF of the `*3` (C) allele
= **0.7286** (110,806 / 152,086).

| population | `*3` frequency |
|---|---|
| Non-Finnish European | **0.931** |
| Finnish | 0.935 |
| Amish | 0.933 |
| Ashkenazi Jewish | 0.918 |
| Middle Eastern | 0.895 |
| Admixed American | 0.780 |
| East Asian | 0.724 |
| South Asian | 0.700 |
| **African / African-American** | **0.305** |
| HGDP Bantu (South Africa) | 0.188 |

A **three-fold** frequency difference between African and European ancestries, and one of the
steepest clines in the pharmacogenome. This is the observation **Thompson, Kalinowski, Bamshad,
Rocchi, Wojnowski et al. (2004), *Am J Hum Genet* 75:1059–1069, doi:10.1086/426406 (PMID 15492926)**
built the selection argument on: >1,000 individuals from 52 worldwide populations, an excess of rare
variants and a homogeneous group of high-frequency long-range haplotypes outside Africa, `*3`
frequency **significantly correlated with distance from the equator**, and — crucially — the
**unlinked** AGT M235T hypertension variant showing the *same* geographic pattern and correlating in
frequency with CYP3A5 `*1/*3`. Two unlinked salt-handling loci with matching clines is much harder
to explain by drift than one. The region proposed as the sweep target runs from ~40 kb upstream of
the CYP3A4 promoter to CYP3A4 intron 6.

**Other alleles (gnomAD v4):**

| allele | rsID | AFR | NFE | EAS | note |
|---|---|---|---|---|---|
| **CYP3A4\*1B** (−392A>G) | rs2740574 | **0.641** | 0.035 | 0.002 | promoter; in LD with CYP3A5\*1 |
| **CYP3A4\*22** | rs35599367 | 0.009 | **0.049** | 0.000 | intron 6; reduced expression |
| CYP3A distal-regulatory | rs776744 | **0.585** | 0.086 | 0.279 | raises CYP3A5 mRNA 1.39× |

**The competing explanations, unresolved.** Thompson et al. argued **salt and water retention** —
CYP3A5 is expressed in kidney and plausibly affects sodium handling, so losing it at high latitude
(where heat-driven salt conservation is not needed) is tolerable or beneficial. Qiu et al. note that
others have proposed **rickets/vitamin D** instead, since CYP3A4 inactivates vitamin D metabolites
(§3.4) and high-latitude populations are vitamin-D-limited. Roca-Umbert et al. 2019 detect the sweep
on **CYP3A4 and CYP3A43** in Europeans (doi:10.1186/s12862-019-1366-7), not CYP3A5 — which fits the
vitamin D arm better than the salt arm. Both remain live; the honest position is that a sweep at the
CYP3A locus outside Africa is well supported and **its target gene and its phenotype are not
settled**.

A caution worth stating plainly: this literature is entangled with the "slavery hypertension
hypothesis," which is contested on historical and population-genetic grounds. The allele frequencies
are solid data. The adaptive narrative attached to them is not.

### 5.4 Post-translational control — and what it implies

- **ERAD / ubiquitination.** CYP3A4 is ubiquitinated by the ER-membrane E3 ligases **gp78/AMFR** and
  **CHIP**; knocking both down stabilises CYP3A4 and the stabilised protein **remains catalytically
  functional** (*Cancer Biol Ther* 11:549–551, doi:10.4161/cbt.11.6.14834). Degradation is a control
  point, not quality control of misfolded protein.
- **Circadian proteolysis.** *Biochem Pharmacol* 253:118330, doi:10.1016/j.bcp.2026.118330 (PMID
  42567488): in clock-synchronised HepaRG cells, CYP3A4 protein and activity oscillate **in phase
  with each other but with no mRNA rhythm**; expressed from a constitutive CMV promoter, CYP3A4
  protein stability still varies with time of day; **gp78 oscillates antiphase** to CYP3A4, and only
  catalytically active gp78 abolishes the rhythm. **CYP3A4 activity is under post-transcriptional
  circadian control via targeted proteolysis.**
- Also reported: phosphorylation (PKA/PKC sites), heme-adduct suicide inactivation (bergamottin,
  6OOB), and irreversible modification of **Lys257** by a reactive ritonavir metabolite — where
  K257A turned out to be **functionally important for CYP3A4 allosterism** independently of
  inactivation (doi:10.1016/j.abb.2024.110071).

**[INFERENCE] What PTM tells us about non-metabolic roles is mostly negative, and that is useful.**
There is no credible moonlighting function for CYP3A4 — no reported nuclear localisation, no
scaffolding role, no protein–protein signalling function. What the PTM data say instead is that the
cell invests heavily in **titrating how much CYP3A4 exists, minute to minute**, on top of the
transcriptional induction machinery of §5.1. Two independent, energetically expensive control layers
on abundance and zero on localisation implies the enzyme's output is something whose *rate* matters
continuously — consistent with the 4β-OHC/LXR and bile-acid/PXR loops of §3, and inconsistent with a
purely reactive "poison arrives, enzyme responds" model.

### 5.5 CYP3A43 — not a pseudogene, just a bad enzyme

Qiu et al. treated it as pseudogenising. That is half right. CYP3A43 produces mostly aberrant
transcripts in human and chimp and is frameshifted in some rhesus, yet **its coding region evolves
under strong purifying constraint** (low terminal-branch ω in Qiu et al. Fig. 5). And recombinant
CYP3A43 is demonstrably catalytic: it makes 4-hydroxyalprazolam and α-hydroxyalprazolam, just slowly
(*Biomedicines* 10:3022, doi:10.3390/biomedicines10123022), and it metabolises **olanzapine**, with
rs472660 apparently raising its expression and the prostate-cancer-associated CYP3A43.3 variant
showing *increased* activity on some substrates (*Xenobiotica* 52:413–425,
doi:10.1080/00498254.2022.2078751). Its highest expression is **prostate**, not liver (UniProt
Q9HB55). Three CYP3A43 point mutants (L293P, T409R, P340A) redirect alprazolam metabolism to a novel
**5-N-oxide** product and raise 4-hydroxylation 4–6-fold.

**[INFERENCE]** Purifying selection on the coding sequence of a gene that is barely transcribed in
liver, with peak expression in prostate, is a signature of a **tissue-restricted endogenous
function** that nobody has identified. Three single mutations that unlock 4–6-fold activity and a
new reaction say the protein is **one or two substitutions away from being a competent enzyme**. The
most economical reading is that CYP3A43 is not decaying; it is **specialised and lowly expressed**,
and its substrate has not been found. It is the obvious place to look for the CYP3A subfamily's
unknown endogenous chemistry.

---

## 6. The inference — what the pocket is built to bind

Everything below is **[INFERENCE]** unless it carries a citation. It is reasoning from shape and
from selective pressure, deliberately not from a list of known drugs.

### 6.1 Derivation

Six constraints fall out of §§1–5, and together they nearly specify the pocket.

1. **Seven hydroxylation positions on one steroid nucleus (§3.1).** A site that reaches C1, C2, C4,
   C6, C16, C18 and C19 of a rigid tetracycle is not orienting its substrate. It is **containing**
   it. The cavity must be substantially larger than a steroid.
2. **But one substrate in that cavity is catalytically dead (§4.4).** The 1-TST species consumes
   NADPH and makes no product. The cavity is therefore *too* large for one small ligand — its
   productive state is **filled**, whether by a second copy, a co-substrate, or one large ligand.
3. **Both a rigid plate and a coiled C20 chain must reach the ferryl oxygen (§3.3).** Steroid
   6β-hydroxylation and anandamide 11,12-epoxidation cannot be served by the same *static* shape.
   The site must be **deformable on the ligand's timescale**, and the deformable part must be the
   lid, because the floor carries the chemistry (§1.3).
4. **The entrance is in the membrane (§4.5).** The F′–G′ face is bilayer-embedded; DHEA-S4 sits on
   it; ligands dock on the surface, migrate into the channel, then enter the cavity through a
   one-way F–F′ gate that opens on a multi-microsecond timescale. **Substrates arrive laterally from
   the lipid phase, not from water.**
5. **The iron is a cheap, saturated anchor (§4.6).** Almost any sp²-N heterocycle reaches it.
   Discrimination lives in the F–F′ contacts and the Ser119 H-bond, not at the metal.
6. **Anionic heads are tolerated, but only at the periphery (§4.4).** Every DHEA-S sulfate in 8GK3
   H-bonds a **backbone amide or a surface lysine/arginine** (Gly481 backbone, Lys224, Arg106,
   Arg243), never a buried salt bridge. The cavity proper is hydrophobic; charge is handled at the
   rim.

### 6.2 The physicochemical envelope

Putting those together, CYP3A4's pocket is built to bind:

- **Size.** Comfortable from ~180 Da (caffeine, metformin) to ~1200 Da (bromoergocryptine,
  cyclosporin), with the **productive optimum around 300–600 Da for a single ligand**, or two
  180–350 Da ligands stacked. Below ~250 Da, expect either multiple occupancy or uncoupling.
- **Lipophilicity.** Strongly biased lipophilic — **logP ~2–6**, with a large apolar surface. The
  membrane-lateral entrance means **partitioning into the bilayer is part of the binding event**, so
  effective affinity is a product of membrane partition coefficient and pocket affinity, not pocket
  affinity alone.
- **Shape.** Either (a) a **rigid fused polycyclic plate** — steroid, macrolide macrocycle, ergoline,
  triazolobenzodiazepine, cannabinoid — or (b) a **long flexible chain that can coil**. Both work;
  what fails is a small rigid sphere with no handle.
- **Polarity distribution.** **Apolar core, polar terminus.** One or two H-bond donors/acceptors
  placed *distally* from the site of oxidation, to be caught by Ser119, Arg105/Arg106, Arg212,
  Thr224, Arg372 or a backbone amide at the rim. Anandamide's ethanolamide, DHEA-S's sulfate,
  midazolam's imidazole N (→Ser119) and fluorol's polar group (→Thr224) are all the same motif.
- **Aromaticity.** Positively valued — the Phe roof means **π-stacking is a primary recognition
  mode**, and caffeine's three sites are all stacking-dominated.
- **Basicity.** Neutral-to-weakly-basic favoured; strong acids are disfavoured *inside* but tolerable
  if the acid sits at the rim.
- **Site of metabolism.** Sterically accessible, electronically activated where possible — allylic
  positions, benzylic positions, N-methyls, sp³ C–H β to a carbonyl — at **4.0–5.5 Å from the iron**
  when the rest of the molecule is anchored (5TE8: C1–Fe 4.4 Å; 8GK3: C16–Fe 5.8 Å).

### 6.3 What makes a good substrate versus a bad one

**Good:** fills the cavity, presents an apolar face to the Phe roof, keeps its polar atoms at the
rim, and has a metabolically soft position that happens to sit near the iron when so anchored. Note
what is *not* on that list — a specific pharmacophore. There isn't one.

**Bad, and why, in four distinct failure modes:**

- **Too small and rigid** (metformin, 5G5J) — binds, doesn't fill, uncouples. Metformin is
  essentially a CYP3A4 non-substrate *despite* binding.
- **Too polar or too charged overall** — never partitions into the membrane, so never reaches the
  lateral entrance. This is why CYP3A4 barely touches highly polar drug space even when the
  molecules are the right size.
- **Too rigid and too large** — the fluorene analogue in the ritonavir series was the weakest binder
  of its whole series precisely because it "is too large to allow unrestricted access"
  (doi:10.1111/cbdd.70043). Size alone is fine; **size without a way to thread the channel** is not.
- **Wrong hinge geometry** — cobicistat's morpholine clashes with the F–F′ connector and disorders
  it (doi:10.1016/j.abb.2024.110071). The lid is the discriminator.

### 6.4 The chemical space we should accept but have never assayed

This is the prediction. Each entry follows from a specific constraint above.

1. **Sulfated and phosphorylated lipophilic metabolites.** From §6.1(6) and 8GK3: the pocket
   demonstrably accommodates a **sulfate ester** on a steroid, held by a backbone amide. CYP3A4 is
   routinely assayed against parent drugs, almost never against **phase-II conjugates**. Steroid
   sulfates beyond DHEA-S, sulfated bile acids, lysophosphatidic acids, sphingosine-1-phosphate and
   sulfated plant glycosides are all in envelope. **Prediction: CYP3A4 oxidises anionic
   lipid-conjugate space far more broadly than the literature records, because nobody looks —
   conjugates are treated as endpoints, not substrates.**
2. **Endocannabinoid- and oxylipin-adjacent chemistry.** CYP3A4 already makes three anandamide
   epoxides and EPA/DHA epoxides (§3.3). **Prediction: 2-AG, N-acyl amino acids, N-acyl taurines,
   oleoylethanolamide, resolvin/protectin precursors and the endovanilloids are substrates.** These
   are exactly "apolar chain + polar head", and they are essentially unassayed against CYP3A4
   because they belong to a different literature.
3. **Macrocycles and beyond-Rule-of-5 space.** Erythromycin (2J0D) and bromoergocryptine (3UA1)
   establish that the expanded conformation takes ~700–1000 Da macrocycles. **Prediction: modern
   bRo5 therapeutics — macrocyclic peptides, molecular glues, PROTACs, cyclic depsipeptides — are
   CYP3A4 substrates at a much higher rate than Ro5 chemistry, and, because two ligands fit, they
   will show pronounced atypical kinetics.** ARV-471 already appears in CYP3A4 DDI work
   (doi:10.1002/cpt.70000). This is the most commercially consequential item on the list.
4. **Covalent / cryptic-electrophile natural products.** Bergamottin's suicide mechanism (6OOB) is
   a template: bind in the lipophilic pocket, get oxidised, alkylate the heme. **Prediction: the
   furanocoumarin, furanoterpenoid, methylenedioxyphenyl and terminal-acetylene classes contain many
   more unrecognised mechanism-based CYP3A4 inactivators than the ~dozen currently known**, and the
   assay that would find them (time-dependent inhibition across natural-product libraries) is rarely
   run.
5. **Two-ligand chemical space — the space with no assay at all.** §4.4's structures and kinetics
   say the productive species is often binary. **Prediction: there exist pairs (A, B), neither
   individually a good CYP3A4 substrate, where B's presence converts A from uncoupled to
   productive** — a specific, testable, currently unmeasured class of food–drug and drug–drug
   interaction. The standard DDI panel measures inhibition, not **reciprocal activation of a
   non-substrate**. Caffeine's simultaneous occupancy of all three stations (8SO1/8SO2) makes it the
   obvious probe partner, and it is in almost everyone's blood.
6. **The CYP3A43 substrate.** From §5.5: a prostate-enriched paralog under purifying selection with
   no known physiological substrate. **Prediction: its substrate is a prostate/androgen-axis lipid or
   steroid** — something in the DHT, oxysterol or steroid-sulfate neighbourhood — and finding it
   would identify a CYP3A function invisible in liver.
7. **[Weaker inference] Sterol side-chain space generally.** CYP3A4 oxidises cholesterol at C22,
   C24, C25 and C26 (§3.2). Phytosterols, ergosterol, oxysterol drugs, cholesterol-tethered
   conjugates and bile-alcohol intermediates are all the same shape. Essentially nothing here has
   been assayed against CYP3A4 as a *substrate* rather than as an inhibitor.

### 6.5 The one-paragraph version

CYP3A4's pocket was not built to bind drugs and was probably not built to bind plants. It was built
to **hold a rigid, membrane-partitioned, mostly apolar molecule of 300–600 Da loosely enough to
oxidise it at several different positions, with its polar end tethered at the rim and its entrance
opening sideways into the lipid bilayer.** Steroids, oxysterols, bile acids, retinoids and
polyunsaturated acyl chains satisfy that description, and so, entirely incidentally, does most of
the small-molecule pharmacopoeia. Primate selection then tuned not breadth but **regioselectivity**,
at residues 437/478/479 on the cavity floor, while the hominoid CYP3A7 burst tuned **plasticity** at
the F/G roof and the membrane anchor. The promiscuity is not the adaptation. **The promiscuity is
the shadow of a large, deformable, membrane-facing steroid-and-lipid site, and the adaptation is
what happens on the floor beneath it.**

---

## 7. Consequences for this repo, and what would falsify the above

**Consequences for the CYP3A4 structure track:**

1. **Multiple occupancy is a live failure mode for single-ligand co-folding** (§4.4). OpenADMET has
   already said several of their ligands fit **multiple mutually exclusive conformations**. If the
   deposited density is multiply occupied or partially occupied and the prediction places one ligand
   centrally, LDDT-PLI penalises a physically defensible pose. Directly testable: does per-ligand
   LDDT-PLI correlate with ligand MW or with cavity-fill fraction?
2. **The floor is conserved and the roof is not** (§§1.3, 4.3, 4.6). This is a structural argument
   for why anchor-local features failed (memory: *anchor-local-features-fail*) and why the residual
   signal must live in substituent placement under the F/G roof — the part that is remodelled, the
   part OpenADMET named as the reason co-folding fails, and the part with the least conserved
   geometry.
3. **Apo structures are weak priors for holo** (§4.2). The 2004 literature disagreed about the apo
   cavity for a reason, and the disagreement was resolved by adding ligands, not by better crystals.
4. **CYP3A5 and CYP3A7 structures are usable negatives** (§4.3). CYP3A5's site is taller and
   narrower; CYP3A7's is rigid. Six plus two structures, 84–93% sequence identity to CYP3A4, and
   **measurably different plasticity**. That is a ready-made, leak-free test of whether a scorer is
   sensitive to lid geometry rather than to sequence identity.
5. **A truncated soluble construct is missing two of the six constraints** (§§4.5, 4.7). No membrane
   means no lateral entrance and no F′–G′ surface site; no CPR means no allosteric modulation from
   the reductase. Co-folding cannot represent either, which bounds how well it can ever do here.

**What would falsify the central inference of §6:**

- If a systematic assay of anionic lipid conjugates (prediction 1) or endocannabinoid-adjacent
  lipids (prediction 2) found CYP3A4 turnover at or below background, the "apolar core / polar rim /
  lateral entrance" model is wrong and the pocket is more water-facing than I claim.
- If two-ligand activation (prediction 5) cannot be demonstrated for any non-substrate/partner pair,
  then the 1-TST uncoupling result is about testosterone specifically and not about cavity volume,
  and constraint 2 of §6.1 collapses.
- If a well-powered selection scan localised the non-African CYP3A sweep to a *non-coding* element
  with no effect on any CYP3A coding sequence, the §1.3 "selection acted on regioselectivity" story
  would survive but the §5.3 phenotype debate would be settled against both salt and rickets.

---

## Appendix — sources, grouped

**Phylogeny and selection.** Qiu et al. 2008, *Pharmacogenet Genomics* 18:53–66,
doi:10.1097/fpc.0b013e3282f313f8, PMID 18216722 (the central paper). Yan & Cai 2010, *PLoS ONE*
5:e14276, doi:10.1371/journal.pone.0014276. *Xenobiotica* 50:1406–1412,
doi:10.1080/00498254.2020.1785580, PMID 32558606. Roca-Umbert et al. 2019, *BMC Evol Biol* 19:39,
doi:10.1186/s12862-019-1366-7, PMID 30704392. Thompson et al. 2004, *Am J Hum Genet* 75:1059–1069,
doi:10.1086/426406, PMID 15492926.

**Diet hypothesis.** Gonzalez & Nebert 1990, *Trends Genet* 6:182–186,
doi:10.1016/0168-9525(90)90174-5, PMID 2196721. *Nat Prod Bioprospect* 3:1–7,
doi:10.1007/s13659-013-0004-0. Johnson et al. 2018, *Nat Genet* 50:1102–1111,
doi:10.1038/s41588-018-0153-5, PMID 29967444. *Animals* 12:2821, doi:10.3390/ani12202821,
PMID 36290207. *Mol Phylogenet Evol* 217:108550, doi:10.1016/j.ympev.2026.108550, PMID 41580058.
*Genome Biol Evol* 13:evab261, doi:10.1093/gbe/evab261, PMID 34850870. *Insects* 16:184,
doi:10.3390/insects16020184. *Comp Biochem Physiol C* 280:109870, doi:10.1016/j.cbpc.2024.109870,
PMID 38428625.

**Endogenous substrates.** UniProt **P08684**, **P20815**, **P24462**, **Q9HB55**
(catalytic-activity blocks). *Steroids* 233–234:109819, doi:10.1016/j.steroids.2026.109819,
PMID 42251959 (4β-OHC). *FEBS J* 293:660–676, doi:10.1111/febs.70277, PMID 41046353
(CYP3A4 I301T / VDDR3). *Cureus* 15:e49976, doi:10.7759/cureus.49976. *Biomedicines* 10:2686,
doi:10.3390/biomedicines10112686, PMID 36359206 (CYP3A in health and disease). *J Endocr Soc*
9:bvaf119, doi:10.1210/jendso/bvaf119. *Front Immunol* 16:1692684,
doi:10.3389/fimmu.2025.1692684.

**Structure and plasticity.** Yano et al. 2004, doi:10.1074/jbc.C400293200, PMID 15258162 (1TQN).
Williams et al. 2004, *Science* 305:683–686, doi:10.1126/science.1099736, PMID 15256616
(1W0E/F/G). Ekroos & Sjögren 2006, *PNAS* 103:13682–13687, doi:10.1073/pnas.0603236103,
PMID 16954191 (2J0D/2V0M). Sevrioukova & Poulos 2017, *PNAS* 114:486–491,
doi:10.1073/pnas.1616198114, PMID 28031486 (5TE8). *JBC* 299:105117,
doi:10.1016/j.jbc.2023.105117, PMID 37524132 (caffeine, 8SO1/8SO2). *JBC* 299:104993,
doi:10.1016/j.jbc.2023.104993, PMID 37392852 (CYP3A7 + 4×DHEA-S, 8GK3). *IJMS* 22:5831,
doi:10.3390/ijms22115831, PMID 34072457 (7MK8, CYP3A7 rigidity). *Mol Pharmacol* 93:14–24,
doi:10.1124/mol.117.109744, PMID 29093019 (5VEU, CYP3A5 cavity). *IJMS* 23:12591,
doi:10.3390/ijms232012591, PMID 36293445 (8DYC, channel site). *J Inorg Biochem* 244:112211,
doi:10.1016/j.jinorgbio.2023.112211, PMID 37080138 (HDX-MS). *Commun Chem* 9:17,
doi:10.1038/s42004-025-01815-5, PMID 41513816 (F–F′ gate MD). *Arch Biochem Biophys* 758:110071,
doi:10.1016/j.abb.2024.110071, PMID 38909836 (cobicistat, 9BBB). *Chem Biol Drug Des* 105:e70043,
doi:10.1111/cbdd.70043, PMID 39792691. *DMD* 54:100210, doi:10.1016/j.dmd.2025.100210 (reductase as
allosteric modulator).

**Cooperativity.** *JBC* 282:7066–7076, doi:10.1074/jbc.m609589200, PMID 17213193 (1/2/3 TST).
*JBC* 286:5540–5545, doi:10.1074/jbc.M110.182055, PMID 21177853. *Arch Biochem Biophys*
488:146–152, doi:10.1016/j.abb.2009.06.013, PMID 19560436. *Arch Biochem Biophys* 519:91–102,
doi:10.1016/j.abb.2011.12.017, PMID 22245335.

**Regulation, zonation, PTM, population genetics.** *Pharmacogenet Genomics* 30:107–116,
doi:10.1097/fpc.0000000000000402, PMID 32301865 (4C/CRISPR, R1–R4). *Clin Transl Sci* 15:2720–2731,
doi:10.1111/cts.13398, PMID 36045613 (shared DRR). *Biochem Pharmacol* 253:118330,
doi:10.1016/j.bcp.2026.118330, PMID 42567488 (gp78 circadian). *Cancer Biol Ther* 11:549–551,
doi:10.4161/cbt.11.6.14834, PMID 21270532 (gp78/CHIP). *Nature* 653:1148–1157,
doi:10.1038/s41586-026-10377-y. *Nat Metab* 8:741–756, doi:10.1038/s42255-026-01459-2. *Annu Rev
Pathol* 21:185–212, doi:10.1146/annurev-pathmechdis-042624-091820. *CPT PSP* 13:1570–1581,
doi:10.1002/psp4.13192, PMID 38923249. **gnomAD v4** (rs776746, rs2740574, rs35599367, rs776744;
queried 2026-09-20). **Ensembl** GRCh38 gene coordinates (queried 2026-09-20). **RCSB** structure
inventory by UniProt accession (queried 2026-09-20).

**CYP3A43.** *Biomedicines* 10:3022, doi:10.3390/biomedicines10123022, PMID 36551778. *Xenobiotica*
52:413–425, doi:10.1080/00498254.2022.2078751, PMID 35582917. *J Neural Transm* 122:29–34,
doi:10.1007/s00702-014-1298-8, PMID 25150845.
