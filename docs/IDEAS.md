# Idea log — cross-domain sweep for the CYP structure track

A standing list, swept daily. The rule that makes it useful: **one priority at a time gets
all the resources**, and an idea only leaves this file when it is either the active
priority or measured and closed. Ideas are cheap; half-finished workstreams are not.

Each entry carries a **kill criterion** written *before* the work starts, because the
noise floor here is +0.0138 (FINDING 007) and "promising" is a meaningless word against it.

Status: `open` · `active` · `measured-dead` · `shipped`

---

## Standing priority

**T5 — structure generation** (BU task, our critical path too). Everything else waits.

---

## A. Seeded by Amit, 2026-09-16

| # | idea | why it might work | kill criterion | status |
|---|---|---|---|---|
| A1 | **Model the protein in pieces, stitch back** | The F/G loop is the plastic region and the reason co-folding fails here (OpenADMET's own claim). Folding the rigid core and the loop separately, then splicing, decouples two problems the model currently solves jointly and badly. | Spliced structures must beat whole-protein co-folding on LDDT-PLI by > +0.020 on the 87-ligand set, and produce chemically valid junctions. | open |
| A2 | **Treat it as a fusion protein: two domains + linker** | Same decoupling as A1 but uses machinery the co-folders already have (they handle multi-domain proteins and linkers natively), so no custom splicing code. | Must beat A1's simpler splice, or A1 wins on cost. | open |
| A3 | **Treat the ligand as a very small protein; use PPI models** | PPI models are trained on a different objective (interface complementarity) than ligand docking (pose RMSD). A Type II CYP inhibitor coordinating Fe *is* an interface problem. | Needs a ligand→pseudo-peptide encoding that preserves the coordinating nitrogen's geometry. Dead if the encoding can't place a donor atom to 0.5 Å. | open |
| A4 | **Treat the complex as a long DNA/sequence object** | Sequence models scale differently and might capture the pocket as long-range context. Speculative. | Dead unless a sequence-only baseline reaches LDDT-PLI > 0.40, i.e. non-trivial. | open |
| A5 | **Fine-tune on synthetic conformers from the model itself** | AF3-style recycling: the model's own high-confidence poses become training data. We have 16.5k poses and ground truth for many. | Must beat the un-tuned engine on held-out TARGETS, not held-out ligands. Prior: PXR fine-tuning returned −0.0020 over 12 days. | open |
| A6 | **Fine-tune on fine-tuned poses (iterated recycling)** | Second-order version of A5. Only worth it if A5 shows any gain at all. | Gated on A5. | open |
| A7 | **Med-chem priors as literal coordinate edits** | The core bet: encode what a medicinal chemist knows (a triazole nitrogen points at the iron; halogens avoid the porphyrin plane) as an explicit geometric correction to predicted coordinates. This is the "biological prior baked in" idea. | Must improve LDDT-PLI post-hoc on poses it edits, and not degrade the ones it leaves alone. | open |
| A8 | **Custom structural scoring function** | Already the repo's thesis (FINDING 011 shipped one). Extends to a CYP-specific physics term. | Tier-0/1/2 measured; tier-1 and tier-2 are dead (FINDING 014). Cross-engine agreement shipped at +0.0381. | partly shipped |

## B. Cross-domain sweep — 2026-09-16

The brief: raid finance, health, cyber and art for method shapes, not metaphors.

| # | source domain | transplanted idea | why it is not just a metaphor | kill criterion | status |
|---|---|---|---|---|---|
| B1 | **Finance — portfolio variance** | Select on the agreement-to-variance ratio rather than agreement alone. | Must beat −z(xeng) by > +0.020 held out. | **measured-dead 2026-09-16** |
| B2 | **Finance — market-making / adverse selection** | Poses where the reference engines *disagree most* are where selection is most valuable. Spend more sampling budget only on high-disagreement ligands instead of uniformly. | Gain tracks pool catastrophe rate (FINDING 012), so budget should follow disagreement. Makes drop-day compute allocation adaptive. | Must beat uniform allocation at matched total job count. | open |
| B3 | **Cyber — differential fuzzing** | Two implementations that should agree, don't → the disagreement localises the bug. Feed the *same* ligand through engines with deliberately perturbed inputs (protonation, tautomer, chirality flip) and treat disagreement as a pose-quality signal. | This is cross-engine agreement generalised from "different engine" to "different input". Input perturbation is free; engine diversity is not. | Must produce a selector beating +0.0381. | open |
| B4 | **Cyber — canary / honeypot** | Insert known-answer decoy ligands into every batch. If the engine places a decoy wrongly, distrust that whole batch. | We have 87 crystal ground truths; salting them into production batches gives per-batch quality telemetry for free. Would have caught the determinism collapse days earlier. | Ships if it detects a seeded regression; it is monitoring, not a scorer. | open |
| B5 | **Art — pentimenti / underpainting** | Painters leave earlier versions beneath. Diffusion co-folders discard intermediate denoising states. The *trajectory* to a pose may say more about its reliability than the endpoint. | `num_steps` changes the trajectory and we showed it changes the answer (FINDING 016). Intermediate states may be retrievable. | Dead if OpenProtein exposes no intermediates — check first, costs one API read. | open |
| B6 | **Art — figure drawing / gesture-first** | Artists block in the gesture before detail. Predict the ligand's *axis and orientation* first, then its substituents, rather than all atoms at once. | FINDING 001: the failure is orientation, not location. A two-stage predictor matches the observed failure mode exactly. | Must beat one-stage on the orientation component specifically. | open |
| B7 | **Health — triage / differential diagnosis** | Don't score all poses equally; first classify the *failure mode* (wrong face, wrong rotamer, right pose), then apply a mode-specific scorer. | The catastrophe detector already implies a bimodal population. A mode classifier could route poses to different scorers. | Must beat the single scorer at matched information. | open |
| B8 | **Health — survival analysis / censoring** | Failed folds are *censored observations*, not missing data. T5's `completed` flag makes this explicit. Modelling which complexes fail is itself predictive of pocket flexibility. | Completion rate is a T5 deliverable already; treating it as a signal rather than a nuisance is free. | Must predict failure better than chance from ligand properties alone. | open |
| B9 | **Finance — backtest overfitting / deflated Sharpe** | With ~30 features tested against a +0.0138 noise floor, some "winner" is guaranteed by chance. Deflate the reported gain by the number of trials, as quant finance does for Sharpe ratios. | We have exactly the multiple-comparisons problem that field formalised, and FINDING 007's null is the raw material. | Not a gain — a correction. Ships as a reporting change. | open |

---

## Method-hopping validation ladder (Amit's spec, 2026-09-16)

The protocol every idea above must pass, in order:

1. **Optimise on an analogous target** with real data and matched train/test chemical
   similarity to the real CYP split.
2. **Hop to a second target** — if the method dies, it was target-specific.
3. **Hop to CYP** and read the post-hoc interpretation.

Status: the P450 superfamily set (FINDING 012, 87 proteins / 491 pairs / crystal ground
truth) is step 1 and 2 already built. What is **missing** is the chemical-similarity
matching: our splits are leave-one-TARGET-out, not similarity-matched to the challenge's
train/test split. See `docs/TODO_similarity_matched_splits.md`.


---

## Measured and closed

### B1 — portfolio variance (finance). **DEAD, 2026-09-16.**

489 pairs / 87 proteins / 6,314 poses, frozen sweep reference, zero new jobs.

| variant | selected | gain |
|---|---|---|
| **mean — the incumbent** | 0.7317 | **+0.0343** |
| median | 0.7325 | +0.0351 |
| mean − 0.5·std (risk-*seeking*) | 0.7321 | +0.0347 |
| mean + 0.5·std | 0.7319 | +0.0345 |
| mean / std (Sharpe-like) | 0.7308 | +0.0334 |
| mean + 1.0·std | 0.7299 | +0.0325 |
| max (worst-case / minimax) | 0.7281 | +0.0308 |
| **std alone** | 0.6291 | **−0.0682** |

The criterion was +0.0543 (incumbent +0.020). The best variant reaches +0.0351, beating
the incumbent by **+0.0008** — two orders below the noise floor.

**Why it fails, which is the part worth keeping.** `std alone` is strongly *negative*:
spread across reference poses is not a quality signal here at all, so there is no risk
dimension for a risk-adjusted score to trade against. Every variant that mixes in std
lands within ±0.002 of the plain mean, which is what you see when the added term carries
no information rather than the wrong information.

Two small orderings, both inside the noise and stated only as direction: median ≥ mean ≥
max. Robustness to a single outlier reference may help a little; committing to the
worst case hurts a little. Neither is worth a change.

The analogy failed at the level of the mapping, not the arithmetic. In a portfolio,
variance is *risk borne by the holder*; here, variance across reference poses is just
disagreement among opinions, and the catastrophe detector already extracts what that
disagreement is worth through the mean.
