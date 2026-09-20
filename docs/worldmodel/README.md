# The shared-embedding world model — plan, status, and kill criteria

**Started 2026-09-20.** Three asks, one object. Kept here so the chain of reasoning is
followable and so no agent's research is lost to a context window.

---

## The unifying idea

A JEPA (Joint Embedding Predictive Architecture) predicts *representations* of masked
parts rather than reconstructing inputs. Pointed at interaction data it becomes:

> given the embedding of entity A and a context, predict the embedding of the entity B
> that binds it.

That is a **generative prior over binders in embedding space**. It is the same object
ask 3 wants — "reason backwards from what CYP3A4 evolved to do, to what it should bind" —
approached from the data side instead of the causal side. Where the learned prior and the
causal derivation agree, we have a candidate. Where they disagree is the interesting part.

It is also a strict generalisation of proteochemometrics: PCM concatenates a drug
embedding with a target embedding and regresses affinity. A JEPA instead learns the
*map* between the two spaces, so it can be asked "what should bind here?" and not only
"do these two bind?".

---

## Honest constraints, recorded up front

| constraint | value | consequence |
|---|---|---|
| GPU | **1× H200**, contended (3 of 4 held by another user for 3–9 days) | no from-scratch pretraining; frozen foundation embeddings + light predictor head |
| Local disk | C: ~34 GB, D: ~20 GB | nothing bulk lands locally, ever |
| Cluster scratch | `/scratch` 298 TB | all corpora live here |
| Cold storage | OneDrive ~4.9 TB via `cypstruct.storage` | finished corpora, never the hot path |
| Interim deadline | **2026-09-24** (4 days) | this program contributes nothing to it; interim ships from the corrected pool |
| Final deadline | 2026-11-03 | realistic target for a proof of concept, not a trained world model |

**The part that cannot be crammed is the data.** Compute can be bought later; a corpus
assembled badly has to be reassembled. So reconnaissance and acquisition run first and
continuously, and modelling waits on them.

---

## The three asks

### Ask 1 — the world model
Two arms, deliberately opposed:
- **Implicit.** Give it everything — DNA, RNA, protein, ligand, every affinity readout —
  and see what structure it finds unprompted.
- **Explicit.** Train/fine-tune only on binding interfaces structurally similar to
  CYP3A4's, and see whether a narrow, well-chosen prior beats a broad one.

The comparison *is* the experiment. If broad beats narrow, biological priors are better
learned than imposed; if narrow wins, the opposite. Either answer is worth having.

### Ask 2 — causal decomposition of CYP3A4
Folds across apo/holo and isoforms; the phylogenetic tree and what each branch was
selected *for*; ancestral and current binding partners; PTMs as evidence of use-cases;
liver cell-type and transcriptional context; then reasoning forward from shape and
sequence to what *should* bind, including chemistry we have no assay for.

### Ask 3 — the dashboard
Real, interactive, continuously updated: the phylogenetic tree, the domains, the
structures, the embedding space as it forms. Merged into the Vercel dashboard.

---

## Kill criteria (pre-registered, per this repo's rules)

1. **Embedding transfer.** If frozen foundation embeddings cannot beat a Morgan-fingerprint
   + one-hot-target PCM baseline on a held-out affinity task, the shared space is not
   carrying biology and the program stops at the baseline.
2. **Partner prediction.** If the JEPA head cannot retrieve true partners above chance on
   held-out interactions (AUC > 0.6 on a difficulty-matched negative set), the predictive
   framing is wrong.
3. **Difficulty matching.** Per FINDING 022, any benchmark whose baseline is already
   >0.80 LDDT-PLI or >85% sub-2 Å is *uninformative* and must not be used to claim a gain.
4. **Cost.** If acquisition exceeds ~10 TB on `/scratch` before a single predictive result,
   stop and prune.

---

## Status

| wave | scope | state |
|---|---|---|
| 1 | model + data reconnaissance, 6 parallel agents | **launched 2026-09-20** |
| 2 | acquisition to `/scratch`, deduplication, splits | pending wave 1 |
| 3 | embedding extraction, JEPA head, baselines | pending wave 2 |
| 2b | CYP3A4 causal decomposition | **launched with wave 1** |

Findings land in `docs/worldmodel/*.md`, one file per domain, written by the agent that
did the work. This file is the index and the argument; those are the evidence.
