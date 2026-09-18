# FINDING 017 — the challenge test set is interpolation, not extrapolation

**Status: stands.** Measured 2026-09-18 on the released blinded test set (750 compounds)
against the full inhibition training set (4,905 compounds), ECFP4 / Morgan radius 2, 2048
bits, nearest-neighbour Tanimoto.

## The measurement

For each test compound, its similarity to the **nearest** training compound. Nearest
rather than mean, because what makes a test compound easy is having one close analogue in
train; a mean over 4,905 compounds washes that out and reports everything as equally novel.

| statistic | value |
|---|---|
| min | 0.318 |
| p5 | 0.500 |
| p25 | 0.543 |
| **median** | **0.587** |
| p75 | 0.639 |
| p95 | 0.741 |
| max | 0.952 |
| mean | 0.598 |

```
  0.30-0.35      1
  0.35-0.40      0
  0.40-0.45      2
  0.45-0.50     29  ###
  0.50-0.55    176  ######################
  0.55-0.60    212  ##########################
  0.60-0.65    182  ######################
  0.65-0.70     71  ########
  0.70-0.75     41  #####
  0.75-0.80     25  ###
  0.80-0.85      8  #
  0.85-0.90      2
  0.95-1.00      1
```

**Truly novel chemistry (NN < 0.4): 3 compounds, 0.4% of the test set.**
**Analog-dense (NN ≥ 0.7): 10.3%.**
**Everything else — 89% — sits in a narrow band between 0.45 and 0.70.**

## What it means

**This is an interpolation task.** Essentially every test compound has a reasonably close
relative in train. Nobody is being asked to extrapolate to new chemical space; they are
being asked to be accurate in the neighbourhood of what they were shown. That shapes what
is worth building:

- Methods that lean on **local similarity** — nearest-neighbour reasoning, analogue
  transfer, local models — are better matched to this task than their reputation suggests.
- The usual worry about **scaffold-novel generalisation does not apply here**, so effort
  spent hardening against it is effort misdirected.
- For the **structure track**: the test ligands are chemically close to ligands we already
  have folded and scored. Structural features learned on the training chemistry should
  transfer, which is a materially more optimistic setting than our own validation implies.

## Why this changes how to read our own numbers

Every generalisation figure in this repo — FINDING 012's +0.0357 over 81 proteins most of
all — comes from **leave-one-TARGET-out**: hold out an entire protein, select on the rest.
That is a legitimate test of cross-protein transfer, and it is **harder than the task the
challenge actually sets**, which holds out *compounds* inside a narrow similarity band on
proteins we have seen.

So the honest reading flips from cautious to slightly optimistic: our held-out number is
measured under a stricter regime than drop day will impose. It is **not** a licence to
expect more — the catastrophe-detector mechanism (FINDING 012) says the gain depends on
pool quality, not on chemical novelty — but it does remove a specific worry that was
implicitly priced in.

## A correction worth recording

My first read of the percentiles was that p5 landing on **exactly 0.5000** looked like a
constructed floor — the organisers filtering the test set to NN ≥ 0.5. It is not: the
minimum is 0.318 and 32 compounds sit below 0.5. Seventeen compounds happen to have NN
exactly 0.5, which is a common value for small-bit-count Tanimoto ratios, and the
percentile landed on that stack.

The instinct was right — a suspiciously round boundary usually *is* a filter, and this
repo has a standing note that too-clean numbers are the tell. The check cost one command
and it said otherwise, which is the point of running it rather than reporting the hunch.

## What it is used for

`scripts/structure/similarity_matched_splits.py` reproduces this distribution on the P450
ligand set, so a method can be validated across the same chemical gap the challenge
imposes rather than across whole proteins. That is steps 1–2 of the method-hopping ladder:
optimise on an analogous target under matched conditions, then hop.

---

## Does the selector care about chemical novelty? Mostly no, and it is strongest where
## the challenge lives

Each P450 ligand's nearest-neighbour Tanimoto to the *other* ligands in our set gives it a
novelty score on the same scale as the measurement above. Selector gain, stratified:

| ligand novelty (NN to rest) | pairs | random | selected | gain |
|---|---|---|---|---|
| most novel, NN < 0.45 | 134 | 0.6684 | 0.6911 | **+0.0227** |
| **0.45–0.60 — the challenge's band** | 82 | 0.7361 | 0.7788 | **+0.0427** |
| 0.60–0.80 | 163 | 0.7286 | 0.7686 | +0.0400 |
| least novel, NN ≥ 0.80 | 110 | 0.6399 | 0.6912 | **+0.0513** |
| all | 489 | 0.6934 | 0.7317 | +0.0383 |

**The number that matters: in the 0.45–0.60 band the challenge actually occupies, the
selector gains +0.0427 — above its own overall +0.0383.** The task's chemistry sits where
this method works slightly better than average, not worse.

**It never collapses.** Even on the most novel ligands it returns +0.0227, above the
+0.0138 noise floor. There is no regime in this data where the selector stops working.

### The confound, stated rather than hidden

Gain does trend upward with familiarity (+0.0227 → +0.0513), which invites the conclusion
that the selector "recognises" chemistry. **The random baselines make that unsafe.** The
least-novel band has the *lowest* baseline (0.6399) — the worst pools — and FINDING 012
says the catastrophe detector pays most exactly where pools are worst. So the apparent
novelty trend is at least partly pool quality wearing a different label, and these n do not
separate them.

What survives the confound is the useful part: the challenge band gains +0.0427 despite
having the **easiest** pools of any band (baseline 0.7361). That direction cannot be
explained by pool quality, since easy pools should depress the gain. It is the one cell in
the table where the two explanations disagree, and the answer favours the selector.

### Why this is what a catastrophe detector should look like

The mechanism keys on disagreement between independent predictions of the *same* molecule,
not on having seen that molecule's scaffold before. There is no learned chemistry in it to
degrade — no fitted parameters at all. A method that leaned on chemical familiarity would
fall off a cliff below NN 0.45; this one loses about half its gain and keeps working.

**For a fitted scorer the conclusion inverts**: anything with learned parameters *should*
be validated on the matched splits, because it has something to overfit that this selector
does not. That is what `p450_matched_splits.json` is for.
