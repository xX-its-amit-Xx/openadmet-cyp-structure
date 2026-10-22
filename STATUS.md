# STATUS — where this repo stands

Last updated: 2026-10-22 (Devin session bfc60e76b2d640c4966654e630406c50)

Start here. Read `CLAUDE.md` first, then `docs/README.md` (the findings index),
then this file's "Pending" section.

## What is done

- Selector shipped: `cypstruct.xengine.select()` — cross-engine agreement in the
  heme frame (FINDING 011/012). val87: selected **0.6164**, oracle **0.6975**,
  random **0.5769**. It beats random; native confidence scores do not.
- Drop-day pipeline ran end to end on the 20 blinded test ligands
  (`docs/DROP_RUN_2026-09-23.md`). Two candidate zips built:
  - `submissions/02_drop.zip` — ungated build (md5 `b8dc9d5c…`)
  - `submissions/03_drop_pbfilter.zip` — **recommended**: removes the 3 poses the
    organisers' PoseBusters gate would zero (validated +0.0650 on val87;
    pre-announcement estimate reverts to ~0.54–0.58).

## Pending — the actual next moves

1. **HUMAN ACTION — submit `submissions/03_drop_pbfilter.zip`** to the HF portal.
   Never submitted; the choice between 02 and 03 was deliberately left to the user.
   Caveats for that call are in `DROP_RUN_2026-09-23.md` §12.5.
2. **Check `scripts/ops/readiness.py`** for upstream changes (dataset growth,
   config flips) before any new run; re-read the Space config if it flips.
3. **T6 scoring-function handoff** is live — brief at
   `docs/HANDOFF_scoring_function.md`; poses/labels on HuggingFace
   (`xX-its-amit-Xx/cyp-cofold-poses`). Integrate results when they land.

## Open threads (not blocking submission)

- Final leaderboard **2026-11-03** — if a second drop is allowed, reference depth
  ≥ 8 is worth +0.0076 (FINDING 041); regenerate the drop reference at depth ≥ 8.
- `P20815` remains an unexplained selector failure (README "Superseded" section) —
  frozen until ≥20 proteins have known outcomes.
- `submit` cannot resubmit failed batches (FINDING 038): failed batches stay marked
  claimed — manual cleanup needed if a run partially fails.
- Secrets: `OPENPROTEIN_PASS` must come from env; no defaults (post-public-repo
  cleanup, commit `3509b56`).

## Resume checklist for the next agent

1. Read `CLAUDE.md` → `docs/README.md` → this file.
2. `git status` should be clean; everything through `3509b56` is pushed.
3. `python scripts/ops/readiness.py` → exit 2 = upstream changed, re-check the Space config.
4. The measurement contract: noise floor **+0.0138**, honest gating = beat the
   incumbent on leave-one-cluster-out; report pool oracle next to every selection
   number. See `CLAUDE.md` "Rules".
