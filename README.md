# OpenADMET CYP3A4 — Structure Prediction Track

Entry for the structure-prediction track of the OpenADMET CYP blind challenge
(announced 2026-09-08; interim leaderboard **2026-09-24**).

Read `CLAUDE.md` first — it carries the measured facts this repo is built on, and the
storage constraints that will otherwise break the machine.

## The short version

OpenADMET says co-folding does poorly on CYP3A4. Their own analysis says *why*: the heme
and backbone are placed correctly, and it is the **ligand orientation** that is wrong,
including near-180° flips. Only 57.8% of Boltz-2 poses reach BiSyRMSD < 2 Å.

Parsing all 116 deposited CYP3A4 entries shows **83 of 101 unique ligands coordinate the
heme iron directly**, almost always through nitrogen, at 2.20 Å median and 171° from the
proximal thiolate. `cypstruct.chem.coordinating_atoms` predicts the right donor atom for
**82 of those 83**.

For the dominant binding mode, then, naming the donor atom determines the orientation.
That is the lever: steer the co-folder with a forced Fe→donor contact, and score poses
with a physics function that knows what a P450 sixth coordination site is.

## Quick start

```bash
pip install -e .
python -m cypstruct.targets                                   # sequences + RCSB inventory
python scripts/structure/build_reference_set.py               # measure all deposited CYP3A4
python scripts/structure/validate_donor_prediction.py         # check the steering mechanism
python scripts/cofold/modal_boltz.py plan --csv data/processed/validation_ligands.csv --tag v1
python scripts/ops/watchdog.py                                # spend, runaways, disk
```

`docs/RUNBOOK.md` says what to launch next and what each step must prove.
`docs/QM_SCORER_DESIGN.md` is the physics scorer.
`docs/FINETUNE_PLAN.md` is fine-tuning, and the case for doing it second.
