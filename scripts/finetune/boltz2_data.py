"""A Boltz-2 training datamodule, because the released one is Boltz-1's.

This is the gap that blocked the fine-tune. `boltz/data/module/trainingv2.py` looks like
the Boltz-2 trainer and is not: it imports `BoltzFeaturizer` and the v1 `Structure`, and
`Boltz2Featurizer` is referenced by exactly one file in the whole package -
`inferencev2.py`. Feeding v1 features to `Boltz2.training_step` does not raise; it trains
on the wrong tensors. So the training path here is built from the *inference* v2 path,
which is the one the released weights were exported against, with `training=True` and a
cropper added.

Three things this does differently from boltz's own trainer, all because N is ~300 and
not ~300,000:

1. An epoch is the record list, once, in order. No `Sampler`, no `samples_per_epoch`.
   Cluster-balanced sampling matters when clusters are unbalanced across millions of
   entries; here it would just hide which structures were seen.
2. The crop is centred on the LIGAND chain, not on a random chain or interface. The
   thing being fine-tuned is ligand placement (FINDING 001: the backbone is already
   ~1.0 A; the error is in the ligand), so a crop that drops the ligand trains nothing
   we care about.
3. A failed item raises instead of silently falling back to item 0. Boltz's loaders
   return `self.__getitem__(0)` on any exception, which at this N would let one bad
   record become a large fraction of an epoch without a single log line.

Import from the trainer:

    from boltz2_data import Boltz2FinetuneDataModule
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytorch_lightning as pl
import torch
from torch.utils.data import DataLoader

from boltz.data import const
from boltz.data.feature.featurizerv2 import Boltz2Featurizer
from boltz.data.mol import load_canonicals, load_molecules
from boltz.data.pad import pad_to_max
from boltz.data.tokenize.boltz2 import Boltz2Tokenizer
from boltz.data.types import MSA, Input, Manifest, StructureV2


def collate(data: list[dict]) -> dict:
    """Boltz's collate, with the exempt-key list replaced by a type check.

    Boltz stacks every value and exempts six keys by NAME - the ragged symmetry and
    all-atom entries. The v2 featurizer emits at least one ragged key that is not on
    that list, so the name-based version dies in the worker with `'list' object has no
    attribute 'shape'`. Checking whether the values are actually Tensors is the same
    rule without the hardcoded list, and it cannot go stale against a new feature.
    """
    collated = {}
    for key in data[0]:
        values = [d[key] for d in data]
        if not all(isinstance(v, torch.Tensor) for v in values):
            collated[key] = values          # ragged: hand through as a list, as boltz does
            continue
        shape = values[0].shape
        if all(v.shape == shape for v in values):
            collated[key] = torch.stack(values, dim=0)
        else:
            collated[key], _ = pad_to_max(values, 0)
    return collated


def load_input(record, target_dir: Path, msa_dir: Path) -> Input:
    """StructureV2 plus whatever MSAs the record's chains name.

    Deliberately not imported from inferencev2: that one also takes constraints,
    templates and affinity paths we do not have, and its signature has changed between
    boltz releases. This is the whole of what we need and it is three lines.
    """
    structure = StructureV2.load(target_dir / "structures" / f"{record.id}.npz")
    msas = {}
    for chain in record.chains:
        msa_id = chain.msa_id
        if msa_id != -1 and msa_id != "":
            msas[chain.chain_id] = MSA(**np.load(msa_dir / f"{msa_id}.npz"))
    # `record` is not optional in practice: Boltz2Tokenizer reads record.affinity, so an
    # Input built without it fails at tokenize time with an AttributeError on None.
    return Input(structure, msas, record=record)


def pick_ligand_chain(structure) -> int | None:
    """The asym_id of the ligand to centre the crop on.

    HEM is excluded: every one of our targets has it, it is never the thing being
    predicted, and a crop centred on it is a crop centred on the middle of the protein.
    Among the remaining non-polymer chains we take the largest, which is the drug-like
    ligand rather than a buffer ion or cryoprotectant.
    """
    nonpoly = const.chain_type_ids["NONPOLYMER"]
    best, best_atoms = None, -1
    for i, ch in enumerate(structure.chains):
        if int(ch["mol_type"]) != nonpoly:
            continue
        name = str(ch["name"])
        res = structure.residues[int(ch["res_idx"]): int(ch["res_idx"]) + int(ch["res_num"])]
        res_names = {str(r["name"]) for r in res}
        if res_names & {"HEM", "HEC", "HEB", "SRM"}:
            continue
        if int(ch["atom_num"]) > best_atoms:
            best, best_atoms = i, int(ch["atom_num"])
        _ = name
    return best


class Boltz2FinetuneDataset(torch.utils.data.Dataset):
    """One item per record, tokenized and featurized exactly as inference does."""

    def __init__(
        self,
        manifest: Manifest,
        target_dir: Path,
        msa_dir: Path,
        mol_dir: Path,
        max_tokens: int = 512,
        max_atoms: int = 4096,
        max_seqs: int = 1024,
        training: bool = True,
        seed: int = 0,
    ) -> None:
        super().__init__()
        self.records = manifest.records
        self.target_dir = Path(target_dir)
        self.msa_dir = Path(msa_dir)
        self.mol_dir = Path(mol_dir)
        self.max_tokens = max_tokens
        self.max_atoms = max_atoms
        self.max_seqs = max_seqs
        self.training = training
        self.seed = seed
        self.tokenizer = Boltz2Tokenizer()
        self.featurizer = Boltz2Featurizer()
        self.canonicals = load_canonicals(self.mol_dir)
        from boltz.data.crop.boltz import BoltzCropper
        self.cropper = BoltzCropper()

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, idx: int) -> dict:
        record = self.records[idx]
        # Validation must be identical every epoch or the curve measures the sampler.
        # Training draws from torch's seeded global RNG, so it varies per epoch and per
        # worker while staying reproducible from `seed_everything`.
        if self.training:
            seed = int(torch.randint(0, 2**31 - 1, (1,)).item())
        else:
            seed = self.seed + idx
        # TWO generators on one seed, because boltz's two consumers want different APIs:
        # the cropper calls `random.randint(n)` (legacy RandomState) and the featurizer
        # expects a Generator. Handing a Generator to the cropper raises AttributeError;
        # handing RandomState to the featurizer silently takes different code paths.
        crop_rng = np.random.RandomState(seed)  # noqa: NPY002
        feat_rng = np.random.default_rng(seed)

        input_data = load_input(record, self.target_dir, self.msa_dir)
        tokenized = self.tokenizer.tokenize(input_data)

        chain_id = pick_ligand_chain(input_data.structure)
        tokenized = self.cropper.crop(
            tokenized,
            max_atoms=self.max_atoms,
            max_tokens=self.max_tokens,
            random=crop_rng,
            chain_id=chain_id,
        )
        if len(tokenized.tokens) == 0:
            msg = f"{record.id}: crop produced no tokens"
            raise ValueError(msg)

        molecules = dict(self.canonicals)
        molecules.update(input_data.extra_mols or {})
        names = set(tokenized.tokens["res_name"].tolist()) - set(molecules)
        molecules.update(load_molecules(self.mol_dir, names))

        feats = self.featurizer.process(
            tokenized,
            molecules=molecules,
            random=feat_rng,
            training=self.training,
            max_atoms=self.max_atoms,
            max_tokens=self.max_tokens,
            max_seqs=self.max_seqs,
            pad_to_max_seqs=True,
            # msa_sampling only during training: at validation the same subsample every
            # time is what makes two epochs comparable.
            msa_sampling=self.training,
            single_sequence_prop=0.0,
            compute_frames=True,
            compute_constraint_features=True,
            # Symmetry-corrected ground truth. Boltz2.training_step routes the true
            # coordinates through minimum_lddt_symmetry_coords, which reads
            # chain_symmetries / amino_acids_symmetries / ligand_symmetries out of the
            # batch; without them every batch raises and training_step returns None, so
            # the run completes having learned nothing. It is also the right choice on
            # the merits: a symmetric ligand matched atom-for-atom is the exact trap that
            # halves LDDT-PLI elsewhere in this repo.
            compute_symmetries=True,
        )
        # training_step's own error handler prints batch['pdb_id']; without it a real
        # error is replaced by a KeyError and the actual cause is never printed.
        feats["pdb_id"] = torch.tensor([idx], dtype=torch.long)
        return feats


class Boltz2FinetuneDataModule(pl.LightningDataModule):
    """Train/validate over one arm's manifests."""

    def __init__(
        self,
        arm_dir: Path,
        msa_dir: Path,
        mol_dir: Path,
        batch_size: int = 1,
        num_workers: int = 2,
        max_tokens: int = 512,
        max_atoms: int = 4096,
        max_seqs: int = 1024,
    ) -> None:
        super().__init__()
        self.arm_dir = Path(arm_dir)
        self.msa_dir = Path(msa_dir)
        self.mol_dir = Path(mol_dir)
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.kw = dict(max_tokens=max_tokens, max_atoms=max_atoms, max_seqs=max_seqs)

    def setup(self, stage: str | None = None) -> None:
        self.train_set = Boltz2FinetuneDataset(
            Manifest.load(self.arm_dir / "manifest_train.json"),
            self.arm_dir, self.msa_dir, self.mol_dir, training=True, **self.kw)
        self.val_set = Boltz2FinetuneDataset(
            Manifest.load(self.arm_dir / "manifest_test.json"),
            self.arm_dir, self.msa_dir, self.mol_dir, training=False, **self.kw)

    def train_dataloader(self) -> DataLoader:
        return DataLoader(self.train_set, batch_size=self.batch_size, shuffle=True,
                          num_workers=self.num_workers, collate_fn=collate,
                          pin_memory=False, persistent_workers=self.num_workers > 0)

    def val_dataloader(self) -> DataLoader:
        return DataLoader(self.val_set, batch_size=self.batch_size, shuffle=False,
                          num_workers=self.num_workers, collate_fn=collate,
                          pin_memory=False, persistent_workers=self.num_workers > 0)
