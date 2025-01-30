import argparse
from dbtk.nn import layers
from dbtk.nn.models import ContrastivePretrainingModel
from dnabert import DnaBertModel, DnaBertPretrainingModel
from dnadb import fasta
import lightning as L
import numpy as np
from pathlib import Path
import pandas as pd
import pickle
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Any, Callable, Optional, Tuple, Union

from sfd2.models import SampleEncoder, GcmsEncoder

def define_arguments(parser: argparse.ArgumentParser):
	group = parser.add_argument_group("DNABERT Model Config")
	group.add_argument("--dnabert_min_length", type=int, default=35)
	group.add_argument("--dnabert_max_length", type=int, default=250)
	group.add_argument("--dnabert_embed_dim", type=int, default=256)
	group.add_argument("--dnabert_num_heads", type=int, default=8)
	group.add_argument("--dnabert_stack", type=int, default=8)
	group.add_argument("--dnabert_kmer", type=int, default=3)
	group.add_argument("--dnabert_state_dict", type=Path, required=True)

	group = parser.add_argument_group("HTS Sample Model Config")
	group.add_argument("--hts_embed_dim", type=int, default=256)
	group.add_argument("--hts_num_heads", type=int, default=8)
	group.add_argument("--hts_stack", type=int, default=8)
	group.add_argument("--hts_sample_size", type=int, default=1000)
	group.add_argument("--hts_chunk_size", type=int, default=256)

	group = parser.add_argument_group("GCMS Sample Model Config")
	group.add_argument("--gcms_embed_dim", type=int, default=256)
	group.add_argument("--gcms_num_heads", type=int, default=8)
	group.add_argument("--gcms_stack", type=int, default=8)
	group.add_argument("--gcms_min_length")
	group.add_argument("--gcms_patch_size", type=int, default=20)

	group = parser.add_argument_group("Data")
	group.add_argument("--sequences", type=Path, required=True)
	group.add_argument("--sequence_mappings", type=Path, required=True)
	group.add_argument("--metadata", type=Path, required=True)
	group.add_argument("--leave_out_index", type=int, required=True)

	group = parser.add_argument_group("Training")
	group.add_argument("--max_steps", type=int, default=30000)
	group.add_argument("--strategy", type=str, default="auto")
	group.add_argument("--num_nodes", type=int, default=1)
	group.add_argument("--devices", type=int, default=1)


def load_dnabert_model(config: argparse.Namespace):
	dnabert = DnaBertPretrainingModel(
		DnaBertModel(
			layers.TransformerEncoder(
				layers.TransformerEncoderBlock(
					layers.RelativeMultiHeadAttention(
						embed_dim=config.dnabert_embed_dim,
						num_heads=config.dnabert_num_heads,
						max_length=config.dnabert_max_length,
					),
					feedforward_dim=config.dnabert_embed_dim,
				),
				num_layers=config.dnabert_stack
			),
			kmer=config.dnabert_kmer,
			kmer_stride=1
		)
	)
	dnabert.load_state_dict(torch.load(config.dnabert_state_dict))
	return dnabert.base


def create_model(config: argparse.Namespace):
	dna_encoder = load_dnabert_model(config)
	sample_encoder = SampleEncoder(
		dna_encoder,
		config.hts_embed_dim,
		num_heads=config.hts_num_heads,
		stack=config.hts_stack)
	gcms_encoder = GcmsEncoder(
		config.gcms_patch_size,
		config.gcms_embed_dim,
		num_heads=config.gcms_num_heads,
		stack=config.gcms_stack)
	return ContrastivePretrainingModel(
		sample_encoder,
		gcms_encoder,
		config.hts_embed_dim,
		config.gcms_embed_dim,
	)


def load_data(config: argparse.Namespace):
	sequences = fasta.FastaDb(config.sequences)
	mappings = sequences.mappings(config.sequence_mappings)
	with open("./data/gcms.pkl", "rb") as f:
		gcms_samples = pickle.load(f)
	gcms_groups = {}
	for f in Path("./data").glob("**/*.mzML"):
		group = f.name.split('-')[0].lower()
		gcms_groups[group] = f
	mappings_dict = {m.name.split('_')[0]: m for m in mappings}
	metadata = pd.read_csv(config.metadata)
	metadata = metadata[metadata["swab.origin"] == "snake"]
	hts_samples = {}
	for _, row in metadata.iterrows():
		group = row["exp.call"]
		if group not in gcms_samples or row["group"] not in mappings_dict:
			continue
		hts_samples[group] = mappings_dict[row["group"]]
	for k in set(gcms_samples.keys()) - set(hts_samples.keys()):
		del gcms_samples[k]
	return hts_samples, gcms_samples


def create_data_loader(config, hts_samples, gcms_samples, sample_names, model):
	tokenizer = model.encoder_a.dna_encoder.tokenizer
	vocabulary = model.encoder_a.dna_encoder.vocabulary
	sample_size = 1000
	min_gcms_length = 15000
	max_gcms_length = 20000
	rng = np.random.default_rng()

	def collate(sample_names):
		# Draw sequences
		def process_entry(entry):
			sequence = entry.sequence
			length = torch.randint(config.dnabert_min_length, min(len(sequence), config.dnabert_max_length) + 1, (1,)).item()
			offset = torch.randint(0, len(sequence) - length + 1, (1,)).item()
			token_ids = torch.full((config.dnabert_min_length - config.dnabert_kmer + 1,), vocabulary["[PAD]"])
			token_ids = torch.tensor(list(vocabulary(tokenizer(sequence[offset:offset+length].encode()))))
			return F.pad(token_ids, (0, (config.dnabert_min_length - config.dnabert_kmer + 1) - token_ids.shape[0]))
		sequences = torch.stack([torch.stack([process_entry(x) for x in hts_samples[name].sample(sample_size, rng)]) for name in sample_names])
		# Draw chromatograms
		retention_times = torch.zeros((len(sample_names), max_gcms_length))
		signals = torch.zeros((len(sample_names), max_gcms_length))
		mask = torch.ones_like(signals, dtype=torch.bool)
		for i, name in enumerate(sample_names):
			retention_time, signal = torch.from_numpy(gcms_samples[name]).T
			# Augment time by adding random shift noise
			dx = (retention_time[-1] - retention_time[0]) / len(retention_time)
			retention_offset = torch.normal(torch.tensor(0.0), dx)
			retention_time = retention_time + retention_offset
			# Augment signal by adding random shift noise
			...
			# Randomly trim
			trimmed_length = torch.randint(min_gcms_length, signal.size(0), (1,)).item()
			offset = torch.randint(0, signal.size(0) - trimmed_length, (1,)).item()
			retention_time = retention_time[offset:offset+trimmed_length]
			signal = signal[offset:offset+trimmed_length]
			signal = signal / torch.max(signal)
			# Insert into batch
			retention_times[i, :trimmed_length] = retention_time
			signals[i, :trimmed_length] = signal
			mask[i, :trimmed_length] = False

		return sequences, (retention_times, signals, mask)

	sampler = torch.utils.data.RandomSampler(sample_names, replacement=False, num_samples=len(sample_names))
	train_loader = torch.utils.data.DataLoader(sample_names, batch_size=len(sample_names), collate_fn=collate, sampler=sampler)
	return train_loader


def main(config: argparse.Namespace):
	# Data
	hts_samples, gcms_samples = load_data(config)
	sample_names = list(hts_samples.keys())
	sample_names.sort()
	leave_out_sample = sample_names[config.leave_out_index]
	sample_names.remove(leave_out_sample)
	print("Leaving out", leave_out_sample)

	# Model
	model = create_model(config)

	# Data Loader
	train_loader = create_data_loader(config, hts_samples, gcms_samples, sample_names, model)

	# Trainer
	trainer = L.Trainer(
		max_steps=config.max_steps,
		logger=L.pytorch.loggers.CSVLogger("logs", name=f"{leave_out_sample}"),
		log_every_n_steps=50,
		strategy=config.strategy,
		num_nodes=config.num_nodes,
		devices=config.devices
	)
	trainer.fit(model, train_loader)

	Path("./models").mkdir(exist_ok=True)
	torch.save(model, f"models/{leave_out_sample}.pt")


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	define_arguments(parser)
	config = parser.parse_args()
	main(config)
