from dbtk.nn import layers
import lightning as L
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.checkpoint import checkpoint
from typing import Optional, Tuple

class SampleEncoder(L.LightningModule):
    def __init__(self, dna_encoder: nn.Module, embed_dim: int, stack: int, num_heads: int, chunk_size: int = 256):
        super().__init__()
        self.dna_encoder = dna_encoder
        self.embed_dim = embed_dim
        self.stack = stack
        self.num_heads = num_heads
        self.chunk_size = chunk_size
        self.transformers = layers.TransformerEncoder(
            layers.TransformerEncoderBlock(
                layers.MultiHeadAttention(
                    embed_dim=self.embed_dim,
                    num_heads=self.num_heads
                ),
                feedforward_dim=self.embed_dim
            ),
            num_layers=self.stack
        )
        self.class_token = nn.Parameter(torch.randn(self.embed_dim))

    def forward(self, inputs, attention_head_mask = None, return_attention_weights = False):
        if not torch.is_floating_point(inputs):
            sequence_embeddings = []
            inputs_flat = inputs.flatten(0, -2)
            for i in range(0, inputs_flat.shape[0], self.chunk_size):
                sequence_embeddings.append(
                    checkpoint(
                        self.dna_encoder,
                        inputs_flat[i:i+self.chunk_size],
                        use_reentrant=False
                    )[0]
                )
            sequence_embeddings = torch.concat(sequence_embeddings, 0).unflatten(0, inputs.shape[:-1])
        else:
            sequence_embeddings = inputs
        class_tokens = self.class_token.expand(*inputs.shape[:-2], 1, -1)
        out = torch.cat((class_tokens, sequence_embeddings), -2)
        out, attention_weights = self.transformers(out, attention_head_mask=attention_head_mask, return_attention_weights=True)

        class_tokens = out.select(-2, 0)
        sequence_embeddings = out.narrow(-2, 1, out.shape[-2] - 1)

        if return_attention_weights:
            return class_tokens, attention_weights
        return class_tokens

class GcmsEncoder(L.LightningModule):
    def __init__(self, patch_size: int, embed_dim: int, stack: int, num_heads: int):
        super().__init__()
        self.patch_size = patch_size
        self.embed_dim = embed_dim
        self.stack = stack
        self.num_heads = num_heads

        self.class_token = nn.Parameter(torch.randn(size=(1, 1, self.embed_dim)))

        self.retention_conv = nn.Conv1d(
            self.embed_dim,
            self.embed_dim,
            kernel_size=self.patch_size,
            stride=self.patch_size)
        self.signal_conv = nn.Conv1d(
            1,
            self.embed_dim,
            kernel_size=self.patch_size,
            stride=self.patch_size)

        self.transformers = layers.TransformerEncoder(
            layers.TransformerEncoderBlock(
                layers.MultiHeadAttention(
                    embed_dim=self.embed_dim,
                    num_heads=self.num_heads
                ),
                feedforward_dim=embed_dim
            ),
            num_layers=self.stack
        )

    def _retention_time_encoding(self, position):
        if self.embed_dim % 2 != 0:
            raise ValueError("Cannot use sin/cos positional encoding with "
                             "odd dim (got dim={:d})".format(self.embed_dim))
        length = position.size(-1)
        pe = torch.zeros(position.size(0), length, self.embed_dim).to(position.device)
        position = position.unsqueeze(-1)
        div_term = torch.exp((torch.arange(0, self.embed_dim, 2, dtype=torch.float) *
                             -(np.log(10000.0) / self.embed_dim)))
        div_term = div_term.to(position.device)
        pe[:, :, 0::2] = torch.sin(position * div_term)
        pe[:, :, 1::2] = torch.cos(position * div_term)
        pe[torch.where((position == 0.0).squeeze(-1))] = 0
        return pe

    def forward(
        self,
        inputs: Tuple[torch.FloatTensor, torch.FloatTensor, Optional[torch.BoolTensor]],
        attention_head_mask = None,
        return_attention_weights = False
    ):
        retention_times, signals, mask = inputs
        retention_times = self._retention_time_encoding(retention_times)
        retention_times = self.retention_conv(retention_times.transpose(-1, -2)).transpose(-1, -2)
        signals = self.signal_conv(signals.unsqueeze(-1).transpose(-1, -2)).transpose(-1, -2)
        out = retention_times + signals

        mask = mask.narrow(-1, 0, mask.shape[-1] // self.patch_size * self.patch_size)
        mask = F.pad(torch.any(mask.unflatten(-1, (-1, 20)), -1), (1, 0))

        # Prepend class tokens
        class_tokens = self.class_token.expand(*out.shape[:-2], 1, -1)
        out = torch.cat([class_tokens, out], dim=1)
        out, attention_weights = self.transformers(
            out,
            src_key_padding_mask=mask,
            attention_head_mask=attention_head_mask,
            return_attention_weights=True)

        class_tokens, tokens = out[:,0,:], out[:,1:,:]
        if return_attention_weights:
            return class_tokens, attention_weights
        return class_tokens


class SampleGcmsContrastiveTrainer(L.LightningModule):
    def __init__(self, sample_encoder: nn.Module, gcms_encoder: nn.Module, max_temp: float = 100.0):
        super().__init__()
        self.sample_encoder = sample_encoder
        self.gcms_encoder = gcms_encoder
        self.max_temp = max_temp
        self.w_a = nn.Linear(sample_encoder.embed_dim, sample_encoder.embed_dim, bias=False)
        self.w_b = self.w_a # shared projections
        # self.w_b = nn.Linear(sample_encoder.embed_dim, sample_encoder.embed_dim, bias=False)
        # self.t = nn.Parameter(torch.tensor(1.0))
        self.t = torch.tensor(0.0)

    def forward(self, x, y):
        a_f = self.sample_encoder(x)
        b_f = self.gcms_encoder(*y)
        a_e = F.normalize(a_f, p=2, dim=1)
        b_e = F.normalize(b_f, p=2, dim=1)
        logits = torch.tensordot(a_e, b_e.T, 1)
        return a_f, b_f, logits

    def _step(self, batch):
        x, y = batch
        labels = torch.arange(x.size(0)).to(self.device)
        a_f, b_f, logits = self(x, y)
        loss_a = F.cross_entropy(logits * torch.exp(self.t), labels)
        loss_b = F.cross_entropy(logits.T * torch.exp(self.t), labels)
        loss = (loss_a + loss_b) / 2
        loss = loss + F.mse_loss(a_f, b_f)
         accuracy = torch.sum(
            (torch.argmax(logits, dim=-1) == labels).int() + (torch.argmax(logits, dim=-2) == labels).int()
        ) / labels.shape[-1] / 2.0
        return loss, accuracy

    def training_step(self, batch, batch_index):
        loss, accuracy = self._step(batch)
        self.log("train/loss", loss, on_step=True, on_epoch=False)
        self.log("train/accuracy", accuracy, on_step=True, on_epoch=False)
        return loss

    def validation_step(self, batch, batch_index):
        loss, accuracy = self._step(batch)
        self.log("validation/loss", loss, on_step=True, on_epoch=False)
        self.log("validation/accuracy", accuracy, on_step=True, on_epoch=False)
        return loss

    def on_train_batch_end(self, *args, **kwargs):
        self.t.data = torch.clamp(self.t, 0.0, self.max_temp)

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=1e-4)
        return optimizer