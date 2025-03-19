import numpy as np
import torch.nn as nn
from torch.nn import Transformer

from DeepCodon.config.config import d_model, n_heads, n_layers
from DeepCodon.src.train.datasets import *

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class PositionalEncoding(nn.Module):
    def __init__(self, d_model, dropout=0.1, max_len=5000):
        super(PositionalEncoding, self).__init__()
        self.dropout = nn.Dropout(p=dropout)
        pos_table = np.array(
            [
                [pos / np.power(10000, 2 * i / d_model) for i in range(d_model)]
                if pos != 0
                else np.zeros(d_model)
                for pos in range(max_len)
            ]
        )
        pos_table[1:, 0::2] = np.sin(pos_table[1:, 0::2])
        pos_table[1:, 1::2] = np.cos(pos_table[1:, 1::2])
        self.pos_table = torch.FloatTensor(
            pos_table
        ).cuda()  # enc_inputs: [seq_len, d_model]

    def forward(self, enc_inputs):  # enc_inputs: [batch_size, seq_len, d_model]
        enc_inputs += self.pos_table[: enc_inputs.size(1), :]
        return self.dropout(enc_inputs.cuda())


class TransformerTranslator(nn.Module):
    def __init__(self):
        super(TransformerTranslator, self).__init__()
        self.device = device
        self.src_emb = nn.Embedding(src_vocab_size, d_model)
        self.tgt_emb = nn.Embedding(tgt_vocab_size, d_model)
        self.transformer = Transformer(
            d_model=d_model,
            nhead=n_heads,
            num_encoder_layers=n_layers,
            num_decoder_layers=n_layers,
            batch_first=True,
        )
        self.projection = nn.Linear(d_model, tgt_vocab_size, bias=True)
        self.pos_emb = PositionalEncoding(d_model)

    def forward(self, src, tgt, opmask=None):
        src_padding_mask = src == 0
        tgt_padding_mask = tgt == 0
        # src_mask = nn.Transformer.generate_square_subsequent_mask(src.size(1)).bool().to(self.device)
        tgt_mask = (
            nn.Transformer.generate_square_subsequent_mask(tgt.size(1))
            .bool()
            .to(self.device)
        )
        src = self.src_emb(src)
        src = self.pos_emb(src)
        tgt = self.tgt_emb(tgt)
        tgt = self.pos_emb(tgt)
        output = self.transformer(
            src,
            tgt,
            tgt_mask=tgt_mask,
            src_key_padding_mask=src_padding_mask,
            tgt_key_padding_mask=tgt_padding_mask,
        )
        if opmask is not None:
            output = self.projection(output) + opmask[:, : output.size(1)].log()
        else:
            output = self.projection(output)
        return output.view(-1, output.shape[-1])
