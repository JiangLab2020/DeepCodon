import torch

from DeepCodon.config.config import *
from DeepCodon.src.train.datasets import codon_vocab


def greedy_algorithm(model, device, src, opmask_inputs):
    model.eval()
    tgt = torch.tensor([[codon_vocab["<S>"]]]).to(device)

    for i in range(pclen + 1):
        out = model(src, tgt, opmask_inputs)
        predict = torch.argmax(out, dim=1)
        y = predict[-1].unsqueeze(0).unsqueeze(0)
        tgt = torch.concat([tgt, y], dim=1)
        if y == codon_vocab["<E>"]:
            break
    return tgt
