import torch
import torch.nn.functional as F

from DeepCodon.config.config import pclen
from DeepCodon.src.train.datasets import codon_vocab


def TopKSampler(inputs, k=6, fill=0):
    # Perform top-k operation on each row
    top_k_values, top_k_indices = torch.topk(inputs, k=2, dim=-1)
    new_probabilities = torch.full_like(inputs, fill)
    new_probabilities.scatter_(-1, top_k_indices, top_k_values)
    return new_probabilities


def top_k_algorithm(model, device, src):
    model.eval()
    # torch.manual_seed(42)
    tgt = torch.tensor([[codon_vocab["<S>"]]]).to(device)
    for i in range(pclen + 1):
        out = model(src, tgt)
        # softmax
        out = F.softmax(out, dim=-1)
        print(out.shape)
        out = TopKSampler(out)
        samples = torch.stack(
            [
                torch.multinomial(out[i], 1, replacement=True)
                for i in range(out.shape[0])
            ],
            dim=0,
        )
        y = samples[-1].unsqueeze(0)
        tgt = torch.concat([tgt, y], dim=1)
        if y == codon_vocab["<E>"]:
            break
    return tgt
