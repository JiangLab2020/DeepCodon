import torch
import torch.nn.functional as F

from DeepCodon.config.config import pclen
from DeepCodon.src.train.datasets import codon_vocab


def TemperatureSampler(inputs, temperature=0.7):
    # print("TemperatureSampler", temperature)
    out = F.softmax(inputs / temperature, dim=-1)
    return out


def temperature_algorithm(model, device, src):
    model.eval()
    # torch.manual_seed(42)
    tgt = torch.tensor([[codon_vocab["<S>"]]]).to(device)
    for i in range(pclen + 1):
        out = model(src, tgt)
        out = TemperatureSampler(out)
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
