import torch
import torch.nn.functional as F

from DeepCodon.config.config import pclen
from DeepCodon.config.server_config import λTemplate
from DeepCodon.src.generate.temperature import TemperatureSampler
from DeepCodon.src.generate.top_p import TopPSampler
from DeepCodon.src.services.preBuildLogitBias import get_logit_bias
from DeepCodon.src.train.datasets import codon_vocab


def MixtureSampler(inputs, logit_bias, i):
    if logit_bias is not None:
        if i >= 0:
            λ = logit_bias[i]
            if λ > 0:
                inputs = get_logit_bias(inputs, λ * λTemplate)
    # Perform top-p operation on each row
    toppout = TopPSampler(inputs, p=0.8, fill=float("-inf"))
    # Perform temperature operation on each row
    temperatureout = TemperatureSampler(toppout, temperature=0.4)
    # print("\n\n\n-------------------------------")
    # print(topkout[-1])
    # print(toppout[-1])
    # print(temperatureout[-1])
    # print("-------------------------------\n\n\n")
    return temperatureout


def mixture_algorithm(
    model, device, src, opmask_inputs, debugFileOut=None, logit_bias=None
):
    model.eval()
    if debugFileOut is not None:
        f = open(debugFileOut, "w")
    # torch.manual_seed(42)
    tgt = torch.tensor([[codon_vocab["<S>"]]]).to(device)
    for i in range(pclen + 1):
        # Perform transformer calculations
        out = model(src, tgt, opmask_inputs)
        # Create a softmax
        out = F.softmax(out, dim=-1)
        if debugFileOut is not None:
            f.write(str(out.tolist()) + "\n")
        out = MixtureSampler(out, logit_bias, i=i - 1)
        # Generate random integer index
        # Perform random sampling
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
    if debugFileOut is not None:
        f.close()
    return tgt
