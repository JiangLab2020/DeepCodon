import torch
import torch.nn.functional as F

from DeepCodon.config.server_config import ϵ
from DeepCodon.src.train.datasets import pc_vocab

property = {
    "<pad>": 1,
    "TTT": 0.57,
    "TTC": 0.43,
    "TTA": 0.15,
    "TTG": 0.12,
    "TCT": 0.11,
    "TCC": 0.11,
    "TCA": 0.15,
    "TCG": 0.16,
    "TAT": 0.53,
    "TAC": 0.47,
    "TAA": 0.64,
    "TAG": 0.0,
    "TGT": 0.42,
    "TGC": 0.58,
    "TGA": 0.36,
    "TGG": 1,
    "CTT": 0.12,
    "CTC": 0.1,
    "CTA": 0.05,
    "CTG": 0.46,
    "CCC": 0.13,
    "CCA": 0.14,
    "CCG": 0.55,
    "CCT": 0.17,
    "CAC": 0.45,
    "CAT": 0.55,
    "CAA": 0.3,
    "CAG": 0.7,
    "CGC": 0.44,
    "CGA": 0.07,
    "CGG": 0.07,
    "CGT": 0.36,
    "ATT": 0.58,
    "ATC": 0.35,
    "ATA": 0.07,
    "ATG": 1,
    "ACC": 0.47,
    "ACA": 0.13,
    "ACG": 0.24,
    "ACT": 0.16,
    "AAC": 0.53,
    "AAT": 0.47,
    "AAA": 0.73,
    "AAG": 0.27,
    "AGC": 0.33,
    "AGT": 0.14,
    "AGA": 0.02,
    "AGG": 0.03,
    "GTT": 0.25,
    "GTC": 0.18,
    "GTA": 0.17,
    "GTG": 0.4,
    "GCC": 0.31,
    "GCA": 0.21,
    "GCG": 0.38,
    "GCT": 0.11,
    "GAC": 0.35,
    "GAT": 0.65,
    "GAA": 0.7,
    "GAG": 0.3,
    "GGC": 0.46,
    "GGA": 0.13,
    "GGG": 0.12,
    "GGT": 0.29,
    "<S>": 1,
    "<E>": 1,
}
reversed_property = {}

for protein_name, codon in pc_vocab.items():
    property_values = [property[c] for c in codon.keys() if c in property]
    min_val = min(property_values)
    max_val = max(property_values)

    for nowcodon, values in codon.items():
        if nowcodon in property:
            norm_value = (property[nowcodon]) / (max_val + ϵ)
            reversed_property[nowcodon] = norm_value
            # reversed_property[nowcodon] = 1 / (property[nowcodon] + ϵ)
        else:
            print(nowcodon)

reversed_property["<E>"] = 0
reversed_property["<S>"] = 0
reversed_property["<pad>"] = 0

# The generated result replaces the original property
# sorted_values = [reversed_property[k] for k in property.keys()]
sorted_values = [reversed_property.get(k, 1) for k in property]
print(sorted_values)


def get_logit_bias(input_tensor, λ):
    # use the zi + λ * pi

    sorted_values_tensor = torch.tensor(sorted_values).to(input_tensor.device)
    mask = input_tensor[-1] > 0

    sorted_values_tensor = torch.where(
        mask,
        sorted_values_tensor,
        torch.tensor(float("-inf"), device=input_tensor.device),
    )
    # add zi and λ * pi
    output = input_tensor.clone()

    output[-1] = F.softmax(output[-1] + λ * sorted_values_tensor, dim=-1)
    # NOTE FIND tensor dimension not match so just modify the last one
    # if λ > 0.4:
    #     print("fai", λ)
    #     print("input_tensor", input_tensor)
    #     print("output", output)
    return output
