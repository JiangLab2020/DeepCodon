import torch
import torch.nn.functional as F

from DeepCodon.config.config import pclen
from DeepCodon.src.train.datasets import codon_vocab


def TopPSampler(inputs, p=0.8, fill=0):
    modified_inputs = inputs.clone()

    # top-p
    if len(inputs.shape) == 2:
        for row_idx, row in enumerate(inputs):
            sorted_probs, sorted_indices = torch.sort(row, descending=True)
            cumulative_probs = torch.cumsum(sorted_probs, dim=-1)
            # Find indexes that exceed the threshold
            topp_indices = (cumulative_probs >= p).nonzero(as_tuple=False)
            if len(topp_indices) > 0:
                topp_index = topp_indices[0][0]
                # Keep the values within the top-p range and replace the rest with 0
                modified_inputs[row_idx][row < sorted_probs[topp_index]] = fill

    elif len(inputs.shape) == 3:
        for sub_tensor_idx, sub_tensor in enumerate(inputs):
            for row_idx, row in enumerate(sub_tensor):
                sorted_probs, sorted_indices = torch.sort(row, descending=True)
                cumulative_probs = torch.cumsum(sorted_probs, dim=-1)
                topp_indices = (cumulative_probs >= p).nonzero(as_tuple=False)
                if len(topp_indices) > 0:
                    topp_index = topp_indices[0][0]
                    modified_inputs[sub_tensor_idx][row_idx][
                        row < sorted_probs[topp_index]
                    ] = fill

    return modified_inputs


def top_p_algorithm(model, device, src):
    model.eval()
    # torch.manual_seed(42)
    tgt = torch.tensor([[codon_vocab["<S>"]]]).to(device)
    for i in range(pclen + 1):
        # transformer
        out = model(src, tgt)
        # softmax
        out = F.softmax(out, dim=-1)
        out = TopPSampler(out)
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
