import torch
import torch.nn.functional as F

from DeepCodon.config.config import *
from DeepCodon.src.train.datasets import codon_vocab


def beam_search(model, device, src, opmask_inputs):
    model.eval()
    # Initialize the starting token for beam search
    start_token = torch.tensor([[codon_vocab["<S>"]]]).to(device)
    # Initialize the beam collection for beam search
    beams = [(start_token, 0.0)]

    for i in range(pclen + 1):
        new_beams = []

        for beam in beams:
            tgt = beam[0]
            # precessing
            out = model(src, tgt, opmask_inputs)
            log_probs = F.softmax(out, dim=1)
            top_probs, top_idx = torch.topk(log_probs[-1], beam_width)

            for j in range(beam_width):
                y = top_idx[j].unsqueeze(0).unsqueeze(0)
                prob = top_probs[j].item()

                # concat
                new_tgt = torch.cat([tgt, y], dim=1)

                # if <eos>，end
                if y == codon_vocab["<E>"]:
                    continue
                new_beams.append((new_tgt, beam[1] + prob))

        # Select the top beam-width beam with the highest probability
        new_beams = sorted(new_beams, key=lambda x: x[1], reverse=True)[:beam_width]
        beams = new_beams

    # TODO
    with open("data/user_data/output/beam.txt", "w") as f:
        for beam in beams:
            f.write(str(beam[0].tolist()) + "\n")
    return beams[0][0]
