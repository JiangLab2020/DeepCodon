import random

from torch.utils.data import Subset

from DeepCodon.config.config import *
from DeepCodon.src.train.datasets import *


def processing(outpath: str) -> str:
    enc_inputs_val = torch.load(enc_inputs_val_path)
    dec_inputs_val = torch.load(dec_inputs_val_path)
    dec_outputs_val = torch.load(dec_outputs_val_path)
    val_data = MyDataSet(enc_inputs_val, dec_inputs_val, dec_outputs_val)

    val_data = MyDataSet(enc_inputs_val, dec_inputs_val, dec_outputs_val)
    total_samples = len(val_data)
    indices = list(range(total_samples))
    random.seed(42)
    random_indices = random.sample(indices, 1000)  # NOTE 1000 samples
    val_data_subset = Subset(val_data, random_indices)

    val_loader = Data.DataLoader(val_data_subset, batch_size=1, shuffle=False)
    with open(outpath, "w", encoding="utf-8") as f:
        for number, (enc_inputs, dec_inputs, dec_outputs) in tqdm(
            enumerate(val_loader)
        ):
            predicted_words = [
                src_idx2word[n.item()] for n in enc_inputs.squeeze(0).squeeze()
            ]
            n = (
                predicted_words.index("<pad>")
                if "<pad>" in predicted_words
                else len(predicted_words)
            )
            predicted_words = [tgt_idx2word[i.item()] for i in dec_inputs.squeeze()][
                1 : n + 1
            ]
            predicted_words_str = "".join(predicted_words)
            f.write(f">{number}\n{predicted_words_str}\n")


if __name__ == "__main__":
    processing(
        PROJECT_ROOT
        / "data"
        / "datasets"
        / "test_data"
        / "compare"
        / "compareAll.fasta",
    )
