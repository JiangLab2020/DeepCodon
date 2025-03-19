import torch
import torch.utils.data as Data
from tqdm import tqdm

from DeepCodon.config.config import minlen, pclen, raw_output_file


class MyDataSet(Data.Dataset):
    def __init__(self, enc_inputs, dec_inputs, dec_outputs):
        super(MyDataSet, self).__init__()
        self.enc_inputs = enc_inputs
        self.dec_inputs = dec_inputs
        self.dec_outputs = dec_outputs

    def __len__(self):
        return self.enc_inputs.shape[0]

    def __getitem__(self, idx):
        return self.enc_inputs[idx], self.dec_inputs[idx], self.dec_outputs[idx]


class MyDataSetuser(Data.Dataset):
    def __init__(self, record_ids, enc_inputs, inputseqs, opmask_inputs):
        super(MyDataSetuser, self).__init__()
        self.record_ids = record_ids
        self.enc_inputs = enc_inputs
        self.inputseqs = inputseqs
        self.opmask_inputs = opmask_inputs

    def __len__(self):
        return self.enc_inputs.shape[0]

    def __getitem__(self, idx):
        return (
            self.record_ids[idx],
            self.enc_inputs[idx],
            self.inputseqs[idx],
            self.opmask_inputs[idx],
        )


protein_vocab = {
    "<pad>": 0,
    "G": 1,
    "A": 2,
    "V": 3,
    "L": 4,
    "I": 5,
    "P": 6,
    "F": 7,
    "Y": 8,
    "W": 9,
    "S": 10,
    "T": 11,
    "C": 12,
    "M": 13,
    "N": 14,
    "Q": 15,
    "D": 16,
    "E": 17,
    "K": 18,
    "R": 19,
    "H": 20,
    "*": 21,
}

codon_vocab = {
    "<pad>": 0,
    "TTT": 1,
    "TTC": 2,
    "TTA": 3,
    "TTG": 4,
    "TCT": 5,
    "TCC": 6,
    "TCA": 7,
    "TCG": 8,
    "TAT": 9,
    "TAC": 10,
    "TAA": 11,
    "TAG": 12,
    "TGT": 13,
    "TGC": 14,
    "TGA": 15,
    "TGG": 16,
    "CTT": 17,
    "CTC": 18,
    "CTA": 19,
    "CTG": 20,
    "CCC": 21,
    "CCA": 22,
    "CCG": 23,
    "CCT": 24,
    "CAC": 25,
    "CAT": 26,
    "CAA": 27,
    "CAG": 28,
    "CGC": 29,
    "CGA": 30,
    "CGG": 31,
    "CGT": 32,
    "ATT": 33,
    "ATC": 34,
    "ATA": 35,
    "ATG": 36,
    "ACC": 37,
    "ACA": 38,
    "ACG": 39,
    "ACT": 40,
    "AAC": 41,
    "AAT": 42,
    "AAA": 43,
    "AAG": 44,
    "AGC": 45,
    "AGT": 46,
    "AGA": 47,
    "AGG": 48,
    "GTT": 49,
    "GTC": 50,
    "GTA": 51,
    "GTG": 52,
    "GCC": 53,
    "GCA": 54,
    "GCG": 55,
    "GCT": 56,
    "GAC": 57,
    "GAT": 58,
    "GAA": 59,
    "GAG": 60,
    "GGC": 61,
    "GGA": 62,
    "GGG": 63,
    "GGT": 64,
    "<S>": 65,
    "<E>": 66,
}

pc_vocab = {
    "F": {
        "TTT": 1,
        "TTC": 2,
    },
    "L": {
        "TTA": 1,
        "TTG": 2,
        "CTT": 3,
        "CTC": 4,
        "CTA": 5,
        "CTG": 6,
    },
    "S": {
        "TCT": 1,
        "TCC": 2,
        "TCA": 3,
        "TCG": 4,
        "AGC": 5,
        "AGT": 6,
    },
    "Y": {
        "TAT": 1,
        "TAC": 2,
    },
    "*": {
        "TAA": 1,
        "TAG": 2,
        "TGA": 3,
    },
    "C": {
        "TGT": 1,
        "TGC": 2,
    },
    "W": {
        "TGG": 1,
    },
    "P": {
        "CCC": 1,
        "CCA": 2,
        "CCG": 3,
        "CCT": 4,
    },
    "H": {
        "CAC": 1,
        "CAT": 2,
    },
    "Q": {
        "CAA": 1,
        "CAG": 2,
    },
    "R": {
        "CGC": 1,
        "CGA": 2,
        "CGG": 3,
        "CGT": 4,
        "AGA": 5,
        "AGG": 6,
    },
    "I": {
        "ATT": 1,
        "ATC": 2,
        "ATA": 3,
    },
    "M": {
        "ATG": 1,
    },
    "T": {
        "ACC": 1,
        "ACA": 2,
        "ACG": 3,
        "ACT": 4,
    },
    "N": {
        "AAC": 1,
        "AAT": 2,
    },
    "K": {
        "AAA": 1,
        "AAG": 2,
    },
    "V": {
        "GTT": 1,
        "GTC": 2,
        "GTA": 3,
        "GTG": 4,
    },
    "A": {
        "GCC": 1,
        "GCA": 2,
        "GCG": 3,
        "GCT": 4,
    },
    "D": {
        "GAC": 1,
        "GAT": 2,
    },
    "E": {
        "GAA": 1,
        "GAG": 2,
    },
    "G": {
        "GGC": 1,
        "GGA": 2,
        "GGG": 3,
        "GGT": 4,
    },
    "<pad>": {
        "<pad>": 0,
    },
}
src_vocab_size = max(len(protein_vocab), len(codon_vocab))
tgt_vocab_size = max(len(protein_vocab), len(codon_vocab))
src_idx2word = {protein_vocab[key]: key for key in protein_vocab}
tgt_idx2word = {codon_vocab[key]: key for key in codon_vocab}


def making_opmask(protein):
    out = []
    for p in protein:
        cache = [0] * len(codon_vocab)
        probkey = pc_vocab[p].keys()
        problist = [codon_vocab[key] for key in probkey]
        for prob in problist:
            cache[int(prob)] = 1
        out.append(cache)

    out = out + [[0] * 66 + [1]] * ((pclen + 1) - len(protein))
    return out


def pad_list(onelist):
    onelist.extend((pclen - len(onelist)) * [0])
    return onelist


def encoding_codon(codon):
    return pad_list(
        [codon_vocab[codon[num : num + 3]] for num in range(0, len(codon), 3)]
    )


def encoding_protein(protein):
    return pad_list([protein_vocab[n] for n in protein])


def make_data(filepath=raw_output_file):
    with open(file=filepath, mode="r") as f:
        print(f"Opening {filepath}")
        protein_inputs, codon_inputs, tai_inputs = [], [], []
        for line, i in tqdm(enumerate(f)):
            name, pseq, cseq, tai = i.strip().split("\t")
            if len(pseq) >= minlen and len(pseq) <= pclen:
                try:
                    proteininput = encoding_protein(pseq)
                    codoninput = encoding_codon(cseq)

                    protein_inputs.append(proteininput)
                    codon_inputs.append(codoninput)
                    tai_inputs.append([float(tai)])
                except:
                    tqdm.write("Found an error in the %dth data" % line)
        tqdm.write("finish coding")
        batch_size = torch.LongTensor(protein_inputs).size(0)
        return (
            torch.cat(
                (
                    torch.LongTensor(protein_inputs),
                    torch.LongTensor(torch.full((batch_size, 1), fill_value=0)),
                ),
                dim=1,
            ),
            torch.cat(
                (
                    torch.LongTensor(torch.full((batch_size, 1), fill_value=65)),
                    torch.LongTensor(codon_inputs),
                ),
                dim=1,
            ),
            torch.cat(
                (
                    torch.LongTensor(codon_inputs),
                    torch.LongTensor(torch.full((batch_size, 1), fill_value=66)),
                ),
                dim=1,
            ),
        )


def getproteinlist(id):
    name = src_idx2word[id]
    outname = pc_vocab[name].keys()
    return [codon_vocab[i] for i in outname]


def make_data_user(filepath=raw_output_file):
    with open(file=filepath, mode="r") as f:
        print(f"Opening {filepath}")
        record_ids, protein_inputs, inputseqs, opmask_inputs = [], [], [], []
        for line, i in tqdm(enumerate(f)):
            name, pseq, iseq = i.strip().split("\t")
            if len(pseq) >= minlen and len(pseq) <= pclen:
                try:
                    record_ids.append(name)
                    proteininput = encoding_protein(pseq)
                    opmask_input = making_opmask(pseq)
                    protein_inputs.append(proteininput)
                    inputseqs.append(iseq)
                    opmask_inputs.append(opmask_input)
                except:
                    tqdm.write("Found an error in the %dth data" % line)
        tqdm.write("finish coding")
        batch_size = torch.LongTensor(protein_inputs).size(0)
        return (
            record_ids,
            torch.cat(
                (
                    torch.LongTensor(protein_inputs),
                    torch.LongTensor(torch.full((batch_size, 1), fill_value=0)),
                ),
                dim=1,
            ),
            inputseqs,
            torch.LongTensor(opmask_inputs),
        )
