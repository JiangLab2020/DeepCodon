import os

from DeepCodon.src.tests.cai.cai import processing as calculate_cai
from DeepCodon.src.tools.generate.corr import isprotein, isproteins
from DeepCodon.src.tools.generate.gc import calculate_gc
from DeepCodon.src.tools.generate.mfe import calculate_mfe
from DeepCodon.src.train.datasets import *

current_path = os.getcwd()
refer_sheet = os.path.join(current_path, "src", "tests", "cai", "ecoli.txt")


def index2content(enc_inputs, predict, inputseqs, record_ids, whichdata):
    # Convert the output to
    predicted_words = [src_idx2word[n.item()] for n in enc_inputs.squeeze(0).squeeze()]
    n = (
        predicted_words.index("<pad>")
        if "<pad>" in predicted_words
        else len(predicted_words)
    )
    predicted_words = [tgt_idx2word[i.item()] for i in predict.squeeze()][1 : n + 1]
    predicted_words_str = "".join(predicted_words)
    if whichdata == "codon":
        try:
            istrue = isprotein(inputseqs[0], predicted_words_str)
        except:
            istrue = "error"
        output_str = f"{record_ids[0]},{inputseqs[0]},{predicted_words_str},{calculate_mfe(inputseqs[0])},{calculate_mfe(predicted_words_str)},{calculate_gc(inputseqs[0])},{calculate_gc(predicted_words_str)},{calculate_cai(infile=refer_sheet, inseq=inputseqs[0])},{calculate_cai(infile=refer_sheet, inseq=predicted_words_str)},{istrue}"
        return output_str
    elif whichdata == "protein":
        try:
            istrue = isproteins(inputseqs[0], predicted_words_str)
        except:
            istrue = "error"
        output_str = f"{record_ids[0]},{inputseqs[0]},{predicted_words_str},{calculate_mfe(predicted_words_str)},{calculate_gc(predicted_words_str)},{calculate_cai(infile=refer_sheet, inseq=predicted_words_str)},{istrue}"
        return output_str


def index2test(enc_inputs, dec_outputs, predict):
    predicted_words = [src_idx2word[n.item()] for n in enc_inputs.squeeze(0).squeeze()]
    n = (
        predicted_words.index("<pad>")
        if "<pad>" in predicted_words
        else len(predicted_words)
    )
    predicted_words = [tgt_idx2word[i.item()] for i in predict.squeeze()][0:n]
    predicted_words_str = "".join(predicted_words)
    dec_words = [tgt_idx2word[i.item()] for i in dec_outputs.squeeze()][0:n]
    dec_words_str = "".join(dec_words)
    output_str = f"{dec_words_str},{predicted_words_str},{calculate_gc(dec_words_str)},{calculate_gc(predicted_words_str)}"
    return output_str
