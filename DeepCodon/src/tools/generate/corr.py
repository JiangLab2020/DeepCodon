import torch
from Bio import SeqIO
from Bio.Seq import Seq

from DeepCodon.src.train.datasets import getproteinlist


def knowcodon(s):
    if s == "UUU" or s == "UUC" or s == "TTT" or s == "TTC" or s == "T\nT\nT":
        out = "苯丙氨酸 Phe/F"
    elif s == "UUA" or s == "UUG" or s == "TTA" or s == "TTG":
        out = "亮氨酸 Leu/L"
    elif (
        s == "UCU"
        or s == "UCC"
        or s == "UCA"
        or s == "UCG"
        or s == "TCT"
        or s == "TCC"
        or s == "TCA"
        or s == "TCG"
    ):
        out = "丝氨酸 Ser/S"
    elif s == "UAU" or s == "UAC" or s == "TAT" or s == "TAC":
        out = "酪氨酸 Tyr/Y"
    elif s == "UAA" or s == "UAG" or s == "TAA" or s == "TAG":
        out = "终止 Ter/end"
    elif s == "UGU" or s == "UGC" or s == "TGT" or s == "TGC":
        out = "半胱氨酸 Cys/C"
    elif s == "UGA" or s == "TGA":
        out = "终止 Ter/end"
    elif s == "UGG" or s == "TGG":
        out = "色氨酸 Trp/W"
    elif (
        s == "CUU"
        or s == "CUC"
        or s == "CUA"
        or s == "CUG"
        or s == "CTT"
        or s == "CTC"
        or s == "CTA"
        or s == "CTG"
    ):
        out = "亮氨酸 Leu/L"
    elif s == "CCU" or s == "CCC" or s == "CCA" or s == "CCG" or s == "CCT":
        out = "脯氨酸 Pro/P"
    elif s == "CAU" or s == "CAC" or s == "CAT":
        out = "组氨酸 His/H"
    elif s == "CAA" or s == "CAG":
        out = "谷氨酰胺 Gln/Q"
    elif s == "CGU" or s == "CGC" or s == "CGA" or s == "CGG" or s == "CGT":
        out = "精氨酸 Arg/R"
    elif (
        s == "AUU" or s == "AUC" or s == "AUA" or s == "ATT" or s == "ATC" or s == "ATA"
    ):
        out = "异亮氨酸 IIe/I"
    elif s == "AUG" or s == "ATG":
        out = "甲硫氨酸 Met/M"
    elif s == "ACU" or s == "ACC" or s == "ACA" or s == "ACG" or s == "ACT":
        out = "苏氨酸 Thr/T"
    elif s == "AAU" or s == "AAC" or s == "AAT":
        out = "天冬酰胺 Asn/N"
    elif s == "AAA" or s == "AAG":
        out = "赖氨酸 Lys/K"
    elif s == "AGU" or s == "AGC" or s == "AGT":
        out = "丝氨酸 Ser/S"
    elif s == "AGA" or s == "AGG":
        out = "精氨酸 Arg/R"
    elif (
        s == "GUU"
        or s == "GUC"
        or s == "GUA"
        or s == "GUG"
        or s == "GTT"
        or s == "GTC"
        or s == "GTA"
        or s == "GTG"
    ):
        out = "缬氨酸 val/V"
    elif s == "GCU" or s == "GCC" or s == "GCA" or s == "GCG" or s == "GCT":
        out = "丙氨酸 Ala/A"
    elif s == "GAU" or s == "GAC" or s == "GAT":
        out = "天冬酰胺 Asp/D"
    elif s == "GAA" or s == "GAG":
        out = "谷氨酸 Glu/E"
    elif s == "GGU" or s == "GGC" or s == "GGA" or s == "GGG" or s == "GGT":
        out = "甘氨酸 Gly/G"
    else:
        out = "不能识别"
    return out


def isprotein(old, new):
    oldList = [old[i : i + 3] for i in range(0, len(old), 3)]
    newList = [new[i : i + 3] for i in range(0, len(new), 3)]

    for i in range(max(len(oldList), len(newList))):
        try:
            p1 = knowcodon(oldList[i])
            p2 = knowcodon(newList[i])
            if p1 != p2:
                return "error"
        except:
            return "error"
    return "right"


def isproteins(old, new):
    protein_sequence = str(Seq(new).translate())
    try:
        if old != protein_sequence:
            return "error"
    except:
        return "error"
    return "right"


def isone(src, tgt, beam_width):
    src = src.item()
    proteinlist = getproteinlist(src)
    tgt = torch.index_select(tgt, 0, torch.tensor(proteinlist).cuda())
    top_probs = tgt
    top_idx = torch.LongTensor(proteinlist).cuda()

    while len(top_idx) < beam_width:
        padding_length = beam_width - len(top_idx)
        top_probs = torch.cat((top_probs, top_probs[:padding_length]))
        top_idx = torch.cat((top_idx, top_idx[:padding_length]))
    return top_probs, top_idx
