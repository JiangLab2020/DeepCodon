import pandas as pd
from Bio import SeqIO

from DeepCodon.src.attached_code.config import *

if SLIDE:
    records = list(SeqIO.parse("data/9_mapping/output_sequences_slide.fasta", "fasta"))
    dfoutpath = "data/9_mapping/dfRIGHT_sorted_slide.csv"
    tsvpath = "data/9_mapping/matches_slide.tsv"
else:
    records = list(SeqIO.parse("data/9_mapping/output_sequences.fasta", "fasta"))
    dfoutpath = "data/9_mapping/dfRIGHT_sorted.csv"
    tsvpath = "data/9_mapping/matches.tsv"

dfin = pd.DataFrame(
    {"id": [rec.id for rec in records], "sequence": [str(rec.seq) for rec in records]}
)
dfin["sequence_length"] = dfin["sequence"].apply(len)
df = pd.read_csv(
    tsvpath,
    sep="\t",
    names=[
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ],
)
dfout = pd.merge(df, dfin, left_on="qseqid", right_on="id")
dfRIGHT = dfout[dfout["length"] == dfout["sequence_length"]]
dfRIGHT_sorted = dfRIGHT.sort_values(by="pident", ascending=False).drop_duplicates(
    subset="qseqid", keep="first"
)
print("Data coverage percentage", len(dfRIGHT_sorted) / len(dfin))
print("Covered", len(dfRIGHT_sorted))
print("Total length", len(dfin))
dfRIGHT_sorted[["qseqid", "sseqid"]].to_csv(dfoutpath, index=False)
