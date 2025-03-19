import pandas as pd
from Bio import SeqIO

from DeepCodon.src.attached_code.config import *

if SLIDE:
    fastadir = "data/9_mapping/output_sequences_slide.fasta"
    founddir = "data/9_mapping/found_files_slide.csv"
else:
    fastadir = "data/9_mapping/output_sequences.fasta"
    founddir = "data/9_mapping/found_files.csv"

records = list(SeqIO.parse(fastadir, "fasta"))
df_fasta = pd.DataFrame(
    {"id": [rec.id for rec in records], "sequence": [str(rec.seq) for rec in records]}
)

df_found = pd.read_csv(founddir)
df_found["name"] = (
    "AF-" + df_found["sseqid"].str.split("_").str[-1] + "-F1-model_v4.cif"
)

df = pd.merge(df_fasta, df_found, left_on="id", right_on="qseqid", how="right")
df = df.rename(columns={"sequence": "protein"})
df[["name", "protein"]].to_csv("data/9_mapping/add_seq.csv", index=False)
