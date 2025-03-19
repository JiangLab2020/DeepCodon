import os
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from DeepCodon.src.attached_code.config import *

if SLIDE:
    data_dir = "data/8_genusOut_slide"
    output_fasta = "data/9_mapping/output_sequences_slide.fasta"
else:
    data_dir = "data/8_genusOut"
    output_fasta = "data/9_mapping/output_sequences.fasta"


def read_csv_file(filepath):
    return pd.read_csv(filepath)


csv_files = [
    os.path.join(data_dir, f) for f in os.listdir(data_dir) if f.endswith(".csv")
]

with ProcessPoolExecutor() as executor:
    all_data = pd.concat(
        tqdm(
            executor.map(read_csv_file, csv_files),
            total=len(csv_files),
            desc="reading csv files",
        ),
        ignore_index=True,
    )
fasta_records = [
    SeqRecord(Seq(protein), id=str(id_), description="")
    for id_, protein in tqdm(
        zip(all_data["id"], all_data["protein"]),
        total=all_data.shape[0],
        desc="fasta records",
    )
]

with open(output_fasta, "w") as fasta_file:
    SeqIO.write(fasta_records, fasta_file, "fasta")

print(f"done :{output_fasta}")
