# Used for joint filtering of datasets by Cai and GC
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import CodonAdaptationIndex
from tqdm import tqdm

from DeepCodon.config.config import (
    PROJECT_ROOT,
    minlen,
    pclen,
    raw_input_file,
    raw_output_file,
)
from DeepCodon.src.tools.generate.gc import calculate_gc

ref_file = PROJECT_ROOT / "EPCD" / "finalData" / "final_tai.fasta"
caiseqs = SeqIO.parse(ref_file, "fasta")

cai = CodonAdaptationIndex(sequences=caiseqs)

input_file = raw_input_file
output_file = raw_output_file

with open(output_file, "w") as tsv_file:
    for record in tqdm(SeqIO.parse(input_file, "fasta")):
        dna_sequence = str(record.seq)
        try:
            tai_info = str(record.description).split(" ")[-1][5:-1]
        except:
            tai_info = "0"
        protein_sequence = str(Seq(dna_sequence).translate())
        if len(protein_sequence) <= pclen and len(protein_sequence) >= minlen:
            if 0.67 > calculate_gc(dna_sequence) > 0.41:
                if 0.89 > cai.calculate(dna_sequence) > 0.50:
                    tsv_file.write(
                        f"{record.id}\t{protein_sequence}\t{dna_sequence}\t{tai_info}\n"
                    )
