from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

from DeepCodon.config.config import minlen, pclen, raw_input_file, raw_output_file

input_file = raw_input_file
output_file = raw_output_file

numless = 0
nummore = 0
with open(output_file, "w") as tsv_file:
    for record in tqdm(SeqIO.parse(input_file, "fasta")):
        dna_sequence = str(record.seq)
        try:
            tai_info = str(record.description).split(" ")[-1][5:-1]
        except:
            tai_info = "0"
        protein_sequence = str(Seq(dna_sequence).translate())
        if len(protein_sequence) <= pclen and len(protein_sequence) >= minlen:
            tsv_file.write(
                f"{record.id}\t{protein_sequence}\t{dna_sequence}\t{tai_info}\n"
            )
        else:
            if len(protein_sequence) < minlen:
                numless += 1
            else:
                nummore += 1
print(f"numless: {numless}, nummore: {nummore}")
