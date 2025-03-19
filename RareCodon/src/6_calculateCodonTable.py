import os
import subprocess
from collections import defaultdict

from Bio import SeqIO
from tqdm import tqdm

CONTINUE = False
if CONTINUE:
    subprocess.run("rm -rf ./data/7_rareTable/*", shell=True)
os.makedirs("./data/7_rareTable", exist_ok=True)
codons = [a + b + c for a in "ATGC" for b in "ATGC" for c in "ATGC"]


def count_codons(sequence):
    """Count the number of occurrences of 64 codons in a sequence"""
    codon_count = defaultdict(int)

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i : i + 3]
        if codon in codons:
            codon_count[codon] += 1
    return codon_count


def process_single_fasta(fasta_file, directory):
    file_path = os.path.join(directory, fasta_file)
    output_file = f"./data/7_rareTable/{fasta_file.split('.')[0]}.txt"

    if os.path.exists(output_file):
        return

    codon_count = defaultdict(int)

    for record in SeqIO.parse(file_path, "fasta"):
        seq = str(record.seq).upper()
        file_codon_count = count_codons(seq)
        for codon, count in file_codon_count.items():
            codon_count[codon] += count

    with open(output_file, "w") as f:
        for codon, count in codon_count.items():
            f.write(f"{codon}\t{count}\n")


def process_fasta_files(directory):
    """Parallel processing of all fna files in the directory"""
    fasta_files = [f for f in os.listdir(directory) if f.endswith(".fna")]
    for file in tqdm(fasta_files):
        process_single_fasta(file, directory)


if __name__ == "__main__":
    process_fasta_files("./data/5_genusData")
