import os
import subprocess
from concurrent.futures import ProcessPoolExecutor

from Bio import SeqIO
from tqdm import tqdm

DEBUG = False


input_folder = "./data/3_allFasta"
output_folder = "./data/4_filterSeq"
subprocess.run(f"rm -rf {output_folder}", shell=True)


start_codon = "ATG"
stop_codons = ["TAA", "TAG", "TGA"]
valid_nucleotides = set("ATGC")


def process_file(input_file_path, protein_output_folder, codon_output_folder):
    with open(protein_output_folder, "w") as Poutput_handle:
        with open(codon_output_folder, "w") as Coutput_handle:
            for record in SeqIO.parse(input_file_path, "fasta"):
                dna_seq = record.seq
                if set(dna_seq.upper()) <= valid_nucleotides:
                    if dna_seq[:3] == start_codon and dna_seq[-3:] in stop_codons:
                        if len(dna_seq) % 3 == 0:
                            protein_seq = dna_seq.translate()
                            if protein_seq.count("*") == 1:
                                if 50 <= len(protein_seq) <= 1000:
                                    SeqIO.write(record, Coutput_handle, "fasta")
                                    record.seq = protein_seq
                                    SeqIO.write(record, Poutput_handle, "fasta")


def main(input_folder, protein_output_folder, codon_output_folder, debug=False):
    if not os.path.exists(protein_output_folder):
        os.makedirs(protein_output_folder)
    if not os.path.exists(codon_output_folder):
        os.makedirs(codon_output_folder)

    with ProcessPoolExecutor() as executor:
        futures = []
        for filename in tqdm(os.listdir(input_folder)):
            if filename.endswith(".fna"):
                input_file_path = os.path.join(input_folder, filename)
                protein_output_file_path = os.path.join(
                    protein_output_folder, filename.replace(".fna", ".faa")
                )
                codon_output_file_path = os.path.join(codon_output_folder, filename)
                futures.append(
                    executor.submit(
                        process_file,
                        input_file_path,
                        protein_output_file_path,
                        codon_output_file_path,
                    )
                )
                if DEBUG:
                    if len(futures) > 2:
                        break
        for future in tqdm(futures):
            future.result()

    print("all done")


if __name__ == "__main__":
    main(input_folder, output_folder, output_folder, debug=DEBUG)
