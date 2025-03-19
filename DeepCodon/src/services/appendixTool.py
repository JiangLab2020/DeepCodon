from Bio import SeqIO
from Bio.Seq import Seq

from DeepCodon.src.services.preBuildLogitBias import property


def check_and_convert_to_aa(sequence):
    # If the string contains other characters besides A, T, C, G, continue
    aa_sequence = sequence
    if set(sequence) <= {"A", "T", "C", "G", "a", "t", "c", "g"}:
        seq_obj = Seq(sequence)
        aa_sequence = str(seq_obj.translate())  # Translated as amino acid sequence
    # Check if the termination codon exists
    if "*" not in aa_sequence:
        aa_sequence += "*"  # If there is no termination codon, add*

    return aa_sequence


def convert_fasta_to_aa(input_file, output_file):
    # Convert the nucleotide sequence in the input fasta file to an amino acid sequence
    with open(output_file, "w") as out_file:
        for record in SeqIO.parse(input_file, "fasta"):
            # Obtain the original sequence
            original_seq = str(record.seq)
            # Check and convert to amino acid sequence
            aa_seq = check_and_convert_to_aa(original_seq)
            # Write the results back to the output file
            record.seq = Seq(aa_seq)
            SeqIO.write(record, out_file, "fasta")
    print(f"Converted fasta file to amino acid sequence: {output_file}")
    return output_file


def assert_all_codon_rare(file):
    # token_dict = {query: [DatabaseToken, DatabaseSeq, QueryToken, SmoothToken]}
    token_dict = {}
    filtered_codons = [
        k for k, v in property.items() if v < 0.3 and k not in {"<pad>", "<S>", "<E>"}
    ]
    with open(file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq = str(record.seq)
            id = record.id
            tokens = ""
            for i in range(0, len(seq), 3):
                codon = seq[i : i + 3]
                # TODO
                if codon in filtered_codons:
                    tokens += "1"
                else:
                    tokens += "0"
            token_dict[id] = ["", "", tokens]
    return token_dict
