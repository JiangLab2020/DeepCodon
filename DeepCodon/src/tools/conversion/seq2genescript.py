from Bio import SeqIO

input_file = "data/test_data/teach/compare1000.fasta"
output_prefix = "data/test_data/teach/geneScript/output_"


file_count = 1
seq_count = 0
max_seqs_per_file = 100


output_file = f"{output_prefix}{file_count}.fasta"
output_handle = open(output_file, "w")


for seq_record in SeqIO.parse(input_file, "fasta"):
    if seq_count >= max_seqs_per_file:
        output_handle.close()
        file_count += 1
        seq_count = 0
        output_file = f"{output_prefix}{file_count}.fasta"
        output_handle = open(output_file, "w")

    SeqIO.write(seq_record, output_handle, "fasta")
    seq_count += 1
output_handle.close()
