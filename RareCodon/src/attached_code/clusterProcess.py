import os
import subprocess
from collections import defaultdict

from Bio import SeqIO


def read_fasta(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.description] = str(record.seq)
    return sequences


def process(input_file_path, work_path):
    clusters = defaultdict(list)
    with open(f"{work_path}/tmp/DB_clu.tsv", "r") as f:
        for line in f:
            cluster_id, seq_id = line.strip().split("\t")
            cluster_id = cluster_id[1:-1]
            seq_id = seq_id[1:-1]
            clusters[cluster_id].append(seq_id)

    filtered_clusters = {k: v for k, v in clusters.items() if len(v) >= 5}

    sequences = read_fasta(input_file_path)
    output_dir = f"{work_path}/clustered_sequences"
    os.makedirs(output_dir, exist_ok=True)

    for num, cluster_id in enumerate(filtered_clusters.keys()):
        with open(os.path.join(output_dir, f"cluster_{num}.fasta"), "w") as f:
            for seq_id in filtered_clusters[cluster_id]:
                f.write(f">{seq_id}\n{sequences[seq_id]}\n")
        subprocess.run(
            f"mafft {output_dir}/cluster_{num}.fasta > {output_dir}/cluster_{num}_aligned.fasta",
            shell=True,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        # subprocess.run(f"rm {output_dir}/cluster_{num}.fasta", shell=True,check=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
