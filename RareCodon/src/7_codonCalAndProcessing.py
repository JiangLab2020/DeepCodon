import concurrent.futures
import glob
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

import pandas as pd
from Bio import SeqIO

from DeepCodon.src.attached_code.codonAndProtein import *
from DeepCodon.src.attached_code.config import *
from DeepCodon.src.attached_code.misc_code import *

# setting
seq_len_tolerance = 0.2  # protein seq offset tolerance
protein_freq_tolerance = 0.8  # min protein levenstein similarity
min_seq_tolerance = 5  # minimum number of sequences in a cluster
sliding_window = 6

# logs
CLUSTER_LOG_FILE = "./log/cluster_process.log"
CLUSTER_LOG_FILE_J = "./log/cluster_process_J.log"
if SLIDE:
    OUT_PATH = "./data/8_genusOut_slide"
    LOG_FILE = "./log/find_cluster_process_slide.log"
    LOG_FILE_J = "./log/find_cluster_process_slide_J.log"
else:
    OUT_PATH = "./data/8_genusOut"
    LOG_FILE = "./log/find_cluster_process.log"
    LOG_FILE_J = "./log/find_cluster_process_J.log"

todoList = list()
with open("log/cluster_process.log", "r") as f:
    for data in f:
        data = data.strip()
        filename = (
            data.split(" ")[-1]
            .replace("5_genusData", "6_clusterOut")
            .replace(".faa", "")
        )
        todoList.append(filename)


def process_codonTable(genusname):
    # 0 Initialize codon frequency table
    for key1 in pc_vocab:
        for key2 in pc_vocab[key1]:
            pc_vocab[key1][key2] = 0
    # 1 Read password sub table
    with open(f"data/7_rareTable/{genusname}.txt", "r") as f:
        for line in f:
            line = line.strip()
            codon, num = line.split("\t")
            # Input the results into the dictionary
            pc_vocab[cp_vocab[codon]][codon] = int(num)
    # 2 Normalize CFD using the maximum value
    cfd_pc_vocab = {}

    for key1, inner_dict in pc_vocab.items():
        max_value = max(inner_dict.values())
        new_inner_dict = {}
        for key2, value in inner_dict.items():
            new_inner_dict[key2] = value / max_value if max_value != 0 else 0
        cfd_pc_vocab[key1] = new_inner_dict
    # 3 Calculate mask or nomask using CFD
    mask_pc_vocab = {}
    mask_c_vocab = {}
    for key1, inner_dict in cfd_pc_vocab.items():
        new_inner_dict = {}
        for key2, value in inner_dict.items():
            if value >= 0.3:
                new_inner_dict[key2] = 0
                mask_c_vocab[key2] = 0
            else:
                new_inner_dict[key2] = 1
                mask_c_vocab[key2] = 1
        mask_pc_vocab[key1] = new_inner_dict
    return mask_c_vocab


def process_align(rare_table, align_file, df_genus):
    """
    Processing MSA files
    MaskTable: Mask Table
    align_file:  Pending Files
    genusname:  genus name
    """
    records = list(SeqIO.parse(align_file, "fasta"))
    df_align = pd.DataFrame(
        {
            "id": [rec.id for rec in records],
            "protein": [str(rec.seq).replace("*", "") for rec in records],
            "len": [len(rec.seq) for rec in records],
        }
    )
    df_align = pd.merge(df_align, df_genus, on="id", how="inner")
    # Select conservative rare codons (with a length within a certain range, num greater than 5 columns, covering 80% of amino acids)
    # 1. Remove sequences that exceed the mode range
    mode_len = df_align["len"].mode()[0]
    tolerance = seq_len_tolerance * mode_len
    df_align = df_align[
        (df_align["len"] >= mode_len - tolerance)
        & (df_align["len"] <= mode_len + tolerance)
    ]
    #  2 80%
    df_align = filter_by_levenshtein_similarity(
        df_align, "protein", protein_freq_tolerance
    )
    # 3 <5
    if len(df_align) < min_seq_tolerance:
        return None
    df_align = give_codonGap(dfname=df_align, inputcol="codon", outputcol="codon_gap")

    df_align["raretable"] = df_align["codon_gap"].apply(
        codon_to_vector, args=(rare_table,)
    )

    df_align["raretable"] = rare_codon_addAND(df_align, "raretable")
    # df_align["raretable_slide"] = vector_with_slidingWindow(
    # df_align["raretable"].iloc[0], sliding_window
    # )
    return df_align


def is_already_processed(file_path, log_file):
    pattern = r"data/6_clusterOut/([A-Za-z]+)"
    with open(log_file, "r") as f:
        lines = f.readlines()
        matches = re.findall(pattern, str(lines))
    return file_path.split("/")[-1] in matches


def process_align_wrapper(align_file, rare_table, df_genus):
    return process_align(
        rare_table=rare_table, align_file=align_file, df_genus=df_genus
    )


def process_pipeline(item):
    genus_seq_dict = dict()
    if not os.path.exists(item):
        print(f"Path does not exist: {item}")
    cluster_tree_path = os.path.join(item, "clustered_sequences")
    genusname = item.split("/")[-1]
    mask_c_vocab = process_codonTable(genusname)
    if not os.path.isdir(cluster_tree_path):
        print(f"ClusterTree folder does not exist: {cluster_tree_path}")
    else:
        align_files = glob.glob(os.path.join(cluster_tree_path, "*aligned.fasta"))
        if align_files:
            records_genus = list(
                SeqIO.parse(f"data/5_genusData/{genusname}.fna", "fasta")
            )
            df_genus = pd.DataFrame(
                {
                    "id": [rec.id for rec in records_genus],
                    "codon": [str(rec.seq)[:-3] for rec in records_genus],
                    "genus": genusname,
                }
            )
            with ThreadPoolExecutor() as executor:
                future_to_file = {
                    executor.submit(
                        process_align_wrapper, align_file, mask_c_vocab, df_genus
                    ): align_file
                    for align_file in align_files
                }
                for future in as_completed(future_to_file):
                    align_file = future_to_file[future]
                    try:
                        dfout = future.result()
                        if dfout is not None:
                            pname, pseq, pmask = get_best_seq(dfout)
                            genus_seq_dict[pname] = [pseq, pmask]
                    except Exception as exc:
                        print(f"Error processing {align_file}: {exc}")
        else:
            print(f"No .fasta files found in: {cluster_tree_path}")
    df_genus_out = pd.DataFrame(
        {
            "id": list(genus_seq_dict.keys()),
            "protein": [rec[0] for rec in genus_seq_dict.values()],
            "mask": [rec[1] for rec in genus_seq_dict.values()],
            "genus": genusname,
        }
    )
    df_genus_out.to_csv(os.path.join(OUT_PATH, f"{genusname}.csv"), index=False)
    return len(align_files)


def process_file(file_path):
    try:
        if not os.path.exists(LOG_FILE):
            open(LOG_FILE, "w").close()
        if not os.path.exists(LOG_FILE_J):
            open(LOG_FILE_J, "w").close()
        if not is_already_processed(file_path, LOG_FILE):
            num_cluster_out = process_pipeline(file_path)
            with open(LOG_FILE, "a") as log_f:
                log_f.write(
                    f"{datetime.now()} - done: {file_path}\tcluster:{num_cluster_out}\n"
                )
        else:
            with open(LOG_FILE_J, "a") as log_f_j:
                log_f_j.write(f"{datetime.now()} - skip: {file_path}\n")
    except Exception as e:
        print(f"Error processing {item}: {e}")
        log_f.write(f"{datetime.now()} - error: {file_path}\n")


if __name__ == "__main__":
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     futures = [executor.submit(process_file, item) for item in todoList]
    #     for future in concurrent.futures.as_completed(futures):
    #         try:
    #             result = future.result()
    #             print(f"Result: {result}")
    #         except Exception as exc:
    #             print(f"Task failed with exception: {exc}")

    with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor:
        executor.map(process_file, todoList)

    # for item in todoList:
    #     process_file(item)
