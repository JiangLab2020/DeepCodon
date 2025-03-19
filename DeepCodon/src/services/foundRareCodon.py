import sqlite3
import subprocess
from io import StringIO

import numpy as np
import pandas as pd
from Bio import SeqIO

from DeepCodon.config.server_config import *


def getSeqsLength(input_file: str) -> int:
    """
    Obtain the sequence length to determine if it meets the criteria
    """
    seq_dict = {
        record.id: len(record.seq) for record in SeqIO.parse(input_file, "fasta")
    }
    return seq_dict


def blastp4RareCodon(nowUuid: str, input_file: str) -> str:
    #  BLASTP
    blastp_cline = [
        "blastp",
        "-query",
        input_file,  # query sequence
        "-db",
        blastp_dataset_path,  # Local BLAST database
        "-num_threads",
        "10",
        "-outfmt",
        "6",  # Select XML format output
        "-out",
        blastp_output_path / f"{nowUuid}.tsv",  # output file
    ]
    result = subprocess.run(
        blastp_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    return result.stdout.decode("utf-8"), result.stderr.decode("utf-8")


def getBlastpResult(nowUuid: str, input_file: str) -> str:
    """
    For multiple sequences, it is necessary to take the maximum value separately and then determine whether the conditions are met before returning the token
    """
    df_blastp = pd.read_csv(
        blastp_output_path / f"{nowUuid}.tsv",
        sep="\t",
        header=None,
        names=[
            "query",
            "subject",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ],
    )
    # Store the length of each sequence
    seqLen = getSeqsLength(input_file)
    df_blastp["query_len"] = df_blastp["query"].map(seqLen)
    df_blastp["coverage"] = df_blastp["length"] / df_blastp["query_len"]
    df_output = df_blastp[df_blastp["coverage"] > coverage_threshold]
    df_output = df_output[df_output["pident"] > pident_threshold * 100]

    # Take the maximum value
    df_max = df_output.groupby("query").apply(
        lambda g: g.nlargest(1, ["pident", "coverage"])
    )
    # Convert the result to a dictionary, as the maximum value has already been taken above, there is only one line here
    # Save the entire row of data corresponding to each query as a list
    result_dict = (
        df_max.set_index("query").apply(lambda row: row.tolist(), axis=1).to_dict()
    )
    return result_dict


def searchDB(input_query: dict) -> dict:
    """
    Used to search the database and return tokens.
    return:  Return the dictionary, where key is the original index and value is the list of labels queried.
    """
    # Connect to database
    conn = sqlite3.connect(rareCodon_db_path)
    cursor = conn.cursor()

    result_dict = {}  # Store query results

    try:
        for key, value in input_query.items():
            if not value or len(value) == 0:
                result_dict[key] = []  # Empty value returns an empty list
                continue

            query_value = value[0]  # Take the first element of value

            # Query the wp2UniRef table to obtain sseqid
            cursor.execute(
                "SELECT link FROM wp2UniRef WHERE qseqid = ?", (query_value,)
            )
            link_results = cursor.fetchall()
            labels = []  # Store the queried label values
            for (link,) in link_results:
                # Query the token table to obtain the label
                cursor.execute("SELECT label FROM token WHERE link = ?", (link,))
                label_results = cursor.fetchall()
                # If the query result is empty, skip it
                if not label_results:
                    continue
                # Extract labels and add them to the list
                labels.extend(label[0] for label in label_results)
                # Query the fasta table to obtain the sequence corresponding to sseqid
                cursor.execute(
                    "SELECT sequence FROM fasta WHERE seqid = ?", (query_value,)
                )
                sequence_results = cursor.fetchall()
                labels.extend(seq[0] for seq in sequence_results)
            # If the query result is empty, skip it
            if len(labels) == 0:
                continue
            result_dict[key] = labels  # store the query result
    except Exception as e:
        print(f"database error: {e}")

    finally:
        # shutdown
        conn.close()
    return result_dict


def doAlign(seqI: str, seqD: str) -> str:
    """
    do alignment
    param seqI: input sequence
    param seqD: database sequence
    return: alignment result
    """
    # NOTE: MAFFT is required to be installed on the server
    mafft_cline = [
        "mafft",
        "--auto",
        "--anysymbol",
        "-",
    ]
    # do MAFFT
    result = subprocess.run(
        mafft_cline,
        input=f">seqI\n{seqI}\n>seqD\n{seqD}",
        text=True,
        capture_output=True,
    )
    return result.stdout


def AlignSeqs(inseqPath: str, inToken: dict) -> str:
    """
    Align the input sequence with the database sequence and return the token.
    """
    # NOTE Only keep tokens of equal length
    inToken = {k: v for k, v in inToken.items() if len(v[1]) == len(v[0])}
    # Read input sequence
    inseq_dict = {}
    for record in SeqIO.parse(inseqPath, "fasta"):
        inseq_dict[record.id] = record.seq
    for key, value in inToken.items():
        if not value or len(value) == 0:
            continue
        alignOut = doAlign(seqI=inseq_dict[key], seqD=value[1])
        print(alignOut)  # TODO
        fasta_io = StringIO(alignOut)
        records = list(SeqIO.parse(fasta_io, "fasta"))
        tokenD = ""
        num = 0
        for aa in records[1].seq:
            if aa == "-":
                tokenD += "0"
            else:
                tokenD += value[0][num]
                num += 1
        print("tokenD:", tokenD)  # TODO
        tokenI = ""
        num = 0
        for aa in records[0].seq:
            if aa == "-":
                pass
            else:
                tokenI += tokenD[num]
            num += 1
        print("tokenI:", tokenI)  # TODO
        inToken[key].append(tokenI)
    return inToken


def gaussian_window(size: int, sigma: float = 1.0):
    """
    Generate a Gaussian distribution weight window with the highest center value and gradually decreasing on both sides.
    : param size:  window size
    : param sigma:  Control the width of the distribution, the larger the attenuation, the slower the attenuation
    : return:  Normalized weight array
    """
    center = size // 2
    x = np.arange(size) - center
    weights = np.exp(-0.5 * (x / sigma) ** 2)  # Gaussian formula
    return weights / weights.sum()  # Normalize to ensure that the total sum is 1


def smooth_binary_sequence_with_decay(
    token_str: str, window_size: int = 5, sigma: float = 1.0
) -> list:
    """
    Perform sliding window processing on strings with the shape of '00000 100001', causing center 1 to decay towards both sides and considering the impact of overlap.
    : param token_str:  A string containing only '0' and '1'
    : param window_size:  Sliding window size
    : param sigma:  Control the degree of Gaussian attenuation
    : return:  A decimal list representing smoothed values
    """
    token_ints = np.array(
        [int(c) for c in token_str], dtype=float
    )  # Convert to numerical array
    smoothed = np.zeros_like(token_ints, dtype=float)  # Initialize smooth array
    weights = gaussian_window(window_size, sigma)  # Calculate Gaussian window weights
    half_win = window_size // 2

    # Traverse every 1 in the sequence and apply Gaussian weights
    for i, val in enumerate(token_ints):
        if val == 1:
            start = max(0, i - half_win)
            end = min(len(token_ints), i + half_win + 1)

            w_start = half_win - (
                i - start
            )  # Calculate the starting point of the corresponding weight window
            w_end = half_win + (end - i)

            smoothed[start:end] += weights[w_start:w_end]  # Accumulate Gaussian weights

    return smoothed.tolist()


def smoon(inToken: dict, userWindowSize=None) -> dict:
    """
    Used to perform sliding window calculation and return a smooth value of 0-1
    """
    for key, value in inToken.items():
        if userWindowSize is not None:
            smooth_window_size = userWindowSize
        else:
            from DeepCodon.config.server_config import smooth_window_size
        smoothToken = smooth_binary_sequence_with_decay(
            value[2], window_size=smooth_window_size, sigma=1.0
        )
        inToken[key].append(smoothToken)

    return inToken
