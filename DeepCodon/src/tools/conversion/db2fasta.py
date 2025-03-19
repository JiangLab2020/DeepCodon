import sqlite3

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from DeepCodon.config.config import *
from DeepCodon.config.server_config import *

conn = sqlite3.connect(rareCodon_db_path)

query = """
SELECT fasta.seqid, fasta.sequence, token.label
FROM wp2UniRef
INNER JOIN token ON wp2UniRef.link = token.link
INNER JOIN fasta ON wp2UniRef.qseqid = fasta.seqid
WHERE token.label LIKE '%1%'
"""


df = pd.read_sql_query(query, conn)


conn.close()


df.to_csv(compare_file_path / "rareDB.csv", index=False)

df_sampled = df.sample(n=100, random_state=42)


fasta_records = [
    SeqRecord(Seq(row["sequence"]), id=row["seqid"], description="")
    for _, row in df_sampled.iterrows()
]

fasta_path = compare_file_path / "rareDB_sampled.fasta"
with open(fasta_path, "w") as fasta_file:
    SeqIO.write(fasta_records, fasta_file, "fasta")
