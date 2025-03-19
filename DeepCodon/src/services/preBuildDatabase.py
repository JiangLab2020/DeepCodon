import sqlite3

import pandas as pd
from Bio import SeqIO

from DeepCodon.config.server_config import (
    blastp_dataset_path,
    rareCodon_db_path,
    rareCodon_path,
    wp2UniRef_path,
)

# read csv
df = pd.read_csv(wp2UniRef_path)
df["link"] = df["sseqid"].apply(lambda x: x.split("_")[1])

conn = sqlite3.connect(rareCodon_db_path)
cursor = conn.cursor()


cursor.execute("""
CREATE TABLE IF NOT EXISTS wp2UniRef (
    qseqid TEXT,
    sseqid TEXT,
    link TEXT
);
""")


df.to_sql("wp2UniRef", conn, if_exists="replace", index=False)

# Create an index to accelerate queries
cursor.execute("CREATE INDEX IF NOT EXISTS idx_qseqid ON wp2UniRef(qseqid);")
cursor.execute("CREATE INDEX IF NOT EXISTS idx_sseqid ON wp2UniRef(sseqid);")
cursor.execute("CREATE INDEX IF NOT EXISTS idx_link ON wp2UniRef(link);")

# ----------------------------
df = pd.read_csv(rareCodon_path)
df["link"] = df["sequence"].apply(lambda x: x.split("-")[1])
cursor.execute("""
CREATE TABLE IF NOT EXISTS token (
    link TEXT,
    sequence TEXT,
    type TEXT,
    chain TEXT,
    label TEXT,
    stage TEXT
);
""")


df.to_sql("token", conn, if_exists="replace", index=False)

cursor.execute("CREATE INDEX IF NOT EXISTS idx_sequence ON token(link);")

cursor.execute("""
CREATE TABLE IF NOT EXISTS fasta (
    seqid TEXT PRIMARY KEY,
    sequence TEXT
);
""")

fasta_data = []
for record in SeqIO.parse(blastp_dataset_path, "fasta"):
    fasta_data.append((record.id, str(record.seq)))


cursor.executemany("INSERT INTO fasta (seqid, sequence) VALUES (?, ?)", fasta_data)

cursor.execute("CREATE INDEX IF NOT EXISTS idx_seqid ON fasta(seqid);")


conn.commit()
conn.close()

print("Data import completed and index created.")
