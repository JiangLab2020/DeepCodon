import os
import subprocess
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from tqdm import tqdm

df_grouped = pd.read_pickle("misc/df_grouped.pkl")
df_groupeds = df_grouped.groupby("genus")
subprocess.run("rm -rf ./data/5_genusData", shell=True)
if not os.path.exists("./data/5_genusData"):
    os.makedirs("./data/5_genusData")


def append_to_fna_faa(group):
    name, group = group
    for index, row in group.iterrows():
        fna_file_path = os.path.join("./data/5_genusData", f"{name}.fna")
        faa_file_path = os.path.join("./data/5_genusData", f"{name}.faa")
        file_in_fna = os.path.join("./data/4_filterSeq", row["fileName"][:-3])
        file_in_faa = file_in_fna.replace("fna", "faa")

        with open(file_in_fna, "r") as fileIn:
            with open(fna_file_path, "a") as fna:
                fna.write(fileIn.read())

        with open(file_in_faa, "r") as fileIn:
            with open(faa_file_path, "a") as faa:
                faa.write(fileIn.read())


with ProcessPoolExecutor() as executor:
    list(tqdm(executor.map(append_to_fna_faa, df_groupeds), total=len(df_groupeds)))
