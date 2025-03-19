import gzip
import os
import shutil
import subprocess

from tqdm import tqdm

source_dir = "./data/2_allSeq"
target_dir = "./data/3_allFasta"
subprocess.run(f"rm -rf {target_dir}", shell=True)

os.makedirs(target_dir, exist_ok=True)


for file_name in tqdm(os.listdir(source_dir)):
    if file_name.endswith(".gz"):
        source_file = os.path.join(source_dir, file_name)
        target_file = os.path.join(target_dir, file_name[:-3])

        try:
            with gzip.open(source_file, "rb") as f_in:
                with open(target_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        except EOFError:
            print(f"Error: {source_file} is incomplete or corrupted.")
        except Exception as e:
            print(f"An error occurred with {source_file}: {e}")
