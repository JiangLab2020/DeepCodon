import concurrent.futures
import os
import subprocess
from datetime import datetime

from DeepCodon.src.attached_code.clusterProcess import process

CONTINUE = False
if CONTINUE:
    subprocess.run("rm -rf ./data/6_clusterOut/*", shell=True)

INPUT_PATH = "./data/5_genusData"
OUTPUT_PATH = "./data/6_clusterOut"
LOG_FILE = "./log/cluster_process.log"
LOG_FILE_J = "./log/cluster_process_J.log"


def cluster_pipeline(input_path):
    input_name = os.path.basename(input_path)[:-4]
    subprocess.run(
        f"rm -rf {OUTPUT_PATH}/{input_name}",
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        f"mkdir -p {OUTPUT_PATH}/{input_name}/tmp",
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        f"mmseqs createdb {input_path} {OUTPUT_PATH}/{input_name}/tmp/DB",
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        f"mmseqs cluster {OUTPUT_PATH}/{input_name}/tmp/DB {OUTPUT_PATH}/{input_name}/tmp/DB_clu {OUTPUT_PATH}/{input_name}/tmp  --min-seq-id 0.9 -c 0.8",
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        f"mmseqs createtsv {OUTPUT_PATH}/{input_name}/tmp/DB {OUTPUT_PATH}/{input_name}/tmp/DB {OUTPUT_PATH}/{input_name}/tmp/DB_clu {OUTPUT_PATH}/{input_name}/tmp/DB_clu.tsv  --full-header",
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    process(input_file_path=input_path, work_path=f"{OUTPUT_PATH}/{input_name}")


def is_already_processed(file_path, log_file):
    with open(log_file, "r") as f:
        return file_path in f.read()


def process_file(file_path):
    if not os.path.exists(LOG_FILE):
        open(LOG_FILE, "w").close()
    if not os.path.exists(LOG_FILE_J):
        open(LOG_FILE_J, "w").close()
    if not is_already_processed(file_path, LOG_FILE):
        cluster_pipeline(file_path)
        with open(LOG_FILE, "a") as log_f:
            log_f.write(f"{datetime.now()} - have done: {file_path}\n")
    else:
        with open(LOG_FILE_J, "a") as log_f_j:
            log_f_j.write(f"{datetime.now()} - skip: {file_path}\n")


if __name__ == "__main__":
    files = [
        os.path.join(INPUT_PATH, f)
        for f in os.listdir(INPUT_PATH)
        if os.path.isfile(os.path.join(INPUT_PATH, f)) and f.endswith(".faa")
    ]
    with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
        executor.map(process_file, files)
