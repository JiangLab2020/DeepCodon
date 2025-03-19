import argparse
import logging
import uuid

from DeepCodon.config.server_config import *
from DeepCodon.src.generate.processing import modelProcess
from DeepCodon.src.services.appendixTool import convert_fasta_to_aa
from DeepCodon.src.services.foundRareCodon import (
    AlignSeqs,
    blastp4RareCodon,
    getBlastpResult,
    searchDB,
    smoon,
)

# -----------------args----------------------
parser = argparse.ArgumentParser(description="Run deepcodon")

parser.add_argument("-i", "--input_file", type=str, help="Input file", required=True)
parser.add_argument("-o", "--output_file", type=str, help="Output file", required=False)
args = parser.parse_args()

nowUuid = uuid.uuid4().hex
# -----------------log----------------------
logging.basicConfig(
    filename=server_log_path / f"{nowUuid}.log",
    level=logging.INFO,
    encoding="utf-8",
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y.%m.%d %H:%M:%S",
)

logging.info(f"Input file: {args.input_file}")
logging.info(f"Output file: {args.output_file}")
logging.info("🔔🔔🔔 Start 🔔🔔🔔")
# -------------run-----------------
# 1.  Used to execute blastp to find similar sequences and return tokens
# 1.1 Execute blastp
# NOTE token dict = {query: [DatabaseToken,DatabaseSeq,QueryToken,SmoothToken]}
logging.info("🔔🔔 RareCodon start 🔔🔔")
user_input_file = convert_fasta_to_aa(args.input_file, str(args.input_file) + "_aa")
bout, berr = blastp4RareCodon(nowUuid, user_input_file)
logging.info("Blastp done")
logging.info(f"Blastp out: {bout}")
logging.info(f"Blastp err: {berr}")
# 1.2 Obtaining Results
result_dict = getBlastpResult(nowUuid, user_input_file)
logging.info("Get blastp result done")
logging.info(f"Result dict: {result_dict}")
# 1.3 Find the corresponding sequence and align the sequence
outToken = searchDB(input_query=result_dict)
logging.info("Search DB done")
logging.info(f"Out token: {outToken}")
# 1.4 Alignment sequence returns Aligntoken
alignToken = AlignSeqs(user_input_file, outToken)
logging.info("Align seqs done")
logging.info(f"Align token: {alignToken}")
# 1.5 Perform sliding window calculation and return a smooth value of 0-1
smoothToken = smoon(alignToken)
logging.info("Smooth token done")
logging.info(f"Smooth token: {smoothToken}")
logging.info("🔔🔔 RareCodon done 🔔🔔")
for k, v in smoothToken.items():
    if len(v[-1]) != len(v[-2]):
        print("error", k)
# 2.  Perform codon optimization
logging.info("🔔🔔 Model process start 🔔🔔")
doneFilePath = modelProcess(
    nowUuid=nowUuid,
    fastafile=user_input_file,
    smoothTokenDict=smoothToken,
    generatenum=1,
    output_path=args.output_file,
    logging=logging,
)
logging.info("🔔🔔 Model process done 🔔🔔")


# 3.  Return data
logging.info(f"🔔🔔 {doneFilePath} plz check out 🔔🔔")
logging.info("🔔🔔🔔 All done 🔔🔔🔔")
