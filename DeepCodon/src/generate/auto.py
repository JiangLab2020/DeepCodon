# for generate no rarecodon
import argparse

from accelerate import Accelerator
from Bio import SeqIO
from Bio.Seq import Seq

from DeepCodon.config.config import *
from DeepCodon.src.generate.beam import beam_search
from DeepCodon.src.generate.convert import index2content
from DeepCodon.src.generate.greedy import greedy_algorithm
from DeepCodon.src.generate.mixture import mixture_algorithm
from DeepCodon.src.generate.temperature import temperature_algorithm
from DeepCodon.src.generate.timestamp import get_timestamp
from DeepCodon.src.generate.top_k import top_k_algorithm
from DeepCodon.src.generate.top_p import top_p_algorithm
from DeepCodon.src.tests.tai.centerCont import *
from DeepCodon.src.train.datasets import *
from DeepCodon.src.train.model import TransformerTranslator

accelerator = Accelerator()
device = accelerator.device

parser = argparse.ArgumentParser(description="1.0")

parser.add_argument("-f", "--file", type=str, help="file name")
parser.add_argument("-n", "--number", type=int, help="number of iterations")

args = parser.parse_args()
nowtime = get_timestamp()
inputdata = args.file.split("/")[-1].split(".")[0]
print("Reading file", inputdata)
with open(f"data/user_data/cache/{nowtime}.tsv", "w") as tsv_file:
    for record in tqdm(SeqIO.parse(args.file, "fasta")):
        # dna or protein
        inputseq = str(record.seq)
        if "M" not in inputseq:
            whichdata = "codon"
            if len(inputseq) % 3 != 0:
                continue
            else:
                protein_sequence = str(Seq(inputseq).translate())
        else:
            whichdata = "protein"
            protein_sequence = inputseq
        if len(protein_sequence) <= pclen and len(protein_sequence) >= minlen:
            if args.number > 1:
                for i in range(args.number):
                    tsv_file.write(f"{record.id}_{i}\t{protein_sequence}\t{inputseq}\n")
            else:
                tsv_file.write(f"{record.id}\t{protein_sequence}\t{inputseq}\n")

record_ids, enc_inputs, inputseqs, opmask_inputs = make_data_user(
    f"data/user_data/cache/{nowtime}.tsv"
)
user_data = MyDataSetuser(record_ids, enc_inputs, inputseqs, opmask_inputs)
user_loader = Data.DataLoader(user_data, shuffle=False)

model = TransformerTranslator()
model, eval_dataloader = accelerator.prepare(model, user_loader)
model = accelerator.unwrap_model(model)
model.load_state_dict(torch.load(whichmodel2generate))

translations = []
for number, (record_ids, enc_inputs, inputseqs, opmask_inputs) in tqdm(
    enumerate(user_loader)
):
    if number < maxgenerateseq:
        # output
        enc_inputs = enc_inputs.cuda()
        opmask_inputs = opmask_inputs.cuda()
        if generationmethod == "greedy":
            predict = greedy_algorithm(model, device, enc_inputs, opmask_inputs)
        elif generationmethod == "beam":
            predict = beam_search(model, device, enc_inputs, opmask_inputs)
        elif generationmethod == "top_p":
            predict = top_p_algorithm(model, device, enc_inputs)
        elif generationmethod == "top_k":
            predict = top_k_algorithm(model, device, enc_inputs)
        elif generationmethod == "temperature":
            predict = temperature_algorithm(model, device, enc_inputs)
        elif generationmethod == "mix":
            predict = mixture_algorithm(
                model,
                device,
                enc_inputs,
                opmask_inputs,
                f"data/user_data/output/debug/{number}",
            )
        else:
            break
        output_str = index2content(
            enc_inputs, predict, inputseqs, record_ids, whichdata
        )
        translations.append(output_str)
    else:
        break

print("finish genereate")
with open(
    f"data/user_data/output/{inputdata}_{nowtime}_{ifcase}_{generationmethod}.csv",
    "w",
    encoding="utf-8",
) as f:
    if whichdata == "codon":
        f.write(
            "id,inputseq,outputseq,inputMFE,outputMFE,inputGC,outputGC,inputCAI,outputCAI,istrue\n"
        )
    else:
        f.write("id,inputseq,outputseq,outputMFE,outputGC,outputCAI,istrue\n")
    for translation in translations:
        f.write(translation + "\n")

if whichdata == "codon":
    process = ProcessData(
        filepath=f"data/user_data/output/{inputdata}_{nowtime}_{ifcase}_{generationmethod}.csv",
        inputCol=1,
        outputCol=2,
    )
    process.processing(
        outpath=f"data/user_data/cache/{inputdata}_{nowtime}_{ifcase}_{generationmethod}.csv"
    )
    dfnotai = pd.read_csv(
        f"data/user_data/output/{inputdata}_{nowtime}_{ifcase}_{generationmethod}.csv"
    )
    dftai = pd.read_csv(
        f"data/user_data/cache/{inputdata}_{nowtime}_{ifcase}_{generationmethod}.csv",
    )
    dftai = dftai.T
    dftai.reset_index(drop=True, inplace=True)
    dftai.drop(index=0, inplace=True)
    dftai.reset_index(drop=True, inplace=True)
    dftai.columns = ["inputTAI", "outputTAI"]
    df = pd.concat([dfnotai, dftai], axis=1)
    df.to_csv(
        f"data/user_data/output/{inputdata}_{nowtime}_{ifcase}_{generationmethod}.csv",
        index=False,
    )
