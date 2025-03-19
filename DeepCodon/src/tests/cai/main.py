import os

from DeepCodon.src.tests.cai.cai import processing

current_path = os.getcwd()

whichmodel = "idt1000"
whichdata = "val"
refer_sheet = os.path.join(current_path, "src", "tests", "cai", "ecoli.txt")

input_path = os.path.join(
    current_path,
    "data",
    "test_data",
    "teach",
    f"{whichmodel}_{whichdata}.csv",
)
output_path = os.path.join(
    current_path, "data", "test_data", "cai", f"{whichmodel}_{whichdata}.csv"
)


with open(output_path, "w") as f:
    with open(input_path, "r") as fileinput:
        for line in fileinput:
            if line.startswith("id"):
                continue
            if whichdata == "val":
                id, inputseq, outputseq, *_ = line.strip().split(",")
                try:
                    inputcai = processing(infile=refer_sheet, inseq=inputseq)
                    outputcai = processing(infile=refer_sheet, inseq=outputseq)
                except:
                    print(f"Error: {inputseq}, {outputseq}")
                f.write(f"{id},{inputseq},{outputseq},{inputcai},{outputcai}\n")
            else:
                id, inputc, inputgc = line.strip().split(",")
                try:
                    inputcai = processing(infile=refer_sheet, inseq=inputc)
                except:
                    print(f"Error: {inputc}")
                f.write(f"{id},{inputc},{inputcai}\n")
