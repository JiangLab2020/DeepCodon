from DeepCodon.src.tools.generate.gc import calculate_gc

whichmodel = "genescript1000"
whichdata = "val"

fw = open(f"data/test_data/teach/{whichmodel}_{whichdata}_o.csv", "w", encoding="utf-8")
f = open(
    f"data/test_data/teach/{whichmodel}_{whichdata}.csv",
    "r",
    encoding="utf-8",
)


for line in f:
    line = line.strip()
    id, inputseq, outputseq, *_ = line.split(",")
    fw.write(
        f"{id},{inputseq},{outputseq},{calculate_gc(inputseq)},{calculate_gc(outputseq)}\n"
    )
print("finish")
