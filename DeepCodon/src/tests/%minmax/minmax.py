from CodonData import download_codon_frequencies_from_kazusa
from CodonEvaluation import get_min_max_percentage

whichmodel = "test"
whichdata = "val"

codonfreq = download_codon_frequencies_from_kazusa(taxonomy_id=83333)
print(codonfreq)


def process(line):
    line = line.strip()
    if whichdata == "all":
        (id, inputseq, inputgc) = line.split(",")
        exit()
    else:
        (id, inputseq, outseq, inputgc, outputgc) = line.split(",")
        minmax_input = get_min_max_percentage(inputseq, codonfreq)
        minmax_output = get_min_max_percentage(outseq, codonfreq)
        print("len inputseq:", len(inputseq))
        print("len inputminmax:", len(minmax_input))
        return id, inputseq, outseq, minmax_input, minmax_output


with open(
    f"data/test_data/minmax/{whichmodel}_{whichdata}.csv", "w", encoding="utf-8"
) as fw:
    with open(
        f"data/test_data/teach/{whichmodel}_{whichdata}.csv",
        "r",
        encoding="utf-8",
    ) as f:
        for line in f:
            if whichdata == "all":
                id, inputseq, minmax = process(line)
                fw.write(f"{id}@{inputseq}@{minmax}\n")
            else:
                id, inputseq, outseq, minmax_input, minmax_output = process(line)
                fw.write(f"{id}@{inputseq}@{outseq}@{minmax_input}@{minmax_output}\n")
print("finish")
