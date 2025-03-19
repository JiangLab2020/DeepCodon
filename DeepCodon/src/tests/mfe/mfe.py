from multiprocessing import Pool

import ViennaRNA
from tqdm import tqdm

whichmodel = "codonop"
whichdata = "val"


def process(line):
    line = line.strip()
    if whichdata == "all":
        (id, inputseq, inputgc) = line.split(",")
        rna_structure_input = ViennaRNA.fold(inputseq)
        return (
            id,
            inputseq,
            rna_structure_input,
        )
    else:
        (id, inputseq, outseq, inputgc, outputgc) = line.split(",")
        rna_structure_input = ViennaRNA.fold(inputseq)
        rna_structure_output = ViennaRNA.fold(outseq)
        return id, inputseq, outseq, rna_structure_input, rna_structure_output


with open(
    f"data/test_data/mfe/{whichmodel}_{whichdata}.csv", "w", encoding="utf-8"
) as fw:
    with open(
        f"data/test_data/teach/{whichmodel}_{whichdata}.csv",
        "r",
        encoding="utf-8",
    ) as f:
        with Pool() as pool:
            results = []
            for line in tqdm(f):
                results.append(
                    pool.apply_async(
                        process,
                        args=(line,),
                    )
                )
            pool.close()
            pool.join()
            for result in results:
                if whichdata == "all":
                    id, inputseq, rna_structure_input = result.get()
                    fw.write(
                        f"{id},{inputseq},{rna_structure_input[1]},{rna_structure_input[0]}\n"
                    )
                else:
                    id, inputseq, outseq, rna_structure_input, rna_structure_output = (
                        result.get()
                    )
                    fw.write(
                        f"{id},{inputseq},{outseq},{rna_structure_input[1]},{rna_structure_input[0]},{rna_structure_output[1]},{rna_structure_output[0]}\n"
                    )
print("finish")
