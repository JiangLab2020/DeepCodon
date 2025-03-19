import os

from centerCont import projectPath, refername

from DeepCodon.config.config import whichtai


def inner_taicontrol(CachePath: str):
    innerpath = f"{CachePath}/tRNAscan.out"
    outputpath = f"{CachePath}/tRNAscan.trna"
    codon_list = [
        "TTT",
        "TTC",
        "TTA",
        "TTG",
        "TCT",
        "TCC",
        "TCA",
        "TCG",
        "TAT",
        "TAC",
        "TAT",
        "TAG",
        "TGT",
        "TGC",
        "TGA",
        "TGG",
        "CTT",
        "CTC",
        "CTA",
        "CTG",
        "CCT",
        "CCC",
        "CCA",
        "CCG",
        "CAT",
        "CAC",
        "CAA",
        "CAG",
        "CGT",
        "CGC",
        "CGA",
        "CGG",
        "ATT",
        "ATC",
        "ATA",
        "ATG",
        "ACT",
        "ACC",
        "ACA",
        "ACG",
        "AAT",
        "AAC",
        "AAA",
        "AAG",
        "AGT",
        "AGC",
        "AGA",
        "AGG",
        "GTT",
        "GTC",
        "GTA",
        "GTG",
        "GCT",
        "GCC",
        "GCA",
        "GCG",
        "GAT",
        "GAC",
        "GAA",
        "GAG",
        "GGT",
        "GGC",
        "GGA",
        "GGG",
    ]
    anticodon_comp_list = []
    for codon in codon_list:
        anticodon = (
            codon.replace("A", "1")
            .replace("T", "2")
            .replace("C", "3")
            .replace("G", "4")
            .replace("1", "T")
            .replace("2", "A")
            .replace("3", "G")
            .replace("4", "C")
        )
        anticodon_comp_list.append(anticodon[::-1])

    fh = open(innerpath)
    trna_copy_num_dict = {}
    fh.readline()
    fh.readline()
    fh.readline()
    for line in fh:
        trna_type = line.split("\t")[5]
        note = line.split("\t")[9].strip()
        if trna_type != "NNN" and note != "pseudo":
            if trna_type not in trna_copy_num_dict.keys():
                trna_copy_num_dict[trna_type] = 1
            else:
                trna_copy_num_dict[trna_type] += 1
    fh.close()

    fh1 = open(outputpath, "w")
    for anticodon_comp in anticodon_comp_list:
        if anticodon_comp not in trna_copy_num_dict.keys():
            fh1.write(str(0) + "\n")
        else:
            fh1.write(str(trna_copy_num_dict[anticodon_comp]) + "\n")
    fh1.close()


CachePath = "data/bintest/-1"
inputcds = "data/raw_data/final_tai_top0.1h.fasta"
# inputcds = "data/bintest/cds_from_genomic.fna"
downloadPath = "data/download_data"
inputgenomic = f"{downloadPath}/{refername}_genomic.fna"
print(f"choiced model {whichtai}")
os.system(f"rm -rf {CachePath}/*")
os.system(f"perl {projectPath}/codonM {inputcds} {CachePath}/output_codonM.m")


# os.system("cp data/bintest/eschColi_K_12_MG1655-tRNAs.out data/bintest/-1/tRNAscan.out")
# os.system("cp data/bintest/author.trna data/bintest/-1/tRNAscan.trna")
os.system(f"tRNAscan-SE -G -o {CachePath}/tRNAscan.out {inputgenomic}")
inner_taicontrol(CachePath=CachePath)
os.system(f"Rscript {projectPath}/run_tai.R {CachePath}")
