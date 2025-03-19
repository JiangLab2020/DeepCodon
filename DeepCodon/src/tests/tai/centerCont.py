import os
import re
import time
from multiprocessing import Pool

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from DeepCodon.config.config import *
from DeepCodon.config.config import whichtai
from DeepCodon.config.server_config import tai_cache_path
from DeepCodon.src.train.datasets import *

# 1.  Obtain dataset
# 2.  Download the corresponding CDS and genomic files
# 3.  Operating Data
# -------------------param----------------
innerCachePath = tai_cache_path
downloadPath = str(downloadPath)
projectPath = "src/tests/tai"
refername = "GCF_000005845.2_ASM584v2"
referlink = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"
)
# --------------------------------------


class ProcessData:
    def __init__(self, filepath, inputCol, outputCol):
        self.dfi, self.dfo = self.init_data(filepath, inputCol, outputCol)
        self.__len__ = len(self.dfi)
        self.link = self.init_link()
        self.whichtai = whichtai

    def init_data(self, filepath, inputCol, outputCol):
        df = pd.read_csv(filepath, header=None, index_col=None)
        return df[inputCol], df[outputCol]

    def init_link(self):
        tqdm.write("start init link")
        link = dict()
        with open(raw_input_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequenceSeq = str(record.seq)
                sequenceDes = record.description
                link[sequenceSeq] = sequenceDes
        tqdm.write("finish init")
        return link

    def find_link(self, inputstr):
        matches = re.findall(r'(?<=http=")(.*?)(?=")', self.link[inputstr])
        return matches[0]

    def download(self, type, inputpath, inputname):
        if type == "_cds_from_genomic":
            os.system(
                f"wget -P {downloadPath} "
                + inputpath
                + "/"
                + inputname
                + "_cds_from_genomic.fna.gz"
            )
            #
            os.system(
                "gunzip " + downloadPath + "/" + inputname + "_cds_from_genomic.fna.gz"
            )
            #
            os.system(
                "rm " + downloadPath + "/" + inputname + "_cds_from_genomic.fna.gz"
            )
        elif type == "_genomic":
            os.system(
                f"wget -P {downloadPath} "
                + inputpath
                + "/"
                + inputname
                + "_genomic.fna.gz"
            )
            #
            os.system("gunzip " + downloadPath + "/" + inputname + "_genomic.fna.gz")
            #
            os.system("rm " + downloadPath + "/" + inputname + "_genomic.fna.gz")

    def checkAndDownload(self, index: int = None, inputlink: str = None):
        """
        Check if there is any downloaded data, if not, download it
        Parameters:
        index (int):  Index of elements
        inputlink (str):  input link
        """
        if whichtai == "self":
            inputname = inputlink.split("/")[-1]
            if not os.path.exists(
                downloadPath + "/" + inputname + "_cds_from_genomic.fna"
            ):
                self.download(
                    type="_cds_from_genomic", inputpath=inputlink, inputname=inputname
                )
            else:
                tqdm.write(f"{index} cds_from_genomic exists")
            if not os.path.exists(downloadPath + "/" + inputname + "_genomic.fna"):
                self.download(type="_genomic", inputpath=inputlink, inputname=inputname)
            else:
                tqdm.write(f"{index} genomic exists")
        elif whichtai == "ecoli":
            if not os.path.exists(downloadPath + "/" + refername + "_genomic.fna"):
                self.download(type="_genomic", inputpath=referlink, inputname=refername)
            else:
                tqdm.write("refer genomic exists")
        else:
            raise ValueError("whichtai must be self or ecoli")

    def taicontrol(
        self,
        inputcds: str,
        inputgenomic: str,
        innerpath: str,
    ):
        """
        Calculate Tai value
        Parameters:
        Inputcds (str): CDS file path
        Inputgenomic (str): genomic file path
        """
        os.system(f"perl {projectPath}/codonM {inputcds} {innerpath}/output_codonM.m")
        os.system(f"tRNAscan-SE -G -o {innerpath}/tRNAscan.out {inputgenomic}")
        self.inner_taicontrol(CachePath=innerpath)
        os.system(f"Rscript {projectPath}/run_tai.R {innerpath}")

    def inner_taicontrol(self, CachePath: str):
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

    def inner_file_control(
        self,
        filepath: str,
        modelinput: str,
        modeloutput: str,
        outpath,
    ):
        """
        Used to modify the order of fasta files for easy operation
        The sequence order is the input and output of the model, and other sequences are included

        Parameters:
        filepath (str):  File Path
        modelinput (str):  Model input
        modeloutput (str):  Model output
        outpath (str):  Output Path
        """
        filepathname = filepath.split("/")[-1]
        with open(filepath, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            with open(f"{outpath}/{filepathname}", "w") as outhandle:
                outhandle.write(f">input\n{modelinput}\n")
                outhandle.write(f">output\n{modeloutput}\n")
                for record in records:
                    if str(record.seq) == modelinput or str(record.seq) == modeloutput:
                        pass
                    else:
                        outhandle.write(f">{record.id}\n{str(record.seq)}\n")
        return f"{outpath}/{filepathname}"

    def ecoli_file_control(self, cachePath):
        outputfasta = f"{cachePath}/eocli.fasta"
        with open(outputfasta, "w") as handle:
            for (index1, linei), (index2, lineo) in zip(
                self.dfi.items(), self.dfo.items()
            ):
                handle.write(f">input{index1}\n{linei}\n")
                handle.write(f">output{index1}\n{lineo}\n")
        return outputfasta

    def process_item(self, index1, linei, lineo, innerCachePath):
        outData = []
        outputlink = self.find_link(linei)
        name = outputlink.split("/")[-1]
        self.checkAndDownload(index1, outputlink)
        CachePath = f"{innerCachePath}/{name}"
        if not os.path.exists(CachePath):
            os.makedirs(CachePath)
        print("CachePath", CachePath)
        # os.system(f"rm -rf {CachePath}/*")
        try:
            inputCdsPath = self.inner_file_control(
                filepath=f"{downloadPath}/{name}_cds_from_genomic.fna",
                modelinput=linei,
                modeloutput=lineo,
                outpath=CachePath,
            )
            # run
            self.taicontrol(
                inputcds=inputCdsPath,
                inputgenomic=f"{downloadPath}/{name}_genomic.fna",
                innerpath=CachePath,
            )
            # read
            with open(f"{CachePath}/output.tai", "r") as handle:
                outData = [line.strip() for line in handle.readlines()]
            return (f"{index1},{outputlink}", outData)
        except Exception as e:
            with open(f"{projectPath}/error", "a") as efile:
                efile.write(f"{e},{index1},{outputlink}\n")

    def process_item_ecoli(self):
        self.checkAndDownload()
        if not os.path.exists(innerCachePath):
            os.makedirs(innerCachePath)
        print("CachePath", innerCachePath)
        # os.system(f"rm -rf {innerCachePath}/*")
        CachePath = f"{innerCachePath}/{(time.time())}"
        if not os.path.exists(CachePath):
            os.makedirs(CachePath)
        theecoliinputcds = self.ecoli_file_control(CachePath)
        self.taicontrol(
            inputcds=theecoliinputcds,
            inputgenomic=f"{downloadPath}/{refername}_genomic.fna",
            innerpath=CachePath,
        )

        with open(f"{CachePath}/output.tai", "r") as handle:
            outData = [line.strip() for line in handle.readlines()]
            outDict = {
                i: [outData[i * 2], outData[i * 2 + 1]]
                for i in range(len(outData) // 2)
            }
        return outDict

    def processing(self, outpath):
        if self.whichtai == "self":
            outDataDic = dict()
            with Pool() as pool:
                results = []
                for (index1, linei), (index2, lineo) in zip(
                    self.dfi.items(), self.dfo.items()
                ):
                    results.append(
                        pool.apply_async(
                            self.process_item,
                            args=(index1, linei, lineo, f"{innerCachePath}/{index1}"),
                        )
                    )
                pool.close()
                pool.join()

                for result in results:
                    try:
                        index, outData = result.get()
                        outDataDic[index] = outData
                    except Exception as e:
                        print(e)

            dfout = pd.DataFrame.from_dict(outDataDic, orient="index")
            dfout = dfout.transpose()
            dfout.to_csv(outpath, index=False)
        elif whichtai == "ecoli":
            outDataDic = self.process_item_ecoli()
            dfout = pd.DataFrame.from_dict(outDataDic, orient="index")
            dfout = dfout.transpose()
            dfout.to_csv(outpath, index=False)
        else:
            raise ValueError("whichtai must be self or ecoli")
