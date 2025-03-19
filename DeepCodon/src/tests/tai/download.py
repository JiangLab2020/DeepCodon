import os
import re

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from DeepCodon.config.config import *
from DeepCodon.src.train.datasets import *


class ProcessData:
    def __init__(self, filepath, inputCol, outputCol):
        self.dfi, self.dfo = self.init_data(filepath, inputCol, outputCol)
        self.__len__ = len(self.dfi)
        self.link = self.init_link()

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

    def init_data(self, filepath, inputCol, outputCol):
        df = pd.read_csv(filepath, header=None, index_col=None)
        return df[inputCol], df[outputCol]

    def checkAndDownload(self, index: int, inputlink: str):
        inputname = inputlink.split("/")[-1]
        if not os.path.exists(downloadPath + "/" + inputname + "_cds_from_genomic.fna"):
            print(f"{index} cds_from_genomic not exists")
            self.download(
                type="_cds_from_genomic", inputpath=inputlink, inputname=inputname
            )
        else:
            tqdm.write(f"{index} cds_from_genomic exists")
        if not os.path.exists(downloadPath + "/" + inputname + "_genomic.fna"):
            self.download(type="_genomic", inputpath=inputlink, inputname=inputname)
        else:
            tqdm.write(f"{index} genomic exists")

    def download(self, type, inputpath, inputname):
        with open("data/download_url/download.txt", "a") as dfile:
            if type == "_cds_from_genomic":
                dfile.write(
                    inputpath + "/" + inputname + "_cds_from_genomic.fna.gz" + "\n"
                )
            elif type == "_genomic":
                dfile.write(inputpath + "/" + inputname + "_genomic.fna.gz" + "\n")

    def find_link(self, inputstr):
        matches = re.findall(r'(?<=http=")(.*?)(?=")', self.link[inputstr])
        return matches[0]

    def processing(self, outpath):
        for (index1, linei), (index2, lineo) in zip(self.dfi.items(), self.dfo.items()):
            outputlink = self.find_link(linei)
            self.checkAndDownload(index1, outputlink)


whichmodel = "mymodel"
whichdata = "val"

process = ProcessData(
    filepath=f"data/test_data/teach/{whichmodel}_{whichdata}.csv",
    inputCol=1,
    outputCol=2,
)
process.processing(outpath=f"data/test_data/tai/{whichmodel}_{whichdata}.csv")
