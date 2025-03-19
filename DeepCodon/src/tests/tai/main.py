from DeepCodon.src.tests.tai.centerCont import *

whichmodel = "idt1000"
whichdata = "val"


process = ProcessData(
    filepath=f"data/test_data/teach/{whichmodel}_{whichdata}.csv",
    inputCol=1,
    outputCol=2,
)
process.processing(outpath=f"data/test_data/tai/{whichmodel}_{whichdata}.csv")
