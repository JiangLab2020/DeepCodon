import datetime
from pathlib import Path

# tensorboard
now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
tensorboardUpdateTime = 1
batchUpdateTime = 100

# path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
checkpointPath = PROJECT_ROOT / "DeepCodon" / "data" / "my_checkpoint"
tensorboardPath = PROJECT_ROOT / "DeepCodon" / "data" / "log" / "tensorboard" / now
modelSavePath = PROJECT_ROOT / "DeepCodon" / "model"
pretrainModelPath = PROJECT_ROOT / "DeepCodon" / "model" / "pretrainingModel"
fineTuneModelPath = PROJECT_ROOT / "DeepCodon" / "model" / "fine-tuningModel"
# dataset
downloadPath = PROJECT_ROOT / "DeepCodon" / "data" / "datasets" / "download_data"
raw_input_file = PROJECT_ROOT / "EPCD" / "finalData"
raw_output_file = PROJECT_ROOT / "DeepCodon" / "data" / "datasets" / "raw_data"
traindataPath = PROJECT_ROOT / "DeepCodon" / "data" / "datasets" / "train_data"
# -------------training-----------------
# checkpoint
useCheckPoint = True
checkpointPathUpdateTime = 1

# model
batch = 6
minlen = 50
pclen = 1000
d_model = 256
d_ff = 64
n_layers = 6
n_heads = 8
maxEpoch = 20000

# optimizer
learningRate = 0.4
decayFactor = 0.9
epsilon = 1e-7
weightDecay = 1e-5

# data
datachoice = "normal"
if datachoice == "normal":
    raw_input_file = raw_input_file / "final_tai.fasta"
    raw_output_file = raw_output_file / "normal.tsv"
    modelSavePath = pretrainModelPath
else:
    raw_input_file = raw_input_file / "final_tai_top0.1.fasta"
    raw_output_file = raw_output_file / "hightTAI.tsv"
    modelSavePath = fineTuneModelPath

enc_inputs_train_path = traindataPath / f"enc_inputs_train_{datachoice}.pth"
enc_inputs_test_path = traindataPath / f"enc_inputs_test_{datachoice}.pth"
enc_inputs_val_path = traindataPath / f"enc_inputs_val_{datachoice}.pth"

dec_inputs_train_path = traindataPath / f"dec_inputs_train_{datachoice}.pth"
dec_inputs_test_path = traindataPath / f"dec_inputs_test_{datachoice}.pth"
dec_inputs_val_path = traindataPath / f"dec_inputs_val_{datachoice}.pth"

dec_outputs_train_path = traindataPath / f"dec_outputs_train_{datachoice}.pth"
dec_outputs_test_path = traindataPath / f"dec_outputs_test_{datachoice}.pth"
dec_outputs_val_path = traindataPath / f"dec_outpusts_val_{datachoice}.pth"


# -------------generating-----------------
maxModelnumber = 2000
# setting 650 680 10
# whichmodel2generate = "model/epoch1/pytorch_model.bin"
# whichmodel2generate = "model/epoch600/pytorch_model.bin"
whichmodel2generate = fineTuneModelPath / "epoch680/pytorch_model.bin"
maxgenerateseq = 1000000000

# mix greedy beam
generationmethod = "mix"
beam_width = 5
ifcase = "nocase"

# tai ecoli self refer
whichtai = "ecoli"

# analysis
compare_file_path = (
    PROJECT_ROOT / "DeepCodon" / "data" / "datasets" / "test_data" / "compare"
)
