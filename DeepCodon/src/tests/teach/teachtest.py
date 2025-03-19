import random

from accelerate import Accelerator
from torch.utils.data import Subset

from DeepCodon.config.config import *
from DeepCodon.src.generate.convert import index2test
from DeepCodon.src.generate.timestamp import get_timestamp
from DeepCodon.src.train.datasets import *
from DeepCodon.src.train.model import TransformerTranslator

# Only used for generating tests
whichmodel = "mymodel"
whichdata = "val"

accelerator = Accelerator()
device = accelerator.device
nowtime = get_timestamp()


enc_inputs_val = torch.load(enc_inputs_val_path)
dec_inputs_val = torch.load(dec_inputs_val_path)
dec_outputs_val = torch.load(dec_outputs_val_path)

print("model path", whichmodel2generate)
print(
    "which data select", enc_inputs_val_path, dec_inputs_val_path, dec_outputs_val_path
)
val_data = MyDataSet(enc_inputs_val, dec_inputs_val, dec_outputs_val)

total_samples = len(val_data)
indices = list(range(total_samples))
random.seed(42)
random_indices = random.sample(indices, 1000)
val_data_subset = Subset(val_data, random_indices)
val_loader = Data.DataLoader(val_data_subset, batch_size=1, shuffle=False)

generate_num = "all"

model = TransformerTranslator()
model, eval_dataloader = accelerator.prepare(model, val_loader)
model = accelerator.unwrap_model(model)
model.load_state_dict(torch.load(whichmodel2generate))

translations = []
for number, (enc_inputs, dec_inputs, dec_outputs) in tqdm(enumerate(val_loader)):
    if generate_num == "all":
        model.eval()
        enc_inputs = enc_inputs.cuda()
        dec_inputs = dec_inputs.cuda()
        dec_outputs = dec_outputs.cuda()
        out = model(enc_inputs, dec_inputs)
        predict = torch.argmax(out, dim=1)
        output_str = index2test(enc_inputs, dec_outputs, predict)
        translations.append(f"{number},{output_str}")
    else:
        if number < generate_num:
            model.eval()
            enc_inputs = enc_inputs.cuda()
            dec_inputs = dec_inputs.cuda()
            dec_outputs = dec_outputs.cuda()
            out = model(enc_inputs, dec_inputs)
            predict = torch.argmax(out, dim=1)
            output_str = index2test(enc_inputs, dec_outputs, predict)
            translations.append(f"{number},{output_str}")
        else:
            break

print("finish genereate")
time = 0
with open(
    f"data/test_data/teach/{whichmodel}_{whichdata}.csv",
    "w",
    encoding="utf-8",
) as f:
    for translation in translations:
        num, inputseq, outputseq, inputgc, outputgc = translation.split(",")
        f.write(translation + "\n")
print("finish write")
