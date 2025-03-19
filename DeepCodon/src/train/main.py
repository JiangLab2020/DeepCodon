import os

import setproctitle
import torch
import torch.nn as nn
import torch.optim as optim
from accelerate import Accelerator, DistributedDataParallelKwargs
from torch.utils.tensorboard import SummaryWriter

from DeepCodon.config.config import *
from DeepCodon.src.train.datasets import Data, MyDataSet
from DeepCodon.src.train.evaluate import evaluate
from DeepCodon.src.train.model import TransformerTranslator
from DeepCodon.src.train.model_save import epoch_log, kmodel, read_epoch

ddp_kwargs = DistributedDataParallelKwargs(find_unused_parameters=True)
accelerator = Accelerator(kwargs_handlers=[ddp_kwargs])
device = accelerator.device

if accelerator.is_main_process:
    log_dir = tensorboardPath
    writer = SummaryWriter(log_dir)
setproctitle.setproctitle("starting train")


enc_inputs_train = torch.load(enc_inputs_train_path)
dec_inputs_train = torch.load(dec_inputs_train_path)
dec_outputs_train = torch.load(dec_outputs_train_path)
enc_inputs_test = torch.load(enc_inputs_test_path)
dec_inputs_test = torch.load(dec_inputs_test_path)
dec_outputs_test = torch.load(dec_outputs_test_path)

train_data = MyDataSet(enc_inputs_train, dec_inputs_train, dec_outputs_train)
test_data = MyDataSet(enc_inputs_test, dec_inputs_test, dec_outputs_test)

train_loader = Data.DataLoader(train_data, batch_size=batch, shuffle=True)
test_loader = Data.DataLoader(test_data, batch_size=batch, shuffle=True)
accelerator.print("Data read completed")

model = TransformerTranslator()

criterion = nn.CrossEntropyLoss()
optimizer = optim.Adadelta(
    model.parameters(),
    lr=learningRate,
    rho=decayFactor,
    eps=epsilon,
    weight_decay=weightDecay,
)
model, optimizer, train_loader, test_loader = accelerator.prepare(
    model, optimizer, train_loader, test_loader
)

if useCheckPoint:
    accelerator.print("Start reading checkpoints")
    try:
        oldepoch = read_epoch(folder_path=checkpointPath)
        accelerator.load_state(input_dir=checkpointPath)
    except FileNotFoundError:
        accelerator.print("The folder does not exist")
    except Exception as e:
        accelerator.print(f"An error has occurred: {str(e)}")
else:
    oldepoch = 0
    accelerator.print("Skip reading checkpoints")

for epoch in range(oldepoch, oldepoch + maxEpoch):
    for Step, (enc_inputs, dec_inputs, dec_outputs) in enumerate(train_loader):
        enc_inputs, dec_inputs, dec_outputs = (
            enc_inputs.cuda(),
            dec_inputs.cuda(),
            dec_outputs.cuda(),
        )
        outputs = model(enc_inputs, dec_inputs)
        loss = criterion(outputs, dec_outputs.view(-1))
        accelerator.backward(loss)
        optimizer.step()
        optimizer.zero_grad()
        # batch loss
        if (Step % batchUpdateTime) == 0:
            if accelerator.is_main_process:
                writer.add_scalar(
                    "batch_loss", float(loss.item()), epoch * len(train_loader) + Step
                )
    # checkpoint
    if useCheckPoint and (epoch % checkpointPathUpdateTime) == 0:
        accelerator.wait_for_everyone()
        accelerator.save_state(output_dir=checkpointPath)
        if accelerator.is_main_process:
            epoch_log(epoch, folder_path=checkpointPath)
    # tensorboard
    if (epoch % tensorboardUpdateTime) == 0:
        train_loss, train_acc = evaluate(model, train_loader, criterion, epoch)
        test_loss, test_acc = evaluate(model, test_loader, criterion, epoch)
        accelerator.print(
            "Epoch:",
            epoch,
            "batch_loss =",
            "{:.6f}".format(loss),
            "train_loss =",
            "{:.6f}".format(train_loss),
            "train_acc =",
            "{:.6f}".format(train_acc),
            "test_loss =",
            "{:.6f}".format(test_loss),
            "test_acc =",
            "{:.6f}".format(test_acc),
        )
        if accelerator.is_main_process:
            writer.add_scalar("train_loss", float(train_loss), epoch)
            writer.add_scalar("train_acc", float(train_acc), epoch)
            writer.add_scalar("test_loss", float(test_loss), epoch)
            writer.add_scalar("test_acc", float(test_acc), epoch)
        setproctitle.setproctitle(f"TaskProgress {int(train_acc * 100)}%")
        # deletemodel
        if accelerator.is_main_process:
            kmodel(k=maxModelnumber, folder_path=modelSavePath)
        accelerator.wait_for_everyone()
        accelerator.save_model(
            model=model, save_directory=f"{str(modelSavePath)}/epoch{epoch}"
        )

accelerator.wait_for_everyone()
accelerator.save_model(model=model, save_directory=f"{str(modelSavePath)}/final_model")
print("finish")
if accelerator.is_main_process:
    writer.close()
