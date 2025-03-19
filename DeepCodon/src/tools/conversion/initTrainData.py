from sklearn.model_selection import train_test_split

from DeepCodon.config.config import (
    dec_inputs_test_path,
    dec_inputs_train_path,
    dec_inputs_val_path,
    dec_outputs_test_path,
    dec_outputs_train_path,
    dec_outputs_val_path,
    enc_inputs_test_path,
    enc_inputs_train_path,
    enc_inputs_val_path,
)
from DeepCodon.src.train.datasets import *


enc_inputs, dec_inputs, dec_outputs = make_data()
print(f"have read {len(enc_inputs)} data")

(
    enc_inputs_train,
    enc_inputs_test,
    dec_inputs_train,
    dec_inputs_test,
    dec_outputs_train,
    dec_outputs_test,
) = train_test_split(
    enc_inputs, dec_inputs, dec_outputs, test_size=0.2, random_state=666
)
(
    enc_inputs_test,
    enc_inputs_val,
    dec_inputs_test,
    dec_inputs_val,
    dec_outputs_test,
    dec_outputs_val,
) = train_test_split(
    enc_inputs_test,
    dec_inputs_test,
    dec_outputs_test,
    test_size=0.5,
    random_state=666,
)

torch.save(enc_inputs_train, enc_inputs_train_path)
torch.save(enc_inputs_test, enc_inputs_test_path)
torch.save(enc_inputs_val, enc_inputs_val_path)

torch.save(dec_inputs_train, dec_inputs_train_path)
torch.save(dec_inputs_test, dec_inputs_test_path)
torch.save(dec_inputs_val, dec_inputs_val_path)

torch.save(dec_outputs_train, dec_outputs_train_path)
torch.save(dec_outputs_test, dec_outputs_test_path)
torch.save(dec_outputs_val, dec_outputs_val_path)
