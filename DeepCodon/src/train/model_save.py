import os
import re
import shutil


def kmodel(k, folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    files = os.listdir(folder_path)
    model_files = [
        (f, int(re.search(r"\d+", f).group())) for f in files if re.search(r"\d+", f)
    ]
    model_files.sort(key=lambda x: x[1], reverse=True)
    try:
        for f, _ in model_files[k:]:
            folder_to_delete = os.path.join(folder_path, f)
            shutil.rmtree(folder_to_delete)
    except FileNotFoundError:
        print(f"The folder {folder_to_delete} does not exist and cannot be deleted.")
    except Exception as e:
        print(f"An error has occurred: {str(e)}")


def epoch_log(epoch, folder_path):
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
    file_path = os.path.join(folder_path, "now_epoch")
    with open(file_path, "w") as file:
        file.write(str(epoch))


def read_epoch(folder_path):
    file_path = os.path.join(folder_path, "now_epoch")
    try:
        with open(file_path, "r") as file:
            oldepoch = file.readline()
    except:
        oldepoch = 0
    return int(oldepoch)
