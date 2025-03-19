import pynvml


def get_available_gpus(min_free_memory=1024):  # 1GB = 1024MB
    """
    Retrieve all available GPUs and check if there are any GPUs with at least min_free_cemory MB of available video memory.

    : param min_free_memory:  Minimum free video memory required for task execution (in MB)
    : return:  List of available GPUs (sorted by remaining memory)
    """
    pynvml.nvmlInit()
    device_count = pynvml.nvmlDeviceGetCount()

    available_gpus = []

    for i in range(device_count):
        handle = pynvml.nvmlDeviceGetHandleByIndex(i)
        mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)

        free_mem_mb = mem_info.free / 1024**2  # MB
        if free_mem_mb >= min_free_memory:
            available_gpus.append((i, free_mem_mb))

    pynvml.nvmlShutdown()

    # sort
    available_gpus.sort(key=lambda x: x[1], reverse=True)

    return [gpu[0] for gpu in available_gpus]


def assign_gpu(min_free_memory=1024):
    """
    Select the most suitable GPU and return its ID
    : param min_free_memory:  Minimum free video memory required for task execution (in MB)
    : return:  Selected GPU ID or None (if no available GPU)
    """
    available_gpus = get_available_gpus(min_free_memory)
    return (
        available_gpus[0] if available_gpus else None
    )  # Choose the GPU with the most remaining video mem


if __name__ == "__main__":
    selected_gpu = assign_gpu(2048)
    if selected_gpu is not None:
        print(f"to GPU: {selected_gpu}")
    else:
        print("no GPU")
