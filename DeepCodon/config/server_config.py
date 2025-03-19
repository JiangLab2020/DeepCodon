from DeepCodon.config.config import PROJECT_ROOT

# log
server_log_path = PROJECT_ROOT / "DeepCodon" / "data" / "userData" / "log"


# cache
blastp_output_path = (
    PROJECT_ROOT / "DeepCodon" / "data" / "userData" / "cache" / "blastp_output"
)
model_cache_path = (
    PROJECT_ROOT / "DeepCodon" / "data" / "userData" / "cache" / "model_cache"
)

# output
defeat_output_path = PROJECT_ROOT / "DeepCodon" / "data" / "userData" / "output"

# tools
blastp_dataset_path = (
    PROJECT_ROOT
    / "DeepCodon"
    / "data"
    / "datasets"
    / "rareCodonData"
    / "output_sequences.fasta"
)

# datasets
wp2UniRef_path = (
    PROJECT_ROOT
    / "DeepCodon"
    / "data"
    / "datasets"
    / "rareCodonData"
    / "found_files.csv"
)
rareCodon_path = (
    PROJECT_ROOT / "DeepCodon" / "data" / "datasets" / "rareCodonData" / "4train.csv"
)
rareCodon_db_path = (
    PROJECT_ROOT / "DeepCodon" / "data" / "datasets" / "rareCodonData" / "rareCodon.db"
)
tai_cache_path = (
    PROJECT_ROOT / "DeepCodon" / "data" / "userData" / "cache" / "tai_cache"
)

# variables
coverage_threshold = 0.8  # Sequence coverage threshold
pident_threshold = 0.8  # Sequence similarity threshold
smooth_window_size = 5  # Sliding window size

# logit_bias
use_logit_bias = False
ϵ = 1e-4
λTemplate = 10
