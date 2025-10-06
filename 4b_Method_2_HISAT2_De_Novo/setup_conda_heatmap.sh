#!/bin/bash

set -euo pipefail

ENV_NAME="${CONDA_DEFAULT_ENV:-GEA_env}"
YAML_FILE="$(dirname "${BASH_SOURCE[0]}")/Heatmap_ENV_R.yml"

# Check conda
if ! command -v conda >/dev/null; then
    echo "Error: Conda not found" >&2
    exit 1
fi

# Install mamba if needed
command -v mamba >/dev/null || conda install -c conda-forge mamba -y

# Create environment if it doesn't exist
conda env list | grep -q "^$ENV_NAME " || mamba create -n "$ENV_NAME" python=3.11 -y

# Install packages
if [[ -f "$YAML_FILE" ]]; then
    mamba env update --name "$ENV_NAME" --file "$YAML_FILE" --prune
else
    mamba install -n "$ENV_NAME" -c conda-forge -c bioconda -y \
        r-base=4.2 r-tidyverse r-pheatmap r-rcolorbrewer r-viridis r-dplyr r-tibble r-biocmanager
fi

# Install DESeq2
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"
Rscript -e "
if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')
BiocManager::install('DESeq2', ask=FALSE, update=FALSE)
"

echo "Setup complete. Use: conda activate $ENV_NAME"
