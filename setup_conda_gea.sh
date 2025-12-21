#!/bin/bash

# Conda Environment Setup for GEA Pipeline
set -euo pipefail

#ENV_NAME="GEA_ENV"
ENV_NAME="gea"
UPDATE_MODE=true
ENV_RESTART_MODE=false

# Parse arguments
for arg in "$@"; do
    case "$arg" in
        --update) UPDATE_MODE=true ;;
        --restart) ENV_RESTART_MODE=true ;;
    esac
done

# Check if mamba is available
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "Using mamba for faster installation"
else
    CONDA_CMD="conda"
    echo "Mamba not found, using conda"
fi

# Create or activate environment
if conda env list | grep -q "^${ENV_NAME}\s"; then
    echo "Environment '${ENV_NAME}' exists"
    if [[ "$ENV_RESTART_MODE" == true ]]; then
        echo "Removing existing environment '${ENV_NAME}'..."
        conda env remove -n ${ENV_NAME} -y
        echo "Recreating environment '${ENV_NAME}'..."
        conda create -n ${ENV_NAME} -y
    elif [[ "$UPDATE_MODE" == true ]]; then
        echo "Updating packages in '${ENV_NAME}'..."
        ${CONDA_CMD} update -n ${ENV_NAME} -c conda-forge -c bioconda --all -y
        echo "Update complete"
        exit 0
    fi
else
    echo "Creating environment '${ENV_NAME}'..."
    conda create -n ${ENV_NAME} -y
fi

echo "Installing packages into '${ENV_NAME}'..."
${CONDA_CMD} install -n ${ENV_NAME} -c conda-forge -c bioconda \
	aria2 parallel-fastq-dump "sra-tools>=3.0" entrez-direct \
	hisat2 stringtie samtools bowtie2 rsem salmon star trim-galore "cutadapt>=4.1" \
	fastqc multiqc \
	parallel \
	r-wgcna r-dynamictreecut r-fastcluster -y || {
    echo "Mamba installation failed, falling back to conda..."
    conda install -n ${ENV_NAME} -c conda-forge -c bioconda \
        aria2 parallel-fastq-dump "sra-tools>=3.0" entrez-direct \
        hisat2 stringtie samtools bowtie2 rsem salmon star trim-galore "cutadapt>=4.1" \
        fastqc multiqc \
        parallel \
        r-wgcna r-dynamictreecut r-fastcluster -y
}

# Configure SRA tools
echo "Configuring SRA tools..."
conda run -n ${ENV_NAME} vdb-config --prefetch-to-cwd

echo "Installation complete. Activate with: conda activate ${ENV_NAME}"
