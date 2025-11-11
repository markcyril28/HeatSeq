#!/bin/bash

# Conda Environment Setup for GEA Pipeline
set -euo pipefail

ENV_NAME="GEA_ENV"
UPDATE_MODE=false

# Parse arguments
if [[ $# -gt 0 ]] && [[ "$1" == "--update" ]]; then
    UPDATE_MODE=true
fi

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
    if [[ "$UPDATE_MODE" == true ]]; then
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
	aria2 parallel-fastq-dump sra-tools \
	hisat2 stringtie samtools bowtie2 rsem salmon trinity trim-galore \
	fastqc multiqc \
	parallel -y || {
    echo "Mamba installation failed, falling back to conda..."
    conda install -n ${ENV_NAME} -c conda-forge -c bioconda \
        aria2 parallel-fastq-dump sra-tools \
        hisat2 stringtie samtools bowtie2 rsem salmon trinity trim-galore \
        fastqc multiqc \
        parallel -y
}

echo "Installation complete. Activate with: conda activate ${ENV_NAME}"
