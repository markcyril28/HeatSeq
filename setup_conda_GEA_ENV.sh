#!/bin/bash

# Conda Environment Setup for GEA Pipeline
set -euo pipefail

echo "Current conda environment: $CONDA_DEFAULT_ENV"
read -p "Make sure you are installing into the correct conda environment. Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Installation cancelled."
    exit 1
fi

mamba install -y -c conda-forge -c bioconda -c defaults python=3.11 \
    sra-tools trim-galore hisat2 samtools stringtie \
    trinity parallel wget curl