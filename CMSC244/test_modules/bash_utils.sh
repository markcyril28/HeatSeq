#!/bin/bash
# ==============================================================================
# BASH UTILITY FUNCTIONS
# ==============================================================================
# Helper functions for the main run_algorithm.sh script.
# ==============================================================================

print_header() {
    echo ""
    echo "========================================================================"
    echo "$1"
    echo "========================================================================"
}

check_inputs() {
    # Check if required input files exist
    # Args: $1=FASTQ_R1, $2=FASTQ_R2, $3=REFERENCE
    local fastq_r1="$1"
    local fastq_r2="$2"
    local reference="$3"
    
    print_header "CHECKING INPUT FILES"
    
    local missing=0
    
    if [[ ! -f "$fastq_r1" ]]; then
        echo "ERROR: Missing FASTQ R1: $fastq_r1"
        missing=1
    else
        echo "Found: $fastq_r1"
    fi
    
    if [[ ! -f "$fastq_r2" ]]; then
        echo "WARNING: Missing FASTQ R2: $fastq_r2 (single-end mode)"
    else
        echo "Found: $fastq_r2"
    fi
    
    if [[ ! -f "$reference" ]]; then
        echo "ERROR: Missing reference: $reference"
        missing=1
    else
        echo "Found: $reference"
    fi
    
    if [[ $missing -eq 1 ]]; then
        echo "Please ensure all required input files are present."
        return 1
    fi
    return 0
}

create_directories() {
    # Create output directories
    # Args: $@=list of directories to create
    print_header "CREATING OUTPUT DIRECTORIES"
    
    for dir in "$@"; do
        mkdir -p "$dir"
        echo "Created: $dir"
    done
}

activate_conda_env() {
    # Activate conda environment, run setup if needed
    # Args: $1=env_name, $2=setup_script
    local env_name="$1"
    local setup_script="$2"
    
    print_header "ACTIVATING CONDA ENVIRONMENT"
    
    eval "$(conda shell.bash hook)"
    
    if conda activate "$env_name" 2>/dev/null; then
        echo "Activated conda environment: $env_name"
    else
        echo "WARNING: Could not activate $env_name"
        echo "Running setup script..."
        bash "$setup_script"
        conda activate "$env_name"
    fi
    
    python --version
}
