#!/bin/bash

#===============================================================================
# HISAT2 De Novo Pipeline - Heatmap Wrapper Script
#===============================================================================
# Switches to Heatmap_ENV conda environment and runs StringTie + heatmap generation
# Usage: bash 0_Heatmap_Wrapper.sh
#===============================================================================

# Configuration
readonly REQUIRED_ENV="Heatmap_ENV"
readonly LOGS_DIR="logs"
readonly LOG_FILE="$LOGS_DIR/heatmap_wrapper_$(date +"%Y%m%d_%H%M%S").log"

# Logging functions
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1" | tee -a "$LOG_FILE"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" | tee -a "$LOG_FILE" >&2
}

# Error handling
set -euo pipefail
trap 'log_error "Script failed at line $LINENO"' ERR

# Main execution
main() {
    mkdir -p "$LOGS_DIR"
    log_message "Starting HISAT2 De Novo Heatmap Wrapper"
    log_message "Current working directory: $(pwd)"
    log_message "Contents of current directory:"
    ls -la | head -10 | while read line; do log_message "  $line"; done
    
    # Initialize conda/mamba for hook setup
    if command -v mamba &> /dev/null; then
        eval "$(mamba shell bash hook)"
        CONDA_CMD="mamba"
    elif command -v conda &> /dev/null; then
        eval "$(conda shell.bash hook)"
        CONDA_CMD="conda"
    else
        log_error "Neither mamba nor conda found in PATH"
        exit 1
    fi
    
    # Setup environment
    if [[ -f "setup_mamba_heatmap.sh" ]]; then
        log_message "Running setup_mamba_heatmap.sh"
        #bash setup_mamba_heatmap.sh
    else
        log_message "setup_mamba_heatmap.sh not found, skipping setup"
    fi
    
    if [[ -f "install_heatmap_libraries.R" ]]; then
        log_message "Installing R libraries"
        #Rscript install_heatmap_libraries.R
    else
        log_message "install_heatmap_libraries.R not found, skipping R setup"
    fi
    
    # Activate environment
    log_message "Activating environment: $REQUIRED_ENV"
    source activate "$REQUIRED_ENV"
    
    if [[ "$CONDA_DEFAULT_ENV" != "$REQUIRED_ENV" ]]; then
        log_error "Failed to activate $REQUIRED_ENV environment"
        exit 1
    fi
    
    # Execute scripts
    log_message "Executing StringTie matrix builder..."
    
    # Ensure we're in the correct directory
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    cd "$SCRIPT_DIR"
    log_message "Changed to script directory: $(pwd)"
    
    if [[ ! -f "a_stringtie_matrix_builder_m2_v4.sh" ]]; then
        log_error "Script not found: a_stringtie_matrix_builder_m2_v4.sh"
        log_error "Available files:"
        ls -la *.sh | while read line; do log_error "  $line"; done
        exit 1
    fi
    
    bash a_stringtie_matrix_builder_m2_v4.sh
    
    log_message "Generating heatmaps..."
    
    if [[ ! -f "b_make_heatmap_of_matrices_v4.R" ]]; then
        log_error "Script not found: b_make_heatmap_of_matrices_v4.R"
        exit 1
    fi
    
    Rscript --verbose b_make_heatmap_of_matrices_v4.R
    
    log_message "Pipeline completed successfully!"
}

main "$@"