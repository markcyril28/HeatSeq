#!/bin/bash

# Conda Environment Setup for GEA Pipeline
set -euo pipefail

ENV_NAME="GEA_HEATMAP_ENV"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YAML_FILE="$SCRIPT_DIR/GEA_ENV_environment.yml"

# Logging
log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
log_step() { log "===== $* ====="; }

check_conda() {
    log_step "Checking conda installation"
    command -v conda >/dev/null 2>&1 || { log "ERROR: Conda not found. Install Miniconda/Anaconda first"; exit 1; }
    log "Found: $(conda --version)"
    
    # Install mamba if missing
    command -v mamba >/dev/null 2>&1 || { log "Installing mamba..."; conda install -c conda-forge mamba -y; }
    
    # Check YAML file
    [[ -f "$YAML_FILE" ]] || { log "ERROR: YAML file not found: $YAML_FILE"; exit 1; }
    log "Environment file: $YAML_FILE"
}

setup_environment() {
    log_step "Setting up GEA_ENV conda environment"
    
    if conda env list | grep -q "^$ENV_NAME "; then
        log "Environment '$ENV_NAME' exists"
        read -p "Update it? (y/N): " update_env
        [[ "${update_env,,}" == "y" ]] && mamba env update --name "$ENV_NAME" --file "$YAML_FILE"
    else
        log "Creating environment '$ENV_NAME'..."
        mamba env create --file "$YAML_FILE"
    fi
}

check_tool() {
    local cmd=$1
    if [[ "$cmd" == "trinity" ]]; then
        # Check both trinity and Trinity (case sensitive)
        if command -v trinity >/dev/null 2>&1; then
            log "✓ trinity: $(trinity --version 2>&1 | head -1 || echo "version unknown")"
            return 0
        elif command -v Trinity >/dev/null 2>&1; then
            log "✓ Trinity: $(Trinity --version 2>&1 | head -1 || echo "version unknown")"
            return 0
        fi
    else
        if command -v "$cmd" >/dev/null 2>&1; then
            log "✓ $cmd: $("$cmd" --version 2>&1 | head -1 || echo "version unknown")"
            return 0
        fi
    fi
    log "✗ $cmd: NOT FOUND"
    return 1
}

verify_tools() {
    log_step "Verifying tools in $ENV_NAME"
    eval "$(conda shell.bash hook)" && conda activate "$ENV_NAME"
    
    local tools=(fasterq-dump trim_galore hisat2 samtools stringtie trinity)
    local missing=()
    
    for tool in "${tools[@]}"; do
        check_tool "$tool" || missing+=("$tool")
    done
    
    # Install missing tools
    if [[ ${#missing[@]} -gt 0 ]]; then
        log "Installing missing tools: ${missing[*]}"
        for tool in "${missing[@]}"; do
            [[ "$tool" == "trinity" ]] && mamba install -n "$ENV_NAME" -c bioconda trinity -y
        done
        
        # Re-verify
        local still_missing=()
        for tool in "${missing[@]}"; do
            check_tool "$tool" || still_missing+=("$tool")
        done
        
        [[ ${#still_missing[@]} -eq 0 ]] && log "All tools installed!" || log "Still missing: ${still_missing[*]}"
    else
        log "All tools verified!"
    fi
}

show_usage() {
    log_step "Setup Complete"
    cat << EOF
To use: conda activate $ENV_NAME
Tools: fasterq-dump, trim-galore, hisat2, samtools, stringtie, trinity
EOF
}

main() {
    log_step "GEA_ENV Setup"
    check_conda
    setup_environment
    verify_tools
    show_usage
}

main "$@"
