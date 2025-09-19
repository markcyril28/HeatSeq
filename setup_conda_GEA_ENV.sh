#!/bin/bash

# ===============================================
# Conda Environment Setup for GEA Pipeline
# ===============================================
# This script creates and sets up the GEA_ENV conda environment
# with all required bioinformatics tools

set -euo pipefail

# Configuration
ENV_NAME="GEA_ENV"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YAML_FILE="$SCRIPT_DIR/GEA_ENV_environment.yml"

# Logging functions
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { printf '[%s] %s\n' "$(timestamp)" "$*"; }
log_info() { log "INFO: $*"; }
log_warn() { log "WARN: $*"; }
log_error() { log "ERROR: $*"; }
log_step() { log "===== $* ====="; }

# ===============================================
# Functions
# ===============================================

check_conda() {
    log_step "Checking conda installation"
    if ! command -v conda >/dev/null 2>&1; then
        log_error "Conda is not installed or not in PATH"
        log_error "Please install Miniconda or Anaconda first"
        exit 1
    fi
    
    conda_version=$(conda --version)
    log_info "Found: $conda_version"
    
    # Check for mamba and install if not available
    if ! command -v mamba >/dev/null 2>&1; then
        log_info "Installing mamba for faster dependency resolution..."
        conda install -c conda-forge mamba -y
    fi
    
    # Check if YAML file exists
    if [[ ! -f "$YAML_FILE" ]]; then
        log_error "YAML environment file not found: $YAML_FILE"
        exit 1
    fi
    log_info "Found environment file: $YAML_FILE"
}

setup_environment() {
    log_step "Setting up GEA_ENV conda environment"
    
    # Check if environment exists
    if conda env list | grep -q "^$ENV_NAME "; then
        log_info "Environment '$ENV_NAME' already exists"
        read -p "Do you want to update it? (y/N): " update_env
        if [[ "${update_env,,}" == "y" ]]; then
            log_info "Updating environment '$ENV_NAME' from $YAML_FILE..."
            mamba env update --name "$ENV_NAME" --file "$YAML_FILE"
        fi
    else
        log_info "Creating environment '$ENV_NAME' from $YAML_FILE..."
        mamba env create --file "$YAML_FILE"
    fi
    
    log_info "Environment '$ENV_NAME' setup completed"
}

verify_tools() {
    log_step "Verifying required tools in $ENV_NAME"
    
    # Activate environment for verification
    eval "$(conda shell.bash hook)"
    conda activate "$ENV_NAME"
    
    local required_cmds=(fasterq-dump trim_galore hisat2 samtools stringtie trinity)
    local missing_tools=()
    
    for cmd in "${required_cmds[@]}"; do
        if command -v "$cmd" >/dev/null 2>&1; then
            version_info=$("$cmd" --version 2>&1 | head -1 || echo "version unknown")
            log_info "✓ $cmd: $version_info"
        else
            log_error "✗ $cmd: NOT FOUND"
            missing_tools+=("$cmd")
        fi
    done
    
    # Try to install missing tools
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        log_warn "Attempting to install missing tools: ${missing_tools[*]}"
        
        for tool in "${missing_tools[@]}"; do
            case "$tool" in
                trinity)
                    log_info "Installing trinity..."
                    mamba install -n "$ENV_NAME" -c bioconda trinity -y || {
                        log_error "Failed to install trinity via conda"
                        log_info "Trinity may need to be installed manually"
                        continue
                    }
                    ;;
                *)
                    log_warn "Don't know how to install $tool automatically"
                    ;;
            esac
        done
        
        # Re-verify after installation attempts
        log_info "Re-verifying tools after installation attempts..."
        local still_missing=()
        for cmd in "${missing_tools[@]}"; do
            if command -v "$cmd" >/dev/null 2>&1; then
                version_info=$("$cmd" --version 2>&1 | head -1 || echo "version unknown")
                log_info "✓ $cmd: $version_info (installed successfully)"
            else
                log_error "✗ $cmd: STILL NOT FOUND"
                still_missing+=("$cmd")
            fi
        done
        
        if [[ ${#still_missing[@]} -eq 0 ]]; then
            log_info "All required tools verified successfully after installation!"
            return 0
        else
            log_error "Still missing tools: ${still_missing[*]}"
            log_info "Manual installation may be required for these tools"
            return 1
        fi
    else
        log_info "All required tools verified successfully!"
        return 0
    fi
}

show_usage_instructions() {
    log_step "Usage Instructions"
    
    cat << EOF

GEA_ENV environment has been setup successfully!

To use the environment:
    1. Activate the environment: conda activate $ENV_NAME
    2. Run your GEA pipeline scripts

Environment location: $(conda info --envs | grep "$ENV_NAME" | awk '{print $2}')

Installed tools:
- sra-tools (fasterq-dump)
- trim-galore
- hisat2
- samtools
- stringtie
- trinity

EOF
}

# ===============================================
# Main Execution
# ===============================================

main() {
    log_step "Starting GEA_ENV environment setup"
    
    # Check prerequisites
    check_conda
    
    # Setup environment
    setup_environment
    
    # Verify tools
    if verify_tools; then
        show_usage_instructions
        log_step "Setup completed successfully!"
        exit 0
    else
        log_warn "Setup completed with some missing tools."
        log_info "You can continue using the environment, but some functionality may be limited."
        show_usage_instructions
        exit 0
    fi
} 

# Legacy function for backward compatibility
require_tools() {
    log_warn "require_tools() function is deprecated. Use main setup script instead."
    main "$@"
}

# Run main function
main "$@"
