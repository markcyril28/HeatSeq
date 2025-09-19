#!/bin/bash

# ===============================================
# Conda Package Installation for GEA Pipeline
# ===============================================
# This script installs required packages into your current conda environment
# for both bash and R scripts:
# - a_stringtie_method_2.2.sh
# - b_heatmap_DeSeq2_v2.R

set -euo pipefail

# Configuration
CURRENT_ENV="${CONDA_DEFAULT_ENV:-Heatmap_ENV}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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
    log_step "Checking conda installation and current environment"
    if ! command -v conda >/dev/null 2>&1; then
        log_error "Conda is not installed or not in PATH"
        log_error "Please install Miniconda or Anaconda first"
        log_error "Download from: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
    
    conda_version=$(conda --version)
    log_info "Found: $conda_version"
    
    if [[ "$CURRENT_ENV" == "base" ]]; then
        log_warn "You are in the base environment"
        read -p "Do you want to continue installing packages in base? (y/N): " continue_base
        if [[ "${continue_base,,}" != "y" ]]; then
            log_info "Please activate your preferred environment and run this script again"
            exit 0
        fi
    fi
    
    log_info "Current environment: $CURRENT_ENV"
}

install_bioinformatics_tools() {
    log_step "Installing bioinformatics tools in $CURRENT_ENV"
    
    log_info "Installing core bioinformatics tools..."
    conda install -c conda-forge -c bioconda -y \
        stringtie=2.2.1
    
    log_info "Bioinformatics tools installed successfully"
}

install_r_packages() {
    log_step "Installing R and required packages in $CURRENT_ENV"
    
    # Install R and basic packages
    log_info "Installing R base and essential packages..."
    conda install -c conda-forge -y \
        r-base=4.3.1 \
        r-essentials \
        r-tidyverse \
        r-devtools \
        r-remotes
    
    # Install Bioconductor packages
    log_info "Installing Bioconductor packages..."
    conda install -c bioconda -c conda-forge -y \
        bioconductor-deseq2 \
        r-pheatmap \
        r-rcolorbrewer \
        r-viridis \
        r-ggplot2 \
        r-dplyr \
        r-readr
    
    log_info "R packages installed successfully"
}

install_additional_tools() {
    log_step "Installing additional utilities in $CURRENT_ENV"
    
    conda install -c conda-forge -y \
        time \
        parallel \
        dos2unix \
        wget \
        curl
    
    log_info "Additional tools installed successfully"
}

verify_installation() {
    log_step "Verifying installation in $CURRENT_ENV"
    
    # Check bioinformatics tools
    local tools=("stringtie")
    local missing_tools=()
    
    for tool in "${tools[@]}"; do
        if command -v "$tool" >/dev/null 2>&1; then
            version_info=$("$tool" --version 2>&1 | head -1 || echo "version unknown")
            log_info "✓ $tool: $version_info"
        else
            log_error "✗ $tool: NOT FOUND"
            missing_tools+=("$tool")
        fi
    done
    
    # Check R packages
    log_info "Checking R packages..."
    if command -v Rscript >/dev/null 2>&1; then
        log_info "✓ R: $(R --version | head -1)"
        
        # Test R packages
        Rscript -e "
        packages <- c('tidyverse', 'DESeq2', 'pheatmap', 'RColorBrewer', 'viridis')
        for (pkg in packages) {
            if (requireNamespace(pkg, quietly = TRUE)) {
                cat('✓', pkg, 'is available\n')
            } else {
                cat('✗', pkg, 'is missing\n')
            }
        }" 2>/dev/null || log_warn "Could not verify R packages"
    else
        log_error "✗ R/Rscript: NOT FOUND"
        missing_tools+=("R")
    fi
    
    if [[ ${#missing_tools[@]} -eq 0 ]]; then
        log_info "All tools verified successfully!"
        return 0
    else
        log_error "Missing tools: ${missing_tools[*]}"
        return 1
    fi
}

show_usage_instructions() {
    log_step "Usage Instructions"
    
    cat << EOF

Packages have been installed successfully in environment '$CURRENT_ENV'!

To run the pipeline scripts (in your current environment):
    1. Run bash script: bash a_stringtie_method_2.2.sh
    2. Run R script: Rscript b_heatmap_DeSeq2_v2.R

If you switch environments, you'll need to run this setup script again
in the new environment.

Current environment location: $(conda info --envs | grep "$CURRENT_ENV" | awk '{print $2}')

EOF
}

# ===============================================
# Main Execution
# ===============================================

main() {
    log_step "Starting package installation for GEA Pipeline in $CURRENT_ENV"
    
    # Check prerequisites
    check_conda
    
    # Install packages in current environment
    install_bioinformatics_tools
    install_r_packages
    install_additional_tools
    
    # Verify everything works
    if verify_installation; then
        show_usage_instructions
        log_step "Setup completed successfully!"
        exit 0
    else
        log_error "Setup completed with errors. Please check the installation."
        exit 1
    fi
}

# Run main function
main "$@"