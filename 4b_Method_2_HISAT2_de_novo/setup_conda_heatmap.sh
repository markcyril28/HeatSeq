#!/bin/bash

# ===============================================
# Conda R Package Installation for Heatmap Pipeline
# ===============================================
# This script installs required R packages using a YAML environment file
# for the R script: b_heatmap_DeSeq2_v2.R

set -euo pipefail

# Configuration
DEFAULT_ENV="${CONDA_DEFAULT_ENV:-Heatmap_ENV}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YAML_FILE="$SCRIPT_DIR/Heatmap_ENV_R.yml"

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
    
    # Check for mamba and install if not available (install in base env only)
    if ! command -v mamba >/dev/null 2>&1; then
        log_info "Installing mamba in base environment for faster dependency resolution..."
        conda install -c conda-forge mamba -y
    fi
    
    log_info "Target environment: $DEFAULT_ENV"
    
    # Check if environment exists, create if necessary
    if ! conda env list | grep -q "^$DEFAULT_ENV "; then
        log_info "Environment '$DEFAULT_ENV' does not exist. Creating it with Python 3.11..."
        mamba create -n "$DEFAULT_ENV" python=3.11 -y
        log_info "Environment '$DEFAULT_ENV' created successfully"
    else
        log_info "Environment '$DEFAULT_ENV' already exists"
    fi
    
    # Verify the target environment exists before proceeding
    if ! conda env list | grep -q "^$DEFAULT_ENV "; then
        log_error "Failed to create or find environment '$DEFAULT_ENV'"
        exit 1
    fi
    
    # Check if YAML file exists
    if [[ ! -f "$YAML_FILE" ]]; then
        log_error "YAML environment file not found: $YAML_FILE"
        exit 1
    fi
    log_info "Found environment file: $YAML_FILE"
}

install_r_packages() {
    log_step "Installing R packages using YAML environment file"
    
    # Double-check target environment exists
    if ! conda env list | grep -q "^$DEFAULT_ENV "; then
        log_error "Target environment '$DEFAULT_ENV' does not exist!"
        exit 1
    fi
    
    # Activate the target environment before installation
    log_info "Activating environment '$DEFAULT_ENV' for package installation..."
    eval "$(conda shell.bash hook)"
    conda activate "$DEFAULT_ENV"
    
    # Verify we're in the correct environment
    if [[ "${CONDA_DEFAULT_ENV}" != "$DEFAULT_ENV" ]]; then
        log_error "Failed to activate environment '$DEFAULT_ENV'. Currently in: ${CONDA_DEFAULT_ENV}"
        exit 1
    fi
    
    log_info "Successfully activated environment '$DEFAULT_ENV'"
    log_info "Installing packages from $YAML_FILE strictly into environment '$DEFAULT_ENV'..."
    
    # Use explicit --name flag to ensure installation in correct environment
    mamba env update --name "$DEFAULT_ENV" --file "$YAML_FILE" --prune || {
        log_warn "YAML installation failed, trying alternative approach..."
        
        # Install packages individually if YAML fails
        log_info "Installing R packages individually..."
        mamba install -n "$DEFAULT_ENV" -c conda-forge -c bioconda -y \
            r-base=4.2 \
            r-tidyverse \
            r-pheatmap \
            r-rcolorbrewer \
            r-viridis \
            r-dplyr \
            r-tibble \
            r-biocmanager || {
            log_error "Failed to install basic R packages"
            exit 1
        }
    }
    
    # Install DESeq2 via BiocManager with simplified approach
    log_info "Installing DESeq2 and other required R libraries via BiocManager..."
    Rscript -e "
    # Install BiocManager if not available
    if (!requireNamespace('BiocManager', quietly = TRUE)) {
        install.packages('BiocManager', repos='https://cran.r-project.org/')
    }

    # Install CRAN packages
    pkgs <- c('pheatmap', 'RColorBrewer', 'dplyr', 'tibble')
    for (pkg in pkgs) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            install.packages(pkg, repos='https://cran.r-project.org/')
        }
    }

    # Install DESeq2 without specifying Bioconductor version
    cat('Installing DESeq2...\n')
    BiocManager::install('DESeq2', ask = FALSE, update = FALSE)

    # Verify installation
    libs <- c('BiocManager', 'DESeq2', 'pheatmap', 'RColorBrewer', 'dplyr', 'tibble')
    for (lib in libs) {
        if (requireNamespace(lib, quietly = TRUE)) {
            cat('✓', lib, 'successfully installed\n')
        } else {
            cat('✗', lib, 'installation failed\n')
            quit(status = 1)
        }
    }
    " 2>&1 | tee /tmp/deseq2_install.log || {
        log_warn "BiocManager installation failed. Log:"
        cat /tmp/deseq2_install.log
        log_info "Trying conda approach as fallback..."
        mamba install -n "$DEFAULT_ENV" -c bioconda bioconductor-deseq2 -y || {
            log_error "All DESeq2 installation methods failed"
            log_error "You may need to install DESeq2 manually after setup completes"
            log_info "Manual installation command:"
            log_info "conda activate $DEFAULT_ENV"
            log_info "Rscript -e \"BiocManager::install('DESeq2')\""
            # Don't exit here, continue with other packages
        }
    }
    
    # Verify installation was targeted correctly
    log_info "Verifying packages were installed in '$DEFAULT_ENV'..."
    
    # List packages in the target environment to confirm
    if mamba list -n "$DEFAULT_ENV" | grep -q "r-base"; then
        log_info "✓ R packages successfully installed in '$DEFAULT_ENV'"
    else
        log_error "✗ R packages not found in '$DEFAULT_ENV' - installation may have failed"
        exit 1
    fi
    
    log_info "R packages installed successfully in environment '$DEFAULT_ENV'"
}

verify_installation() {
    log_step "Verifying R packages installation in $DEFAULT_ENV"
    
    # Verify we're targeting the correct environment
    if ! conda env list | grep -q "^$DEFAULT_ENV "; then
        log_error "Target environment '$DEFAULT_ENV' does not exist!"
        return 1
    fi
    
    # Activate environment for verification
    eval "$(conda shell.bash hook)"
    conda activate "$DEFAULT_ENV"
    
    # Double-check we're in the correct environment
    if [[ "${CONDA_DEFAULT_ENV}" != "$DEFAULT_ENV" ]]; then
        log_error "Failed to activate environment '$DEFAULT_ENV'. Currently in: ${CONDA_DEFAULT_ENV}"
        return 1
    fi
    
    log_info "Successfully activated environment '$DEFAULT_ENV'"
    
    # Check R packages
    log_info "Checking R packages..."
    if command -v Rscript >/dev/null 2>&1; then
        log_info "✓ R: $(R --version | head -1)"
        
        # Test R packages
        Rscript -e "
        packages <- c('tidyverse', 'DESeq2', 'pheatmap', 'RColorBrewer', 'viridis', 'ggplot2', 'dplyr', 'readr')
        missing_packages <- c()
        for (pkg in packages) {
            if (requireNamespace(pkg, quietly = TRUE)) {
                cat('✓', pkg, 'is available\n')
            } else {
                cat('✗', pkg, 'is missing\n')
                missing_packages <- c(missing_packages, pkg)
            }
        }
        if (length(missing_packages) > 0) {
            cat('Missing packages:', paste(missing_packages, collapse=', '), '\n')
            quit(status = 1)
        }
        " || {
            log_error "Some R packages are missing"
            return 1
        }
    else
        log_error "✗ R/Rscript: NOT FOUND"
        return 1
    fi
    
    log_info "All R packages verified successfully!"
    return 0
}

show_usage_instructions() {
    log_step "Usage Instructions"
    
    cat << EOF

R packages have been installed successfully in environment '$DEFAULT_ENV'!

To run the R heatmap script:
    1. Activate the environment: conda activate $DEFAULT_ENV
    2. Run the script: Rscript b_heatmap_DeSeq2_v2.R

Current environment location: $(conda info --envs | grep "$DEFAULT_ENV" | awk '{print $2}')

Installed R packages:
- tidyverse (includes ggplot2, dplyr, readr)
- DESeq2 (Bioconductor)
- pheatmap
- RColorBrewer
- viridis

EOF
}

# ===============================================
# Main Execution
# ===============================================

main() {
    log_step "Starting R package installation for Heatmap Pipeline in $DEFAULT_ENV"
    
    # Store current environment to restore later if needed
    ORIGINAL_ENV="${CONDA_DEFAULT_ENV:-base}"
    
    # Check prerequisites
    check_conda
    
    # Install R packages using YAML file (this will activate the target env)
    install_r_packages
    
    # Verify everything works (verification also activates the env)
    if verify_installation; then
        show_usage_instructions
        log_step "Setup completed successfully!"
        log_info "Environment '$DEFAULT_ENV' is now active and ready to use"
        exit 0
    else
        log_error "Setup completed with errors. Please check the installation."
        exit 1
    fi
}

# Run main function
main "$@"