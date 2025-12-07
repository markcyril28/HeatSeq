#!/bin/bash
# ==============================================================================
# CMSC244 Conda Environment Setup
# ==============================================================================
# Creates conda environment for running alignment algorithm tests
# Default: mamba, falls back to conda if mamba fails
# ==============================================================================

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

ENV_NAME="cmsc"
PYTHON_VERSION="3.10"

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

print_header() {
    echo ""
    echo "========================================================================"
    echo "$1"
    echo "========================================================================"
}

# Determine package manager (mamba preferred, fallback to conda)
get_pkg_manager() {
    if command -v mamba &> /dev/null; then
        echo "mamba"
    else
        echo "conda"
    fi
}

# ==============================================================================
# MAIN SETUP
# ==============================================================================

main() {
    print_header "CMSC244 CONDA ENVIRONMENT SETUP"
    
    # Initialize conda
    eval "$(conda shell.bash hook)"
    
    # Get package manager
    PKG_MGR=$(get_pkg_manager)
    echo "Using package manager: $PKG_MGR"
    
    # Check if environment exists
    if conda env list | grep -q "^${ENV_NAME} "; then
        echo "Environment '$ENV_NAME' already exists."
        read -p "Do you want to recreate it? (y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            echo "Removing existing environment..."
            conda env remove -n "$ENV_NAME" -y
        else
            echo "Keeping existing environment."
            echo "To activate: conda activate $ENV_NAME"
            exit 0
        fi
    fi
    
    # Create environment
    print_header "CREATING CONDA ENVIRONMENT"
    
    if [[ "$PKG_MGR" == "mamba" ]]; then
        echo "Creating environment with mamba..."
        if ! mamba create -n "$ENV_NAME" python="$PYTHON_VERSION" -y; then
            echo "Mamba failed, falling back to conda..."
            PKG_MGR="conda"
            conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y
        fi
    else
        echo "Creating environment with conda..."
        conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y
    fi
    
    # Activate environment
    conda activate "$ENV_NAME"
    
    # Install packages
    print_header "INSTALLING PACKAGES"
    
    # Core scientific packages
    echo "Installing numpy..."
    if ! $PKG_MGR install -n "$ENV_NAME" numpy -y 2>/dev/null; then
        pip install numpy
    fi
    
    # Matplotlib for graphs
    echo "Installing matplotlib..."
    if ! $PKG_MGR install -n "$ENV_NAME" matplotlib -y 2>/dev/null; then
        pip install matplotlib
    fi
    
    # Optional: pandas for data analysis
    echo "Installing pandas..."
    if ! $PKG_MGR install -n "$ENV_NAME" pandas -y 2>/dev/null; then
        pip install pandas
    fi
    
    # Rust toolchain for Rust implementations
    print_header "INSTALLING RUST TOOLCHAIN"
    echo "Installing Rust/Cargo from conda-forge..."
    if ! $PKG_MGR install -n "$ENV_NAME" rust -c conda-forge -y 2>/dev/null; then
        echo "Conda Rust install failed, trying rustup..."
        if ! command -v cargo &> /dev/null; then
            curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
            source "$HOME/.cargo/env"
        fi
    fi
    
    # Verify installation
    print_header "VERIFYING INSTALLATION"
    
    echo "Python version:"
    python --version
    
    echo ""
    echo "Installed packages:"
    pip list | grep -E "numpy|matplotlib|pandas" || true
    
    # Test imports
    echo ""
    echo "Testing imports..."
    python -c "import numpy; print('numpy:', numpy.__version__)"
    python -c "import matplotlib; print('matplotlib:', matplotlib.__version__)" || echo "matplotlib not available"
    
    # Check Rust
    echo ""
    echo "Rust toolchain:"
    if command -v cargo &> /dev/null; then
        cargo --version
        rustc --version
        
        # Build Rust aligners
        print_header "BUILDING RUST ALIGNERS"
        SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        if [[ -d "$SCRIPT_DIR/test_modules/rust_aligners" ]]; then
            cd "$SCRIPT_DIR/test_modules/rust_aligners"
            echo "Building release binaries..."
            cargo build --release
            echo "Rust aligners built successfully."
            ls -la target/release/ | grep -E "hisat2|bowtie2|salmon" || true
            cd "$SCRIPT_DIR"
        else
            echo "Rust aligners directory not found."
        fi
    else
        echo "Rust not available (optional for Rust tests)"
    fi
    
    print_header "SETUP COMPLETE"
    echo ""
    echo "Environment: $ENV_NAME"
    echo "Python: $PYTHON_VERSION"
    echo ""
    echo "To activate: conda activate $ENV_NAME"
    echo "To run tests: bash run_test.sh"
}

# ==============================================================================
# RUN
# ==============================================================================

main "$@"
