#!/bin/bash
# ==============================================================================
# GPU PREPARATION SCRIPT FOR POST-PROCESSING
# ==============================================================================
# Prepares the GPU environment specifically for post-processing modules
# including R torch GPU acceleration, WGCNA, and matrix computations
#
# R torch uses system CUDA (e.g., /usr/local/cuda) for GPU acceleration
#
# Usage:
#   ./prep_GPU_post-proc.sh           # Standard GPU setup
#   ./prep_GPU_post-proc.sh --check   # Check GPU status only
#   ./prep_GPU_post-proc.sh --install # Full installation with R packages
#
# Supports: NVIDIA GPUs (CUDA), WSL2
# ==============================================================================

set -e

# ==============================================================================
# CONFIGURATION
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POST_PROC_DIR="$SCRIPT_DIR"
LOG_DIR="${LOG_DIR:-$SCRIPT_DIR/logs}"
CONDA_ENV="${CONDA_ENV:-gea}"

# Create log directory with proper error handling
if ! mkdir -p "$LOG_DIR" 2>/dev/null; then
    # Fallback to /tmp if local directory fails
    LOG_DIR="/tmp/heatseq_logs"
    mkdir -p "$LOG_DIR" 2>/dev/null || true
fi

LOG_FILE="$LOG_DIR/gpu_post_proc_$(date +%Y%m%d_%H%M%S).log"

# Ensure log file is writable (create empty file to test)
if ! touch "$LOG_FILE" 2>/dev/null; then
    LOG_FILE="/dev/null"  # Disable logging if we can't write
fi

# ==============================================================================
# LOGGING FUNCTIONS
# ==============================================================================

_log() {
    local level="$1"
    shift
    local msg="$*"
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    echo "[$level] $msg"
    echo "[$timestamp] [$level] $msg" >> "$LOG_FILE" 2>/dev/null || true
}

log_info()  { _log "INFO" "$*"; }
log_warn()  { _log "WARN" "$*"; }
log_error() { _log "ERROR" "$*" >&2; }
log_step()  { echo ""; echo "=== $* ==="; _log "STEP" "$*"; }

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

is_wsl() {
    [[ -f /proc/version ]] && grep -qi "microsoft\|wsl" /proc/version 2>/dev/null
}

check_command() {
    command -v "$1" &>/dev/null
}

check_root() {
    if [[ $EUID -eq 0 ]]; then
        log_error "This script should not be run as root"
        exit 1
    fi
}

# ==============================================================================
# GPU DETECTION
# ==============================================================================

detect_gpu() {
    GPU_AVAILABLE="false"
    GPU_COUNT=0
    GPU_MEMORY_MB=0
    CUDA_VERSION=""
    CUDA_READY="false"

    if check_command nvidia-smi; then
        local gpu_info
        gpu_info=$(nvidia-smi --query-gpu=count,memory.total --format=csv,noheader,nounits 2>/dev/null) || true
        if [[ -n "$gpu_info" ]]; then
            GPU_AVAILABLE="true"
            GPU_COUNT=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | wc -l)
            GPU_MEMORY_MB=$(echo "$gpu_info" | head -1 | cut -d',' -f2 | tr -d ' ')
            CUDA_VERSION=$(nvidia-smi 2>/dev/null | grep -oP 'CUDA Version: \K[0-9.]+' | head -1) || true

            if check_command nvcc; then
                CUDA_READY="true"
            fi
        fi
    fi

    export GPU_AVAILABLE GPU_COUNT GPU_MEMORY_MB CUDA_VERSION CUDA_READY
}

# ==============================================================================
# GPU STATUS CHECK
# ==============================================================================

check_gpu_status() {
    log_step "GPU Status Check"
    
    detect_gpu
    
    if [[ "$GPU_AVAILABLE" == "true" ]]; then
        log_info "GPU Available: YES"
        log_info "GPU Count: $GPU_COUNT"
        log_info "GPU Memory: ${GPU_MEMORY_MB} MB"
        log_info "CUDA Version (Driver): ${CUDA_VERSION:-Not detected}"
        log_info "CUDA Toolkit: $([ "$CUDA_READY" == "true" ] && echo 'Installed' || echo 'Not installed')"
        
        # Check conda CUDA version if available
        local conda_cuda_ver=""
        if [[ -n "${CONDA_PREFIX:-}" ]] && [[ -f "$CONDA_PREFIX/bin/nvcc" ]]; then
            conda_cuda_ver=$("$CONDA_PREFIX/bin/nvcc" --version 2>/dev/null | grep -oP 'release \K[0-9.]+' | head -1) || true
        elif check_command nvcc; then
            conda_cuda_ver=$(nvcc --version 2>/dev/null | grep -oP 'release \K[0-9.]+' | head -1) || true
        fi
        if [[ -n "$conda_cuda_ver" ]]; then
            log_info "CUDA Version (Conda): $conda_cuda_ver"
        else
            log_info "CUDA Version (Conda): Not installed"
        fi
        
        echo ""
        log_info "GPU Details:"
        nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader 2>/dev/null | while read line; do
            log_info "  $line"
        done
    else
        log_warn "No NVIDIA GPU detected"
        if is_wsl; then
            log_info "Running in WSL2 - ensure NVIDIA drivers are installed on Windows host"
            log_info "Download from: https://www.nvidia.com/Download/index.aspx"
        fi
    fi
    
    # Check WSL-specific settings
    if is_wsl; then
        echo ""
        log_info "WSL2 Environment Detected"
        if [[ -e /dev/dxg ]]; then
            log_info "  DirectX GPU interface: Available"
        else
            log_warn "  DirectX GPU interface: Not available"
        fi
    fi
}

# ==============================================================================
# CONDA/MAMBA INITIALIZATION
# ==============================================================================

init_conda() {
    log_step "Initializing Conda Environment"
    
    if [[ -z "${CONDA_EXE:-}" ]]; then
        if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
            source "$HOME/miniconda3/etc/profile.d/conda.sh"
        elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
            source "$HOME/anaconda3/etc/profile.d/conda.sh"
        elif [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
            source "/opt/conda/etc/profile.d/conda.sh"
        else
            log_error "Cannot find conda installation"
            log_info "Please install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
            return 1
        fi
    fi
    
    eval "$(conda shell.bash hook)"
    
    # Check for mamba
    if check_command mamba; then
        PKG_MGR="mamba"
        log_info "Using mamba for faster package management"
    else
        PKG_MGR="conda"
        log_info "Using conda (install mamba for faster operations)"
    fi
    
    export PKG_MGR
}

# ==============================================================================
# CUDA ENVIRONMENT SETUP
# ==============================================================================

setup_cuda_env() {
    log_step "Configuring CUDA Environment"
    
    # Find CUDA installation
    local cuda_home=""
    for path in /usr/local/cuda /usr/local/cuda-* /opt/cuda; do
        if [[ -d "$path" ]]; then
            cuda_home="$path"
            break
        fi
    done
    
    if [[ -n "$cuda_home" ]]; then
        export CUDA_HOME="$cuda_home"
        export PATH="$CUDA_HOME/bin:$PATH"
        export LD_LIBRARY_PATH="$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
        log_info "CUDA_HOME set to: $CUDA_HOME"
    else
        log_info "No system CUDA found (conda CUDA will be used)"
    fi
    
    # WSL2-specific: Set up library paths
    if is_wsl; then
        if [[ -d "/usr/lib/wsl/lib" ]]; then
            export LD_LIBRARY_PATH="/usr/lib/wsl/lib:${LD_LIBRARY_PATH:-}"
            log_info "Added WSL2 GPU library path"
        fi
    fi
}

# ==============================================================================
# GPU PACKAGES INSTALLATION
# ==============================================================================

install_gpu_packages() {
    log_step "Configuring GPU Environment for R torch"
    
    init_conda || return 1
    
    # Check if environment exists
    if ! conda env list | grep -q "^${CONDA_ENV}\s"; then
        log_error "Conda environment '$CONDA_ENV' not found"
        log_info "Run first: cd $POST_PROC_DIR && bash setup_conda_post_proc.sh"
        return 1
    fi
    
    log_info "Activating environment: $CONDA_ENV"
    conda activate "$CONDA_ENV"
    
    log_info "R torch will use system CUDA ${CUDA_VERSION:-unknown}"
    
    # Configure environment for R torch to use system CUDA
    log_info "Configuring environment for R torch CUDA..."
    mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
    cat > "$CONDA_PREFIX/etc/conda/activate.d/cuda_env.sh" << 'EOF'
#!/bin/bash
# CUDA Configuration for R torch
# R torch uses system CUDA via these environment variables

# System CUDA for R torch
if [[ -d "/usr/local/cuda" ]]; then
    export CUDA_HOME="/usr/local/cuda"
    export CUDA_PATH="/usr/local/cuda"
    export PATH="/usr/local/cuda/bin:$PATH"
    export LD_LIBRARY_PATH="/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}"
fi

# WSL2 GPU library support
if [[ -d "/usr/lib/wsl/lib" ]]; then
    export LD_LIBRARY_PATH="/usr/lib/wsl/lib:${LD_LIBRARY_PATH:-}"
fi

# Conda libs
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
EOF
    chmod +x "$CONDA_PREFIX/etc/conda/activate.d/cuda_env.sh"
    
    # Source it now
    source "$CONDA_PREFIX/etc/conda/activate.d/cuda_env.sh"
    
    log_info "GPU environment configured for R torch"
    log_info "R torch will use system CUDA ${CUDA_VERSION:-unknown}"
}

# ==============================================================================
# R GPU PACKAGES
# ==============================================================================

install_r_gpu_packages() {
    log_step "Installing R GPU Packages"
    
    init_conda || return 1
    conda activate "$CONDA_ENV"
    
    # Check if R is available
    if ! check_command Rscript; then
        log_error "R is not available in the conda environment"
        return 1
    fi
    
    # Check system CUDA version (prefer system CUDA for R torch)
    local cuda_ver=""
    
    # First check system CUDA
    if [[ -x "/usr/local/cuda/bin/nvcc" ]]; then
        cuda_ver=$(/usr/local/cuda/bin/nvcc --version 2>/dev/null | grep -oP 'release \K[0-9]+\.[0-9]+' | head -1) || true
        log_info "System CUDA version: ${cuda_ver:-not found}"
    fi
    
    # Fallback to PATH nvcc
    if [[ -z "$cuda_ver" ]] && check_command nvcc; then
        cuda_ver=$(nvcc --version 2>/dev/null | grep -oP 'release \K[0-9]+\.[0-9]+' | head -1) || true
        log_info "CUDA toolkit version: ${cuda_ver:-not found}"
    fi
    
    # R torch supports CUDA - use system CUDA
    local torch_cuda_supported=false
    if [[ -n "$cuda_ver" ]]; then
        torch_cuda_supported=true
        log_info "Using system CUDA $cuda_ver for R torch"
    else
        log_warn "No CUDA toolkit found - R torch will use CPU backend"
    fi
    
    log_info "Installing R torch package..."
    
    # R torch only has pre-built binaries for specific CUDA versions
    # CUDA 12.6 has the latest libtorch (2.7.1) - compatible with CUDA 13.x drivers
    local R_TORCH_CUDA="12.6"
    log_info "R torch will use CUDA $R_TORCH_CUDA binaries (compatible with driver ${cuda_ver:-unknown})"
    
    # Set environment variables BEFORE running R
    # These must be set in the shell environment, not just in R
    export CUDA="$R_TORCH_CUDA"
    export CUDA_VERSION="$R_TORCH_CUDA"
    export TORCH_CUDA_VERSION="$R_TORCH_CUDA"
    
    # Ensure CUDA libraries are in LD_LIBRARY_PATH
    if [[ -d "/usr/local/cuda/lib64" ]]; then
        export LD_LIBRARY_PATH="/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}"
    fi
    if [[ -d "/usr/lib/wsl/lib" ]]; then
        export LD_LIBRARY_PATH="/usr/lib/wsl/lib:${LD_LIBRARY_PATH:-}"
    fi
    
    log_info "Set CUDA_VERSION=$R_TORCH_CUDA for R torch installation"
    log_info "LD_LIBRARY_PATH includes CUDA libs"
    
    # First, remove any existing torch installation to force clean CUDA install
    log_info "Removing existing R torch installation (if any)..."
    
    # Remove torch cache directories from shell (more reliable than R)
    local torch_home="${HOME}/.torch"
    local torch_cache="${HOME}/.cache/torch"
    
    if [[ -d "$torch_home" ]]; then
        log_info "Removing $torch_home..."
        rm -rf "$torch_home"
    fi
    
    if [[ -d "$torch_cache" ]]; then
        log_info "Removing $torch_cache..."
        rm -rf "$torch_cache"
    fi
    
    # Remove torch deps from R library paths
    Rscript --no-save -e '
# Remove torch deps from all R library paths
lib_paths <- .libPaths()
for (lib_path in lib_paths) {
    torch_pkg_path <- file.path(lib_path, "torch")
    if (dir.exists(torch_pkg_path)) {
        deps_path <- file.path(torch_pkg_path, "deps")
        if (dir.exists(deps_path)) {
            cat("Removing torch deps at:", deps_path, "\n")
            unlink(deps_path, recursive = TRUE)
        }
        # Also remove lantern folder
        lantern_path <- file.path(torch_pkg_path, "lantern")
        if (dir.exists(lantern_path)) {
            cat("Removing lantern at:", lantern_path, "\n")
            unlink(lantern_path, recursive = TRUE)
        }
    }
}
cat("R library cleanup complete\n")
' 2>&1 || true
    
    log_info "Cleanup complete"

    # Now install torch with CUDA
    log_info "Installing R torch with CUDA $R_TORCH_CUDA backend..."
    
    # Create an R script file to ensure environment is properly set
    local r_install_script=$(mktemp /tmp/install_torch_cuda.R.XXXXXX)
    cat > "$r_install_script" << 'RSCRIPT'
# Force CUDA version for libtorch download
# CUDA 12.6 has the latest libtorch 2.7.1 available
cuda_ver <- "12.6"

cat("=== R torch CUDA Installation ===\n")
cat("Target CUDA version:", cuda_ver, "\n\n")

# Set ALL environment variables that torch might check
Sys.setenv(CUDA = cuda_ver)
Sys.setenv(CUDA_VERSION = cuda_ver)
Sys.setenv(TORCH_CUDA_VERSION = cuda_ver)
Sys.setenv(TORCH_INSTALL = "cuda")  # Critical: tells torch to install CUDA version

cat("Environment check:\n")
cat("  CUDA:", Sys.getenv("CUDA"), "\n")
cat("  TORCH_INSTALL:", Sys.getenv("TORCH_INSTALL"), "\n")

# Increase timeout for large download (~2GB)
options(timeout = 1200)

# Install torch package if needed
if (!requireNamespace("torch", quietly = TRUE)) {
    cat("\nInstalling torch R package from CRAN...\n")
    install.packages("torch", repos = "https://cloud.r-project.org")
}

cat("\nInstalling libtorch with CUDA backend...\n")
cat("This downloads ~2GB, please be patient...\n\n")

# Use install_torch with explicit parameters
# The 'type' parameter is deprecated in newer versions, use 'backend' instead
tryCatch({
    # Try newer API first (backend parameter)
    torch::install_torch(
        backend = "cuda",
        reinstall = TRUE
    )
    cat("\n=== libtorch CUDA installation completed ===\n")
}, error = function(e) {
    cat("First attempt error:", conditionMessage(e), "\n")
    cat("\nTrying with type parameter...\n")
    
    tryCatch({
        torch::install_torch(
            type = "cuda",
            reinstall = TRUE
        )
        cat("\n=== libtorch CUDA installation completed ===\n")
    }, error = function(e2) {
        cat("Second attempt error:", conditionMessage(e2), "\n")
    })
})

cat("\n=== Installation Complete ===\n")
cat("IMPORTANT: CUDA will only work after restarting R!\n")
cat("Run verification in a NEW R session.\n")
RSCRIPT

    # Run the R script with environment variables set
    # TORCH_INSTALL=cuda is the key variable for newer torch versions
    CUDA="$R_TORCH_CUDA" \
    CUDA_VERSION="$R_TORCH_CUDA" \
    TORCH_CUDA_VERSION="$R_TORCH_CUDA" \
    TORCH_INSTALL="cuda" \
    Rscript --no-save "$r_install_script" 2>&1 || log_warn "R torch installation had some issues"
    
    rm -f "$r_install_script"
    
    # Verify in a FRESH R session (critical for CUDA detection)
    log_info "Verifying CUDA in fresh R session..."
    CUDA="$R_TORCH_CUDA" \
    TORCH_INSTALL="cuda" \
    Rscript --no-save -e '
library(torch)
cat("=== Fresh Session Verification ===\n")
cuda_ok <- tryCatch(cuda_is_available(), error = function(e) FALSE)
cat("CUDA available:", cuda_ok, "\n")
if (cuda_ok) {
    cat("CUDA devices:", cuda_device_count(), "\n")
    cat("GPU:", cuda_get_device_name(0), "\n")
    x <- torch_randn(10, 10, device = "cuda")
    cat("GPU tensor test: PASSED\n")
} else {
    # Check what was installed
    torch_home <- file.path(Sys.getenv("HOME"), ".torch")
    cat("\nDiagnostics:\n")
    cat("TORCH_HOME:", torch_home, "\n")
    if (dir.exists(torch_home)) {
        cuda_files <- list.files(torch_home, pattern = "cuda|nvrtc|cublas", recursive = TRUE)
        cat("CUDA libs found:", length(cuda_files), "\n")
        if (length(cuda_files) == 0) {
            cat("ERROR: CPU-only libtorch was installed!\n")
        }
    }
}
' 2>&1 || true

    # Check for GPU-accelerated WGCNA support
    log_info "Verifying WGCNA installation..."
    Rscript --no-save -e '
if (requireNamespace("WGCNA", quietly = TRUE)) {
    cat("WGCNA: Installed\n")
    # Enable multi-threading for correlation calculations
    WGCNA::enableWGCNAThreads(nThreads = NULL)
    cat("WGCNA threading enabled\n")
} else {
    cat("WGCNA: Not installed\n")
}
' 2>&1 || true

    log_info "R packages setup complete"
    if [[ "$torch_cuda_supported" == "true" ]]; then
        log_info "R torch configured to use system CUDA $cuda_ver"
    else
        log_info "R torch using CPU backend"
    fi
}

# ==============================================================================
# VERIFICATION
# ==============================================================================

verify_gpu_setup() {
    log_step "Verifying GPU Setup for Post-Processing"
    
    local all_ok=true
    
    # Check nvidia-smi
    if check_command nvidia-smi && nvidia-smi &>/dev/null; then
        log_info "✓ nvidia-smi: Working"
    else
        log_warn "✗ nvidia-smi: Not available"
        all_ok=false
    fi
    
    # Check conda environment first and use conda CUDA for testing
    init_conda 2>/dev/null || true
    if conda env list 2>/dev/null | grep -q "^${CONDA_ENV}\s"; then
        log_info "✓ Conda environment '$CONDA_ENV': Exists"
        
        # Activate conda environment to use conda CUDA
        conda activate "$CONDA_ENV" 2>/dev/null || true
        
        # Source conda CUDA environment if available
        if [[ -f "$CONDA_PREFIX/etc/conda/activate.d/cuda_env.sh" ]]; then
            source "$CONDA_PREFIX/etc/conda/activate.d/cuda_env.sh"
            log_info "✓ Conda CUDA environment: Loaded"
        fi
        
        # Check system CUDA (for R torch)
        if [[ -x "/usr/local/cuda/bin/nvcc" ]]; then
            local sys_cuda_ver=$(/usr/local/cuda/bin/nvcc --version 2>/dev/null | grep release | awk '{print $6}' | tr -d ',')
            log_info "✓ System CUDA: $sys_cuda_ver (for R torch)"
        elif check_command nvcc; then
            local sys_cuda_ver=$(nvcc --version 2>/dev/null | grep release | awk '{print $6}' | tr -d ',')
            log_info "✓ CUDA toolkit: $sys_cuda_ver"
        else
            log_info "○ System CUDA: Not found (R torch will use CPU)"
        fi
        
        # Check R torch
        if Rscript -e 'cat(requireNamespace("torch", quietly=TRUE))' 2>/dev/null | grep -q "TRUE"; then
            local cuda_status=$(Rscript -e 'if(requireNamespace("torch",quietly=TRUE)){cat(torch::cuda_is_available())}' 2>/dev/null || echo "FALSE")
            if [[ "$cuda_status" == "TRUE" ]]; then
                log_info "✓ R torch CUDA: Available"
            else
                log_info "○ R torch: Installed (CPU mode)"
            fi
        else
            log_info "○ R torch: Not installed"
        fi
    else
        log_warn "✗ Conda environment '$CONDA_ENV': Not found"
        all_ok=false
        
        # Fallback: Check system CUDA if no conda env
        if check_command nvcc; then
            log_info "✓ System CUDA toolkit: $(nvcc --version 2>/dev/null | grep release | awk '{print $6}' | tr -d ',')"
        else
            log_info "○ CUDA toolkit: Not in PATH"
        fi
    fi
    
    echo ""
    if $all_ok; then
        log_info "GPU setup verification: PASSED"
    else
        log_warn "GPU setup verification: PARTIAL (some components missing)"
        log_info "Run: $0 --install  to install missing components"
    fi
}

# ==============================================================================
# USAGE
# ==============================================================================

show_usage() {
    cat << EOF
GPU Preparation Script for Post-Processing

Usage: $0 [OPTIONS]

Options:
  --check     Check GPU status and availability only
  --install   Full installation: GPU packages + R packages
  --verify    Verify GPU setup for post-processing
  --help      Show this help message

Without options, performs GPU check and environment setup.

Examples:
  $0                  # Check GPU and setup environment
  $0 --check          # Quick GPU status check
  $0 --install        # Install all GPU packages
  $0 --verify         # Verify complete setup

Environment Variables:
  CONDA_ENV           Conda environment name (default: gea)
  LOG_DIR             Log directory (default: ./logs)

Requirements:
  - NVIDIA GPU with drivers installed
  - For WSL2: NVIDIA drivers on Windows host
  - Conda/Miniconda installed
  - 'gea' conda environment (from setup_conda_post_proc.sh)

EOF
}

# ==============================================================================
# MAIN
# ==============================================================================

main() {
    log_step "GPU Preparation for Post-Processing"
    log_info "Log file: $LOG_FILE"
    log_info "Script directory: $SCRIPT_DIR"
    
    check_root
    
    case "${1:-}" in
        --check)
            check_gpu_status
            ;;
        --install)
            check_gpu_status
            setup_cuda_env
            install_gpu_packages
            install_r_gpu_packages
            verify_gpu_setup
            ;;
        --verify)
            check_gpu_status
            verify_gpu_setup
            ;;
        --help|-h)
            show_usage
            exit 0
            ;;
        "")
            # Default: check and setup environment
            check_gpu_status
            setup_cuda_env
            verify_gpu_setup
            
            echo ""
            log_step "Next Steps"
            if [[ "$GPU_AVAILABLE" == "true" ]]; then
                log_info "GPU detected! To install GPU packages:"
                log_info "  $0 --install"
                echo ""
                log_info "Or run post-processing with GPU enabled:"
                log_info "  cd $POST_PROC_DIR"
                log_info "  ENABLE_GPU=TRUE bash run_all_post_processing.sh"
            else
                log_warn "No GPU detected. Post-processing will use CPU."
                log_info "For GPU acceleration, ensure NVIDIA drivers are installed."
            fi
            ;;
        *)
            log_error "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
    
    echo ""
    log_info "Log saved to: $LOG_FILE"
}

# Run main
main "$@"
