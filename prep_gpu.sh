#!/bin/bash
# ==============================================================================
# GPU PREPARATION SCRIPT
# ==============================================================================
# Prepares system for GPU monitoring and CUDA acceleration
# Supports: NVIDIA, AMD, Intel GPUs
# ==============================================================================

set -e

# ==============================================================================
# CONFIGURATION
# ==============================================================================

NVIDIA_DRIVER_VERSION="${NVIDIA_DRIVER_VERSION:-535}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ==============================================================================
# SOURCE GPU UTILITIES
# ==============================================================================

if [[ -f "$SCRIPT_DIR/modules/logging/gpu_utils.sh" ]]; then
	source "$SCRIPT_DIR/modules/logging/gpu_utils.sh"
else
	echo "[ERROR] gpu_utils.sh not found at $SCRIPT_DIR/modules/logging/gpu_utils.sh"
	exit 1
fi

# ==============================================================================
# MAIN
# ==============================================================================

gpu_prep_main "$@"
