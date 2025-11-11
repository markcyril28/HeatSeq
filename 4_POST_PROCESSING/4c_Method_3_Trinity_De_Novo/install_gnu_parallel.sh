#!/bin/bash
# Install GNU Parallel for Trinity pipeline speedup

set -euo pipefail

echo "=== GNU Parallel Installation Helper ==="
echo ""

# Check if parallel is already installed
if command -v parallel >/dev/null 2>&1; then
    version=$(parallel --version 2>&1 | head -n1)
    echo "✓ GNU Parallel is already installed: $version"
    exit 0
fi

echo "GNU Parallel not found. Installing..."
echo ""

# Detect OS
if [[ -f /etc/os-release ]]; then
    . /etc/os-release
    OS=$ID
elif [[ "$OSTYPE" == "darwin"* ]]; then
    OS="macos"
else
    OS="unknown"
fi

case "$OS" in
    ubuntu|debian)
        echo "Detected Ubuntu/Debian system"
        echo "Running: sudo apt-get update && sudo apt-get install -y parallel"
        sudo apt-get update && sudo apt-get install -y parallel
        ;;
    centos|rhel|fedora)
        echo "Detected RHEL/CentOS/Fedora system"
        echo "Running: sudo yum install -y parallel"
        sudo yum install -y parallel
        ;;
    macos)
        echo "Detected macOS system"
        if command -v brew >/dev/null 2>&1; then
            echo "Running: brew install parallel"
            brew install parallel
        else
            echo "Error: Homebrew not found. Install from https://brew.sh"
            exit 1
        fi
        ;;
    *)
        echo "Unknown OS. Please install GNU Parallel manually:"
        echo "  Ubuntu/Debian: sudo apt-get install parallel"
        echo "  RHEL/CentOS:   sudo yum install parallel"
        echo "  macOS:         brew install parallel"
        exit 1
        ;;
esac

# Verify installation
if command -v parallel >/dev/null 2>&1; then
    version=$(parallel --version 2>&1 | head -n1)
    echo ""
    echo "✓ GNU Parallel installed successfully: $version"
    echo ""
    echo "Performance improvement: 2-4x speedup for Salmon quantification"
else
    echo "✗ Installation failed"
    exit 1
fi
