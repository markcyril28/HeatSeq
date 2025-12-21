#!/bin/bash
# Delete all conda environments except base
set -euo pipefail

conda env list | grep -v "^#" | grep -v "^base\s" | awk '{print $1}' | grep -v "^$" | while read env; do
    echo "Removing environment: $env"
    conda env remove -n "$env" -y
done

echo "All non-base environments removed."
