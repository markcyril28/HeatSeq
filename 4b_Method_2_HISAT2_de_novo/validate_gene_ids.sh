#!/bin/bash

# ===============================================
# TSV MERGE VALIDATION SCRIPT
# ===============================================
# This script validates the results of TSV merging operations
# and provides detailed reports on gene matching and data integrity.
# 
# Usage: bash validate_gene_ids.sh <target_group_path> <version>
# 
# Parameters:
#   target_group_path : Path to the target gene group directory
#   version          : Version tag (v1 or v2) to validate
# ===============================================

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <target_group_path> <version>"
    echo ""
    echo "Parameters:"
    echo "  target_group_path : Path to the target gene group directory"
    echo "  version          : Version tag (v1 or v2) to validate"
    echo ""
    echo "Example:"
    echo "  $0 5_stringtie/a_Method2_RESULTS/TEST v1"
    exit 1
}

# Check arguments
if [[ $# -ne 2 ]]; then
    usage
fi

TARGET_GROUP_PATH="$1"
VERSION="$2"

# Validate version
if [[ "$VERSION" != "v1" && "$VERSION" != "v2" ]]; then
    echo "Error: Version must be 'v1' or 'v2'"
    exit 1
fi

# Validate directory
if [[ ! -d "$TARGET_GROUP_PATH" ]]; then
    echo "Error: Target group path does not exist: $TARGET_GROUP_PATH"
    exit 1
fi

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log() { echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $*"; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $*"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $*"; }
error() { echo -e "${RED}[ERROR]${NC} $*"; }

# Validate gene ID consistency
validate_gene_consistency() {
    local files=("$@")
    local tmpdir
    tmpdir=$(mktemp -d)
    trap "rm -rf '$tmpdir'" EXIT
    
    log "Validating gene ID consistency across ${#files[@]} files..."
    
    # Extract gene IDs from each file
    local file_count=0
    for file in "${files[@]}"; do
        ((file_count++))
        local srr_id
        srr_id=$(basename "$file" | sed -n 's/^\(SRR[0-9]\+\)_.*/\1/p')
        tail -n +2 "$file" | cut -f1 | sort > "$tmpdir/genes_${file_count}_${srr_id}.txt"
        echo "  File $file_count: $srr_id ($(wc -l < "$tmpdir/genes_${file_count}_${srr_id}.txt") genes)"
    done
    
    # Compare gene sets
    local first_file="$tmpdir/genes_1_*.txt"
    first_file=$(ls $first_file | head -1)
    local first_count
    first_count=$(wc -l < "$first_file")
    
    local consistent=true
    for ((i=2; i<=file_count; i++)); do
        local current_file
        current_file=$(ls "$tmpdir/genes_${i}_"*.txt | head -1)
        local current_count
        current_count=$(wc -l < "$current_file")
        
        if [[ $current_count -ne $first_count ]]; then
            warning "  Gene count mismatch: File $i has $current_count genes vs $first_count in first file"
            consistent=false
        fi
        
        local common_count
        common_count=$(comm -12 "$first_file" "$current_file" | wc -l)
        local match_percentage
        match_percentage=$(( (common_count * 100) / first_count ))
        
        if [[ $match_percentage -lt 95 ]]; then
            warning "  Gene ID match: File $i has only ${match_percentage}% common genes"
            consistent=false
        else
            log "  Gene ID match: File $i has ${match_percentage}% common genes"
        fi
    done
    
    if [[ "$consistent" == true ]]; then
        success "Gene ID consistency validation passed"
    else
        warning "Gene ID consistency issues detected"
    fi
    
    return 0
}

# Check for expected gene ID patterns
validate_gene_patterns() {
    local files=("$@")
    
    log "Validating gene ID patterns for version $VERSION..."
    
    local expected_pattern
    case "$VERSION" in
        "v1") expected_pattern="STRG\." ;;
        "v2") expected_pattern="MSTRG\." ;;
    esac
    
    for file in "${files[@]}"; do
        local srr_id
        srr_id=$(basename "$file" | sed -n 's/^\(SRR[0-9]\+\)_.*/\1/p')
        
        local total_genes
        total_genes=$(tail -n +2 "$file" | wc -l)
        
        local pattern_matches
        pattern_matches=$(tail -n +2 "$file" | cut -f1 | grep -c "^$expected_pattern" || true)
        
        local match_percentage
        if [[ $total_genes -gt 0 ]]; then
            match_percentage=$(( (pattern_matches * 100) / total_genes ))
        else
            match_percentage=0
        fi
        
        if [[ $match_percentage -ge 90 ]]; then
            success "  $srr_id: ${match_percentage}% genes match expected pattern ($expected_pattern)"
        else
            warning "  $srr_id: Only ${match_percentage}% genes match expected pattern ($expected_pattern)"
        fi
    done
}

# Generate summary report
generate_summary_report() {
    local files=("$@")
    local group_name
    group_name=$(basename "$TARGET_GROUP_PATH")
    
    echo ""
    log "=== VALIDATION SUMMARY REPORT ==="
    log "Group: $group_name"
    log "Version: $VERSION"
    log "Files validated: ${#files[@]}"
    echo ""
    
    log "File Details:"
    for file in "${files[@]}"; do
        local srr_id
        srr_id=$(basename "$file" | sed -n 's/^\(SRR[0-9]\+\)_.*/\1/p')
        local gene_count
        gene_count=$(tail -n +2 "$file" | wc -l)
        local file_size
        file_size=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null || echo "unknown")
        
        printf "  %-12s: %6d genes, %10s bytes\n" "$srr_id" "$gene_count" "$file_size"
    done
    
    echo ""
    log "Recommendations:"
    if [[ ${#files[@]} -gt 1 ]]; then
        log "1. Files are ready for matrix generation"
        log "2. Consider running visualization scripts"
        log "3. Backup validated files before further processing"
    else
        warning "Only 1 file found - consider adding more samples"
    fi
}

# Main function
main() {
    log "Starting validation for version $VERSION in group $(basename "$TARGET_GROUP_PATH")"
    
    # Find all TSV files for the specified version
    local files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(find "$TARGET_GROUP_PATH" -type f -name "*gene_abundances_de_novo_${VERSION}.tsv" -print0)
    
    if [[ ${#files[@]} -eq 0 ]]; then
        warning "No files found for version $VERSION in $TARGET_GROUP_PATH"
        exit 0
    fi
    
    log "Found ${#files[@]} files for validation"
    
    # Run validations
    validate_gene_patterns "${files[@]}"
    
    if [[ ${#files[@]} -gt 1 ]]; then
        validate_gene_consistency "${files[@]}"
    else
        log "Skipping consistency check (only 1 file found)"
    fi
    
    # Generate summary
    generate_summary_report "${files[@]}"
    
    success "Validation completed"
}

# Execute main function
main "$@"
