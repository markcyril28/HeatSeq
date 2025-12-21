#!/bin/bash

# ==============================================================================
# CLEANUP UTILITIES
# ==============================================================================
# Functions for cleaning up intermediate files and directories
# ==============================================================================

set -euo pipefail

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -f "$SCRIPT_DIR/logging_utils.sh" ]] && source "$SCRIPT_DIR/logging_utils.sh"
[[ -f "$SCRIPT_DIR/config.sh" ]] && source "$SCRIPT_DIR/config.sh"

# Delete trimmed FASTQ files for specified SRR IDs
delete_trimmed_fastq_by_srr_list() {
	local srr_list=("$@")
	local trim_root="${TRIM_DIR_ROOT:-2_TRIMMED_SRR}"
	local deleted_count=0 skipped_count=0

	log_step "Deleting trimmed FASTQ files for ${#srr_list[@]} SRR IDs"

	for srr_id in "${srr_list[@]}"; do
		local srr_dir="${trim_root}/${srr_id}"
		if [[ -d "$srr_dir" ]]; then
			local file_count=$(find "$srr_dir" -type f \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) 2>/dev/null | wc -l)
			if [[ $file_count -gt 0 ]]; then
				rm -rf "$srr_dir"
				log_info "Deleted: $srr_dir ($file_count files)"
				((deleted_count++))
			else
				log_warn "No FASTQ files in: $srr_dir"
				((skipped_count++))
			fi
		else
			log_warn "Directory not found: $srr_dir"
			((skipped_count++))
		fi
	done
	log_info "Cleanup complete: $deleted_count directories deleted, $skipped_count skipped"
}

# Delete BAM files to save disk space
delete_bam_files() {
	local root_dir="${1:-.}"
	log_step "Deleting BAM files in $root_dir"
	find "$root_dir" -type f \( -name "*.bam" -o -name "*.bam.bai" \) -delete 2>/dev/null || true
	log_info "BAM cleanup complete"
}

# Delete SAM files
delete_sam_files() {
	local root_dir="${1:-.}"
	log_step "Deleting SAM files in $root_dir"
	find "$root_dir" -type f -name "*.sam" -delete 2>/dev/null || true
	log_info "SAM cleanup complete"
}
