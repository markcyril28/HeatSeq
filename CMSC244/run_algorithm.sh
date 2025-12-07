#!/bin/bash
# ==============================================================================
# CMSC244 - Alignment Algorithm Test Runner
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ==============================================================================
# CONFIGURATION (TOGGLES)
# ==============================================================================

CLEAR_LOGS=true       # Clear logs folder before running
CLEAR_OUTPUTS=true    # Clear output folders before running
RUN_TEST=true         # Run test mode (small subset for algorithm analysis)
RUN_FULL=true         # Run full mode (all reads in input file)
NUM_THREADS=32         # Number of threads for alignment

# ==============================================================================
# INPUT/OUTPUT PATHS
# ==============================================================================

FASTQ_R1="test_inputs/SRR3884686_1_val_1.fq.gz"
FASTQ_R2="test_inputs/SRR3884686_2_val_2.fq.gz"
REFERENCE="test_inputs/Eggplant_V4.1_transcripts.function.fa"
#REFERENCE="test_inputs/All_Smel_Genes.fasta"

OUTPUT_HISAT2="Outputs/HISAT2"
OUTPUT_BOWTIE="Outputs/Bowtie"
OUTPUT_SALMON="Outputs/Salmon_Saf"

CONDA_ENV="cmsc"

# ==============================================================================
# SOURCE UTILITIES
# ==============================================================================

source test_modules/logging_utils.sh
source test_modules/bash_utils.sh

# ==============================================================================
# MAIN
# ==============================================================================

main() {
    print_header "CMSC244 ALIGNMENT ALGORITHM TEST"
    echo "Date: $(date)"
    echo "Clear logs: $CLEAR_LOGS | Clear outputs: $CLEAR_OUTPUTS"
    echo "Run Test: $RUN_TEST | Run Full: $RUN_FULL | Threads: $NUM_THREADS"
    
    # Clear folders if toggled
    if [[ "$CLEAR_LOGS" == "true" ]] && [[ -d "logs" ]]; then
        rm -rf logs/* 2>/dev/null || true
        echo "Cleared: logs/"
    fi
    
    if [[ "$CLEAR_OUTPUTS" == "true" ]]; then
        for dir in "$OUTPUT_HISAT2" "$OUTPUT_BOWTIE" "$OUTPUT_SALMON"; do
            [[ -d "$dir" ]] && rm -rf "$dir"/* 2>/dev/null || true
        done
        rm -f Outputs/*.png Outputs/*.txt Outputs/*.csv 2>/dev/null || true
        echo "Cleared: Outputs/"
    fi
    
    setup_logging false
    check_inputs "$FASTQ_R1" "$FASTQ_R2" "$REFERENCE" || exit 1
    create_directories "$OUTPUT_HISAT2" "$OUTPUT_BOWTIE" "$OUTPUT_SALMON" "logs"
    
    # Run tests
    if [[ "$RUN_TEST" == "true" ]]; then
        activate_conda_env "$CONDA_ENV" "setup_cmsc244.sh"
        print_header "RUNNING TEST MODE"
        run_with_space_time_log --input "test_inputs" --output "Outputs" \
            python test_modules/run_alignment_tests.py --mode test --threads "$NUM_THREADS"
    fi
    
    if [[ "$RUN_FULL" == "true" ]]; then
        activate_conda_env "$CONDA_ENV" "setup_cmsc244.sh"
        print_header "RUNNING FULL MODE"
        run_with_space_time_log --input "test_inputs" --output "Outputs" \
            python test_modules/run_alignment_tests.py --mode full --threads "$NUM_THREADS"
    fi
    
    print_header "COMPLETE"
    echo "Results: $OUTPUT_HISAT2, $OUTPUT_BOWTIE, $OUTPUT_SALMON"
    echo "Logs: logs/"
}

main
