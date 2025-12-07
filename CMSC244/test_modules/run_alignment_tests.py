#!/usr/bin/env python3
"""
Alignment Test Runner
=====================
Runs all three alignment algorithms (HISAT2, Bowtie2, Salmon) with actual FASTQ input
and generates complexity reports.

Modes:
  - test: Small subset of reads for algorithm analysis (default)
  - full: Process all reads in the input file

Only basic Python: variables, loops, conditionals, functions.
"""

import os
import sys
import time

# Add test_modules to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fastq_utils import read_fastq, read_fasta, write_alignments_to_sam
from complexity_analysis import (
    create_complexity_tracker, add_measurement, measure_memory_usage,
    generate_full_report, generate_combined_comparison
)
from shared_utils import reverse_complement

# Import alignment algorithms
from hisat2_alignment import hisat2_align
from bowtie2_alignment import bowtie2_align
from salmon_saf_alignment import salmon_quantify

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input files
FASTQ_R1 = "test_inputs/SRR3884686/SRR3884686_1_val_1.fq.gz"
FASTQ_R2 = "test_inputs/SRR3884686/SRR3884686_2_val_2.fq.gz"

# Reference options (change as needed)
REFERENCE_FASTA = "test_inputs/All_Smel_Genes.fasta"

# Test mode settings (small subset for algorithm analysis)
TEST_SIZES = [10, 50, 100, 200, 500]
TEST_REF_LIMIT = 10000      # Truncate reference for test mode
TEST_TRANSCRIPT_LIMIT = 10  # Number of transcripts for Salmon test

# Full mode settings
FULL_BATCH_SIZE = 1000      # Process reads in batches for full mode

# Threading configuration (can be overridden via --threads)
NUM_THREADS = 4

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def ensure_dir(directory):
    """Create directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)


def print_header(text):
    """Print a formatted header."""
    print("")
    print("=" * 70)
    print(text)
    print("=" * 70)


def print_subheader(text):
    """Print a formatted subheader."""
    print("")
    print("-" * 50)
    print(text)
    print("-" * 50)


# =============================================================================
# HISAT2 TEST RUNNER
# =============================================================================

def run_hisat2_test(reads, reference, ref_name, output_dir):
    """
    Run HISAT2 alignment test on a set of reads.
    
    Returns:
        Tuple of (alignments, runtime, memory_mb)
    """
    print("  Running HISAT2 alignment on", len(reads), "reads...")
    
    # Measure initial memory
    mem_before = measure_memory_usage()
    
    # Run alignment
    start_time = time.time()
    
    alignments = []
    aligned_count = 0
    
    for i, read in enumerate(reads):
        seq = read['sequence']
        
        # Try forward alignment
        alns = hisat2_align(seq, reference, max_mismatches=2)
        
        # If no alignment, try reverse complement
        if len(alns) == 0:
            rc_seq = reverse_complement(seq)
            alns = hisat2_align(rc_seq, reference, max_mismatches=2)
        
        if len(alns) > 0:
            best = alns[0]
            aligned_count += 1
            alignments.append({
                'read_id': read['id'],
                'position': best['position'],
                'cigar': best['cigar'],
                'mapq': min(60, best['score']),
                'sequence': seq,
                'quality': read.get('quality', '*'),
                'unmapped': False
            })
        else:
            alignments.append({
                'read_id': read['id'],
                'sequence': seq,
                'quality': read.get('quality', '*'),
                'unmapped': True
            })
        
        # Progress indicator (suppress individual read output)
        if (i + 1) % 50 == 0:
            print("    Processed", i + 1, "/", len(reads), "reads...")
    
    end_time = time.time()
    runtime = end_time - start_time
    
    # Measure memory
    mem_after = measure_memory_usage()
    memory_used = max(0, mem_after - mem_before)
    
    print("  Aligned:", aligned_count, "/", len(reads), "reads")
    print("  Runtime:", round(runtime, 4), "seconds")
    
    return alignments, runtime, memory_used


# =============================================================================
# BOWTIE2 TEST RUNNER
# =============================================================================

def run_bowtie2_test(reads, reference, ref_name, output_dir):
    """
    Run Bowtie2 alignment test on a set of reads.
    
    Returns:
        Tuple of (alignments, runtime, memory_mb)
    """
    print("  Running Bowtie2 alignment on", len(reads), "reads...")
    
    mem_before = measure_memory_usage()
    start_time = time.time()
    
    alignments = []
    aligned_count = 0
    
    for i, read in enumerate(reads):
        seq = read['sequence']
        
        # Run local mode alignment
        alns = bowtie2_align(seq, reference, mode="local", seed_len=15)
        
        if len(alns) > 0:
            best = alns[0]
            aligned_count += 1
            alignments.append({
                'read_id': read['id'],
                'position': best['position'],
                'cigar': best['cigar'],
                'mapq': best.get('mapq', 60),
                'sequence': seq,
                'quality': read.get('quality', '*'),
                'unmapped': False
            })
        else:
            alignments.append({
                'read_id': read['id'],
                'sequence': seq,
                'quality': read.get('quality', '*'),
                'unmapped': True
            })
        
        if (i + 1) % 50 == 0:
            print("    Processed", i + 1, "/", len(reads), "reads...")
    
    end_time = time.time()
    runtime = end_time - start_time
    
    mem_after = measure_memory_usage()
    memory_used = max(0, mem_after - mem_before)
    
    print("  Aligned:", aligned_count, "/", len(reads), "reads")
    print("  Runtime:", round(runtime, 4), "seconds")
    
    return alignments, runtime, memory_used


# =============================================================================
# SALMON TEST RUNNER
# =============================================================================

def run_salmon_test(reads, transcripts, output_dir):
    """
    Run Salmon quantification test on a set of reads.
    
    Returns:
        Tuple of (tpm_results, runtime, memory_mb)
    """
    print("  Running Salmon quantification on", len(reads), "reads...")
    
    mem_before = measure_memory_usage()
    start_time = time.time()
    
    # Extract just sequences for Salmon
    read_sequences = []
    for read in reads:
        read_sequences.append(read['sequence'])
    
    # Run quantification
    tpm = salmon_quantify(read_sequences, transcripts, k=15, validate=True)
    
    end_time = time.time()
    runtime = end_time - start_time
    
    mem_after = measure_memory_usage()
    memory_used = max(0, mem_after - mem_before)
    
    # Count transcripts with expression
    expressed = sum(1 for t in tpm.values() if t > 0)
    
    print("  Transcripts with expression:", expressed, "/", len(tpm))
    print("  Runtime:", round(runtime, 4), "seconds")
    
    return tpm, runtime, memory_used


# =============================================================================
# MAIN TEST FUNCTION
# =============================================================================

def get_output_dirs(mode):
    """Get output directories based on run mode."""
    base = "test" if mode == "test" else "full"
    return {
        'hisat2': "Outputs/HISAT2/python_" + base,
        'bowtie': "Outputs/Bowtie/python_" + base,
        'salmon': "Outputs/Salmon_Saf/python_" + base,
        'combined': "Outputs/python_" + base + "_comparison"
    }


def run_test_mode(transcripts, reference, ref_name):
    """Run test mode: small subset for algorithm analysis."""
    dirs = get_output_dirs("test")
    
    for d in dirs.values():
        ensure_dir(d)
    
    # Truncate reference for test mode
    if len(reference) > TEST_REF_LIMIT:
        reference = reference[:TEST_REF_LIMIT]
        print("Truncated reference to", TEST_REF_LIMIT, "bp for test mode")
    
    # Create complexity trackers
    hisat2_tracker = create_complexity_tracker()
    hisat2_tracker['algorithm'] = 'HISAT2'
    bowtie2_tracker = create_complexity_tracker()
    bowtie2_tracker['algorithm'] = 'Bowtie2'
    salmon_tracker = create_complexity_tracker()
    salmon_tracker['algorithm'] = 'Salmon'
    
    # Prepare test transcripts for Salmon
    test_transcripts = {}
    count = 0
    for tid in transcripts:
        if count >= TEST_TRANSCRIPT_LIMIT:
            break
        test_transcripts[tid] = transcripts[tid][:5000]
        count += 1
    
    # Run tests for each size
    for num_reads in TEST_SIZES:
        print_subheader("Testing with " + str(num_reads) + " reads")
        
        reads = read_fastq(FASTQ_R1, max_reads=num_reads)
        print("Loaded", len(reads), "reads")
        
        if len(reads) == 0:
            continue
        
        # HISAT2
        print_subheader("HISAT2 Test")
        alns, runtime, memory = run_hisat2_test(reads, reference, ref_name, dirs['hisat2'])
        add_measurement(hisat2_tracker, len(reads), runtime, memory, 
                       len(reads) * len(reference), str(len(reads)) + " reads")
        sam_file = os.path.join(dirs['hisat2'], "alignments_" + str(num_reads) + ".sam")
        write_alignments_to_sam(sam_file, alns, ref_name, len(reference))
        
        # Bowtie2
        print_subheader("Bowtie2 Test")
        alns, runtime, memory = run_bowtie2_test(reads, reference, ref_name, dirs['bowtie'])
        add_measurement(bowtie2_tracker, len(reads), runtime, memory,
                       len(reads) * len(reference), str(len(reads)) + " reads")
        sam_file = os.path.join(dirs['bowtie'], "alignments_" + str(num_reads) + ".sam")
        write_alignments_to_sam(sam_file, alns, ref_name, len(reference))
        
        # Salmon
        print_subheader("Salmon Test")
        tpm, runtime, memory = run_salmon_test(reads, test_transcripts, dirs['salmon'])
        add_measurement(salmon_tracker, len(reads), runtime, memory,
                       len(reads) * len(test_transcripts), str(len(reads)) + " reads")
        tpm_file = os.path.join(dirs['salmon'], "quant_" + str(num_reads) + ".tsv")
        with open(tpm_file, 'w') as f:
            f.write("transcript_id\tTPM\n")
            for tid, val in sorted(tpm.items(), key=lambda x: x[1], reverse=True):
                f.write(tid + "\t" + str(round(val, 4)) + "\n")
    
    # Generate reports
    print_header("GENERATING COMPLEXITY REPORTS")
    generate_full_report(hisat2_tracker, dirs['hisat2'])
    generate_full_report(bowtie2_tracker, dirs['bowtie'])
    generate_full_report(salmon_tracker, dirs['salmon'])
    generate_combined_comparison([hisat2_tracker, bowtie2_tracker, salmon_tracker], dirs['combined'])
    
    return dirs


def run_full_mode(transcripts, reference, ref_name):
    """Run full mode: process all reads."""
    dirs = get_output_dirs("full")
    
    for d in dirs.values():
        ensure_dir(d)
    
    print_header("FULL RUN MODE - Processing all reads")
    
    # Load all reads
    print("Loading all reads from FASTQ...")
    reads = read_fastq(FASTQ_R1, max_reads=None)
    total_reads = len(reads)
    print("Total reads:", total_reads)
    
    if total_reads == 0:
        print("ERROR: No reads loaded")
        return dirs
    
    # HISAT2 full run
    print_subheader("HISAT2 Full Run")
    start_time = time.time()
    all_alns = []
    aligned_count = 0
    
    for i, read in enumerate(reads):
        seq = read['sequence']
        alns = hisat2_align(seq, reference, max_mismatches=2)
        if len(alns) == 0:
            alns = hisat2_align(reverse_complement(seq), reference, max_mismatches=2)
        
        if len(alns) > 0:
            best = alns[0]
            aligned_count += 1
            all_alns.append({'read_id': read['id'], 'position': best['position'],
                           'cigar': best['cigar'], 'mapq': min(60, best['score']),
                           'sequence': seq, 'quality': read.get('quality', '*'), 'unmapped': False})
        else:
            all_alns.append({'read_id': read['id'], 'sequence': seq,
                           'quality': read.get('quality', '*'), 'unmapped': True})
        
        if (i + 1) % 1000 == 0:
            print("  Processed", i + 1, "/", total_reads, "reads...")
    
    runtime = time.time() - start_time
    print("HISAT2 aligned:", aligned_count, "/", total_reads, "in", round(runtime, 2), "s")
    sam_file = os.path.join(dirs['hisat2'], "alignments_full.sam")
    write_alignments_to_sam(sam_file, all_alns, ref_name, len(reference))
    
    # Bowtie2 full run
    print_subheader("Bowtie2 Full Run")
    start_time = time.time()
    all_alns = []
    aligned_count = 0
    
    for i, read in enumerate(reads):
        seq = read['sequence']
        alns = bowtie2_align(seq, reference, mode="local", seed_len=15)
        
        if len(alns) > 0:
            best = alns[0]
            aligned_count += 1
            all_alns.append({'read_id': read['id'], 'position': best['position'],
                           'cigar': best['cigar'], 'mapq': best.get('mapq', 60),
                           'sequence': seq, 'quality': read.get('quality', '*'), 'unmapped': False})
        else:
            all_alns.append({'read_id': read['id'], 'sequence': seq,
                           'quality': read.get('quality', '*'), 'unmapped': True})
        
        if (i + 1) % 1000 == 0:
            print("  Processed", i + 1, "/", total_reads, "reads...")
    
    runtime = time.time() - start_time
    print("Bowtie2 aligned:", aligned_count, "/", total_reads, "in", round(runtime, 2), "s")
    sam_file = os.path.join(dirs['bowtie'], "alignments_full.sam")
    write_alignments_to_sam(sam_file, all_alns, ref_name, len(reference))
    
    # Salmon full run
    print_subheader("Salmon Full Run")
    start_time = time.time()
    read_seqs = [r['sequence'] for r in reads]
    tpm = salmon_quantify(read_seqs, transcripts, k=15, validate=True)
    runtime = time.time() - start_time
    
    expressed = sum(1 for t in tpm.values() if t > 0)
    print("Salmon quantified:", expressed, "/", len(tpm), "transcripts in", round(runtime, 2), "s")
    tpm_file = os.path.join(dirs['salmon'], "quant_full.tsv")
    with open(tpm_file, 'w') as f:
        f.write("transcript_id\tTPM\n")
        for tid, val in sorted(tpm.items(), key=lambda x: x[1], reverse=True):
            f.write(tid + "\t" + str(round(val, 4)) + "\n")
    
    return dirs


def run_all_tests(mode="test"):
    """Run all alignment tests based on mode."""
    
    print_header("ALIGNMENT ALGORITHM TESTING")
    print("Mode:", mode)
    print("Input FASTQ:", FASTQ_R1)
    print("Reference:", REFERENCE_FASTA)
    
    ensure_dir("logs")
    
    # Check input files
    if not os.path.exists(FASTQ_R1):
        print("ERROR: FASTQ file not found:", FASTQ_R1)
        return
    if not os.path.exists(REFERENCE_FASTA):
        print("ERROR: Reference FASTA not found:", REFERENCE_FASTA)
        return
    
    # Load reference
    print_subheader("Loading Reference")
    transcripts = read_fasta(REFERENCE_FASTA)
    print("Loaded", len(transcripts), "sequences")
    
    ref_name = list(transcripts.keys())[0]
    reference = transcripts[ref_name]
    print("Using reference:", ref_name, "length:", len(reference))
    
    # Run based on mode
    if mode == "full":
        dirs = run_full_mode(transcripts, reference, ref_name)
    else:
        dirs = run_test_mode(transcripts, reference, ref_name)
    
    # Summary
    print_header("TEST COMPLETE")
    print("Mode:", mode)
    print("Output directories:")
    for k, v in dirs.items():
        print(" ", k + ":", v)


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    # Parse command line arguments
    mode = "test"
    for i, arg in enumerate(sys.argv):
        if arg == "--mode" and i + 1 < len(sys.argv):
            mode = sys.argv[i + 1]
        if arg == "--threads" and i + 1 < len(sys.argv):
            NUM_THREADS = int(sys.argv[i + 1])
    
    print("Using", NUM_THREADS, "threads")
    run_all_tests(mode)
