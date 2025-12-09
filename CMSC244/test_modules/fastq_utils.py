#!/usr/bin/env python3
"""
FASTQ Utilities for Alignment Testing
=====================================
Provides functions to read FASTQ files (gzipped or plain) and FASTA reference files.
"""

import gzip
import os

# =============================================================================
# STEP 1: READ FASTQ FILES
# =============================================================================

def read_fastq(filepath, max_reads=None):
    """
    Read sequences from a FASTQ file (supports .gz compression).
    
    Args:
        filepath: Path to FASTQ file (.fq, .fastq, .fq.gz, .fastq.gz)
        max_reads: Maximum number of reads to return (None = all)
    
    Returns:
        List of dictionaries with 'id', 'sequence', 'quality' keys
    """
    reads = []
    
    # Determine if file is gzipped
    is_gzipped = filepath.endswith('.gz')
    
    # Open file appropriately
    if is_gzipped:
        file_handle = gzip.open(filepath, 'rt')
    else:
        file_handle = open(filepath, 'r')
    
    try:
        line_num = 0
        current_read = {}
        
        for line in file_handle:
            line = line.strip()
            position = line_num % 4
            
            if position == 0:
                # Header line (starts with @)
                current_read = {'id': line[1:]}  # Remove @ symbol
            elif position == 1:
                # Sequence line
                current_read['sequence'] = line
            elif position == 2:
                # Plus line (ignored)
                pass
            elif position == 3:
                # Quality line
                current_read['quality'] = line
                reads.append(current_read)
                
                # Check if we've reached max reads
                if max_reads is not None and len(reads) >= max_reads:
                    break
            
            line_num += 1
    
    finally:
        file_handle.close()
    
    return reads


def count_fastq_reads(filepath):
    """
    Count total number of reads in a FASTQ file without loading all into memory.
    """
    is_gzipped = filepath.endswith('.gz')
    
    if is_gzipped:
        file_handle = gzip.open(filepath, 'rt')
    else:
        file_handle = open(filepath, 'r')
    
    count = 0
    try:
        for i, _ in enumerate(file_handle):
            pass
        count = (i + 1) // 4
    finally:
        file_handle.close()
    
    return count


def read_fastq_paired(filepath1, filepath2, max_reads=None):
    """
    Read paired-end FASTQ files.
    
    Returns:
        List of tuples (read1, read2) where each read is a dict
    """
    reads1 = read_fastq(filepath1, max_reads)
    reads2 = read_fastq(filepath2, max_reads)
    
    # Pair them up
    paired = []
    min_len = min(len(reads1), len(reads2))
    
    for i in range(min_len):
        paired.append((reads1[i], reads2[i]))
    
    return paired


# =============================================================================
# STEP 2: READ FASTA REFERENCE FILES
# =============================================================================

def read_fasta(filepath):
    """
    Read sequences from a FASTA file.
    
    Returns:
        Dictionary mapping sequence_id -> sequence
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    # Handle gzipped files
    is_gzipped = filepath.endswith('.gz')
    if is_gzipped:
        file_handle = gzip.open(filepath, 'rt')
    else:
        file_handle = open(filepath, 'r')
    
    try:
        for line in file_handle:
            line = line.strip()
            
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                
                # Start new sequence
                current_id = line[1:].split()[0]  # Get first word after >
                current_seq = []
            else:
                # Sequence line
                current_seq.append(line.upper())
        
        # Save last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    
    finally:
        file_handle.close()
    
    return sequences


def read_fasta_single(filepath):
    """
    Read a FASTA file and return concatenated sequence (for single reference).
    
    Returns:
        Tuple of (sequence_id, sequence_string)
    """
    sequences = read_fasta(filepath)
    
    if len(sequences) == 0:
        return None, ""
    
    # Return first sequence
    first_id = list(sequences.keys())[0]
    return first_id, sequences[first_id]


# =============================================================================
# STEP 3: WRITE OUTPUT FILES
# =============================================================================

def write_sam_header(output_file, reference_name, reference_length):
    """
    Write SAM format header.
    """
    with open(output_file, 'w') as f:
        f.write("@HD\tVN:1.6\tSO:unsorted\n")
        f.write("@SQ\tSN:" + reference_name + "\tLN:" + str(reference_length) + "\n")
        f.write("@PG\tID:test_aligner\tPN:test_aligner\tVN:1.0\n")


def write_sam_alignment(output_file, read_id, flag, ref_name, position, mapq, cigar, sequence, quality):
    """
    Append a SAM alignment record to file.
    """
    with open(output_file, 'a') as f:
        # SAM fields: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
        line = "\t".join([
            read_id,
            str(flag),
            ref_name,
            str(position + 1),  # SAM uses 1-based positions
            str(mapq),
            cigar,
            "*",  # RNEXT
            "0",  # PNEXT
            "0",  # TLEN
            sequence,
            quality
        ])
        f.write(line + "\n")


def write_alignments_to_sam(output_file, alignments, ref_name, ref_length):
    """
    Write all alignments to SAM file.
    
    Args:
        alignments: List of dicts with alignment info
        ref_name: Reference sequence name
        ref_length: Reference sequence length
    """
    write_sam_header(output_file, ref_name, ref_length)
    
    for aln in alignments:
        flag = 0
        if aln.get('unmapped', False):
            flag = 4
        
        write_sam_alignment(
            output_file,
            aln.get('read_id', 'unknown'),
            flag,
            ref_name if not aln.get('unmapped', False) else "*",
            aln.get('position', 0),
            aln.get('mapq', 255),
            aln.get('cigar', '*'),
            aln.get('sequence', '*'),
            aln.get('quality', '*')
        )


# =============================================================================
# STEP 4: UTILITY FUNCTIONS
# =============================================================================

def get_file_size_mb(filepath):
    """Get file size in megabytes."""
    if os.path.exists(filepath):
        size_bytes = os.path.getsize(filepath)
        return size_bytes / (1024 * 1024)
    return 0


def count_reads_in_fastq(filepath):
    """Count total reads in a FASTQ file without loading all into memory."""
    count = 0
    
    is_gzipped = filepath.endswith('.gz')
    if is_gzipped:
        file_handle = gzip.open(filepath, 'rt')
    else:
        file_handle = open(filepath, 'r')
    
    try:
        for line in file_handle:
            count += 1
    finally:
        file_handle.close()
    
    return count // 4  # 4 lines per read in FASTQ


def count_sequences_in_fasta(filepath):
    """Count sequences in a FASTA file."""
    count = 0
    
    is_gzipped = filepath.endswith('.gz')
    if is_gzipped:
        file_handle = gzip.open(filepath, 'rt')
    else:
        file_handle = open(filepath, 'r')
    
    try:
        for line in file_handle:
            if line.startswith('>'):
                count += 1
    finally:
        file_handle.close()
    
    return count


# =============================================================================
# STEP 5: EXAMPLE/TEST
# =============================================================================

if __name__ == "__main__":
    print("FASTQ Utilities Test")
    print("=" * 50)
    
    # Test paths
    test_fastq1 = "test_inputs/SRR3884686_1_val_1.fq.gz"
    test_fastq2 = "test_inputs/SRR3884686_2_val_2.fq.gz"
    test_fasta = "test_inputs/All_Smel_Genes.fasta"
    
    # Test FASTQ reading
    if os.path.exists(test_fastq1):
        print("\nReading FASTQ file:", test_fastq1)
        reads = read_fastq(test_fastq1, max_reads=5)
        print("  Read", len(reads), "reads")
        if len(reads) > 0:
            print("  First read ID:", reads[0]['id'])
            print("  First read length:", len(reads[0]['sequence']))
    
    # Test FASTA reading
    if os.path.exists(test_fasta):
        print("\nReading FASTA file:", test_fasta)
        sequences = read_fasta(test_fasta)
        print("  Found", len(sequences), "sequences")
        if len(sequences) > 0:
            first_id = list(sequences.keys())[0]
            print("  First sequence ID:", first_id)
            print("  First sequence length:", len(sequences[first_id]))
