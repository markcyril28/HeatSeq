#!/usr/bin/env python3
"""
Bowtie2 Alignment Algorithm Implementation (Refactored)
======================================================
Uses FM-Index (Burrows-Wheeler Transform) with dynamic programming for gapped alignment.
Implements: BWT, FM-Index, Smith-Waterman local alignment, Needleman-Wunsch global alignment.

Only basic Python: variables, loops, conditionals, functions.
Uses shared_utils for common BWT/FM-Index operations.
"""

from shared_utils import (
    build_suffix_array, build_bwt, build_c_table, build_occ_table,
    fm_backward_search, reverse_complement
)

# =============================================================================
# FM-INDEX PATTERN LOCATION
# =============================================================================

def find_pattern(pattern, suffix_array, c_table, occ):
    """Find all occurrences of pattern in reference."""
    n = len(suffix_array)
    lo, hi = fm_backward_search(pattern, c_table, occ, n)
    if lo >= hi or lo < 0:
        return []
    return sorted([suffix_array[i] for i in range(lo, hi)])


# =============================================================================
# SMITH-WATERMAN LOCAL ALIGNMENT
# =============================================================================

def compress_cigar(ops):
    """Compress CIGAR operations into standard format."""
    if len(ops) == 0:
        return ""
    cigar = ""
    current_op = ops[0]
    count = 1
    for i in range(1, len(ops)):
        if ops[i] == current_op:
            count += 1
        else:
            cigar += str(count) + current_op
            current_op = ops[i]
            count = 1
    cigar += str(count) + current_op
    return cigar


def smith_waterman(query, target, match_score=2, mismatch_penalty=-4, gap_extend=-1):
    """
    Smith-Waterman local alignment algorithm.
    Returns: (score, cigar_string, query_start, target_start)
    """
    m, n = len(query), len(target)
    
    # Initialize matrices
    H = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[0] * (n + 1) for _ in range(m + 1)]
    max_score, max_i, max_j = 0, 0, 0
    
    # Fill matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = H[i-1][j-1] + (match_score if query[i-1] == target[j-1] else mismatch_penalty)
            delete = H[i-1][j] + gap_extend
            insert = H[i][j-1] + gap_extend
            
            H[i][j] = max(0, match, delete, insert)
            
            if H[i][j] == match:
                traceback[i][j] = 1
            elif H[i][j] == delete:
                traceback[i][j] = 3
            elif H[i][j] == insert:
                traceback[i][j] = 2
            
            if H[i][j] > max_score:
                max_score, max_i, max_j = H[i][j], i, j
    
    # Traceback
    cigar_ops = []
    i, j = max_i, max_j
    while i > 0 and j > 0 and H[i][j] > 0:
        if traceback[i][j] == 1:
            cigar_ops.append("M")
            i, j = i - 1, j - 1
        elif traceback[i][j] == 3:
            cigar_ops.append("I")
            i -= 1
        elif traceback[i][j] == 2:
            cigar_ops.append("D")
            j -= 1
        else:
            break
    
    cigar_ops.reverse()
    return max_score, compress_cigar(cigar_ops), i, j


# =============================================================================
# NEEDLEMAN-WUNSCH GLOBAL ALIGNMENT
# =============================================================================

def needleman_wunsch(query, target, match_score=2, mismatch_penalty=-4, gap_open=-6, gap_extend=-1):
    """
    Needleman-Wunsch global alignment algorithm.
    Returns: (score, cigar_string)
    """
    m, n = len(query), len(target)
    NEG_INF = -999999
    
    # Initialize matrices
    H = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    traceback = [[0] * (n + 1) for _ in range(m + 1)]
    
    H[0][0] = 0
    for j in range(1, n + 1):
        H[0][j] = gap_open + j * gap_extend
    for i in range(1, m + 1):
        H[i][0] = gap_open + i * gap_extend
    
    # Fill matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = H[i-1][j-1] + (match_score if query[i-1] == target[j-1] else mismatch_penalty)
            up = H[i-1][j] + gap_extend
            left = H[i][j-1] + gap_extend
            
            H[i][j] = max(diag, up, left)
            
            if H[i][j] == diag:
                traceback[i][j] = 1
            elif H[i][j] == up:
                traceback[i][j] = 3
            else:
                traceback[i][j] = 2
    
    # Traceback
    cigar_ops = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback[i][j] == 1:
            cigar_ops.append("M")
            i, j = i - 1, j - 1
        elif j > 0 and (i == 0 or traceback[i][j] == 2):
            cigar_ops.append("D")
            j -= 1
        else:
            cigar_ops.append("I")
            i -= 1
    
    cigar_ops.reverse()
    return H[m][n], compress_cigar(cigar_ops)


# =============================================================================
# SEED EXTRACTION
# =============================================================================

def extract_seeds(read, seed_len, seed_interval):
    """Extract seeds from read at regular intervals."""
    seeds = []
    offset = 0
    while offset <= len(read) - seed_len:
        seeds.append((read[offset:offset + seed_len], offset))
        offset += seed_interval
    return seeds


# =============================================================================
# MAIN BOWTIE2 ALIGNMENT FUNCTION
# =============================================================================

def bowtie2_align(read, reference, mode="local", seed_len=22, seed_interval=15, max_mismatches=2):
    """
    Main Bowtie2 alignment function.
    mode: "local" (Smith-Waterman) or "end-to-end" (Needleman-Wunsch)
    """
    ref_with_term = reference + "$"
    
    print("Building FM-index...")
    suffix_array = build_suffix_array(ref_with_term)
    bwt = build_bwt(ref_with_term, suffix_array)
    c_table, _ = build_c_table(bwt)
    occ = build_occ_table(bwt)
    
    print("Extracting seeds and finding hits...")
    seeds = extract_seeds(read, seed_len, seed_interval)
    candidate_positions = []
    
    for seed, offset in seeds:
        positions = find_pattern(seed, suffix_array, c_table, occ)
        for pos in positions:
            read_start = pos - offset
            if 0 <= read_start <= len(reference) - len(read):
                if read_start not in candidate_positions:
                    candidate_positions.append(read_start)
    
    print("Extending candidates...")
    alignments = []
    
    for pos in candidate_positions:
        ref_region = reference[pos:pos + len(read) + 20]
        
        if mode == "local":
            score, cigar, q_start, t_start = smith_waterman(read, ref_region)
            if score > 0:
                alignments.append({"position": pos + t_start, "cigar": cigar, "score": score, "mode": "local"})
        else:
            ref_segment = reference[pos:pos + len(read)]
            if len(ref_segment) == len(read):
                score, cigar = needleman_wunsch(read, ref_segment)
                alignments.append({"position": pos, "cigar": cigar, "score": score, "mode": "end-to-end"})
    
    alignments.sort(key=lambda x: x["score"], reverse=True)
    
    for i, aln in enumerate(alignments):
        aln["mapq"] = min(60, aln["score"]) if i == 0 else max(0, 60 - 10 * i)
    
    return alignments


# =============================================================================
# PAIRED-END ALIGNMENT
# =============================================================================

def bowtie2_align_paired(read1, read2, reference, min_insert=0, max_insert=500):
    """Align paired-end reads with insert size constraints."""
    alignments1 = bowtie2_align(read1, reference, mode="local")
    read2_rc = reverse_complement(read2)
    alignments2 = bowtie2_align(read2_rc, reference, mode="local")
    
    concordant_pairs = []
    for aln1 in alignments1:
        for aln2 in alignments2:
            insert_size = abs(aln2["position"] - aln1["position"]) + len(read1)
            if min_insert <= insert_size <= max_insert:
                concordant_pairs.append({
                    "read1": aln1, "read2": aln2,
                    "insert_size": insert_size, "concordant": True,
                    "pair_score": aln1["score"] + aln2["score"]
                })
    
    if concordant_pairs:
        concordant_pairs.sort(key=lambda x: x["pair_score"], reverse=True)
        return concordant_pairs[0]
    
    return {
        "read1": alignments1[0] if alignments1 else None,
        "read2": alignments2[0] if alignments2 else None,
        "concordant": False
    }


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    reference = "ATCGATCGATCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAG"
    
    test_reads = ["ATCGATCGATCG", "ATCGATCGATCA", "TAGCTAGCTAGC"]
    
    print("=" * 60)
    print("BOWTIE2 ALIGNMENT ALGORITHM DEMONSTRATION")
    print("=" * 60)
    print("Reference length:", len(reference))
    
    for read in test_reads:
        print("-" * 50)
        print("Read:", read)
        
        print("LOCAL mode:")
        alignments = bowtie2_align(read, reference, mode="local", seed_len=8)
        for i, aln in enumerate(alignments[:3]):
            print("  Alignment", i+1, "- Pos:", aln["position"], "Score:", aln["score"], "MAPQ:", aln["mapq"])
        
        print("END-TO-END mode:")
        alignments = bowtie2_align(read, reference, mode="end-to-end", seed_len=8)
        for i, aln in enumerate(alignments[:3]):
            print("  Alignment", i+1, "- Pos:", aln["position"], "Score:", aln["score"])
