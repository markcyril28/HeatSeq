#!/usr/bin/env python3
"""
HISAT2 Alignment Algorithm Implementation (Refactored)
======================================================
Uses Hierarchical Graph FM Index (HGFM) built on Burrows-Wheeler Transform (BWT).
Implements: BWT, FM-Index, backward search, and splice-aware alignment.

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

def locate_pattern(pattern, suffix_array, c_table, occ):
    """Find all positions of pattern in reference."""
    n = len(suffix_array)
    top, bottom = fm_backward_search(pattern, c_table, occ, n)
    if top >= bottom or top < 0:
        return []
    positions = []
    for i in range(top, bottom):
        positions.append(suffix_array[i])
    return sorted(positions)


# =============================================================================
# SEED-AND-EXTEND WITH MISMATCHES
# =============================================================================

def count_mismatches(seq1, seq2):
    """Count number of mismatches between two sequences."""
    mismatches = 0
    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] != seq2[i]:
            mismatches += 1
    return mismatches


def seed_and_extend(read, reference, suffix_array, c_table, occ, max_mismatches):
    """
    Seed-and-extend strategy for approximate matching.
    1. Extract seeds (short exact matches)
    2. Find seed positions using FM-index
    3. Extend seeds and count mismatches
    """
    alignments = []
    read_len = len(read)
    ref_len = len(reference)
    
    seed_len = max(8, read_len // (max_mismatches + 1))
    seed_positions = list(range(0, read_len - seed_len + 1, seed_len))
    checked_positions = {}
    
    for seed_offset in seed_positions:
        seed = read[seed_offset:seed_offset + seed_len]
        hit_positions = locate_pattern(seed, suffix_array, c_table, occ)
        
        for hit_pos in hit_positions:
            read_start = hit_pos - seed_offset
            if read_start < 0 or read_start + read_len > ref_len:
                continue
            if read_start in checked_positions:
                continue
            checked_positions[read_start] = True
            
            ref_segment = reference[read_start:read_start + read_len]
            mm_count = count_mismatches(read, ref_segment)
            
            if mm_count <= max_mismatches:
                alignments.append({
                    "position": read_start,
                    "cigar": str(read_len) + "M",
                    "mismatches": mm_count,
                    "score": read_len - mm_count,
                    "spliced": False
                })
    return alignments


# =============================================================================
# SPLICE-AWARE ALIGNMENT (HISAT2 SPECIFIC)
# =============================================================================

def check_canonical_splice_site(reference, donor_pos, acceptor_pos):
    """Check for canonical GT-AG splice site signals."""
    if donor_pos + 2 > len(reference) or acceptor_pos < 2:
        return False
    donor = reference[donor_pos:donor_pos + 2]
    acceptor = reference[acceptor_pos - 2:acceptor_pos]
    return donor == "GT" and acceptor == "AG"


def spliced_alignment(read, reference, suffix_array, c_table, occ):
    """
    Attempt spliced alignment for reads spanning introns.
    1. Split read at different positions
    2. Align segments independently
    3. Check for valid intron (GT-AG splice sites)
    """
    alignments = []
    read_len = len(read)
    min_anchor = 8
    max_intron = 500000
    min_intron = 50
    
    for split_pos in range(min_anchor, read_len - min_anchor):
        left_segment = read[:split_pos]
        right_segment = read[split_pos:]
        
        left_positions = locate_pattern(left_segment, suffix_array, c_table, occ)
        right_positions = locate_pattern(right_segment, suffix_array, c_table, occ)
        
        for left_pos in left_positions:
            left_end = left_pos + len(left_segment)
            for right_pos in right_positions:
                intron_length = right_pos - left_end
                
                if min_intron <= intron_length <= max_intron:
                    if check_canonical_splice_site(reference, left_end, right_pos):
                        cigar = str(len(left_segment)) + "M" + str(intron_length) + "N" + str(len(right_segment)) + "M"
                        alignments.append({
                            "position": left_pos,
                            "cigar": cigar,
                            "mismatches": 0,
                            "score": read_len,
                            "spliced": True,
                            "intron_start": left_end,
                            "intron_end": right_pos
                        })
    return alignments


# =============================================================================
# MAIN HISAT2 ALIGNMENT FUNCTION
# =============================================================================

def hisat2_align(read, reference, max_mismatches=2):
    """
    Main HISAT2 alignment function.
    1. Build FM-index from reference
    2. Try exact matching
    3. Try approximate matching with mismatches
    4. Try spliced alignment
    """
    ref_with_term = reference + "$"
    
    print("Building FM-index...")
    suffix_array = build_suffix_array(ref_with_term)
    bwt = build_bwt(ref_with_term, suffix_array)
    c_table, _ = build_c_table(bwt)
    occ = build_occ_table(bwt)
    
    print("Searching for exact matches...")
    exact_positions = locate_pattern(read, suffix_array, c_table, occ)
    
    alignments = []
    for pos in exact_positions:
        alignments.append({
            "position": pos,
            "cigar": str(len(read)) + "M",
            "mismatches": 0,
            "score": len(read),
            "spliced": False
        })
    
    if len(alignments) == 0:
        print("Trying approximate matching...")
        alignments = seed_and_extend(read, reference, suffix_array, c_table, occ, max_mismatches)
    
    if len(alignments) == 0:
        print("Trying spliced alignment...")
        alignments = spliced_alignment(read, reference, suffix_array, c_table, occ)
    
    alignments.sort(key=lambda x: x["score"], reverse=True)
    return alignments


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    reference = "ATCGATCGATCGATCGATCGATCGATCG" + "GT" + "N" * 100 + "AG" + "TAGCTAGCTAGCTAGCTAGCTAGCTAG"
    
    test_reads = [
        "ATCGATCGATCG",
        "ATCGATCGATCA",
        "TAGCTAGCTAGC",
    ]
    
    print("=" * 60)
    print("HISAT2 ALIGNMENT ALGORITHM DEMONSTRATION")
    print("=" * 60)
    print("Reference length:", len(reference))
    
    for read in test_reads:
        print("-" * 50)
        print("Read:", read)
        alignments = hisat2_align(read, reference, max_mismatches=2)
        
        if len(alignments) > 0:
            print("Alignments found:", len(alignments))
            for i, aln in enumerate(alignments):
                print("  Alignment", i + 1, "- Pos:", aln["position"], "CIGAR:", aln["cigar"], "Score:", aln["score"])
        else:
            print("No alignments found")
