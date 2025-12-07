#!/usr/bin/env python3
"""
Shared Utilities for Alignment Algorithms
==========================================
Common functions used by HISAT2, Bowtie2, and Salmon aligners.
Implements: BWT, FM-Index construction, sequence utilities.

Only basic Python: variables, loops, conditionals, functions.
"""

# =============================================================================
# DNA SEQUENCE UTILITIES
# =============================================================================

def reverse_complement(seq):
    """Compute reverse complement of DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    result = ""
    for base in reversed(seq):
        result += complement.get(base, "N")
    return result


# =============================================================================
# BURROWS-WHEELER TRANSFORM (BWT) CONSTRUCTION
# =============================================================================

def build_suffix_array(text):
    """
    Build suffix array by sorting all suffixes of text.
    Suffix array stores starting positions of sorted suffixes.
    Time: O(n^2 log n), Space: O(n)
    """
    n = len(text)
    suffixes = []
    for i in range(n):
        suffixes.append((text[i:], i))
    suffixes.sort(key=lambda x: x[0])
    suffix_array = []
    for suffix, index in suffixes:
        suffix_array.append(index)
    return suffix_array


def build_bwt(text, suffix_array):
    """
    Construct BWT from suffix array.
    BWT[i] = character before the i-th suffix in sorted order.
    """
    n = len(text)
    bwt = ""
    for i in range(n):
        sa_index = suffix_array[i]
        if sa_index == 0:
            bwt += text[n - 1]
        else:
            bwt += text[sa_index - 1]
    return bwt


# =============================================================================
# FM-INDEX CONSTRUCTION
# =============================================================================

def build_c_table(bwt):
    """
    Build C table: cumulative count of characters lexicographically smaller.
    C[c] = starting position of character c in first column of BWT matrix.
    """
    counts = {}
    for c in bwt:
        if c not in counts:
            counts[c] = 0
        counts[c] += 1
    
    sorted_chars = sorted(counts.keys())
    c_table = {}
    total = 0
    for c in sorted_chars:
        c_table[c] = total
        total += counts[c]
    return c_table, counts


def build_occ_table(bwt):
    """
    Build occurrence table for rank queries.
    occ[c][i] = number of occurrences of c in bwt[0:i].
    """
    n = len(bwt)
    alphabet = []
    for c in bwt:
        if c not in alphabet:
            alphabet.append(c)
    
    occ = {}
    for c in alphabet:
        occ[c] = [0] * (n + 1)
    
    for i in range(n):
        for c in alphabet:
            occ[c][i + 1] = occ[c][i]
        curr_char = bwt[i]
        occ[curr_char][i + 1] += 1
    return occ


def fm_backward_search(pattern, c_table, occ, bwt_len):
    """
    Perform backward search on FM-index.
    Returns (top, bottom) range of suffix array, or (-1, -1) if not found.
    """
    top = 0
    bottom = bwt_len
    
    for i in range(len(pattern) - 1, -1, -1):
        c = pattern[i]
        if c not in c_table:
            return -1, -1
        top = c_table[c] + (occ[c][top] if top > 0 else 0)
        bottom = c_table[c] + occ[c][bottom]
        if top >= bottom:
            return -1, -1
    return top, bottom


# =============================================================================
# K-MER UTILITIES (FOR SALMON)
# =============================================================================

def hash_kmer(kmer):
    """
    Hash a k-mer using polynomial rolling hash.
    Maps A=0, C=1, G=2, T=3.
    """
    base_map = {"A": 0, "C": 1, "G": 2, "T": 3}
    h = 0
    for c in kmer:
        h = h * 4 + base_map.get(c, 0)
    return h


def get_minimizers(sequence, k, w):
    """
    Extract minimizers from sequence.
    Minimizer = smallest k-mer hash in each window of size w.
    """
    if len(sequence) < k:
        return []
    
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        kmers.append((hash_kmer(kmer), i))
    
    minimizers = []
    last_minimizer = None
    
    for i in range(len(kmers) - w + 1):
        window = kmers[i:i + w]
        min_kmer = min(window, key=lambda x: x[0])
        if last_minimizer is None or min_kmer != last_minimizer:
            minimizers.append(min_kmer)
            last_minimizer = min_kmer
    return minimizers


# =============================================================================
# ALIGNMENT SCORING
# =============================================================================

def score_match(a, b, match=2, mismatch=-4):
    """Simple scoring for base pair comparison."""
    return match if a == b else mismatch


def compute_alignment_score(seq1, seq2, start_pos, match=2, mismatch=-4, gap=-6):
    """
    Compute simple alignment score between two sequences.
    Used for validation in multiple aligners.
    """
    score = 0
    min_len = min(len(seq1), len(seq2))
    for i in range(min_len):
        if seq1[i] == seq2[i]:
            score += match
        else:
            score += mismatch
    return score
