#!/usr/bin/env python3
"""
Salmon Selective Alignment (SAF) Algorithm Implementation (Refactored)
=====================================================================
Uses quasi-mapping with minimizers, selective alignment validation,
and EM algorithm for transcript quantification.

"""

import math
from shared_utils import reverse_complement, hash_kmer, get_minimizers

# =============================================================================
# BUILD SALMON INDEX
# =============================================================================

def build_salmon_index(transcripts, k=31):
    """
    Build quasi-index from transcripts.
    Index maps minimizer hashes to (transcript_id, position, orientation).
    """
    index = {}
    transcript_lengths = {}
    
    for tid, seq in transcripts.items():
        transcript_lengths[tid] = len(seq)
        
        # Index forward strand
        for h, pos in get_minimizers(seq, k, 10):
            if h not in index:
                index[h] = []
            index[h].append((tid, pos, "+"))
        
        # Index reverse complement
        rc_seq = reverse_complement(seq)
        for h, pos in get_minimizers(rc_seq, k, 10):
            if h not in index:
                index[h] = []
            index[h].append((tid, len(seq) - pos - k, "-"))
    
    return index, transcript_lengths


def add_decoys_to_index(index, decoy_sequences, k=31):
    """Add decoy sequences to index for filtering spurious mappings."""
    decoy_ids = []
    for i, seq in enumerate(decoy_sequences):
        decoy_id = "decoy_" + str(i)
        decoy_ids.append(decoy_id)
        for h, pos in get_minimizers(seq, k, 10):
            if h not in index:
                index[h] = []
            index[h].append((decoy_id, pos, "+"))
    return decoy_ids


# =============================================================================
# QUASI-MAPPING
# =============================================================================

def quasi_map(read, index, transcript_lengths, k=31, min_hits=3):
    """
    Perform quasi-mapping of a read to transcripts, then
    Returns list of candidate mappings based on minimizer hits.
    """
    minimizers = get_minimizers(read, k, 10)
    hit_counts = {}
    
    for read_hash, read_pos in minimizers:
        if read_hash in index:
            for tid, ref_pos, orientation in index[read_hash]:
                key = (tid, orientation)
                if key not in hit_counts:
                    hit_counts[key] = []
                hit_counts[key].append((read_pos, ref_pos))
    
    mappings = []
    for (tid, orientation), hits in hit_counts.items():
        if len(hits) >= min_hits and not tid.startswith("decoy_"):
            offsets = sorted([ref_pos - read_pos for read_pos, ref_pos in hits])
            position = offsets[len(offsets) // 2]
            coverage = len(hits) / max(1, len(minimizers))
            
            mappings.append({
                "transcript_id": tid,
                "position": position,
                "orientation": orientation,
                "num_hits": len(hits),
                "coverage": coverage
            })
    
    mappings.sort(key=lambda x: x["coverage"], reverse=True)
    return mappings


# =============================================================================
# SELECTIVE ALIGNMENT (VALIDATION)
# =============================================================================

def banded_alignment(query, target, bandwidth=10, match=2, mismatch=-4, gap=-2):
    """Banded semi-global alignment for fast scoring."""
    m, n = len(query), len(target)
    NEG_INF = -999999
    
    dp = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    for j in range(n + 1):
        dp[0][j] = 0
    
    for i in range(1, m + 1):
        j_start = max(1, i - bandwidth)
        j_end = min(n + 1, i + bandwidth)
        
        for j in range(j_start, j_end):
            diag = dp[i-1][j-1] + (match if query[i-1] == target[j-1] else mismatch)
            up = dp[i-1][j] + gap if dp[i-1][j] > NEG_INF else NEG_INF
            left = dp[i][j-1] + gap if dp[i][j-1] > NEG_INF else NEG_INF
            dp[i][j] = max(diag, up, left)
    
    return max(dp[m])


def selective_align(read, mapping, transcripts):
    """Validate quasi-mapping using alignment score."""
    tid = mapping["transcript_id"]
    if tid not in transcripts:
        return None
    
    transcript = transcripts[tid]
    start = max(0, mapping["position"])
    end = min(len(transcript), mapping["position"] + len(read) + 20)
    ref_region = transcript[start:end]
    
    query = reverse_complement(read) if mapping["orientation"] == "-" else read
    score = banded_alignment(query, ref_region)
    
    if score > 0:
        return {
            "transcript_id": tid,
            "position": start,
            "score": score,
            "orientation": mapping["orientation"],
            "cigar": str(len(read)) + "M"
        }
    return None


# =============================================================================
# EXPECTATION-MAXIMIZATION (EM) FOR QUANTIFICATION
# =============================================================================

def em_quantify(alignments_per_read, transcript_lengths, max_iter=1000, tolerance=1e-8):
    """
    Run EM algorithm to estimate transcript abundances.
    Returns: dict mapping transcript_id -> TPM
    """
    transcripts = list(transcript_lengths.keys())
    n_transcripts = len(transcripts)
    
    if n_transcripts == 0:
        return {}
    
    tid_to_idx = {tid: i for i, tid in enumerate(transcripts)}
    theta = [1.0 / n_transcripts] * n_transcripts
    
    # Build read-to-transcript mapping
    read_mappings = []
    for read_alns in alignments_per_read:
        mapping = {}
        for aln in read_alns:
            if aln and aln["transcript_id"] in tid_to_idx:
                mapping[tid_to_idx[aln["transcript_id"]]] = aln["score"]
        if mapping:
            read_mappings.append(mapping)
    
    if not read_mappings:
        return {tid: 0.0 for tid in transcripts}
    
    # EM iterations
    for iteration in range(max_iter):
        theta_old = theta[:]
        expected_counts = [0.0] * n_transcripts
        
        for mapping in read_mappings:
            probs = [0.0] * n_transcripts
            total = 0.0
            
            for idx, score in mapping.items():
                prob = theta[idx] * math.exp(score / 10.0)
                probs[idx] = prob
                total += prob
            
            if total > 0:
                for idx in mapping:
                    expected_counts[idx] += probs[idx] / total
        
        # M-step with effective length normalization
        fragment_length = 150
        effective_lengths = [max(1, transcript_lengths[tid] - fragment_length + 1) for tid in transcripts]
        
        theta_sum = 0.0
        for i in range(n_transcripts):
            theta[i] = expected_counts[i] / effective_lengths[i]
            theta_sum += theta[i]
        
        if theta_sum > 0:
            theta = [t / theta_sum for t in theta]
        
        # Check convergence
        if max(abs(theta[i] - theta_old[i]) for i in range(n_transcripts)) < tolerance:
            print("EM converged after", iteration + 1, "iterations")
            break
    
    return {transcripts[i]: theta[i] * 1e6 for i in range(n_transcripts)}


# =============================================================================
# MAIN SALMON QUANTIFICATION FUNCTION
# =============================================================================

def salmon_quantify(reads, transcripts, k=31, validate=True, decoy_sequences=None):
    """
    Main Salmon quantification function.
    Returns: dict mapping transcript_id -> TPM
    """
    print("Building Salmon index...")
    index, transcript_lengths = build_salmon_index(transcripts, k)
    
    if decoy_sequences:
        print("Adding decoys to index...")
        add_decoys_to_index(index, decoy_sequences, k)
    
    print("Mapping reads...")
    all_alignments = []
    
    for i, read in enumerate(reads):
        mappings = quasi_map(read, index, transcript_lengths, k)
        
        if validate and mappings:
            validated = []
            for mapping in mappings[:10]:
                aln = selective_align(read, mapping, transcripts)
                if aln:
                    validated.append(aln)
            all_alignments.append(validated)
        else:
            pseudo_alns = [{
                "transcript_id": m["transcript_id"],
                "position": m["position"],
                "score": m["coverage"] * 100,
                "orientation": m["orientation"],
                "cigar": str(len(read)) + "M"
            } for m in mappings]
            all_alignments.append(pseudo_alns)
        
        if (i + 1) % 100 == 0:
            print("  Processed", i + 1, "reads...")
    
    print("Running EM quantification...")
    return em_quantify(all_alignments, transcript_lengths)


# =============================================================================
# PAIRED-END QUANTIFICATION
# =============================================================================

def salmon_quantify_paired(reads1, reads2, transcripts, k=31, min_frag=100, max_frag=500):
    """Quantify from paired-end reads with fragment length constraints."""
    print("Building Salmon index...")
    index, transcript_lengths = build_salmon_index(transcripts, k)
    
    print("Mapping paired-end reads...")
    all_alignments = []
    
    for i in range(len(reads1)):
        r1, r2_rc = reads1[i], reverse_complement(reads2[i])
        mappings1 = quasi_map(r1, index, transcript_lengths, k)
        mappings2 = quasi_map(r2_rc, index, transcript_lengths, k)
        
        concordant = []
        for m1 in mappings1:
            for m2 in mappings2:
                if m1["transcript_id"] == m2["transcript_id"]:
                    frag_len = abs(m2["position"] - m1["position"]) + len(r1)
                    if min_frag <= frag_len <= max_frag:
                        aln1 = selective_align(r1, m1, transcripts)
                        aln2 = selective_align(r2_rc, m2, transcripts)
                        if aln1 and aln2:
                            concordant.append({
                                "transcript_id": m1["transcript_id"],
                                "position": min(m1["position"], m2["position"]),
                                "score": aln1["score"] + aln2["score"],
                                "orientation": m1["orientation"],
                                "cigar": aln1["cigar"] + "..." + aln2["cigar"]
                            })
        all_alignments.append(concordant)
    
    print("Running EM quantification...")
    return em_quantify(all_alignments, transcript_lengths)


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    transcripts = {
        "transcript_A": "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "transcript_B": "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
        "transcript_C": "GCTAGCTAGCATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAG",
    }
    
    decoys = ["NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"]
    
    test_reads = [
        "ATCGATCGATCGATCGATCGATCG",
        "ATCGATCGATCGATCGATCGATCG",
        "TAGCTAGCTAGCTAGCTAGCTAGC",
        "GCTAGCTAGCATCGATCGATCGAT",
        "GCTAGCTAGCATCGATCGATCGAT",
        "GCTAGCTAGCATCGATCGATCGAT",
    ]
    
    print("=" * 60)
    print("SALMON SELECTIVE ALIGNMENT QUANTIFICATION DEMONSTRATION")
    print("=" * 60)
    print("Transcripts:", len(transcripts), "Reads:", len(test_reads))
    
    tpm = salmon_quantify(test_reads, transcripts, k=15, validate=True, decoy_sequences=decoys)
    
    print("\nTRANSCRIPT QUANTIFICATION RESULTS (TPM)")
    print("-" * 40)
    for tid, abundance in sorted(tpm.items(), key=lambda x: x[1], reverse=True):
        print("  " + tid + ": " + str(round(abundance, 2)) + " TPM")
    print("Total TPM:", round(sum(tpm.values()), 2))
