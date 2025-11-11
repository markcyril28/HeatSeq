#!/bin/bash
################################################################################
# Script: fix_gtf_for_stringtie.sh
# Description: Fixes GTF file for StringTie compatibility by:
#              1. Adding missing gene features (StringTie requires them)
#              2. Fixing chromosome/sequence name mismatch between FASTA and GTF
#              3. Creating proper mapping from Chr01-Chr12 to gene sequences
################################################################################

set -euo pipefail

# Input files
INPUT_GTF="All_Smel_Genes_Full_Name_reformatted.gtf"
INPUT_FASTA="All_Smel_Genes.fasta"
OUTPUT_GTF="All_Smel_Genes_Full_Name_FIXED.gtf"
MAPPING_FILE="chr_to_sequence_mapping.tsv"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== GTF Fix Script for StringTie ===${NC}"
echo ""

# Check input files exist
if [[ ! -f "$INPUT_GTF" ]]; then
    echo -e "${RED}ERROR: GTF file not found: $INPUT_GTF${NC}"
    exit 1
fi

if [[ ! -f "$INPUT_FASTA" ]]; then
    echo -e "${RED}ERROR: FASTA file not found: $INPUT_FASTA${NC}"
    exit 1
fi

echo -e "${YELLOW}Step 1: Creating chromosome to sequence mapping...${NC}"
# Extract mapping from FASTA headers
# Format: >SmelDMP_01.730     SMEL4.1_01g000730.1.01  gene=SMEL4.1_01g000730.1
awk '/^>/ {
    # Extract sequence name (first field after >)
    seq_name = $1
    gsub(/>/, "", seq_name)
    
    # Extract gene_id
    for(i=2; i<=NF; i++) {
        if($i ~ /^gene=/) {
            gene_id = $i
            gsub(/gene=/, "", gene_id)
            print gene_id "\t" seq_name
            break
        }
    }
}' "$INPUT_FASTA" > "$MAPPING_FILE"

echo -e "${GREEN}✓ Created mapping file: $MAPPING_FILE${NC}"
echo "  Total mappings: $(wc -l < "$MAPPING_FILE")"
echo ""

echo -e "${YELLOW}Step 2: Fixing GTF - Adding gene features and updating references...${NC}"

# Create fixed GTF with:
# 1. Gene features added
# 2. Chromosome names replaced with sequence names
awk -v mapping_file="$MAPPING_FILE" '
BEGIN {
    OFS="\t"
    
    # Load gene_id to sequence_name mapping
    while((getline < mapping_file) > 0) {
        gene_to_seq[$1] = $2
    }
    close(mapping_file)
    
    print "# Fixed GTF for StringTie - Generated on", strftime("%Y-%m-%d %H:%M:%S")
    print "# Added gene features and fixed sequence references"
}

# Process GTF lines
{
    # Extract gene_id from attributes
    if(match($0, /gene_id "([^"]+)"/, arr)) {
        current_gene_id = arr[1]
        
        # Check if we have a mapping for this gene
        if(current_gene_id in gene_to_seq) {
            new_ref = gene_to_seq[current_gene_id]
        } else {
            new_ref = $1  # Keep original if no mapping
        }
        
        # If this is a transcript line, add gene feature first
        if($3 == "transcript") {
            # Only add gene if different from previous gene
            if(current_gene_id != prev_gene_id) {
                # Print gene feature
                print new_ref, $2, "gene", $4, $5, $6, $7, $8, "gene_id \"" current_gene_id "\";"
                prev_gene_id = current_gene_id
            }
        }
        
        # Print original line with fixed reference
        $1 = new_ref
        print
    } else {
        # No gene_id found, print as-is (shouldn'\''t happen in valid GTF)
        print
    }
}' "$INPUT_GTF" > "$OUTPUT_GTF"

echo -e "${GREEN}✓ Created fixed GTF: $OUTPUT_GTF${NC}"
echo ""

echo -e "${YELLOW}Step 3: Validating output...${NC}"

# Count features
echo "Feature counts in original GTF:"
cut -f3 "$INPUT_GTF" | sort | uniq -c | awk '{printf "  %-15s: %d\n", $2, $1}'
echo ""

echo "Feature counts in fixed GTF (excluding comments):"
grep -v "^#" "$OUTPUT_GTF" | cut -f3 | sort | uniq -c | awk '{printf "  %-15s: %d\n", $2, $1}'
echo ""

# Count unique sequences
echo "Reference sequences in fixed GTF:"
grep -v "^#" "$OUTPUT_GTF" | cut -f1 | sort -u | head -20
echo ""

echo -e "${GREEN}=== DONE ===${NC}"
echo ""
echo -e "${YELLOW}Next steps:${NC}"
echo "1. Review the fixed GTF: $OUTPUT_GTF"
echo "2. Update your pipeline to use: $OUTPUT_GTF"
echo "3. Ensure FASTA file matches: $INPUT_FASTA"
echo "4. Rebuild HISAT2 index with the fixed GTF"
echo ""
echo -e "${YELLOW}Sample usage in pipeline:${NC}"
echo "  GTF=\"0_INPUT_FASTAs/$OUTPUT_GTF\""
echo "  FASTA=\"0_INPUT_FASTAs/$INPUT_FASTA\""
echo ""
