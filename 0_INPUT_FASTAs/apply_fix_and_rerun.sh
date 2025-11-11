#!/bin/bash
################################################################################
# Quick Fix Application Script
# Run this to apply the fix and rerun Method 1 pipeline
################################################################################

set -euo pipefail

echo "=========================================="
echo "Method 1 HISAT2 Reference-Guided - FIX"
echo "=========================================="
echo ""

# Configuration
FIXED_GTF="0_INPUT_FASTAs/All_Smel_Genes_Full_Name_FIXED.gtf"
FASTA="0_INPUT_FASTAs/All_Smel_Genes.fasta"
METHOD1_DIR="4_OUTPUTS/4a_Method_1_HISAT2_Ref_Guided"

# Verify fixed GTF exists
if [[ ! -f "$FIXED_GTF" ]]; then
    echo "ERROR: Fixed GTF not found: $FIXED_GTF"
    echo "Run fix_gtf_for_stringtie.sh first!"
    exit 1
fi

echo "✓ Fixed GTF found: $FIXED_GTF"
echo ""

# Prompt user
read -p "This will DELETE old Method 1 results and rebuild. Continue? (y/N): " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 0
fi

echo ""
echo "Step 1: Removing old HISAT2 index (built with wrong GTF)..."
rm -rf "$METHOD1_DIR/1_HISAT2_Ref_Guided_Index/All_Smel_Genes"*
echo "✓ Removed old index"
echo ""

echo "Step 2: Removing old StringTie results..."
rm -rf "$METHOD1_DIR/5_stringtie_WD/a_Method_1_RAW_RESULTs/All_Smel_Genes_copy"
echo "✓ Removed old quantifications"
echo ""

echo "Step 3: Running Method 1 with FIXED GTF..."
echo "Command: ./b_GEA_script_v12.sh --method 1 --fasta \"$FASTA\" --gtf \"$FIXED_GTF\""
echo ""
read -p "Execute now? (y/N): " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    ./b_GEA_script_v12.sh --method 1 --fasta "$FASTA" --gtf "$FIXED_GTF"
else
    echo ""
    echo "Skipped execution. Run manually when ready:"
    echo "  ./b_GEA_script_v12.sh --method 1 --fasta \"$FASTA\" --gtf \"$FIXED_GTF\""
fi

echo ""
echo "=========================================="
echo "DONE"
echo "=========================================="
