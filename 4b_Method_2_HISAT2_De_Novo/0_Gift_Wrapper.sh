#!/bin/bash 

: << 'DESCRIPTIONS'
Run this script in sequence:
1. bash a_stringtie_method_2.1.sh
2. Rscript b_heatmap_DeSeq2_v2.R
DESCRIPTIONS

#bash a_stringtie_method_2.2.sh

Rscript --verbose b_heatmap_METHOD2_RESULTS_matrices.R
