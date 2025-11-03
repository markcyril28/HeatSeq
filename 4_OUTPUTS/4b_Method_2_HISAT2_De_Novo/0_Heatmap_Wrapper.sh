#!/bin/bash

#===============================================================================
# HISAT2 De Novo Pipeline - Heatmap Wrapper Script
#===============================================================================
# Switches to Heatmap_ENV conda environment and runs StringTie + heatmap generation
# Usage: bash 0_Heatmap_Wrapper.sh
#===============================================================================

# Activate the Heatmap_ENV conda environment
#conda activate Heatmap_ENV

#bash a_stringtie_matrix_builder_m2_v4.sh

#Rscript --no-save b_make_heatmap_of_matrices_v4.R

Rscript --no-save c_make_heatmap_with_CV_of_matrices_v4.R

Rscript --no-save d_make_BarGraph_of_matrices_v4.R
