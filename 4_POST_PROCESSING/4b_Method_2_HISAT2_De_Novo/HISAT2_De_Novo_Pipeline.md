# HISAT2 De Novo Pipeline - Description

This folder contains scripts and resources for the HISAT2 De Novo gene expression analysis pipeline (Method 2). Below is an overview of the workflow and main steps:

---

## Overview

This pipeline processes RNA-seq data by assembling transcripts de novo (without a reference genome) and quantifying gene expression across samples. It is designed for eggplant transcriptomics but is adaptable to other organisms.

---

## Workflow Steps

1. **Input Preparation**

   - Raw RNA-seq FASTQ files for each sample (SRR IDs) are stored in `../../1_RAW_SRR/` and trimmed reads in `../../2_TRIMMED_SRR/`.
   - Gene reference FASTA files are in `../../0_INPUT_FASTAs/`.
2. **Environment Setup**

   - Use `setup_mamba_heatmap.sh` to install required bioinformatics tools and R packages for downstream analysis.
3. **Alignment & Assembly**

   - Align trimmed reads to de novo assembled transcriptome using HISAT2.
   - Assemble transcripts and quantify expression with StringTie.
   - Scripts: `b_GEA_script_v10_Method_2_HISAT2_De_Novo.sh`, `a_stringtie_matrix_builder_m2_v4.sh`.
4. **Matrix Generation**

   - Generate count matrices (coverage, FPKM, TPM) for each gene group and sample set using `a_stringtie_matrix_builder_m2_v4.sh`.
   - Output: `5_stringtie_WD/b_Method_2_COUNT_MATRICES/`
5. **Visualization & Analysis**

   - **Heatmaps**: `b_make_heatmap_of_matrices_v4.R` and `c_make_heatmap_with_CV_of_matrices_v4.R` create heatmaps and coefficient of variation plots.
   - **Bar Graphs**: `d_make_bargraph_of_matrices_v4.R` generates bar graphs for individual gene expression across developmental stages.
   - Output: `6_Heatmap_Visualizations/`, `7_Heatmap_with_Covariability_Visualizations/`, `8_Bar_Graph_Visualizations/`
6. **Wrapper Script**

   - `0_Heatmap_Wrapper.sh` automates environment activation and sequential execution of matrix building and visualization scripts.

---

## Output Organization

- Results are organized by gene group, normalization scheme, and visualization type.
- All steps are logged for reproducibility.

---

## Execution

- Run the main pipeline via `../../a_GEA_run.sh` or directly use the wrapper script in this folder.

---

## Summary Diagram

```
Raw Data → Trim → HISAT2 Alignment → StringTie Assembly/Quantification → Matrix Generation → Visualization
```

---

For more details, see individual scripts and subfolders.
