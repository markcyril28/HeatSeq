# Gene Expression Analysis (GEA) Pipeline - Overview

This pipeline processes RNA-seq data for gene expression analysis in eggplant, supporting multiple quantification methods. Below is a high-level description of the workflow and main steps:

---

## 1. Input Data Preparation

- **Raw RNA-seq Data**: FASTQ files for various samples (SRR IDs) are stored in `1_RAW_SRR/`.
- **Reference Files**: Gene FASTA files and annotation files are in `0_INPUT_FASTAs/`.

## 2. Environment Setup

- **Conda/Mamba Environment**: Scripts like `setup_GEA_env.sh` and `setup_heatmap.sh` install required bioinformatics tools (HISAT2, StringTie, Salmon, Bowtie2, RSEM, FastQC, MultiQC, etc.) and R packages for downstream analysis.

## 3. Data Preprocessing

- **Download & Trim**: Raw reads are downloaded and trimmed for quality using `parallel-fastq-dump` and `trim-galore`. Results are stored in `2_TRIMMED_SRR/`.
- **Quality Control**: FastQC and MultiQC generate QC reports in `3_FastQC/`.

## 4. Alignment & Quantification (Multiple Methods)

- **Method 1: HISAT2 Reference Guided**
  - Align reads to reference genome using HISAT2.
  - Assemble transcripts and quantify with StringTie.
  - Output: `4a_Method_1_HISAT2_Ref_Guided/`
- **Method 2: HISAT2 De Novo (Main)**
  - Align reads to de novo assembled transcriptome.
  - Quantify with StringTie.
  - Output: `4b_Method_2_HISAT2_De_Novo/`
- **Method 3: Trinity De Novo**
  - Assemble transcripts with Trinity, quantify with StringTie.
  - Output: `4c_Method_3_Trinity_De_Novo/`
- **Method 4: Salmon SAF Quantification**
  - Quantify expression using Salmon.
  - Output: `4d_Method_4_Salmon_Saf_Quantification/`
- **Method 5: Bowtie2 + RSEM Quantification**
  - Quantify expression using Bowtie2 and RSEM.
  - Output: `4e_Method_5_Bowtie2_Quantification/`

## 5. Matrix Generation

- **StringTie Matrix Builder**: Scripts like `a_stringtie_matrix_builder_m2_v4.sh` generate count matrices (coverage, FPKM, TPM) for each gene group and sample set. Output is in `5_stringtie_WD/b_Method_2_COUNT_MATRICES/`.

## 6. Visualization & Statistical Analysis

- **Heatmap Generation**: R scripts (`b_make_heatmap_of_matrices_v4.R`, `c_make_heatmap_with_CV_of_matrices_v4.R`) create heatmaps of gene expression and coefficient of variation.
- **Bar Graph Generation**: `d_make_bargraph_of_matrices_v4.R` generates bar graphs for individual gene expression across developmental stages, with multiple normalization schemes. Output is in `8_Bar_Graph_Visualizations/`.

## 7. Output & Organization

- All results are organized by method, gene group, normalization scheme, and visualization type.
- Logs and intermediate files are stored in `logs/` and method-specific subfolders.

---

## Pipeline Execution

- The main pipeline is run via `a_GEA_run.sh`, which calls `b_GEA_script_v10_Method_2_HISAT2_De_Novo.sh` and orchestrates all steps above.
- Wrapper scripts (e.g., `0_Heatmap_Wrapper.sh`) activate environments and run visualization scripts in sequence.

---

## Summary Diagram

```
Raw Data → QC → Alignment/Assembly → Quantification → Matrix Generation → Visualization
         (FastQC)   (HISAT2/Trinity/Salmon/Bowtie2)   (StringTie/RSEM)   (R scripts)
```

---

## References

- See individual method folders and scripts for detailed configuration and options.
- All steps are logged for reproducibility and troubleshooting.

---

This pipeline is modular and supports extension to new gene groups, normalization schemes, and visualization types as needed.
