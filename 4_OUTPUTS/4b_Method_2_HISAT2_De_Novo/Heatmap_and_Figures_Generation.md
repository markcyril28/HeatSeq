# Complete Pipeline Description: StringTie to R Visualization

## Overview

This document describes the complete bioinformatics pipeline for gene expression analysis from RNA-seq data, starting from StringTie quantification through visualization in R. The pipeline is designed for **Method 2: HISAT2 De Novo** approach but can be adapted for other methods.

---

## Pipeline Architecture

```
StringTie Quantification → Matrix Generation → Normalization → Visualization
        (Bash)                 (Python/Bash)        (R)            (R)
```

---

## Phase 1: StringTie Gene Expression Quantification

### Input
- **Aligned BAM files** from HISAT2 alignment (in `4_HISAT2_WD/`)
- **Reference FASTA files** for gene groups (in `../../0_INPUT_FASTAs/`)
- **Sample IDs** (SRR accessions) from RNA-seq projects

### Process
StringTie is run on each aligned BAM file to:
1. Assemble transcripts de novo (without reference genome)
2. Quantify gene expression for each sample
3. Generate abundance files with expression metrics

### Output Files
**Location:** `5_stringtie_WD/a_Method_2_RAW_RESULTs/{Gene_Group}/{SRR_ID}/`

**File format:** `{SRR}_{Gene_Group}_gene_abundances_de_novo.tsv`

**Example:** `SRR3884597_SmelGRFs_gene_abundances_de_novo.tsv`

**TSV Structure:**
```
Gene ID   Gene Name   Length   Coverage   FPKM      TPM       Num_Reads
Gene_001  GRF1        1200     45.2       12.5      8.3       542
Gene_002  GRF2        980      32.1       9.8       6.5       315
...
```

**Key columns:**
- **Gene Name** (Column 3): Gene identifier from reference FASTA
- **Coverage** (Column 7): Average read depth across gene
- **FPKM** (Column 8): Fragments Per Kilobase Million (length-normalized)
- **TPM** (Column 9): Transcripts Per Million (library-normalized)

---

## Phase 2: Count Matrix Generation

### Script: `a_stringtie_matrix_builder_m2_v4.sh`

This Bash script consolidates individual StringTie abundance files into organized count matrices.

### Input Configuration

**Gene Groups to Process:**
```bash
Gene_Groups_Boilerplates=(
    "Best_Cell_Cycle_Associated_Control_Genes"
    "Best_Control_Genes"
    "SmelDMPs"
    "SmelGIFs"
    "SmelGRFs"
    "Selected_GRF_GIF_Genes"
)
```

**Sample Lists:**
```bash
SRR_LIST_COMBINED=(
    SRR3884685   # Radicles
    SRR3884677   # Cotyledons
    SRR3884675   # Roots
    SRR3884690   # Stems
    SRR3884689   # Leaves_vegetative
    SRR3884684   # Senescent_leaves
    SRR3884686   # Buds_0.7cm
    SRR3884687   # Opened_Buds
    SRR3884597   # Flowers
    # ... additional samples
)
```

**SRR to Organ Mapping:**
```bash
declare -A SRR_TO_ORGAN=(
    ["SRR3884675"]="Roots_1"
    ["SRR3884690"]="Stems_1"
    ["SRR3884689"]="Leaves_1"
    ["SRR3884597"]="Flowers_1"
    # ... more mappings
)
```

### Process Flow

1. **For each gene group:**
   - Load reference TSV file from `5_stringtie_WD/0_Ref_Boilerplate_TSVs/{Gene_Group}.tsv`
   - Extract list of gene names from reference

2. **For each sample (SRR):**
   - Locate abundance file: `{Gene_Group}/{SRR}/{SRR}_{Gene_Group}_gene_abundances_de_novo.tsv`
   - Extract relevant columns (gene name + count metric)

3. **For each count type (coverage, FPKM, TPM):**
   - Build matrix using `matrix_builder.py`
   - Create two versions: SRR labels and Organ labels

### Matrix Builder Script: `matrix_builder.py`

**Purpose:** Merge individual sample files into a unified expression matrix

**Algorithm:**
```python
1. Read gene names from reference file (rows)
2. For each sample file:
   - Parse gene_name → count mapping
3. For each gene in reference:
   - Look up count in each sample
   - Insert "0" if gene is missing
4. Output tab-separated matrix
```

**Output Format:**
```
GeneName    SRR3884597  SRR3884686  SRR3884690  ...
GRF1        12.5        8.3         15.2        ...
GRF2        9.8         6.5         11.4        ...
GIF1        5.2         4.1         7.8         ...
```

### Output Matrices

**Location:** `5_stringtie_WD/b_Method_2_COUNT_MATRICES/{Gene_Group}/`

**File naming pattern:**
```
{Gene_Group}_{count_type}_counts_{gene_type}_{label_type}[_from_{Master_Reference}].tsv

Examples:
- SmelGRFs_coverage_counts_geneName_SRR_from_All_Smel_Genes.tsv
- SmelGRFs_fpkm_counts_geneName_Organ_from_All_Smel_Genes.tsv
- SmelGRFs_tpm_counts_geneName_SRR.tsv
```

**Matrix types generated:**
- **Count types:** coverage, FPKM, TPM (3 types)
- **Gene types:** geneID, geneName (2 types)
- **Label types:** SRR, Organ (2 types)
- **Total:** 3 × 2 × 2 = **12 matrices per gene group**

---

## Phase 3: Heatmap Visualization

### Script: `b_make_heatmap_of_matrices_v4.R`

This R script generates comprehensive heatmaps with multiple normalization schemes.

### Libraries Used
```r
library(ComplexHeatmap)  # Advanced heatmap generation
library(circlize)        # Color palettes
library(RColorBrewer)    # Color schemes
library(dplyr)           # Data manipulation
library(tibble)          # Data frames
library(grid)            # Graphics parameters
```

### Normalization Schemes

The script applies **5 normalization methods** to the raw count matrices:

#### 1. **Raw Normalization**
```r
preprocess_for_raw(data_matrix)
```
- No transformation applied
- Direct count values
- Use case: Visualizing absolute expression levels

#### 2. **Raw Normalized** (identical to Raw)
```r
preprocess_for_raw_normalized(data_matrix)
```
- Same as raw (for comparison purposes)

#### 3. **Count-Type Normalized**
```r
preprocess_for_count_type_normalized(data_matrix, count_type)
```
**For Coverage counts:**
```r
lib_sizes <- colSums(data_matrix)  # Total reads per sample
CPM <- counts / (lib_sizes / 1e6)  # Counts Per Million
log2_CPM <- log2(CPM + 1)          # Log transformation
```

**For FPKM/TPM counts:**
```r
log2_expression <- log2(counts + 0.1)  # Direct log transformation
```

- **Purpose:** Library size normalization, log-scale for visualization
- **Use case:** Comparing expression across samples with different sequencing depths

#### 4. **Z-score Normalized**
```r
preprocess_for_zscore(data_matrix, count_type)
```
**Process:**
1. Apply count-type normalization first
2. Calculate Z-scores row-wise (per gene):
```r
z_score = (x - mean(x)) / sd(x)
```
- **Purpose:** Standardize expression levels to show relative patterns
- **Use case:** Identifying genes with consistent vs. variable expression patterns

#### 5. **CPM (Counts Per Million) Normalized**
```r
preprocess_for_cpm(data_matrix)
```
```r
CPM <- (counts / lib_sizes) × 1,000,000
log2_CPM <- log2(CPM + 1)
```
- **Purpose:** Normalize for sequencing depth differences
- **Use case:** Cross-sample comparison of expression levels

### Color Scheme

**Eggplant-themed gradient** (white → purple):
```r
eggplant_colors <- colorRamp2(
  seq(min, max, length = 8),
  c("#FFFFFF", "#F3E5F5", "#CE93D8", "#AB47BC",
    "#8E24AA", "#6A1B9A", "#4A148C", "#2F1B69")
)
```

### Heatmap Features

**Matrix orientation:**
- **Rows:** Genes (y-axis, left side)
- **Columns:** Samples/Organs (x-axis, top)

**Visual elements:**
- Hierarchical clustering disabled (keeps developmental order)
- Gene names displayed on left (if < 50 genes)
- Sample names displayed on top (45° rotation)
- White grid lines separating cells
- Legend on right side

### Output Files

**Location:** `6_Heatmap_Visualizations/{Gene_Group}/{Normalization_Type}/`

**Subdirectories:**
```
├── raw/
├── raw_normalized/
├── Count-Type_Normalized/
├── Z-score_Normalized/
└── CPM_Normalized/
```

**File types:**
1. **PNG heatmaps:** `{matrix_name}_{normalization}_heatmap.png`
2. **TSV matrices:** `{matrix_name}_{normalization}_heatmap.tsv` (same data in table format)

**Example output:**
```
SmelGRFs_coverage_counts_geneName_Organ_from_All_Smel_Genes_zscore_normalized_heatmap.png
SmelGRFs_coverage_counts_geneName_Organ_from_All_Smel_Genes_zscore_normalized_heatmap.tsv
```

**Total files per matrix:** 5 normalizations × 2 file types = **10 files**

---

## Phase 4: Coefficient of Variation (CV) Heatmaps

### Script: `c_make_heatmap_with_CV_of_matrices_v4.R`

This script generates specialized heatmaps with **coefficient of variation** annotations to identify genes with variable expression patterns.

### Additional Normalization: Z-score Scaled to [0, 10]

```r
preprocess_for_zscore_scaled_to_ten(data_matrix, count_type)
```

**Process:**
1. Apply count-type normalization
2. Calculate Z-scores
3. Scale entire matrix to [0, 10] range:
```r
scaled_value = ((z_score - min) / (max - min)) × 10
```

**Total normalization schemes:** Raw, Count-Type, Z-score, CPM, **Z-score-scaled-to-10** (5 schemes)

### Coefficient of Variation Calculation

**Formula:**
```r
CV = standard_deviation / mean
```

**Interpretation:**
- **High CV (>1.0):** Gene shows high variability across developmental stages
- **Low CV (<0.3):** Gene shows consistent expression (housekeeping-like)
- **Medium CV (0.3-1.0):** Moderate stage-specific regulation

### Heatmap Features

**Color scheme** (white → yellow → purple):
```r
heatmap_colors <- colorRamp2(
  seq(-2, 2, length = 10),
  c("#FFFFFF", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", 
    "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
)
```

**Layout:**
- Left: Gene names
- Center: Expression heatmap (Z-score normalized)
- Right: CV values (numerical annotation)
- Bottom: Developmental stages (sample names)
- Hierarchical clustering: **Enabled** (shows expression pattern similarity)

### Sorting Options

Each matrix generates **2 versions**:

#### 1. **Sorted by Organ (Developmental Order)**
- Samples arranged by developmental stage
- Location: `sorted_by_organ/`
- Purpose: Follow temporal expression patterns

#### 2. **Sorted by Expression (Intensity Order)**
- Samples sorted by mean expression (low → high)
- Location: `sorted_by_expression/`
- Purpose: Identify highest/lowest expressing stages

### Label Types

Each sorting version generates **2 label variants**:

1. **SRR labels:** `SRR3884597`, `SRR3884686`, etc.
2. **Organ labels:** `Flowers_1`, `Buds_1`, `Roots_1`, etc.

### Output Structure

**Location:** `7_Heatmap_with_Covariability_Visualizations/{Gene_Group}/{Normalization}/`

```
{Gene_Group}/
├── raw/
│   ├── sorted_by_organ/
│   │   ├── {matrix}_raw_SRR_cv_heatmap.png
│   │   ├── {matrix}_raw_SRR_cv_heatmap_CV_data.tsv
│   │   ├── {matrix}_raw_Organ_cv_heatmap.png
│   │   └── {matrix}_raw_Organ_cv_heatmap_CV_data.tsv
│   └── sorted_by_expression/
│       ├── {matrix}_raw_SRR_cv_heatmap.png
│       └── ... (same 4 files)
├── count_type_normalized/
│   ├── sorted_by_organ/
│   └── sorted_by_expression/
├── zscore/
├── cpm/
└── zscore_scaled_to_ten/
```

**Files per matrix:**
- 5 normalizations × 2 sortings × 2 labels × 2 file types = **40 files**

---

## Phase 5: Bar Graph Visualization

### Script: `d_make_BarGraph_of_matrices_v4.R`

This script generates individual bar graphs for **each gene**, showing expression across developmental stages.

### Libraries Used
```r
library(ggplot2)     # Bar graph generation
library(dplyr)       # Data manipulation
library(RColorBrewer) # Color palettes
library(grid)        # Graphics
library(scales)      # Axis formatting
library(gridExtra)   # Multi-plot layouts
```

### Bar Graph Features

**Each gene generates 4 bar graphs:**

1. **Sorted by Organ + SRR labels**
2. **Sorted by Organ + Organ labels**
3. **Sorted by Expression + SRR labels**
4. **Sorted by Expression + Organ labels**

### Color Coding

**Expression-based color gradient** (white → violet):
```r
custom_colors <- c(
  "#FFFFFF",  # White (lowest)
  "#FFEB3B",  # Yellow
  "#CDDC39",  # Lime
  "#8BC34A",  # Light green
  "#4CAF50",  # Green
  "#009688",  # Teal
  "#00BCD4",  # Cyan
  "#3F51B5",  # Blue
  "#673AB7",  # Purple
  "#4A148C"   # Violet (highest)
)
```

**Color mapping:**
- Bar color intensity = expression level
- Highest expression → deepest violet
- Lowest expression → white

### Y-axis Labels (by normalization)

```r
y_label <- switch(normalization_scheme,
  "raw" = "Raw counts",
  "raw_normalized" = "Raw counts",
  "count_type_normalized" = "Log2(CPM + 1)" or "Log2(Expression + 0.1)",
  "zscore" = "Z-score",
  "cpm" = "Log2(CPM + 1)"
)
```

### Graph Layout

**Visual elements:**
- X-axis: Developmental stages (45° rotation)
- Y-axis: Expression value (labeled by normalization)
- Title: Gene name
- Bars: Color-coded by expression intensity
- Border: Black outline on each bar
- Background: White with minimal gridlines

**Dimensions:** 8" × 6" (width × height) at 300 DPI

### Output Structure

**Location:** `8_Bar_Graph_Visualizations/{Gene_Group}/{Normalization}/`

```
{Gene_Group}/
├── raw/
│   ├── sorted_by_organ/
│   │   ├── GRF1_raw_bargraph.png
│   │   ├── GRF2_raw_bargraph.png
│   │   └── ... (one per gene)
│   └── sorted_by_expression/
│       └── ... (same genes)
├── raw_normalized/
├── count_type_normalized/
├── zscore/
└── cpm/
```

**Files per gene group:**
- N genes × 5 normalizations × 2 sortings × 2 labels = **20N files**
- Example: 10 genes → 200 bar graph files

---

## Phase 6: Wrapper Script Execution

### Script: `0_Heatmap_Wrapper.sh`

This Bash script orchestrates the entire visualization pipeline.

### Execution Sequence

```bash
#!/bin/bash

# 1. Activate conda environment (if needed)
# conda activate Heatmap_ENV

# 2. Generate count matrices (optional - usually pre-run)
# bash a_stringtie_matrix_builder_m2_v4.sh

# 3. Generate standard heatmaps (optional)
# Rscript --no-save b_make_heatmap_of_matrices_v4.R

# 4. Generate CV heatmaps
Rscript --no-save c_make_heatmap_with_CV_of_matrices_v4.R

# 5. Generate bar graphs
Rscript --no-save d_make_BarGraph_of_matrices_v4.R
```

### Usage

```bash
cd 4b_Method_2_HISAT2_De_Novo
bash 0_Heatmap_Wrapper.sh
```

**Runtime:** Depends on number of gene groups and samples (typically 5-30 minutes)

---

## Summary Statistics

### Files Generated Per Gene Group

| Phase | Script | Files Generated | File Types |
|-------|--------|----------------|------------|
| Matrix Building | `a_stringtie_matrix_builder_m2_v4.sh` | 12 | TSV matrices |
| Standard Heatmaps | `b_make_heatmap_of_matrices_v4.R` | 120 | PNG + TSV (10 per matrix) |
| CV Heatmaps | `c_make_heatmap_with_CV_of_matrices_v4.R` | 480 | PNG + TSV (40 per matrix) |
| Bar Graphs | `d_make_BarGraph_of_matrices_v4.R` | 20 × N | PNG (N = gene count) |

**Total for 1 gene group with 10 genes:** ~850 files

---

## Key Configuration Parameters

### Gene Groups
```bash
FASTA_GROUPS=(
    "Best_Cell_Cycle_Associated_Control_Genes"
    "Best_Control_Genes"
    "SmelDMPs"
    "SmelGIFs"
    "SmelGRFs"
    "Selected_GRF_GIF_Genes"
)
```

### Count Metrics
```bash
COUNT_TYPES=("coverage" "fpkm" "tpm")
```

### Normalization Methods
```r
NORMALIZATION_SCHEMES <- c(
  "raw",
  "raw_normalized",
  "count_type_normalized",
  "zscore",
  "cpm",
  "zscore_scaled_to_ten"  # CV heatmaps only
)
```

### Samples Configuration
- **SRR IDs:** RNA-seq accession numbers
- **Organ labels:** Mapped developmental stages
- **Projects:** PRJNA328564, SAMN28540077, SAMN28540068

---

## Data Flow Diagram

```
StringTie Abundance Files
         ↓
[a_stringtie_matrix_builder_m2_v4.sh + matrix_builder.py]
         ↓
Count Matrices (12 per gene group)
         ↓
    ┌────┴────┬────────────┐
    ↓         ↓            ↓
[Standard] [CV Heatmaps] [Bar Graphs]
Heatmaps    (Clustered)   (Individual)
    ↓         ↓            ↓
  PNG/TSV   PNG/TSV      PNG files
```

---

## Output Directory Structure

```
4b_Method_2_HISAT2_De_Novo/
├── 5_stringtie_WD/
│   ├── 0_Ref_Boilerplate_TSVs/           # Reference gene lists
│   ├── a_Method_2_RAW_RESULTs/           # StringTie abundance files
│   │   └── {Gene_Group}/{SRR}/{SRR}_{Gene_Group}_gene_abundances_de_novo.tsv
│   └── b_Method_2_COUNT_MATRICES/        # Generated matrices
│       └── {Gene_Group}/{Gene_Group}_{count_type}_counts_{gene_type}_{label_type}.tsv
├── 6_Heatmap_Visualizations/             # Standard heatmaps
│   └── {Gene_Group}/{Normalization}/{matrix}_heatmap.png
├── 7_Heatmap_with_Covariability_Visualizations/  # CV heatmaps
│   └── {Gene_Group}/{Normalization}/{Sorting}/{matrix}_cv_heatmap.png
└── 8_Bar_Graph_Visualizations/           # Individual gene bar graphs
    └── {Gene_Group}/{Normalization}/{Sorting}/{Gene}_bargraph.png
```

---

## Quality Control & Validation

### Matrix Builder Validation
- Checks for duplicate gene names (makes unique if found)
- Replaces missing values with "0"
- Logs processing statistics (samples found, matrices generated)

### R Script Validation
- Filters out matrices with < 2 genes
- Handles NA and infinite values (sets to 0)
- Removes matrices with no variance (all identical values)
- Reports success/failure rates in summary

### Log Files
- Matrix builder: `5_stringtie_WD/b_Method_2_COUNT_MATRICES/logs/matrix_builder_{timestamp}.log`
- R scripts: Console output with summary statistics

---

## Troubleshooting

### Common Issues

1. **Missing input files**
   - Check StringTie abundance files exist in `5_stringtie_WD/a_Method_2_RAW_RESULTs/`
   - Verify reference TSV files in `5_stringtie_WD/0_Ref_Boilerplate_TSVs/`

2. **R package errors**
   - Install missing packages: `install_heatmap_libraries.R`
   - Activate correct conda environment: `conda activate Heatmap_ENV`

3. **Memory issues**
   - Reduce number of gene groups processed simultaneously
   - Filter high-variance genes (top 1000) for large datasets

4. **Visualization errors**
   - Ensure at least 2 genes in each group
   - Check for zero-variance data (all samples identical)
   - Verify color palette compatibility

---

## Best Practices

### For Matrix Generation
1. Use master reference when querying multiple gene groups (`QUERY_AGAINST_MASTER_REFERENCE=TRUE`)
2. Maintain consistent SRR-to-organ mapping across analyses
3. Keep log files for reproducibility

### For Visualization
1. Use Z-score normalization for pattern detection
2. Use count-type normalization for absolute expression comparison
3. Generate both organ-sorted and expression-sorted versions
4. Review CV values to identify housekeeping vs. regulated genes

### For Publication
1. Use high-resolution outputs (300 DPI PNG files)
2. Include both SRR and organ labels for supplementary materials
3. Save TSV matrices alongside visualizations for data sharing
4. Document normalization methods in figure legends

---

## Citations & Dependencies

### Bioinformatics Tools
- **HISAT2:** Graph-based alignment
- **StringTie:** Transcript assembly and quantification
- **Python 3:** Matrix building script

### R Packages
- **ComplexHeatmap:** Advanced heatmap visualization
- **ggplot2:** Grammar of graphics plotting
- **dplyr/tibble:** Data manipulation
- **circlize:** Color palettes and circular layouts
- **RColorBrewer:** Pre-defined color schemes

---

## Version Information

- **Matrix Builder:** v4
- **Heatmap Scripts:** v4
- **Bar Graph Script:** v4
- **Pipeline Method:** Method 2 (HISAT2 De Novo)

**Last Updated:** 2025

---

## Contact & Support

For questions or issues with this pipeline, refer to:
- Main pipeline documentation: `../../CMSC244_Reporting_Alignment_Methods_in_GEA.md`
- Methods comparison: `../../GEA_Methods_Comparison_Detailed.md`
- Setup instructions: `setup_mamba_heatmap.sh`

---

**End of Pipeline Description**
