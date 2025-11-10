#!/usr/bin/env Rscript

# ===============================================
# BAR GRAPH GENERATION FOR METHOD 5 (BOWTIE2/RSEM)
# ===============================================
# Generates individual bar graphs for each gene
# from DESeq2/tximport-normalized RSEM data

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("modules_method_5_bowtie2/2_utility_functions.R")

# ===============================================
# CONFIGURATION
# ===============================================

MATRICES_DIR <- file.path("4_matrices")
# Output directly to consolidated directory
BAR_GRAPH_OUT_DIR <- "ALL_HEATMAPS_CONSOLIDATED/Bar_Graphs"
MASTER_REFERENCE <- "All_Smel_Genes"

# Read gene groups from wrapper script
GENE_GROUPS_FILE <- "modules_method_5_bowtie2/.gene_groups_temp.txt"
if (file.exists(GENE_GROUPS_FILE)) {
  GENE_GROUPS <- readLines(GENE_GROUPS_FILE)
  GENE_GROUPS <- GENE_GROUPS[nzchar(GENE_GROUPS)]  # Remove empty lines
  cat("Gene groups from wrapper:", paste(GENE_GROUPS, collapse = ", "), "\n")
} else {
  # Fallback to default gene groups if file not found
  GENE_GROUPS <- c(
    "Selected_GRF_GIF_Genes_vAll_GIFs",
    "Selected_GRF_GIF_Genes_vTwo_GIFs"
  )
  cat("WARNING: Gene groups file not found, using defaults\n")
}

# Read overwrite setting from wrapper script
OVERWRITE_FILE <- "modules_method_5_bowtie2/.overwrite_temp.txt"
if (file.exists(OVERWRITE_FILE)) {
  overwrite_setting <- tolower(trimws(readLines(OVERWRITE_FILE, n = 1)))
  OVERWRITE_EXISTING <- overwrite_setting == "true"
  cat("Overwrite existing files:", OVERWRITE_EXISTING, "\n")
} else {
  OVERWRITE_EXISTING <- TRUE  # Default to overwrite
  cat("WARNING: Overwrite setting not found, defaulting to TRUE\n")
}

COUNT_TYPES <- c("expected_count")
GENE_TYPES <- c("geneID", "geneName")
LABEL_TYPES <- c("SRR", "Organ")

# Processing levels (gene-level and isoform-level)
PROCESSING_LEVELS <- c("gene_level", "isoform_level")

# Normalization schemes
NORM_SCHEMES <- c(
  "count_type_normalized",
  "zscore",
  "zscore_scaled_to_ten"
)

SAMPLE_LABELS <- c(
    # Roots
    "SRR3884675" = "Roots",      # PRJNA328564
    #"SRR20722229" = "Roots_2",     # SAMN28540077
    #"SRR31755282" = "Roots_3",     # SAMN28540068

    # Stems
    "SRR3884690" = "Stems",      # PRJNA328564
    #"SRR20722227" = "Stems_2",     # SAMN28540077
    #"SRR20722384" = "Stems_3",     # SAMN28540068

    # Leaves
    "SRR3884689" = "Leaves",     # PRJNA328564
    #"SRR20722230" = "Leaves_2",    # SAMN28540077
    #"SRR20722386" = "Leaves_3",    # SAMN28540068
    "SRR3884684" = "Senescent_leaves", # PRJNA328564

    # Buds
    "SRR3884686" = "Buds",       # PRJNA328564
    #"SRR21010466" = "Buds_2",      # SAMN28540077
    #"SRR20722297" = "Buds_3",      # SAMN28540068

    # Opened Buds
    "SRR3884687" = "Opened_Buds", # PRJNA328564

    # Flowers
    "SRR3884597" = "Flowers",    # PRJNA328564
    #"SRR20722234" = "Flowers_2",   # SAMN28540077
    #"SRR23909863" = "Flowers_3",   # SAMN28540068

    # Fruits
    "SRR3884631" = "Fruits",     # PRJNA328564
    #"SRR2072232" = "Fruits_2",     # SAMN28540077
    #"SRR20722387" = "Fruits_3",    # SAMN28540068
    "SRR3884608" = "Fruits_1cm",   # PRJNA328564
    "SRR3884620" = "Fruits_Stage_1", # PRJNA328564
    "SRR3884642" = "Fruits_Skin_Stage_2", # PRJNA328564
    "SRR3884653" = "Fruits_Flesh_Stage_2", # PRJNA328564
    "SRR3884664" = "Fruits_Calyx_Stage_2", # PRJNA328564
    "SRR3884680" = "Fruits_Skin_Stage_3", # PRJNA328564
    "SRR3884681" = "Fruits_Flesh_Stage_3", # PRJNA328564
    "SRR3884678" = "Fruits_peduncle", # PRJNA328564

    # Other organs
    "SRR3884685" = "Radicles",     # PRJNA328564
    "SRR3884677" = "Cotyledons",   # PRJNA328564
    "SRR3884679" = "Pistils"       # PRJNA328564
)


dir.create(BAR_GRAPH_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ===============================================
# MAIN PROCESSING
# ===============================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("BAR GRAPH GENERATION - METHOD 5 (BOWTIE2/RSEM)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Configuration:\n")
cat("  • Overwrite existing files:", OVERWRITE_EXISTING, "\n")
cat("  • Gene groups:", paste(GENE_GROUPS, collapse = ", "), "\n\n")

total_graphs <- 0
successful_graphs <- 0

for (gene_group in GENE_GROUPS) {
  cat("Processing gene group:", gene_group, "\n")

  gene_group_output_dir <- file.path(BAR_GRAPH_OUT_DIR, gene_group)
  if (dir.exists(gene_group_output_dir)) unlink(gene_group_output_dir, recursive = TRUE)
  dir.create(gene_group_output_dir, recursive = TRUE, showWarnings = FALSE)

  for (processing_level in PROCESSING_LEVELS) {
    for (count_type in COUNT_TYPES) {
      for (gene_type in GENE_TYPES) {
        for (label_type in LABEL_TYPES) {

          # Input file path (now includes processing_level)
          input_file <- file.path(MATRICES_DIR, MASTER_REFERENCE, processing_level, gene_group,
                                 paste0(gene_group, "_", count_type, "_", gene_type,
                                       "_from_", MASTER_REFERENCE, "_", processing_level, ".tsv"))

          if (!file.exists(input_file)) {
            cat("  Skipping (file not found):", processing_level, "|", basename(input_file), "\n")
            next
          }

          # Read raw data matrix
          raw_data_matrix <- read_count_matrix(input_file)
          if (is.null(raw_data_matrix) || nrow(raw_data_matrix) < 1) {
            cat("  Skipping (insufficient data):", processing_level, "|", basename(input_file), "\n")
            next
          }

          cat("  Processing:", processing_level, "|", count_type, "|", gene_type, "|", label_type, "\n")

          # Generate bar graphs for each normalization scheme
          for (norm_scheme in NORM_SCHEMES) {

            # Apply normalization
            data_normalized <- apply_normalization(raw_data_matrix, norm_scheme, count_type)

            if (is.null(data_normalized)) {
              cat("    Failed normalization:", norm_scheme, "\n")
              next
            }

            # Generate versions with different sorting
            sorting_options <- list(
              list(sort = FALSE, sort_name = "Sorted_by_Organ"),
              list(sort = TRUE, sort_name = "Sorted_by_Expression")
            )

            for (sorting in sorting_options) {

              # Create output directory structure (include processing_level)
              version_dir <- file.path(gene_group_output_dir, processing_level, count_type, gene_type,
                                      norm_scheme, sorting$sort_name)
              dir.create(version_dir, recursive = TRUE, showWarnings = FALSE)

              # Generate bar graphs
              n_graphs <- generate_gene_bargraphs(
                data_matrix = data_normalized,
                output_dir = version_dir,
                title_prefix = gene_group,
                count_type = count_type,
                label_type = label_type,
                normalization_scheme = norm_scheme,
                sort_by_expression = sorting$sort,
                overwrite = OVERWRITE_EXISTING
              )

              total_graphs <- total_graphs + nrow(data_normalized)
              successful_graphs <- successful_graphs + n_graphs

              if (n_graphs > 0) {
                cat("    Generated", n_graphs, "bar graphs for", processing_level, "|", norm_scheme, "|", sorting$sort_name, "\n")
              }
            }
          }
        }
      }
    }
  }
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY:", successful_graphs, "/", total_graphs, "bar graphs generated\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

if (successful_graphs > 0) {
  cat("Output directory:", BAR_GRAPH_OUT_DIR, "\n")
  cat("Bar graphs organized by:\n")
  cat("  - Gene group\n")
  cat("  - Count type (expected_count)\n")
  cat("  - Gene identifier type (geneID/geneName)\n")
  cat("  - Normalization scheme\n")
  cat("  - Sorting (by organ/by expression)\n")
  cat("\nEach gene has its own individual bar graph showing expression across samples.\n")
}
