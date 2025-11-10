#!/usr/bin/env Rscript

# ===============================================
# HEATMAP GENERATION FOR METHOD 5 (BOWTIE2/RSEM)
# ===============================================
# Generates heatmaps from DESeq2/tximport-normalized RSEM data
# Multiple normalization schemes and output formats

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

source("modules_method_5_bowtie2/2_utility_functions.R")

# ===============================================
# CONFIGURATION
# ===============================================

LEGEND_POSITION <- "bottom"
BASE_DIR <- getwd()
MATRICES_DIR <- file.path("4_matrices")
# Output directly to consolidated directory
HEATMAP_OUT_DIR <- "ALL_HEATMAPS_CONSOLIDATED/Basic_Heatmaps"
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

# RSEM provides expected_count, TPM, and FPKM
# After DESeq2/tximport processing we get normalized counts
COUNT_TYPES <- c("expected_count")
GENE_TYPES <- c("geneID", "geneName")
LABEL_TYPES <- c("SRR", "Organ")

# Processing levels (gene-level and isoform-level)
PROCESSING_LEVELS <- c("gene_level", "isoform_level")

# Normalization schemes
NORM_SCHEMES <- c(
  "count_type_normalized",  # Log2 transformation
  "zscore",                 # Z-score normalization
  "zscore_scaled_to_ten"    # Z-score scaled to [0,10]
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


dir.create(HEATMAP_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ===============================================
# MAIN PROCESSING
# ===============================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("BASIC HEATMAP GENERATION - METHOD 5 (BOWTIE2/RSEM)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Configuration:\n")
cat("  • Overwrite existing files:", OVERWRITE_EXISTING, "\n")
cat("  • Gene groups:", paste(GENE_GROUPS, collapse = ", "), "\n\n")

total_heatmaps <- 0
successful_heatmaps <- 0
skipped_heatmaps <- 0

for (gene_group in GENE_GROUPS) {
  cat("Processing gene group:", gene_group, "\n")

  gene_group_output_dir <- file.path(HEATMAP_OUT_DIR, gene_group)
  if (dir.exists(gene_group_output_dir)) unlink(gene_group_output_dir, recursive = TRUE)
  dir.create(gene_group_output_dir, recursive = TRUE, showWarnings = FALSE)

  for (processing_level in PROCESSING_LEVELS) {
    for (count_type in COUNT_TYPES) {
      for (gene_type in GENE_TYPES) {
        for (label_type in LABEL_TYPES) {

          # Input file path from DESeq2 normalized counts (now includes processing_level)
          input_file <- file.path(MATRICES_DIR, MASTER_REFERENCE, processing_level, gene_group,
                                 paste0(gene_group, "_", count_type, "_", gene_type,
                                       "_from_", MASTER_REFERENCE, "_", processing_level, ".tsv"))

          if (!file.exists(input_file)) {
            cat("  Skipping (file not found):", processing_level, "|", basename(input_file), "\n")
            next
          }

          # Read matrix
          raw_data_matrix <- read_count_matrix(input_file)
          if (is.null(raw_data_matrix)) {
            cat("  Skipping (failed to read):", processing_level, "|", basename(input_file), "\n")
            next
          }
          if (nrow(raw_data_matrix) < 2) {
            cat("  Skipping (need ≥2 genes for heatmap, found", nrow(raw_data_matrix), "):", processing_level, "|", basename(input_file), "\n")
            next
          }

          cat("  Processing:", processing_level, "|", count_type, "|", gene_type, "|", label_type, "- Found", nrow(raw_data_matrix), "genes\n")

          # Generate heatmaps for each normalization scheme
          for (norm_scheme in NORM_SCHEMES) {

            # Apply normalization
            data_normalized <- apply_normalization(raw_data_matrix, norm_scheme, count_type)

            if (is.null(data_normalized)) {
              cat("    Failed normalization:", norm_scheme, "\n")
              next
            }

            # Determine normalization display name
            norm_display <- switch(norm_scheme,
              "count_type_normalized" = "Count-Type_Normalized",
              "zscore" = "Z-score_Normalized",
              "zscore_scaled_to_ten" = "Z-score_Scaled_to_Ten",
              norm_scheme
            )

            # Generate all versions (transposed/sorted combinations)
            sorting_options <- list(
              list(sort = FALSE, sort_name = "Sorted_by_Organ"),
              list(sort = TRUE, sort_name = "Sorted_by_Expression")
            )

            orientation_options <- list(
              list(transpose = FALSE, orient_name = "Original_RowGene_ColumnOrgan"),
              list(transpose = TRUE, orient_name = "Transposed_RowOrgan_ColumnGene")
            )

            for (orient in orientation_options) {
              for (sorting in sorting_options) {

              # Create output directory structure (include processing_level)
              version_dir <- file.path(gene_group_output_dir, processing_level, count_type, gene_type,
                                      norm_scheme, orient$orient_name, sorting$sort_name)
              dir.create(version_dir, recursive = TRUE, showWarnings = FALSE)

              # Build title and output filename (include processing_level)
              title_base <- paste0(gene_group, "_", count_type, "_", gene_type, "_",
                                  label_type, "_from_", MASTER_REFERENCE, "_", processing_level, "_", norm_scheme)

              output_path <- file.path(version_dir, paste0(title_base, ".png"))
              
              total_heatmaps <- total_heatmaps + 1

                # Check if file exists and skip if not overwriting
                if (!OVERWRITE_EXISTING && file.exists(output_path)) {
                  cat("      Skipping (already exists):", basename(output_path), "\n")
                  skipped_heatmaps <- skipped_heatmaps + 1
                  next
                }

                # Generate heatmap
                success <- generate_heatmap_violet(
                  data_matrix = data_normalized,
                  output_path = output_path,
                  title = title_base,
                  count_type = count_type,
                  label_type = label_type,
                  normalization_type = norm_display,
                  transpose = orient$transpose,
                  sort_by_expression = sorting$sort
                )

                if (success) {
                  successful_heatmaps <- successful_heatmaps + 1
                } else {
                  cat("      Failed:", basename(output_path), "\n")
                }
              }
            }
          }
        }
      }
    }
  }
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY:", successful_heatmaps, "/", total_heatmaps, "heatmaps generated")
if (skipped_heatmaps > 0) {
  cat(" (", skipped_heatmaps, " skipped)\n", sep = "")
} else {
  cat("\n")
}
cat(paste(rep("=", 60), collapse = ""), "\n\n")

if (successful_heatmaps > 0) {
  cat("Output directory:", HEATMAP_OUT_DIR, "\n")
  cat("Heatmaps organized by:\n")
  cat("  - Gene group\n")
  cat("  - Count type (expected_count)\n")
  cat("  - Gene identifier type (geneID/geneName)\n")
  cat("  - Normalization scheme\n")
  cat("  - Orientation (original/transposed)\n")
  cat("  - Sorting (by organ/by expression)\n")
}
