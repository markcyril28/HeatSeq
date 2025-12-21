# ===============================================
# SHARED CONFIGURATION FOR METHOD 5 HEATMAP GENERATION
# ===============================================
# Centralized configuration to avoid duplication across scripts

# ===============================================
# CONSTANTS (Single source of truth for directories)
# ===============================================

MATRICES_DIR <- file.path("6_matrices")
CONSOLIDATED_BASE_DIR <- "7_Figures_Outputs"
#MASTER_REFERENCE <- "All_Smel_Genes"
MASTER_REFERENCE <- "Eggplant_V4.1_transcripts.function"

# Output subdirectories - must match bash configuration
HEATMAP_SUBDIR <- "I_Basic_Heatmap"
CV_HEATMAP_SUBDIR <- "II_Heatmap_with_CV"
BAR_GRAPH_SUBDIR <- "III_Bar_Graphs"
WGCNA_SUBDIR <- "IV_Coexpression_WGCNA"
DEA_SUBDIR <- "V_Differential_Expression"
GSEA_SUBDIR <- "VI_Gene_Set_Enrichment"
DIM_REDUCTION_SUBDIR <- "VII_Dimensionality_Reduction"
CORRELATION_SUBDIR <- "VIII_Sample_Correlation"
TISSUE_SPEC_SUBDIR <- "IX_Tissue_Specificity"

COUNT_TYPES <- c("expected_count")
GENE_TYPES <- c(
  #"geneID", 
  "geneName"
)

LABEL_TYPES <- c(
  #"SRR", 
  "Organ"
)

PROCESSING_LEVELS <- c(
  #"gene_level", 
  "isoform_level"
)

NORM_SCHEMES <- c(
  #"count_type_normalized",
  #"zscore",
  "zscore_scaled_to_ten"
)

# Sample labels mapping (SRR ID -> Organ name)
SAMPLE_LABELS <- c(
  "SRR3884675" = "Root",
  "SRR20722229" = "Root_2",
  "SRR31755282" = "Root_3",
  "SRR3884690" = "Stem",
  "SRR20722227" = "Stem_2",
  "SRR20722384" = "Stem_3",
  "SRR20722233" = "Leaf_buds_2",
  "SRR20722296" = "Leaf_buds_3",
  "SRR20722383" = "Young_Leaves",
  "SRR3884689" = "Mature_Leaves",
  "SRR20722230" = "Mature_Leaves_2",
  "SRR20722386" = "Mature_Leaves_3",
  #"SRR3884678" = "Fruits_peduncle",
  #"SRR20722228" = "Sepals_2",
  #"SRR20722385" = "Sepals_3",
  "SRR3884686" = "Flower_Buds",
  "SRR4243802" = "Flower_Buds_2",
  "SRR20722297" = "Flower_Buds_3",
  "SRR3884687" = "Opened_Buds",
  "SRR3884597" = "Flowers",
  "SRR20722234" = "Flowers_2",
  "SRR3884608" = "Young_Fruits",
  "SRR20722226" = "Young_Fruits_2",
  "SRR3884631" = "Mature_Fruits",
  "SRR20722232" = "Mature_Fruits_2",
  "SRR20722387" = "Mature_Fruits_3",
  "SRR3884685" = "Radicles",
  "SRR3884677" = "Cotyledons",
  "SRR3884679" = "Pistils"
)

# Sample IDs
SAMPLE_IDS <- names(SAMPLE_LABELS)

# ===============================================
# HELPER FUNCTIONS
# ===============================================

# Read configuration from wrapper script files
read_config_file <- function(file_path, default_value, is_boolean = FALSE) {
  if (!file.exists(file_path)) {
    return(default_value)
  }
  
  value <- trimws(readLines(file_path, warn = FALSE))
  
  if (is_boolean) {
    return(tolower(value[1]) == "true")
  } else {
    return(value[nzchar(value)])  # Remove empty lines
  }
}

# Load gene groups and overwrite settings
load_runtime_config <- function() {
  gene_groups <- read_config_file(
    "B_M5_modules/.gene_groups_temp.txt",
    default_value = c("Selected_GRF_GIF_Genes_vAll_GIFs", 
                      "Selected_GRF_GIF_Genes_vTwo_GIFs")
  )
  
  master_reference <- read_config_file(
    "B_M5_modules/.master_reference_temp.txt",
    default_value = "Eggplant_V4.1_transcripts.function"
  )
  # Use first element if multiple lines
  if (length(master_reference) > 1) {
    master_reference <- master_reference[1]
  }
  
  overwrite <- read_config_file(
    "B_M5_modules/.overwrite_temp.txt",
    default_value = TRUE,
    is_boolean = TRUE
  )
  
  list(
    gene_groups = gene_groups,
    master_reference = master_reference,
    overwrite_existing = overwrite
  )
}

# Print separator line
print_separator <- function(char = "=", width = 60) {
  cat("\n", paste(rep(char, width), collapse = ""), "\n")
}

# Print configuration summary
print_config_summary <- function(title, config) {
  print_separator()
  cat(title, "\n")
  print_separator()
  cat("\nConfiguration:\n")
  cat("  • Overwrite existing files:", config$overwrite_existing, "\n")
  cat("  • Gene groups:", paste(config$gene_groups, collapse = ", "), "\n\n")
}

# Print summary statistics
print_summary <- function(successful, total, skipped = NULL) {
  print_separator()
  cat("SUMMARY:", successful, "/", total, "items generated")
  if (!is.null(skipped) && skipped > 0) {
    cat(" (", skipped, " skipped)\n", sep = "")
  } else {
    cat("\n")
  }
  print_separator()
}

# Build input file path
build_input_path <- function(gene_group, processing_level, count_type, gene_type) {
  # If gene_group is same as MASTER_REFERENCE, use root level (no subdirectory)
  if (gene_group == MASTER_REFERENCE) {
    file.path(
      MATRICES_DIR, MASTER_REFERENCE, processing_level,
      paste0(MASTER_REFERENCE, "_", count_type, "_", gene_type,
             "_from_", MASTER_REFERENCE, "_", processing_level, ".tsv")
    )
  } else {
    file.path(
      MATRICES_DIR, MASTER_REFERENCE, processing_level, gene_group,
      paste0(gene_group, "_", count_type, "_", gene_type,
             "_from_", MASTER_REFERENCE, "_", processing_level, ".tsv")
    )
  }
}

# Build title base
build_title_base <- function(gene_group, count_type, gene_type, label_type, processing_level, norm_scheme) {
  paste0(gene_group, "_", count_type, "_", gene_type, "_",
         label_type, "_from_", MASTER_REFERENCE, "_", processing_level, "_", norm_scheme)
}

# Validate and read matrix
validate_and_read_matrix <- function(input_file, min_rows = 2) {
  if (!file.exists(input_file)) {
    return(list(success = FALSE, reason = "file not found"))
  }
  
  matrix_data <- read_count_matrix(input_file)
  
  if (is.null(matrix_data)) {
    return(list(success = FALSE, reason = "failed to read"))
  }
  
  if (nrow(matrix_data) < min_rows) {
    return(list(success = FALSE, reason = paste0("need ≥", min_rows, " rows, found ", nrow(matrix_data))))
  }
  
  list(success = TRUE, data = matrix_data, n_genes = nrow(matrix_data))
}

# Check if should skip (file exists and not overwriting)
should_skip_existing <- function(output_path, overwrite) {
  !overwrite && file.exists(output_path)
}

# Create output directory safely
ensure_output_dir <- function(dir_path, clean = FALSE) {
  if (clean && dir.exists(dir_path)) {
    unlink(dir_path, recursive = TRUE)
  }
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# Print final output summary (DRY helper)
print_output_summary <- function(output_dir, organization_items = NULL) {
  if (!is.null(output_dir)) {
    cat("Output directory:", output_dir, "\n")
  }
  if (!is.null(organization_items) && length(organization_items) > 0) {
    cat("Outputs organized by:\n")
    for (item in organization_items) {
      cat("  -", item, "\n")
    }
  }
}
