# ===============================================
# SHARED CONFIGURATION FOR METHOD 5 HEATMAP GENERATION
# ===============================================
# Centralized configuration to avoid duplication across scripts

# ===============================================
# CONSTANTS
# ===============================================

MATRICES_DIR <- file.path("6_matrices")
CONSOLIDATED_BASE_DIR <- "7_Heatmap_Outputs"
MASTER_REFERENCE <- "All_Smel_Genes"

# Output subdirectories - must match bash configuration
HEATMAP_SUBDIR <- "I_Basic_Heatmap"
CV_HEATMAP_SUBDIR <- "II_Heatmap_with_CV"
BAR_GRAPH_SUBDIR <- "III_Bar_Graphs"

COUNT_TYPES <- c("expected_count")
GENE_TYPES <- c(
  #"geneID", 
  "geneName"
)
LABEL_TYPES <- c("SRR", "Organ")
PROCESSING_LEVELS <- c("gene_level", "isoform_level")

NORM_SCHEMES <- c(
  #"count_type_normalized",
  #"zscore",
  "zscore_scaled_to_ten"
)

# Sample labels mapping
SAMPLE_LABELS <- c(
  # Roots
  "SRR3884675" = "Roots",
  # Stems
  "SRR3884690" = "Stems",
  # Leaves
  "SRR3884689" = "Leaves",
  "SRR3884684" = "Senescent_leaves",
  # Buds
  "SRR3884686" = "Buds",
  "SRR3884687" = "Opened_Buds",
  # Flowers
  "SRR3884597" = "Flowers",
  # Fruits
  "SRR3884631" = "Fruits",
  "SRR3884608" = "Fruits_1cm",
  "SRR3884620" = "Fruits_Stage_1",
  "SRR3884642" = "Fruits_Skin_Stage_2",
  "SRR3884653" = "Fruits_Flesh_Stage_2",
  "SRR3884664" = "Fruits_Calyx_Stage_2",
  "SRR3884680" = "Fruits_Skin_Stage_3",
  "SRR3884681" = "Fruits_Flesh_Stage_3",
  "SRR3884678" = "Fruits_peduncle",
  # Other organs
  "SRR3884685" = "Radicles",
  "SRR3884677" = "Cotyledons",
  "SRR3884679" = "Pistils"
)

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
    "b_modules_for_Method_5/.gene_groups_temp.txt",
    default_value = c("Selected_GRF_GIF_Genes_vAll_GIFs", 
                      "Selected_GRF_GIF_Genes_vTwo_GIFs")
  )
  
  overwrite <- read_config_file(
    "b_modules_for_Method_5/.overwrite_temp.txt",
    default_value = TRUE,
    is_boolean = TRUE
  )
  
  list(
    gene_groups = gene_groups,
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
  file.path(
    MATRICES_DIR, MASTER_REFERENCE, processing_level, gene_group,
    paste0(gene_group, "_", count_type, "_", gene_type,
           "_from_", MASTER_REFERENCE, "_", processing_level, ".tsv")
  )
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
