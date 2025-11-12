# ===============================================
# SHARED CONFIGURATION FOR METHOD 4 HEATMAP GENERATION
# ===============================================
# Purpose: Centralized configuration to maintain DRY principle
# Contains:
#   - Constants for directory paths and processing parameters
#   - Sample-to-organ label mappings
#   - Helper functions for configuration and validation
#   - Reusable utility functions for path building and logging
# ===============================================

# ===============================================
# DIRECTORY AND FILE CONSTANTS
# ===============================================

MATRICES_DIR <- file.path("6_matrices_from_Salmon")  # tximport output location
CONSOLIDATED_BASE_DIR <- "7_Heatmap_Outputs"          # All visualizations
MASTER_REFERENCE <- "All_Smel_Genes"                  # Reference dataset name

# ===============================================
# PROCESSING PARAMETERS
# ===============================================

# Count type from Salmon (expected_count is standard)
COUNT_TYPES <- c("expected_count")

# Gene identifier types (ID vs human-readable name)
GENE_TYPES <- c(
  #"geneID", 
  "geneName"
)

# Sample label types (SRR accession vs organ name)
LABEL_TYPES <- c("SRR", "Organ")

# Processing granularity (gene-level summary vs isoform-level detail)
PROCESSING_LEVELS <- c("gene_level", "isoform_level")

# Normalization schemes for visualization
NORM_SCHEMES <- c(
  #"count_type_normalized",    # Log2 transformation
  #"zscore",                   # Z-score standardization
  "zscore_scaled_to_ten"      # Z-score scaled to [0,10] range
)

# ===============================================
# SAMPLE-TO-ORGAN LABEL MAPPING
# ===============================================
# Maps SRR accessions to biological organ/tissue names
# Used for generating human-readable heatmap labels
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
# HELPER FUNCTIONS FOR CONFIGURATION MANAGEMENT
# ===============================================

# Read configuration from temporary files written by bash wrapper
# Args:
#   file_path: Path to config file
#   default_value: Fallback if file doesn't exist
#   is_boolean: Whether to parse as TRUE/FALSE
# Returns: Configuration value(s)
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

# Load runtime configuration from temporary files
# Returns: List with gene_groups and overwrite_existing settings
load_runtime_config <- function() {
  gene_groups <- read_config_file(
    "b_modules_for_Method_4/.gene_groups_temp.txt",
    default_value = c("Selected_GRF_GIF_Genes_vAll_GIFs", 
                      "Selected_GRF_GIF_Genes_vTwo_GIFs")
  )
  
  overwrite <- read_config_file(
    "b_modules_for_Method_4/.overwrite_temp.txt",
    default_value = TRUE,
    is_boolean = TRUE
  )
  
  list(
    gene_groups = gene_groups,
    overwrite_existing = overwrite
  )
}

# ===============================================
# OUTPUT FORMATTING UTILITIES
# ===============================================

# Print formatted separator line
# Args:
#   char: Character to use for line (default "=")
#   width: Line width in characters (default 60)
print_separator <- function(char = "=", width = 60) {
  cat("\n", paste(rep(char, width), collapse = ""), "\n")
}

# Print configuration summary header
# Args:
#   title: Analysis title to display
#   config: Configuration list from load_runtime_config()
print_config_summary <- function(title, config) {
  print_separator()
  cat(title, "\n")
  print_separator()
  cat("\nConfiguration:\n")
  cat("  • Overwrite existing files:", config$overwrite_existing, "\n")
  cat("  • Gene groups:", paste(config$gene_groups, collapse = ", "), "\n\n")
}

# Print completion summary with statistics
# Args:
#   successful: Number of successful operations
#   total: Total operations attempted
#   skipped: Number of skipped operations (optional)
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

# ===============================================
# PATH CONSTRUCTION UTILITIES
# ===============================================

# Build input file path for count matrix
# Args: gene_group, processing_level, count_type, gene_type
# Returns: Full path to TSV matrix file
build_input_path <- function(gene_group, processing_level, count_type, gene_type) {
  file.path(
    MATRICES_DIR, MASTER_REFERENCE, processing_level, gene_group,
    paste0(gene_group, "_", count_type, "_", gene_type,
           "_from_", MASTER_REFERENCE, "_", processing_level, ".tsv")
  )
}

# Build standardized title for output files
# Args: gene_group, count_type, gene_type, label_type, processing_level, norm_scheme
# Returns: Descriptive title string
build_title_base <- function(gene_group, count_type, gene_type, label_type, processing_level, norm_scheme) {
  paste0(gene_group, "_", count_type, "_", gene_type, "_",
         label_type, "_from_", MASTER_REFERENCE, "_", processing_level, "_", norm_scheme)
}

# ===============================================
# DATA VALIDATION UTILITIES
# ===============================================

# Validate and read count matrix with error handling
# Args:
#   input_file: Path to matrix file
#   min_rows: Minimum required rows (default 2)
# Returns: List with success status and data/reason
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

# Check if file should be skipped due to existing output
# Args:
#   output_path: Path to check
#   overwrite: Whether overwrite mode is enabled
# Returns: TRUE if should skip, FALSE if should process
should_skip_existing <- function(output_path, overwrite) {
  !overwrite && file.exists(output_path)
}

# ===============================================
# DIRECTORY MANAGEMENT
# ===============================================

# Create output directory safely with optional cleaning
# Args:
#   dir_path: Directory to create
#   clean: Whether to remove existing directory first (default FALSE)
ensure_output_dir <- function(dir_path, clean = FALSE) {
  if (clean && dir.exists(dir_path)) {
    unlink(dir_path, recursive = TRUE)
  }
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}
