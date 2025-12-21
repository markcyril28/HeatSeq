# ---
# Title: "Inputs, Configurations, and Outputs"
# ---

# --- Input Files ---
count_matrix_file <- "5_stringtie_WD/a_Method_1_RAW_RESULTs/transcript_count_matrix.csv"
sample_metadata_file <- "sample_metadata.txt"
gene_groups_dir <- "gene_groups"

# --- Output Directories ---
base_output_dir <- "6_DeSeq2_Method_1_Ref_Guided"
basic_heatmap_dir <- file.path(base_output_dir, "1_Basic_Heatmap")
cv_heatmap_dir <- file.path(base_output_dir, "2_Heatmap_with_CV")
bar_graphs_dir <- file.path(base_output_dir, "3_Bar_Graphs")

# --- Create Output Directories ---
dir.create(basic_heatmap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cv_heatmap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bar_graphs_dir, recursive = TRUE, showWarnings = FALSE)

# --- Analysis Configurations ---
normalization_schemes <- c("normalized", "zscore") # Add more if needed
