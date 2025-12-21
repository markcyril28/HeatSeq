#!/usr/bin/env Rscript

# ===============================================
# WGCNA COEXPRESSION ANALYSIS FOR METHOD 5 (BOWTIE2/RSEM)
# ===============================================
# Performs weighted gene co-expression network analysis using WGCNA
# Identifies modules of co-expressed genes and hub genes

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
  library(fastcluster)
  library(ggplot2)
  library(RColorBrewer)
  library(igraph)
})

# Enable multi-threading for WGCNA
allowWGCNAThreads()

source("B_M5_modules/0_shared_config.R")
source("B_M5_modules/2_utility_functions.R")

# ===============================================
# CONFIGURATION
# ===============================================

WGCNA_SUBDIR <- "IV_Coexpression_WGCNA"
WGCNA_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, WGCNA_SUBDIR)

# WGCNA parameters
SOFT_POWER_RANGE <- 1:20           # Range for soft threshold testing. Higher = stronger emphasis on high correlations
MIN_MODULE_SIZE <- 15              # Min genes/module. ↑ = fewer but more robust modules; ↓ = more but smaller modules
MERGE_CUT_HEIGHT <- 0.15           # Module merge threshold (0-1). ↑ = more merging (fewer modules); ↓ = distinct modules
NETWORK_TYPE <- "signed"           # "signed" (considers direction: +/-) or "unsigned" (magnitude only)
MIN_GENES_FOR_WGCNA <- 20          # Min genes for analysis. ↑ = stringent quality; ↓ = more permissive

# Hub gene and coexpression parameters
N_HUB_GENES <- 3                   # Top hubs/module. ↓ = most connected only; ↑ = include secondary hubs
N_COEXPRESSED_GENES <- 10          # Coexpressed genes/query. ↓ = strongest associations; ↑ = broader network
N_NETWORK_GENES <- 30              # Max genes in network plot. ↓ = focused visualization; ↑ = comprehensive but cluttered
COR_THRESHOLD <- 0.9               # Edge correlation cutoff (0-1). ↑ = only strongest; ↓ = include weaker connections

# Load runtime configuration
config <- load_runtime_config()
ensure_output_dir(WGCNA_OUT_DIR)

# FIGURE GENERATION CONTROL (set to FALSE to skip)
GENERATE_FIGURES <- list(
  soft_threshold_plot = FALSE,           # Soft power selection diagnostic plots
  module_dendrogram = FALSE,             # Gene clustering tree with module colors
  eigengene_adjacency = FALSE,           # Module-module relationship heatmap
  module_sizes_barplot = FALSE,          # Number of genes per module
  correlation_network = TRUE            # Top hub genes network graph
)

# ===============================================
# WGCNA CORE FUNCTIONS
# ===============================================

# Pick optimal soft-thresholding power
pick_soft_threshold <- function(data_matrix, output_dir, gene_group) {
  powers <- SOFT_POWER_RANGE
  sft <- pickSoftThreshold(data_matrix, powerVector = powers, 
                           networkType = NETWORK_TYPE, verbose = 0)
  
  # Plot scale-free topology fit
  if (GENERATE_FIGURES$soft_threshold_plot) {
    png(file.path(output_dir, paste0(gene_group, "_soft_threshold.png")), 
        width = 1000, height = 500, res = 100)
    par(mfrow = c(1, 2))
    
    # Scale-free topology fit index
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
         xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (R^2)",
         type = "n", main = paste0(gene_group, " - Scale Independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
         labels = powers, col = "red")
    abline(h = 0.80, col = "red")
    
    # Mean connectivity
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
         type = "n", main = paste0(gene_group, " - Mean Connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red")
    
    dev.off()
  }
  
  # Select power: first to reach R^2 > 0.80, or highest R^2 if none reach threshold
  fit_index <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  power_selected <- which(fit_index > 0.80)[1]
  if (is.na(power_selected)) {
    power_selected <- which.max(fit_index)
  }
  
  return(powers[power_selected])
}

# Build network and detect modules
build_network_and_detect_modules <- function(data_matrix, soft_power, output_dir, gene_group) {
  # One-step network construction and module detection
  net <- blockwiseModules(
    data_matrix, 
    power = soft_power,
    networkType = NETWORK_TYPE,
    TOMType = "signed",
    minModuleSize = MIN_MODULE_SIZE,
    reassignThreshold = 0,
    mergeCutHeight = MERGE_CUT_HEIGHT,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 0
  )
  
  # Convert numeric labels to colors
  module_colors <- labels2colors(net$colors)
  
  # Plot dendrogram with module colors
  if (GENERATE_FIGURES$module_dendrogram) {
    png(file.path(output_dir, paste0(gene_group, "_module_dendrogram.png")), 
        width = 1200, height = 600, res = 100)
    plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                        "Module colors", dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = paste0(gene_group, " - Gene Dendrogram and Module Colors"))
    dev.off()
  }
  
  return(list(net = net, module_colors = module_colors))
}

# Calculate module eigengenes and create heatmap
create_module_eigengene_heatmap <- function(data_matrix, module_colors, output_dir, gene_group) {
  ME_list <- moduleEigengenes(data_matrix, colors = module_colors)
  MEs <- ME_list$eigengenes
  
  if (ncol(MEs) < 2) {
    cat("    Only one module detected, skipping eigengene heatmap\n")
    return(MEs)
  }
  
  # Eigengene adjacency heatmap
  if (GENERATE_FIGURES$eigengene_adjacency) {
    png(file.path(output_dir, paste0(gene_group, "_eigengene_adjacency.png")), 
        width = 800, height = 800, res = 100)
    plotEigengeneNetworks(MEs, "", marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2),
                          plotDendrograms = TRUE, xLabelsAngle = 90)
    dev.off()
  }
  
  return(MEs)
}

# Calculate and export module membership (kME) and connectivity
calculate_module_metrics <- function(data_matrix, module_colors, MEs, output_dir, gene_group) {
  # Module membership (kME)
  kME <- signedKME(data_matrix, MEs, outputColumnName = "kME")
  
  # For large datasets, skip full connectivity calculation to avoid memory issues
  # Use kME as the primary connectivity metric instead
  gene_info <- data.frame(
    Gene = colnames(data_matrix),
    Module = module_colors,
    stringsAsFactors = FALSE
  )
  gene_info <- cbind(gene_info, kME)
  
  # Sort by module then by module membership
  modules <- unique(gene_info$Module)
  for (mod in modules) {
    kME_col <- paste0("kME", mod)
    if (kME_col %in% colnames(gene_info)) {
      gene_info <- gene_info[order(gene_info$Module, -gene_info[[kME_col]]), ]
      break
    }
  }
  
  # Save gene module assignments
  write.table(gene_info, 
              file = file.path(output_dir, paste0(gene_group, "_gene_module_assignments.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(gene_info)
}

# Identify hub genes per module
identify_hub_genes <- function(gene_info, n_top = N_HUB_GENES, output_dir, gene_group) {
  modules <- unique(gene_info$Module)
  modules <- modules[modules != "grey"]
  
  hub_genes <- data.frame()
  
  for (mod in modules) {
    mod_genes <- gene_info[gene_info$Module == mod, ]
    kME_col <- paste0("kME", mod)
    
    if (kME_col %in% colnames(mod_genes)) {
      mod_genes <- mod_genes[order(-mod_genes[[kME_col]]), ]
      top_n <- min(n_top, nrow(mod_genes))
      hub_genes <- rbind(hub_genes, mod_genes[1:top_n, ])
    }
  }
  
  # Save hub genes
  write.table(hub_genes, 
              file = file.path(output_dir, paste0(gene_group, "_hub_genes.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(hub_genes)
}

# Create module summary barplot
create_module_summary_plot <- function(module_colors, output_dir, gene_group) {
  if (!GENERATE_FIGURES$module_sizes_barplot) return(invisible(NULL))
  
  module_counts <- table(module_colors)
  module_df <- data.frame(Module = names(module_counts), 
                          Count = as.numeric(module_counts))
  module_df <- module_df[order(-module_df$Count), ]
  module_df$Module <- factor(module_df$Module, levels = module_df$Module)
  
  p <- ggplot(module_df, aes(x = Module, y = Count, fill = Module)) +
    geom_bar(stat = "identity") +
    scale_fill_identity() +
    theme_minimal() +
    labs(title = paste0(gene_group, " - Genes per Module"),
         x = "Module", y = "Number of Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave(file.path(output_dir, paste0(gene_group, "_module_sizes.png")), 
         p, width = 10, height = 6, dpi = 100)
}

# Create gene correlation network graph
create_correlation_network <- function(data_matrix, module_colors, soft_power, 
                                        output_dir, gene_group, 
                                        top_n_genes = N_NETWORK_GENES, 
                                        cor_threshold = COR_THRESHOLD) {
  if (!GENERATE_FIGURES$correlation_network) return(invisible(NULL))
  
  # Limit to top hub genes per module to keep network manageable
  n_genes <- ncol(data_matrix)
  if (n_genes > top_n_genes) {
    # Calculate connectivity and select top genes
    adj <- adjacency(data_matrix, power = soft_power, type = NETWORK_TYPE)
    connectivity <- rowSums(adj) - 1
    top_idx <- order(connectivity, decreasing = TRUE)[1:top_n_genes]
    data_subset <- data_matrix[, top_idx]
    colors_subset <- module_colors[top_idx]
  } else {
    data_subset <- data_matrix
    colors_subset <- module_colors
  }
  
  # Calculate correlation matrix
  cor_mat <- cor(data_subset, use = "pairwise.complete.obs")
  diag(cor_mat) <- 0
  
  # Apply threshold to create adjacency
  adj_mat <- abs(cor_mat)
  adj_mat[adj_mat < cor_threshold] <- 0
  
  # Build igraph object
  graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Remove isolated nodes
  isolated <- which(degree(graph) == 0)
  if (length(isolated) > 0) {
    graph <- delete_vertices(graph, isolated)
    colors_subset <- colors_subset[-isolated]
  }
  
  if (vcount(graph) < 3) {
    cat("    Skipping network plot: too few connected genes\n")
    return(NULL)
  }
  
  # Set vertex attributes
  V(graph)$color <- colors_subset
  V(graph)$size <- scales::rescale(degree(graph), to = c(5, 15))
  V(graph)$label.cex <- 0.6
  
  # Set edge attributes based on correlation strength
  E(graph)$width <- scales::rescale(E(graph)$weight, to = c(0.5, 3))
  E(graph)$color <- adjustcolor("grey50", alpha.f = 0.5)
  
  # Plot network
  png(file.path(output_dir, paste0(gene_group, "_correlation_network.png")), 
      width = 1200, height = 1200, res = 120)
  
  set.seed(42)
  layout <- layout_with_fr(graph)
  
  plot(graph, layout = layout, vertex.frame.color = NA,
       edge.curved = 0.1, main = paste0(gene_group, " - Gene Correlation Network"),
       sub = paste0("Top ", vcount(graph), " genes, r >= ", cor_threshold))
  
  # Add legend for modules
  unique_colors <- unique(colors_subset)
  legend("bottomleft", legend = unique_colors, fill = unique_colors, 
         title = "Module", cex = 0.7, bty = "n")
  
  dev.off()
  
  return(graph)
}

# Generate comprehensive results summary file
generate_results_summary <- function(output_dir, gene_group, processing_level, count_type, gene_type,
                                     n_genes, n_modules, soft_power, hub_genes, gene_info) {
  summary_file <- file.path(output_dir, paste0(gene_group, "_WGCNA_RESULTS_SUMMARY.txt"))
  
  sink(summary_file)
  
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("WGCNA COEXPRESSION ANALYSIS RESULTS\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  cat("ANALYSIS INFORMATION\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("Gene Group:", gene_group, "\n")
  cat("Processing Level:", processing_level, "\n")
  cat("Count Type:", count_type, "\n")
  cat("Gene Type:", gene_type, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("PARAMETERS USED\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("Network Type:", NETWORK_TYPE, "\n")
  cat("Soft Threshold Power:", soft_power, "\n")
  cat("Min Module Size:", MIN_MODULE_SIZE, "genes\n")
  cat("Merge Cut Height:", MERGE_CUT_HEIGHT, "\n")
  cat("Correlation Threshold:", COR_THRESHOLD, "\n")
  cat("Number of Hub Genes per Module:", N_HUB_GENES, "\n")
  cat("Coexpressed Genes per Query:", N_COEXPRESSED_GENES, "\n\n")
  
  cat("RESULTS SUMMARY\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("Total Genes Analyzed:", n_genes, "\n")
  cat("Number of Modules Detected:", n_modules, "(excluding grey/unassigned)\n")
  
  module_counts <- table(gene_info$Module)
  grey_count <- ifelse("grey" %in% names(module_counts), module_counts["grey"], 0)
  cat("Unassigned Genes (grey module):", grey_count, sprintf("(%.1f%%)\n", 100*grey_count/n_genes))
  cat("Assigned Genes:", n_genes - grey_count, sprintf("(%.1f%%)\n\n", 100*(n_genes-grey_count)/n_genes))
  
  cat("MODULE COMPOSITION\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  modules <- names(module_counts)[names(module_counts) != "grey"]
  for (mod in modules) {
    count <- module_counts[mod]
    pct <- 100 * count / n_genes
    cat(sprintf("  %s: %d genes (%.1f%%)\n", mod, count, pct))
  }
  cat("\n")
  
  cat("HUB GENES (Top", N_HUB_GENES, "per module)\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  if (nrow(hub_genes) > 0) {
    current_module <- ""
    for (i in 1:nrow(hub_genes)) {
      if (hub_genes$Module[i] != current_module) {
        current_module <- hub_genes$Module[i]
        cat("\n", toupper(current_module), "MODULE:\n")
      }
      kME_col <- paste0("kME", current_module)
      kME_val <- ifelse(kME_col %in% colnames(hub_genes), sprintf("%.3f", hub_genes[i, kME_col]), "N/A")
      cat(sprintf("  %d. %s (kME = %s)\n", sum(hub_genes$Module[1:i] == current_module), hub_genes$Gene[i], kME_val))
    }
  } else {
    cat("  No hub genes identified\n")
  }
  cat("\n")
  
  cat("OUTPUT FILES\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("1. *_soft_threshold.png - Soft power selection plots\n")
  cat("2. *_module_dendrogram.png - Gene clustering with module colors\n")
  cat("3. *_eigengene_adjacency.png - Module relationship heatmap\n")
  cat("4. *_module_sizes.png - Genes per module barplot\n")
  cat("5. *_correlation_network.png - Top", N_NETWORK_GENES, "genes network (r >=", COR_THRESHOLD, ")\n")
  cat("6. *_gene_module_assignments.tsv - All genes with modules and kME values\n")
  cat("7. *_hub_genes.tsv - Top", N_HUB_GENES, "hubs per module\n")
  cat("8. *_module_eigengenes.tsv - Module representative expression profiles\n")
  cat("9. coexpression_*/* - Top", N_COEXPRESSED_GENES, "coexpressed genes per query\n\n")
  
  cat("INTERPRETATION\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("• Modules = co-expressed gene groups with similar patterns across samples\n")
  cat("• Hub genes = highly connected, potentially important regulators\n")
  cat("• kME (module membership): gene-module correlation (>0.8 strong, >0.7 moderate)\n")
  cat("• Grey module = unassigned genes not fitting any module\n\n")
  
  cat("NEXT STEPS\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("1. Perform GO/KEGG enrichment on modules\n")
  cat("2. Correlate module eigengenes with sample traits\n")
  cat("3. Validate hub genes experimentally\n")
  cat("4. Compare modules across conditions\n\n")
  
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  sink()
  
  cat("    Generated results summary:", basename(summary_file), "\n")
}

# Extract coexpressed genes for specific gene of interest
extract_coexpressed_genes <- function(gene_info, gene_of_interest, top_n = N_COEXPRESSED_GENES, output_dir, gene_group) {
  if (!gene_of_interest %in% gene_info$Gene) {
    cat("    Gene", gene_of_interest, "not found\n")
    return(NULL)
  }
  
  gene_module <- gene_info$Module[gene_info$Gene == gene_of_interest]
  
  if (gene_module == "grey") {
    cat("    Gene", gene_of_interest, "is unassigned (grey module)\n")
    return(NULL)
  }
  
  # Get all genes in same module
  module_genes <- gene_info[gene_info$Module == gene_module, ]
  kME_col <- paste0("kME", gene_module)
  
  # Sort by module membership
  if (kME_col %in% colnames(module_genes)) {
    module_genes <- module_genes[order(-module_genes[[kME_col]]), ]
  }
  
  # Get top coexpressed genes (excluding gene of interest itself)
  coexpressed <- module_genes[module_genes$Gene != gene_of_interest, ]
  coexpressed <- head(coexpressed, top_n)
  
  # Save results
  output_file <- file.path(output_dir, paste0(gene_group, "_coexpressed_with_", gene_of_interest, ".tsv"))
  write.table(coexpressed, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("    Found", nrow(coexpressed), "coexpressed genes with", gene_of_interest, 
      "in", gene_module, "module\n")
  
  return(coexpressed)
}

# ===============================================
# MAIN PROCESSING FUNCTION
# ===============================================

process_wgcna <- function(gene_group, processing_level, count_type, gene_type) {
  # Build input path
  input_file <- build_input_path(gene_group, processing_level, count_type, gene_type)
  
  # Validate and read matrix
  validation <- validate_and_read_matrix(input_file, min_rows = MIN_GENES_FOR_WGCNA)
  
  if (!validation$success) {
    cat("  Skipping WGCNA (", validation$reason, "):", basename(input_file), "\n", sep = "")
    return(list(success = FALSE))
  }
  
  cat("  Processing:", processing_level, "|", count_type, "|", gene_type, 
      "- Found", validation$n_genes, "genes\n")
  
  # WGCNA requires samples as rows, genes as columns
  data_matrix <- t(validation$data)
  
  # Check for sufficient samples
  if (nrow(data_matrix) < 3) {
    cat("    Skipping: insufficient samples (need >= 3, have", nrow(data_matrix), ")\n")
    return(list(success = FALSE))
  }
  
  # Create output directory
  output_dir <- file.path(WGCNA_OUT_DIR, gene_group, processing_level, count_type, gene_type)
  ensure_output_dir(output_dir)
  
  # Check for good genes (remove genes with zero variance)
  gsg <- goodSamplesGenes(data_matrix, verbose = 0)
  if (!gsg$allOK) {
    data_matrix <- data_matrix[gsg$goodSamples, gsg$goodGenes]
    cat("    Removed", sum(!gsg$goodGenes), "genes with zero variance\n")
  }
  
  if (ncol(data_matrix) < MIN_GENES_FOR_WGCNA) {
    cat("    Skipping: insufficient genes after filtering (need >=", MIN_GENES_FOR_WGCNA, ", have", ncol(data_matrix), ")\n")
    return(list(success = FALSE))
  }
  
  # Step 1: Pick soft threshold
  cat("    Picking soft threshold...\n")
  soft_power <- pick_soft_threshold(data_matrix, output_dir, gene_group)
  cat("    Selected soft power:", soft_power, "\n")
  
  # Step 2: Build network and detect modules
  cat("    Building network and detecting modules...\n")
  network_result <- build_network_and_detect_modules(data_matrix, soft_power, output_dir, gene_group)
  
  n_modules <- length(unique(network_result$module_colors)) - 1  # Exclude grey
  cat("    Detected", n_modules, "modules\n")
  
  # Step 3: Calculate module eigengenes
  cat("    Calculating module eigengenes...\n")
  MEs <- create_module_eigengene_heatmap(data_matrix, network_result$module_colors, output_dir, gene_group)
  
  # Step 4: Calculate module metrics
  cat("    Calculating module metrics...\n")
  gene_info <- calculate_module_metrics(data_matrix, network_result$module_colors, MEs, output_dir, gene_group)
  
  # Step 5: Identify hub genes
  cat("    Identifying hub genes...\n")
  hub_genes <- identify_hub_genes(gene_info, output_dir = output_dir, gene_group = gene_group)
  
  # Step 6: Create summary plots
  cat("    Creating summary plots...\n")
  create_module_summary_plot(network_result$module_colors, output_dir, gene_group)
  
  # Step 7: Create correlation network graph
  cat("    Creating correlation network graph...\n")
  create_correlation_network(data_matrix, network_result$module_colors, soft_power, 
                             output_dir, gene_group)
  
  # Step 8: Extract coexpressed genes for each gene group
  for (gene_grp in config$gene_groups) {
    cat("    Extracting coexpressed genes for:", gene_grp, "\n")
    
    # Load genes of interest from gene group file
    gene_group_file <- file.path("A_GeneGroups_InputList", 
                                  paste0("for_", gene_group), paste0(gene_grp, ".txt"))
    
    if (!file.exists(gene_group_file)) {
      cat("      Warning: Gene group file not found:", gene_group_file, "\n")
      next
    }
    
    genes_of_interest <- readLines(gene_group_file)
    genes_of_interest <- genes_of_interest[genes_of_interest != ""]
    
    cat("      Found", length(genes_of_interest), "genes of interest\n")
    
    # Create gene group output directory
    goi_output_dir <- file.path(output_dir, paste0("coexpression_", gene_grp))
    ensure_output_dir(goi_output_dir)
    
    for (goi in genes_of_interest) {
      extract_coexpressed_genes(gene_info, goi, output_dir = goi_output_dir, gene_group = gene_grp)
    }
  }
  
  # Save eigengene expression matrix (for external use)
  write.table(MEs, file = file.path(output_dir, paste0(gene_group, "_module_eigengenes.tsv")),
              sep = "\t", row.names = TRUE, quote = FALSE)
  
  # Generate comprehensive results summary
  cat("    Generating results summary...\n")
  generate_results_summary(
    output_dir = output_dir,
    gene_group = gene_group,
    processing_level = processing_level,
    count_type = count_type,
    gene_type = gene_type,
    n_genes = ncol(data_matrix),
    n_modules = n_modules,
    soft_power = soft_power,
    hub_genes = hub_genes,
    gene_info = gene_info
  )
  
  return(list(success = TRUE, n_modules = n_modules, n_genes = ncol(data_matrix)))
}

# ===============================================
# MAIN EXECUTION
# ===============================================

print_separator()
cat("WGCNA COEXPRESSION ANALYSIS - METHOD 5 (BOWTIE2/RSEM)\n")
print_separator()
cat("\nConfiguration:\n")
cat("  • Overwrite existing files:", config$overwrite_existing, "\n")
cat("  • Master reference:", config$master_reference, "\n")
cat("  • Genes of interest groups:", paste(config$gene_groups, collapse = ", "), "\n")
cat("  • Network type:", NETWORK_TYPE, "\n")
cat("  • Min module size:", MIN_MODULE_SIZE, "\n\n")

# Run WGCNA on full transcriptome (MASTER_REFERENCE), extract coexpression for GENE_GROUPS
total_processed <- 0
total_successful <- 0

# Use master_reference for full transcriptome analysis
master_ref <- config$master_reference

cat("\nRunning WGCNA on full transcriptome:", master_ref, "\n")
cat("Will extract coexpression partners for:", paste(config$gene_groups, collapse = ", "), "\n\n")

for (processing_level in PROCESSING_LEVELS) {
  for (count_type in COUNT_TYPES) {
    for (gene_type in GENE_TYPES) {
      total_processed <- total_processed + 1
      
      result <- tryCatch({
        process_wgcna(master_ref, processing_level, count_type, gene_type)
      }, error = function(e) {
        cat("    Error:", e$message, "\n")
        list(success = FALSE)
      })
      
      if (result$success) {
        total_successful <- total_successful + 1
      }
    }
  }
}

# Print summary
print_separator()
cat("SUMMARY:", total_successful, "/", total_processed, "analyses completed\n")
print_separator()

if (total_successful > 0) {
  print_output_summary(WGCNA_OUT_DIR, c(
    "Soft threshold selection plot", "Module dendrogram with colors",
    "Module eigengene adjacency heatmap", "Gene module assignments (TSV)",
    "Hub genes per module (TSV)", "Module eigengenes matrix (TSV)", 
    "Module size barplot", "Gene correlation network graph"
  ))
}
