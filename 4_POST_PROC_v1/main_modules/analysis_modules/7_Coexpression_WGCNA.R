#!/usr/bin/env Rscript

# ===============================================
# WGCNA COEXPRESSION ANALYSIS MODULE
# ===============================================
# Weighted gene co-expression network analysis

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
  library(fastcluster)
  library(ggplot2)
  library(RColorBrewer)
  library(igraph)
})

allowWGCNAThreads()

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# CONFIGURATION
# ===============================================

WGCNA_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$WGCNA)

# WGCNA parameters (optimized for accuracy with 24GB RAM)
SOFT_POWER_RANGE <- 1:30           # Extended range for better fit
# MIN_MODULE_SIZE: Set dynamically based on gene count (default 30, min 5)
# For small gene groups (<100), use smaller module size
MIN_MODULE_SIZE_DEFAULT <- 30
MIN_MODULE_SIZE_SMALL <- 5
MERGE_CUT_HEIGHT <- 0.15           # Standard merge threshold
NETWORK_TYPE <- "signed"           # "signed" preserves biological interpretation (more accurate)
DEEP_SPLIT <- 2                    # Module detection sensitivity (0-4, higher = more modules)
PAM_STAGE <- TRUE                  # PAM for more accurate module assignment
# Note: MIN_GENES_WGCNA is defined in 0_shared_config.R (default: 20)

# Coexpression parameters (increased for comprehensive analysis)
N_HUB_GENES <- 5                   # Top hub genes per module
N_COEXPRESSED_GENES <- 100         # More coexpressed genes per query
N_NETWORK_GENES <- 50              # Larger network visualization
COR_THRESHOLD <- 0.95              # Slightly lower threshold for more connections

# Figure toggles (all enabled for comprehensive output)
GENERATE_WGCNA_FIGURES <- list(
  soft_threshold_plot = TRUE,      # Scale-free topology fit
  module_dendrogram = TRUE,        # Gene clustering dendrogram
  eigengene_adjacency = TRUE,      # Module relationships
  correlation_network = TRUE       # Network visualization
)

# ===============================================
# WGCNA CORE FUNCTIONS
# ===============================================

pick_soft_threshold <- function(data_matrix, output_dir, gene_group) {
  powers <- SOFT_POWER_RANGE
  sft <- pickSoftThreshold(data_matrix, powerVector = powers, 
                           networkType = NETWORK_TYPE, verbose = 0)
  
  if (GENERATE_WGCNA_FIGURES$soft_threshold_plot) {
    png(file.path(output_dir, paste0(gene_group, "_soft_threshold.png")), 
        width = 1000, height = 500, res = 100)
    par(mfrow = c(1, 2))
    fit_index <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
    plot(sft$fitIndices[,1], fit_index,
         xlab = "Soft Threshold", ylab = "R^2",
         type = "n", main = "Scale Independence")
    text(sft$fitIndices[,1], fit_index, labels = powers, col = "red")
    abline(h = 0.80, col = "red")
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab = "Soft Threshold", ylab = "Mean Connectivity",
         type = "n", main = "Mean Connectivity")
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red")
    dev.off()
  }
  
  fit_index <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  power_selected <- which(fit_index > 0.80)[1]
  if (is.na(power_selected)) power_selected <- which.max(fit_index)
  
  return(powers[power_selected])
}

build_network_and_detect_modules <- function(data_matrix, soft_power, output_dir, gene_group) {
  # Dynamic module size based on gene count
  n_genes <- ncol(data_matrix)  # Genes are columns after transpose for WGCNA
  min_module_size <- if (n_genes < 100) MIN_MODULE_SIZE_SMALL else MIN_MODULE_SIZE_DEFAULT
  
  # Use PAM and deepSplit for more accurate module detection
  deep_split_val <- if (exists("DEEP_SPLIT")) DEEP_SPLIT else 2
  pam_stage_val <- if (exists("PAM_STAGE")) PAM_STAGE else TRUE
  
  net <- blockwiseModules(
    data_matrix, 
    power = soft_power,
    networkType = NETWORK_TYPE,
    TOMType = "signed",
    minModuleSize = min_module_size,
    reassignThreshold = 0,
    mergeCutHeight = MERGE_CUT_HEIGHT,
    deepSplit = deep_split_val,        # Sensitivity for module detection (0-4)
    pamStage = pam_stage_val,          # PAM for accurate gene assignment
    pamRespectsDendro = TRUE,          # PAM respects dendrogram structure
    numericLabels = TRUE,
    saveTOMs = TRUE,                   # Save TOM for downstream analysis
    saveTOMFileBase = file.path(output_dir, paste0(gene_group, "_TOM")),
    verbose = 1                        # Show progress
  )
  
  module_colors <- labels2colors(net$colors)
  
  if (GENERATE_WGCNA_FIGURES$module_dendrogram) {
    png(file.path(output_dir, paste0(gene_group, "_module_dendrogram.png")), 
        width = 1200, height = 600, res = 100)
    plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                        "Module colors", dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = paste0(gene_group, " - Module Dendrogram"))
    dev.off()
  }
  
  return(list(net = net, module_colors = module_colors))
}

calculate_module_eigengenes <- function(data_matrix, module_colors, output_dir, gene_group) {
  ME_list <- moduleEigengenes(data_matrix, colors = module_colors)
  MEs <- ME_list$eigengenes
  
  if (ncol(MEs) >= 2 && GENERATE_WGCNA_FIGURES$eigengene_adjacency) {
    png(file.path(output_dir, paste0(gene_group, "_eigengene_adjacency.png")), 
        width = 800, height = 800, res = 100)
    plotEigengeneNetworks(MEs, "", marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2),
                          plotDendrograms = TRUE, xLabelsAngle = 90)
    dev.off()
  }
  
  return(MEs)
}

identify_hub_genes <- function(gene_info, n_top = N_HUB_GENES) {
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
  
  return(hub_genes)
}

create_correlation_network <- function(data_matrix, query_genes, output_dir, gene_group) {
  if (!GENERATE_WGCNA_FIGURES$correlation_network) return(NULL)
  
  # Compute correlation matrix (use GPU if available)
  cor_matrix <- gpu_cor(t(data_matrix))
  
  # Filter to query genes and their top correlated genes
  matched <- query_genes[query_genes %in% rownames(cor_matrix)]
  if (length(matched) == 0) return(NULL)
  
  # Get top correlated genes for each query gene
  all_genes <- matched
  for (qg in matched) {
    cors <- cor_matrix[qg, ]
    cors <- cors[!names(cors) %in% matched]
    top_cors <- names(sort(abs(cors), decreasing = TRUE)[1:min(N_NETWORK_GENES, length(cors))])
    all_genes <- unique(c(all_genes, top_cors))
  }
  
  # Create adjacency for network
  sub_cor <- cor_matrix[all_genes, all_genes]
  adj <- abs(sub_cor)
  adj[adj < COR_THRESHOLD] <- 0
  diag(adj) <- 0
  
  # Create igraph network
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
  V(g)$is_query <- V(g)$name %in% matched
  
  # Plot
  png(file.path(output_dir, paste0(gene_group, "_correlation_network.png")),
      width = 1200, height = 1000, res = 100)
  
  colors <- ifelse(V(g)$is_query, "#B2182B", "#2166AC")
  sizes <- ifelse(V(g)$is_query, 12, 8)
  
  plot(g, 
       vertex.color = colors,
       vertex.size = sizes,
       vertex.label.cex = 0.6,
       edge.width = E(g)$weight * 2,
       main = paste0(gene_group, " - Correlation Network"))
  
  dev.off()
  
  return(g)
}

# ===============================================
# MAIN WGCNA FUNCTION
# ===============================================

run_wgcna <- function(config = NULL, matrices_dir = MATRICES_DIR) {
  if (is.null(config)) config <- load_runtime_config()
  ensure_output_dir(WGCNA_OUT_DIR)
  
  print_config_summary("WGCNA COEXPRESSION ANALYSIS", config)
  
  successful <- 0
  total <- 0
  
  for (gene_group in config$gene_groups) {
    cat("Processing:", gene_group, "\n")
    total <- total + 1
    
    output_dir <- file.path(WGCNA_OUT_DIR, gene_group)
    ensure_output_dir(output_dir)
    
    # Read matrix
    input_file <- build_input_path(gene_group, PROCESSING_LEVELS[1], 
                                   COUNT_TYPES[1], GENE_TYPES[1],
                                   matrices_dir, config$master_reference)
    
    validation <- validate_and_read_matrix(input_file, MIN_GENES_WGCNA)
    if (!validation$success) {
      cat("  Skipped:", validation$reason, "\n")
      next
    }
    
    # WGCNA requires variance-stabilized data and variance filtering
    # Apply log2(count+1) transformation for variance stabilization
    data_log <- log2(validation$data + 1)
    
    # Filter genes with low variance (uninformative for network construction)
    gene_vars <- apply(data_log, 1, var, na.rm = TRUE)
    var_threshold <- quantile(gene_vars, 0.25, na.rm = TRUE)  # Remove bottom 25%
    data_filtered <- data_log[gene_vars > var_threshold, , drop = FALSE]
    
    if (nrow(data_filtered) < MIN_GENES_WGCNA) {
      cat("  Skipped: Too few genes after variance filtering\n")
      next
    }
    cat("  Genes after variance filtering:", nrow(data_filtered), "\n")
    
    data_matrix <- t(data_filtered)  # WGCNA expects samples as rows
    
    # Run WGCNA pipeline
    soft_power <- pick_soft_threshold(data_matrix, output_dir, gene_group)
    cat("  Soft power:", soft_power, "\n")
    
    network <- build_network_and_detect_modules(data_matrix, soft_power, output_dir, gene_group)
    cat("  Modules detected:", length(unique(network$module_colors)), "\n")
    
    MEs <- calculate_module_eigengenes(data_matrix, network$module_colors, output_dir, gene_group)
    
    # Export gene-module assignments
    kME <- signedKME(data_matrix, MEs, outputColumnName = "kME")
    gene_info <- data.frame(
      Gene = colnames(data_matrix),
      Module = network$module_colors,
      stringsAsFactors = FALSE
    )
    gene_info <- cbind(gene_info, kME)
    
    write.table(gene_info, file.path(output_dir, paste0(gene_group, "_module_assignments.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Hub genes
    hubs <- identify_hub_genes(gene_info)
    write.table(hubs, file.path(output_dir, paste0(gene_group, "_hub_genes.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Correlation network
    create_correlation_network(t(data_matrix), colnames(data_matrix), output_dir, gene_group)
    
    successful <- successful + 1
    cat("  Complete\n")
  }
  
  print_summary(successful, total)
}

if (!interactive() && identical(environment(), globalenv())) {
  run_wgcna()
}
