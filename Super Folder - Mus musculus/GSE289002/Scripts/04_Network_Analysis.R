# Load required libraries
library(WGCNA)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)  # Add this for gene name conversion

# Set working directory
setwd("/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE289002")

# Create results directory
dir.create("Outputs/Network_Analysis", recursive = TRUE, showWarnings = FALSE)

# Load processed data - Using correct RDS files
normalized_counts <- readRDS("Outputs/Enhanced_QC/Processed_Data/dge_normalized_final.rds")
sample_metadata <- read.csv("Outputs/Enhanced_QC/Processed_Data/final_metadata.csv", row.names = 1)

# Extract expression matrix from DGEList object
if (class(normalized_counts) == "DGEList") {
  normalized_counts <- normalized_counts$counts
}

cat("Data loaded successfully!\n")
cat("Genes:", nrow(normalized_counts), "\n")
cat("Samples:", ncol(normalized_counts), "\n\n")

# ============================================================================
# GENE NAME CONVERSION SETUP
# ============================================================================

cat("Setting up gene name conversion...\n")

# Function to convert Ensembl IDs to gene names
convert_ensembl_to_gene_names <- function(ensembl_ids) {
  # Remove any non-Ensembl IDs (like __alignment_not_unique, etc.)
  valid_ensembl <- ensembl_ids[grepl("^ENSMUSG", ensembl_ids)]
  
  if (length(valid_ensembl) == 0) {
    # If no valid Ensembl IDs, return original IDs
    gene_map <- data.frame(
      ensembl_gene_id = ensembl_ids,
      external_gene_name = ensembl_ids,
      stringsAsFactors = FALSE
    )
    return(gene_map)
  }
  
  tryCatch({
    # Connect to Ensembl biomart
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    # Get gene names
    gene_info <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "description"),
      filters = "ensembl_gene_id",
      values = valid_ensembl,
      mart = mart
    )
    
    # Create a mapping for all input IDs
    gene_map <- data.frame(
      ensembl_gene_id = ensembl_ids,
      external_gene_name = ensembl_ids,  # Default to original ID
      stringsAsFactors = FALSE
    )
    
    # Update with gene names where available
    for (i in 1:nrow(gene_info)) {
      idx <- which(gene_map$ensembl_gene_id == gene_info$ensembl_gene_id[i])
      if (length(idx) > 0 && gene_info$external_gene_name[i] != "") {
        gene_map$external_gene_name[idx] <- gene_info$external_gene_name[i]
      }
    }
    
    return(gene_map)
    
  }, error = function(e) {
    cat("Warning: Could not connect to biomaRt. Using Ensembl IDs.\n")
    gene_map <- data.frame(
      ensembl_gene_id = ensembl_ids,
      external_gene_name = ensembl_ids,
      stringsAsFactors = FALSE
    )
    return(gene_map)
  })
}

# Get gene name mapping for all genes
all_gene_ids <- rownames(normalized_counts)
gene_name_map <- convert_ensembl_to_gene_names(all_gene_ids)

# Create lookup function
get_gene_name <- function(ensembl_id) {
  idx <- match(ensembl_id, gene_name_map$ensembl_gene_id)
  if (is.na(idx)) {
    return(ensembl_id)
  } else {
    return(gene_name_map$external_gene_name[idx])
  }
}

cat("Gene name conversion setup completed!\n")
cat("Converted", sum(gene_name_map$ensembl_gene_id != gene_name_map$external_gene_name), 
    "out of", nrow(gene_name_map), "gene IDs\n\n")

# ============================================================================
# 1. CO-EXPRESSION NETWORK ANALYSIS USING WGCNA
# ============================================================================

cat("Starting WGCNA co-expression network analysis...\n")

# Prepare data for WGCNA (genes as columns)
wgcna_data <- t(normalized_counts)

# Ensure data is numeric (not integer) to avoid REAL() error
wgcna_data <- apply(wgcna_data, 2, as.numeric)

# Check for genes with too many missing values
gsg <- goodSamplesGenes(wgcna_data, verbose = 3)
if (!gsg$allOK) {
  wgcna_data <- wgcna_data[gsg$goodSamples, gsg$goodGenes]
}

cat("After quality filtering:\n")
cat("Samples:", nrow(wgcna_data), "\n")
cat("Genes:", ncol(wgcna_data), "\n\n")

# Filter genes by variance to reduce dataset size (keep top 5000 most variable genes)
cat("Filtering for most variable genes to reduce computational load...\n")
gene_vars <- apply(wgcna_data, 2, var, na.rm = TRUE)
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:5000])
wgcna_data <- wgcna_data[, top_var_genes]

cat("After variance filtering:\n")
cat("Samples:", nrow(wgcna_data), "\n")
cat("Genes:", ncol(wgcna_data), "\n\n")

# Choose soft threshold power
cat("Selecting soft threshold power...\n")
powers <- c(c(1:10), seq(from = 12, to=20, by=2))

# Enable multi-threading for WGCNA
enableWGCNAThreads()

sft <- pickSoftThreshold(wgcna_data, powerVector = powers, verbose = 5)

# Plot scale independence and mean connectivity
pdf("Outputs/Network_Analysis/soft_threshold_selection.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Select optimal power
soft_power <- sft$powerEstimate
if (is.na(soft_power)) {
  # Find the first power where R^2 > 0.8
  r2_threshold <- 0.8
  suitable_powers <- which(-sign(sft$fitIndices[,3])*sft$fitIndices[,2] > r2_threshold)
  if (length(suitable_powers) > 0) {
    soft_power <- powers[suitable_powers[1]]
  } else {
    soft_power <- 6  # Default fallback
  }
}

cat("Selected soft threshold power:", soft_power, "\n\n")

# Construct co-expression network and detect modules with reduced parameters for large datasets
cat("Constructing co-expression network and detecting modules...\n")

# Use one-step network construction for better memory management
net <- blockwiseModules(wgcna_data, 
                       power = soft_power,
                       TOMType = "unsigned", 
                       minModuleSize = 50,  # Increased from 30
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, 
                       saveTOMFileBase = "Outputs/Network_Analysis/TOM",
                       maxBlockSize = 3000,  # Reduced block size
                       verbose = 3)

# Module colors - IMPORTANT: Use the correct gene names
module_colors <- labels2colors(net$colors)
names(module_colors) <- colnames(wgcna_data)

cat("Module detection completed!\n")
cat("Number of modules:", length(unique(module_colors)) - 1, "\n") # Excluding grey
print(table(module_colors))

# Plot dendrogram with module colors
pdf("Outputs/Network_Analysis/module_dendrogram.pdf", width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                   "Module colors", dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05)
dev.off()

# ============================================================================
# 2. HUB GENE IDENTIFICATION - WITH GENE NAME CONVERSION
# ============================================================================

cat("Identifying hub genes...\n")

# Calculate module membership and gene significance
module_membership <- as.data.frame(cor(wgcna_data, net$MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(module_membership), nrow(wgcna_data)))

# Create the full adjacency matrix manually for connectivity calculations
cat("Creating adjacency matrix for connectivity calculations...\n")
adjacency_full <- adjacency(wgcna_data, power = soft_power)

# Calculate connectivity within modules using the full adjacency matrix
connectivity <- intramodularConnectivity(adjacency_full, module_colors)

cat("Connectivity calculation completed.\n")

# FIXED: Identify hub genes with proper name preservation
hub_genes_list <- list()
all_modules <- unique(module_colors)
cat("Processing", length(all_modules), "modules\n")

for (module in all_modules) {
  cat("Processing module:", module, "\n")
  
  if (module != "grey") {  # grey module contains unassigned genes
    module_genes <- names(module_colors)[module_colors == module]
    cat("Module", module, "has", length(module_genes), "genes\n")
    
    if (length(module_genes) > 0) {  # Check if module has genes
      # FIXED METHOD: Use setNames to preserve gene names
      module_connectivity <- setNames(connectivity[module_genes, "kWithin"], module_genes)
      
      cat("Connectivity range for", module, ":", range(module_connectivity), "\n")
      
      # Calculate number of hub genes based on module size
      if (length(module_connectivity) >= 1000) {
        percent_hubs <- 0.02  # 2% for very large modules
      } else if (length(module_connectivity) >= 500) {
        percent_hubs <- 0.03  # 3% for large modules
      } else if (length(module_connectivity) >= 100) {
        percent_hubs <- 0.05  # 5% for medium-large modules
      } else if (length(module_connectivity) >= 50) {
        percent_hubs <- 0.10  # 10% for medium modules
      } else if (length(module_connectivity) >= 20) {
        percent_hubs <- 0.15  # 15% for small modules
      } else {
        percent_hubs <- 0.20  # 20% for very small modules
      }
      
      # Calculate number of hubs based on percentage
      n_hubs_percent <- ceiling(length(module_connectivity) * percent_hubs)
      
      # Set absolute minimums based on module size
      if (length(module_connectivity) >= 100) {
        min_hubs <- 5
      } else if (length(module_connectivity) >= 50) {
        min_hubs <- 3
      } else if (length(module_connectivity) >= 10) {
        min_hubs <- 2
      } else {
        min_hubs <- 1
      }
      
      # Take the maximum of percentage-based and minimum required
      n_hubs <- max(min_hubs, n_hubs_percent)
      
      # Ensure we don't exceed total genes in module
      n_hubs <- min(n_hubs, length(module_connectivity))
      
      # Ensure at least 1 hub gene
      n_hubs <- max(1, n_hubs)
      
      cat("Calculated n_hubs for", module, ":", n_hubs, "(from", length(module_connectivity), "genes)\n")
      
      # Select hub genes - this will now work because names are preserved
      sorted_connectivity <- sort(module_connectivity, decreasing = TRUE)
      hub_genes <- names(sorted_connectivity[1:n_hubs])
      
      hub_genes_list[[module]] <- hub_genes
      cat("Selected", length(hub_genes), "hub genes for module", module, "\n")
      
      # Convert to gene names for display
      hub_gene_names <- sapply(hub_genes[1:min(3, length(hub_genes))], get_gene_name)
      cat("Top 3 hub genes:", paste(hub_gene_names, collapse = ", "), "\n")
    }
  }
}

# Create hub genes summary WITH GENE NAMES
if (length(hub_genes_list) > 0 && sum(sapply(hub_genes_list, length)) > 0) {
  all_hub_ensembl <- unlist(hub_genes_list)
  hub_genes_summary <- data.frame(
    Ensembl_ID = all_hub_ensembl,
    Gene_Name = sapply(all_hub_ensembl, get_gene_name),
    Module = rep(names(hub_genes_list), sapply(hub_genes_list, length)),
    Connectivity = connectivity[all_hub_ensembl, "kWithin"],
    stringsAsFactors = FALSE
  )
  rownames(hub_genes_summary) <- NULL
} else {
  cat("Warning: No hub genes identified. Creating empty summary.\n")
  hub_genes_summary <- data.frame(
    Ensembl_ID = character(0),
    Gene_Name = character(0),
    Module = character(0),
    Connectivity = numeric(0)
  )
}

cat("Hub genes identified:", nrow(hub_genes_summary), "\n\n")

# ============================================================================
# 3. MODULE CHARACTERIZATION
# ============================================================================

cat("Characterizing modules...\n")

# Calculate module eigengenes
MEs <- moduleEigengenes(wgcna_data, module_colors)$eigengenes
ME_dissim <- 1 - cor(MEs)

# Cluster module eigengenes
ME_tree <- hclust(as.dist(ME_dissim), method = "average")

# Plot module relationships
pdf("Outputs/Network_Analysis/module_eigengene_clustering.pdf", width = 12, height = 9)
plot(ME_tree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col = "red")
dev.off()

# Module-trait relationships using processed metadata
sample_info <- sample_metadata

# Check and clean metadata for model matrix creation
cat("Preparing metadata for trait analysis...\n")

# Remove columns with only one unique value or all missing values
valid_cols <- sapply(sample_info, function(x) {
  if (is.factor(x) || is.character(x)) {
    length(unique(x[!is.na(x)])) > 1
  } else {
    var(x, na.rm = TRUE) > 0 && !all(is.na(x))
  }
})

sample_info_clean <- sample_info[, valid_cols, drop = FALSE]

if (ncol(sample_info_clean) > 0) {
  # Convert categorical traits to numeric
  trait_data <- model.matrix(~ . - 1, data = sample_info_clean)
  
  # Calculate correlations
  module_trait_cor <- cor(MEs, trait_data, use = "p")
  module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(trait_data))
  
  # Plot heatmap
  pdf("Outputs/Network_Analysis/module_trait_relationships.pdf", width = 10, height = 8)
  textMatrix <- paste(signif(module_trait_cor, 2), "\n(",
                     signif(module_trait_pvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(module_trait_cor)
  
  labeledHeatmap(Matrix = module_trait_cor,
                xLabels = colnames(trait_data),
                yLabels = names(MEs),
                ySymbols = names(MEs),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix,
                setStdMargins = FALSE,
                cex.text = 0.5,
                zlim = c(-1,1),
                main = paste("Module-trait relationships"))
  dev.off()
} else {
  cat("Warning: No valid traits found for correlation analysis.\n")
  module_trait_cor <- NULL
  module_trait_pvalue <- NULL
}

# ============================================================================
# 4. NETWORK VISUALIZATION - WITH GENE NAMES
# ============================================================================

cat("Creating network visualizations...\n")

# Convert to igraph object for visualization
TOM_full <- TOMsimilarity(adjacency_full)

# Create a simplified network for visualization (using top connections)
# Select top 1% of TOM values to avoid overly dense networks
tom_threshold <- quantile(TOM_full[upper.tri(TOM_full)], 0.99)
TOM_viz <- TOM_full
TOM_viz[TOM_viz < tom_threshold] <- 0

# Only create network if there are connections
if (sum(TOM_viz > 0) > 0) {
  coexpr_network <- graph_from_adjacency_matrix(TOM_viz, mode = "undirected", 
                                               weighted = TRUE, diag = FALSE)
  
  # Add module information to vertices
  if (vcount(coexpr_network) > 0) {
    # Get vertex names that exist in the network
    vertex_names <- V(coexpr_network)$name
    
    # Assign module colors for vertices that exist
    V(coexpr_network)$module <- module_colors[vertex_names]
    
    # Assign hub status
    if (nrow(hub_genes_summary) > 0 && length(hub_genes_summary$Ensembl_ID) > 0) {
      all_hub_genes <- hub_genes_summary$Ensembl_ID
      hub_status <- ifelse(vertex_names %in% all_hub_genes, "Hub", "Regular")
    } else {
      hub_status <- rep("Regular", length(vertex_names))
    }
    
    # Validate length before assignment
    if (length(hub_status) == length(vertex_names)) {
      V(coexpr_network)$hub <- hub_status
      cat("Network created with", vcount(coexpr_network), "vertices and", ecount(coexpr_network), "edges\n")
      cat("Hub genes in network:", sum(hub_status == "Hub"), "\n")
    } else {
      cat("Warning: Length mismatch in hub assignment. Using all Regular.\n")
      V(coexpr_network)$hub <- rep("Regular", length(vertex_names))
    }
  }
} else {
  cat("Warning: No connections above threshold for network visualization.\n")
  coexpr_network <- NULL
}

# Plot network for largest modules (limit to 2 to avoid memory issues)
largest_modules <- names(sort(table(module_colors), decreasing = TRUE))[1:2]
largest_modules <- largest_modules[largest_modules != "grey"]

for (module in largest_modules) {
  module_genes <- names(module_colors)[module_colors == module]
  
  if (length(module_genes) > 1) {  # Only process if module has multiple genes
    # Create subnetwork for this module
    module_adj <- TOM_full[module_genes, module_genes]
    module_threshold <- quantile(module_adj[upper.tri(module_adj)], 0.95)
    module_adj[module_adj < module_threshold] <- 0
    
    # Only create network if there are connections
    if (sum(module_adj > 0) > 0) {
      module_network <- graph_from_adjacency_matrix(module_adj, mode = "undirected", 
                                                   weighted = TRUE, diag = FALSE)
      
      # Only plot if network has edges
      if (ecount(module_network) > 0) {
        pdf(paste0("Outputs/Network_Analysis/network_", module, "_module.pdf"), width = 12, height = 12)
        
        # Set vertex attributes with proper error handling
        module_hub_genes <- if (module %in% names(hub_genes_list) && length(hub_genes_list[[module]]) > 0) {
          hub_genes_list[[module]]
        } else {
          character(0)
        }
        
        module_vertex_names <- V(module_network)$name
        module_hub_status <- ifelse(module_vertex_names %in% module_hub_genes, "Hub", "Regular")
        
        # Ensure proper assignment with length validation
        if (length(module_hub_status) == length(module_vertex_names)) {
          V(module_network)$hub_status <- module_hub_status
        } else {
          cat("Warning: Length mismatch in module", module, "hub assignment.\n")
          V(module_network)$hub_status <- rep("Regular", length(module_vertex_names))
        }
        
        # Convert hub gene labels to gene names for display
        hub_gene_labels <- ifelse(V(module_network)$hub_status == "Hub", 
                                 sapply(V(module_network)$name, get_gene_name), 
                                 NA)
        
        plot(module_network,
             vertex.color = module,
             vertex.size = ifelse(V(module_network)$hub_status == "Hub", 8, 4),
             vertex.label = hub_gene_labels,
             vertex.label.cex = 0.6,
             vertex.label.color = "black",
             edge.width = E(module_network)$weight * 3,
             edge.color = "grey70",
             layout = layout_with_fr,
             main = paste("Co-expression Network -", module, "Module"))
        
        # Add legend
        legend("topright", legend = c("Hub Gene", "Regular Gene"), 
               pch = 21, pt.bg = module, pt.cex = c(2, 1), bty = "n")
        
        dev.off()
        cat("Network plot saved for", module, "module\n")
      }
    }
  }
}

# ============================================================================
# 5. SAVE RESULTS - WITH GENE NAMES
# ============================================================================

cat("Saving results...\n")

# Save all network analysis results
save(net, module_colors, hub_genes_list, hub_genes_summary, gene_name_map,
     MEs, module_trait_cor, coexpr_network, adjacency_full, TOM_full,
     file = "Outputs/Network_Analysis/network_analysis_results.RData")

# Export hub genes WITH GENE NAMES
write.csv(hub_genes_summary, "Outputs/Network_Analysis/hub_genes.csv", row.names = FALSE)

# Export module assignments WITH GENE NAMES
module_assignments <- data.frame(
  Ensembl_ID = names(module_colors),
  Gene_Name = sapply(names(module_colors), get_gene_name),
  Module = module_colors,
  stringsAsFactors = FALSE
)
write.csv(module_assignments, "Outputs/Network_Analysis/module_assignments.csv", row.names = FALSE)

# Export gene name mapping
write.csv(gene_name_map, "Outputs/Network_Analysis/gene_name_mapping.csv", row.names = FALSE)

# Export module eigengenes
write.csv(MEs, "Outputs/Network_Analysis/module_eigengenes.csv", row.names = TRUE)

# Export module-trait correlations if available
if (!is.null(module_trait_cor)) {
  write.csv(module_trait_cor, "Outputs/Network_Analysis/module_trait_correlations.csv", row.names = TRUE)
}

# Create comprehensive summary report WITH GENE NAMES
sink("Outputs/Network_Analysis/network_analysis_summary.txt")
cat("================================================================\n")
cat("NETWORK ANALYSIS SUMMARY REPORT\n")
cat("================================================================\n\n")
cat("Analysis Date:", as.character(Sys.time()), "\n\n")

cat("INPUT DATA:\n")
cat("- Original genes:", nrow(normalized_counts), "\n")
cat("- Genes analyzed (top 5000 most variable):", ncol(wgcna_data), "\n")
cat("- Samples analyzed:", nrow(wgcna_data), "\n")
cat("- Gene names converted:", sum(gene_name_map$ensembl_gene_id != gene_name_map$external_gene_name), 
    "out of", nrow(gene_name_map), "\n\n")

cat("WGCNA RESULTS:\n")
cat("- Soft threshold power used:", soft_power, "\n")
cat("- Total modules detected:", length(unique(module_colors)) - 1, "\n")
cat("- Hub genes identified:", nrow(hub_genes_summary), "\n\n")

cat("MODULE BREAKDOWN:\n")
module_table <- table(module_colors)
for (i in 1:length(module_table)) {
  cat(names(module_table)[i], "module:", module_table[i], "genes\n")
}

if (nrow(hub_genes_summary) > 0) {
  cat("\nHUB GENES BY MODULE (with gene names):\n")
  for (module in names(hub_genes_list)) {
    cat(module, "module:", length(hub_genes_list[[module]]), "hub genes\n")
    # Show top 5 hub genes with names
    module_hubs <- hub_genes_summary[hub_genes_summary$Module == module, ]
    top_hubs <- head(module_hubs[order(module_hubs$Connectivity, decreasing = TRUE), ], 5)
    for (j in 1:nrow(top_hubs)) {
      cat("  ", j, ". ", top_hubs$Gene_Name[j], " (", top_hubs$Ensembl_ID[j], 
          ") - Connectivity: ", round(top_hubs$Connectivity[j], 2), "\n", sep="")
    }
  }
} else {
  cat("\nNo hub genes identified with current thresholds.\n")
}

cat("\nFILES GENERATED:\n")
cat("- network_analysis_results.RData: Complete analysis object\n")
cat("- hub_genes.csv: Hub gene list with gene names and connectivity scores\n") 
cat("- module_assignments.csv: Gene-to-module assignments with gene names\n")
cat("- gene_name_mapping.csv: Ensembl ID to gene name mapping\n")
cat("- module_eigengenes.csv: Module eigengene values\n")
if (!is.null(module_trait_cor)) {
  cat("- module_trait_correlations.csv: Module-trait correlation matrix\n")
}
cat("- soft_threshold_selection.pdf: WGCNA parameter selection\n")
cat("- module_dendrogram.pdf: Module clustering dendrogram\n")
cat("- module_eigengene_clustering.pdf: Module relationships\n")
if (!is.null(module_trait_cor)) {
  cat("- module_trait_relationships.pdf: Module-trait correlations heatmap\n")
}
cat("- network_[module]_module.pdf: Individual module network plots with gene names\n")

sink()

cat("Network analysis completed successfully!\n")
cat("Results saved in Outputs/Network_Analysis/\n")
cat("Key outputs:\n")
cat("- Hub genes identified:", nrow(hub_genes_summary), "\n")
cat("- Modules detected:", length(unique(module_colors)) - 1, "\n")
cat("- Gene names converted:", sum(gene_name_map$ensembl_gene_id != gene_name_map$external_gene_name), 
    "out of", nrow(gene_name_map), "\n")
cat("- Network files and visualizations saved with gene names\n")

# Disable multi-threading
disableWGCNAThreads()