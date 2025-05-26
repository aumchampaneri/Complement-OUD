# Cell type annotation using original GSE207128 paper markers with integration compatibility

# Set working directory
setwd("/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128/")
input_dir <- "Outputs/Processed Data"
output_dir <- "Outputs/Processed Data/Marker Based Outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(patchwork)

print("=== CELL TYPE ANNOTATION WORKFLOW ===")
print("Using original GSE207128 paper markers with integration compatibility")

# =============================================================================
# 1. LOAD BEST AVAILABLE INTEGRATED DATA
# =============================================================================

print("Loading best available Seurat object...")

# Try to load the best available processed object (in order of preference)
if (file.exists(file.path(input_dir, "harmony_integrated_seurat.rds"))) {
  seurat_obj <- readRDS(file.path(input_dir, "harmony_integrated_seurat.rds"))
  integration_method <- "Harmony"
  print("✓ Loaded Harmony-integrated data")
} else if (file.exists(file.path(input_dir, "batch_corrected_seurat.rds"))) {
  seurat_obj <- readRDS(file.path(input_dir, "batch_corrected_seurat.rds"))
  integration_method <- "Batch Regression"
  print("✓ Loaded batch-corrected data")
} else if (file.exists(file.path(input_dir, "fastmnn_integrated_seurat.rds"))) {
  seurat_obj <- readRDS(file.path(input_dir, "fastmnn_integrated_seurat.rds"))
  integration_method <- "FastMNN"
  print("✓ Loaded FastMNN-integrated data")
} else if (file.exists(file.path(input_dir, "processed_seurat.rds"))) {
  seurat_obj <- readRDS(file.path(input_dir, "processed_seurat.rds"))
  integration_method <- "Standard Processing"
  print("✓ Loaded standard processed data")
} else {
  stop("ERROR: No processed Seurat object found! Please run data processing first.")
}

# Print basic info
print(paste("Number of cells:", ncol(seurat_obj)))
print(paste("Number of genes:", nrow(seurat_obj)))
print(paste("Integration method:", integration_method))

# Ensure RNA assay is default for marker analysis
DefaultAssay(seurat_obj) <- "RNA"

# Handle Seurat v5 objects - join layers if needed
if (inherits(seurat_obj[["RNA"]], "Assay5")) {
  print("Detected Seurat v5 object. Joining layers...")
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
}

# Check if clusters exist, if not create them
if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  print("No existing clusters found. Running clustering on integrated data...")
  
  # Determine which reduction to use based on integration method
  if ("harmony" %in% names(seurat_obj@reductions)) {
    reduction_for_clustering <- "harmony"
    print("Using Harmony reduction for clustering")
  } else if ("fastmnn" %in% names(seurat_obj@reductions)) {
    reduction_for_clustering <- "fastmnn"
    print("Using FastMNN reduction for clustering")
  } else {
    reduction_for_clustering <- "pca"
    print("Using PCA reduction for clustering")
  }
  
  # Run clustering
  seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_for_clustering, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
  print(paste("Created", length(unique(seurat_obj$seurat_clusters)), "clusters"))
} else {
  print(paste("Using existing", length(unique(seurat_obj$seurat_clusters)), "clusters"))
}

# =============================================================================
# 2. DEFINE ORIGINAL PAPER MARKERS
# =============================================================================

print("Setting up marker genes from original GSE207128 paper...")

# EXACT markers from the original paper
cell_type_markers <- list(
  ASC = "Gja1",           # Astrocytes
  OPC = "Pdgfra",         # Oligodendrocyte Precursor Cells
  MG = "Tmem119",         # Microglia
  EC = "Cldn5",           # Endothelial Cells
  NEUR = "Syt1",          # Neurons
  OLG = c("Cldn11", "Mobp"), # Oligodendrocytes
  MAC = "Pf4",            # Macrophages
  PC = "Vtn",             # Pericytes
  NFOLG = "Enpp6",        # Newly Formed Oligodendrocytes
  VSMC = "Acta2",         # Vascular Smooth Muscle Cells
  DC = c("Cd74", "Cd209a"), # Dendritic Cells
  EPC = "Ccdc153",        # Ependymal Cells
  NSC = "Thbs4",          # Neural Stem Cells
  ARP = "Cd44",           # Arachnoid Barrier-like Cells
  NRP = c("Top2a", "Cdk1"), # Neural Restricted Progenitors
  Tcells = "Cd3d",        # T cells
  NEUT = "S100a9",        # Neutrophils
  EPIC = "Ttr",           # Epithelial Cells (contamination)
  VLMC = "Slc6a13",       # Vascular Leptomeningeal Cells (contamination)
  ABC = "Slc47a1"         # Arachnoid Barrier Cells (contamination)
)

# Additional neuron subtype markers for refinement
neuron_subtype_markers <- list(
  Excitatory = c("Slc17a7", "Camk2a", "Grin2a"),
  Inhibitory = c("Gad1", "Gad2", "Slc32a1"),
  PV = "Pvalb",
  SST = "Sst", 
  VIP = "Vip",
  CCK = "Cck"
)

# Astrocyte subtype markers
astrocyte_subtype_markers <- list(
  Protoplasmic = "Aqp4",
  Fibrous = "Gfap",
  Reactive = c("Lcn2", "Serpina3n")
)

print(paste("Defined", length(cell_type_markers), "cell type markers"))

# =============================================================================
# 3. FIND CLUSTER MARKERS
# =============================================================================

print("Finding cluster-specific markers...")
print("This may take a while for large datasets...")

# Find markers for all clusters (following the paper's approach)
all_markers <- FindAllMarkers(seurat_obj, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25,
                             verbose = FALSE)

# Get top markers per cluster
top_markers_per_cluster <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# Save all markers
write.csv(all_markers, file.path(output_dir, "all_cluster_markers.csv"), row.names = FALSE)
write.csv(top_markers_per_cluster, file.path(output_dir, "top_cluster_markers.csv"), row.names = FALSE)

print(paste("Found", nrow(all_markers), "significant cluster markers"))

# =============================================================================
# 4. CHECK MARKER AVAILABILITY
# =============================================================================

print("Checking marker gene availability...")
all_marker_genes <- unique(unlist(cell_type_markers))
available_markers <- all_marker_genes[all_marker_genes %in% rownames(seurat_obj)]
missing_markers <- all_marker_genes[!all_marker_genes %in% rownames(seurat_obj)]

print(paste("Available markers:", length(available_markers), "out of", length(all_marker_genes)))
if (length(missing_markers) > 0) {
  print("Missing markers:")
  print(missing_markers)
  
  # Save missing markers info
  missing_df <- data.frame(
    marker = missing_markers,
    cell_type = sapply(missing_markers, function(x) {
      names(cell_type_markers)[sapply(cell_type_markers, function(y) x %in% y)]
    })
  )
  write.csv(missing_df, file.path(output_dir, "missing_markers.csv"), row.names = FALSE)
}

# =============================================================================
# 5. CALCULATE MARKER EXPRESSION (ORIGINAL PAPER METHOD)
# =============================================================================

print("Calculating marker expression scores using original paper method...")
marker_expression <- data.frame(row.names = levels(factor(seurat_obj$seurat_clusters)))

for (cell_type in names(cell_type_markers)) {
  markers <- cell_type_markers[[cell_type]]
  available_markers_for_type <- markers[markers %in% rownames(seurat_obj)]
  
  if (length(available_markers_for_type) > 0) {
    # Calculate average expression per cluster (EXACT ORIGINAL METHOD)
    avg_exp <- AverageExpression(seurat_obj, 
                                features = available_markers_for_type, 
                                group.by = "seurat_clusters")$RNA
    
    if (length(available_markers_for_type) == 1) {
      marker_expression[[cell_type]] <- avg_exp[available_markers_for_type, ]
    } else {
      marker_expression[[cell_type]] <- colMeans(avg_exp[available_markers_for_type, ])
    }
  } else {
    marker_expression[[cell_type]] <- 0
    print(paste("Warning: No markers available for", cell_type))
  }
}

# =============================================================================
# 6. ASSIGN CELL TYPES (ORIGINAL PAPER METHOD)
# =============================================================================

print("Assigning cell types using original paper criteria...")
cluster_assignments <- data.frame(
  cluster = rownames(marker_expression),
  predicted_type = apply(marker_expression, 1, function(x) {
    if (max(x) > 0) {
      names(x)[which.max(x)]
    } else {
      "Unknown"
    }
  }),
  max_score = apply(marker_expression, 1, max),
  stringsAsFactors = FALSE
)

# Add confidence threshold (following paper's approach)
confidence_threshold <- 0.5  # As used in original paper
cluster_assignments$confident <- cluster_assignments$max_score > confidence_threshold

# Add second-best score for ambiguity assessment
cluster_assignments$second_best_score <- apply(marker_expression, 1, function(x) {
  if (length(x) > 1) {
    sort(x, decreasing = TRUE)[2]
  } else {
    0
  }
})

cluster_assignments$score_difference <- cluster_assignments$max_score - cluster_assignments$second_best_score

# Add top cluster markers for manual review
cluster_assignments$top_markers <- sapply(cluster_assignments$cluster, function(x) {
  top_genes <- top_markers_per_cluster %>% 
    filter(cluster == x) %>% 
    pull(gene) %>% 
    head(5)
  paste(top_genes, collapse = ", ")
})

# =============================================================================
# 7. ADD ANNOTATIONS TO SEURAT OBJECT
# =============================================================================

print("Adding annotations to Seurat object...")

# Create cell type assignments for individual cells
seurat_obj$predicted_celltype <- cluster_assignments$predicted_type[
  match(seurat_obj$seurat_clusters, cluster_assignments$cluster)]
seurat_obj$celltype_confidence <- cluster_assignments$max_score[
  match(seurat_obj$seurat_clusters, cluster_assignments$cluster)]
seurat_obj$confident_prediction <- cluster_assignments$confident[
  match(seurat_obj$seurat_clusters, cluster_assignments$cluster)]

# Filter out contaminating cell types (EXACTLY as in original paper)
contaminating_types <- c("EPIC", "VLMC", "ABC")
seurat_obj$is_contaminant <- seurat_obj$predicted_celltype %in% contaminating_types

print("Cell type distribution:")
print(table(seurat_obj$predicted_celltype))

# =============================================================================
# 8. CREATE COMPREHENSIVE VISUALIZATIONS
# =============================================================================

print("Creating comprehensive visualizations...")

# 1. UMAP with predicted cell types
p1 <- DimPlot(seurat_obj, 
              group.by = "predicted_celltype", 
              label = TRUE, 
              pt.size = 0.5,
              repel = TRUE) +
  ggtitle(paste("Predicted Cell Types (", integration_method, ")", sep = "")) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 3, override.aes = list(size = 3)))

# 2. Confidence scoring
p2 <- FeaturePlot(seurat_obj, 
                  features = "celltype_confidence", 
                  pt.size = 0.5) +
  ggtitle("Cell Type Prediction Confidence") +
  scale_color_viridis_c(name = "Confidence")

# 3. Contamination identification
p3 <- DimPlot(seurat_obj, 
              group.by = "is_contaminant", 
              pt.size = 0.5,
              cols = c("FALSE" = "lightgray", "TRUE" = "red")) +
  ggtitle("Contaminating Cell Types (To Filter)") +
  labs(color = "Contaminant")

# 4. Original clusters
p4 <- DimPlot(seurat_obj, 
              group.by = "seurat_clusters", 
              label = TRUE, 
              pt.size = 0.5) +
  ggtitle("Seurat Clusters") +
  theme(legend.position = "none")

# 5. Sample distribution (if available)
if ("sample_id" %in% colnames(seurat_obj@meta.data)) {
  p5 <- DimPlot(seurat_obj, 
                group.by = "sample_id", 
                pt.size = 0.3) +
    ggtitle("Sample Distribution") +
    theme(legend.position = "bottom")
  
  # 6. Condition distribution (if available)
  if ("condition" %in% colnames(seurat_obj@meta.data)) {
    p6 <- DimPlot(seurat_obj, 
                  group.by = "condition", 
                  pt.size = 0.3) +
      ggtitle("Condition Distribution") +
      theme(legend.position = "bottom")
    
    # Combined overview plot
    overview_plot <- (p1 | p2) / (p3 | p4) / (p5 | p6)
    plot_height <- 24
  } else {
    overview_plot <- (p1 | p2) / (p3 | p4) / p5
    plot_height <- 18
  }
} else {
  overview_plot <- (p1 | p2) / (p3 | p4)
  plot_height <- 12
}

# Save individual plots
ggsave(file.path(output_dir, "predicted_celltypes_UMAP.pdf"), p1, width = 14, height = 10)
ggsave(file.path(output_dir, "confidence_scores.pdf"), p2, width = 10, height = 8)
ggsave(file.path(output_dir, "contaminants_to_filter.pdf"), p3, width = 10, height = 8)
ggsave(file.path(output_dir, "original_clusters.pdf"), p4, width = 10, height = 8)

# Save combined overview
ggsave(file.path(output_dir, "annotation_overview.pdf"), overview_plot, width = 16, height = plot_height)

# =============================================================================
# 9. MARKER EXPRESSION HEATMAP
# =============================================================================

print("Creating enhanced marker expression heatmap...")
pdf(file.path(output_dir, "marker_expression_heatmap.pdf"), width = 14, height = 10)

# Add cluster annotations for the heatmap
cluster_annotation <- data.frame(
  Cluster = rownames(marker_expression),
  Predicted_Type = cluster_assignments$predicted_type,
  Confidence = ifelse(cluster_assignments$confident, "High", "Low"),
  Integration = integration_method
)
rownames(cluster_annotation) <- cluster_annotation$Cluster

# Create color scheme
n_types <- length(unique(cluster_annotation$Predicted_Type))
cell_type_colors <- rainbow(n_types)
names(cell_type_colors) <- unique(cluster_annotation$Predicted_Type)

annotation_colors <- list(
  Predicted_Type = cell_type_colors,
  Confidence = c("High" = "darkgreen", "Low" = "orange"),
  Integration = c("Harmony" = "blue", "Batch Regression" = "purple", 
                 "FastMNN" = "darkred", "Standard Processing" = "gray")
)

pheatmap(as.matrix(marker_expression), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         scale = "column",
         annotation_row = cluster_annotation[, c("Predicted_Type", "Confidence", "Integration")],
         annotation_colors = annotation_colors,
         main = paste("Average Marker Expression by Cluster\n(", integration_method, " Integration)", sep = ""),
         angle_col = 45,
         fontsize = 10)
dev.off()

# =============================================================================
# 10. KEY MARKER FEATURE PLOTS
# =============================================================================

print("Creating feature plots for key markers...")
key_markers <- c("Gja1", "Syt1", "Tmem119", "Pdgfra", "Cldn5", "Mobp", "Pf4", "Cldn11", "Acta2")
available_key_markers <- key_markers[key_markers %in% rownames(seurat_obj)]

if (length(available_key_markers) > 0) {
  p_features <- FeaturePlot(seurat_obj, 
                           features = available_key_markers, 
                           ncol = 3, 
                           pt.size = 0.3)
  ggsave(file.path(output_dir, "key_markers_featurePlot.pdf"), p_features, 
         width = 15, height = ceiling(length(available_key_markers)/3) * 4)
}

# =============================================================================
# 11. SAVE DATA AND SUMMARY TABLES
# =============================================================================

print("Creating summary tables and saving data...")

# Cell type frequency tables
celltype_freq <- table(seurat_obj$predicted_celltype)
cluster_celltype_table <- table(seurat_obj$seurat_clusters, seurat_obj$predicted_celltype)

# Confidence summary by cell type
confidence_summary <- seurat_obj@meta.data %>%
  group_by(predicted_celltype) %>%
  summarize(
    n_cells = n(),
    mean_confidence = mean(celltype_confidence),
    sd_confidence = sd(celltype_confidence),
    confident_cells = sum(confident_prediction),
    confident_percent = round(sum(confident_prediction)/n()*100, 1)
  )

# Save tables
write.csv(as.data.frame(celltype_freq), file.path(output_dir, "celltype_frequencies.csv"))
write.csv(cluster_celltype_table, file.path(output_dir, "cluster_celltype_assignments.csv"))
write.csv(cluster_assignments, file.path(output_dir, "cluster_annotations.csv"), row.names = FALSE)
write.csv(marker_expression, file.path(output_dir, "marker_expression_by_cluster.csv"))
write.csv(confidence_summary, file.path(output_dir, "confidence_summary_by_celltype.csv"), row.names = FALSE)

# Summary statistics
n_cells_total <- ncol(seurat_obj)
n_contaminants <- sum(seurat_obj$is_contaminant)
n_confident <- sum(seurat_obj$confident_prediction)
n_celltypes <- length(unique(seurat_obj$predicted_celltype[!seurat_obj$is_contaminant]))

summary_report <- data.frame(
  Metric = c("Total Cells", "Contaminant Cells", "Clean Cells", 
             "Confident Predictions", "Cell Types Identified", 
             "Contamination Rate", "Confidence Rate", "Integration Method"),
  Value = c(n_cells_total, n_contaminants, n_cells_total - n_contaminants,
            n_confident, n_celltypes,
            paste0(round(n_contaminants/n_cells_total*100, 1), "%"),
            paste0(round(n_confident/n_cells_total*100, 1), "%"),
            integration_method)
)

write.csv(summary_report, file.path(output_dir, "analysis_summary.csv"), row.names = FALSE)

# =============================================================================
# 12. CREATE FILTERED DATASET
# =============================================================================

print("Creating filtered dataset (removing contaminants)...")
seurat_filtered <- subset(seurat_obj, subset = is_contaminant == FALSE)

print("Saving annotated objects...")
# Save both original and filtered objects
saveRDS(seurat_obj, file.path(output_dir, "seurat_obj_annotated.rds"))
saveRDS(seurat_filtered, file.path(output_dir, "seurat_obj_filtered.rds"))

# Also save to main processed data directory for easy access
saveRDS(seurat_obj, file.path(input_dir, "annotated_seurat.rds"))
saveRDS(seurat_filtered, file.path(input_dir, "filtered_annotated_seurat.rds"))

# =============================================================================
# 13. FINAL SUMMARY REPORT
# =============================================================================

print("\n=== CELL TYPE ANNOTATION COMPLETED ===")
print(paste("Integration method used:", integration_method))
print(paste("Results saved in:", output_dir))

print("\nSummary:")
print(paste("- Total cells analyzed:", n_cells_total))
print(paste("- Contaminant cells identified:", n_contaminants, 
            "(", round(n_contaminants/n_cells_total*100, 1), "%)"))
print(paste("- Clean cells retained:", n_cells_total - n_contaminants))
print(paste("- Cell types identified:", n_celltypes))
print(paste("- Confident predictions:", n_confident, 
            "(", round(n_confident/n_cells_total*100, 1), "%)"))

print("\nCell type frequencies (excluding contaminants):")
clean_freq <- table(seurat_filtered$predicted_celltype)
print(sort(clean_freq, decreasing = TRUE))

print("\nCluster to cell type mapping:")
cluster_summary <- cluster_assignments %>%
  arrange(desc(max_score)) %>%
  select(cluster, predicted_type, max_score, confident, top_markers)
print(cluster_summary)

print("\nFiles created:")
print("- seurat_obj_annotated.rds (with contaminants)")
print("- seurat_obj_filtered.rds (contaminants removed)")
print("- annotated_seurat.rds (main directory)")
print("- filtered_annotated_seurat.rds (main directory)")
print("- annotation_overview.pdf")
print("- marker_expression_heatmap.pdf")
print("- key_markers_featurePlot.pdf")
print("- cluster_annotations.csv")
print("- all_cluster_markers.csv")
print("- analysis_summary.csv")

print("\nNext steps:")
print("1. Review annotation plots for accuracy")
print("2. Use filtered_annotated_seurat.rds for downstream analysis")
print("3. Proceed with differential expression analysis")
print("4. Investigate condition-specific cell type changes")

# Clean up
gc()

print("\n=== CELL TYPE ANNOTATION WORKFLOW COMPLETE ===")
