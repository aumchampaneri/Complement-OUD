# 03_Annotations.R - Production-ready marker-based brain cell annotation

# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(biomaRt)

#' Run marker-based cell type annotation on a Seurat object
#'
#' @param seurat_obj Seurat object with normalized data
#' @param markers Named list of marker genes for each cell type
#' @param output_dir Directory to save results
#' @param ctrl_size Number of control genes for scoring (default: 50)
#' @param score_threshold Minimum gap between top scores to avoid "Uncertain" label (default: 0.1)
#' @param seed Random seed for reproducibility (default: 123)
#' @param overwrite Whether to overwrite existing score columns (default: FALSE)
#' @return List containing: annotated Seurat object, composition tables, and score matrix
run_annotation <- function(seurat_obj,
                          markers,
                          output_dir = "annotation_results",
                          ctrl_size = 50,
                          score_threshold = 0.1,
                          seed = 123,
                          overwrite = FALSE) {

  # Set seed for reproducibility
  set.seed(seed)

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Function to log messages with timestamps
  log_info <- function(msg) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  }

  # Safety check for existing score columns
  existing <- grep("_score$", colnames(seurat_obj@meta.data), value=TRUE)
  if(length(existing) > 0 && !overwrite) {
    stop("Found existing score columns: ", paste(existing, collapse=", "),
         ". Use overwrite=TRUE to replace them or rename existing columns.")
  }

  #' Check marker genes present in dataset and log warnings for missing ones
  #'
  #' @param cell_type Name of the cell type
  #' @param markers Vector of marker gene names
  #' @param seurat_obj Seurat object
  #' @return Vector of markers present in the dataset
  check_markers <- function(cell_type, markers, seurat_obj) {
    present <- intersect(markers, rownames(seurat_obj))
    missing <- setdiff(markers, rownames(seurat_obj))

    if(length(missing) > 0) {
      log_info(sprintf("Warning: %d/%d markers missing for %s: %s",
                       length(missing), length(markers), cell_type,
                       paste(missing, collapse=", ")))
    }
    return(present)
  }

  #' Add module score and rename the column properly
  #'
  #' @param seurat_obj Seurat object
  #' @param cell_type Cell type name (used for column naming)
  #' @param markers Markers to score
  #' @param ctrl Number of control genes
  #' @param seed Random seed
  #' @return Updated Seurat object
  add_score <- function(seurat_obj, cell_type, markers, ctrl, seed) {
    seurat_obj <- Seurat::AddModuleScore(
      seurat_obj,
      features = list(markers),
      name = cell_type,
      ctrl = ctrl,
      seed = seed
    )

    # Rename the score column (remove trailing "1")
    old_name <- paste0(cell_type, "1")
    new_name <- paste0(cell_type, "_score")
    seurat_obj[[new_name]] <- seurat_obj[[old_name]]
    seurat_obj[[old_name]] <- NULL

    return(seurat_obj)
  }

  # Process marker genes and generate scores
  log_info("Beginning cell type annotation with marker genes")

  # Plot canonical markers before scoring to validate data
  canonical_markers <- c("Gfap", "Mbp", "P2ry12", "Snap25", "Cldn5")
  present_canonical <- intersect(canonical_markers, rownames(seurat_obj))

  if(length(present_canonical) > 0) {
    log_info(paste("Creating QC violin plots for", length(present_canonical), "canonical markers"))
    tryCatch({
      vln_pre <- Seurat::VlnPlot(seurat_obj, features = present_canonical,
                         group.by = "seurat_clusters", pt.size = 0, ncol = 3)
      ggplot2::ggsave(file.path(output_dir, "canonical_markers_vln.png"), vln_pre,
             width = min(15, length(present_canonical) * 3), height = 7)
    }, error = function(e) {
      warning("QC plot creation failed: ", e$message)
    })
  }

  # Score cells for each marker set
  log_info("Scoring cells with marker genes...")
  for(cell_type in names(markers)) {
    # Check markers and get present ones
    present_markers <- check_markers(cell_type, markers[[cell_type]], seurat_obj)

    if(length(present_markers) > 0) {
      log_info(sprintf("Scoring %s with %d markers", cell_type, length(present_markers)))
      seurat_obj <- add_score(seurat_obj, cell_type, present_markers, ctrl_size, seed)
    } else {
      log_info(sprintf("WARNING: No markers found for %s - adding zero score", cell_type))
      # Add placeholder score
      seurat_obj[[paste0(cell_type, "_score")]] <- 0
    }
  }

  # Get score columns and assign cell types
  log_info("Assigning cell types based on scores...")
  score_cols <- grep("_score$", colnames(seurat_obj@meta.data), value = TRUE)

  # Create a score dataframe for analysis
  scores_df <- as.matrix(seurat_obj@meta.data[, score_cols])

  # Calculate top types and score gaps for calling (vectorized approach)
  cell_types <- gsub("_score$", "", score_cols)
  top_indices <- max.col(scores_df, ties.method = "first")
  top_types <- cell_types[top_indices]

  # Vectorized score gap calculation
  if(requireNamespace("matrixStats", quietly = TRUE)) {
    # Efficient approach using matrixStats for large datasets
    ord <- t(apply(scores_df, 1, function(x) order(x, decreasing = TRUE)[1:2]))
    score_gaps <- scores_df[cbind(1:nrow(scores_df), ord[,1])] -
                  scores_df[cbind(1:nrow(scores_df), ord[,2])]
  } else {
    # Fallback to apply approach
    score_gaps <- apply(scores_df, 1, function(row) {
      if(length(row) <= 1) return(1)
      sorted <- sort(row, decreasing = TRUE)
      return(sorted[1] - sorted[2])
    })
  }

  # Assign cell types with threshold for ambiguous cells
  seurat_obj$predicted_celltype <- top_types
  seurat_obj$score_gap <- score_gaps
  seurat_obj$predicted_celltype[score_gaps < score_threshold] <- "Uncertain"

  # Make predicted_celltype a factor with consistent ordering
  seurat_obj$predicted_celltype <- factor(
    seurat_obj$predicted_celltype,
    levels = c(cell_types, "Uncertain")
  )

  # Create cluster-celltype composition tables
  log_info("Creating cluster composition tables...")
  cluster_type_table <- table(seurat_obj$seurat_clusters, seurat_obj$predicted_celltype)

  # Calculate percentages
  cluster_type_pct <- prop.table(cluster_type_table, margin = 1) * 100
  cluster_type_df <- as.data.frame.matrix(cluster_type_pct)
  cluster_type_df$Cluster <- rownames(cluster_type_df)

  # Save composition tables
  write.csv(cluster_type_table, file.path(output_dir, "cluster_celltype_counts.csv"))
  write.csv(cluster_type_pct, file.path(output_dir, "cluster_celltype_percent.csv"))

  # Create heat map of cluster composition
  tryCatch({
    if(requireNamespace("pheatmap", quietly = TRUE)) {
      pheatmap::pheatmap(cluster_type_pct,
                        filename = file.path(output_dir, "cluster_celltype_heatmap.png"),
                        display_numbers = TRUE,
                        number_format = "%.1f",
                        fontsize_number = 7,
                        width = 10,
                        height = 8)
    }
  }, error = function(e) {
    warning("Heatmap creation failed: ", e$message)
  })

  # Visualizations
  log_info("Creating visualizations...")
  tryCatch({
    p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "predicted_celltype",
                  label = TRUE, repel = TRUE) +
          ggplot2::ggtitle("Predicted Cell Types") +
          ggplot2::theme(legend.position = "right")

    p2 <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters",
                  label = TRUE, repel = TRUE) +
          ggplot2::ggtitle("Original Clusters")

    combined_plot <- p1 + p2
    ggplot2::ggsave(file.path(output_dir, "brain_cell_annotations.png"), combined_plot, width = 14, height = 8)
  }, error = function(e) {
    warning("UMAP plot creation failed: ", e$message)
  })

  # Plot marker scores with consistent cutoffs
  tryCatch({
    score_plots <- Seurat::FeaturePlot(seurat_obj,
                               features = score_cols,
                               ncol = 3,
                               min.cutoff = "q10",
                               max.cutoff = "q90",
                               order = TRUE)
    ggplot2::ggsave(file.path(output_dir, "marker_scores.png"), score_plots,
           width = 15, height = ceiling(length(score_cols)/3) * 5)
  }, error = function(e) {
    warning("Score plot creation failed: ", e$message)
  })

  # Create dotplot of top markers per type
  tryCatch({
    top_markers <- lapply(names(markers), function(ct) {
      present <- intersect(markers[[ct]], rownames(seurat_obj))
      if(length(present) > 0) head(present, 3) else NULL
    })

    top_markers_flat <- unlist(top_markers)
    if(length(top_markers_flat) > 0) {
      dotplot <- Seurat::DotPlot(seurat_obj, features = top_markers_flat,
                        group.by = "predicted_celltype") +
                ggplot2::RotatedAxis() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      ggplot2::ggsave(file.path(output_dir, "marker_dotplot.png"), dotplot, width = 12, height = 8)
    }
  }, error = function(e) {
    warning("Dotplot creation failed: ", e$message)
  })

  # Save results
  log_info("Saving annotated Seurat object...")
  saveRDS(seurat_obj, file.path(output_dir, "seurat_annotated.rds"))

  # Export cell annotations as CSV
  cell_annotations <- seurat_obj@meta.data %>%
    dplyr::select(seurat_clusters, predicted_celltype, score_gap, dplyr::all_of(score_cols))
  write.csv(cell_annotations, file.path(output_dir, "cell_annotations.csv"), row.names = TRUE)

  # Save session info for reproducibility
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

  # Save structured session info if devtools is available
  if(requireNamespace("devtools", quietly = TRUE)) {
    saveRDS(devtools::session_info(), file.path(output_dir, "session_info.rds"))
  }

  log_info(paste0("Cell annotation complete! Results saved to: ", output_dir))

  # Return a structured result
  return(list(
    seurat = seurat_obj,
    tables = list(
      counts = cluster_type_table,
      percent = cluster_type_pct
    ),
    scores = scores_df,
    markers = markers,
    parameters = list(
      ctrl_size = ctrl_size,
      score_threshold = score_threshold,
      seed = seed
    )
  ))
}

# Define brain-specific cell markers for mouse
brain_markers <- list(
  "Neurons" = c("Rbfox3", "Snap25", "Syt1", "Tubb3"),
  "Excitatory_neurons" = c("Slc17a7", "Camk2a", "Nrgn", "Satb2"),
  "Inhibitory_neurons" = c("Gad1", "Gad2", "Slc32a1", "Dlx1"),
  "Astrocytes" = c("Gfap", "Aldh1l1", "Aqp4", "Slc1a3"),
  "Oligodendrocytes" = c("Mbp", "Mag", "Mog", "Plp1"),
  "OPC" = c("Pdgfra", "Cspg4", "Sox10", "Olig1"),
  "Microglia" = c("Cx3cr1", "P2ry12", "Tmem119", "Csf1r"),
  "Endothelial" = c("Cldn5", "Pecam1", "Flt1", "Tek"),
  "Pericytes" = c("Pdgfrb", "Anpep", "Rgs5", "Notch3")
)

# Define paths
harmony_dir <- "/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE207128/harmony_results"
output_dir <- file.path(harmony_dir, "cell_annotation")

# Load Harmony integrated Seurat object
seurat_harmony <- readRDS(file.path(harmony_dir, "seurat_harmony_integrated.rds"))

# Check if gene IDs are in Ensembl format and convert markers if needed
gene_ids <- head(rownames(seurat_harmony), 5)
cat("Gene ID format in dataset:", gene_ids, "\n")

if(any(grepl("^ENSMUS", gene_ids))) {
  cat("Dataset uses Ensembl IDs. Converting markers from gene symbols...\n")

  # Set up biomaRt connection
  mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

  # Collect all unique markers
  all_markers <- unique(unlist(brain_markers))

  # Get mappings for all markers at once
  mapping <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "external_gene_name",
    values = all_markers,
    mart = mart
  )

  # Convert brain_markers list to use Ensembl IDs
  brain_markers_ensembl <- lapply(brain_markers, function(genes) {
    ensembl_ids <- mapping$ensembl_gene_id[mapping$external_gene_name %in% genes]
    if(length(ensembl_ids) == 0) {
      cat("No Ensembl IDs found for cell type markers:", paste(genes, collapse=", "), "\n")
    } else {
      cat("Found", length(ensembl_ids), "Ensembl IDs for markers\n")
    }
    return(ensembl_ids)
  })

  # Replace original markers with Ensembl IDs
  brain_markers <- brain_markers_ensembl
}

# Run the annotation pipeline
annotation_results <- run_annotation(
  seurat_obj = seurat_harmony,
  markers = brain_markers,
  output_dir = output_dir,
  ctrl_size = 50,
  score_threshold = 0.1,
  seed = 123,
  overwrite = TRUE  # Set to TRUE to replace any existing score columns
)

# The results are now available as a structured list:
# - annotation_results$seurat: Annotated Seurat object
# - annotation_results$tables$counts: Raw count table of cell types per cluster
# - annotation_results$tables$percent: Percentage composition table
# - annotation_results$scores: Raw score matrix for all markers