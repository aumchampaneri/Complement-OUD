# ==============================================================================
# Advanced Visualizations for Neuroinflammatory OUD Analysis
# ==============================================================================
# Purpose: Create publication-ready advanced plots and interactive visualizations
# Features: Chord diagrams, sankey plots, volcano plots, heatmaps, network graphs
# Input: Results from comprehensive neuroinflammatory analysis
# Output: High-impact visualizations for publication
# ==============================================================================

# Configuration
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens"
NEUROINFLAMM_DIR <- file.path(BASE_DIR, "Results", "Neuroinflammatory_Analysis")
OUTPUTS_BASE <- file.path(BASE_DIR, "Outputs", "Neuroinflammatory_Analysis")
VIZ_OUTPUTS <- file.path(OUTPUTS_BASE, "Advanced_Visualizations")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(plotly)
  library(circlize)      # Chord diagrams
  library(networkD3)     # Sankey plots
  library(ComplexHeatmap) # Advanced heatmaps
  library(ggraph)        # Network graphs
  library(igraph)
  library(viridis)       # Color palettes
  library(patchwork)     # Plot composition
  library(ggridges)      # Ridge plots
  library(ggbeeswarm)    # Beeswarm plots
  library(tidyr)
  library(RColorBrewer)
  library(cowplot)
  library(pheatmap)
  library(tibble)        # Add this for column_to_rownames and rownames_to_column
})

# ==============================================================================
# 1. INTERACTIVE VOLCANO PLOTS (FIXED)
# ==============================================================================

#' Create interactive volcano plots with hover information (with fallback)
create_interactive_volcano_plots <- function(comprehensive_results) {
  cat("\n=== Creating Volcano Plots ===\n")
  
  volcano_dir <- file.path(VIZ_OUTPUTS, "Volcano_Plots")
  if (!dir.exists(volcano_dir)) dir.create(volcano_dir, recursive = TRUE)
  
  for (dataset in names(comprehensive_results$enrichment_results)) {
    for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
      
      # Get results if they exist
      if ("GO_BP" %in% names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
        go_results <- comprehensive_results$enrichment_results[[dataset]][[method]]$GO_BP
        
        if (!is.null(go_results) && nrow(go_results@result) > 0) {
          pathway_data <- go_results@result
          
          # Prepare volcano data
          volcano_data <- pathway_data %>%
            mutate(
              neg_log10_pval = -log10(pvalue),
              neg_log10_padj = -log10(p.adjust),
              is_significant = p.adjust < 0.05,
              size_by_count = as.numeric(Count)
            ) %>%
            arrange(desc(is_significant), desc(neg_log10_padj))
          
          # Try interactive plot first, fallback to static
          tryCatch({
            # Create interactive volcano plot
            p_volcano <- plot_ly(
              data = volcano_data,
              x = ~size_by_count,
              y = ~neg_log10_padj,
              color = ~is_significant,
              colors = c("FALSE" = "lightgray", "TRUE" = "red"),
              size = ~size_by_count,
              sizes = c(20, 200),
              text = ~paste("Pathway:", Description, "<br>Count:", Count),
              hovertemplate = "%{text}<extra></extra>",
              type = "scatter",
              mode = "markers"
            ) %>%
            layout(
              title = paste("Interactive Volcano Plot:", dataset, "-", method),
              xaxis = list(title = "Gene Count"),
              yaxis = list(title = "-log10(Adjusted P-value)"),
              showlegend = TRUE
            )
            
            # Save interactive plot (without selfcontained to avoid pandoc issue)
            volcano_file <- file.path(volcano_dir, paste0(dataset, "_", method, "_volcano.html"))
            htmlwidgets::saveWidget(p_volcano, volcano_file, selfcontained = FALSE)
            cat("✓ Interactive volcano saved:", volcano_file, "\n")
            
          }, error = function(e) {
            cat("⚠ Interactive plot failed, creating static version:", e$message, "\n")
            
            # Create static volcano plot as fallback
            p_static <- ggplot(volcano_data, aes(x = size_by_count, y = neg_log10_padj)) +
              geom_point(aes(color = is_significant, size = size_by_count), alpha = 0.7) +
              scale_color_manual(values = c("FALSE" = "lightgray", "TRUE" = "red"), 
                               name = "Significant") +
              scale_size_continuous(name = "Gene Count", range = c(1, 8)) +
              labs(
                title = paste("Volcano Plot:", dataset, "-", method),
                x = "Gene Count",
                y = "-log10(Adjusted P-value)"
              ) +
              theme_minimal()
            
            static_file <- file.path(volcano_dir, paste0(dataset, "_", method, "_volcano_static.png"))
            ggsave(static_file, p_static, width = 10, height = 8, dpi = 300)
            cat("✓ Static volcano saved:", static_file, "\n")
          })
        }
      }
    }
  }
}

# ==============================================================================
# 2. PATHWAY CHORD DIAGRAMS (FIXED)
# ==============================================================================

#' Create chord diagrams showing pathway-dataset relationships
create_pathway_chord_diagrams <- function(comprehensive_results) {
  cat("\n=== Creating Pathway Chord Diagrams ===\n")
  
  chord_dir <- file.path(VIZ_OUTPUTS, "Chord_Diagrams")
  if (!dir.exists(chord_dir)) dir.create(chord_dir, recursive = TRUE)
  
  # Extract pathway data for chord diagram
  pathway_data <- data.frame()
  
  for (dataset in names(comprehensive_results$enrichment_results)) {
    for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
      for (db in names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
        result_obj <- comprehensive_results$enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          sig_pathways <- result_obj@result[result_obj@result$p.adjust < 0.05, ]
          
          if (nrow(sig_pathways) > 0) {
            pathway_data <- rbind(pathway_data, data.frame(
              Dataset = dataset,
              Database = db,
              Pathway = sig_pathways$Description,
              Count = as.numeric(sig_pathways$Count),
              NegLogP = -log10(sig_pathways$p.adjust),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  if (nrow(pathway_data) > 0) {
    # Get top pathways for chord diagram
    top_pathways <- pathway_data %>%
      group_by(Pathway) %>%
      summarise(
        total_count = sum(Count),
        max_neglogp = max(NegLogP),
        n_datasets = n_distinct(Dataset),
        .groups = 'drop'
      ) %>%
      arrange(desc(max_neglogp)) %>%
      slice_head(n = 20)
    
    # Create matrix for chord diagram
    chord_data <- pathway_data %>%
      filter(Pathway %in% top_pathways$Pathway) %>%
      group_by(Dataset, Database) %>%
      summarise(total_pathways = n(), .groups = 'drop')
    
    # Convert to matrix using base R approach
    chord_matrix <- chord_data %>%
      pivot_wider(names_from = Database, values_from = total_pathways, values_fill = 0)
    
    # Set row names manually
    dataset_names <- chord_matrix$Dataset
    chord_matrix$Dataset <- NULL
    chord_matrix <- as.matrix(chord_matrix)
    rownames(chord_matrix) <- dataset_names
    
    # Only create chord diagram if we have enough data
    if (nrow(chord_matrix) >= 2 && ncol(chord_matrix) >= 2) {
      # Create chord diagram
      png(file.path(chord_dir, "pathway_chord_diagram.png"), 
          width = 12, height = 12, units = "in", res = 300)
      
      # Set colors
      grid_colors <- rainbow(nrow(chord_matrix) + ncol(chord_matrix))
      names(grid_colors) <- c(rownames(chord_matrix), colnames(chord_matrix))
      
      # Create chord diagram
      circos.clear()
      circos.par(start.degree = 90, clock.wise = FALSE)
      chordDiagram(
        chord_matrix,
        grid.col = grid_colors,
        transparency = 0.3,
        annotationTrack = "grid",
        preAllocateTracks = 1
      )
      
      # Add labels
      circos.trackPlotRegion(
        track.index = 1, 
        panel.fun = function(x, y) {
          xlim = get.cell.meta.data("xlim")
          ylim = get.cell.meta.data("ylim")
          sector.name = get.cell.meta.data("sector.index")
          circos.text(mean(xlim), ylim[1] + 0.1, sector.name, 
                     facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
        },
        bg.border = NA
      )
      
      title("Pathway Database Relationships Across Datasets")
      dev.off()
      circos.clear()
      
      cat("✓ Chord diagram saved\n")
    } else {
      cat("⚠ Insufficient data for chord diagram\n")
    }
  }
}

# ==============================================================================
# 3. SANKEY FLOW DIAGRAMS (FIXED)
# ==============================================================================

#' Create Sankey diagrams showing pathway flow across methods (with fallback)
create_pathway_sankey_diagrams <- function(comprehensive_results) {
  cat("\n=== Creating Pathway Flow Diagrams ===\n")
  
  sankey_dir <- file.path(VIZ_OUTPUTS, "Sankey_Diagrams")
  if (!dir.exists(sankey_dir)) dir.create(sankey_dir, recursive = TRUE)
  
  for (dataset in names(comprehensive_results$enrichment_results)) {
    # Extract pathway overlaps between methods
    method_pathways <- list()
    
    for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
      if ("GO_BP" %in% names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
        go_results <- comprehensive_results$enrichment_results[[dataset]][[method]]$GO_BP
        
        if (!is.null(go_results) && nrow(go_results@result) > 0) {
          sig_pathways <- go_results@result[go_results@result$p.adjust < 0.05, ]
          method_pathways[[method]] <- sig_pathways$Description[1:min(20, nrow(sig_pathways))]
        }
      }
    }
    
    if (length(method_pathways) >= 2) {
      # Create simple overlap visualization instead of complex Sankey
      overlap_data <- data.frame()
      methods <- names(method_pathways)
      
      for (i in 1:(length(methods)-1)) {
        for (j in (i+1):length(methods)) {
          method1 <- methods[i]
          method2 <- methods[j]
          
          # Find common pathways
          common_pathways <- intersect(method_pathways[[method1]], method_pathways[[method2]])
          
          overlap_data <- rbind(overlap_data, data.frame(
            Method1 = method1,
            Method2 = method2,
            Overlap = length(common_pathways),
            Method1_Total = length(method_pathways[[method1]]),
            Method2_Total = length(method_pathways[[method2]])
          ))
        }
      }
      
      if (nrow(overlap_data) > 0) {
        # Create overlap bar plot
        p_overlap <- ggplot(overlap_data, aes(x = paste(Method1, "vs", Method2), y = Overlap)) +
          geom_col(fill = "steelblue", alpha = 0.7) +
          geom_text(aes(label = Overlap), vjust = -0.5) +
          labs(
            title = paste("Pathway Overlap Between Methods -", dataset),
            x = "Method Comparison",
            y = "Number of Shared Pathways"
          ) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        overlap_file <- file.path(sankey_dir, paste0(dataset, "_method_overlap.png"))
        ggsave(overlap_file, p_overlap, width = 10, height = 6, dpi = 300)
        cat("✓ Method overlap plot saved:", overlap_file, "\n")
      }
    }
  }
}

# ==============================================================================
# 4. ADVANCED HEATMAPS WITH CLUSTERING (FIXED)
# ==============================================================================

#' Create advanced heatmaps with hierarchical clustering (with fallback)
create_advanced_heatmaps <- function(comprehensive_results) {
  cat("\n=== Creating Advanced Heatmaps ===\n")
  
  heatmap_dir <- file.path(VIZ_OUTPUTS, "Advanced_Heatmaps")
  if (!dir.exists(heatmap_dir)) dir.create(heatmap_dir, recursive = TRUE)
  
  # Collect pathway data
  pathway_matrix_data <- data.frame()
  
  for (dataset in names(comprehensive_results$enrichment_results)) {
    for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
      for (db in names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
        result_obj <- comprehensive_results$enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          sig_pathways <- result_obj@result[result_obj@result$p.adjust < 0.05, ]
          
          if (nrow(sig_pathways) > 0) {
            pathway_matrix_data <- rbind(pathway_matrix_data, data.frame(
              Dataset = dataset,
              Method = method,
              Database = db,
              Pathway = sig_pathways$Description,
              NegLogP = -log10(sig_pathways$p.adjust),
              Count = as.numeric(sig_pathways$Count),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  if (nrow(pathway_matrix_data) > 0) {
    # Create pathway significance matrix
    top_pathways <- pathway_matrix_data %>%
      group_by(Pathway) %>%
      summarise(max_neglogp = max(NegLogP), .groups = 'drop') %>%
      arrange(desc(max_neglogp)) %>%
      slice_head(n = 50)
    
    heatmap_data <- pathway_matrix_data %>%
      filter(Pathway %in% top_pathways$Pathway) %>%
      unite("Analysis", Dataset, Method, Database, sep = "_") %>%
      dplyr::select(Analysis, Pathway, NegLogP) %>%  # Use dplyr:: explicitly
      pivot_wider(names_from = Analysis, values_from = NegLogP, values_fill = 0)
    
    # Set row names manually using base R
    pathway_names <- heatmap_data$Pathway
    heatmap_data$Pathway <- NULL
    heatmap_data <- as.matrix(heatmap_data)
    rownames(heatmap_data) <- pathway_names
    
    # Only proceed if we have sufficient data
    if (nrow(heatmap_data) >= 5 && ncol(heatmap_data) >= 2) {
      # Try ComplexHeatmap first, fallback to pheatmap
      tryCatch({
        # Create advanced heatmap with ComplexHeatmap
        png(file.path(heatmap_dir, "pathway_significance_heatmap_complex.png"),
            width = 16, height = 12, units = "in", res = 300)
        
        # Color function
        col_fun <- colorRamp2(c(0, 2, 5, 10), c("white", "yellow", "orange", "red"))
        
        # Create heatmap
        ht <- Heatmap(
          heatmap_data,
          name = "-log10(p.adj)",
          col = col_fun,
          clustering_distance_rows = "euclidean",
          clustering_distance_columns = "euclidean",
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 10),
          column_names_rot = 45,
          show_row_names = TRUE,
          show_column_names = TRUE,
          heatmap_legend_param = list(title = "-log10(Adj P-value)")
        )
        
        draw(ht)
        dev.off()
        cat("✓ Complex heatmap saved\n")
        
      }, error = function(e) {
        cat("⚠ ComplexHeatmap failed, using pheatmap:", e$message, "\n")
        
        # Fallback to pheatmap
        png(file.path(heatmap_dir, "pathway_significance_heatmap_simple.png"),
            width = 16, height = 12, units = "in", res = 300)
        
        pheatmap(
          heatmap_data,
          color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
          clustering_distance_rows = "euclidean",
          clustering_distance_cols = "euclidean",
          fontsize_row = 8,
          fontsize_col = 10,
          angle_col = 45,
          main = "Pathway Significance Heatmap"
        )
        
        dev.off()
        cat("✓ Simple heatmap saved\n")
      })
    } else {
      cat("⚠ Insufficient data for heatmap (need ≥5 pathways, ≥2 analyses)\n")
    }
  }
}

# ==============================================================================
# 5. NETWORK GRAPH VISUALIZATIONS
# ==============================================================================

#' Create pathway network graphs with ggraph
create_pathway_network_graphs <- function(comprehensive_results) {
  cat("\n=== Creating Pathway Network Graphs ===\n")
  
  network_dir <- file.path(VIZ_OUTPUTS, "Network_Graphs")
  if (!dir.exists(network_dir)) dir.create(network_dir, recursive = TRUE)
  
  # Extract pathway relationships
  pathway_data <- data.frame()
  
  for (dataset in names(comprehensive_results$enrichment_results)) {
    for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
      if ("GO_BP" %in% names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
        go_results <- comprehensive_results$enrichment_results[[dataset]][[method]]$GO_BP
        
        if (!is.null(go_results) && nrow(go_results@result) > 0) {
          sig_pathways <- go_results@result[go_results@result$p.adjust < 0.05, ]
          
          if (nrow(sig_pathways) > 0) {
            pathway_data <- rbind(pathway_data, data.frame(
              Dataset = dataset,
              Method = method,
              Pathway = sig_pathways$Description,
              NegLogP = -log10(sig_pathways$p.adjust),
              Count = as.numeric(sig_pathways$Count),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  if (nrow(pathway_data) > 0) {
    # Create network edges based on shared pathways
    top_pathways <- pathway_data %>%
      group_by(Pathway) %>%
      summarise(
        datasets = paste(unique(Dataset), collapse = ","),
        max_neglogp = max(NegLogP),
        .groups = 'drop'
      ) %>%
      arrange(desc(max_neglogp)) %>%
      slice_head(n = 30)
    
    # Create nodes and edges
    nodes <- data.frame(
      id = 1:nrow(top_pathways),
      label = top_pathways$Pathway,
      significance = top_pathways$max_neglogp,
      size = log10(top_pathways$max_neglogp + 1) * 3
    )
    
    # Create edges based on pathway similarity (simplified)
    edges <- data.frame()
    for (i in 1:(nrow(nodes)-1)) {
      for (j in (i+1):nrow(nodes)) {
        # Simple similarity based on shared words
        pathway1_words <- tolower(strsplit(nodes$label[i], " ")[[1]])
        pathway2_words <- tolower(strsplit(nodes$label[j], " ")[[1]])
        shared_words <- intersect(pathway1_words, pathway2_words)
        
        if (length(shared_words) >= 2) {  # At least 2 shared words
          edges <- rbind(edges, data.frame(
            from = i,
            to = j,
            weight = length(shared_words)
          ))
        }
      }
    }
    
    if (nrow(edges) > 0) {
      # Create igraph object
      g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
      
      # Create network plot with ggraph
      p_network <- ggraph(g, layout = "fr") +
        geom_edge_link(aes(alpha = weight), color = "gray60", show.legend = FALSE) +
        geom_node_point(aes(size = significance, color = significance)) +
        geom_node_text(aes(label = label), size = 2.5, repel = TRUE, max.overlaps = 20) +
        scale_color_viridis_c(name = "-log10(p.adj)") +
        scale_size_continuous(name = "Significance", range = c(3, 10)) +
        theme_graph() +
        labs(title = "Pathway Network Graph",
             subtitle = "Nodes: Pathways, Edges: Shared terminology")
        
      ggsave(file.path(network_dir, "pathway_network_graph.png"),
             p_network, width = 16, height = 12, dpi = 300)
      
      cat("✓ Network graph saved\n")
    }
  }
}

# ==============================================================================
# 6. RIDGE PLOTS FOR EXPRESSION DISTRIBUTIONS (ENHANCED)
# ==============================================================================

#' Create ridge plots showing expression distributions (enhanced fallback)
create_expression_ridge_plots <- function(comprehensive_results) {
  cat("\n=== Creating Expression Ridge Plots ===\n")
  
  ridge_dir <- file.path(VIZ_OUTPUTS, "Ridge_Plots")
  if (!dir.exists(ridge_dir)) dir.create(ridge_dir, recursive = TRUE)
  
  # Check for expression data in comprehensive results
  if (!is.null(comprehensive_results$expression_data)) {
    for (dataset_name in names(comprehensive_results$expression_data)) {
      expr_data <- comprehensive_results$expression_data[[dataset_name]]
      
      if (!is.null(expr_data$log_cpm)) {
        # Sample a subset of genes for visualization
        n_genes_to_plot <- 20
        gene_variances <- apply(expr_data$log_cpm, 1, var, na.rm = TRUE)
        top_var_genes <- names(sort(gene_variances, decreasing = TRUE))[1:n_genes_to_plot]
        
        # Prepare data for ridge plot using base R approach
        ridge_data <- expr_data$log_cpm[top_var_genes, ] %>%
          as.data.frame()
        
        # Add gene names manually
        ridge_data$Gene <- rownames(ridge_data)
        
        # Reshape to long format
        ridge_data <- ridge_data %>%
          pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")
        
        # Add metadata - check column names
        metadata_cols <- colnames(expr_data$metadata)
        if ("Sample_ID" %in% metadata_cols) {
          sample_col <- "Sample_ID"
        } else if ("Sample" %in% metadata_cols) {
          sample_col <- "Sample"
        } else {
          sample_col <- colnames(expr_data$metadata)[1]  # Use first column as fallback
        }
        
        if ("Region" %in% metadata_cols) {
          region_col <- "Region"
        } else if ("Brain_Region" %in% metadata_cols) {
          region_col <- "Brain_Region"
        } else {
          # Find a column that might be region-related
          region_col <- metadata_cols[grepl("region", tolower(metadata_cols))][1]
          if (is.na(region_col)) {
            region_col <- metadata_cols[2]  # Use second column as fallback
          }
        }
        
        # Create metadata mapping
        metadata_subset <- expr_data$metadata[, c(sample_col, region_col), drop = FALSE]
        colnames(metadata_subset) <- c("Sample", "Region")
        
        # Join with ridge data
        ridge_data <- ridge_data %>%
          left_join(metadata_subset, by = "Sample")
        
        # Only create plot if we have Region data
        if (!"Region" %in% colnames(ridge_data) || all(is.na(ridge_data$Region))) {
          cat("⚠ No region information found for", dataset_name, "\n")
          next
        }
        
        # Create ridge plot
        p_ridge <- ggplot(ridge_data, aes(x = Expression, y = Gene, fill = Region)) +
          geom_density_ridges(alpha = 0.7, scale = 0.9) +
          scale_fill_viridis_d() +
          labs(
            title = paste("Expression Distribution Ridge Plot -", dataset_name),
            subtitle = "Top 20 most variable genes",
            x = "Log2-CPM Expression",
            y = "Gene"
          ) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 8))
        
        ggsave(file.path(ridge_dir, paste0(dataset_name, "_expression_ridges.png")),
               p_ridge, width = 12, height = 10, dpi = 300)
        
        cat("✓ Ridge plot saved for", dataset_name, "\n")
      }
    }
  } else {
    cat("⚠ No expression data found in comprehensive results\n")
    cat("💡 Creating pathway distribution plot instead...\n")
    
    # Alternative: Create pathway distribution plot
    tryCatch({
      # Extract pathway counts by dataset and method
      pathway_counts <- data.frame()
      
      for (dataset in names(comprehensive_results$enrichment_results)) {
        for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
          for (db in names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
            result_obj <- comprehensive_results$enrichment_results[[dataset]][[method]][[db]]
            
            if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
              sig_pathways <- result_obj@result[result_obj@result$p.adjust < 0.05, ]
              
              if (nrow(sig_pathways) > 0) {
                pathway_counts <- rbind(pathway_counts, data.frame(
                  Dataset = dataset,
                  Method = method,
                  Database = db,
                  Count = nrow(sig_pathways),
                  NegLogP_Mean = mean(-log10(sig_pathways$p.adjust), na.rm = TRUE)
                ))
              }
            }
          }
        }
      }
      
      if (nrow(pathway_counts) > 0) {
        # Create pathway distribution plot
        p_pathway_dist <- ggplot(pathway_counts, aes(x = Dataset, y = Count, fill = Method)) +
          geom_col(position = "dodge", alpha = 0.8) +
          facet_wrap(~Database, scales = "free_y") +
          scale_fill_viridis_d() +
          labs(
            title = "Pathway Count Distribution",
            subtitle = "Significant pathways by dataset, method, and database",
            x = "Dataset",
            y = "Number of Significant Pathways"
          ) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggsave(file.path(ridge_dir, "pathway_distribution_alternative.png"),
               p_pathway_dist, width = 14, height = 10, dpi = 300)
        
        cat("✓ Alternative pathway distribution plot saved\n")
      }
      
    }, error = function(e) {
      cat("⚠ Alternative plot also failed:", e$message, "\n")
    })
  }
}

# ==============================================================================
# 7. COMBINED DASHBOARD PLOT (FIXED)
# ==============================================================================

#' Create a combined dashboard with multiple visualizations
create_combined_dashboard <- function(comprehensive_results) {
  cat("\n=== Creating Combined Dashboard ===\n")
  
  dashboard_dir <- file.path(VIZ_OUTPUTS, "Combined_Dashboard")
  if (!dir.exists(dashboard_dir)) dir.create(dashboard_dir, recursive = TRUE)
  
  # Create multiple plot components
  plots <- list()
  
  # 1. Method comparison bar plot
  summary_data <- data.frame()
  for (dataset in names(comprehensive_results$enrichment_results)) {
    for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
      total_sig <- 0
      for (db in names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
        result_obj <- comprehensive_results$enrichment_results[[dataset]][[method]][[db]]
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          total_sig <- total_sig + sum(result_obj@result$p.adjust < 0.05, na.rm = TRUE)
        }
      }
      summary_data <- rbind(summary_data, data.frame(
        Dataset = dataset,
        Method = method,
        Significant_Pathways = total_sig
      ))
    }
  }
  
  p1 <- ggplot(summary_data, aes(x = Method, y = Significant_Pathways, fill = Dataset)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_viridis_d() +
    labs(title = "A. Significant Pathways by Method", x = "Method", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 2. Database comparison pie chart
  db_summary <- data.frame()
  for (dataset in names(comprehensive_results$enrichment_results)) {
    for (method in names(comprehensive_results$enrichment_results[[dataset]])) {
      for (db in names(comprehensive_results$enrichment_results[[dataset]][[method]])) {
        result_obj <- comprehensive_results$enrichment_results[[dataset]][[method]][[db]]
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          sig_count <- sum(result_obj@result$p.adjust < 0.05, na.rm = TRUE)
          db_summary <- rbind(db_summary, data.frame(
            Database = db,
            Count = sig_count
          ))
        }
      }
    }
  }
  
  db_totals <- db_summary %>%
    group_by(Database) %>%
    summarise(Total = sum(Count), .groups = 'drop')
  
  p2 <- ggplot(db_totals, aes(x = "", y = Total, fill = Database)) +
    geom_col() +
    coord_polar("y", start = 0) +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    labs(title = "B. Pathways by Database") +
    theme_void()
  
  # 3. Technology comparison (FIXED)
  tech_data <- summary_data %>%
    mutate(Technology = ifelse(grepl("RNA|microarray", Method, ignore.case = TRUE), "RNA-Seq/microarray",
                            ifelse(grepl("Methylation", Method, ignore.case = TRUE), "Methylation",
                                   ifelse(grepl("Proteomics", Method, ignore.case = TRUE), "Proteomics",
                                          "Other")))) %>%
    group_by(Technology) %>%
    summarise(Total_Significant = sum(Significant_Pathways), .groups = 'drop')
  
  p3 <- ggplot(tech_data, aes(x = Technology, y = Total_Significant, fill = Technology)) +
    geom_col(show.legend = FALSE) +
    labs(title = "C. Total Significant Pathways by Technology", x = "Technology", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 4. Combined heatmap (top 10 pathways per dataset)
  combined_heatmap_data <- pathway_matrix_data %>%
    group_by(Pathway) %>%
    summarise(mean_neglogp = mean(NegLogP), .groups = 'drop') %>%
    arrange(desc(mean_neglogp)) %>%
    slice_head(n = 10) %>%
    inner_join(pathway_matrix_data, by = "Pathway")
  
  if (nrow(combined_heatmap_data) > 0) {
    p4 <- ggplot(combined_heatmap_data, aes(x = reorder(Pathway, -mean_neglogp), y = Dataset_Method_Database)) +
      geom_tile(aes(fill = NegLogP), color = "white") +
      scale_fill_viridis_c(option = "C") +
      labs(title = "D. Heatmap of Top 10 Pathways (Combined)", x = "Pathway", y = "Dataset_Method_Database") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = 8))
  } else {
    p4 <- NULL
  }
  
  # 5. Expression boxplots (if data available)
  if (!is.null(comprehensive_results$expression_data)) {
    expr_boxplot_data <- comprehensive_results$expression_data %>%
      bind_rows(.id = "Dataset") %>%
      pivot_longer(-c(Dataset, Sample), names_to = "Gene", values_to = "Expression") %>%
      group_by(Dataset, Gene) %>%
      summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = 'drop')
    
    top_expr_genes <- expr_boxplot_data %>%
      group_by(Gene) %>%
      summarise(Mean_Expression = mean(Mean_Expression), .groups = 'drop') %>%
      arrange(desc(Mean_Expression)) %>%
      slice_head(n = 10) %>%
      pull(Gene)
    
    p5 <- ggplot(expr_boxplot_data %>% filter(Gene %in% top_expr_genes), aes(x = Gene, y = Mean_Expression, fill = Dataset)) +
      geom_col(position = "dodge") +
      labs(title = "E. Mean Expression of Top Genes", x = "Gene", y = "Mean Expression") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    p5 <- NULL
  }
  
  # Combine all plots into a dashboard layout
  all_plots <- (p1 + p2 + p3 + p4 + p5) + plot_layout(ncol = 2)
  
  # Save dashboard
  ggsave(file.path(dashboard_dir, "combined_dashboard.png"), all_plots, width = 16, height = 12, dpi = 300)
  cat("✓ Combined dashboard saved\n")
}
