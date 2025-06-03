# ==============================================================================
# Network Integration and Technology Comparison Analysis
# ==============================================================================
# Purpose: Advanced network analysis and cross-technology validation
# Combines: TF networks, pathway networks, technology comparisons, meta-analysis
# Input: Results from 02_Neuroinflammatory_Analysis.R
# Output: Network plots and integration analysis
# ==============================================================================

# Configuration
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens"
NEUROINFLAMM_DIR <- file.path(BASE_DIR, "Results", "Neuroinflammatory_Analysis")
OUTPUTS_BASE <- file.path(BASE_DIR, "Outputs", "Neuroinflammatory_Analysis")

# Output directories
NETWORK_OUTPUTS <- file.path(OUTPUTS_BASE, "Network_Analysis")
PATHWAY_OUTPUTS <- file.path(OUTPUTS_BASE, "Pathway_Enrichment")
DECOUPLER_OUTPUTS <- file.path(OUTPUTS_BASE, "TF_Analysis")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(igraph)
  library(corrplot)
  library(tidyr)
  library(RColorBrewer)
})

# ==============================================================================
# ENHANCED TF ACTIVITY NETWORK ANALYSIS
# ==============================================================================

#' Create TF activity correlation networks
create_tf_correlation_networks <- function(tf_results_list) {
  cat("\n=== Creating TF Activity Correlation Networks ===\n")
  
  # Create TF output directory
  tf_network_dir <- file.path(NETWORK_OUTPUTS, "TF_Networks")
  if (!dir.exists(tf_network_dir)) dir.create(tf_network_dir, recursive = TRUE)
  
  for (dataset in names(tf_results_list)) {
    tryCatch({
      tf_data <- tf_results_list[[dataset]]
      
      if (!is.null(tf_data) && nrow(tf_data) > 5) {
        cat(sprintf("Processing TF correlations for %s...\n", dataset))
        
        # Create TF correlation matrix
        tf_matrix <- tf_data %>%
          select_if(is.numeric) %>%
          as.matrix()
        
        if (ncol(tf_matrix) > 2) {
          # Calculate correlations
          tf_cor <- cor(tf_matrix, use = "complete.obs")
          
          # Create correlation network plot
          png(file.path(tf_network_dir, paste0(dataset, "_TF_correlation_network.png")),
              width = 10, height = 8, units = "in", res = 300)
          
          corrplot(tf_cor, 
                   method = "circle",
                   type = "upper",
                   order = "hclust",
                   tl.cex = 0.8,
                   tl.col = "black",
                   title = paste("TF Activity Correlations:", dataset),
                   mar = c(0,0,2,0))
          
          dev.off()
          cat("✓ TF correlation network saved for", dataset, "\n")
        } else {
          cat("⚠ Insufficient TF data for correlation analysis in", dataset, "\n")
        }
      } else {
        cat("⚠ No TF data available for", dataset, "\n")
      }
      
    }, error = function(e) {
      cat("✗ TF network analysis failed for", dataset, ":", e$message, "\n")
    })
  }
}

# ==============================================================================
# PATHWAY INTERACTION NETWORKS
# ==============================================================================

#' Create pathway interaction networks
create_pathway_interaction_networks <- function(enrichment_results) {
  cat("\n=== Creating Pathway Interaction Networks ===\n")
  
  pathway_network_dir <- file.path(NETWORK_OUTPUTS, "Pathway_Networks")
  if (!dir.exists(pathway_network_dir)) dir.create(pathway_network_dir, recursive = TRUE)
  
  # Extract all significant pathways
  pathway_data <- data.frame()
  
  for (dataset in names(enrichment_results)) {
    for (method in names(enrichment_results[[dataset]])) {
      for (db in names(enrichment_results[[dataset]][[method]])) {
        result_obj <- enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          sig_pathways <- result_obj@result[result_obj@result$p.adjust < 0.05, ]
          
          if (nrow(sig_pathways) > 0) {
            pathway_data <- rbind(pathway_data, data.frame(
              Dataset = dataset,
              Method = method,
              Database = db,
              Pathway = sig_pathways$Description,
              NegLogP = -log10(sig_pathways$p.adjust),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  if (nrow(pathway_data) > 0) {
    tryCatch({
      # Create pathway co-occurrence network - FIX: Use proper data aggregation
      pathway_summary <- pathway_data %>%
        group_by(Pathway, Dataset) %>%
        summarise(AvgNegLogP = mean(NegLogP, na.rm = TRUE), .groups = 'drop')
      
      # Check if we have enough data for matrix creation
      if (nrow(pathway_summary) < 5) {
        cat("⚠ Insufficient pathway data for network analysis (", nrow(pathway_summary), "rows)\n")
        return(FALSE)
      }
      
      # Create pathway matrix using base R approach
      unique_pathways <- unique(pathway_summary$Pathway)
      unique_datasets <- unique(pathway_summary$Dataset)
      
      pathway_matrix <- matrix(0, 
                              nrow = length(unique_pathways), 
                              ncol = length(unique_datasets),
                              dimnames = list(unique_pathways, unique_datasets))
      
      # Fill matrix manually
      for (i in 1:nrow(pathway_summary)) {
        pathway <- pathway_summary$Pathway[i]
        dataset <- pathway_summary$Dataset[i]
        value <- pathway_summary$AvgNegLogP[i]
        pathway_matrix[pathway, dataset] <- value
      }
      
      # Only proceed if we have multiple datasets
      if (ncol(pathway_matrix) < 2) {
        cat("⚠ Need multiple datasets for pathway correlation analysis\n")
        
        # Create simple pathway frequency plot instead
        pathway_counts <- pathway_data %>%
          count(Pathway, sort = TRUE) %>%
          slice_head(n = 20)
        
        p_freq <- ggplot(pathway_counts, aes(x = reorder(Pathway, n), y = n)) +
          geom_col(fill = "steelblue", alpha = 0.7) +
          coord_flip() +
          labs(title = "Most Frequent Significant Pathways",
               x = "Pathway", y = "Frequency") +
          theme_minimal()
        
        ggsave(file.path(pathway_network_dir, "pathway_frequency_plot.png"),
               p_freq, width = 12, height = 8, dpi = 300)
        
        cat("✓ Pathway frequency plot saved\n")
        return(TRUE)
      }
      
      # Calculate pathway similarities (transpose for pathway-pathway correlation)
      pathway_cor <- cor(t(pathway_matrix), use = "complete.obs")
      
      # Create network from high correlations
      threshold <- 0.5
      pathway_cor[abs(pathway_cor) < threshold] <- 0
      diag(pathway_cor) <- 0
      
      # Check if we have any connections
      if (sum(abs(pathway_cor) > 0) == 0) {
        cat("⚠ No significant correlations found above threshold", threshold, "\n")
        
        # Lower threshold and try again
        threshold <- 0.3
        pathway_cor <- cor(t(pathway_matrix), use = "complete.obs")
        pathway_cor[abs(pathway_cor) < threshold] <- 0
        diag(pathway_cor) <- 0
      }
      
      if (sum(abs(pathway_cor) > 0) > 0) {
        # Convert to igraph - FIX: Use absolute values for weights to avoid negative weights
        pathway_graph <- graph_from_adjacency_matrix(abs(pathway_cor), 
                                                     mode = "undirected", 
                                                     weighted = TRUE)
        
        # Only plot if we have edges
        if (ecount(pathway_graph) > 0) {
          # Plot network with better layout handling
          png(file.path(pathway_network_dir, "pathway_interaction_network.png"),
              width = 12, height = 10, units = "in", res = 300)
          
          # Use a different layout that doesn't require positive weights
          tryCatch({
            # Try Fruchterman-Reingold with absolute weights (already done above)
            layout_coords <- layout_with_fr(pathway_graph)
          }, error = function(e) {
            # Fallback to a simpler layout
            cat("⚠ FR layout failed, using spring layout\n")
            layout_coords <- layout_with_dh(pathway_graph)
          })
          
          plot(pathway_graph,
               vertex.size = 8,
               vertex.label.cex = 0.6,
               vertex.label.color = "black",
               edge.width = E(pathway_graph)$weight * 3,
               edge.color = "gray60",
               layout = layout_coords,
               main = "Pathway Interaction Network\n(Based on Cross-Dataset Correlations)")
          
          dev.off()
          cat("✓ Pathway interaction network saved\n")
        } else {
          cat("⚠ No pathway connections to plot\n")
        }
      } else {
        cat("⚠ No significant pathway correlations found\n")
      }
      
    }, error = function(e) {
      cat("✗ Pathway network creation failed:", e$message, "\n")
      
      # Create fallback pathway frequency plot
      tryCatch({
        pathway_counts <- pathway_data %>%
          count(Pathway, sort = TRUE) %>%
          slice_head(n = 20)
        
        p_freq <- ggplot(pathway_counts, aes(x = reorder(Pathway, n), y = n)) +
          geom_col(fill = "steelblue", alpha = 0.7) +
          coord_flip() +
          labs(title = "Most Frequent Significant Pathways",
               subtitle = "Fallback visualization due to network creation failure",
               x = "Pathway", y = "Frequency") +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 8))
        
        ggsave(file.path(pathway_network_dir, "pathway_frequency_plot.png"),
               p_freq, width = 12, height = 8, dpi = 300)
        
        cat("✓ Fallback pathway frequency plot saved\n")
        
      }, error = function(e2) {
        cat("✗ Fallback visualization also failed:", e2$message, "\n")
      })
      
      return(FALSE)
    })
  } else {
    cat("⚠ No pathway data available for network analysis\n")
    return(FALSE)
  }
  
  return(TRUE)
}

# ==============================================================================
# TECHNOLOGY COMPARISON ANALYSIS
# ==============================================================================

#' Systematic comparison of bulk vs snRNA-seq findings
create_technology_comparison_analysis <- function(enrichment_results, expression_data = NULL) {
  cat("\n=== Technology Comparison Analysis ===\n")
  
  comparison_dir <- file.path(NETWORK_OUTPUTS, "Technology_Comparison")
  if (!dir.exists(comparison_dir)) dir.create(comparison_dir, recursive = TRUE)
  
  # 1. Pathway enrichment comparison
  tryCatch({
    pathway_summary <- data.frame()
    
    for (dataset in names(enrichment_results)) {
      for (method in names(enrichment_results[[dataset]])) {
        total_sig <- 0
        
        for (db in names(enrichment_results[[dataset]][[method]])) {
          result_obj <- enrichment_results[[dataset]][[method]][[db]]
          
          if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
            sig_count <- sum(result_obj@result$p.adjust < 0.05, na.rm = TRUE)
            total_sig <- total_sig + sig_count
          }
        }
        
        pathway_summary <- rbind(pathway_summary, data.frame(
          Dataset = dataset,
          Method = method,
          Technology = ifelse(grepl("174409", dataset), "Bulk RNA-seq", "snRNA-seq"),
          Significant_Pathways = total_sig,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    if (nrow(pathway_summary) > 0) {
      # Create comparison plot
      p_tech_comp <- ggplot(pathway_summary, aes(x = Method, y = Significant_Pathways, fill = Technology)) +
        geom_col(position = "dodge", alpha = 0.8, color = "black", size = 0.3) +
        geom_text(aes(label = Significant_Pathways), 
                  position = position_dodge(width = 0.9), 
                  vjust = -0.5, size = 3) +
        scale_fill_manual(values = c("Bulk RNA-seq" = "#FF6B6B", "snRNA-seq" = "#4ECDC4")) +
        labs(title = "Technology Comparison: Pathway Enrichment",
             subtitle = "Significant pathways across bulk RNA-seq vs snRNA-seq",
             x = "Analysis Method",
             y = "Number of Significant Pathways",
             fill = "Technology") +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
      
      ggsave(file.path(comparison_dir, "technology_pathway_comparison.png"),
             p_tech_comp, width = 10, height = 6, dpi = 300)
      
      cat("✓ Technology comparison plot saved\n")
      
      # Save summary table
      write.csv(pathway_summary, file.path(comparison_dir, "technology_comparison_summary.csv"), 
                row.names = FALSE)
      cat("✓ Technology comparison summary saved\n")
    }
    
  }, error = function(e) {
    cat("⚠ Technology comparison analysis failed:", e$message, "\n")
  })
}

# ==============================================================================
# CROSS-DATASET GENE OVERLAP ANALYSIS
# ==============================================================================

#' Analyze gene overlap patterns across datasets
analyze_cross_dataset_gene_overlap <- function(enrichment_results) {
  cat("\n=== Cross-Dataset Gene Overlap Analysis ===\n")
  
  overlap_dir <- file.path(NETWORK_OUTPUTS, "Gene_Overlap")
  if (!dir.exists(overlap_dir)) dir.create(overlap_dir, recursive = TRUE)
  
  # Extract all genes from pathway results
  dataset_genes <- list()
  
  for (dataset in names(enrichment_results)) {
    all_genes <- c()
    
    for (method in names(enrichment_results[[dataset]])) {
      for (db in names(enrichment_results[[dataset]][[method]])) {
        result_obj <- enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          # Extract genes from significant pathways
          sig_pathways <- result_obj@result[result_obj@result$p.adjust < 0.05, ]
          
          if (nrow(sig_pathways) > 0 && "geneID" %in% colnames(sig_pathways)) {
            pathway_genes <- unique(unlist(strsplit(sig_pathways$geneID, "/")))
            all_genes <- c(all_genes, pathway_genes)
          }
        }
      }
    }
    
    if (length(all_genes) > 0) {
      dataset_genes[[dataset]] <- unique(all_genes)
      cat(sprintf("%s: %d unique genes\n", dataset, length(unique(all_genes))))
    }
  }
  
  # Calculate overlap if we have data from multiple datasets
  if (length(dataset_genes) >= 2) {
    dataset_names <- names(dataset_genes)
    overlap_count <- length(intersect(dataset_genes[[1]], dataset_genes[[2]]))
    
    cat(sprintf("Gene overlap between datasets: %d genes\n", overlap_count))
    
    # Create overlap visualization
    if (overlap_count > 0) {
      overlap_data <- data.frame(
        Dataset1 = length(dataset_genes[[1]]),
        Dataset2 = length(dataset_genes[[2]]),
        Overlap = overlap_count,
        Dataset1_Name = dataset_names[1],
        Dataset2_Name = dataset_names[2]
      )
      
      write.csv(overlap_data, file.path(overlap_dir, "gene_overlap_summary.csv"), 
                row.names = FALSE)
      cat("✓ Gene overlap summary saved\n")
    }
  }
  
  return(dataset_genes)
}

# ==============================================================================
# NETWORK INTEGRATION SUMMARY
# ==============================================================================

#' Create comprehensive network integration summary
create_network_integration_summary <- function(enrichment_results, dataset_genes) {
  cat("\n=== Creating Network Integration Summary ===\n")
  
  summary_dir <- file.path(NETWORK_OUTPUTS, "Integration_Summary")
  if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)
  
  # Collect summary statistics
  summary_stats <- list()
  
  for (dataset in names(enrichment_results)) {
    dataset_summary <- list(
      dataset = dataset,
      methods = length(enrichment_results[[dataset]]),
      total_pathways = 0,
      significant_pathways = 0
    )
    
    for (method in names(enrichment_results[[dataset]])) {
      for (db in names(enrichment_results[[dataset]][[method]])) {
        result_obj <- enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          dataset_summary$total_pathways <- dataset_summary$total_pathways + nrow(result_obj@result)
          dataset_summary$significant_pathways <- dataset_summary$significant_pathways + 
            sum(result_obj@result$p.adjust < 0.05, na.rm = TRUE)
        }
      }
    }
    
    summary_stats[[dataset]] <- dataset_summary
  }
  
  # Convert to data frame
  summary_df <- do.call(rbind, lapply(summary_stats, data.frame))
  
  # Save summary
  write.csv(summary_df, file.path(summary_dir, "network_integration_summary.csv"), 
            row.names = FALSE)
  
  cat("✓ Network integration summary saved\n")
  cat("\n=== Integration Summary ===\n")
  print(summary_df)
  
  return(summary_df)
}

# ==============================================================================
# MAIN NETWORK INTEGRATION PIPELINE
# ==============================================================================

#' Run comprehensive network integration analysis
run_network_integration_analysis <- function() {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("NETWORK INTEGRATION AND TECHNOLOGY COMPARISON\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  tryCatch({
    # Create output directories
    dirs_to_create <- c(
      file.path(NETWORK_OUTPUTS, "TF_Networks"),
      file.path(NETWORK_OUTPUTS, "Pathway_Networks"),
      file.path(NETWORK_OUTPUTS, "Technology_Comparison"),
      file.path(NETWORK_OUTPUTS, "Gene_Overlap"),
      file.path(NETWORK_OUTPUTS, "Integration_Summary")
    )
    
    for (dir in dirs_to_create) {
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
      }
    }
    
    # Load comprehensive results
    results_file <- file.path(NEUROINFLAMM_DIR, "comprehensive_neuroinflammatory_analysis_enhanced.rds")
    
    if (!file.exists(results_file)) {
      stop("Comprehensive results file not found. Please run 02_Neuroinflammatory_Analysis.R first.")
    }
    
    comprehensive_results <- readRDS(results_file)
    
    # Validate data structure
    if (is.null(comprehensive_results$enrichment_results)) {
      stop("No enrichment results found in comprehensive analysis")
    }
    
    # Run network analyses
    cat("\n🔗 Creating TF correlation networks...\n")
    if (!is.null(comprehensive_results$tf_results)) {
      create_tf_correlation_networks(comprehensive_results$tf_results)
    } else {
      cat("⚠ No TF results available for correlation analysis\n")
    }
    
    cat("\n🕸️ Creating pathway interaction networks...\n")
    create_pathway_interaction_networks(comprehensive_results$enrichment_results)
    
    cat("\n⚖️ Running technology comparison analysis...\n")
    create_technology_comparison_analysis(comprehensive_results$enrichment_results, 
                                        comprehensive_results$expression_data)
    
    cat("\n🧬 Analyzing cross-dataset gene overlap...\n")
    dataset_genes <- analyze_cross_dataset_gene_overlap(comprehensive_results$enrichment_results)
    
    cat("\n📊 Creating integration summary...\n")
    summary_df <- create_network_integration_summary(comprehensive_results$enrichment_results, 
                                                   dataset_genes)
    
    # List created files
    cat("\n📁 Generated Files:\n")
    all_files <- list.files(NETWORK_OUTPUTS, recursive = TRUE, pattern = "\\.(png|csv)$")
    for (file in all_files) {
      cat("  ✓", file.path(NETWORK_OUTPUTS, file), "\n")
    }
    
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("✓ Network integration analysis complete!\n")
    cat("🔗 TF correlation networks created\n")
    cat("🕸️ Pathway interaction networks generated\n")
    cat("⚖️ Technology comparison completed\n")
    cat("📊 Integration summary available\n")
    cat("📁 All outputs saved to:", NETWORK_OUTPUTS, "\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    
    return(list(
      summary = summary_df,
      gene_overlap = dataset_genes,
      output_dir = NETWORK_OUTPUTS
    ))
    
  }, error = function(e) {
    cat("\n💥 ERROR in network integration analysis:", e$message, "\n")
    cat("🔧 Please check that 02_Neuroinflammatory_Analysis.R has been run successfully\n")
    stop(e)
  })
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if (!exists("SOURCED")) {
  network_results <- run_network_integration_analysis()
}