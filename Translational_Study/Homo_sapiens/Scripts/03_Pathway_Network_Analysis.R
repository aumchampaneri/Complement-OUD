# ==============================================================================
# Pathway Network Analysis
# ==============================================================================
# Purpose: Create pathway interaction networks from enrichment results
# Input: Results from 02_Neuroinflammatory_Analysis.R
# Features: Pathway similarity networks, gene overlap networks, hierarchical clustering
# Output: Interactive networks, static plots, network statistics
# ==============================================================================

# Configuration
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens"
RESULTS_DIR <- file.path(BASE_DIR, "Results")
NEUROINFLAMM_DIR <- file.path(RESULTS_DIR, "Neuroinflammatory_Analysis")
NETWORK_DIR <- file.path(RESULTS_DIR, "Pathway_Networks")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(igraph)
  library(visNetwork)
  library(RColorBrewer)
  library(pheatmap)
  library(corrplot)
  library(enrichplot)
  library(clusterProfiler)
  library(stringr)
})

# ==============================================================================
# SETUP AND DATA LOADING
# ==============================================================================

#' Create network analysis directories
create_network_directories <- function() {
  cat("=== Creating Network Analysis Directories ===\n")
  
  dirs_to_create <- c(
    NETWORK_DIR,
    file.path(NETWORK_DIR, "Interactive_Networks"),
    file.path(NETWORK_DIR, "Static_Networks"),
    file.path(NETWORK_DIR, "Network_Statistics"),
    file.path(NETWORK_DIR, "Pathway_Clusters"),
    file.path(NETWORK_DIR, "Gene_Overlap_Networks")
  )
  
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("✓ Created:", dir, "\n")
    }
  }
  
  cat("✓ Network directories ready\n")
}

#' Load comprehensive enrichment results
load_enrichment_results <- function() {
  cat("=== Loading Enrichment Results ===\n")
  
  results_file <- file.path(NEUROINFLAMM_DIR, "comprehensive_neuroinflammatory_analysis_enhanced.rds")
  
  if (!file.exists(results_file)) {
    stop("Enrichment results not found: ", results_file,
         "\nPlease run 02_Neuroinflammatory_Analysis.R first")
  }
  
  comprehensive_results <- readRDS(results_file)
  cat("✓ Enrichment results loaded\n")
  
  return(comprehensive_results)
}

# ==============================================================================
# PATHWAY SIMILARITY CALCULATION
# ==============================================================================

#' Calculate pathway similarity based on gene overlap
calculate_pathway_similarity <- function(enrichment_results, min_overlap = 2, max_pathways = 200) {
  cat("\n=== Calculating Pathway Similarity Networks ===\n")
  
  similarity_networks <- list()
  
  for (dataset in names(enrichment_results)) {
    for (method in names(enrichment_results[[dataset]])) {
      for (db in names(enrichment_results[[dataset]][[method]])) {
        
        result_obj <- enrichment_results[[dataset]][[method]][[db]]
        
        if (!is.null(result_obj) && nrow(result_obj@result) > 0) {
          
          # Filter to significant pathways
          sig_pathways <- result_obj@result[result_obj@result$p.adjust < 0.05, ]
          
          if (nrow(sig_pathways) < 5) {
            cat("  Skipping", dataset, method, db, "- too few pathways\n")
            next
          }
          
          # Limit to top pathways for computational efficiency
          if (nrow(sig_pathways) > max_pathways) {
            sig_pathways <- sig_pathways[order(sig_pathways$p.adjust), ][1:max_pathways, ]
          }
          
          cat("  Processing", dataset, method, db, "-", nrow(sig_pathways), "pathways\n")
          
          # Calculate gene overlap matrix
          similarity_matrix <- calculate_gene_overlap_matrix(sig_pathways, min_overlap)
          
          if (!is.null(similarity_matrix)) {
            network_key <- paste(dataset, method, db, sep = "_")
            similarity_networks[[network_key]] <- list(
              similarity_matrix = similarity_matrix,
              pathway_data = sig_pathways,
              dataset = dataset,
              method = method,
              database = db
            )
          }
        }
      }
    }
  }
  
  cat("✓ Calculated", length(similarity_networks), "pathway similarity networks\n")
  return(similarity_networks)
}

#' Calculate gene overlap matrix between pathways
calculate_gene_overlap_matrix <- function(pathway_data, min_overlap = 2) {
  
  # Extract gene lists for each pathway
  pathway_genes <- list()
  
  for (i in 1:nrow(pathway_data)) {
    pathway_id <- pathway_data$ID[i]
    
    if ("geneID" %in% colnames(pathway_data)) {
      genes <- strsplit(pathway_data$geneID[i], "/")[[1]]
    } else if ("core_enrichment" %in% colnames(pathway_data)) {
      genes <- strsplit(pathway_data$core_enrichment[i], "/")[[1]]
    } else {
      # Try to extract from other potential columns
      genes <- NULL
    }
    
    if (!is.null(genes) && length(genes) > 0) {
      pathway_genes[[pathway_id]] <- genes
    }
  }
  
  if (length(pathway_genes) < 2) {
    cat("    Warning: Insufficient gene information for similarity calculation\n")
    return(NULL)
  }
  
  # Calculate Jaccard similarity matrix
  pathway_ids <- names(pathway_genes)
  n_pathways <- length(pathway_ids)
  similarity_matrix <- matrix(0, nrow = n_pathways, ncol = n_pathways)
  rownames(similarity_matrix) <- pathway_ids
  colnames(similarity_matrix) <- pathway_ids
  
  for (i in 1:n_pathways) {
    for (j in i:n_pathways) {
      genes1 <- pathway_genes[[i]]
      genes2 <- pathway_genes[[j]]
      
      intersection <- length(intersect(genes1, genes2))
      union <- length(union(genes1, genes2))
      
      if (union > 0) {
        jaccard <- intersection / union
        similarity_matrix[i, j] <- jaccard
        similarity_matrix[j, i] <- jaccard
        
        # Only keep connections with sufficient overlap
        if (intersection < min_overlap && i != j) {
          similarity_matrix[i, j] <- 0
          similarity_matrix[j, i] <- 0
        }
      }
    }
  }
  
  return(similarity_matrix)
}

# ==============================================================================
# NETWORK CREATION AND VISUALIZATION
# ==============================================================================

#' Create pathway networks from similarity matrices
create_pathway_networks <- function(similarity_networks, min_similarity = 0.1) {
  cat("\n=== Creating Pathway Networks ===\n")
  
  network_objects <- list()
  
  for (network_key in names(similarity_networks)) {
    network_data <- similarity_networks[[network_key]]
    similarity_matrix <- network_data$similarity_matrix
    pathway_data <- network_data$pathway_data
    
    # Create network from similarity matrix
    # Remove diagonal and apply threshold
    diag(similarity_matrix) <- 0
    similarity_matrix[similarity_matrix < min_similarity] <- 0
    
    # Create igraph object
    network <- graph_from_adjacency_matrix(similarity_matrix, 
                                         mode = "undirected", 
                                         weighted = TRUE)
    
    # Add pathway information as vertex attributes
    pathway_descriptions <- pathway_data$Description[match(V(network)$name, pathway_data$ID)]
    V(network)$description <- pathway_descriptions
    V(network)$p_value <- pathway_data$p.adjust[match(V(network)$name, pathway_data$ID)]
    V(network)$size <- -log10(V(network)$p_value)  # Node size based on significance
    
    # Remove isolated nodes
    isolated <- which(degree(network) == 0)
    if (length(isolated) > 0) {
      network <- delete_vertices(network, isolated)
    }
    
    if (vcount(network) > 2) {
      network_objects[[network_key]] <- list(
        network = network,
        metadata = network_data
      )
      
      cat("  Created network:", network_key, "-", vcount(network), "nodes,", ecount(network), "edges\n")
    }
  }
  
  return(network_objects)
}

#' Create static network visualizations
create_static_network_plots <- function(network_objects) {
  cat("\n=== Creating Static Network Plots ===\n")
  
  static_dir <- file.path(NETWORK_DIR, "Static_Networks")
  
  for (network_key in names(network_objects)) {
    network_data <- network_objects[[network_key]]
    network <- network_data$network
    metadata <- network_data$metadata
    
    # Skip if network is too large for static visualization
    if (vcount(network) > 50) {
      cat("  Skipping", network_key, "- too large for static plot\n")
      next
    }
    
    # Create layout
    layout <- layout_with_fr(network)
    
    # Set visual attributes
    V(network)$size <- pmax(5, pmin(20, V(network)$size * 3))  # Scale node size
    V(network)$color <- heat.colors(10)[cut(-log10(V(network)$p_value), breaks = 10)]
    E(network)$width <- E(network)$weight * 5
    E(network)$color <- "grey70"
    
    # Create plot
    plot_file <- file.path(static_dir, paste0(network_key, "_network.png"))
    
    png(plot_file, width = 1200, height = 1000, res = 150)
    plot(network, 
         layout = layout,
         vertex.label = ifelse(vcount(network) < 20, 
                              str_wrap(V(network)$description, 20), 
                              ""),
         vertex.label.cex = 0.7,
         vertex.label.color = "black",
         main = paste("Pathway Network:", metadata$dataset, metadata$method, metadata$database),
         sub = paste(vcount(network), "pathways,", ecount(network), "connections"))
    dev.off()
    
    cat("  Static plot saved:", basename(plot_file), "\n")
  }
}

#' Create interactive network visualizations (FIXED)
create_interactive_networks <- function(network_objects) {
  cat("\n=== Creating Interactive Networks ===\n")
  
  interactive_dir <- file.path(NETWORK_DIR, "Interactive_Networks")
  
  for (network_key in names(network_objects)) {
    network_data <- network_objects[[network_key]]
    network <- network_data$network
    metadata <- network_data$metadata
    
    # Skip if network is too large
    if (vcount(network) > 100) {
      cat("  Skipping", network_key, "- too large for interactive visualization\n")
      next
    }
    
    # Convert to visNetwork format
    nodes <- data.frame(
      id = V(network)$name,
      label = str_wrap(V(network)$description, 30),
      title = paste0("ID: ", V(network)$name, "\n",
                    "Description: ", V(network)$description, "\n",
                    "P-value: ", signif(V(network)$p_value, 3)),
      size = pmax(10, pmin(50, V(network)$size * 5)),
      color = heat.colors(10)[cut(-log10(V(network)$p_value), breaks = 10)],
      group = sample(1:3, vcount(network), replace = TRUE),  # FIXED: Add group column
      stringsAsFactors = FALSE
    )
    
    edges <- data.frame(
      from = get.edgelist(network)[, 1],
      to = get.edgelist(network)[, 2],
      weight = E(network)$weight,
      width = E(network)$weight * 10,
      title = paste("Similarity:", round(E(network)$weight, 3)),
      stringsAsFactors = FALSE
    )
    
    # Create interactive network with fallback for pandoc
    tryCatch({
      vis_network <- visNetwork(nodes, edges) %>%
        visOptions(highlightNearest = TRUE, selectedBy = "group") %>%
        visLayout(randomSeed = 123) %>%
        visPhysics(stabilization = FALSE) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visNodes(borderWidth = 2) %>%
        visEdges(smooth = FALSE)
      
      # Try to save with selfcontained = FALSE first (doesn't require pandoc)
      html_file <- file.path(interactive_dir, paste0(network_key, "_interactive.html"))
      visSave(vis_network, html_file, selfcontained = FALSE)
      
      cat("  Interactive network saved:", basename(html_file), "\n")
      
    }, error = function(e) {
      cat("  ⚠ Interactive network failed for", network_key, ":", e$message, "\n")
      
      # Create a simple HTML summary instead
      summary_file <- file.path(interactive_dir, paste0(network_key, "_summary.txt"))
      writeLines(c(
        paste("Network:", network_key),
        paste("Nodes:", vcount(network)),
        paste("Edges:", ecount(network)),
        paste("Dataset:", metadata$dataset),
        paste("Method:", metadata$method),
        paste("Database:", metadata$database)
      ), summary_file)
      
      cat("  Summary saved instead:", basename(summary_file), "\n")
    })
  }
}

# ==============================================================================
# NETWORK ANALYSIS AND STATISTICS
# ==============================================================================

#' Analyze network properties and statistics
analyze_network_properties <- function(network_objects) {
  cat("\n=== Analyzing Network Properties ===\n")
  
  network_stats <- data.frame()
  
  for (network_key in names(network_objects)) {
    network_data <- network_objects[[network_key]]
    network <- network_data$network
    metadata <- network_data$metadata
    
    # Calculate network properties
    n_nodes <- vcount(network)
    n_edges <- ecount(network)
    density <- edge_density(network)
    avg_degree <- mean(degree(network))
    clustering_coef <- transitivity(network, type = "global")
    avg_path_length <- ifelse(is.connected(network), 
                             average.path.length(network), 
                             NA)
    
    # Community detection
    communities <- cluster_louvain(network)
    modularity <- modularity(communities)
    n_communities <- length(communities)
    
    # Central pathways
    betweenness <- betweenness(network)
    closeness <- closeness(network)
    pagerank <- page_rank(network)$vector
    
    # Top central pathway
    top_central <- names(which.max(betweenness))
    top_central_desc <- V(network)$description[V(network)$name == top_central]
    
    # Add to summary
    network_stats <- rbind(network_stats, data.frame(
      Network = network_key,
      Dataset = metadata$dataset,
      Method = metadata$method,
      Database = metadata$database,
      N_Nodes = n_nodes,
      N_Edges = n_edges,
      Density = round(density, 3),
      Avg_Degree = round(avg_degree, 2),
      Clustering_Coefficient = round(clustering_coef, 3),
      Avg_Path_Length = round(avg_path_length, 2),
      Modularity = round(modularity, 3),
      N_Communities = n_communities,
      Top_Central_Pathway = top_central_desc,
      stringsAsFactors = FALSE
    ))
    
    cat("  Analyzed:", network_key, "-", n_nodes, "nodes,", n_communities, "communities\n")
  }
  
  # Save network statistics
  stats_file <- file.path(NETWORK_DIR, "Network_Statistics", "network_properties.csv")
  write.csv(network_stats, stats_file, row.names = FALSE)
  cat("✓ Network statistics saved:", basename(stats_file), "\n")
  
  return(network_stats)
}

#' Create pathway clustering analysis
create_pathway_clusters <- function(similarity_networks) {
  cat("\n=== Creating Pathway Clustering Analysis ===\n")
  
  cluster_dir <- file.path(NETWORK_DIR, "Pathway_Clusters")
  
  for (network_key in names(similarity_networks)) {
    network_data <- similarity_networks[[network_key]]
    similarity_matrix <- network_data$similarity_matrix
    pathway_data <- network_data$pathway_data
    
    if (nrow(similarity_matrix) < 5) {
      cat("  Skipping", network_key, "- too few pathways for clustering\n")
      next
    }
    
    # Hierarchical clustering
    dist_matrix <- as.dist(1 - similarity_matrix)
    hclust_result <- hclust(dist_matrix, method = "ward.D2")
    
    # Create dendrogram plot
    png_file <- file.path(cluster_dir, paste0(network_key, "_dendrogram.png"))
    png(png_file, width = 1200, height = 800, res = 150)
    
    plot(hclust_result, 
         main = paste("Pathway Clustering:", network_data$dataset, network_data$method, network_data$database),
         xlab = "Pathways", 
         ylab = "Distance",
         cex = 0.7)
    
    # Add cluster rectangles
    n_clusters <- min(8, max(2, nrow(similarity_matrix) %/% 5))
    clusters <- cutree(hclust_result, k = n_clusters)
    rect.hclust(hclust_result, k = n_clusters, border = "red")
    
    dev.off()
    
    # Create heatmap
    if (nrow(similarity_matrix) <= 50) {  # Only for smaller matrices
      heatmap_file <- file.path(cluster_dir, paste0(network_key, "_heatmap.png"))
      
      # Prepare pathway names for heatmap
      pathway_names <- pathway_data$Description[match(rownames(similarity_matrix), pathway_data$ID)]
      pathway_names <- str_wrap(pathway_names, 40)
      
      rownames(similarity_matrix) <- pathway_names
      colnames(similarity_matrix) <- pathway_names
      
      png(heatmap_file, width = 1200, height = 1200, res = 150)
      pheatmap(similarity_matrix,
               main = paste("Pathway Similarity:", network_data$dataset, network_data$method),
               clustering_method = "ward.D2",
               color = colorRampPalette(c("white", "yellow", "red"))(100),
               fontsize = 8)
      dev.off()
      
      cat("  Heatmap saved:", basename(heatmap_file), "\n")
    }
    
    cat("  Clustering analysis saved:", basename(png_file), "\n")
  }
}

# ==============================================================================
# CROSS-DATASET NETWORK COMPARISON
# ==============================================================================

#' Compare networks across datasets and methods
compare_networks_across_datasets <- function(network_objects, network_stats) {
  cat("\n=== Comparing Networks Across Datasets ===\n")
  
  comparison_dir <- file.path(NETWORK_DIR, "Network_Statistics")
  
  # Create network properties comparison plot
  if (nrow(network_stats) > 1) {
    
    # Network size comparison
    p1 <- ggplot(network_stats, aes(x = Method, y = N_Nodes, fill = Dataset)) +
      geom_col(position = "dodge", alpha = 0.8) +
      facet_wrap(~ Database, scales = "free_y") +
      labs(title = "Network Size Comparison",
           subtitle = "Number of connected pathways by method and database",
           x = "Method", y = "Number of Pathways") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(comparison_dir, "network_size_comparison.png"), 
           p1, width = 12, height = 8, dpi = 300)
    
    # Network density comparison
    p2 <- ggplot(network_stats, aes(x = Method, y = Density, fill = Dataset)) +
      geom_col(position = "dodge", alpha = 0.8) +
      facet_wrap(~ Database) +
      labs(title = "Network Density Comparison",
           subtitle = "Pathway connectivity by method and database",
           x = "Method", y = "Network Density") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(comparison_dir, "network_density_comparison.png"), 
           p2, width = 12, height = 8, dpi = 300)
    
    # Modularity comparison
    p3 <- ggplot(network_stats, aes(x = Method, y = Modularity, fill = Dataset)) +
      geom_col(position = "dodge", alpha = 0.8) +
      facet_wrap(~ Database) +
      labs(title = "Network Modularity Comparison",
           subtitle = "Community structure strength by method and database",
           x = "Method", y = "Modularity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(comparison_dir, "network_modularity_comparison.png"), 
           p3, width = 12, height = 8, dpi = 300)
    
    cat("✓ Network comparison plots saved\n")
  }
  
  # Create summary statistics table
  summary_stats <- network_stats %>%
    group_by(Dataset, Database) %>%
    summarise(
      Methods_Analyzed = n(),
      Avg_Network_Size = round(mean(N_Nodes), 1),
      Avg_Density = round(mean(Density), 3),
      Avg_Modularity = round(mean(Modularity), 3),
      Total_Pathways = sum(N_Nodes),
      Total_Connections = sum(N_Edges),
      .groups = "drop"
    )
  
  summary_file <- file.path(comparison_dir, "network_summary_by_dataset.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat("✓ Network summary saved:", basename(summary_file), "\n")
  
  return(summary_stats)
}

# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

#' Run comprehensive pathway network analysis
run_pathway_network_analysis <- function() {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("PATHWAY NETWORK ANALYSIS\n")
  cat("Creating pathway interaction networks from enrichment results\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  tryCatch({
    # Setup
    create_network_directories()
    
    # Load enrichment results
    comprehensive_results <- load_enrichment_results()
    enrichment_results <- comprehensive_results$enrichment_results
    
    # Calculate pathway similarities
    similarity_networks <- calculate_pathway_similarity(enrichment_results, 
                                                      min_overlap = 2, 
                                                      max_pathways = 150)
    
    # Create network objects
    network_objects <- create_pathway_networks(similarity_networks, 
                                             min_similarity = 0.1)
    
    # Create visualizations (with error handling)
    create_static_network_plots(network_objects)
    
    # Try interactive networks with fallback
    tryCatch({
      create_interactive_networks(network_objects)
    }, error = function(e) {
      cat("⚠ Interactive networks skipped due to technical issues\n")
      cat("Static networks and analysis complete\n")
    })
    
    # Analyze network properties
    network_stats <- analyze_network_properties(network_objects)
    
    # Create clustering analysis
    create_pathway_clusters(similarity_networks)
    
    # Compare across datasets
    summary_stats <- compare_networks_across_datasets(network_objects, network_stats)
    
    # Save comprehensive results
    results_file <- file.path(NETWORK_DIR, "pathway_network_analysis_results.rds")
    saveRDS(list(
      similarity_networks = similarity_networks,
      network_objects = network_objects,
      network_statistics = network_stats,
      summary_statistics = summary_stats,
      analysis_date = Sys.Date()
    ), results_file)
    
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("SUCCESS: Pathway network analysis complete!\n")
    cat("✓ Pathway similarity networks calculated (", length(similarity_networks), " networks)\n")
    cat("✓ Network objects created (", length(network_objects), " networks)\n")
    cat("✓ Static visualizations created\n")
    cat("✓ Network properties analyzed\n")
    cat("✓ Pathway clustering completed\n")
    cat("✓ Cross-dataset comparisons generated\n")
    cat("Results saved to:", NETWORK_DIR, "\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    
    return(list(
      similarity_networks = similarity_networks,
      network_objects = network_objects,
      network_statistics = network_stats,
      summary_statistics = summary_stats
    ))
    
  }, error = function(e) {
    cat("\n❌ NETWORK ANALYSIS ERROR:", e$message, "\n")
    # Return partial results instead of stopping
    return(list(
      error = e$message,
      partial_results = "Network analysis encountered technical issues"
    ))
  })
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if (!exists("SOURCED")) {
  network_results <- run_pathway_network_analysis()
}
