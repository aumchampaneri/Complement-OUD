# Load required libraries
library(STRINGdb)
library(ggplot2)
library(VennDiagram)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(biomaRt)
library(gridExtra)

# Set timeout to 300 seconds
options(timeout = 300)

# Define input and output folder paths
input_folder <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Folder - GSE225158/GSE225158/DESeq2 outputs'
output_folder <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/Super Folder - GSE225158/GSE225158/STRING outputs'

# Create output directory if it doesn't exist
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Error handling wrapper
run_safely <- function(expr, default = NULL, msg = "Error in operation") {
  tryCatch(
    expr,
    error = function(e) {
      message(msg, ": ", e$message)
      default
    }
  )
}

# Step 1: Load DESeq2 output CSVs
deg1 <- run_safely(read.csv(file.path(input_folder, "deseq2_results_M_OUD_vs_M_None.csv")),
                  msg = "Error loading condition 1 file")
deg2 <- run_safely(read.csv(file.path(input_folder, "deseq2_results_F_OUD_vs_F_None.csv")),
                  msg = "Error loading condition 2 file")

# Step 2: Filter for significant genes
sig1 <- run_safely(subset(deg1, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1),
                  msg = "Error filtering condition 1 data")
sig2 <- run_safely(subset(deg2, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1),
                  msg = "Error filtering condition 2 data")

# Extract gene names
genes1 <- data.frame(gene = sig1$gene)
genes2 <- data.frame(gene = sig2$gene)

# Step 3: Initialize STRINGdb
string_db <- run_safely(
  STRINGdb$new(version = "11.5", species = 9606, score_threshold = 700, input_directory = ""),
  msg = "Error initializing STRING database"
)

# Step 4: Map gene names to STRING IDs
mapped1 <- run_safely(string_db$map(genes1, "gene", removeUnmappedRows = FALSE),
                     msg = "Error mapping condition 1 genes")
mapped2 <- run_safely(string_db$map(genes2, "gene", removeUnmappedRows = FALSE),
                     msg = "Error mapping condition 2 genes")

# Save unmapped genes for debugging
unmapped_genes1 <- mapped1[is.na(mapped1$STRING_id), ]
unmapped_genes2 <- mapped2[is.na(mapped2$STRING_id), ]
run_safely(write.csv(unmapped_genes1, file.path(output_folder, "unmapped_genes_condition1.csv"), row.names = FALSE))
run_safely(write.csv(unmapped_genes2, file.path(output_folder, "unmapped_genes_condition2.csv"), row.names = FALSE))

# Remove unmapped rows
mapped1 <- mapped1[!is.na(mapped1$STRING_id), ]
mapped2 <- mapped2[!is.na(mapped2$STRING_id), ]

# Save mapped genes for debugging
run_safely(write.csv(mapped1, file.path(output_folder, "mapped_genes_condition1.csv"), row.names = FALSE))
run_safely(write.csv(mapped2, file.path(output_folder, "mapped_genes_condition2.csv"), row.names = FALSE))

# Step 5: Get enrichment results
enrichment1 <- run_safely(string_db$get_enrichment(mapped1$STRING_id), msg = "Error getting enrichment for condition 1")
enrichment2 <- run_safely(string_db$get_enrichment(mapped2$STRING_id), msg = "Error getting enrichment for condition 2")

# Save enrichment results
run_safely(write.csv(enrichment1, file.path(output_folder, "enrichment1.csv"), row.names = FALSE))
run_safely(write.csv(enrichment2, file.path(output_folder, "enrichment2.csv"), row.names = FALSE))

# Step 6: Enrichment Comparison (Shared/Unique Pathways)
e1 <- run_safely(enrichment1[, c("term", "description", "fdr", "category")],
                msg = "Error extracting columns from enrichment1")
e2 <- run_safely(enrichment2[, c("term", "description", "fdr", "category")],
                msg = "Error extracting columns from enrichment2")

if (!is.null(e1) && !is.null(e2)) {
  e1$source <- "Condition_1"
  e2$source <- "Condition_2"

  merged <- run_safely(merge(e1, e2, by = "term", all = TRUE, suffixes = c("_1", "_2")),
                      msg = "Error merging enrichment data")

  if (!is.null(merged)) {
    shared <- merged[!is.na(merged$fdr_1) & !is.na(merged$fdr_2), ]
    unique_1 <- merged[!is.na(merged$fdr_1) & is.na(merged$fdr_2), ]
    unique_2 <- merged[is.na(merged$fdr_1) & !is.na(merged$fdr_2), ]

    run_safely(write.csv(shared, file.path(output_folder, "shared_enriched_pathways.csv"), row.names = FALSE))
    run_safely(write.csv(unique_1, file.path(output_folder, "unique_condition1_pathways.csv"), row.names = FALSE))
    run_safely(write.csv(unique_2, file.path(output_folder, "unique_condition2_pathways.csv"), row.names = FALSE))

    # Step 7: Venn Diagram of Enriched Terms
    run_safely({
      venn.plot <- draw.pairwise.venn(
        area1 = length(e1$term),
        area2 = length(e2$term),
        cross.area = length(intersect(e1$term, e2$term)),
        category = c("Condition 1", "Condition 2"),
        fill = c("skyblue", "pink1"),
        lty = "blank",
        cex = 2,
        cat.cex = 2,
        cat.pos = c(-20, 20),
        cat.dist = 0.05
      )
      pdf(file.path(output_folder, "venn_diagram.pdf"))
      grid.newpage()
      grid.draw(venn.plot)
      dev.off()
    }, msg = "Error creating Venn diagram")

    # Step 8: Interactive Plotly Visualization of Shared Pathways
    run_safely({
      if (nrow(shared) > 0) {
        shared$logFDR_1 <- -log10(shared$fdr_1)
        shared$logFDR_2 <- -log10(shared$fdr_2)

        p <- plot_ly(
          data = shared,
          x = ~logFDR_1,
          y = ~logFDR_2,
          type = 'scatter',
          mode = 'markers',
          text = ~paste(description_1, "<br>Category:", category_1,
                        "<br>FDR Condition 1:", signif(fdr_1, 3),
                        "<br>FDR Condition 2:", signif(fdr_2, 3)),
          hoverinfo = 'text',
          marker = list(color = 'royalblue', size = 10)
        ) %>%
          layout(
            title = "Shared Enriched Pathways (Interactive)",
            xaxis = list(title = "-log10 FDR (Condition 1)"),
            yaxis = list(title = "-log10 FDR (Condition 2)")
          )

        saveWidget(p, file.path(output_folder, "shared_pathways_plot.html"))
      }
    }, msg = "Error creating interactive plot")
  }
}

# Step 9: STRING Network Clustering and Individual Cluster Processing
clusters1 <- run_safely(string_db$get_clusters(mapped1$STRING_id), msg = "Error getting clusters for condition 1")
clusters2 <- run_safely(string_db$get_clusters(mapped2$STRING_id), msg = "Error getting clusters for condition 2")

run_safely(write.csv(clusters1, file.path(output_folder, "clusters_condition1.csv"), row.names = FALSE))
run_safely(write.csv(clusters2, file.path(output_folder, "clusters_condition2.csv"), row.names = FALSE))

# Individual processing - Condition 1
if (!is.null(clusters1) && !is.null(mapped1)) {
  run_safely({
    mapped1_clusters <- merge(mapped1, clusters1, by.x = "STRING_id", by.y = "stringId")

    # Create summary file for condition 1
    summary_file <- file.path(output_folder, "cluster_summary_condition1.txt")
    cat("Clusters in Condition 1:\n", file = summary_file)

    for (c in unique(mapped1_clusters$cluster)) {
      cluster_genes <- mapped1_clusters[mapped1_clusters$cluster == c, ]

      # Add to summary file
      cat(paste0("Cluster ", c, ": ", nrow(cluster_genes), " genes\n"),
          file = summary_file, append = TRUE)
      cat(paste0(cluster_genes$preferred_name, collapse = ", "),
          file = summary_file, append = TRUE)
      cat("\n\n", file = summary_file, append = TRUE)

      # Only process clusters with at least 3 genes
      if (nrow(cluster_genes) >= 3) {
        # Get enrichment for this specific cluster
        result <- string_db$get_enrichment(cluster_genes$STRING_id)

        # Only save if we got results
        if (is.data.frame(result) && nrow(result) > 0) {
          result$cluster <- c
          write.csv(result,
                    file.path(output_folder, paste0("cluster_", c, "_enrichment_cond1.csv")),
                    row.names = FALSE)
        }
      }
    }
  }, msg = "Error processing clusters for condition 1")
}

# Individual processing - Condition 2
if (!is.null(clusters2) && !is.null(mapped2)) {
  run_safely({
    mapped2_clusters <- merge(mapped2, clusters2, by.x = "STRING_id", by.y = "stringId")

    # Create summary file for condition 2
    summary_file <- file.path(output_folder, "cluster_summary_condition2.txt")
    cat("Clusters in Condition 2:\n", file = summary_file)

    for (c in unique(mapped2_clusters$cluster)) {
      cluster_genes <- mapped2_clusters[mapped2_clusters$cluster == c, ]

      # Add to summary file
      cat(paste0("Cluster ", c, ": ", nrow(cluster_genes), " genes\n"),
          file = summary_file, append = TRUE)
      cat(paste0(cluster_genes$preferred_name, collapse = ", "),
          file = summary_file, append = TRUE)
      cat("\n\n", file = summary_file, append = TRUE)

      # Only process clusters with at least 3 genes
      if (nrow(cluster_genes) >= 3) {
        # Get enrichment for this specific cluster
        result <- string_db$get_enrichment(cluster_genes$STRING_id)

        # Only save if we got results
        if (is.data.frame(result) && nrow(result) > 0) {
          result$cluster <- c
          write.csv(result,
                    file.path(output_folder, paste0("cluster_", c, "_enrichment_cond2.csv")),
                    row.names = FALSE)
        }
      }
    }
  }, msg = "Error processing clusters for condition 2")
}

cat("Analysis completed successfully.\n")