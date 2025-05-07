# Load required libraries
library(decoupleR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
library(msigdbr)
library(pheatmap)
library(limma)
library(tibble)

# Check if required packages are installed and install if needed
required_packages <- c("decoupleR", "dplyr", "tidyr", "ggplot2", "ComplexHeatmap",
                      "circlize", "org.Hs.eg.db", "msigdbr", "pheatmap", "limma", "tibble")

missing_packages <- required_packages[!requireNamespace(required_packages, quietly = TRUE)]
if(length(missing_packages) > 0) {
  if(!require(BiocManager, quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(missing_packages)
}

# Set seed for reproducibility
set.seed(123)

# Create output directory
output_dir <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/Glial Analysis/decoupleR_results/'
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load edgeR results
norm_counts <- read.csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/Glial Analysis/edgeR Glial outputs/log_cpm_counts.csv', row.names=1)
meta <- read.csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/meta.csv', row.names=1)
meta_glial <- meta[colnames(norm_counts),]

# Create sample groups for visualization
meta_glial$Sex_Dx_OUD <- factor(paste(meta_glial$Sex, meta_glial$Dx_OUD, sep="_"))

# Load edgeR contrast results
contrast_files <- list.files('/Users/aumchampaneri/PycharmProjects/Complement-OUD/Glial Analysis/edgeR Glial outputs/',
                            pattern = "edgeR_results_", full.names = TRUE)
contrast_results <- list()

for(file in contrast_files) {
  contrast_name <- gsub(".*edgeR_results_|\\.csv", "", file)
  contrast_results[[contrast_name]] <- read.csv(file)
}

# 1. Get prior knowledge resources

# a. PROGENy for pathway activities
progeny_human <- get_progeny(organism = "human", top = 500)

# Check PROGENy column structure and adapt as needed
cat("PROGENy columns:", paste(colnames(progeny_human), collapse=", "), "\n")

# Fix PROGENy column names
if ("weight" %in% colnames(progeny_human)) {
  # Use explicit column name reference instead of bare name
  progeny_human <- progeny_human %>%
    dplyr::rename("mor" = "weight")
}

# Add likelihood column if missing
if (!"likelihood" %in% colnames(progeny_human)) {
  progeny_human <- progeny_human %>%
    mutate(likelihood = 1)
}

# b. DoRothEA for transcription factor activities
dorothea_human <- get_dorothea(organism = "human", levels = c("A", "B", "C"))

# Check DoRothEA column structure and adapt as needed
cat("DoRothEA columns:", paste(colnames(dorothea_human), collapse=", "), "\n")

# Same fix for DoRothEA
if ("weight" %in% colnames(dorothea_human)) {
  dorothea_human <- dorothea_human %>%
    dplyr::rename("mor" = "weight")
}

# Add likelihood column if missing for DoRothEA
if (!"likelihood" %in% colnames(dorothea_human)) {
  dorothea_human <- dorothea_human %>%
    mutate(likelihood = 1)
}

# c. Custom MSigDB gene sets focused on relevant pathways
msigdb_collections <- list(
  hallmark = msigdbr(species = "Homo sapiens", collection = "H")
)

# Get all C2 collection and filter by name pattern instead of using subcollections
c2_all <- msigdbr(species = "Homo sapiens", collection = "C2")
msigdb_collections$kegg <- c2_all %>% filter(grepl("^KEGG", gs_name))
msigdb_collections$reactome <- c2_all %>% filter(grepl("^REACTOME", gs_name))

# Get GO:BP collection
msigdb_collections$go_bp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")

# Filter for pathways of interest based on your previous ORA script
pathways_of_interest <- c(
  # Neuroinflammation & Glial Activation
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNF_ALPHA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_TNF_SIGNALING_PATHWAY",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "REACTOME_INNATE_IMMUNE_SYSTEM",
  "REACTOME_INTERFERON_SIGNALING",

  # Intracellular Signaling
  "HALLMARK_COMPLEMENT",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "KEGG_PI3K_AKT_SIGNALING_PATHWAY",
  "KEGG_APOPTOSIS",
  "REACTOME_APOPTOSIS",

  # Neurotransmission
  "KEGG_DOPAMINERGIC_SYNAPSE",
  "KEGG_SEROTONERGIC_SYNAPSE",
  "KEGG_GABAERGIC_SYNAPSE",
  "KEGG_GLUTAMATERGIC_SYNAPSE",
  "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
  "REACTOME_NEUROTRANSMITTER_RECEPTOR_BINDING_AND_DOWNSTREAM_TRANSMISSION_IN_THE_POSTSYNAPTIC_CELL",
  "REACTOME_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING",

  # Bonus
  "HALLMARK_HYPOXIA"
)

# Create custom gene set collection focused on these
custom_gs <- bind_rows(lapply(names(msigdb_collections), function(coll) {
  df <- msigdb_collections[[coll]]
  df %>%
    filter(grepl(paste(pathways_of_interest, collapse="|"), gs_name, ignore.case=TRUE) |
           grepl("complement|inflam|microglia|astro|glia|oligo|immune|cytokine|neuro", gs_name, ignore.case=TRUE)) %>%
    mutate(collection = coll)
}))

# Format for decoupleR
custom_net <- custom_gs %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(source = gs_name, target = gene_symbol) %>%
  mutate(mor = 1) %>%
  distinct(source, target, .keep_all = TRUE)  # Remove duplicate edges

# 2. Run decoupleR methods on each contrast
all_methods_results <- list()

for(contrast_name in names(contrast_results)) {
  # Extract stats from contrast results
  degs <- contrast_results[[contrast_name]]

  # Debug: Check column names
  cat("Columns in edgeR results:", paste(colnames(degs), collapse=", "), "\n")

  # Check if the expected columns exist
  if(!"gene" %in% colnames(degs)) {
    # Try to find gene name column
    if("genes" %in% colnames(degs)) {
      degs$gene <- degs$genes
    } else if("ID" %in% colnames(degs)) {
      degs$gene <- degs$ID
    } else {
      cat("ERROR: Cannot find gene identifier column in contrast results\n")
      next
    }
  }

  if(!"logFC" %in% colnames(degs)) {
    if("log2FoldChange" %in% colnames(degs)) {
      degs$logFC <- degs$log2FoldChange
    } else {
      cat("ERROR: Cannot find logFC column in contrast results\n")
      next
    }
  }

  stat_table <- degs %>%
    dplyr::select(gene, logFC) %>%
    column_to_rownames("gene")

  # Prepare inputs
  mat <- as.matrix(stat_table)

  # Run multiple methods on pathway gene sets
  pathway_activities <- run_wmean(mat, network = custom_net, .source = "source",
                                 .target = "target", .mor = "mor", times = 1000)

  # Run multiple methods on PROGENy pathways
  progeny_activities <- run_mlm(mat, network = progeny_human, .source = "source",
                                .target = "target", .mor = "mor", .likelihood = "likelihood")

  # Run multiple methods on DoRothEA transcription factors
  tf_activities <- run_ulm(mat, network = dorothea_human, .source = "source",
                          .target = "target", .mor = "mor", .likelihood = "likelihood")

  # Store results
  all_methods_results[[contrast_name]] <- list(
    pathways = pathway_activities,
    progeny = progeny_activities,
    tfs = tf_activities
  )

  # Save individual results
  write.csv(pathway_activities, file = paste0(output_dir, "pathways_", contrast_name, ".csv"), row.names = FALSE)
  write.csv(progeny_activities, file = paste0(output_dir, "progeny_", contrast_name, ".csv"), row.names = FALSE)
  write.csv(tf_activities, file = paste0(output_dir, "tfs_", contrast_name, ".csv"), row.names = FALSE)
}

# 3. Generate visualizations for each contrast
for(contrast_name in names(all_methods_results)) {
  results <- all_methods_results[[contrast_name]]

  # Create contrast-specific visualization directory
  contrast_viz_dir <- paste0(output_dir, contrast_name, "/")
  dir.create(contrast_viz_dir, showWarnings = FALSE)

  # Top pathways heatmap - with error checking
  top_pathways <- results$pathways %>%
    filter(statistic == "norm_wmean") %>%
    top_n(20, abs(score)) %>%
    arrange(desc(abs(score)))

  if(nrow(top_pathways) > 0) {
    p1 <- ggplot(top_pathways, aes(x = "Activity", y = reorder(source, score), fill = score)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      theme_minimal() +
      labs(title = paste0("Top Pathway Activities - ", contrast_name),
           x = "", y = "",
           subtitle = "Normalized weighted means method") +
      theme(axis.text.y = element_text(size = 10))

    ggsave(paste0(contrast_viz_dir, "top_pathways.pdf"), p1, width = 10, height = 8)
  } else {
    cat("No pathway data for contrast", contrast_name, "\n")
  }

  # Apply similar checks to other visualizations
  if(nrow(results$progeny) > 0) {
    p2 <- ggplot(results$progeny, aes(x = statistic, y = source, fill = score)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      theme_minimal() +
      labs(title = paste0("PROGENy Pathway Activities - ", contrast_name),
           x = "Method", y = "Pathway") +
      theme(axis.text.y = element_text(size = 10))

    ggsave(paste0(contrast_viz_dir, "progeny_pathways.pdf"), p2, width = 8, height = 6)
  }

  # Top TFs - with error checking
  # Check TF results structure first
  cat("TF results columns for", contrast_name, ":", paste(colnames(results$tfs), collapse=", "), "\n")

  # Use the correct statistic column value
  tf_stat_value <- if("norm_ulm" %in% unique(results$tfs$statistic)) "norm_ulm" else "ulm"

  top_tfs <- results$tfs %>%
    filter(statistic == tf_stat_value) %>%
    top_n(30, abs(score)) %>%
    arrange(desc(abs(score)))

  if(nrow(top_tfs) > 0) {
    p3 <- ggplot(top_tfs, aes(x = "Activity", y = reorder(source, score), fill = score)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      theme_minimal() +
      labs(title = paste0("Transcription Factor Activities - ", contrast_name),
           x = "", y = "",
           subtitle = paste0("Using ", tf_stat_value, " statistic")) +
      theme(axis.text.y = element_text(size = 10))

    ggsave(paste0(contrast_viz_dir, "top_tfs.pdf"), p3, width = 10, height = 10)
  }
}

# 4. Generate integrated sample-level pathway activity scores for visualization
cat("Running sample-level pathway activities. This may take time...\n")

# Prepare counts for scoring - try with transposed orientation (samples as columns)
mat_counts <- as.matrix(norm_counts)

# Get overlapping genes
gene_symbols_in_data <- rownames(mat_counts)
network_genes <- unique(custom_net$target)
overlap_genes <- intersect(gene_symbols_in_data, network_genes)

# Print diagnostic information
cat("Expression data has", length(gene_symbols_in_data), "genes\n")
cat("Network has", length(network_genes), "genes\n")
cat("Overlap:", length(overlap_genes), "genes\n")

# Filter matrix and network to only include overlapping genes
net_filtered <- custom_net[custom_net$target %in% overlap_genes, ]

# Count genes per pathway
pathway_counts <- table(net_filtered$source)
cat("After filtering, pathways have between", min(pathway_counts), "and", max(pathway_counts), "genes\n")
cat("Number of pathways with ≥1 gene:", sum(pathway_counts >= 1), "\n")
cat("Number of pathways with ≥2 genes:", sum(pathway_counts >= 2), "\n")
cat("Number of pathways with ≥3 genes:", sum(pathway_counts >= 3), "\n")

tryCatch({
  # First attempt with full network
  sample_activities <- decoupleR::run_wmean(
    mat = mat_counts,
    network = net_filtered,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = 1
  )

  # Save results
  write.csv(sample_activities,
           file = paste0(output_dir, "sample_pathway_activities.csv"),
           row.names = FALSE)

  # Create heatmap of sample-level pathway activities
  # Get top variable pathways
  top_var_paths <- sample_activities %>%
    dplyr::filter(statistic == "norm_wmean") %>%
    dplyr::group_by(source) %>%
    dplyr::summarize(variance = var(score)) %>%
    dplyr::top_n(30, variance) %>%
    dplyr::pull(source)

  # Extract scores for heatmap
  path_matrix <- sample_activities %>%
    dplyr::filter(statistic == "norm_wmean" & source %in% top_var_paths) %>%
    tidyr::pivot_wider(id_cols = "source", names_from = "condition", values_from = "score") %>%
    tibble::column_to_rownames("source") %>%
    as.matrix()

  # Create annotation for samples
  sample_anno <- meta_glial %>%
    dplyr::select(Sex, Dx_OUD) %>%
    as.data.frame()

  # Create heatmap
  pdf(paste0(output_dir, "sample_pathway_heatmap.pdf"), width = 12, height = 10)
  pheatmap(path_matrix,
          annotation_col = sample_anno,
          show_colnames = FALSE,
          main = "Sample-level Pathway Activities",
          color = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()

  cat("Sample-level pathway analysis complete! Results saved.\n")
}, error = function(e) {
  cat("ERROR in sample-level pathway analysis:", conditionMessage(e), "\n")

  # If the above fails, try with minimal network
  tryCatch({
    # Create a minimal network with only the top 50 pathways
    top_pathways <- names(sort(pathway_counts, decreasing = TRUE)[1:50])
    minimal_net <- net_filtered[net_filtered$source %in% top_pathways, ]

    # Try with smaller feature set
    sample_activities <- decoupleR::run_wmean(
      mat = mat_counts,
      network = minimal_net,
      .source = "source",
      .target = "target",
      .mor = "mor"
    )

    write.csv(sample_activities,
              file = paste0(output_dir, "sample_pathway_activities_minimal.csv"),
              row.names = FALSE)
    cat("Sample-level analysis with minimal network succeeded!\n")
  }, error = function(e2) {
    cat("Second attempt failed:", conditionMessage(e2), "\n")
    cat("Proceeding with contrast-level results only.\n")
  })
})

# Add visualization for the successful minimal network analysis - WITH FIXES
tryCatch({
  # Check if minimal network results exist
  if(file.exists(paste0(output_dir, "sample_pathway_activities_minimal.csv"))) {
    # Load the saved sample activities
    sample_activities <- read.csv(paste0(output_dir, "sample_pathway_activities_minimal.csv"))

    # Simple summary plots that won't fail with the 'gpar' error
    cat("Creating robust summary visualizations...\n")

    # 1. Simple barplot of most variable pathways
    if(nrow(sample_activities) > 0) {
      # Get pathway variability
      path_var <- sample_activities %>%
        dplyr::filter(statistic == "norm_wmean") %>%
        dplyr::group_by(source) %>%
        dplyr::summarize(
          mean_score = mean(score, na.rm = TRUE),
          var_score = var(score, na.rm = TRUE)
        ) %>%
        dplyr::arrange(desc(var_score)) %>%
        dplyr::slice_head(n = 20)

      if(nrow(path_var) > 0) {
        # Create bar plot of most variable pathways
        p1 <- ggplot(path_var, aes(x = reorder(source, var_score), y = var_score)) +
          geom_col(fill = "steelblue") +
          coord_flip() +
          labs(title = "Most Variable Pathways Across Samples",
               x = "Pathway", y = "Variance") +
          theme_minimal()

        ggsave(paste0(output_dir, "pathway_variance_barplot.pdf"), p1, width = 10, height = 8)

        # 2. Mean pathway activity
        p2 <- ggplot(path_var, aes(x = reorder(source, mean_score), y = mean_score)) +
          geom_col(fill = ifelse(path_var$mean_score >= 0, "red", "blue")) +
          coord_flip() +
          labs(title = "Mean Pathway Activity",
               x = "Pathway", y = "Mean Activity Score") +
          theme_minimal()

        ggsave(paste0(output_dir, "pathway_mean_barplot.pdf"), p2, width = 10, height = 8)

        # 3. Create a faceted boxplot of pathway activities by condition
        path_data <- sample_activities %>%
          dplyr::filter(statistic == "norm_wmean" & source %in% path_var$source)

        # Add metadata for faceting
        meta_minimal <- data.frame(
          condition = colnames(norm_counts),
          Sex = meta_glial$Sex,
          Dx_OUD = meta_glial$Dx_OUD
        )

        path_data <- path_data %>%
          dplyr::left_join(meta_minimal, by = "condition")

        # Boxplot of top pathways by experimental group
        p3 <- ggplot(path_data, aes(x = Dx_OUD, y = score, fill = Sex)) +
          geom_boxplot() +
          facet_wrap(~source, scales = "free_y", ncol = 4) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(size = 8)
          ) +
          labs(title = "Pathway Activities by Condition", y = "Activity Score")

        ggsave(paste0(output_dir, "pathway_boxplots.pdf"), p3, width = 12, height = 10)

        cat("Successfully created alternative visualizations!\n")
      } else {
        cat("Warning: No variable pathways found for visualization\n")
      }
    } else {
      cat("Warning: No data available in sample_activities\n")
    }
  } else {
    cat("No minimal network results file found, skipping visualization\n")
  }
}, error = function(e) {
  cat("Error creating visualizations:", conditionMessage(e), "\n")

  # Final fallback to ensure at least one visualization works
  tryCatch({
    sample_act <- read.csv(paste0(output_dir, "sample_pathway_activities_minimal.csv"))

    # Count number of samples per pathway
    pathway_counts <- sample_act %>%
      dplyr::filter(statistic == "norm_wmean") %>%
      dplyr::group_by(source) %>%
      dplyr::summarize(count = n())

    # Simple dotplot - extremely robust
    p_simple <- ggplot(pathway_counts, aes(x = count, y = reorder(source, count))) +
      geom_point(size = 3) +
      labs(title = "Number of Samples per Pathway",
           x = "Count", y = "Pathway") +
      theme_minimal()

    ggsave(paste0(output_dir, "absolute_minimal_plot.pdf"), p_simple, width = 8, height = 10)
    cat("Created absolute minimal plot as last resort\n")
  }, error = function(e2) {
    cat("All visualization attempts failed\n")
  })
})

# 5. Generate summary report of key findings for each contrast
sink(paste0(output_dir, "summary_report.txt"))

cat("# decoupleR Analysis Results Summary\n\n")
cat("Analysis performed on glial cell data from OUD study\n\n")

for(contrast_name in names(all_methods_results)) {
  cat("## Contrast:", contrast_name, "\n\n")

  # Get pathway results and check column names
  pathways_result <- all_methods_results[[contrast_name]]$pathways
  cat("Available columns in pathway results:", paste(colnames(pathways_result), collapse=", "), "\n\n")

  # Check if padj exists, otherwise use p_value
  p_value_col <- if("padj" %in% colnames(pathways_result)) "padj" else "p_value"

  # Top upregulated pathways
  cat("### Top Upregulated Pathways:\n")
  top_up <- pathways_result %>%
    filter(statistic == "norm_wmean") %>%
    arrange(desc(score)) %>%
    head(10)

  for(i in 1:nrow(top_up)) {
    cat(sprintf("%d. %s (score = %.3f, p-value = %.3e)\n",
                i, top_up$source[i], top_up$score[i], top_up[[p_value_col]][i]))
  }
  cat("\n")

  # Top downregulated pathways
  cat("### Top Downregulated Pathways:\n")
  top_down <- pathways_result %>%
    filter(statistic == "norm_wmean") %>%
    arrange(score) %>%
    head(10)

  for(i in 1:nrow(top_down)) {
    cat(sprintf("%d. %s (score = %.3f, p-value = %.3e)\n",
                i, top_down$source[i], top_down$score[i], top_down[[p_value_col]][i]))
  }
  cat("\n")

  # Top significant TFs
  cat("### Key Transcription Factors:\n")

  # Check TF results and get columns
  tf_result <- all_methods_results[[contrast_name]]$tfs

  if(nrow(tf_result) > 0) {
    cat("Available columns in TF results:", paste(colnames(tf_result), collapse=", "), "\n\n")

    # Check which statistic value is present
    tf_stat_value <- if("norm_ulm" %in% unique(tf_result$statistic)) "norm_ulm" else "ulm"
    cat("Using statistic value:", tf_stat_value, "\n\n")

    # Check if padj exists, otherwise use p_value for TFs
    tf_p_value_col <- if("padj" %in% colnames(tf_result)) "padj" else "p_value"

    top_tfs <- tf_result %>%
      filter(statistic == tf_stat_value) %>%
      arrange(desc(abs(score))) %>%
      head(10)

    for(i in 1:nrow(top_tfs)) {
      cat(sprintf("%d. %s (score = %.3f, p-value = %.3e)\n",
                  i, top_tfs$source[i], top_tfs$score[i], top_tfs[[tf_p_value_col]][i]))
    }
  } else {
    cat("No transcription factor results available for this contrast.\n")
  }
  cat("\n\n")
}

sink()

cat("decoupleR analysis complete! Results saved to:", output_dir, "\n")