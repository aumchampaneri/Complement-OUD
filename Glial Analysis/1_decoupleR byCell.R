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
library(stringr)

# Check if required packages are installed and install if needed
required_packages <- c("decoupleR", "dplyr", "tidyr", "ggplot2", "ComplexHeatmap",
                      "circlize", "org.Hs.eg.db", "msigdbr", "pheatmap", "limma", "tibble", "stringr")

missing_packages <- required_packages[!requireNamespace(required_packages, quietly = TRUE)]
if(length(missing_packages) > 0) {
  if(!require(BiocManager, quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(missing_packages)
}

# Set seed for reproducibility
set.seed(123)

# Main output directory
base_output_dir <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/Glial Analysis/decoupleR_results/'
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

# Load metadata once
meta <- read.csv('/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/meta.csv', row.names=1)

# Get list of cell types from the directory structure
# Handle spaces in directory path properly with normalizePath
celltype_dir <- normalizePath('/Users/aumchampaneri/PycharmProjects/Complement-OUD/Glial Analysis/edgeR Glial Per_Celltype/', mustWork = FALSE)
cell_types <- list.dirs(celltype_dir, full.names = FALSE, recursive = FALSE)

# Remove any empty strings or hidden directories
cell_types <- cell_types[cell_types != "" & !startsWith(cell_types, ".")]

# If no cell types found, provide error message with instructions
if(length(cell_types) == 0) {
  cat("No cell type directories found. Please check the path:", celltype_dir, "\n")
  cat("The script expects a directory structure like:\n")
  cat("/path/to/Per_Celltype/Astrocytes/edgeR_results_*.csv\n")
  cat("/path/to/Per_Celltype/Microglia/edgeR_results_*.csv\n")
  cat("etc.\n")

  # Try alternative location
  alt_dir <- '/Users/aumchampaneri/PycharmProjects/Complement-OUD/Glial Analysis/edgeR Per_Celltype'
  if(dir.exists(alt_dir)) {
    celltype_dir <- alt_dir
    cat("Found cell types in alternative location:", alt_dir, "\n")
    cell_types <- list.dirs(celltype_dir, full.names = FALSE, recursive = FALSE)
    cell_types <- cell_types[cell_types != "" & !startsWith(cell_types, ".")]
  }

  # If still no cell types found, check one directory level up and filter
  if(length(cell_types) == 0) {
    parent_dir <- dirname(celltype_dir)
    all_dirs <- list.dirs(parent_dir, full.names = FALSE, recursive = FALSE)
    cell_candidates <- grep("Per_Celltype|Celltype", all_dirs, value = TRUE)

    if(length(cell_candidates) > 0) {
      for(cand in cell_candidates) {
        potential_path <- file.path(parent_dir, cand)
        potential_types <- list.dirs(potential_path, full.names = FALSE, recursive = FALSE)
        if(length(potential_types) > 0) {
          celltype_dir <- potential_path
          cell_types <- potential_types[potential_types != "" & !startsWith(potential_types, ".")]
          cat("Found cell types in:", potential_path, "\n")
          break
        }
      }
    }
  }

  # If all automatic detection fails, use manual fallback
  if(length(cell_types) == 0) {
    # Manually specify cell types
    cell_types <- c("Astrocytes", "Microglia", "Oligodendrocytes", "OPCs")
    cat("Using manually specified cell types:", paste(cell_types, collapse=", "), "\n")
  }
}

cat("Processing the following cell types:", paste(cell_types, collapse=", "), "\n\n")

# Get resources once to avoid repetition
# a. PROGENy for pathway activities
progeny_human <- get_progeny(organism = "human", top = 500)

# Fix PROGENy column names
if ("weight" %in% colnames(progeny_human)) {
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

# Add likelihood column to custom_net if missing
if(!"likelihood" %in% colnames(custom_net)) {
  custom_net <- custom_net %>%
    mutate(likelihood = 1)
}

# Prepare network for sample-level pathway analysis
minimal_net <- custom_net %>%
  filter(grepl("COMPLEMENT|MICROGLIA|ASTRO|OLIGO|INFLAM|IMMUNE|NEURON", source, ignore.case = TRUE))

cat("Created custom networks with", nrow(custom_net), "total edges and",
    nrow(minimal_net), "focused edges for minimal network\n\n")

# Create a summary report file
summary_file <- paste0(base_output_dir, "summary_report.txt")
file.create(summary_file)
writeLines("# Cell Type decoupleR Analysis Summary\n", summary_file)

# Create a combined output file for all TFs/pathways across cell types
combined_file <- paste0(base_output_dir, "combined_results.csv")
combined_df <- data.frame()
first_write <- TRUE

# Process each cell type
for(cell_type in cell_types) {
  cat("\nProcessing", cell_type, "...\n")

  # Create output directory for this cell type
  output_dir <- paste0(base_output_dir, cell_type, "/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Set up contrast metadata
  cell_type_fixed <- gsub(" ", "_", cell_type)  # Fix spaces in cell type name if any

  # Look for edgeR results for this cell type
  result_files <- list.files(
    path = file.path(celltype_dir, cell_type),
    pattern = "edgeR_results_.*\\.csv$",
    full.names = TRUE
  )

  # If no files found, try alternative cell type name formats
  if(length(result_files) == 0) {
    # Try looking for results with alternative naming
    alternative_pattern <- paste0(".*", gsub("s$", "", cell_type), ".*\\.csv$")
    result_files <- list.files(
      path = celltype_dir,
      pattern = alternative_pattern,
      recursive = TRUE,
      full.names = TRUE
    )
  }

  # Skip if no results found for this cell type
  if(length(result_files) == 0) {
    cat("No edgeR results found for", cell_type, "- skipping\n")
    next
  }

  cat("Found", length(result_files), "results file(s) for", cell_type, "\n")

  # Add to summary report
  cat("\n## Cell Type:", cell_type, "\n\n", file = summary_file, append = TRUE)

  # Process each contrast file
  all_methods_results <- list()

  for(file_path in result_files) {
    # Extract contrast name from file name
    contrast_name <- gsub("edgeR_results_|\\.csv$", "", basename(file_path))
    cat("Processing contrast:", contrast_name, "\n")

    # Add to summary report
    cat("## Contrast:", contrast_name, "\n\n", file = summary_file, append = TRUE)

    # Load contrast data
    edger_results <- read.csv(file_path)

    # Debug info about file structure
    cat("Preview of input file columns:", paste(colnames(edger_results), collapse=", "), "\n")

  # Ensure we have proper gene column
  gene_col <- intersect(c("GeneSymbol", "Symbol", "Gene", "gene_symbol", "SYMBOL", "gene", "gene_name", "X"),
                       colnames(edger_results))[1]

  if(is.na(gene_col)) {
    cat("Warning: No standard gene symbol column found in", file_path, "\n")
    # Debug - show available columns
    cat("Available columns:", paste(colnames(edger_results), collapse=", "), "\n")
    cat("First few rows of data:\n")
    print(head(edger_results[, 1:min(5, ncol(edger_results))]))

    # Try first column if it looks like gene symbols
    first_col <- colnames(edger_results)[1]
    if(!is.null(first_col) &&
       all(grepl("^[A-Za-z0-9]", as.character(edger_results[[first_col]])))) {
      gene_col <- first_col
      cat("Using first column as gene symbols:", gene_col, "\n")
    } else {
      cat("Error: Cannot process file without gene symbols\n")
      next
    }
  }

  # Prepare input matrix for decoupleR
  # Use logFC column for decoupleR input - check for different naming conventions
  logfc_col <- intersect(c("logFC", "log2FoldChange", "log2FC"), colnames(edger_results))[1]
  pval_col <- intersect(c("FDR", "padj", "PValue", "pvalue", "P.Value"), colnames(edger_results))[1]

  if(is.na(logfc_col)) {
    cat("Error: No logFC column found in", file_path, "\n")
    next
  }

  # Create input matrix
  input <- edger_results %>%
    dplyr::select(all_of(c(gene_col, logfc_col))) %>%
    dplyr::rename(gene = all_of(gene_col), score = all_of(logfc_col)) %>%
    dplyr::filter(!is.na(score))

  # MOVED: Put diagnostic output AFTER creating the input matrix
  cat("Found", nrow(input), "genes with score values\n")
  cat("Sample genes:", paste(head(input$gene), collapse=", "), "\n")

    # Prepare input matrix for decoupleR
    # Use logFC column for decoupleR input - check for different naming conventions
    logfc_col <- intersect(c("logFC", "log2FoldChange", "log2FC"), colnames(edger_results))[1]
    pval_col <- intersect(c("FDR", "padj", "PValue", "pvalue", "P.Value"), colnames(edger_results))[1]

    if(is.na(logfc_col)) {
      cat("Error: No logFC column found in", file_path, "\n")
      next
    }

    # Create input matrix
    input <- edger_results %>%
      dplyr::select(all_of(c(gene_col, logfc_col))) %>%
      dplyr::rename(gene = all_of(gene_col), score = all_of(logfc_col)) %>%
      dplyr::filter(!is.na(score))

    # Apply standard preprocessing
    rownames(input) <- input$gene
    input <- input[!is.na(input$gene) & input$gene != "",]
    input <- input[!duplicated(input$gene),]

    # Debug info
    cat("Number of genes in input matrix:", nrow(input), "\n")
    cat("First few genes:", paste(head(input$gene), collapse=", "), "\n")

    # Log transformation for better distribution if values are very large or small
    if(max(abs(input$score), na.rm = TRUE) > 100) {
      input$score <- sign(input$score) * log2(abs(input$score) + 1)
      cat("Applied log transformation to large fold changes\n")
    }

    # Create an output object to store all results for this contrast
    contrast_results <- list()

    # Get p-value column for annotation if available
    if(!is.na(pval_col)) {
      p_values <- edger_results %>%
        dplyr::select(all_of(c(gene_col, pval_col))) %>%
        dplyr::rename(gene = all_of(gene_col), p_value = all_of(pval_col))

      rownames(p_values) <- p_values$gene
      p_values <- p_values[rownames(input), "p_value", drop = FALSE]
    }

    # Run PROGENy analysis with fallback to smaller minsize
    cat("Running PROGENy for pathway analysis...\n")
    progeny_success <- FALSE

    # Try with progressively smaller minsize values
    for(try_size in c(5, 3, 2, 1)) {
      tryCatch({
        progeny_activities <- run_wmean(
          mat = input,
          network = progeny_human,
          .source = "source",
          .target = "target",
          .mor = "mor",
          .likelihood = "likelihood",
          times = 100,
          minsize = try_size
        )
        cat("PROGENy analysis successful with minsize =", try_size, "\n")
        progeny_success <- TRUE
        break
      }, error = function(e) {
        cat("PROGENy attempt failed with minsize =", try_size, ":", conditionMessage(e), "\n")
      })
    }

    # FIXED: Remove duplicate check
    if(!progeny_success) {
      cat("ERROR: PROGENy analysis failed. Skipping this contrast.\n")
      next
    }

    # Run DoRothEA analysis with fallback to smaller minsize
    cat("Running DoRothEA for transcription factor analysis...\n")
    dorothea_success <- FALSE

    for(try_size in c(5, 3, 2, 1)) {
      tryCatch({
        dorothea_activities <- run_ulm(
          mat = input,
          network = dorothea_human,
          .source = "source",
          .target = "target",
          .mor = "mor",
          .likelihood = "likelihood",
          minsize = try_size
        )
        cat("DoRothEA analysis successful with minsize =", try_size, "\n")
        dorothea_success <- TRUE
        break
      }, error = function(e) {
        cat("DoRothEA attempt failed with minsize =", try_size, ":", conditionMessage(e), "\n")
      })
    }

    if(!dorothea_success) {
      cat("ERROR: DoRothEA analysis failed. Trying to continue without TF analysis.\n")
      dorothea_activities <- data.frame(source=character(), target=character(), score=numeric(), statistic=character())
    }

    # Run custom gene set analysis with fallback to smaller minsize
    cat("Running custom pathway analysis...\n")
    custom_success <- FALSE

    for(try_size in c(5, 3, 2, 1)) {
      tryCatch({
        custom_activities <- run_wmean(
          mat = input,
          network = custom_net,
          .source = "source",
          .target = "target",
          .mor = "mor",
          .likelihood = "likelihood",
          times = 100,
          minsize = try_size
        )
        cat("Custom pathway analysis successful with minsize =", try_size, "\n")
        custom_success <- TRUE
        break
      }, error = function(e) {
        cat("Custom pathway analysis failed with minsize =", try_size, ":", conditionMessage(e), "\n")
      })
    }

    if(!custom_success) {
      cat("ERROR: Custom pathway analysis failed. Skipping this contrast.\n")
      next
    }

    # Combine results and add condition column
    progeny_activities$condition <- contrast_name
    dorothea_activities$condition <- contrast_name
    custom_activities$condition <- contrast_name

    # Process and store results
    contrast_results$progeny <- progeny_activities
    contrast_results$tfs <- dorothea_activities
    contrast_results$pathways <- custom_activities

    # Store in the main results list
    all_methods_results[[contrast_name]] <- contrast_results

    # Only store normalized activities for reporting
    norm_activities <- progeny_activities %>%
      dplyr::filter(statistic == "norm_wmean") %>%
      dplyr::arrange(desc(abs(score)))

    # Get p-value column from the results
    p_value_col <- "p_value"

    # Add p-value if missing
    if(!"p_value" %in% colnames(custom_activities)) {
      custom_activities$p_value <- NA
    }

    # Add to summary report
    cat("Available columns in pathway results:",
        paste(colnames(custom_activities), collapse = ", "), "\n\n",
        file = summary_file, append = TRUE)

    # Extract and format pathway results for combined output
    pathways_result <- custom_activities %>%
      dplyr::filter(statistic == "norm_wmean") %>%
      dplyr::mutate(
        cell_type = cell_type,
        contrast = contrast_name
      )

    # Add to combined output
    combined_entry <- pathways_result %>%
      dplyr::select(cell_type, contrast, source, score, all_of(p_value_col)) %>%
      dplyr::rename(pathway = source)

    if(nrow(combined_df) == 0) {
      combined_df <- combined_entry
    } else {
      combined_df <- bind_rows(combined_df, combined_entry)
    }

    # Write to combined CSV file
    if(first_write) {
      write.csv(combined_entry, combined_file, row.names = FALSE)
      first_write <- FALSE
    } else {
      write.table(combined_entry, combined_file, sep = ",",
                  append = TRUE, row.names = FALSE, col.names = FALSE)
    }

    # Create output visualizations for each contrast

    # TF Activity Plots
    if(dorothea_success) {
      tf_result <- dorothea_activities %>%
        dplyr::filter(statistic == "ulm") %>%
        dplyr::arrange(desc(abs(score))) %>%
        dplyr::slice_head(n = 30)

      if(nrow(tf_result) > 0) {
        tf_plot <- ggplot(tf_result, aes(x = reorder(source, score), y = score, fill = score > 0)) +
          geom_col() +
          coord_flip() +
          scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                            labels = c("TRUE" = "Up", "FALSE" = "Down"),
                            name = "Direction") +
          labs(title = paste0("Top TFs - ", contrast_name),
              x = "Transcription Factor",
              y = "Activity Score") +
          theme_minimal() +
          theme(
            axis.text.y = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.position = "bottom"
          )

        # Save plot
        ggsave(paste0(output_dir, "tf_activity_", contrast_name, ".pdf"), tf_plot, width = 10, height = 8)
      }
    }

    # Pathway Activity Plots
    pathway_result <- custom_activities %>%
      dplyr::filter(statistic == "norm_wmean") %>%
      dplyr::arrange(desc(abs(score))) %>%
      dplyr::slice_head(n = 30)

    if(nrow(pathway_result) > 0) {
      # Format pathway names for better readability
      pathway_result$source <- gsub("_", " ", pathway_result$source)
      pathway_result$wrapped_source <- str_wrap(pathway_result$source, width = 35)

      # Create plot
      pathway_plot <- ggplot(pathway_result, aes(x = reorder(wrapped_source, score), y = score, fill = score > 0)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                          labels = c("TRUE" = "Up", "FALSE" = "Down"),
                          name = "Direction") +
        labs(title = paste0("Top Pathways - ", contrast_name),
             x = "Pathway",
             y = "Activity Score") +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 9),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom"
        )

      # Save plot
      ggsave(paste0(output_dir, "pathway_activity_", contrast_name, ".pdf"), pathway_plot, width = 12, height = 10)
    }

    # Add top pathways to summary report
    cat("### Top Upregulated Pathways:\n", file = summary_file, append = TRUE)
    top_up <- pathways_result %>%
      filter(statistic == "norm_wmean") %>%
      filter(score > 0) %>%
      top_n(10, score) %>%
      arrange(desc(score))

    for(i in 1:nrow(top_up)) {
      cat("- ", top_up$source[i], ": score =", round(top_up$score[i], 3),
          ", p-value =", format(top_up[[p_value_col]][i], scientific = TRUE, digits = 3), "\n",
          file = summary_file, append = TRUE)
    }
    cat("\n", file = summary_file, append = TRUE)

    # Top downregulated pathways
    cat("### Top Downregulated Pathways:\n", file = summary_file, append = TRUE)
    top_down <- pathways_result %>%
      filter(statistic == "norm_wmean") %>%
      filter(score < 0) %>%
      top_n(10, abs(score)) %>%
      arrange(score)

    for(i in 1:nrow(top_down)) {
      cat("- ", top_down$source[i], ": score =", round(top_down$score[i], 3),
          ", p-value =", format(top_down[[p_value_col]][i], scientific = TRUE, digits = 3), "\n",
          file = summary_file, append = TRUE)
    }
    cat("\n", file = summary_file, append = TRUE)

    # Top significant TFs
    cat("### Key Transcription Factors:\n", file = summary_file, append = TRUE)

    # Check TF results and get columns
    tf_result <- all_methods_results[[contrast_name]]$tfs

    if(!is.null(tf_result) && nrow(tf_result) > 0) {
      cat("Available columns in TF results:",
          paste(colnames(tf_result), collapse = ", "), "\n\n",
          file = summary_file, append = TRUE)

      tf_stat_value <- if("norm_ulm" %in% unique(tf_result$statistic)) "norm_ulm" else "ulm"
      cat("Using statistic value:", tf_stat_value, "\n\n",
          file = summary_file, append = TRUE)

      top_tfs <- tf_result %>%
        filter(statistic == tf_stat_value) %>%
        top_n(15, abs(score)) %>%
        arrange(desc(abs(score)))

      for(i in 1:min(nrow(top_tfs), 15)) {
        direction <- ifelse(top_tfs$score[i] > 0, "up", "down")
        cat(i, ". ", top_tfs$source[i], " (score = ", round(top_tfs$score[i], 3),
            ", ", direction, "regulated)\n",
            file = summary_file, append = TRUE)
      }
    } else {
      cat("No significant transcription factors found for this contrast.\n",
          file = summary_file, append = TRUE)
    }
    cat("\n\n", file = summary_file, append = TRUE)
  }  # End of contrast loop

  # Create heatmap of pathway activities across all contrasts
  if(length(all_methods_results) > 1) {
    # Create a matrix of pathway scores for all contrasts
    pathway_matrix <- lapply(names(all_methods_results), function(contrast) {
      pathway_df <- all_methods_results[[contrast]]$pathways %>%
        filter(statistic == "norm_wmean") %>%
        select(source, score)

      # Set rownames to pathway names for this contrast column
      pathway_scores <- pathway_df$score
      names(pathway_scores) <- pathway_df$source
      return(pathway_scores)
    })

    # Convert list to matrix
    pathway_matrix <- do.call(cbind, pathway_matrix)
    colnames(pathway_matrix) <- names(all_methods_results)

    # Filter for most variable pathways
    if(nrow(pathway_matrix) > 20) {
      path_vars <- apply(pathway_matrix, 1, var, na.rm = TRUE)
      top_paths <- names(sort(path_vars, decreasing = TRUE))[1:20]
      pathway_matrix <- pathway_matrix[top_paths, , drop = FALSE]
    }

    # Create a heatmap
    if(nrow(pathway_matrix) > 0 && ncol(pathway_matrix) > 1) {
      # Fix pathway names for display
      rownames(pathway_matrix) <- gsub("_", " ", rownames(pathway_matrix))

      # Create colormap
      col_fun <- colorRamp2(c(min(pathway_matrix, na.rm = TRUE),
                             0,
                             max(pathway_matrix, na.rm = TRUE)),
                          c("blue", "white", "red"))

# Generate heatmap
      ht <- Heatmap(pathway_matrix,
                   name = "Score",
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   col = col_fun,
                   row_names_gp = gpar(fontsize = 10),
                   column_names_gp = gpar(fontsize = 10),
                   heatmap_legend_param = list(title = "Activity Score"))

      # Save heatmap
      pdf(paste0(output_dir, "pathway_heatmap.pdf"), width = 10, height = 12)
      draw(ht)
      dev.off()
    }
  }

  # Collect sample-level activities for variability analysis if available
  if(exists("sample_activities")) {
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
      # Create bar plot of most variable pathways with wrapped labels
      p1 <- ggplot(path_var, aes(x = reorder(source, var_score), y = var_score)) +
        geom_col(fill = "steelblue") +
        coord_flip() +
        labs(title = paste0("Most Variable Pathways - ", cell_type),
             x = "Pathway", y = "Variance") +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5)
        ) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 40))

      ggsave(paste0(output_dir, "pathway_variance_barplot.pdf"), p1, width = 12, height = 8)
    }
  }

  # Generate sample count analysis if multiple contrasts
  if(length(all_methods_results) > 1) {
    # Count how many samples each pathway appears in the top 50
    pathway_counts <- lapply(names(all_methods_results), function(contrast) {
      paths <- all_methods_results[[contrast]]$pathways %>%
        filter(statistic == "norm_wmean") %>%
        arrange(desc(abs(score))) %>%
        slice_head(n = 50) %>%
        pull(source)

      return(data.frame(source = paths, contrast = contrast))
    })

    # Combine all contrasts
    pathway_counts <- bind_rows(pathway_counts)

    # Count occurrences of each pathway
    pathway_counts <- pathway_counts %>%
      group_by(source) %>%
      summarize(count = n()) %>%
      arrange(desc(count)) %>%
      slice_head(n = 30)

    # Format pathway names
    pathway_counts$source <- gsub("_", " ", pathway_counts$source)
    pathway_counts$wrapped_source <- str_wrap(pathway_counts$source, width = 40)

    # Create count plot
    if(nrow(pathway_counts) > 0) {
      p2 <- ggplot(pathway_counts, aes(x = reorder(wrapped_source, count), y = count)) +
        geom_col(fill = "darkgreen") +
        coord_flip() +
        labs(title = paste0("Pathways by Contrast Count - ", cell_type),
             x = "Pathway", y = "Number of Contrasts") +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5)
        )

      # Save plot
      ggsave(paste0(output_dir, "pathway_counts.pdf"), p2, width = 12, height = 8)

      # Simple dotplot - extremely robust
      p_simple <- ggplot(pathway_counts, aes(x = count, y = reorder(wrapped_source, count))) +
        geom_point(size = 3, color = "darkblue") +
        theme_minimal() +
        labs(title = paste0("Pathways by Sample Count - ", cell_type),
             x = "Number of Contrasts", y = "Pathway") +
        theme(
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5)
        )

      ggsave(paste0(output_dir, "pathway_dotplot.pdf"), p_simple, width = 10, height = 8)
    }
  }

  # Add closing message for this cell type
  cat("decoupleR analysis for", cell_type, "complete! Results saved to:", output_dir, "\n")
}  # End of cell type loop

# Create combined heatmap of top pathways across all cell types if multiple cell types
if(length(cell_types) > 1 && !first_write) {
  # Read the combined results for heatmap creation
  combined_results <- read.csv(combined_file)

  # Get top pathways by absolute score
  top_pathways <- combined_results %>%
    group_by(pathway) %>%
    summarize(mean_abs_score = mean(abs(score), na.rm = TRUE)) %>%
    top_n(30, mean_abs_score) %>%
    pull(pathway)

  # Filter results to these top pathways
  filtered_results <- combined_results %>%
    filter(pathway %in% top_pathways)

  # Create a matrix for heatmap
  matrix_data <- filtered_results %>%
    pivot_wider(
      id_cols = pathway,
      names_from = c(cell_type, contrast),
      values_from = score,
      names_sep = "_"
    )

  # Convert to matrix
  pathway_matrix <- as.matrix(matrix_data[, -1])
  rownames(pathway_matrix) <- matrix_data$pathway

  # Format pathway names
  rownames(pathway_matrix) <- gsub("_", " ", rownames(pathway_matrix))

  # Create heatmap if data exists
  if(nrow(pathway_matrix) > 0 && ncol(pathway_matrix) > 1) {
    col_fun <- colorRamp2(c(min(pathway_matrix, na.rm = TRUE),
                           0,
                           max(pathway_matrix, na.rm = TRUE)),
                        c("blue", "white", "red"))

    # Generate combined heatmap
    ht <- Heatmap(pathway_matrix,
                 name = "Score",
                 cluster_rows = TRUE,
                 cluster_columns = TRUE,
                 show_row_names = TRUE,
                 show_column_names = TRUE,
                 col = col_fun,
                 row_names_gp = gpar(fontsize = 10),
                 column_names_gp = gpar(fontsize = 8),
                 column_names_rot = 45,
                 heatmap_legend_param = list(title = "Activity Score"))

    # Save heatmap
    pdf(paste0(base_output_dir, "combined_pathway_heatmap.pdf"), width = 14, height = 12)
    draw(ht)
    dev.off()

    # Create interactive heatmap version if ComplexHeatmap has interactive capabilities
    # Note: This requires the InteractiveComplexHeatmap package
    if(requireNamespace("InteractiveComplexHeatmap", quietly = TRUE)) {
      InteractiveComplexHeatmap::saveHeatmap(ht, paste0(base_output_dir, "interactive_heatmap.html"))
    }
  }

  # Create bubble plot for combined visualization
  filtered_results <- filtered_results %>%
    mutate(
      significance = -log10(p_value),
      direction = ifelse(score > 0, "Up", "Down"),
      contrast_label = paste(cell_type, contrast, sep = "_")
    )

  # Replace Inf values
  filtered_results$significance[is.infinite(filtered_results$significance)] <-
    max(filtered_results$significance[!is.infinite(filtered_results$significance)], na.rm = TRUE) * 1.1

  # Create bubble plot
  bubble_plot <- ggplot(filtered_results,
                      aes(x = contrast_label, y = pathway,
                          size = abs(score), color = score,
                          alpha = significance)) +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_alpha(range = c(0.3, 1)) +
    scale_size(range = c(1, 8)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = "Pathway Activities Across Cell Types and Contrasts",
      x = "",
      y = "",
      size = "Absolute\nScore",
      color = "Score",
      alpha = "-log10(p-value)"
    )

  # Save bubble plot
  ggsave(paste0(base_output_dir, "combined_bubble_plot.pdf"),
         bubble_plot, width = 14, height = 12)
}

cat("\nAll cell type analyses completed successfully!\n")