# ========================================================================
# ENHANCED CROSS-DATASET COMPLEMENT ANALYSIS
# Comparing Available OUD-Related Datasets with Focus on Complement System
# ========================================================================

// ...existing code until dataset loading...

# ========================================================================
# ENHANCED DATA LOADING AND VALIDATION
# ========================================================================

# Dataset information
datasets <- list(
  GSE151807 = list(
    path = "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE151807",
    tissue = "Nucleus Accumbens",
    treatment = "Morphine vs Saline",
    species = "Mus musculus",
    available = FALSE
  ),
  GSE66351 = list(
    path = "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE66351", 
    tissue = "Striatum",
    treatment = "Cocaine vs Saline",
    species = "Mus musculus",
    available = FALSE
  ),
  GSE239387 = list(
    path = "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE239387",
    tissue = "Nucleus Accumbens", 
    treatment = "Morphine vs Control",
    species = "Mus musculus",
    available = FALSE
  )
)

# Check dataset availability
cat("=== DATASET AVAILABILITY CHECK ===\n")
available_datasets <- list()
de_results <- list()
complement_results <- list()

for(dataset_name in names(datasets)) {
  dataset_info <- datasets[[dataset_name]]
  
  # Check for DE results
  de_file <- file.path(dataset_info$path, "Outputs/02_Differential_Expression", 
                       paste0(dataset_name, "_DE_results.csv"))
  
  # Check for complement results  
  complement_file <- file.path(dataset_info$path, "Outputs/02_Differential_Expression/Data",
                              "complement_DE_results.csv")
  
  if(file.exists(de_file)) {
    cat("✓", dataset_name, "DE results found\n")
    datasets[[dataset_name]]$available <- TRUE
    available_datasets[[dataset_name]] <- dataset_info
    
    # Load DE results
    de_results[[dataset_name]] <- read.csv(de_file)
    cat("  -", nrow(de_results[[dataset_name]]), "genes analyzed\n")
    
    # Load complement results if available
    if(file.exists(complement_file)) {
      complement_results[[dataset_name]] <- read.csv(complement_file)
      cat("  -", nrow(complement_results[[dataset_name]]), "complement genes\n")
    }
    
  } else {
    cat("✗", dataset_name, "DE results not found\n")
    cat("  Expected:", de_file, "\n")
  }
}

cat("\nAvailable datasets:", length(available_datasets), "out of", length(datasets), "\n\n")

# ========================================================================
# SINGLE DATASET ANALYSIS (GSE239387 FOCUS)
# ========================================================================

if("GSE239387" %in% names(available_datasets)) {
  cat("=== DETAILED GSE239387 ANALYSIS ===\n")
  
  gse239387_de <- de_results[["GSE239387"]]
  
  # Basic statistics
  sig_genes <- gse239387_de[gse239387_de$P.Value < 0.05, ]
  up_genes <- sig_genes[sig_genes$logFC > 0, ]
  down_genes <- sig_genes[sig_genes$logFC < 0, ]
  
  cat("Total genes analyzed:", nrow(gse239387_de), "\n")
  cat("Significantly changed genes (p < 0.05):", nrow(sig_genes), "\n")
  cat("  - Upregulated:", nrow(up_genes), "\n")
  cat("  - Downregulated:", nrow(down_genes), "\n")
  
  # Complement gene analysis
  if("GSE239387" %in% names(complement_results)) {
    complement_de <- complement_results[["GSE239387"]]
    sig_complement <- complement_de[complement_de$P.Value < 0.05, ]
    
    cat("\nComplement system analysis:\n")
    cat("Total complement genes:", nrow(complement_de), "\n")
    cat("Significantly changed:", nrow(sig_complement), "\n")
    
    if(nrow(sig_complement) > 0) {
      up_complement <- sig_complement[sig_complement$logFC > 0, ]
      down_complement <- sig_complement[sig_complement$logFC < 0, ]
      
      cat("  - Upregulated complement genes:", nrow(up_complement), "\n")
      cat("  - Downregulated complement genes:", nrow(down_complement), "\n")
      
      # Show top changed complement genes
      cat("\nTop 10 most significantly changed complement genes:\n")
      top_complement <- head(sig_complement[order(sig_complement$P.Value), ], 10)
      for(i in 1:nrow(top_complement)) {
        direction <- ifelse(top_complement$logFC[i] > 0, "↑", "↓")
        cat(sprintf("  %d. %s %s (logFC=%.2f, p=%.2e) [%s]\n", 
                    i, top_complement$Symbol[i], direction,
                    top_complement$logFC[i], top_complement$P.Value[i],
                    top_complement$pathway_class[i]))
      }
      
      # Complement pathway breakdown
      cat("\nComplement pathway analysis:\n")
      pathway_summary <- table(sig_complement$pathway_class)
      for(pathway in names(pathway_summary)) {
        cat("  -", pathway, ":", pathway_summary[[pathway]], "genes\n")
      }
    }
  }
}

# ========================================================================
# PREPARE FOR MULTI-DATASET COMPARISON (WHEN AVAILABLE)  
# ========================================================================

cat("\n=== PREPARING FOR MULTI-DATASET COMPARISON ===\n")

if(length(available_datasets) > 1) {
  cat("Multiple datasets available - performing cross-dataset analysis\n")
  
  # Find common genes across datasets
  common_genes <- Reduce(intersect, lapply(de_results, function(x) x$Symbol))
  cat("Common genes across datasets:", length(common_genes), "\n")
  
  # Cross-dataset effect size comparison
  effect_sizes <- data.frame(Gene = common_genes)
  
  for(dataset in names(de_results)) {
    de_data <- de_results[[dataset]]
    matched_data <- de_data[match(common_genes, de_data$Symbol), ]
    effect_sizes[[paste0(dataset, "_logFC")]] <- matched_data$logFC
    effect_sizes[[paste0(dataset, "_pval")]] <- matched_data$P.Value
  }
  
  # Save cross-dataset comparison
  write.csv(effect_sizes, "Cross_Dataset_Analysis/Gene_Effect_Sizes_Comparison.csv", row.names = FALSE)
  
} else {
  cat("Only one dataset available (GSE239387)\n")
  cat("To enable cross-dataset analysis:\n")
  cat("1. Process GSE151807 differential expression\n")
  cat("2. Process GSE66351 differential expression\n")
  cat("3. Re-run this cross-dataset analysis\n")
}

# ========================================================================
# COMPLEMENT-SPECIFIC META-ANALYSIS FRAMEWORK
# ========================================================================

cat("\n=== COMPLEMENT META-ANALYSIS FRAMEWORK ===\n")

# Prepare complement meta-analysis structure
complement_meta <- list()

for(dataset in names(complement_results)) {
  complement_data <- complement_results[[dataset]]
  
  complement_meta[[dataset]] <- list(
    total_genes = nrow(complement_data),
    significant_genes = nrow(complement_data[complement_data$P.Value < 0.05, ]),
    upregulated = nrow(complement_data[complement_data$P.Value < 0.05 & complement_data$logFC > 0, ]),
    downregulated = nrow(complement_data[complement_data$P.Value < 0.05 & complement_data$logFC < 0, ]),
    tissue = datasets[[dataset]]$tissue,
    treatment = datasets[[dataset]]$treatment
  )
}

# Create complement summary table
if(length(complement_meta) > 0) {
  complement_summary <- data.frame(
    Dataset = names(complement_meta),
    Tissue = sapply(complement_meta, function(x) x$tissue),
    Treatment = sapply(complement_meta, function(x) x$treatment),
    Total_Complement_Genes = sapply(complement_meta, function(x) x$total_genes),
    Significant_Genes = sapply(complement_meta, function(x) x$significant_genes),
    Upregulated = sapply(complement_meta, function(x) x$upregulated),
    Downregulated = sapply(complement_meta, function(x) x$downregulated),
    Percent_Changed = round(sapply(complement_meta, function(x) 
      x$significant_genes / x$total_genes * 100), 1)
  )
  
  write.csv(complement_summary, "Cross_Dataset_Analysis/Complement_Summary_Across_Datasets.csv", row.names = FALSE)
  
  cat("Complement summary across datasets:\n")
  print(complement_summary)
}

# ========================================================================
# DATASET PROCESSING RECOMMENDATIONS
# ========================================================================

cat("\n=== DATASET PROCESSING RECOMMENDATIONS ===\n")

missing_datasets <- setdiff(names(datasets), names(available_datasets))

if(length(missing_datasets) > 0) {
  cat("Missing datasets to process:\n")
  
  for(missing in missing_datasets) {
    cat("\n", missing, ":\n")
    cat("  1. Check if raw data files exist in:", datasets[[missing]]$path, "\n")
    cat("  2. Create differential expression script\n") 
    cat("  3. Run complement gene annotation\n")
    cat("  4. Generate pathway enrichment results\n")
    
    # Check if directory exists
    if(dir.exists(datasets[[missing]]$path)) {
      cat("  ✓ Directory exists\n")
      data_files <- list.files(file.path(datasets[[missing]]$path, "Data"), 
                              pattern = "\\.(txt|csv|tsv)$", full.names = FALSE)
      if(length(data_files) > 0) {
        cat("  ✓ Data files found:", length(data_files), "\n")
        cat("    -", paste(head(data_files, 3), collapse = ", "), 
            ifelse(length(data_files) > 3, "...", ""), "\n")
      } else {
        cat("  ✗ No data files found in Data/ directory\n")
      }
    } else {
      cat("  ✗ Directory does not exist\n")
    }
  }
}

# ========================================================================
# ENHANCED REPORTING AND VISUALIZATION PREPARATION
# ========================================================================

cat("\n=== GENERATING ENHANCED REPORTS ===\n")

# Create detailed summary report
summary_file <- "Cross_Dataset_Analysis/Cross_Dataset_Summary_Report.txt"

cat("", file = summary_file)
cat("========================================================================\n", file = summary_file, append = TRUE)
cat("CROSS-DATASET COMPLEMENT ANALYSIS SUMMARY REPORT\n", file = summary_file, append = TRUE)
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = summary_file, append = TRUE)
cat("========================================================================\n\n", file = summary_file, append = TRUE)

cat("DATASET AVAILABILITY:\n", file = summary_file, append = TRUE)
for(dataset in names(datasets)) {
  status <- ifelse(datasets[[dataset]]$available, "✓ PROCESSED", "✗ MISSING")
  cat(sprintf("- %s: %s (%s, %s)\n", dataset, status, 
              datasets[[dataset]]$tissue, datasets[[dataset]]$treatment), 
      file = summary_file, append = TRUE)
}

if(length(available_datasets) > 0) {
  cat("\nPROCESSED DATASETS SUMMARY:\n", file = summary_file, append = TRUE)
  for(dataset in names(available_datasets)) {
    de_data <- de_results[[dataset]]
    sig_count <- nrow(de_data[de_data$P.Value < 0.05, ])
    cat(sprintf("- %s: %d total genes, %d significant (%.1f%%)\n", 
                dataset, nrow(de_data), sig_count, 
                sig_count/nrow(de_data)*100), file = summary_file, append = TRUE)
  }
}

if(length(complement_results) > 0) {
  cat("\nCOMPLEMENT SYSTEM FINDINGS:\n", file = summary_file, append = TRUE)
  for(dataset in names(complement_results)) {
    comp_data <- complement_results[[dataset]]
    sig_comp <- nrow(comp_data[comp_data$P.Value < 0.05, ])
    cat(sprintf("- %s: %d complement genes, %d significantly changed\n", 
                dataset, nrow(comp_data), sig_comp), file = summary_file, append = TRUE)
  }
}

cat("\nNEXT STEPS:\n", file = summary_file, append = TRUE)
if(length(missing_datasets) > 0) {
  cat("1. Process missing datasets:", paste(missing_datasets, collapse = ", "), "\n", 
      file = summary_file, append = TRUE)
  cat("2. Re-run cross-dataset analysis for comprehensive comparison\n", 
      file = summary_file, append = TRUE)
} else {
  cat("1. All datasets processed - ready for meta-analysis\n", 
      file = summary_file, append = TRUE)
}

cat("3. Generate publication-ready figures\n", file = summary_file, append = TRUE)
cat("4. Perform functional enrichment meta-analysis\n", file = summary_file, append = TRUE)

cat("========================================================================\n", file = summary_file, append = TRUE)

cat("✓ Enhanced cross-dataset analysis complete!\n")
cat("✓ Summary report saved to:", summary_file, "\n")
cat("✓ Results saved to: Cross_Dataset_Analysis/\n")

if(length(missing_datasets) > 0) {
  cat("\n📋 TODO: Process missing datasets:\n")
  for(missing in missing_datasets) {
    cat("   -", missing, "\n")
  }
  cat("\nOnce all datasets are processed, re-run this script for full meta-analysis!\n")
} else {
  cat("\n🎉 All datasets ready for comprehensive meta-analysis!\n")
}
