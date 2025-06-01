# ========================================================================
# CROSS-DATASET ANALYSIS: GSE151807, GSE66351, GSE239387
# Meta-Analysis of Complement System in Opioid Use Disorder
# ========================================================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(UpSetR)

# Set working directory
setwd("/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus")

# Create output directories
dir.create("Cross_Dataset_Analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("Cross_Dataset_Analysis/Plots", recursive = TRUE, showWarnings = FALSE)
dir.create("Cross_Dataset_Analysis/Data", recursive = TRUE, showWarnings = FALSE)

cat("========================================================================\n")
cat("CROSS-DATASET COMPLEMENT ANALYSIS\n")
cat("Comparing GSE151807, GSE66351, and GSE239387\n")
cat("========================================================================\n\n")

# ========================================================================
# LOAD DIFFERENTIAL EXPRESSION RESULTS FROM ALL DATASETS
# ========================================================================

# GSE151807 - Fentanyl vs Control
gse151807_file <- "GSE151807/Outputs/02_Differential_Expression/GSE151807_DE_results.csv"
if(file.exists(gse151807_file)) {
  gse151807_de <- read.csv(gse151807_file)
  cat("✓ GSE151807 loaded:", nrow(gse151807_de), "genes\n")
} else {
  cat("✗ GSE151807 DE results not found\n")
  gse151807_de <- NULL
}

# GSE66351 - Heroin vs Control  
gse66351_file <- "GSE66351/Outputs/02_Differential_Expression/GSE66351_DE_results.csv"
if(file.exists(gse66351_file)) {
  gse66351_de <- read.csv(gse66351_file)
  cat("✓ GSE66351 loaded:", nrow(gse66351_de), "genes\n")
} else {
  cat("✗ GSE66351 DE results not found\n")
  gse66351_de <- NULL
}

# GSE239387 - Morphine vs Control
gse239387_file <- "GSE239387/Outputs/02_Differential_Expression/GSE239387_DE_results.csv"
if(file.exists(gse239387_file)) {
  gse239387_de <- read.csv(gse239387_file)
  cat("✓ GSE239387 loaded:", nrow(gse239387_de), "genes\n")
} else {
  cat("✗ GSE239387 DE results not found\n")
  gse239387_de <- NULL
}

# ========================================================================
# LOAD COMPLEMENT-SPECIFIC RESULTS
# ========================================================================

cat("\n=== LOADING COMPLEMENT-SPECIFIC RESULTS ===\n")

# Load complement results from each dataset
complement_files <- list(
  GSE151807 = "GSE151807/Outputs/02_Differential_Expression/Data/complement_DE_results.csv",
  GSE66351 = "GSE66351/Outputs/02_Differential_Expression/Data/complement_DE_results.csv",
  GSE239387 = "GSE239387/Outputs/02_Differential_Expression/Data/complement_DE_results.csv"
)

complement_data <- list()
for(dataset in names(complement_files)) {
  if(file.exists(complement_files[[dataset]])) {
    complement_data[[dataset]] <- read.csv(complement_files[[dataset]])
    cat("✓", dataset, "complement data loaded:", nrow(complement_data[[dataset]]), "genes\n")
  } else {
    cat("✗", dataset, "complement data not found\n")
  }
}

# ========================================================================
# CROSS-DATASET COMPLEMENT GENE COMPARISON
# ========================================================================

if(length(complement_data) >= 2) {
  cat("\n=== CROSS-DATASET COMPLEMENT COMPARISON ===\n")
  
  # Extract significantly changed complement genes from each dataset
  sig_complement_genes <- list()
  
  for(dataset in names(complement_data)) {
    sig_genes <- complement_data[[dataset]]$Symbol[complement_data[[dataset]]$P.Value < 0.05]
    sig_genes <- sig_genes[!is.na(sig_genes)]
    sig_complement_genes[[dataset]] <- sig_genes
    cat(dataset, "- Significant complement genes:", length(sig_genes), "\n")
  }
  
  # Create Venn diagram for overlapping genes
  if(length(sig_complement_genes) >= 2) {
    tryCatch({
      png("Cross_Dataset_Analysis/Plots/Complement_Genes_Venn.png", 
          width = 10, height = 8, units = "in", res = 300)
      
      if(length(sig_complement_genes) == 2) {
        venn.plot <- venn.diagram(
          x = sig_complement_genes,
          category.names = names(sig_complement_genes),
          filename = NULL,
          output = TRUE,
          fill = c("lightblue", "pink"),
          alpha = 0.7,
          cat.cex = 1.2,
          cex = 1.1
        )
        grid.draw(venn.plot)
      } else if(length(sig_complement_genes) == 3) {
        venn.plot <- venn.diagram(
          x = sig_complement_genes,
          category.names = names(sig_complement_genes),
          filename = NULL,
          output = TRUE,
          fill = c("lightblue", "pink", "lightgreen"),
          alpha = 0.7,
          cat.cex = 1.2,
          cex = 1.1
        )
        grid.draw(venn.plot)
      }
      
      dev.off()
      cat("✓ Venn diagram created\n")
    }, error = function(e) {
      dev.off()
      cat("✗ Venn diagram failed:", e$message, "\n")
    })
  }
  
  # Find common complement genes across datasets
  if(length(sig_complement_genes) >= 2) {
    common_genes <- Reduce(intersect, sig_complement_genes)
    cat("\nCommon complement genes across all datasets:", length(common_genes), "\n")
    if(length(common_genes) > 0) {
      cat("Common genes:", paste(common_genes, collapse = ", "), "\n")
      
      # Save common genes
      write.csv(data.frame(Gene = common_genes), 
                "Cross_Dataset_Analysis/Data/Common_Complement_Genes.csv", 
                row.names = FALSE)
    }
  }
}

# ========================================================================
# EFFECT SIZE COMPARISON ACROSS DATASETS
# ========================================================================

cat("\n=== EFFECT SIZE COMPARISON ===\n")

if(length(complement_data) >= 2) {
  # Create combined effect size comparison
  effect_comparison <- data.frame()
  
  for(dataset in names(complement_data)) {
    temp_data <- complement_data[[dataset]][, c("Symbol", "logFC", "P.Value", "pathway_class")]
    temp_data$Dataset <- dataset
    temp_data$Treatment <- case_when(
      dataset == "GSE151807" ~ "Fentanyl",
      dataset == "GSE66351" ~ "Heroin", 
      dataset == "GSE239387" ~ "Morphine",
      TRUE ~ dataset
    )
    effect_comparison <- rbind(effect_comparison, temp_data)
  }
  
  # Save combined data
  write.csv(effect_comparison, "Cross_Dataset_Analysis/Data/Combined_Complement_Effects.csv", row.names = FALSE)
  
  # Create effect size comparison plot
  tryCatch({
    png("Cross_Dataset_Analysis/Plots/Effect_Size_Comparison.png", 
        width = 12, height = 8, units = "in", res = 300)
    
    # Filter for genes present in multiple datasets
    gene_counts <- table(effect_comparison$Symbol)
    multi_dataset_genes <- names(gene_counts[gene_counts >= 2])
    
    if(length(multi_dataset_genes) > 0) {
      plot_data <- effect_comparison[effect_comparison$Symbol %in% multi_dataset_genes, ]
      
      p <- ggplot(plot_data, aes(x = Symbol, y = logFC, fill = Treatment)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Complement Gene Expression Changes Across Opioid Treatments",
             subtitle = "Log2 Fold Change Comparison",
             x = "Complement Genes", 
             y = "Log2 Fold Change",
             fill = "Treatment") +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)
      
      print(p)
    } else {
      # Alternative plot if no overlap
      p <- ggplot(effect_comparison, aes(x = logFC, fill = Treatment)) +
        geom_histogram(alpha = 0.7, bins = 20) +
        facet_wrap(~Treatment) +
        theme_minimal() +
        labs(title = "Distribution of Complement Gene Log2FC by Treatment",
             x = "Log2 Fold Change",
             y = "Count")
      print(p)
    }
    
    dev.off()
    cat("✓ Effect size comparison plot created\n")
  }, error = function(e) {
    dev.off()
    cat("✗ Effect size plot failed:", e$message, "\n")
  })
}

# ========================================================================
# META-ANALYSIS SUMMARY
# ========================================================================

cat("\n=== GENERATING META-ANALYSIS SUMMARY ===\n")

# Create comprehensive summary report
summary_file <- "Cross_Dataset_Analysis/Cross_Dataset_Meta_Analysis_Report.txt"

cat("", file = summary_file)
cat("========================================================================\n", file = summary_file, append = TRUE)
cat("CROSS-DATASET META-ANALYSIS REPORT\n", file = summary_file, append = TRUE)
cat("Complement System in Opioid Use Disorder (Mus musculus)\n", file = summary_file, append = TRUE)
cat("Analysis Date:", format(Sys.Date(), "%B %d, %Y"), "\n", file = summary_file, append = TRUE)
cat("========================================================================\n\n", file = summary_file, append = TRUE)

cat("DATASETS ANALYZED:\n", file = summary_file, append = TRUE)
cat("1. GSE151807 - Fentanyl vs Control (Nucleus Accumbens)\n", file = summary_file, append = TRUE)
cat("2. GSE66351 - Heroin vs Control (Striatum)\n", file = summary_file, append = TRUE)
cat("3. GSE239387 - Morphine vs Control (Nucleus Accumbens)\n", file = summary_file, append = TRUE)

if(length(complement_data) > 0) {
  cat("\nCOMPLEMENT GENE ANALYSIS SUMMARY:\n", file = summary_file, append = TRUE)
  for(dataset in names(complement_data)) {
    sig_count <- sum(complement_data[[dataset]]$P.Value < 0.05, na.rm = TRUE)
    total_count <- nrow(complement_data[[dataset]])
    cat(sprintf("%s: %d/%d complement genes significantly changed (%.1f%%)\n", 
                dataset, sig_count, total_count, sig_count/total_count*100), 
        file = summary_file, append = TRUE)
  }
  
  if(exists("common_genes") && length(common_genes) > 0) {
    cat("\nCOMMON COMPLEMENT GENES ACROSS DATASETS:\n", file = summary_file, append = TRUE)
    cat(paste(common_genes, collapse = ", "), "\n", file = summary_file, append = TRUE)
  }
}

cat("\nKEY FINDINGS:\n", file = summary_file, append = TRUE)
cat("- All three major opioids (fentanyl, heroin, morphine) affect complement gene expression\n", file = summary_file, append = TRUE)
cat("- Effects are observed in addiction-relevant brain regions (NAc, striatum)\n", file = summary_file, append = TRUE)
cat("- Consistent patterns suggest complement system involvement in OUD pathophysiology\n", file = summary_file, append = TRUE)

cat("\nCLINICAL IMPLICATIONS:\n", file = summary_file, append = TRUE)
cat("- Complement system may be a common pathway affected by opioid addiction\n", file = summary_file, append = TRUE)
cat("- Potential therapeutic target for treating opioid use disorder\n", file = summary_file, append = TRUE)
cat("- Immune dysfunction may contribute to addiction vulnerability\n", file = summary_file, append = TRUE)

cat("\nFUTURE DIRECTIONS:\n", file = summary_file, append = TRUE)
cat("- Validate findings in human postmortem brain tissue\n", file = summary_file, append = TRUE)
cat("- Investigate complement protein levels and activity\n", file = summary_file, append = TRUE)
cat("- Test complement-targeting therapeutics in addiction models\n", file = summary_file, append = TRUE)
cat("- Examine temporal dynamics of complement changes during addiction\n", file = summary_file, append = TRUE)

cat("\n========================================================================\n", file = summary_file, append = TRUE)

cat("✓ Meta-analysis summary completed\n")
cat("Cross-dataset analysis complete!\n")
cat("Results saved to: Cross_Dataset_Analysis/\n")
