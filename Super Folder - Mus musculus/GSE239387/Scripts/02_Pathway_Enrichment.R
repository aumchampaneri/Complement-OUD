# ========================================================================
# GSE239387 Pathway Enrichment Analysis
# KEGG and GO Pathway Analysis of Morphine vs Control
# ========================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

# Set working directory
setwd("/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE239387")

# Create output directory
dir.create("Outputs/03_Pathway_Enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("Outputs/03_Pathway_Enrichment/Plots", recursive = TRUE, showWarnings = FALSE)

# Load differential expression results
de_results <- read.csv("Outputs/02_Differential_Expression/GSE239387_DE_results.csv")

# Prepare gene lists
# Significantly upregulated genes (p < 0.05)
up_genes <- de_results$Symbol[de_results$P.Value < 0.05 & de_results$logFC > 0]
up_genes <- up_genes[!is.na(up_genes)]

# Significantly downregulated genes (p < 0.05)
down_genes <- de_results$Symbol[de_results$P.Value < 0.05 & de_results$logFC < 0]
down_genes <- down_genes[!is.na(down_genes)]

# All significant genes
all_sig_genes <- de_results$Symbol[de_results$P.Value < 0.05]
all_sig_genes <- all_sig_genes[!is.na(all_sig_genes)]

cat("Pathway enrichment analysis:\n")
cat("Upregulated genes:", length(up_genes), "\n")
cat("Downregulated genes:", length(down_genes), "\n")
cat("Total significant genes:", length(all_sig_genes), "\n")

# KEGG pathway enrichment
if(length(all_sig_genes) > 10) {
  tryCatch({
    cat("Converting gene symbols to Entrez IDs...\n")
    
    # Convert gene symbols to Entrez IDs for KEGG analysis with error handling
    entrez_conversion <- bitr(all_sig_genes, 
                             fromType = "SYMBOL", 
                             toType = "ENTREZID", 
                             OrgDb = org.Mm.eg.db,
                             drop = TRUE)
    
    cat("Converted", nrow(entrez_conversion), "out of", length(all_sig_genes), "genes to Entrez IDs\n")
    
    if(nrow(entrez_conversion) > 5) {
      kegg_enrich <- enrichKEGG(gene = entrez_conversion$ENTREZID,
                                organism = 'mmu',
                                pAdjustMethod = 'BH',
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2,
                                minGSSize = 3,
                                maxGSSize = 500)
      
      if(!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
        write.csv(kegg_enrich@result, "Outputs/03_Pathway_Enrichment/KEGG_enrichment_results.csv", row.names = FALSE)
        
        # Plot KEGG results
        png("Outputs/03_Pathway_Enrichment/Plots/KEGG_dotplot.png", width = 12, height = 8, units = "in", res = 300)
        print(dotplot(kegg_enrich, showCategory = 20) + ggtitle("KEGG Pathway Enrichment - Morphine vs Control"))
        dev.off()
        
        cat("✓ KEGG pathway analysis completed\n")
        cat("Significant KEGG pathways found:", nrow(kegg_enrich@result), "\n")
        
        # Print top pathways
        if(nrow(kegg_enrich@result) > 0) {
          cat("\nTop KEGG pathways:\n")
          top_kegg <- head(kegg_enrich@result[order(kegg_enrich@result$pvalue), ], 5)
          for(i in 1:nrow(top_kegg)) {
            cat(sprintf("  %d. %s (p = %.2e)\n", i, top_kegg$Description[i], top_kegg$pvalue[i]))
          }
        }
      } else {
        cat("No significant KEGG pathways found\n")
      }
    } else {
      cat("Not enough genes successfully converted to Entrez IDs for KEGG analysis\n")
    }
  }, error = function(e) {
    cat("KEGG analysis failed:", e$message, "\n")
    cat("This might be due to internet connectivity issues or gene ID conversion problems\n")
  })
} else {
  cat("Not enough significant genes for KEGG analysis\n")
}

# GO enrichment analysis
tryCatch({
  cat("Running GO Biological Process enrichment...\n")
  go_bp <- enrichGO(gene = all_sig_genes,
                    OrgDb = org.Mm.eg.db,
                    keyType = 'SYMBOL',
                    ont = 'BP',
                    pAdjustMethod = 'BH',
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    minGSSize = 3,
                    maxGSSize = 500)
  
  if(!is.null(go_bp) && nrow(go_bp@result) > 0) {
    write.csv(go_bp@result, "Outputs/03_Pathway_Enrichment/GO_BP_enrichment_results.csv", row.names = FALSE)
    
    # Plot GO results
    png("Outputs/03_Pathway_Enrichment/Plots/GO_BP_dotplot.png", width = 12, height = 10, units = "in", res = 300)
    print(dotplot(go_bp, showCategory = 20) + ggtitle("GO Biological Process Enrichment - Morphine vs Control"))
    dev.off()
    
    cat("✓ GO Biological Process analysis completed\n")
    cat("Significant GO BP terms found:", nrow(go_bp@result), "\n")
    
    # Print top GO terms
    if(nrow(go_bp@result) > 0) {
      cat("\nTop GO Biological Process terms:\n")
      top_go <- head(go_bp@result[order(go_bp@result$pvalue), ], 5)
      for(i in 1:nrow(top_go)) {
        cat(sprintf("  %d. %s (p = %.2e)\n", i, top_go$Description[i], top_go$pvalue[i]))
      }
    }
  } else {
    cat("No significant GO BP terms found\n")
  }
}, error = function(e) cat("GO analysis failed:", e$message, "\n"))

# Additional GO analyses
# GO Molecular Function
tryCatch({
  go_mf <- enrichGO(gene = all_sig_genes,
                    OrgDb = org.Mm.eg.db,
                    keyType = 'SYMBOL',
                    ont = 'MF',
                    pAdjustMethod = 'BH',
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)
  
  if(!is.null(go_mf) && nrow(go_mf@result) > 0) {
    write.csv(go_mf@result, "Outputs/03_Pathway_Enrichment/GO_MF_enrichment_results.csv", row.names = FALSE)
    cat("✓ GO Molecular Function analysis completed\n")
  }
}, error = function(e) cat("GO MF analysis failed:", e$message, "\n"))

# GO Cellular Component
tryCatch({
  go_cc <- enrichGO(gene = all_sig_genes,
                    OrgDb = org.Mm.eg.db,
                    keyType = 'SYMBOL',
                    ont = 'CC',
                    pAdjustMethod = 'BH',
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)
  
  if(!is.null(go_cc) && nrow(go_cc@result) > 0) {
    write.csv(go_cc@result, "Outputs/03_Pathway_Enrichment/GO_CC_enrichment_results.csv", row.names = FALSE)
    cat("✓ GO Cellular Component analysis completed\n")
  }
}, error = function(e) cat("GO CC analysis failed:", e$message, "\n"))

# ========================================================================
# COMPLEMENT-SPECIFIC PATHWAY ANALYSIS
# ========================================================================

cat("\n=== COMPLEMENT-SPECIFIC PATHWAY ANALYSIS ===\n")

# Load complement genes if available
if(file.exists("Outputs/02_Differential_Expression/Data/complement_DE_results.csv")) {
  complement_de <- read.csv("Outputs/02_Differential_Expression/Data/complement_DE_results.csv")
  sig_complement <- complement_de[complement_de$P.Value < 0.05, ]
  
  if(nrow(sig_complement) > 3) {
    cat("Analyzing", nrow(sig_complement), "significantly changed complement genes\n")
    
    # Create complement-specific gene list
    complement_genes <- sig_complement$Symbol[!is.na(sig_complement$Symbol)]
    
    # GO analysis for complement genes
    tryCatch({
      complement_go <- enrichGO(gene = complement_genes,
                                OrgDb = org.Mm.eg.db,
                                keyType = 'SYMBOL',
                                ont = 'BP',
                                pAdjustMethod = 'BH',
                                pvalueCutoff = 0.1,  # More lenient for smaller gene set
                                qvalueCutoff = 0.3)
      
      if(!is.null(complement_go) && nrow(complement_go@result) > 0) {
        write.csv(complement_go@result, "Outputs/03_Pathway_Enrichment/Complement_GO_enrichment.csv", row.names = FALSE)
        
        # Plot complement GO results
        png("Outputs/03_Pathway_Enrichment/Plots/Complement_GO_dotplot.png", width = 10, height = 6, units = "in", res = 300)
        print(dotplot(complement_go, showCategory = 15) + 
              ggtitle("GO Enrichment: Complement Genes (Morphine vs Control)"))
        dev.off()
        
        cat("✓ Complement GO analysis completed\n")
        cat("Significant complement GO terms found:", nrow(complement_go@result), "\n")
      } else {
        cat("No significant GO terms found for complement genes\n")
      }
    }, error = function(e) cat("Complement GO analysis failed:", e$message, "\n"))
    
    # Additional complement analysis - compare with immune pathways
    cat("\nAnalyzing complement gene overlap with immune pathways...\n")
    
    # Save detailed complement results (fix dplyr conflict)
    if(nrow(sig_complement) > 0) {
      complement_summary <- sig_complement %>%
        dplyr::arrange(P.Value) %>%
        dplyr::select(Symbol, logFC, P.Value, pathway_class) %>%
        dplyr::mutate(
          Direction = ifelse(logFC > 0, "Upregulated", "Downregulated"),
          Significance = dplyr::case_when(
            P.Value < 0.001 ~ "***",
            P.Value < 0.01 ~ "**", 
            P.Value < 0.05 ~ "*",
            TRUE ~ ""
          )
        )
      
      write.csv(complement_summary, "Outputs/03_Pathway_Enrichment/Complement_Gene_Summary.csv", row.names = FALSE)
      cat("✓ Complement gene summary saved\n")
      
      # Print complement gene summary
      cat("\nSummary of significantly changed complement genes:\n")
      up_complement <- sum(complement_summary$Direction == "Upregulated")
      down_complement <- sum(complement_summary$Direction == "Downregulated")
      cat("  - Upregulated complement genes:", up_complement, "\n")
      cat("  - Downregulated complement genes:", down_complement, "\n")
      
      # Show most significant complement genes
      cat("\nMost significant complement genes:\n")
      top_complement <- head(complement_summary, 5)
      for(i in 1:nrow(top_complement)) {
        cat(sprintf("  %d. %s (%s, logFC = %.2f, p = %.2e)\n", 
                    i, top_complement$Symbol[i], top_complement$Direction[i],
                    top_complement$logFC[i], top_complement$P.Value[i]))
      }
    }
  } else {
    cat("Not enough significantly changed complement genes for pathway analysis\n")
  }
} else {
  cat("Complement DE results not found - run differential expression analysis first\n")
}

# ========================================================================
# DIRECTIONAL PATHWAY ANALYSIS
# ========================================================================

cat("\n=== DIRECTIONAL PATHWAY ANALYSIS ===\n")

# Separate analysis for up and down regulated genes
if(length(up_genes) > 5) {
  cat("Analyzing upregulated genes (n =", length(up_genes), ")\n")
  
  tryCatch({
    up_go <- enrichGO(gene = up_genes,
                      OrgDb = org.Mm.eg.db,
                      keyType = 'SYMBOL',
                      ont = 'BP',
                      pAdjustMethod = 'BH',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
    
    if(!is.null(up_go) && nrow(up_go@result) > 0) {
      write.csv(up_go@result, "Outputs/03_Pathway_Enrichment/Upregulated_GO_enrichment.csv", row.names = FALSE)
      
      png("Outputs/03_Pathway_Enrichment/Plots/Upregulated_GO_dotplot.png", width = 12, height = 8, units = "in", res = 300)
      print(dotplot(up_go, showCategory = 15) + 
            ggtitle("GO Enrichment: Upregulated Genes (Morphine vs Control)"))
      dev.off()
      
      cat("✓ Upregulated genes GO analysis completed\n")
      cat("Significant upregulated GO terms found:", nrow(up_go@result), "\n")
    } else {
      cat("No significant GO terms found for upregulated genes\n")
    }
  }, error = function(e) cat("Upregulated GO analysis failed:", e$message, "\n"))
} else {
  cat("Not enough upregulated genes for pathway analysis\n")
}

if(length(down_genes) > 5) {
  cat("Analyzing downregulated genes (n =", length(down_genes), ")\n")
  
  tryCatch({
    down_go <- enrichGO(gene = down_genes,
                        OrgDb = org.Mm.eg.db,
                        keyType = 'SYMBOL',
                        ont = 'BP',
                        pAdjustMethod = 'BH',
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
    
    if(!is.null(down_go) && nrow(down_go@result) > 0) {
      write.csv(down_go@result, "Outputs/03_Pathway_Enrichment/Downregulated_GO_enrichment.csv", row.names = FALSE)
      
      png("Outputs/03_Pathway_Enrichment/Plots/Downregulated_GO_dotplot.png", width = 12, height = 8, units = "in", res = 300)
      print(dotplot(down_go, showCategory = 15) + 
            ggtitle("GO Enrichment: Downregulated Genes (Morphine vs Control)"))
      dev.off()
      
      cat("✓ Downregulated genes GO analysis completed\n")
      cat("Significant downregulated GO terms found:", nrow(down_go@result), "\n")
    } else {
      cat("No significant GO terms found for downregulated genes\n")
    }
  }, error = function(e) cat("Downregulated GO analysis failed:", e$message, "\n"))
} else {
  cat("Not enough downregulated genes for pathway analysis\n")
}

# ========================================================================
# FUNCTIONAL THEME ANALYSIS
# ========================================================================

cat("\n=== FUNCTIONAL THEME ANALYSIS ===\n")

# Look for specific functional themes related to addiction, neuroplasticity, immune function
addiction_related_terms <- c("addiction", "reward", "dopamine", "neurotransmitter", "synaptic", 
                            "plasticity", "learning", "memory", "behavior")

immune_related_terms <- c("immune", "inflammation", "cytokine", "complement", "innate immunity",
                         "adaptive immunity", "defense response")

if(exists("go_bp") && !is.null(go_bp) && nrow(go_bp@result) > 0) {
  # Search for addiction-related terms
  addiction_terms <- go_bp@result[grepl(paste(addiction_related_terms, collapse = "|"), 
                                       go_bp@result$Description, ignore.case = TRUE), ]
  
  if(nrow(addiction_terms) > 0) {
    cat("Found", nrow(addiction_terms), "addiction/neuroplasticity-related GO terms\n")
    write.csv(addiction_terms, "Outputs/03_Pathway_Enrichment/Addiction_Related_GO_terms.csv", row.names = FALSE)
  }
  
  # Search for immune-related terms
  immune_terms <- go_bp@result[grepl(paste(immune_related_terms, collapse = "|"), 
                                    go_bp@result$Description, ignore.case = TRUE), ]
  
  if(nrow(immune_terms) > 0) {
    cat("Found", nrow(immune_terms), "immune-related GO terms\n")
    write.csv(immune_terms, "Outputs/03_Pathway_Enrichment/Immune_Related_GO_terms.csv", row.names = FALSE)
    
    # Show top immune-related terms
    cat("\nTop immune-related GO terms:\n")
    top_immune <- head(immune_terms[order(immune_terms$pvalue), ], 5)
    for(i in 1:nrow(top_immune)) {
      cat(sprintf("  %d. %s (p = %.2e)\n", i, top_immune$Description[i], top_immune$pvalue[i]))
    }
  }
}

# Specific analysis for morphine addiction pathway from KEGG
if(exists("kegg_enrich") && !is.null(kegg_enrich)) {
  morphine_pathway <- kegg_enrich@result[grepl("morphine|addiction|cocaine|amphetamine", 
                                               kegg_enrich@result$Description, ignore.case = TRUE), ]
  
  if(nrow(morphine_pathway) > 0) {
    cat("\nAddiction-related KEGG pathways found:\n")
    for(i in 1:nrow(morphine_pathway)) {
      cat(sprintf("  %d. %s (p = %.2e, genes = %d)\n", 
                  i, morphine_pathway$Description[i], 
                  morphine_pathway$pvalue[i], morphine_pathway$Count[i]))
    }
    write.csv(morphine_pathway, "Outputs/03_Pathway_Enrichment/Addiction_KEGG_pathways.csv", row.names = FALSE)
  }
}

# ========================================================================
# PATHWAY ANALYSIS SUMMARY AND REPORT
# ========================================================================

# Create comprehensive summary report
summary_file <- "Outputs/03_Pathway_Enrichment/Pathway_Analysis_Summary.txt"

cat("", file = summary_file)
cat("========================================================================\n", file = summary_file, append = TRUE)
cat("GSE239387 PATHWAY ENRICHMENT ANALYSIS SUMMARY\n", file = summary_file, append = TRUE)
cat("Morphine vs Control in Nucleus Accumbens\n", file = summary_file, append = TRUE)
cat("Analysis Date:", format(Sys.Date(), "%B %d, %Y"), "\n", file = summary_file, append = TRUE)
cat("========================================================================\n\n", file = summary_file, append = TRUE)

cat("GENE SET SIZES:\n", file = summary_file, append = TRUE)
cat("Total significant genes (p < 0.05):", length(all_sig_genes), "\n", file = summary_file, append = TRUE)
cat("Upregulated genes:", length(up_genes), "\n", file = summary_file, append = TRUE)
cat("Downregulated genes:", length(down_genes), "\n", file = summary_file, append = TRUE)

if(exists("sig_complement") && nrow(sig_complement) > 0) {
  cat("Significant complement genes:", nrow(sig_complement), "\n", file = summary_file, append = TRUE)
}

cat("\nANALYSES PERFORMED:\n", file = summary_file, append = TRUE)
cat("- KEGG pathway enrichment\n", file = summary_file, append = TRUE)
cat("- GO Biological Process enrichment\n", file = summary_file, append = TRUE)
cat("- GO Molecular Function enrichment\n", file = summary_file, append = TRUE)
cat("- GO Cellular Component enrichment\n", file = summary_file, append = TRUE)
cat("- Directional pathway analysis (up/down genes)\n", file = summary_file, append = TRUE)
cat("- Complement-specific pathway analysis\n", file = summary_file, append = TRUE)
cat("- Functional theme analysis (addiction, immune)\n", file = summary_file, append = TRUE)

cat("\nKEY FINDINGS:\n", file = summary_file, append = TRUE)
cat("- Morphine treatment affects", length(all_sig_genes), "genes significantly\n", file = summary_file, append = TRUE)
cat("- More genes are", ifelse(length(down_genes) > length(up_genes), "DOWNREGULATED", "UPREGULATED"), 
    "than", ifelse(length(down_genes) > length(up_genes), "upregulated", "downregulated"), "\n", file = summary_file, append = TRUE)

if(exists("sig_complement") && nrow(sig_complement) > 0) {
  cat("- Complement system shows significant changes with morphine treatment\n", file = summary_file, append = TRUE)
}

# Add pathway enrichment results summary
if(exists("kegg_enrich") && !is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
  cat("- KEGG pathways significantly enriched:", nrow(kegg_enrich@result), "\n", file = summary_file, append = TRUE)
}

if(exists("go_bp") && !is.null(go_bp) && nrow(go_bp@result) > 0) {
  cat("- GO Biological Process terms enriched:", nrow(go_bp@result), "\n", file = summary_file, append = TRUE)
}

cat("\nOUTPUT FILES GENERATED:\n", file = summary_file, append = TRUE)
output_files <- c(
  "KEGG_enrichment_results.csv",
  "GO_BP_enrichment_results.csv", 
  "GO_MF_enrichment_results.csv",
  "GO_CC_enrichment_results.csv",
  "Upregulated_GO_enrichment.csv",
  "Downregulated_GO_enrichment.csv"
)

if(exists("sig_complement") && nrow(sig_complement) > 0) {
  output_files <- c(output_files, c("Complement_GO_enrichment.csv", "Complement_Gene_Summary.csv"))
}

for(file in output_files) {
  if(file.exists(paste0("Outputs/03_Pathway_Enrichment/", file))) {
    cat("✓", file, "\n", file = summary_file, append = TRUE)
  }
}

cat("\nPLOTS GENERATED:\n", file = summary_file, append = TRUE)
plot_files <- c(
  "KEGG_dotplot.png",
  "GO_BP_dotplot.png",
  "Upregulated_GO_dotplot.png", 
  "Downregulated_GO_dotplot.png",
  "Complement_GO_dotplot.png"
)

for(file in plot_files) {
  if(file.exists(paste0("Outputs/03_Pathway_Enrichment/Plots/", file))) {
    cat("✓", file, "\n", file = summary_file, append = TRUE)
  }
}

cat("\n========================================================================\n", file = summary_file, append = TRUE)
cat("END OF PATHWAY ANALYSIS SUMMARY\n", file = summary_file, append = TRUE)
cat("========================================================================\n", file = summary_file, append = TRUE)

cat("\n✓ Pathway enrichment analysis complete!\n")
cat("Summary saved to:", summary_file, "\n")
cat("All results saved to: Outputs/03_Pathway_Enrichment/\n")
cat("\nReady for cross-dataset analysis once other datasets are processed!\n")
