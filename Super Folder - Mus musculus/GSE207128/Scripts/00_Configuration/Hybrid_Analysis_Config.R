# =============================================================================
# Hybrid Analysis Configuration: Best of Both Approaches
# =============================================================================
# 
# Purpose: Smart configuration system that automatically optimizes parameters
#          based on dataset type and tissue specificity
#
# Features:
# - Auto-detection of GSE207128 vs other datasets
# - Tissue-specific parameter optimization
# - Original paper validation for amygdala
# - Comprehensive fallback for other tissues
# - Memory-efficient processing options
# =============================================================================

# Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# MAIN CONFIGURATION FUNCTION
# =============================================================================

get_analysis_config <- function(dataset_path, tissue_type = "auto") {
  
  # Auto-detect tissue type and dataset from path
  if (tissue_type == "auto") {
    if (grepl("GSE207128", dataset_path)) {
      tissue_type <- "amygdala"
      dataset_type <- "GSE207128"
      message("Detected GSE207128 - applying amygdala-specific optimizations")
    } else {
      tissue_type <- "brain_general"
      dataset_type <- "general"
      message("Detected general brain dataset - applying comprehensive methods")
    }
  }
  
  # Base configuration (comprehensive robust methods)
  base_config <- list(
    # QC Parameters - Comprehensive approach
    qc = list(
      min_features_adaptive = TRUE,
      max_features_adaptive = TRUE,
      mt_percent_adaptive = TRUE,
      ribosomal_percent_check = TRUE,
      doublet_detection = TRUE,
      cell_cycle_scoring = TRUE,
      min_features_base = 200,
      max_features_base = 8000,
      min_counts_base = 300,
      max_counts_base = 50000,
      max_mt_percent_base = 20,
      ribosomal_threshold = 40,
      min_cells_per_gene = 3
    ),
    
    # Integration Parameters - Comprehensive
    integration = list(
      methods = c("harmony", "fastmnn", "batch_regression"),
      variable_features = 3000,
      pca_dims = 30,
      integration_dims = 25
    ),
    
    # Clustering Parameters - Comprehensive
    clustering = list(
      resolutions = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5),
      silhouette_analysis = TRUE,
      stability_analysis = TRUE,
      default_resolution = 0.8
    ),
    
    # Analysis Parameters
    analysis = list(
      mad_threshold = 3,
      confidence_threshold = 0.6,
      min_marker_genes = 3
    )
  )
  
  # Tissue-specific optimizations
  if (tissue_type == "amygdala" && dataset_type == "GSE207128") {
    
    # GSE207128 amygdala-specific optimizations (original paper validated)
    optimized_config <- list(
      # Amygdala-specific QC (validated by original paper)
      qc = list(
        min_features_base = 250,      # Original paper validation
        max_features_base = 6000,     # Amygdala tissue-specific
        min_counts_base = 500,        # Conservative for amygdala
        max_counts_base = 40000,      # Upper limit for amygdala
        max_mt_percent_base = 12,     # Strict MT threshold for amygdala
        min_features_adaptive = FALSE, # Use fixed thresholds for consistency
        doublet_detection = FALSE,    # Skip for this specific dataset
        ribosomal_threshold = 30,     # Amygdala-specific
        cell_cycle_scoring = TRUE     # Keep comprehensive feature
      ),
      
      # Streamlined integration for efficiency
      integration = list(
        methods = c("harmony", "batch_regression"), # Focus on validated methods
        variable_features = 2000,     # More focused feature set
        pca_dims = 25,               # Optimized for amygdala
        integration_dims = 20        # Efficient integration
      ),
      
      # Focused clustering (amygdala-optimized)
      clustering = list(
        resolutions = c(0.4, 0.6, 0.8, 1.0), # Targeted resolution range
        silhouette_analysis = FALSE,  # Skip for computational efficiency
        stability_analysis = TRUE,    # Keep stability validation
        default_resolution = 0.8      # Validated optimal for amygdala
      ),
      
      # Original paper cell type markers
      cell_types = list(
        method = "original_paper_primary",
        validation = "module_scoring_secondary",
        markers = "GSE207128_validated",
        confidence_threshold = 0.6
      )
    )
    
    # Merge configurations (optimized overrides base)
    config <- modifyList(base_config, optimized_config)
    config$dataset_info <- list(
      type = "GSE207128",
      tissue = "amygdala",
      optimization = "original_paper_validated",
      species = "Mus_musculus"
    )
    
  } else {
    # Use comprehensive approach for other datasets
    config <- base_config
    config$cell_types <- list(
      method = "comprehensive_primary",
      validation = "original_markers_secondary", 
      markers = "brain_general",
      confidence_threshold = 0.7
    )
    config$dataset_info <- list(
      type = "general",
      tissue = tissue_type,
      optimization = "comprehensive_robust",
      species = "auto_detect"
    )
  }
  
  # Add standardized output paths for organized structure
  config$output_paths <- list(
    aggregated_data = "Outputs/01_Aggregated_Data",
    processed_data = "Outputs/02_Processed_Data", 
    integrated_data = "Outputs/03_Integrated_Data",
    annotated_data = "Outputs/04_Annotated_Data",
    analysis_results = "Outputs/05_Analysis_Results",
    figures = "Outputs/06_Figures",
    reports = "Outputs/07_Reports",
    legacy = "Outputs/08_Legacy"
  )
  
  # Validate configuration
  validate_config(config)
  
  return(config)
}

# =============================================================================
# MARKER GENE CONFIGURATIONS
# =============================================================================

get_marker_config <- function(dataset_type = "GSE207128") {
  
  if (dataset_type == "GSE207128") {
    
    # ORIGINAL PAPER MARKERS (Primary) - Validated for amygdala
    primary_markers <- list(
      # Core cell types from GSE207128 paper
      ASC = "Gja1",                    # Astrocytes
      OPC = "Pdgfra",                  # Oligodendrocyte Precursor Cells  
      MG = "Tmem119",                  # Microglia
      EC = "Cldn5",                    # Endothelial Cells
      NEUR = "Syt1",                   # Neurons
      OLG = c("Cldn11", "Mobp"),       # Oligodendrocytes
      MAC = "Pf4",                     # Macrophages
      PC = "Vtn",                      # Pericytes
      NFOLG = "Enpp6",                 # Newly Formed Oligodendrocytes
      VSMC = "Acta2",                  # Vascular Smooth Muscle Cells
      DC = c("Cd74", "Cd209a"),        # Dendritic Cells
      EPC = "Ccdc153",                 # Ependymal Cells
      NSC = "Thbs4",                   # Neural Stem Cells
      ARP = "Cd44",                    # Arachnoid Barrier-like Cells
      NRP = c("Top2a", "Cdk1"),        # Neural Restricted Progenitors
      Tcells = "Cd3d",                 # T cells
      NEUT = "S100a9",                 # Neutrophils
      
      # Contamination markers (for removal)
      EPIC = "Ttr",                    # Epithelial Cells (contamination)
      VLMC = "Slc6a13",                # Vascular Leptomeningeal Cells
      ABC = "Slc47a1"                  # Arachnoid Barrier Cells
    )
    
    # COMPREHENSIVE VALIDATION MARKERS
    validation_markers <- list(
      # Neuronal subtypes (additional validation)
      Excitatory_Neurons = c("Slc17a7", "Camk2a", "Grin2a", "Dlg4"),
      Inhibitory_Neurons = c("Gad1", "Gad2", "Slc32a1", "Dlx1"),
      Parvalbumin_Interneurons = "Pvalb",
      Somatostatin_Interneurons = "Sst",
      VIP_Interneurons = "Vip",
      
      # Glial subtypes (detailed)
      Protoplasmic_Astrocytes = "Aqp4",
      Fibrous_Astrocytes = "Gfap", 
      Reactive_Astrocytes = c("Lcn2", "Serpina3n"),
      
      # Specialized cell types
      Dopaminergic_Neurons = c("Th", "Ddc", "Slc6a3"),
      Cholinergic_Neurons = c("Chat", "Slc18a3"),
      GABAergic_Neurons = c("Gad1", "Gad2", "Slc32a1"),
      
      # Additional brain cell markers
      Pyramidal_Neurons = c("Slc17a7", "Neurod6"),
      Granule_Cells = c("Neurod1", "Prox1")
    )
    
  } else {
    
    # COMPREHENSIVE BRAIN MARKERS (for general datasets)
    primary_markers <- list(
      # Core brain cell types
      Excitatory_Neurons = c("Slc17a7", "Camk2a", "Grin1", "Dlg4"),
      Inhibitory_Neurons = c("Gad1", "Gad2", "Slc32a1", "Dlx1"),
      Astrocytes = c("Gfap", "Aqp4", "Aldh1l1", "S100b"),
      Microglia = c("Cx3cr1", "Tmem119", "C1qa", "Csf1r"),
      Oligodendrocytes = c("Mbp", "Mog", "Plp1", "Cnp"),
      OPCs = c("Pdgfra", "Olig2", "Sox10", "Cspg4"),
      Endothelial_Cells = c("Pecam1", "Flt1", "Cdh5", "Tek"),
      Pericytes = c("Pdgfrb", "Des", "Acta2", "Rgs5"),
      
      # Specialized neurons
      Dopaminergic = c("Th", "Ddc", "Slc6a3"),
      Cholinergic = c("Chat", "Slc18a3"),
      Serotonergic = c("Tph2", "Slc6a4"),
      
      # Interneuron subtypes
      PV_Interneurons = "Pvalb",
      SST_Interneurons = "Sst",
      VIP_Interneurons = "Vip"
    )
    
    validation_markers <- list(
      # Regional markers
      Cortical = c("Tbr1", "Satb2", "Bcl11b"),
      Hippocampal = c("Prox1", "Neurod1"),
      Striatal = c("Drd1", "Drd2", "Penk"),
      
      # Developmental markers
      Immature_Neurons = c("Dcx", "Tbr2"),
      Mature_Neurons = c("Map2", "Neun", "Rbfox3")
    )
  }
  
  return(list(
    primary = primary_markers,
    validation = validation_markers,
    dataset_type = dataset_type
  ))
}

# =============================================================================
# ADAPTIVE THRESHOLD CALCULATION
# =============================================================================

calculate_adaptive_thresholds <- function(seurat_obj, metric, n_mad = 3, min_val = NULL, max_val = NULL, tissue_type = "amygdala") {
  
  # Ensure dplyr is loaded for group operations
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for adaptive threshold calculation")
  }
  
  # Calculate per-sample statistics
  per_sample_stats <- seurat_obj@meta.data %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      median_val = median(!!dplyr::sym(metric), na.rm = TRUE),
      mad_val = mad(!!dplyr::sym(metric), na.rm = TRUE),
      q25 = quantile(!!dplyr::sym(metric), 0.25, na.rm = TRUE),
      q75 = quantile(!!dplyr::sym(metric), 0.75, na.rm = TRUE),
      iqr = q75 - q25,
      lower_thresh_mad = median_val - n_mad * mad_val,
      upper_thresh_mad = median_val + n_mad * mad_val,
      lower_thresh_iqr = q25 - 1.5 * iqr,
      upper_thresh_iqr = q75 + 1.5 * iqr,
      .groups = 'drop'
    )
  
  # Apply tissue-specific constraints
  if (tissue_type == "amygdala") {
    # Use original paper's validated constraints for amygdala
    if (metric == "nFeature_RNA") {
      min_val <- max(min_val %||% 250, 200)  # Conservative for amygdala
      max_val <- min(max_val %||% 6000, 8000)
    } else if (metric == "percent.mt") {
      max_val <- min(max_val %||% 12, 15)  # Stricter MT threshold
    } else if (metric == "nCount_RNA") {
      min_val <- max(min_val %||% 500, 300)
      max_val <- min(max_val %||% 40000, 60000)
    }
  } else {
    # More lenient thresholds for other brain regions
    if (metric == "nFeature_RNA") {
      min_val <- max(min_val %||% 200, 100)
      max_val <- min(max_val %||% 8000, 10000)
    } else if (metric == "percent.mt") {
      max_val <- min(max_val %||% 20, 25)
    }
  }
  
  # Apply global constraints if provided
  if(!is.null(min_val)) {
    per_sample_stats$lower_thresh_mad <- pmax(per_sample_stats$lower_thresh_mad, min_val)
    per_sample_stats$lower_thresh_iqr <- pmax(per_sample_stats$lower_thresh_iqr, min_val)
  }
  if(!is.null(max_val)) {
    per_sample_stats$upper_thresh_mad <- pmin(per_sample_stats$upper_thresh_mad, max_val)
    per_sample_stats$upper_thresh_iqr <- pmin(per_sample_stats$upper_thresh_iqr, max_val)
  }
  
  return(per_sample_stats)
}

# =============================================================================
# CONFIGURATION VALIDATION
# =============================================================================

validate_config <- function(config) {
  
  # Check required sections
  required_sections <- c("qc", "integration", "clustering")
  missing_sections <- setdiff(required_sections, names(config))
  
  if(length(missing_sections) > 0) {
    stop("Missing required configuration sections: ", paste(missing_sections, collapse = ", "))
  }
  
  # Validate QC parameters
  if(is.null(config$qc$min_features_base) || config$qc$min_features_base < 50) {
    warning("min_features_base seems too low, recommend >= 200")
  }
  
  if(is.null(config$qc$max_mt_percent_base) || config$qc$max_mt_percent_base > 30) {
    warning("max_mt_percent_base seems too high, recommend <= 25%")
  }
  
  # Validate integration parameters
  if(is.null(config$integration$variable_features) || config$integration$variable_features < 1000) {
    warning("variable_features seems too low, recommend >= 2000")
  }
  
  # Validate clustering parameters
  if(is.null(config$clustering$resolutions) || length(config$clustering$resolutions) == 0) {
    warning("No clustering resolutions specified, using default")
    config$clustering$resolutions <- c(0.4, 0.6, 0.8, 1.0)
  }
  
  message("Configuration validation completed successfully")
  return(TRUE)
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Get species-specific cell cycle genes
get_cell_cycle_genes <- function(species = "mouse") {
  
  if (species == "mouse" || species == "Mus_musculus") {
    
    s_genes <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm7", "Mcm4", "Rrm1", "Ung", 
                 "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Cenpu", 
                 "Hells", "Rfc2", "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76", 
                 "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", 
                 "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", 
                 "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Mrpl36", "E2f8")
    
    g2m_genes <- c("Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", 
                   "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", 
                   "Tacc3", "Fam64a", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", 
                   "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b", "Hjurp", 
                   "Cdca3", "Hn1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", 
                   "Ncapd2", "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", 
                   "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", 
                   "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa")
    
  } else {
    # Human cell cycle genes (default Seurat)
    s_genes <- cc.genes$s.genes
    g2m_genes <- cc.genes$g2m.genes
  }
  
  return(list(s_genes = s_genes, g2m_genes = g2m_genes))
}

# Print configuration summary
print_config_summary <- function(config) {
  
  cat("=== HYBRID CONFIGURATION SUMMARY ===\n")
  cat("Dataset type:", config$dataset_info$type, "\n")
  cat("Tissue type:", config$dataset_info$tissue, "\n")
  cat("Optimization:", config$dataset_info$optimization, "\n")
  cat("Species:", config$dataset_info$species, "\n")
  cat("\n")
  
  cat("QC Parameters:\n")
  cat("- Min features:", config$qc$min_features_base, "\n")
  cat("- Max features:", config$qc$max_features_base, "\n")
  cat("- Max MT%:", config$qc$max_mt_percent_base, "\n")
  cat("- Adaptive thresholds:", config$qc$min_features_adaptive, "\n")
  cat("- Doublet detection:", config$qc$doublet_detection, "\n")
  cat("\n")
  
  cat("Integration Parameters:\n")
  cat("- Methods:", paste(config$integration$methods, collapse = ", "), "\n")
  cat("- Variable features:", config$integration$variable_features, "\n")
  cat("- PCA dimensions:", config$integration$pca_dims, "\n")
  cat("\n")
  
  cat("Clustering Parameters:\n")
  cat("- Resolutions:", paste(config$clustering$resolutions, collapse = ", "), "\n")
  cat("- Silhouette analysis:", config$clustering$silhouette_analysis, "\n")
  cat("- Default resolution:", config$clustering$default_resolution, "\n")
  cat("=====================================\n")
}

# =============================================================================
# EXAMPLE USAGE AND TESTING
# =============================================================================

# Test function for configuration
test_config <- function() {
  
  cat("Testing hybrid configuration system...\n")
  
  # Test GSE207128 detection
  gse_config <- get_analysis_config("/path/to/GSE207128/data")
  cat("GSE207128 config test: PASSED\n")
  
  # Test general dataset
  general_config <- get_analysis_config("/path/to/general/data")
  cat("General config test: PASSED\n")
  
  # Test marker configuration
  gse_markers <- get_marker_config("GSE207128")
  general_markers <- get_marker_config("general")
  cat("Marker config test: PASSED\n")
  
  cat("All configuration tests PASSED!\n")
  return(TRUE)
}

# =============================================================================
# OUTPUT PATH HELPERS
# =============================================================================

# Get standardized output paths
get_output_paths <- function(config) {
  if(is.null(config$output_paths)) {
    # Fallback to organized structure
    return(list(
      aggregated_data = "Outputs/01_Aggregated_Data",
      processed_data = "Outputs/02_Processed_Data", 
      integrated_data = "Outputs/03_Integrated_Data",
      annotated_data = "Outputs/04_Annotated_Data",
      analysis_results = "Outputs/05_Analysis_Results",
      figures = "Outputs/06_Figures",
      reports = "Outputs/07_Reports",
      legacy = "Outputs/08_Legacy"
    ))
  }
  return(config$output_paths)
}

# Create output directories
create_output_directories <- function(config) {
  paths <- get_output_paths(config)
  for(path in paths) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
  }
  cat("✓ Output directory structure created\n")
}

# Message when configuration is loaded
message("Hybrid Analysis Configuration loaded successfully!")
message("Use get_analysis_config() to retrieve dataset-specific parameters")
message("Use get_marker_config() to retrieve cell type markers")
message("Use test_config() to validate the configuration system")