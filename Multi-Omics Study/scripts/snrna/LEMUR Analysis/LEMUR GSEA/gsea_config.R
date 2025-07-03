# =============================================================================
# GSEA Analysis Configuration File
# =============================================================================
#
# Purpose: Centralized configuration for GSEA pathway analysis
# Author: Multi-Omics OUD Study
# Date: 2024
# =============================================================================

# =============================================================================
# PROJECT PATHS
# =============================================================================

# Base project directory
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/scripts/snrna/LEMUR Analysis"

# Input and output paths
PATHS <- list(
  input_file = file.path(BASE_DIR, "outputs/tables/top_de_genes.csv"),
  output_dir = file.path(BASE_DIR, "LEMUR GSEA/outputs"),
  plot_dir = file.path(BASE_DIR, "LEMUR GSEA/outputs/plots"),
  table_dir = file.path(BASE_DIR, "LEMUR GSEA/outputs/tables"),
  report_dir = file.path(BASE_DIR, "LEMUR GSEA/outputs/reports"),
  log_dir = file.path(BASE_DIR, "LEMUR GSEA/outputs/logs")
)

# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================

# Statistical thresholds
STATS <- list(
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.2,
  fdr_threshold = 0.05,
  min_gs_size = 10,
  max_gs_size = 500,
  eps = 1e-10
)

# Gene filtering parameters
GENE_FILTER <- list(
  effect_size_threshold = 0.1,
  min_genes_for_analysis = 5,
  conversion_rate_threshold = 0.5 # Minimum 50% gene ID conversion rate
)

# =============================================================================
# DATABASE CONFIGURATION
# =============================================================================

# MSigDB collections to analyze
MSIGDB_COLLECTIONS <- list(
  hallmark = list(
    collection = "H",
    subcollection = NULL,
    name = "Hallmark Pathways",
    description = "Hallmark gene sets from MSigDB"
  ),
  canonical = list(
    collection = "C2",
    subcollection = "CP",
    name = "Canonical Pathways",
    description = "Canonical pathways from MSigDB"
  ),
  go_bp = list(
    collection = "C5",
    subcollection = "GO:BP",
    name = "GO Biological Process",
    description = "Gene Ontology Biological Process terms"
  ),
  go_cc = list(
    collection = "C5",
    subcollection = "GO:CC",
    name = "GO Cellular Component",
    description = "Gene Ontology Cellular Component terms"
  ),
  go_mf = list(
    collection = "C5",
    subcollection = "GO:MF",
    name = "GO Molecular Function",
    description = "Gene Ontology Molecular Function terms"
  )
)

# Additional pathway databases
PATHWAY_DATABASES <- list(
  kegg = list(
    enabled = TRUE,
    organism = "hsa",
    name = "KEGG Pathways",
    description = "Kyoto Encyclopedia of Genes and Genomes"
  ),
  reactome = list(
    enabled = TRUE,
    name = "Reactome Pathways",
    description = "Reactome pathway database"
  ),
  disease_ontology = list(
    enabled = FALSE, # Currently not working
    ontology = "DO",
    name = "Disease Ontology",
    description = "Disease-associated gene sets"
  )
)

# Enrichr databases (web-based analysis)
ENRICHR_DATABASES <- c(
  "GO_Biological_Process_2023",
  "GO_Molecular_Function_2023",
  "KEGG_2021_Human",
  "Reactome_2022",
  "WikiPathway_2023_Human",
  "BioPlanet_2019"
)

# =============================================================================
# NEUROSCIENCE FILTERING KEYWORDS
# =============================================================================

# Keywords for filtering neuroscience/addiction relevant pathways
NEURO_KEYWORDS <- c(
  # Core neuroscience terms
  "synaptic", "synapse", "neuron", "neural", "neuronal",
  "axon", "dendrite", "spine", "vesicle", "presynaptic", "postsynaptic",

  # Neurotransmitter systems
  "dopamine", "serotonin", "GABA", "glutamate", "acetylcholine",
  "norepinephrine", "noradrenaline", "glycine",

  # Addiction and reward
  "addiction", "reward", "reinforcement", "motivation",
  "craving", "withdrawal", "tolerance", "dependence",

  # Plasticity and learning
  "plasticity", "learning", "memory", "LTP", "LTD",
  "potentiation", "depression",

  # Ion channels and signaling
  "channel", "receptor", "calcium", "potassium", "sodium",
  "voltage", "ligand", "AMPA", "NMDA", "GABA_A", "GABA_B",

  # Neurotransmitter signaling
  "neurotransmitter", "signal", "transduction", "second.messenger",
  "cAMP", "cGMP", "PKA", "PKC",

  # Cellular processes
  "excitatory", "inhibitory", "membrane.potential", "action.potential",
  "calcium.signaling", "protein.kinase"
)

# =============================================================================
# PLOTTING CONFIGURATION
# =============================================================================

# Plot parameters
PLOT_PARAMS <- list(
  # Basic settings
  width = 12,
  height = 8,
  dpi = 300,
  format = "png",

  # Display limits
  max_pathways_plot = 20,
  max_pathway_name_length = 50,
  min_pathways_for_network = 5,

  # Color schemes
  continuous_palette = "viridis",
  categorical_palette = "Set2",
  diverging_palette = "RdBu",

  # Point and text sizes
  point_size_range = c(2, 8),
  text_size_title = 14,
  text_size_subtitle = 12,
  text_size_axis = 10,
  text_size_label = 9
)

# Plot types to generate
PLOT_TYPES <- list(
  # Individual database plots
  enhanced_dotplot = TRUE,
  enhanced_barplot = TRUE,
  volcano_plot = TRUE,
  network_plot = TRUE,

  # Summary plots
  database_summary = TRUE,
  pathway_categories = TRUE,
  top_pathways_combined = TRUE,

  # Advanced plots
  heatmap_summary = FALSE, # Currently disabled
  interactive_plots = FALSE # Future feature
)

# =============================================================================
# ANALYSIS WORKFLOW CONFIGURATION
# =============================================================================

# Analysis steps to run
WORKFLOW <- list(
  # Setup and validation
  validate_setup = TRUE,
  validate_input_data = TRUE,

  # Gene processing
  convert_gene_ids = TRUE,
  filter_gene_list = TRUE,

  # Database preparation
  load_msigdb = TRUE,
  filter_neuro_pathways = TRUE,

  # Analysis methods
  run_gsea = TRUE, # Currently has issues - may skip
  run_ora = TRUE,
  run_kegg = TRUE,
  run_reactome = TRUE,
  run_enrichr = FALSE, # Requires internet, optional

  # Visualization
  create_plots = TRUE,
  create_enhanced_plots = TRUE,

  # Reporting
  generate_reports = TRUE,
  save_workspace = TRUE
)

# =============================================================================
# PACKAGE REQUIREMENTS
# =============================================================================

# Required packages by category
REQUIRED_PACKAGES <- list(
  # Core analysis packages (Bioconductor)
  bioc = c(
    "clusterProfiler", "org.Hs.eg.db", "msigdbr", "enrichplot",
    "DOSE", "ReactomePA", "pathview"
  ),

  # Data manipulation and visualization (CRAN)
  cran = c(
    "ggplot2", "dplyr", "readr", "stringr", "tidyr",
    "RColorBrewer", "viridis", "scales", "ggrepel"
  ),

  # Optional packages
  optional = c(
    "httr", "jsonlite", # For Enrichr web API
    "pheatmap", "VennDiagram", # Additional plotting
    "plotly", "DT", # Interactive features
    "knitr", "rmarkdown" # Report generation
  )
)

# =============================================================================
# PATHWAY CATEGORIES FOR ANALYSIS
# =============================================================================

# Functional categories for pathway grouping
PATHWAY_CATEGORIES <- list(
  "Synaptic Function" = c("synap", "vesicle", "neurotransmitter", "axon", "dendrite"),
  "Addiction/Reward" = c("dopamine", "reward", "addiction", "substance", "reinforcement"),
  "Inflammation" = c("interleukin", "immune", "inflammatory", "cytokine", "interferon"),
  "Cell Adhesion" = c("adhesion", "junction", "cadherin", "cell.cell", "tight.junction"),
  "Ion Channels" = c("channel", "calcium", "potassium", "sodium", "membrane.potential"),
  "Transcription" = c("transcription", "stat", "jak", "gene.expression", "nuclear"),
  "Metabolism" = c("metabolic", "glycol", "lipid", "amino.acid", "energy"),
  "Plasticity" = c("plasticity", "learning", "memory", "LTP", "potentiation"),
  "Apoptosis" = c("apoptosis", "cell.death", "programmed.cell.death", "caspase")
)

# =============================================================================
# OUTPUT FORMATTING
# =============================================================================

# Report generation settings
REPORT_SETTINGS <- list(
  include_methods = TRUE,
  include_plots = TRUE,
  include_gene_lists = TRUE,
  max_pathways_in_report = 50,
  significant_only = TRUE,
  format = "markdown"
)

# File naming conventions
FILE_NAMING <- list(
  timestamp_format = "%Y%m%d_%H%M%S",
  use_timestamps = TRUE,
  prefix = "gsea",
  separator = "_"
)

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

# Logging settings
LOGGING <- list(
  log_level = "INFO", # DEBUG, INFO, WARNING, ERROR
  log_to_file = TRUE,
  log_to_console = TRUE,
  timestamp_logs = TRUE,
  max_log_size_mb = 10
)

# =============================================================================
# HELPER FUNCTION TO VALIDATE CONFIG
# =============================================================================

validate_config <- function() {
  cat("🔧 Validating configuration...\n")

  # Check if base directory exists
  if (!dir.exists(BASE_DIR)) {
    stop("❌ Base directory does not exist: ", BASE_DIR)
  }

  # Check statistical parameters
  if (STATS$pvalue_cutoff <= 0 || STATS$pvalue_cutoff >= 1) {
    stop("❌ Invalid p-value cutoff: ", STATS$pvalue_cutoff)
  }

  # Validate gene filter parameters
  if (GENE_FILTER$effect_size_threshold < 0) {
    stop("❌ Effect size threshold must be non-negative")
  }

  # Check plot parameters
  if (PLOT_PARAMS$width <= 0 || PLOT_PARAMS$height <= 0) {
    stop("❌ Invalid plot dimensions")
  }

  cat("✅ Configuration validation passed\n")
  return(TRUE)
}

# =============================================================================
# CONFIGURATION SUMMARY
# =============================================================================

print_config_summary <- function() {
  cat("📋 GSEA Analysis Configuration Summary\n")
  cat("=====================================\n")
  cat("📁 Base Directory:", BASE_DIR, "\n")
  cat("📊 P-value Cutoff:", STATS$pvalue_cutoff, "\n")
  cat("🧬 Min Gene Set Size:", STATS$min_gs_size, "\n")
  cat("🗃️ MSigDB Collections:", length(MSIGDB_COLLECTIONS), "\n")
  cat("🔍 Neuro Keywords:", length(NEURO_KEYWORDS), "\n")
  cat("🎨 Plot Types:", sum(unlist(PLOT_TYPES)), "enabled\n")
  cat(
    "📦 Required Packages:",
    length(REQUIRED_PACKAGES$bioc) + length(REQUIRED_PACKAGES$cran), "\n"
  )
  cat("=====================================\n")
}

# Automatically validate configuration when loaded
if (interactive()) {
  validate_config()
  print_config_summary()
}
