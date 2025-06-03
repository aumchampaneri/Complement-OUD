# ==============================================================================
# Complete Analysis Pipeline with Visualizations
# ==============================================================================
# Purpose: Run the complete pipeline from data loading to advanced visualizations
# Features: Automated execution, error handling, comprehensive outputs
# ==============================================================================

# Configuration
BASE_DIR <- "/Users/aumchampaneri/Complement-OUD/Translational_Study/Homo_sapiens"
SCRIPTS_DIR <- file.path(BASE_DIR, "Scripts")

# Set flag to prevent individual script execution
SOURCED <- TRUE

# ==============================================================================
# PIPELINE EXECUTION FUNCTIONS
# ==============================================================================

#' Run complete analysis pipeline
run_complete_pipeline <- function() {
  cat("="*80, "\n")
  cat("COMPLETE NEUROINFLAMMATORY OUD ANALYSIS PIPELINE\n")
  cat("="*80, "\n")
  
  start_time <- Sys.time()
  
  tryCatch({
    # Step 1: Load and analyze GSE174409 (Bulk RNA-seq)
    cat("\n🧬 STEP 1: GSE174409 Analysis (Bulk RNA-seq)\n")
    cat("-"*50, "\n")
    source(file.path(SCRIPTS_DIR, "00_GSE174409_Paired-Analysis.R"))
    cat("✓ GSE174409 analysis complete\n")
    
    # Step 2: Load and analyze GSE225158 (snRNA-seq pseudobulk)
    cat("\n🔬 STEP 2: GSE225158 Analysis (snRNA-seq)\n")
    cat("-"*50, "\n")
    
    # First run Python data preparation if needed
    python_script <- file.path(SCRIPTS_DIR, "00_GSE225158_Python-Loading.py")
    if (file.exists(python_script)) {
      cat("Running Python data preparation...\n")
      system(paste("python", python_script))
    }
    
    # Then run R analysis
    source(file.path(SCRIPTS_DIR, "00_GSE225158_Paired-Pseudobulk.R"))
    cat("✓ GSE225158 analysis complete\n")
    
    # Step 3: Cross-dataset integration
    cat("\n🔗 STEP 3: Cross-Dataset Integration\n")
    cat("-"*50, "\n")
    source(file.path(SCRIPTS_DIR, "01_Cross-Dataset_Integration.R"))
    cat("✓ Integration analysis complete\n")
    
    # Step 4: Comprehensive neuroinflammatory analysis
    cat("\n🧠 STEP 4: Neuroinflammatory Pathway Analysis\n")
    cat("-"*50, "\n")
    source(file.path(SCRIPTS_DIR, "02_Neuroinflammatory_Analysis.R"))
    cat("✓ Neuroinflammatory analysis complete\n")
    
    # Step 5: Advanced visualizations
    cat("\n🎨 STEP 5: Advanced Visualizations\n")
    cat("-"*50, "\n")
    source(file.path(SCRIPTS_DIR, "06_Advanced_Visualizations.R"))
    cat("✓ Advanced visualizations complete\n")
    
    # Summary
    end_time <- Sys.time()
    runtime <- difftime(end_time, start_time, units = "mins")
    
    cat("\n", "="*80, "\n")
    cat("🎉 PIPELINE COMPLETE! 🎉\n")
    cat("Total runtime:", round(runtime, 2), "minutes\n")
    cat("="*80, "\n")
    
    # List all outputs
    list_all_outputs()
    
  }, error = function(e) {
    cat("\n❌ PIPELINE ERROR:", e$message, "\n")
    cat("Please check individual scripts for detailed error messages\n")
    stop(e)
  })
}

#' List all generated outputs
list_all_outputs <- function() {
  cat("\n📁 GENERATED OUTPUTS:\n")
  
  outputs_base <- file.path(BASE_DIR, "Outputs")
  results_base <- file.path(BASE_DIR, "Results")
  
  if (dir.exists(outputs_base)) {
    cat("\n🎨 VISUALIZATIONS & FIGURES:\n")
    output_dirs <- list.dirs(outputs_base, recursive = FALSE)
    for (dir in output_dirs) {
      cat("  📊", basename(dir), "\n")
      files <- list.files(dir, recursive = TRUE, pattern = "\\.(png|html|pdf)$")
      if (length(files) > 0) {
        for (file in head(files, 5)) {
          cat("    -", file, "\n")
        }
        if (length(files) > 5) {
          cat("    ... and", length(files) - 5, "more files\n")
        }
      }
    }
  }
  
  if (dir.exists(results_base)) {
    cat("\n📊 DATA & RESULTS:\n")
    result_dirs <- list.dirs(results_base, recursive = FALSE)
    for (dir in result_dirs) {
      cat("  📁", basename(dir), "\n")
      files <- list.files(dir, pattern = "\\.(rds|csv)$")
      if (length(files) > 0) {
        for (file in head(files, 3)) {
          cat("    -", file, "\n")
        }
        if (length(files) > 3) {
          cat("    ... and", length(files) - 3, "more files\n")
        }
      }
    }
  }
  
  cat("\n💡 TIP: Check the HTML volcano plots for interactive pathway exploration!\n")
}

# ==============================================================================
# QUICK ANALYSIS FUNCTIONS
# ==============================================================================

#' Run just the visualization pipeline (if data already exists)
run_visualizations_only <- function() {
  cat("🎨 Running Advanced Visualizations Only\n")
  cat("-"*50, "\n")
  
  SOURCED <- TRUE
  source(file.path(SCRIPTS_DIR, "06_Advanced_Visualizations.R"))
  
  cat("✓ Visualizations complete!\n")
}

#' Run just the neuroinflammatory analysis (to fix expression data)
run_neuroinflammatory_only <- function() {
  cat("🧠 Running Neuroinflammatory Analysis Only\n")
  cat("-"*50, "\n")
  
  SOURCED <- TRUE
  source(file.path(SCRIPTS_DIR, "02_Neuroinflammatory_Analysis.R"))
  
  cat("✓ Neuroinflammatory analysis complete!\n")
}

#' Quick status check of all analyses
check_analysis_status <- function() {
  cat("📋 ANALYSIS STATUS CHECK\n")
  cat("="*40, "\n")
  
  # Check for key result files
  checks <- list(
    "GSE174409 Analysis" = file.path(BASE_DIR, "Results", "GSE174409", "GSE174409_region_analysis_results.rds"),
    "GSE225158 Analysis" = file.path(BASE_DIR, "Results", "GSE225158", "GSE225158_region_analysis_results.rds"),
    "Integration Analysis" = file.path(BASE_DIR, "Results", "Integration", "integration_summary.csv"),
    "Neuroinflammatory Analysis" = file.path(BASE_DIR, "Results", "Neuroinflammatory_Analysis", "comprehensive_neuroinflammatory_analysis_enhanced.rds"),
    "Advanced Visualizations" = file.path(BASE_DIR, "Outputs", "Neuroinflammatory_Analysis", "Advanced_Visualizations", "Combined_Dashboard", "combined_analysis_dashboard.png")
  )
  
  for (name in names(checks)) {
    status <- if (file.exists(checks[[name]])) "✅ Complete" else "❌ Missing"
    cat(sprintf("%-25s: %s\n", name, status))
  }
  
  cat("\n💡 Run run_complete_pipeline() to generate missing components\n")
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if (!exists("SOURCED")) {
  cat("🚀 Starting Complete Analysis Pipeline...\n")
  cat("This will run all analysis steps sequentially.\n")
  cat("Estimated time: 10-15 minutes\n\n")
  
  run_complete_pipeline()
}
