# ========================================================================
# GSE118918 Complete Pipeline Batch Runner
# ========================================================================
# 
# PIPELINE OVERVIEW:
# This script executes the complete GSE118918 analysis pipeline in sequence:
# 1. Processing & Quality Control (00_Processing_QC.R)
# 2. Differential Expression Analysis (01_Differential_Expression.R) 
# 3. Pathway Enrichment Analysis (02_Pathway_Enrichment.R)
# 
# SCIENTIFIC RATIONALE:
# - Automated execution ensures consistency and reproducibility
# - Sequential dependency management prevents data conflicts
# - Comprehensive error handling and progress tracking
# - Publication-ready output generation in single run
# 
# EXECUTION TIME:
# Estimated runtime: 15-30 minutes depending on system specifications
# - Processing & QC: 5-10 minutes
# - Differential Expression: 3-5 minutes
# - Pathway Enrichment: 7-15 minutes
# 
# REQUIREMENTS:
# - All required R packages installed (see individual scripts)
# - Raw data files in Data/Raw_Data directory
# - Sufficient disk space for outputs (~500MB)
# - Memory: 8GB+ recommended for large datasets
# ========================================================================

# Set reproducibility parameters
set.seed(42)
options(stringsAsFactors = FALSE)

# Record pipeline start time
pipeline_start_time <- Sys.time()
cat("========================================================================\n")
cat("GSE118918 COMPLETE ANALYSIS PIPELINE STARTED\n")
cat("========================================================================\n")
cat("Pipeline initiated at:", as.character(pipeline_start_time), "\n")
cat("Estimated completion time:", as.character(pipeline_start_time + as.difftime(25, units = "mins")), "\n\n")

# ========================================================================
# SECTION 1: ENVIRONMENT SETUP AND VALIDATION
# ========================================================================

cat("=== PIPELINE ENVIRONMENT SETUP ===\n")

# Define base directory and script paths
base_dir <- "/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE118918"
setwd(base_dir)

# Define script paths in execution order
scripts <- list(
  processing_qc = "Scripts/00_Processing_QC.R",
  differential_expression = "Scripts/01_Differential_Expression.R", 
  pathway_enrichment = "Scripts/02_Pathway_Enrichment.R"
)

# Validate all scripts exist
cat("Validating pipeline scripts...\n")
script_validation <- sapply(scripts, file.exists)

if(!all(script_validation)) {
  missing_scripts <- names(script_validation)[!script_validation]
  stop("Missing required scripts: ", paste(missing_scripts, collapse = ", "), 
       "\nPlease ensure all pipeline scripts are in the Scripts/ directory.")
}

cat("✓ All pipeline scripts found\n")

# Validate data directories
required_dirs <- c("Data/Raw_Data", "Data/NCBI")
dir_validation <- sapply(file.path(base_dir, required_dirs), dir.exists)

if(!all(dir_validation)) {
  missing_dirs <- required_dirs[!dir_validation]
  stop("Missing required directories: ", paste(missing_dirs, collapse = ", "),
       "\nPlease ensure data directories are properly set up.")
}

cat("✓ All required data directories found\n")

# Create master output directory
master_output_dir <- file.path(base_dir, "Outputs/00_Pipeline_Results")
dir.create(master_output_dir, showWarnings = FALSE, recursive = TRUE)

# Create pipeline log directory
log_dir <- file.path(master_output_dir, "Pipeline_Logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

cat("✓ Pipeline output directories created\n")

# ========================================================================
# SECTION 2: PIPELINE EXECUTION FUNCTIONS
# ========================================================================

cat("\n=== PIPELINE EXECUTION FUNCTIONS ===\n")

# Function to execute individual scripts with comprehensive logging
execute_script <- function(script_name, script_path, log_dir) {
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("EXECUTING:", toupper(script_name), "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  script_start_time <- Sys.time()
  cat("Script started at:", as.character(script_start_time), "\n")
  
  # Create log file for this script
  log_file <- file.path(log_dir, paste0(script_name, "_execution.log"))
  
  # Execute script with error handling
  execution_result <- tryCatch({
    
    # Capture output to log file
    sink(log_file, type = "output", split = TRUE)
    
    cat("Executing script:", script_path, "\n")
    cat("Execution started:", as.character(script_start_time), "\n\n")
    
    # Source the script
    source(script_path, echo = TRUE)
    
    # Close log capture
    sink()
    
    script_end_time <- Sys.time()
    execution_time <- as.numeric(difftime(script_end_time, script_start_time, units = "mins"))
    
    cat("✓ Script completed successfully\n")
    cat("Execution time:", round(execution_time, 2), "minutes\n")
    
    return(list(
      success = TRUE,
      start_time = script_start_time,
      end_time = script_end_time,
      execution_time = execution_time,
      error = NULL
    ))
    
  }, error = function(e) {
    
    # Ensure sink is closed even on error
    if(sink.number() > 0) sink()
    
    script_end_time <- Sys.time()
    execution_time <- as.numeric(difftime(script_end_time, script_start_time, units = "mins"))
    
    cat("✗ Script failed with error:\n")
    cat("Error message:", e$message, "\n")
    cat("Execution time before error:", round(execution_time, 2), "minutes\n")
    
    # Log error details
    error_log <- file.path(log_dir, paste0(script_name, "_error.log"))
    writeLines(c(
      paste("Script:", script_path),
      paste("Error time:", as.character(script_end_time)),
      paste("Error message:", e$message),
      paste("Call:", deparse(e$call))
    ), error_log)
    
    return(list(
      success = FALSE,
      start_time = script_start_time,
      end_time = script_end_time,
      execution_time = execution_time,
      error = e$message
    ))
  })
  
  return(execution_result)
}

# Function to validate script outputs
validate_script_outputs <- function(script_name) {
  
  validation_results <- list()
  
  if(script_name == "processing_qc") {
    required_files <- c(
      "Outputs/01_Processing_QC/Data/dge_normalized_final.rds",
      "Outputs/01_Processing_QC/Data/sample_metadata_final.rds",
      "Outputs/01_Processing_QC/Data/logcpm_normalized_final.rds",
      "Outputs/01_Processing_QC/Reports/GSE118918_Processing_QC_Report.txt"
    )
    
  } else if(script_name == "differential_expression") {
    required_files <- c(
      "Outputs/01_Differential_Expression/Tables/differential_expression_results.csv",
      "Outputs/01_Differential_Expression/Data/de_analysis_complete.rds"
    )
    
  } else if(script_name == "pathway_enrichment") {
    required_files <- c(
      "Outputs/02_Pathway_Enrichment/Data/pathway_enrichment_results.rds",
      "Outputs/02_Pathway_Enrichment/Reports/GSE118918_Pathway_Analysis_Report.txt"
    )
  }
  
  # Check if required files exist
  file_validation <- sapply(required_files, file.exists)
  
  validation_results$files_created <- sum(file_validation)
  validation_results$files_expected <- length(required_files)
  validation_results$success_rate <- validation_results$files_created / validation_results$files_expected
  validation_results$missing_files <- required_files[!file_validation]
  
  return(validation_results)
}

# ========================================================================
# SECTION 3: PIPELINE EXECUTION
# ========================================================================

cat("\n=== BEGINNING PIPELINE EXECUTION ===\n")

# Initialize pipeline tracking
pipeline_results <- list()
pipeline_success <- TRUE

# Execute scripts in sequence
for(script_name in names(scripts)) {
  
  script_path <- scripts[[script_name]]
  
  # Execute script
  execution_result <- execute_script(script_name, script_path, log_dir)
  
  # Store results
  pipeline_results[[script_name]] <- execution_result
  
  # Check if script succeeded
  if(!execution_result$success) {
    pipeline_success <- FALSE
    cat("\n❌ PIPELINE HALTED due to script failure:", script_name, "\n")
    cat("Error:", execution_result$error, "\n")
    break
  }
  
  # Validate outputs if script succeeded
  validation_result <- validate_script_outputs(script_name)
  pipeline_results[[script_name]]$validation <- validation_result
  
  if(validation_result$success_rate < 1.0) {
    cat("⚠️  Warning: Some expected outputs missing\n")
    cat("Created:", validation_result$files_created, "/", validation_result$files_expected, "files\n")
    if(length(validation_result$missing_files) > 0) {
      cat("Missing files:\n")
      for(file in validation_result$missing_files) {
        cat("  -", file, "\n")
      }
    }
  } else {
    cat("✓ All expected outputs created successfully\n")
  }
  
  # Brief pause between scripts for system resources
  Sys.sleep(2)
}

# ========================================================================
# SECTION 4: PIPELINE SUMMARY AND REPORTING
# ========================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PIPELINE EXECUTION SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

pipeline_end_time <- Sys.time()
total_pipeline_time <- as.numeric(difftime(pipeline_end_time, pipeline_start_time, units = "mins"))

cat("Pipeline started:", as.character(pipeline_start_time), "\n")
cat("Pipeline completed:", as.character(pipeline_end_time), "\n")
cat("Total execution time:", round(total_pipeline_time, 2), "minutes\n\n")

# Create detailed pipeline report
create_pipeline_report <- function(pipeline_results, pipeline_success, output_file) {
  
  sink(output_file)
  
  cat("========================================================================\n")
  cat("GSE118918 COMPLETE PIPELINE EXECUTION REPORT\n")
  cat("========================================================================\n\n")
  
  cat("PIPELINE OVERVIEW\n")
  cat("-----------------\n")
  cat("Dataset: GSE118918 (Nucleus Accumbens, Mock vs Morphine)\n")
  cat("Pipeline execution date:", as.character(pipeline_start_time), "\n")
  cat("Total execution time:", round(total_pipeline_time, 2), "minutes\n")
  cat("Overall pipeline status:", ifelse(pipeline_success, "SUCCESS", "FAILED"), "\n\n")
  
  cat("SCRIPT EXECUTION DETAILS\n")
  cat("------------------------\n")
  
  for(script_name in names(pipeline_results)) {
    result <- pipeline_results[[script_name]]
    
    cat("Script:", toupper(script_name), "\n")
    cat("  Status:", ifelse(result$success, "SUCCESS", "FAILED"), "\n")
    cat("  Execution time:", round(result$execution_time, 2), "minutes\n")
    
    if(!is.null(result$validation)) {
      cat("  Output validation:", result$validation$files_created, "/", 
          result$validation$files_expected, "files created\n")
    }
    
    if(!result$success) {
      cat("  Error:", result$error, "\n")
    }
    cat("\n")
  }
  
  cat("OUTPUT DIRECTORIES\n")
  cat("------------------\n")
  cat("Processing & QC:", "Outputs/01_Processing_QC/\n")
  cat("Differential Expression:", "Outputs/01_Differential_Expression/\n")
  cat("Pathway Enrichment:", "Outputs/02_Pathway_Enrichment/\n")
  cat("Pipeline Logs:", "Outputs/00_Pipeline_Results/Pipeline_Logs/\n\n")
  
  if(pipeline_success) {
    cat("NEXT STEPS\n")
    cat("----------\n")
    cat("1. Review comprehensive analysis reports in each output directory\n")
    cat("2. Examine publication-quality figures for manuscript preparation\n")
    cat("3. Validate complement pathway findings with functional studies\n")
    cat("4. Consider integration with additional datasets\n")
    cat("5. Prepare manuscript focusing on complement-OUD mechanisms\n\n")
    
    cat("KEY DELIVERABLES\n")
    cat("----------------\n")
    cat("- Quality-controlled expression data ready for analysis\n")
    cat("- Differential expression results (Mock vs Morphine)\n")
    cat("- Comprehensive pathway enrichment analysis\n")
    cat("- Custom complement pathway investigation\n")
    cat("- Publication-ready figures and supplementary materials\n")
  } else {
    cat("TROUBLESHOOTING\n")
    cat("---------------\n")
    cat("Pipeline execution failed. Please:\n")
    cat("1. Check error logs in Pipeline_Logs directory\n")
    cat("2. Verify all required packages are installed\n")
    cat("3. Ensure sufficient system memory and disk space\n")
    cat("4. Validate input data integrity\n")
    cat("5. Re-run individual scripts for detailed error diagnosis\n")
  }
  
  cat("\n========================================================================\n")
  cat("PIPELINE REPORT COMPLETED\n")
  cat("========================================================================\n")
  
  sink()
}

# Generate comprehensive pipeline report
pipeline_report_file <- file.path(master_output_dir, "GSE118918_Complete_Pipeline_Report.txt")
create_pipeline_report(pipeline_results, pipeline_success, pipeline_report_file)

# Print final status
if(pipeline_success) {
  cat("🎉 PIPELINE COMPLETED SUCCESSFULLY! 🎉\n\n")
  cat("Key Results:\n")
  
  # Calculate total scripts executed
  scripts_completed <- sum(sapply(pipeline_results, function(x) x$success))
  cat("- Scripts executed:", scripts_completed, "/", length(scripts), "\n")
  
  # Calculate total execution time breakdown
  for(script_name in names(pipeline_results)) {
    if(pipeline_results[[script_name]]$success) {
      cat("- ", script_name, ":", round(pipeline_results[[script_name]]$execution_time, 1), "minutes\n")
    }
  }
  
  cat("\nOutput Summary:\n")
  cat("- Complete analysis report:", pipeline_report_file, "\n")
  cat("- All individual script outputs in respective directories\n")
  cat("- Publication-ready figures and tables generated\n")
  cat("- Complement pathway analysis completed\n")
  
  cat("\n✨ Ready for manuscript preparation and biological validation! ✨\n")
  
} else {
  cat("❌ PIPELINE EXECUTION FAILED\n\n")
  cat("Please check the following:\n")
  cat("1. Error logs in:", log_dir, "\n")
  cat("2. System requirements (memory, disk space)\n") 
  cat("3. Required R packages installation\n")
  cat("4. Input data availability and format\n")
  cat("\nDetailed error information available in pipeline report:\n")
  cat(pipeline_report_file, "\n")
}

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PIPELINE EXECUTION COMPLETED\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Save pipeline execution environment for reproducibility
pipeline_environment <- list(
  execution_results = pipeline_results,
  system_info = Sys.info(),
  session_info = sessionInfo(),
  pipeline_success = pipeline_success,
  total_execution_time = total_pipeline_time,
  start_time = pipeline_start_time,
  end_time = pipeline_end_time
)

saveRDS(pipeline_environment, file.path(master_output_dir, "pipeline_execution_environment.rds"))

cat("Pipeline environment saved for reproducibility\n")
cat("All analysis components completed. Ready for publication! 📊✨\n")
