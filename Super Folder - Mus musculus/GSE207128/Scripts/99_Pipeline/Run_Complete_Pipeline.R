# =============================================================================
# Enhanced GSE207128 Analysis Pipeline - Organized Structure v2.0
# =============================================================================

print("=== ENHANCED GSE207128 ANALYSIS PIPELINE ===")
print("Professional Organized Structure - Version 2.0")
print(paste("Started at:", Sys.time()))

# Set working directory
setwd("/Users/aumchampaneri/Complement-OUD/Super Folder - Mus musculus/GSE207128/")

# Create pipeline session directory
pipeline_start <- Sys.time()
pipeline_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
session_dir <- paste0("Outputs/07_Reports/Pipeline_Session_", pipeline_id)
dir.create(session_dir, showWarnings = FALSE, recursive = TRUE)

# Setup logging
log_file <- file.path(session_dir, "pipeline_log.txt")
cat("Pipeline started at:", as.character(pipeline_start), "\n", file = log_file)

cat("=== PIPELINE SESSION INFO ===\n")
cat("Session ID:", pipeline_id, "\n")
cat("Log file:", log_file, "\n")
cat("Results will be saved in organized Outputs/ folders\n")

tryCatch({

  # Step 0: Configuration Check
  print("\n=== STEP 0: CONFIGURATION CHECK ===")
  source("Scripts/00_Configuration/Hybrid_Analysis_Config.R")
  config <- get_analysis_config(getwd())
  print("✓ Configuration loaded successfully")
  cat("Configuration check completed\n", file = log_file, append = TRUE)

  # Step 1: Data Aggregation (if needed)
  if(!file.exists("Outputs/01_Aggregated_Data/combined_seurat.rds")) {
    print("\n=== STEP 1: DATA AGGREGATION ===")
    step1_start <- Sys.time()
    source("Scripts/01_Data_Processing/01_Data_Aggregation.R")
    step1_time <- round(as.numeric(difftime(Sys.time(), step1_start, units = "mins")), 2)
    print(paste("✓ Step 1 completed in", step1_time, "minutes"))
    cat("Step 1 completed in", step1_time, "minutes\n", file = log_file, append = TRUE)
  } else {
    print("\n=== STEP 1: SKIPPED (Data already aggregated) ===")
    cat("Step 1 skipped - data already aggregated\n", file = log_file, append = TRUE)
  }

  # Step 2: QC and Processing
  print("\n=== STEP 2: QC AND PROCESSING ===")
  step2_start <- Sys.time()
  source("Scripts/01_Data_Processing/02_QC_and_Processing.R")
  step2_time <- round(as.numeric(difftime(Sys.time(), step2_start, units = "mins")), 2)
  print(paste("✓ Step 2 completed in", step2_time, "minutes"))
  cat("Step 2 completed in", step2_time, "minutes\n", file = log_file, append = TRUE)

  # Step 3: Integration
  print("\n=== STEP 3: INTEGRATION ===")
  step3_start <- Sys.time()
  source("Scripts/01_Data_Processing/03_Integration.R")
  step3_time <- round(as.numeric(difftime(Sys.time(), step3_start, units = "mins")), 2)
  print(paste("✓ Step 3 completed in", step3_time, "minutes"))
  cat("Step 3 completed in", step3_time, "minutes\n", file = log_file, append = TRUE)

  # Step 4: Cell Type Annotation
  print("\n=== STEP 4: CELL TYPE ANNOTATION ===")
  step4_start <- Sys.time()
  source("Scripts/01_Data_Processing/04_CellType_Annotation.R")
  step4_time <- round(as.numeric(difftime(Sys.time(), step4_start, units = "mins")), 2)
  print(paste("✓ Step 4 completed in", step4_time, "minutes"))
  cat("Step 4 completed in", step4_time, "minutes\n", file = log_file, append = TRUE)

  # Pipeline completion summary
  total_time <- round(as.numeric(difftime(Sys.time(), pipeline_start, units = "mins")), 2)

  print("\n=== PIPELINE COMPLETED SUCCESSFULLY ===")
  print(paste("Total pipeline time:", total_time, "minutes"))
  print(paste("Session ID:", pipeline_id))
  print("📊 Results organized in:")
  print("  - Outputs/01_Aggregated_Data/     (Raw combined data)")
  print("  - Outputs/02_Processed_Data/      (QC and processed data)")
  print("  - Outputs/03_Integrated_Data/     (Integration results)")
  print("  - Outputs/04_Annotated_Data/      (Cell type annotations)")
  print("  - Outputs/07_Reports/             (Pipeline reports)")
  print(paste("Finished at:", Sys.time()))

  # Save completion summary
  completion_summary <- list(
    pipeline_version = "2.0_organized",
    session_id = pipeline_id,
    start_time = pipeline_start,
    end_time = Sys.time(),
    total_duration_minutes = total_time,
    steps_completed = c("Configuration", "Data Aggregation", "QC Processing", "Integration", "Cell Type Annotation"),
    outputs_location = "Organized in Outputs/ subdirectories"
  )

  saveRDS(completion_summary, file.path(session_dir, "pipeline_completion_summary.rds"))
  cat("\nPipeline completed successfully in", total_time, "minutes\n", file = log_file, append = TRUE)

}, error = function(e) {
  print(paste("ERROR in pipeline:", e$message))
  print("Check individual script outputs and log file for details.")
  print(traceback())
  cat("Pipeline error:", e$message, "\n", file = log_file, append = TRUE)

}, finally = {
  cat("Pipeline session ended at:", as.character(Sys.time()), "\n", file = log_file, append = TRUE)
})

print("\n🎉 ENHANCED PIPELINE COMPLETE!")
print("📁 Check organized Outputs/ folders for results")
print(paste("📝 Session log:", log_file))
