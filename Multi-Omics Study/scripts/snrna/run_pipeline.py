#!/usr/bin/env python3
'''
🚀 snRNA-seq Processing Pipeline Runner
Sequential execution with storage cleanup optimization

Pipeline steps:
1. 01b_preprocessing_scvi.py     - scVI preprocessing & batch correction (Python)
2. 02b_cell_type_annotation.py  - Cell type annotation (Python)
3. 03b_differential_expression.py - Differential expression analysis (Python)
4. 04_clinical_analysis.py       - Clinical associations (Python)
5. 05_final_integration.py       - Final results integration (Python)

Features:
- Sequential execution with dependency checking
- Support for both Python (.py) and R (.R/.r) scripts
- Automatic cleanup of intermediate files
- Raw data preservation
- Storage monitoring
- Error handling and recovery
'''

import os
import sys
import shutil
import subprocess
import time
import logging
from pathlib import Path
import pandas as pd

# =======================================
# 📁 CONFIGURATION
# =======================================

# Base paths
BASE_DIR = Path("/Users/aumchampaneri/Complement-OUD/Multi-Omics Study")
SCRIPTS_DIR = BASE_DIR / "scripts/snrna"
DATA_DIR = BASE_DIR / "data"
RAW_DATA_DIR = DATA_DIR / "raw/snrna"
RESULTS_DIR = BASE_DIR / "results"

# Pipeline scripts in execution order
PIPELINE_STEPS = [
    {
        'step': '01b',
        'script': '01b_preprocessing_scvi.py',
        'script_type': 'python',  # Added script type detection
        'description': 'scVI Preprocessing & Batch Correction',
        'output_dirs': ['data/processed/snrna_scvi'],
        'preserve_files': [],  # Files to keep after cleanup
        'required_memory_gb': 32
    },
    {
        'step': '02b', 
        'script': '02b_cell_type_annotation.py',
        'script_type': 'python',
        'description': 'Cell Type Annotation',
        'output_dirs': ['data/processed/snrna_annotated'],
        'preserve_files': ['GSE225158_annotated.h5ad'],  # Keep final annotated data
        'required_memory_gb': 16
    },
    {
        'step': '03b',
        'script': '03b_differential_expression.py',
        'script_type': 'python', 
        'description': 'Differential Expression Analysis',
        'output_dirs': ['results/differential_expression'],
        'preserve_files': ['de_results_summary.csv', 'volcano_plots'],  # Keep DE results
        'required_memory_gb': 20
    },
    {
        'step': '04',
        'script': '04_clinical_analysis.py',
        'script_type': 'python',
        'description': 'Clinical Associations Analysis', 
        'output_dirs': ['results/clinical_analysis'],
        'preserve_files': ['clinical_associations.csv', 'clinical_plots'],  # Keep clinical results
        'required_memory_gb': 12
    },
    {
        'step': '05',
        'script': '05_final_integration.py',
        'script_type': 'python',
        'description': 'Final Results Integration',
        'output_dirs': ['results/final_integration'],
        'preserve_files': ['*'],  # Keep everything from final step
        'required_memory_gb': 8
    }
]

# Storage thresholds (GB)
MIN_FREE_SPACE_GB = 50  # Minimum free space required
CLEANUP_THRESHOLD_GB = 100  # Start cleanup when less than this available

# Script execution configuration
PYTHON_EXECUTABLE = sys.executable  # Use same Python as pipeline runner
R_EXECUTABLE = "Rscript"  # R script executor
R_TIMEOUT = 10800  # 3 hours for R scripts (often slower)
PYTHON_TIMEOUT = 7200  # 2 hours for Python scripts

# =======================================
# 🔧 UTILITY FUNCTIONS
# =======================================

def setup_logging():
    """Setup logging configuration"""
    log_dir = BASE_DIR / "logs"
    log_dir.mkdir(exist_ok=True)
    
    log_file = log_dir / f"pipeline_run_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    return logging.getLogger(__name__)

def detect_script_type(script_path):
    """Auto-detect script type from file extension"""
    suffix = script_path.suffix.lower()
    if suffix == '.py':
        return 'python'
    elif suffix in ['.r', '.R']:
        return 'r'
    else:
        return 'unknown'

def get_directory_size(path):
    """Get directory size in GB"""
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            if os.path.exists(filepath):
                total_size += os.path.getsize(filepath)
    return total_size / (1024**3)  # Convert to GB

def get_disk_usage(path):
    """Get disk usage information"""
    statvfs = os.statvfs(path)
    total_gb = (statvfs.f_blocks * statvfs.f_frsize) / (1024**3)
    free_gb = (statvfs.f_bavail * statvfs.f_frsize) / (1024**3)
    used_gb = total_gb - free_gb
    return total_gb, used_gb, free_gb

def check_storage_space(logger, required_gb=0):
    """Check available storage space"""
    total, used, free = get_disk_usage(BASE_DIR)
    
    logger.info(f"💽 Storage Status:")
    logger.info(f"   • Total: {total:.1f} GB")
    logger.info(f"   • Used: {used:.1f} GB ({used/total*100:.1f}%)")
    logger.info(f"   • Free: {free:.1f} GB ({free/total*100:.1f}%)")
    
    if required_gb > 0:
        logger.info(f"   • Required: {required_gb} GB")
        if free < required_gb + MIN_FREE_SPACE_GB:
            return False, f"Insufficient space: {free:.1f} GB available, {required_gb + MIN_FREE_SPACE_GB} GB needed"
    
    return True, "Storage OK"

def cleanup_previous_step(step_config, logger):
    """Clean up files from previous pipeline step"""
    logger.info(f"🗑️  Cleaning up after step {step_config['step']}...")
    
    cleaned_size = 0
    
    for output_dir in step_config['output_dirs']:
        full_path = BASE_DIR / output_dir
        
        if not full_path.exists():
            continue
            
        # Calculate size before cleanup
        dir_size = get_directory_size(full_path)
        logger.info(f"   • {output_dir}: {dir_size:.2f} GB")
        
        # If preserve_files is ['*'], keep everything
        if step_config['preserve_files'] == ['*']:
            logger.info(f"     ✅ Preserving all files in {output_dir}")
            continue
            
        # Clean up directory while preserving specified files
        if step_config['preserve_files']:
            # Create temp directory for preserved files
            temp_preserve = full_path.parent / f"{full_path.name}_preserve_temp"
            temp_preserve.mkdir(exist_ok=True)
            
            # Move preserved files to temp
            for preserve_pattern in step_config['preserve_files']:
                for file_path in full_path.glob(preserve_pattern):
                    if file_path.is_file():
                        shutil.move(str(file_path), str(temp_preserve / file_path.name))
                    elif file_path.is_dir():
                        shutil.move(str(file_path), str(temp_preserve / file_path.name))
            
            # Remove original directory
            shutil.rmtree(full_path)
            
            # Move temp back to original location
            temp_preserve.rename(full_path)
            
            preserved_size = get_directory_size(full_path)
            cleaned_size += (dir_size - preserved_size)
            logger.info(f"     🗑️  Cleaned: {dir_size - preserved_size:.2f} GB")
            logger.info(f"     ✅ Preserved: {preserved_size:.2f} GB")
        else:
            # Remove entire directory
            shutil.rmtree(full_path)
            cleaned_size += dir_size
            logger.info(f"     🗑️  Removed: {dir_size:.2f} GB")
    
    if cleaned_size > 0:
        logger.info(f"✅ Total cleaned: {cleaned_size:.2f} GB")
    else:
        logger.info("✅ No cleanup needed")

def protect_raw_data(logger):
    """Ensure raw data is protected from cleanup"""
    raw_files = list(RAW_DATA_DIR.glob("**/*"))
    
    logger.info(f"🛡️  Raw data protection check:")
    logger.info(f"   • Raw data directory: {RAW_DATA_DIR}")
    logger.info(f"   • Files protected: {len(raw_files)}")
    
    total_raw_size = get_directory_size(RAW_DATA_DIR)
    logger.info(f"   • Total raw data size: {total_raw_size:.2f} GB")

def check_dependencies(logger):
    """Check if required executables are available"""
    logger.info("🔍 Checking dependencies...")
    
    # Check Python
    try:
        python_version = subprocess.run([PYTHON_EXECUTABLE, '--version'], 
                                      capture_output=True, text=True)
        if python_version.returncode == 0:
            logger.info(f"   ✅ Python: {python_version.stdout.strip()}")
        else:
            logger.error(f"   ❌ Python check failed")
            return False
    except Exception as e:
        logger.error(f"   ❌ Python not found: {e}")
        return False
    
    # Check R (optional - only if R scripts are present)
    r_scripts = [step for step in PIPELINE_STEPS if step.get('script_type') == 'r']
    if r_scripts:
        try:
            r_version = subprocess.run([R_EXECUTABLE, '--version'], 
                                     capture_output=True, text=True)
            if r_version.returncode == 0:
                version_line = r_version.stdout.split('\n')[0]
                logger.info(f"   ✅ R: {version_line}")
            else:
                logger.error(f"   ❌ R check failed")
                return False
        except Exception as e:
            logger.error(f"   ❌ R not found: {e}")
            logger.error(f"   💡 Install R or remove R scripts from pipeline")
            return False
    
    return True

# ...existing code...

def run_script(script_path, step_config, logger):
    """Run a pipeline script with error handling for both Python and R"""
    logger.info(f"🚀 Starting: {step_config['description']}")
    logger.info(f"   • Script: {script_path}")
    logger.info(f"   • Memory requirement: {step_config['required_memory_gb']} GB")
    
    # Determine script type
    script_type = step_config.get('script_type')
    if script_type is None:
        script_type = detect_script_type(script_path)
        logger.info(f"   • Auto-detected type: {script_type}")
    else:
        logger.info(f"   • Script type: {script_type}")
    
    # Prepare execution command
    if script_type == 'python':
        cmd = [PYTHON_EXECUTABLE, script_path.name]
        timeout = PYTHON_TIMEOUT
    elif script_type == 'r':
        cmd = [R_EXECUTABLE, script_path.name]
        timeout = R_TIMEOUT
    else:
        logger.error(f"   ❌ Unsupported script type: {script_type}")
        return False, f"Unsupported script type: {script_type}"
    
    logger.info(f"   • Command: {' '.join(cmd)}")
    logger.info(f"   • Timeout: {timeout/3600:.1f} hours")
    
    start_time = time.time()
    
    try:
        # Change to scripts directory
        original_cwd = os.getcwd()
        os.chdir(SCRIPTS_DIR)
        
        # Run the script
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        
        # Restore working directory
        os.chdir(original_cwd)
        
        duration = time.time() - start_time
        
        if result.returncode == 0:
            logger.info(f"✅ Completed: {step_config['description']} ({duration:.1f}s)")
            
            # Log any important output (last few lines)
            if result.stdout:
                stdout_lines = result.stdout.strip().split('\n')
                if len(stdout_lines) > 0:
                    logger.info(f"   • Output summary:")
                    for line in stdout_lines[-3:]:  # Last 3 lines
                        if line.strip():
                            logger.info(f"     {line.strip()}")
            
            return True, "Success"
        else:
            logger.error(f"❌ Failed: {step_config['description']}")
            logger.error(f"   • Exit code: {result.returncode}")
            logger.error(f"   • Duration: {duration:.1f}s")
            
            # Log stderr output
            if result.stderr:
                logger.error(f"   • Error output:")
                stderr_lines = result.stderr.strip().split('\n')
                for line in stderr_lines[-10:]:  # Last 10 error lines
                    if line.strip():
                        logger.error(f"     {line.strip()}")
            
            # Log stdout as it might contain useful info even on failure
            if result.stdout:
                logger.error(f"   • Standard output:")
                stdout_lines = result.stdout.strip().split('\n')
                for line in stdout_lines[-5:]:  # Last 5 output lines
                    if line.strip():
                        logger.error(f"     {line.strip()}")
            
            return False, result.stderr
            
    except subprocess.TimeoutExpired:
        timeout_hours = timeout / 3600
        logger.error(f"❌ Timeout: {step_config['description']} (>{timeout_hours:.1f} hours)")
        return False, "Timeout"
    except Exception as e:
        logger.error(f"❌ Exception: {step_config['description']}: {e}")
        return False, str(e)
    finally:
        os.chdir(original_cwd)

def validate_pipeline_scripts(logger):
    """Validate that all pipeline scripts exist and have correct types"""
    logger.info("📋 Validating pipeline scripts...")
    
    all_valid = True
    
    for step in PIPELINE_STEPS:
        script_path = SCRIPTS_DIR / step['script']
        
        if not script_path.exists():
            logger.error(f"   ❌ {step['step']}: Script not found - {script_path}")
            all_valid = False
            continue
        
        # Check script type consistency
        detected_type = detect_script_type(script_path)
        declared_type = step.get('script_type', detected_type)
        
        if detected_type != declared_type and detected_type != 'unknown':
            logger.warning(f"   ⚠️  {step['step']}: Type mismatch - declared:{declared_type}, detected:{detected_type}")
        
        logger.info(f"   ✅ {step['step']}: {step['script']} ({declared_type})")
    
    return all_valid

def main():
    """Main pipeline execution"""
    
    # Setup logging
    logger = setup_logging()
    
    logger.info("=" * 70)
    logger.info("🚀 snRNA-seq PROCESSING PIPELINE RUNNER")
    logger.info("   Multi-language support (Python & R)")
    logger.info("=" * 70)
    
    # Validate dependencies and scripts
    if not check_dependencies(logger):
        logger.error("❌ Dependency check failed")
        return False
    
    if not validate_pipeline_scripts(logger):
        logger.error("❌ Script validation failed")
        return False
    
    # Initial setup
    protect_raw_data(logger)
    
    # Check initial storage
    storage_ok, storage_msg = check_storage_space(logger)
    if not storage_ok:
        logger.error(f"❌ {storage_msg}")
        return False
    
    successful_steps = []
    
    try:
        for i, step_config in enumerate(PIPELINE_STEPS):
            logger.info(f"\n{'='*50}")
            logger.info(f"STEP {i+1}/{len(PIPELINE_STEPS)}: {step_config['description']}")
            logger.info(f"{'='*50}")
            
            # Check if script exists
            script_path = SCRIPTS_DIR / step_config['script']
            if not script_path.exists():
                logger.error(f"❌ Script not found: {script_path}")
                break
            
            # Check storage before step
            storage_ok, storage_msg = check_storage_space(
                logger, step_config['required_memory_gb']
            )
            if not storage_ok:
                logger.warning(f"⚠️  {storage_msg}")
                logger.info("🗑️  Attempting emergency cleanup...")
                
                # Emergency cleanup of non-essential files
                for prev_step in successful_steps:
                    cleanup_previous_step(prev_step, logger)
                
                # Recheck storage
                storage_ok, storage_msg = check_storage_space(
                    logger, step_config['required_memory_gb']
                )
                if not storage_ok:
                    logger.error(f"❌ Still insufficient space after cleanup: {storage_msg}")
                    break
            
            # Run the step
            success, error_msg = run_script(script_path, step_config, logger)
            
            if success:
                successful_steps.append(step_config)
                
                # Cleanup previous step (except the last one)
                if i > 0 and i < len(PIPELINE_STEPS) - 1:
                    cleanup_previous_step(PIPELINE_STEPS[i-1], logger)
                
                # Check storage after step
                check_storage_space(logger)
                
            else:
                logger.error(f"❌ Pipeline failed at step {step_config['step']}: {error_msg}")
                logger.error(f"💡 Check the logs above for detailed error information")
                break
        
        # Final summary
        logger.info(f"\n{'='*70}")
        logger.info("📊 PIPELINE EXECUTION SUMMARY")
        logger.info(f"{'='*70}")
        
        logger.info(f"✅ Successful steps: {len(successful_steps)}/{len(PIPELINE_STEPS)}")
        for step in successful_steps:
            script_type = step.get('script_type', 'unknown')
            logger.info(f"   • {step['step']}: {step['description']} ({script_type})")
        
        if len(successful_steps) == len(PIPELINE_STEPS):
            logger.info("🎉 PIPELINE COMPLETED SUCCESSFULLY!")
            
            # Final storage check
            check_storage_space(logger)
            
            return True
        else:
            failed_steps = len(PIPELINE_STEPS) - len(successful_steps)
            logger.warning(f"⚠️  Pipeline incomplete: {failed_steps} steps failed")
            
            # Show next step that would run
            if len(successful_steps) < len(PIPELINE_STEPS):
                next_step = PIPELINE_STEPS[len(successful_steps)]
                logger.info(f"💡 Next step would be: {next_step['step']} - {next_step['description']}")
            
            return False
            
    except KeyboardInterrupt:
        logger.warning("⚠️  Pipeline interrupted by user")
        return False
    except Exception as e:
        logger.error(f"❌ Unexpected error: {e}")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
