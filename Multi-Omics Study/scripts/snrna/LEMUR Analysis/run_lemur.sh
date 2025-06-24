#!/bin/bash

#' @title LEMUR Analysis Execution Script
#' @description Shell script to run LEMUR analysis pipeline
#' @author Generated Analysis Pipeline
#' @date 2024

# Exit on any error
set -e

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Print header
print_header() {
    echo -e "${BLUE}"
    echo "============================================================"
    echo "🧬 LEMUR Analysis Pipeline for GSE225158"
    echo "   Single-cell RNA-seq: OUD vs Control"
    echo "============================================================"
    echo -e "${NC}"
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites..."

    # Check if R is installed
    if ! command -v R &> /dev/null; then
        error "R is not installed or not in PATH"
        exit 1
    fi

    # Check R version
    R_VERSION=$(R --version | head -n1 | sed 's/.*R version \([0-9]\+\.[0-9]\+\).*/\1/')
    log "Found R version: $R_VERSION"

    # Check if input file exists
    INPUT_FILE="$PROJECT_DIR/data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad"
    if [[ ! -f "$INPUT_FILE" ]]; then
        error "Input H5AD file not found: $INPUT_FILE"
        exit 1
    fi

    success "Prerequisites check passed"
}

# Install R packages if needed
install_packages() {
    log "Checking and installing R packages..."

    cat > /tmp/install_packages.R << 'EOF'
# Install required packages
required_packages <- c(
    "SingleCellExperiment", "scran", "scater", "lemur",
    "tidyverse", "ggplot2", "viridis", "patchwork",
    "zellkonverter", "HDF5Array", "uwot", "harmony",
    "BiocNeighbors", "Matrix", "limma", "ComplexHeatmap",
    "ggrepel", "circlize"
)

install_if_missing <- function(pkg) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat("Installing:", pkg, "\n")
        if (pkg %in% rownames(available.packages())) {
            install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
        } else {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager", repos = "https://cloud.r-project.org/")
            }
            BiocManager::install(pkg, dependencies = TRUE)
        }
    } else {
        cat("Already installed:", pkg, "\n")
    }
}

cat("Installing required R packages...\n")
invisible(lapply(required_packages, install_if_missing))
cat("Package installation complete!\n")
EOF

    R --vanilla < /tmp/install_packages.R
    rm /tmp/install_packages.R

    success "R packages installation complete"
}

# Run the main LEMUR analysis
run_analysis() {
    log "Starting LEMUR analysis..."

    # Change to script directory
    cd "$SCRIPT_DIR"

    # Run the R script
    if Rscript run_LEMUR_GSE225158.R; then
        success "LEMUR analysis completed successfully!"
    else
        error "LEMUR analysis failed"
        exit 1
    fi
}

# Generate summary report
generate_summary() {
    log "Generating analysis summary..."

    OUTPUT_DIR="$PROJECT_DIR/Multi-Omics Study/results/snrna/lemur_analysis"

    if [[ -d "$OUTPUT_DIR" ]]; then
        echo -e "${GREEN}"
        echo "============================================================"
        echo "📊 ANALYSIS COMPLETE!"
        echo "============================================================"
        echo -e "${NC}"

        echo "📂 Output directory: $OUTPUT_DIR"
        echo ""
        echo "📋 Generated files:"

        # List key output files
        if [[ -f "$OUTPUT_DIR/tables/lemur_de_results.csv" ]]; then
            echo "  ✅ DE results: tables/lemur_de_results.csv"
        fi

        if [[ -f "$OUTPUT_DIR/tables/top_de_genes.csv" ]]; then
            echo "  ✅ Top DE genes: tables/top_de_genes.csv"
        fi

        if [[ -f "$OUTPUT_DIR/plots/umap_overview.png" ]]; then
            echo "  ✅ UMAP overview: plots/umap_overview.png"
        fi

        if [[ -f "$OUTPUT_DIR/plots/volcano_plot.png" ]]; then
            echo "  ✅ Volcano plot: plots/volcano_plot.png"
        fi

        if [[ -f "$OUTPUT_DIR/reports/analysis_report.txt" ]]; then
            echo "  ✅ Analysis report: reports/analysis_report.txt"
            echo ""
            echo "📋 Quick summary:"
            tail -20 "$OUTPUT_DIR/reports/analysis_report.txt"
        fi

        echo ""
        echo -e "${BLUE}To view results:${NC}"
        echo "  cd $OUTPUT_DIR"
        echo "  open plots/umap_overview.png"
        echo "  open plots/volcano_plot.png"
        echo "  cat reports/analysis_report.txt"

    else
        error "Output directory not found. Analysis may have failed."
        exit 1
    fi
}

# Cleanup function
cleanup() {
    log "Cleaning up temporary files..."
    # Add any cleanup commands here
}

# Main execution function
main() {
    print_header

    # Parse command line arguments
    INSTALL_PACKAGES=false
    SKIP_ANALYSIS=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            --install-packages)
                INSTALL_PACKAGES=true
                shift
                ;;
            --skip-analysis)
                SKIP_ANALYSIS=true
                shift
                ;;
            -h|--help)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Options:"
                echo "  --install-packages    Install R packages before running analysis"
                echo "  --skip-analysis      Skip the main analysis (useful for testing)"
                echo "  -h, --help           Show this help message"
                echo ""
                echo "Examples:"
                echo "  cd 'Multi-Omics Study/scripts/snrna/LEMUR Analysis'"
                echo "  ./run_lemur.sh                                    # Run complete analysis"
                echo "  ./run_lemur.sh --install-packages                # Install packages and run analysis"
                echo "  ./run_lemur.sh --install-packages --skip-analysis # Only install packages"
                exit 0
                ;;
            *)
                error "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done

    # Execute pipeline steps
    check_prerequisites

    if [[ "$INSTALL_PACKAGES" == true ]]; then
        install_packages
    fi

    if [[ "$SKIP_ANALYSIS" == false ]]; then
        run_analysis
        generate_summary
    fi

    cleanup

    success "Pipeline execution complete!"
}

# Set up error handling
trap 'error "Script failed on line $LINENO"' ERR
trap cleanup EXIT

# Run main function with all arguments
main "$@"
