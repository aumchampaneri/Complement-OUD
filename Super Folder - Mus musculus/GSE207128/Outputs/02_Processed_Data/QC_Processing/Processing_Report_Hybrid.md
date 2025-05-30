# Single-Cell RNA-seq Processing Report - HYBRID CONFIGURATION
Generated on: 2025-05-25
Processing time: 10.16 minutes

## Configuration Applied
- Configuration type: Amygdala-optimized
- Filtering method: Fixed thresholds (tissue-optimized)
- Tissue optimization: TRUE

## Processing Parameters
- Random seed: 42
- Minimum features: 250
- Maximum features: 6000
- Maximum mitochondrial %: 12
- Variable features: 2000
- PCA dimensions: 25
- Optimal clustering resolution: 0.8
- Silhouette analysis: FALSE

## Processing Features
- ✅ Hybrid configuration system
- ✅ Tissue-specific optimization
- ✅ Smart QC threshold selection
- ✅ Comprehensive visualization
- ✅ Memory-efficient processing
- ✅ Resolution optimization
- ✅ Cell cycle scoring
- ✅ Robust error handling
-  ❌  Doublet detection

## Results Summary
- Initial cells: 121704
- Final cells: 79477
- Cells removed: 42227 ( 34.7 %)
- Final clusters: 35
- Median genes per cell: 2044
- Median UMIs per cell: 4886
- Median MT percent: 6.79 %

## Cell Type Markers Available
- Cd68, Cx3cr1, Tmem119, Gfap, Mbp, Snap25, Slc17a7, Gad1

## Output Files
### Main Results
- processed_seurat.rds - Final processed Seurat object

### Quality Control Plots
- QC_violins_pre_filtering.png - Pre-filtering QC metrics
- QC_violins_post_filtering.png - Post-filtering QC metrics
- QC_histograms_pre_filtering.png - QC histograms with thresholds
- filtering_comparison.png - Pre vs post filtering comparison

### Analysis Plots
- variable_features.png - Highly variable features
- PCA_analysis.png - PCA visualization
- clustering_resolution_tree.png - Resolution optimization tree (if available)
- resolution_optimization.png - Resolution analysis plot
- UMAP_combined_overview.png - UMAP overview
- cell_cycle_phases.png - Cell cycle analysis
- marker_*.png - Individual marker expression plots
- UMAP_final.png - Final processed UMAP

### Statistics and Data
- filtering_summary_by_sample.csv - Filtering statistics per sample
- resolution_analysis.csv - Clustering resolution analysis
- cell_cycle_summary.csv - Cell cycle phase distribution
- final_processing_metrics.rds - Complete processing metrics

### Adaptive Thresholds
- adaptive_feature_thresholds.csv - Sample-specific feature thresholds
- adaptive_count_thresholds.csv - Sample-specific count thresholds
- adaptive_mt_thresholds.csv - Sample-specific mitochondrial thresholds

## Next Steps
1. Run 2_Integration.R for batch correction
2. Run 3_CellType_Annotation.R for cell type identification
3. Proceed with downstream analysis
