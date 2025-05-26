# Single-Cell RNA-seq Processing Report - STREAMLINED
Generated on: 2025-05-24

## Processing Parameters
- Random seed: 42
- Minimum features: 200
- Maximum features: 8000
- Minimum counts: 500
- Maximum counts: 40000
- Maximum mitochondrial %: 15
- Variable features: 2000
- PCA dimensions: 1-30
- Optimal clustering resolution: 0.1

## Processing Features
- ✅ Comprehensive QC analysis
- ✅ Adaptive threshold calculation
- ✅ Variable feature scaling only
- ✅ Resolution optimization
- ✅ Cell cycle scoring
- ❌ Doublet detection (skipped for performance)

## Results Summary
- Initial cells: 121704
- Final cells: 93985
- Cells removed (QC only): 27719 ( 22.78 %)
- Final clusters: 17
- Median genes per cell: 1889
- Median UMIs per cell: 4417

## Output Files
### Main Results
- combined_seurat_streamlined_complete.rds - Final processed Seurat object

### Plots
- QC_violins_pre_filtering.png - Pre-filtering QC metrics
- QC_violins_post_filtering.png - Post-filtering QC metrics
- QC_histograms_pre_filtering.png - QC histograms with thresholds
- filtering_comparison.png - Pre vs post filtering comparison
- variable_features.png - Highly variable features
- PCA_analysis.png - PCA visualization
- clustering_resolution_tree.png - Resolution optimization tree
- resolution_optimization.png - Resolution analysis plot
- UMAP_combined_overview.png - UMAP overview
- cell_cycle_phases.png - Cell cycle analysis
- marker_*.png - Individual marker expression plots
- UMAP_final.png - Final processed UMAP

### Statistics
- filtering_summary_by_sample.csv - Filtering statistics per sample
- resolution_analysis.csv - Clustering resolution analysis
- cell_cycle_summary.csv - Cell cycle phase distribution
- final_processing_metrics.rds - Complete processing metrics
