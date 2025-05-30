# Single-Cell RNA-seq Integration Report
Generated on: 2025-05-25
Processing time: 15.27 minutes

## Integration Method
- Method used: Batch Regression (FastMNN fallback)
- Reduction: pca
- Dimensions used: 19
- Configuration: comprehensive_robust

## Dataset Information
- Input cells: 79477
- Output cells: 79477
- Variable features: 3000
- Samples: 9
- Conditions: 3

## Clustering Results
- Resolutions tested: 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5
- Default resolution: 0.8
- Final clusters: 34

## Integration Quality
- Sample mixing entropy: 2.888
- Condition mixing entropy: 1.506

## Output Files
### Main Results
-  seurat_integrated.rds  - Integrated Seurat object

### Visualizations
- integration_comparison.png - Before/after integration
- integration_overview.png - Post-integration overview
- UMAP_by_condition.png - Conditions on UMAP
- UMAP_by_clusters.png - Clusters on UMAP
- UMAP_by_sample.png - Samples on UMAP

### Metrics and Data
- integration_summary.rds - Complete integration summary
- integration_metrics.rds - Quality metrics
- sample_mixing_by_cluster.csv - Sample distribution per cluster
- condition_mixing_by_cluster.csv - Condition distribution per cluster

## Next Steps
1. Run 3_CellType_Annotation.R for cell type identification
2. Perform downstream analysis on integrated data
3. Validate integration quality with known markers
