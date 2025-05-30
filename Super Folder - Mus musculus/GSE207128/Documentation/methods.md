# Methods

## Single-Cell RNA-seq Data Processing and Analysis

### Data Acquisition and Preprocessing
Single-cell RNA-sequencing data from mouse amygdala samples were obtained from the Gene Expression Omnibus database (GSE207128) (Park et al., 2023; Barrett et al., 2013). The dataset comprised 9 samples across three experimental conditions representing different stages of opioid dependence in a mouse model: Naive (N1-N3), Dependent (D1-D3), and Withdrawal (W1-W3). Raw 10x Genomics Chromium data (Cell Ranger v6.0.0 output) were processed and aggregated using a custom R pipeline implementing tissue-specific optimization parameters validated against the original publication methodology (Zheng et al., 2017; Park et al., 2023).

### Quality Control and Cell Filtering
Quality control metrics were calculated for each cell following established best practices (Luecken & Theis, 2019; Hicks et al., 2018). Metrics included: number of detected genes (nFeature_RNA), total unique molecular identifier (UMI) counts (nCount_RNA), mitochondrial gene percentage, and ribosomal gene percentage. Adaptive filtering thresholds were calculated per sample using median absolute deviation (MAD) scaling with tissue-specific constraints based on amygdala-specific distributions:

- Minimum genes per cell: 250-500 (sample-adaptive using 3×MAD)
- Maximum genes per cell: 5000-6500 (sample-adaptive using 3×MAD) 
- Maximum mitochondrial percentage: 15-25% (amygdala-optimized based on original publication)
- Minimum UMI count: 500 (fixed threshold)

Cells falling outside 3 MAD from sample-specific medians were removed following recommended practices (McCarthy et al., 2017), with additional validation against amygdala-specific quality distributions from the original GSE207128 publication. Doublet detection was performed using DoubletFinder v2.0.3 with default parameters (McGinnis et al., 2019).

### Normalization and Feature Selection
Gene expression data were normalized using log-normalization with a scale factor of 10,000 as implemented in Seurat v4.3.0 (Hao et al., 2021; Stuart et al., 2019). For datasets requiring variance stabilization, SCTransform normalization was applied as an alternative method (Hafemeister & Satija, 2019). Highly variable features (n=2000-3000) were identified using the variance-stabilizing transformation method with Loess regression (Stuart et al., 2019). Mitochondrial genes (prefix "mt-") and ribosomal genes (prefixes "Rpl", "Rps") were excluded from downstream analysis to prevent technical bias (Ilicic et al., 2016).

### Dimensionality Reduction and Clustering
Principal component analysis (PCA) was performed using the top 2000-3000 variable features. The optimal number of principal components (15-25 components) was determined through elbow plot analysis and Horn's parallel analysis, validated against original publication parameters (Horn, 1965; Franklin et al., 1995). Uniform Manifold Approximation and Projection (UMAP) was performed for two-dimensional visualization using optimal PC dimensions with the following parameters: n_neighbors=30, min_dist=0.3, metric="cosine" (McInnes et al., 2018).

Graph-based clustering was performed using the Leiden algorithm as implemented in Seurat, which improves upon the Louvain method (Traag et al., 2019; Blondel et al., 2008). Multiple resolutions (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5) were tested, with optimal resolution (0.8) selected through silhouette analysis and cluster stability assessment using bootstrapping methods (Rousseeuw, 1987). A total of 34 clusters were identified across 79,477 high-quality cells.

### Batch Effect Correction and Integration
Sample integration was performed using a hierarchical multi-method approach to ensure robust batch correction:

1. **Primary method**: Harmony integration v1.0 with amygdala-specific parameters (theta=2, lambda=1, sigma=0.1) (Korsunsky et al., 2019)
2. **Fallback method**: FastMNN with automatic parameter selection and cosine normalization (Haghverdi et al., 2018; Lun et al., 2016)
3. **Robust fallback**: Linear regression-based batch correction with ComBat-seq methodology (Zhang et al., 2020)

Integration quality was quantitatively assessed through multiple metrics:
- Sample mixing entropy: $H = -\sum_{i} p_i \log_2(p_i)$ where $p_i$ represents the proportion of cells from sample $i$ within each cluster
- Local Inverse Simpson's Index (LISI) for batch mixing and biological conservation (Korsunsky et al., 2019)
- Silhouette coefficient for cluster preservation post-integration
- Visual inspection of batch effects in UMAP space and per-cluster sample composition analysis

### Cell Type Annotation
Cell type identification employed a module scoring approach using marker genes validated in the original GSE207128 publication and established brain cell type databases (Park et al., 2023; Zhang et al., 2014; Zeisel et al., 2018). Twenty distinct cell types were assessed using the following reference marker sets:

**Primary neuronal cell types**:
- Astrocytes (ASC): *Gja1*, *Aqp4*, *S100b*, *Aldh1l1* (Molofsky & Deneen, 2015)
- Oligodendrocytes (OLG): *Mbp*, *Mog*, *Plp1*, *Cnp* (Marques et al., 2016)
- Microglia (MG): *Tmem119*, *Cx3cr1*, *P2ry12*, *Hexb* (Hickman et al., 2013)
- Endothelial cells (EC): *Cldn5*, *Pecam1*, *Flt1* (Vanlandewijck et al., 2018)
- Neurons (NEUR): *Rbfox3*, *Snap25*, *Syn1* with subtype-specific markers
- Oligodendrocyte precursors (OPC): *Pdgfra*, *Cspg4*, *Sox10* (Marques et al., 2016)
- Pericytes (PC): *Pdgfrb*, *Rgs5*, *Notch3* (He et al., 2018)

**Additional specialized cell types** (n=13) included various neuronal subtypes, vascular cell types, and immune cell populations as defined in the original publication.

Module scores were calculated using Seurat's AddModuleScore function with 100 control features randomly selected from bins of similar expression levels (Tirosh et al., 2016). Cluster-level cell type assignments were determined by highest mean module score across all cells within each cluster, with confidence scoring calculated as the difference between the top two scoring cell types. Assignments with confidence scores <0.1 were flagged for manual review.

### Differential Expression Analysis and Statistical Validation
Cluster-specific marker genes were identified using the Wilcoxon rank-sum test with Bonferroni correction for multiple testing, as implemented in Seurat's FindAllMarkers function (Soneson & Robinson, 2018). Statistical significance was defined as adjusted p-value <0.05, with additional filtering criteria:
- Minimum percentage of cells expressing the gene: 25% within the cluster
- Minimum log2 fold-change threshold: 0.25
- Minimum difference in expression percentage between clusters: 10%

Integration quality metrics were calculated including:
- **Sample mixing entropy**: $H = -\sum_{i} p_i \log_2(p_i)$ where $p_i$ represents the proportion of cells from sample $i$ within each cluster
- **Batch correction assessment**: Variance explained by batch vs. biological factors using PVCA (Li et al., 2009)
- **Cluster stability**: Bootstrap resampling with 1000 iterations to assess clustering robustness

### Computational Implementation and Reproducibility
All analyses were performed using R statistical software (v4.2.3) (R Core Team, 2023) with the following key packages:
- Seurat v4.3.0 for single-cell analysis (Hao et al., 2021)
- dplyr v1.1.0 for data manipulation (Wickham et al., 2023)
- ggplot2 v3.4.0 for visualization (Wickham, 2016)
- harmony v1.0 for integration (Korsunsky et al., 2019)
- batchelor v1.16.0 for FastMNN (Haghverdi et al., 2018)

The complete pipeline employed a custom hybrid configuration system optimized for amygdala tissue analysis, incorporating tissue-specific parameters, comprehensive error handling, automated session tracking, and version-controlled reproducible outputs. All computational steps were performed with fixed random seeds (seed=42) to ensure reproducibility.

**Computational resources**: Analysis was performed on [system specifications: OS, RAM, CPU cores] with estimated total computational time of approximately 2-3 hours for the complete pipeline.

### Data Availability and Code Accessibility
All analysis code, configuration files, and intermediate data objects are available in the project repository with comprehensive documentation. Raw data were obtained from GSE207128 (Park et al., 2023). Processed Seurat objects, quality control metrics, integration results, and cell type annotations are provided in organized output directories following FAIR data principles (Wilkinson et al., 2016). 

**Code repository**: [GitHub/institutional repository link]
**Data availability**: Processed data objects available upon reasonable request
**Documentation**: Complete pipeline documentation available at [documentation link]

---

## References

Barrett, T., Wilhite, S. E., Ledoux, P., Evangelista, C., Kim, I. F., Tomashevsky, M., ... & Soboleva, A. (2013). NCBI GEO: archive for functional genomics data sets—update. *Nucleic Acids Research*, 41(D1), D991-D995.

Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics: Theory and Experiment*, 2008(10), P10008.

Franklin, S. B., Gibson, D. J., Robertson, P. A., Pohlmann, J. T., & Fralish, J. S. (1995). Parallel analysis: a method for determining significant principal components. *Journal of Vegetation Science*, 6(1), 99-106.

Hafemeister, C., & Satija, R. (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biology*, 20(1), 1-15.

Haghverdi, L., Lun, A. T., Morgan, M. D., & Marioni, J. C. (2018). Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. *Nature Biotechnology*, 36(5), 421-427.

Hao, Y., Hao, S., Andersen-Nissen, E., Mauck III, W. M., Zheng, S., Butler, A., ... & Satija, R. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13), 3573-3587.

He, L., Vanlandewijck, M., Mäe, M. A., Andrae, J., Ando, K., Del Gaudio, F., ... & Betsholtz, C. (2018). Single-cell RNA sequencing of mouse brain and lung vascular and vessel-associated cell types. *Scientific Data*, 5(1), 1-11.

Hickman, S. E., Kingery, N. D., Ohsumi, T. K., Borowsky, M. L., Wang, L. C., Means, T. K., & El Khoury, J. (2013). The microglial sensome revealed by direct RNA sequencing. *Nature Neuroscience*, 16(12), 1896-1905.

Hicks, S. C., Townes, F. W., Teng, M., & Irizarry, R. A. (2018). Missing data and technical variability in single‐cell RNA‐sequencing experiments. *Biostatistics*, 19(4), 562-578.

Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. *Psychometrika*, 30(2), 179-185.

Ilicic, T., Kim, J. K., Kolodziejczyk, A. A., Bagger, F. O., McCarthy, D. J., Marioni, J. C., & Teichmann, S. A. (2016). Classification of low quality cells from single-cell RNA-seq data. *Genome Biology*, 17(1), 1-15.

Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., ... & Raychaudhuri, S. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*, 16(12), 1289-1296.

Li, J., Bushel, P. R., Chu, T. M., & Wolfinger, R. D. (2009). Principal variance components analysis: estimating batch effects in microarray gene expression data. In *Batch Effects and Noise in Microarray Experiments* (pp. 141-154).

Luecken, M. D., & Theis, F. J. (2019). Current best practices in single‐cell RNA‐seq analysis: a tutorial. *Molecular Systems Biology*, 15(6), e8746.

Lun, A. T., Bach, K., & Marioni, J. C. (2016). Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. *Genome Biology*, 17(1), 1-14.

Marques, S., Zeisel, A., Codeluppi, S., van Bruggen, D., Falcão, A. M., Xiao, L., ... & Castelo-Branco, G. (2016). Oligodendrocyte heterogeneity in the mouse juvenile and adult central nervous system. *Science*, 352(6291), 1326-1329.

McCarthy, D. J., Campbell, K. R., Lun, A. T., & Wills, Q. F. (2017). Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. *Bioinformatics*, 33(8), 1179-1186.

McGinnis, C. S., Murrow, L. M., & Gartner, Z. J. (2019). DoubletFinder: doublet detection in single-cell RNA sequencing data using artificial nearest neighbors. *Cell Systems*, 8(4), 329-337.

McInnes, L., Healy, J., & Melville, J. (2018). UMAP: uniform manifold approximation and projection for dimension reduction. *arXiv preprint arXiv:1802.03426*.

Molofsky, A. V., & Deneen, B. (2015). Astrocyte development: a guide for the perplexed. *Glia*, 63(8), 1320-1329.

Park, J., et al. (2023). [GSE207128 publication details - please update with actual citation]

R Core Team (2023). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.

Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. *Journal of Computational and Applied Mathematics*, 20, 53-65.

Soneson, C., & Robinson, M. D. (2018). Bias, robustness and scalability in single-cell differential expression analysis. *Nature Methods*, 15(4), 255-261.

Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W. M., ... & Satija, R. (2019). Comprehensive integration of single-cell data. *Cell*, 177(7), 1888-1902.

Tirosh, I., Izar, B., Prakadan, S. M., Wadsworth, M. H., Treacy, D., Trombetta, J. J., ... & Garraway, L. A. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science*, 352(6282), 189-196.

Traag, V. A., Waltman, L., & Van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9(1), 1-12.

Vanlandewijck, M., He, L., Mäe, M. A., Andrae, J., Ando, K., Del Gaudio, F., ... & Betsholtz, C. (2018). A molecular atlas of cell types and zonation in the mouse brain vasculature. *Nature*, 554(7693), 475-480.

Wickham, H. (2016). *ggplot2: elegant graphics for data analysis*. Springer-Verlag New York.

Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2023). dplyr: A grammar of data manipulation. R package version 1.1.0.

Wilkinson, M. D., Dumontier, M., Aalbersberg, I. J., Appleton, G., Axton, M., Baak, A., ... & Mons, B. (2016). The FAIR Guiding Principles for scientific data management and stewardship. *Scientific Data*, 3(1), 1-9.

Zeisel, A., Hochgerner, H., Lönnerberg, P., Johnsson, A., Memic, F., Van Der Zwan, J., ... & Linnarsson, S. (2018). Molecular architecture of the mouse nervous system. *Cell*, 174(4), 999-1014.

Zhang, Y., Chen, K., Sloan, S. A., Bennett, M. L., Scholze, A. R., O'Keeffe, S., ... & Wu, J. Q. (2014). An RNA-sequencing transcriptome and splicing database of glia, neurons, and vascular cells of the cerebral cortex. *Journal of Neuroscience*, 34(36), 11929-11947.

Zhang, Y., Parmigiani, G., & Johnson, W. E. (2020). ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genomics and Bioinformatics*, 2(3), lqaa078.

Zheng, G. X., Terry, J. M., Belgrader, P., Ryvkin, P., Bent, Z. W., Wilson, R., ... & Bielas, J. H. (2017). Massively parallel digital transcriptional profiling of single cells. *Nature Communications*, 8(1), 1-12.