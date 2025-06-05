'''
🔬 Differential Expression Analysis - OUD vs Control (edgeR)
GSE225158 - Pseudo-bulk DE using edgeR for statistical rigor

Strategy:
1. Load scVI-annotated data  
2. Aggregate raw counts by donor/condition/cell_type (pseudo-bulk)
3. Run edgeR GLM with proper dispersion estimation
4. Test multiple contrasts with FDR correction
5. Generate comprehensive results

Dependencies:
- R with edgeR package installed
- rpy2 for Python-R integration
'''

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns  # Removed - not used in this script
import os
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

# R integration for edgeR
try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.packages import importr
    R_AVAILABLE = True
    print("✅ R integration available")
except ImportError:
    R_AVAILABLE = False
    print("❌ R integration not available - install rpy2 and R with edgeR")
    exit(1)

# Paths
BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
INPUT_H5AD = f"{BASE_DIR}/data/processed/snrna_scvi/GSE225158_annotated_scvi.h5ad"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/differential_expression_edgeR"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"

# Create output directories
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

def load_and_prepare_data():
    """Load annotated data and prepare for pseudo-bulk analysis"""
    print("📂 LOADING AND PREPARING DATA")
    print("=" * 50)
    
    if not os.path.exists(INPUT_H5AD):
        raise FileNotFoundError(f"Annotated data not found: {INPUT_H5AD}")
    
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"   Loaded: {adata.shape}")
    
    # Check required columns
    required_cols = ['cell_type', 'Dx_OUD']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Look for brain region and sex columns
    region_col = None
    sex_col = None
    
    for col in ['Region', 'region', 'brain_region', 'area']:
        if col in adata.obs.columns:
            region_col = col
            print(f"   ✅ Found region column: {col}")
            break
    
    for col in ['Sex', 'sex', 'gender', 'Gender']:
        if col in adata.obs.columns:
            sex_col = col
            print(f"   ✅ Found sex column: {col}")
            break
    
    # Create pseudo-bulk grouping variable
    if region_col and sex_col:
        adata.obs['pseudobulk_group'] = (
            adata.obs[region_col].astype(str) + '_' + 
            adata.obs[sex_col].astype(str)
        )
        print(f"   ✅ Using region + sex grouping: {adata.obs['pseudobulk_group'].nunique()} groups")
    else:
        raise ValueError("Both region and sex columns required for edgeR analysis")
    
    return adata

def create_pseudobulk(adata):
    """Aggregate single cells to pseudo-bulk samples"""
    print("\n🧮 CREATING PSEUDO-BULK AGGREGATES")
    print("=" * 50)
    
    # Use raw counts for DE analysis
    if adata.raw is not None:
        count_matrix = adata.raw.X
        gene_names = adata.raw.var_names
        print("   Using raw counts from adata.raw")
    else:
        count_matrix = adata.X
        gene_names = adata.var_names
        print("   Using counts from adata.X")
    
    # Convert to dense if sparse
    if sparse.issparse(count_matrix):
        count_matrix = count_matrix.toarray()
    
    # Create aggregation metadata
    agg_metadata = adata.obs[['pseudobulk_group', 'Dx_OUD']].copy()
    agg_metadata['sample_id'] = (
        agg_metadata['pseudobulk_group'].astype(str) + '_' +
        agg_metadata['Dx_OUD'].astype(str)
    )
    
    # Aggregate counts by sample_id
    unique_samples = agg_metadata['sample_id'].unique()
    pseudobulk_counts = np.zeros((len(gene_names), len(unique_samples)))
    sample_metadata = []
    
    print(f"   Aggregating {len(adata)} cells into {len(unique_samples)} pseudo-bulk samples...")
    
    for i, sample_id in enumerate(unique_samples):
        sample_mask = agg_metadata['sample_id'] == sample_id
        sample_cells = np.where(sample_mask)[0]
        
        # Sum counts across all cells in this sample
        pseudobulk_counts[:, i] = np.sum(count_matrix[sample_cells, :], axis=0)
        
        # Store metadata
        sample_info = agg_metadata.loc[sample_mask].iloc[0]
        sample_metadata.append({
            'sample_id': sample_id,
            'pseudobulk_group': sample_info['pseudobulk_group'],
            'condition': sample_info['Dx_OUD'],
            'n_cells': sample_mask.sum()
        })
    
    # Create DataFrame
    sample_ids = [s['sample_id'] for s in sample_metadata]
    pseudobulk_df = pd.DataFrame(
        pseudobulk_counts.T,
        columns=gene_names,
        index=sample_ids
    )
    
    metadata_df = pd.DataFrame(sample_metadata)
    metadata_df.set_index('sample_id', inplace=True)
    
    print(f"   Created pseudo-bulk matrix: {pseudobulk_df.shape}")
    print(f"   Samples per condition: {metadata_df['condition'].value_counts().to_dict()}")
    
    return pseudobulk_df, metadata_df

def run_edger_analysis(counts_df, metadata_df):
    """Run comprehensive edgeR differential expression analysis with subgroup contrasts"""
    print("\n🔬 RUNNING EDGER ANALYSIS")
    print("=" * 50)
    print("   🧠 Testing same contrasts as Python implementation:")
    print("     • OUD vs Control in Putamen")
    print("     • OUD vs Control in Caudate") 
    print("     • OUD vs Control in Males")
    print("     • OUD vs Control in Females")
    print("     • How OUD effects differ between Males vs Females")
    print("     • How OUD effects differ between Putamen vs Caudate")
    print("     • Pooled OUD vs Control (all samples)")
    
    if not R_AVAILABLE:
        raise RuntimeError("R integration not available")
    
    # Check sample sizes
    if len(counts_df) < 4:
        raise ValueError(f"Too few samples ({len(counts_df)})")
    
    condition_counts = metadata_df['condition'].value_counts()
    if any(condition_counts < 2):
        raise ValueError(f"Too few replicates per condition: {condition_counts.to_dict()}")
    
    # Extract region and sex from pseudobulk_group for subgroup analysis
    region_sex_split = metadata_df['pseudobulk_group'].str.split('_', expand=True)
    metadata_df['region'] = region_sex_split[0]
    metadata_df['sex'] = region_sex_split[1]
    
    results = {}
    
    # Use context-based conversion for safer R integration
    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        # Import R packages
        try:
            edger = importr('edgeR')
            base = importr('base')
            print("   ✅ edgeR package loaded successfully")
        except Exception as e:
            raise RuntimeError(f"Could not import edgeR: {e}")
        
        # Prepare base data for R
        counts_clean = counts_df.select_dtypes(include=[np.number]).fillna(0).astype(int)
        ro.r('library(edgeR)')
        ro.globalenv['gene_names'] = ro.StrVector(counts_clean.columns.astype(str))
        
        # Define subgroup contrasts to match Python implementation exactly
        subgroup_contrasts = {
            'OUD_vs_Control_Putamen': {'region': 'Putamen'},
            'OUD_vs_Control_Caudate': {'region': 'Caudate'}, 
            'OUD_vs_Control_Male': {'sex': 'M'},
            'OUD_vs_Control_Female': {'sex': 'F'},
            'OUD_Effect_Male_vs_Female': 'interaction_sex',
            'OUD_Effect_Putamen_vs_Caudate': 'interaction_region',
            'Pooled_OUD_vs_Control': {}  # All samples
        }
        
        for contrast_name, filter_criteria in subgroup_contrasts.items():
            print(f"\n   🔬 Running edgeR for {contrast_name}...")
            
            # Handle interaction contrasts differently
            if contrast_name.startswith('OUD_Effect_'):
                if filter_criteria == 'interaction_sex':
                    # Test if OUD effect differs between males and females using edgeR
                    results[contrast_name] = run_edger_interaction_sex(counts_clean, metadata_df, contrast_name)
                elif filter_criteria == 'interaction_region':
                    # Test if OUD effect differs between regions using edgeR
                    results[contrast_name] = run_edger_interaction_region(counts_clean, metadata_df, contrast_name)
                continue
            
            # Filter samples for this subgroup
            if filter_criteria:  # Subgroup analysis
                mask = pd.Series([True] * len(metadata_df), index=metadata_df.index)
                for col, val in filter_criteria.items():
                    mask = mask & (metadata_df[col] == val)
                
                if mask.sum() < 4:
                    print(f"     ⚠️  Too few samples in {contrast_name} ({mask.sum()})")
                    continue
                
                subgroup_metadata = metadata_df[mask]
                subgroup_counts = counts_clean.loc[mask]
                
                # Check if both conditions present
                subgroup_conditions = subgroup_metadata['condition'].unique()
                if len(subgroup_conditions) < 2:
                    print(f"     ⚠️  Only one condition in {contrast_name}")
                    continue
                
            else:  # Pooled analysis (all samples)
                subgroup_metadata = metadata_df
                subgroup_counts = counts_clean
            
            try:
                # Convert subgroup data to R
                r_counts = ro.conversion.py2rpy(subgroup_counts.T.values)
                r_condition = ro.StrVector(subgroup_metadata['condition'].astype(str).values)
                r_sample_names = ro.StrVector(subgroup_counts.index.astype(str))
                
                # Set up R environment for this contrast
                ro.globalenv['counts_matrix'] = r_counts
                ro.globalenv['condition'] = r_condition
                ro.globalenv['sample_names'] = r_sample_names
                
                # Run edgeR pipeline for this subgroup
                ro.r('''
                # Set up count matrix
                rownames(counts_matrix) <- gene_names
                colnames(counts_matrix) <- sample_names
                
                # Simple design matrix: ~ condition only for subgroup analysis
                design <- model.matrix(~ condition)
                print(paste("Design matrix dimensions:", paste(dim(design), collapse="x")))
                
                # Create DGEList
                dge <- DGEList(counts = counts_matrix)
                print(paste("DGEList created:", nrow(dge), "genes x", ncol(dge), "samples"))
                
                # Filter low-expressed genes (more lenient for subgroups)
                min_count <- ifelse(ncol(dge) < 6, 5, 10)  # Lower threshold for small subgroups
                keep <- filterByExpr(dge, design, min.count = min_count, min.total.count = 10)
                dge <- dge[keep, , keep.lib.sizes = FALSE]
                print(paste("After filtering:", nrow(dge), "genes retained"))
                
                # Normalization
                dge <- calcNormFactors(dge, method = "TMM")
                
                # Estimate dispersions (with appropriate settings for small samples)
                if (ncol(dge) < 6) {
                    # For very small samples, use more conservative dispersion estimation
                    dge <- estimateDisp(dge, design, robust = FALSE)
                } else {
                    dge <- estimateDisp(dge, design, robust = TRUE)
                }
                
                print(paste("Common dispersion:", round(dge$common.dispersion, 4)))
                
                # Fit GLM
                fit <- glmQLFit(dge, design, robust = TRUE)
                
                # Test OUD vs Control (coefficient 2 = conditionOUD)
                qlf <- glmQLFTest(fit, coef = 2)
                current_result <- topTags(qlf, n = Inf, sort.by = "PValue")$table
                
                print(paste("edgeR completed for", length(sample_names), "samples"))
                ''')
                
                # Extract results for this contrast
                results[contrast_name] = extract_edger_results(contrast_name)
                
            except Exception as e:
                print(f"     ❌ Failed edgeR for {contrast_name}: {e}")
                continue
        
        print(f"   ✅ edgeR completed: {len(results)} contrasts analyzed")
        return results

def extract_edger_results(contrast_name):
    """Extract edgeR results from R to Python DataFrame"""
    ro.r('''
    result_genes <- rownames(current_result)
    result_logFC <- current_result$logFC
    result_logCPM <- current_result$logCPM
    result_PValue <- current_result$PValue
    result_FDR <- current_result$FDR
    ''')
    
    # Convert to Python
    genes = list(ro.r('result_genes'))
    logFC = list(ro.r('result_logFC'))
    logCPM = list(ro.r('result_logCPM'))
    pvals = list(ro.r('result_PValue'))
    fdr = list(ro.r('result_FDR'))
    
    results_df = pd.DataFrame({
        'gene': genes,
        'logFC': logFC,
        'logCPM': logCPM,
        'PValue': pvals,
        'FDR': fdr,
        'contrast': contrast_name
    })
    
    print(f"     ✅ {contrast_name}: {len(results_df)} genes")
    return results_df


def run_edger_interaction_sex(counts_clean, metadata_df, contrast_name):
    """Run edgeR interaction analysis to test if OUD effect differs between males and females"""
    print(f"     🔬 Running interaction analysis: {contrast_name}")
    
    try:
        # Convert data to R with interaction design
        r_counts = ro.conversion.py2rpy(counts_clean.T.values)
        r_condition = ro.StrVector(metadata_df['condition'].astype(str).values)
        r_sex = ro.StrVector(metadata_df['sex'].astype(str).values)
        r_sample_names = ro.StrVector(counts_clean.index.astype(str))
        
        # Set up R environment
        ro.globalenv['counts_matrix'] = r_counts
        ro.globalenv['condition'] = r_condition
        ro.globalenv['sex'] = r_sex
        ro.globalenv['sample_names'] = r_sample_names
        
        # Run edgeR interaction analysis
        ro.r('''
        # Set up count matrix
        rownames(counts_matrix) <- gene_names
        colnames(counts_matrix) <- sample_names
        
        # Interaction design matrix: ~ condition + sex + condition:sex
        design <- model.matrix(~ condition * sex)
        print(paste("Interaction design matrix dimensions:", paste(dim(design), collapse="x")))
        
        # Create DGEList
        dge <- DGEList(counts = counts_matrix)
        
        # Filter low-expressed genes
        keep <- filterByExpr(dge, design, min.count = 10)
        dge <- dge[keep, , keep.lib.sizes = FALSE]
        print(paste("After filtering:", nrow(dge), "genes retained"))
        
        # Normalization and dispersion estimation
        dge <- calcNormFactors(dge, method = "TMM")
        dge <- estimateDisp(dge, design, robust = TRUE)
        
        # Fit GLM
        fit <- glmQLFit(dge, design, robust = TRUE)
        
        # Test interaction term (condition:sex)
        # This tests if OUD effect differs between males and females
        qlf <- glmQLFTest(fit, coef = ncol(design))  # Last coefficient is interaction
        current_result <- topTags(qlf, n = Inf, sort.by = "PValue")$table
        
        print(paste("Interaction analysis completed"))
        ''')
        
        # Extract results
        return extract_edger_results(contrast_name)
        
    except Exception as e:
        print(f"     ❌ Failed interaction analysis for {contrast_name}: {e}")
        return pd.DataFrame()


def run_edger_interaction_region(counts_clean, metadata_df, contrast_name):
    """Run edgeR interaction analysis to test if OUD effect differs between regions"""
    print(f"     🔬 Running interaction analysis: {contrast_name}")
    
    try:
        # Convert data to R with interaction design
        r_counts = ro.conversion.py2rpy(counts_clean.T.values)
        r_condition = ro.StrVector(metadata_df['condition'].astype(str).values)
        r_region = ro.StrVector(metadata_df['region'].astype(str).values)
        r_sample_names = ro.StrVector(counts_clean.index.astype(str))
        
        # Set up R environment
        ro.globalenv['counts_matrix'] = r_counts
        ro.globalenv['condition'] = r_condition
        ro.globalenv['region'] = r_region
        ro.globalenv['sample_names'] = r_sample_names
        
        # Run edgeR interaction analysis
        ro.r('''
        # Set up count matrix
        rownames(counts_matrix) <- gene_names
        colnames(counts_matrix) <- sample_names
        
        # Interaction design matrix: ~ condition + region + condition:region
        design <- model.matrix(~ condition * region)
        print(paste("Interaction design matrix dimensions:", paste(dim(design), collapse="x")))
        
        # Create DGEList
        dge <- DGEList(counts = counts_matrix)
        
        # Filter low-expressed genes
        keep <- filterByExpr(dge, design, min.count = 10)
        dge <- dge[keep, , keep.lib.sizes = FALSE]
        print(paste("After filtering:", nrow(dge), "genes retained"))
        
        # Normalization and dispersion estimation
        dge <- calcNormFactors(dge, method = "TMM")
        dge <- estimateDisp(dge, design, robust = TRUE)
        
        # Fit GLM
        fit <- glmQLFit(dge, design, robust = TRUE)
        
        # Test interaction term (condition:region)
        # This tests if OUD effect differs between regions
        qlf <- glmQLFTest(fit, coef = ncol(design))  # Last coefficient is interaction
        current_result <- topTags(qlf, n = Inf, sort.by = "PValue")$table
        
        print(paste("Interaction analysis completed"))
        ''')
        
        # Extract results
        return extract_edger_results(contrast_name)
        
    except Exception as e:
        print(f"     ❌ Failed interaction analysis for {contrast_name}: {e}")
        return pd.DataFrame()


def save_edger_results(de_results_dict, summary_df, pseudobulk_df, metadata_df):
    """Save edgeR differential expression results"""
    print("\n💾 SAVING EDGER RESULTS")
    print("=" * 50)
    
    try:
        # Save individual contrast results
        for contrast_name, results_df in de_results_dict.items():
            if not results_df.empty:
                # Add significance column
                results_df['significant'] = (results_df['FDR'] < 0.05) & (abs(results_df['logFC']) > 0.5)
                
                # Save full results
                output_file = f"{OUTPUT_DIR}/edgeR_{contrast_name}_results.csv"
                results_df.to_csv(output_file, index=False)
                print(f"   📁 {contrast_name}: {output_file}")
                
                # Save significant genes only
                sig_genes = results_df[results_df['significant']]
                if len(sig_genes) > 0:
                    sig_output_file = f"{OUTPUT_DIR}/edgeR_{contrast_name}_significant.csv"
                    sig_genes.to_csv(sig_output_file, index=False)
                    print(f"   📁 {contrast_name} (significant): {sig_output_file}")
        
        # Save summary
        summary_output_file = f"{OUTPUT_DIR}/edgeR_summary.csv"
        summary_df.to_csv(summary_output_file, index=False)
        print(f"   📁 Summary: {summary_output_file}")
        
        # Save combined results
        if de_results_dict:
            combined_df = pd.concat(de_results_dict.values(), ignore_index=True)
            combined_output_file = f"{OUTPUT_DIR}/edgeR_all_contrasts.csv"
            combined_df.to_csv(combined_output_file, index=False)
            print(f"   📁 Combined: {combined_output_file}")
        
        # Save pseudobulk data and metadata for reference
        pseudobulk_output_file = f"{OUTPUT_DIR}/pseudobulk_counts.csv"
        pseudobulk_df.to_csv(pseudobulk_output_file)
        print(f"   📁 Pseudobulk: {pseudobulk_output_file}")
        
        metadata_output_file = f"{OUTPUT_DIR}/pseudobulk_metadata.csv"
        metadata_df.to_csv(metadata_output_file)
        print(f"   📁 Metadata: {metadata_output_file}")
        
        print("   ✅ All results saved successfully!")
        
    except Exception as e:
        print(f"   ❌ Error saving results: {e}")
        raise

def analyze_edger_results(de_results_dict):
    """Analyze edgeR results with proper FDR thresholds - matching Python analysis"""
    print("\n📊 ANALYZING EDGER RESULTS")
    print("=" * 50)
    
    summary_stats = []
    
    for contrast_name, de_results in de_results_dict.items():
        print(f"\n   {contrast_name}:")
        print(f"     • Total genes tested: {len(de_results)}")
        print(f"     • LogFC range: {de_results['logFC'].min():.2f} to {de_results['logFC'].max():.2f}")
        print(f"     • P-value range: {de_results['PValue'].min():.2e} to {de_results['PValue'].max():.2e}")
        print(f"     • FDR range: {de_results['FDR'].min():.2e} to {de_results['FDR'].max():.2e}")
        
        # Define significance thresholds (matching Python implementation)
        de_results['significant_fdr05'] = (de_results['FDR'] < 0.05) & (abs(de_results['logFC']) > 1.0)
        de_results['significant_fdr10'] = (de_results['FDR'] < 0.10) & (abs(de_results['logFC']) > 0.5)
        de_results['significant_nominal'] = (de_results['PValue'] < 0.01) & (abs(de_results['logFC']) > 0.5)
        
        # Use nominal p-values for discovery (matching Python approach)
        de_results['significant_strict'] = de_results['significant_fdr05']
        de_results['significant'] = de_results['significant_nominal']  # Match Python
        de_results['direction'] = np.where(de_results['logFC'] > 0, 'Up', 'Down')
        
        # Count significant genes  
        n_sig_fdr05 = de_results['significant_fdr05'].sum()
        n_sig_fdr10 = de_results['significant_fdr10'].sum()
        n_sig_nominal = de_results['significant_nominal'].sum()
        n_up = ((de_results['significant']) & (de_results['logFC'] > 0)).sum()
        n_down = ((de_results['significant']) & (de_results['logFC'] < 0)).sum()
        
        print(f"     • FDR < 0.05 & |FC| > 1.0: {n_sig_fdr05} genes")
        print(f"     • FDR < 0.10 & |FC| > 0.5: {n_sig_fdr10} genes")
        print(f"     • Significant (nominal p<0.01, |FC|>0.5): {n_sig_nominal} genes ({n_up}↑, {n_down}↓)")
        
        # Show top genes by p-value (matching Python)
        if len(de_results) > 0:
            top_genes = de_results.nsmallest(5, 'PValue')[['gene', 'logFC', 'PValue', 'FDR']]
            print(f"     • Top 5 genes by p-value:")
            for _, row in top_genes.iterrows():
                print(f"       - {row['gene']}: logFC={row['logFC']:.2f}, p={row['PValue']:.2e}, FDR={row['FDR']:.2e}")
        
        summary_stats.append({
            'contrast': contrast_name,
            'total_genes': len(de_results),
            'significant_genes': n_sig_nominal,  # Match Python
            'significant_strict': n_sig_fdr05,
            'significant_fdr05': n_sig_fdr05,
            'significant_fdr10': n_sig_fdr10,
            'upregulated': n_up,
            'downregulated': n_down,
            'pct_significant': n_sig_nominal / len(de_results) * 100
        })
    
    return pd.DataFrame(summary_stats)

def create_edger_plots(de_results_dict, summary_df):
    """Create visualization plots matching Python implementation"""
    print("\n📈 CREATING EDGER VISUALIZATION PLOTS")
    print("=" * 50)
    
    if de_results_dict is None or summary_df is None:
        print("   No results to plot")
        return
    
    # Set matplotlib backend
    plt.switch_backend('Agg')
    
    # Create plots for each contrast (matching Python style)
    for contrast_name, de_results in de_results_dict.items():
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # 1. Overall DE summary
        n_sig = de_results['significant'].sum()
        n_up = ((de_results['significant']) & (de_results['logFC'] > 0)).sum()
        n_down = ((de_results['significant']) & (de_results['logFC'] < 0)).sum()
        n_nonsig = len(de_results) - n_sig
        
        categories = ['Upregulated', 'Downregulated', 'Non-significant']
        counts = [n_up, n_down, n_nonsig]
        colors = ['red', 'blue', 'gray']
        
        axes[0,0].bar(categories, counts, color=colors)
        axes[0,0].set_title(f'edgeR DE Gene Summary - {contrast_name}')
        axes[0,0].set_ylabel('Number of Genes')
        
        # 2. Volcano plot
        scatter_colors = de_results['significant'].map({True: 'red', False: 'gray'})
        axes[0,1].scatter(
            de_results['logFC'], 
            -np.log10(de_results['PValue']),
            c=scatter_colors,
            alpha=0.6, s=20
        )
        axes[0,1].set_xlabel('Log2 Fold Change')
        axes[0,1].set_ylabel('-log10(P-value)')
        axes[0,1].set_title(f'edgeR Volcano Plot - {contrast_name}')
        axes[0,1].axhline(y=-np.log10(0.01), color='black', linestyle='--', alpha=0.5)
        axes[0,1].axvline(x=0.5, color='black', linestyle='--', alpha=0.5)
        axes[0,1].axvline(x=-0.5, color='black', linestyle='--', alpha=0.5)
        
        # 3. Top significant genes
        if n_sig > 0:
            top_genes = de_results[de_results['significant']].nsmallest(min(20, n_sig), 'PValue')
            
            y_pos = np.arange(len(top_genes))
            colors = ['red' if fc > 0 else 'blue' for fc in top_genes['logFC']]
            
            axes[1,0].barh(y_pos, top_genes['logFC'], color=colors)
            axes[1,0].set_yticks(y_pos)
            axes[1,0].set_yticklabels(top_genes['gene'], fontsize=8)
            axes[1,0].set_xlabel('Log2 Fold Change')
            axes[1,0].set_title(f'Top {len(top_genes)} Significant DE Genes')
            axes[1,0].invert_yaxis()
        else:
            axes[1,0].text(0.5, 0.5, 'No significant genes found', 
                          transform=axes[1,0].transAxes, ha='center', va='center')
            axes[1,0].set_title('Top Significant DE Genes')
        
        # 4. P-value distribution
        axes[1,1].hist(de_results['PValue'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        axes[1,1].set_xlabel('P-value')
        axes[1,1].set_ylabel('Frequency')
        axes[1,1].set_title('P-value Distribution')
        axes[1,1].axvline(x=0.01, color='red', linestyle='--', alpha=0.7)
        
        plt.suptitle(f'edgeR Differential Expression - {contrast_name}', fontsize=16)
        plt.tight_layout()
        
        # Save plot
        safe_name = contrast_name.replace(' ', '_').replace('vs', 'vs')
        plt.savefig(f"{PLOTS_DIR}/edgeR_analysis_{safe_name}.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"   edgeR plots saved to: {PLOTS_DIR}/")

def main():
    """Main edgeR workflow"""
    print("🔬 EDGER DIFFERENTIAL EXPRESSION ANALYSIS")
    print("=" * 70)
    print("Using pseudo-bulk aggregation + edgeR GLM (matching Python contrasts)")
    print("=" * 70)
    
    try:
        # Load data
        adata = load_and_prepare_data()
        
        # Create pseudo-bulk
        pseudobulk_df, metadata_df = create_pseudobulk(adata)
        
        # Run edgeR analysis (same contrasts as Python)
        de_results_dict = run_edger_analysis(pseudobulk_df, metadata_df)
        
        # Analyze results (matching Python thresholds)
        summary_df = analyze_edger_results(de_results_dict)
        
        # Create plots
        create_edger_plots(de_results_dict, summary_df)
        
        # Save results
        save_edger_results(de_results_dict, summary_df, pseudobulk_df, metadata_df)
        
        print("\n" + "=" * 70)
        print("✅ EDGER ANALYSIS COMPLETE!")
        if de_results_dict is not None:
            total_sig = sum(df.get('significant', pd.Series([])).sum() for df in de_results_dict.values())
            print(f"📊 Found {total_sig} significant DE genes across {len(de_results_dict)} contrasts")
        print(f"📁 Results saved to: {OUTPUT_DIR}/")
        print("🚀 Ready for comparison with Python results!")
        print("=" * 70)
        
        return de_results_dict, summary_df
        
    except Exception as e:
        print(f"❌ edgeR analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None

if __name__ == "__main__":
    main()
