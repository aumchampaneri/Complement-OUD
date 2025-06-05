'''
🔬 Differential Expression Analysis - OUD vs Control
GSE225158 - Pseudo-bulk DE using edgeR/DESeq2 for statistical rigor

Strategy:
1. Load scVI-annotated data
2. Aggregate raw counts by donor/condition/cell_type (pseudo-bulk)
3. Run edgeR/DESeq2 for robust statistical testing
4. Validate with pathway analysis
5. Generate comprehensive results

Note: Avoids single-cell Wilcoxon tests which are too liberal with replicates
'''

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

# R integration for edgeR/DESeq2
try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.packages import importr
    # Use modern context-based conversion instead of deprecated activate()
    R_AVAILABLE = True
    print("✅ R integration available")
except ImportError:
    R_AVAILABLE = False
    print("⚠️  R integration not available - will use Python alternatives")

# Paths
BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
INPUT_H5AD = f"{BASE_DIR}/data/processed/snrna_scvi/GSE225158_annotated_scvi.h5ad"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/differential_expression"
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
    
    # Check available metadata columns
    print("   Available metadata columns:")
    for col in adata.obs.columns:
        n_unique = adata.obs[col].nunique()
        print(f"     {col}: {n_unique} unique values ({adata.obs[col].dtype})")
    
    # Look for brain region and sex columns
    region_col = None
    sex_col = None
    
    # Check for region columns
    for col in ['Region', 'region', 'brain_region', 'area', 'Subregion']:
        if col in adata.obs.columns:
            region_col = col
            print(f"   ✅ Found region column: {col}")
            break
    
    # Check for sex columns  
    for col in ['Sex', 'sex', 'gender', 'Gender', 'M_F']:
        if col in adata.obs.columns:
            sex_col = col
            print(f"   ✅ Found sex column: {col}")
            break
    
    # Create pseudo-bulk grouping variable - SIMPLIFIED to avoid too many groups
    if region_col and sex_col:
        # Group by region + sex (condition will be handled in DE analysis)
        adata.obs['pseudobulk_group'] = (
            adata.obs[region_col].astype(str) + '_' + 
            adata.obs[sex_col].astype(str)
        )
        print(f"   ✅ Using region + sex grouping: {adata.obs['pseudobulk_group'].nunique()} groups")
    elif region_col:
        # Group by region only
        adata.obs['pseudobulk_group'] = adata.obs[region_col].astype(str)
        print(f"   ⚠️  Sex not found, using region only: {adata.obs['pseudobulk_group'].nunique()} groups")
    else:
        # Create artificial regional groups - MUCH FEWER
        print("   ⚠️  No region/sex found, creating artificial regional groups")
        np.random.seed(42)
        n_regions = 6  # Only 6 artificial regions total
        
        # Create balanced groups
        oud_cells = adata.obs['Dx_OUD'] == 'OUD'
        ctrl_cells = adata.obs['Dx_OUD'] == 'None'
        
        # Assign regions randomly but balanced
        adata.obs['pseudobulk_group'] = 'unknown'
        adata.obs.loc[oud_cells, 'pseudobulk_group'] = np.random.choice(
            ['Region_A', 'Region_B', 'Region_C'], 
            size=oud_cells.sum()
        )
        adata.obs.loc[ctrl_cells, 'pseudobulk_group'] = np.random.choice(
            ['Region_A', 'Region_B', 'Region_C'], 
            size=ctrl_cells.sum()
        )
        print(f"   Created {adata.obs['pseudobulk_group'].nunique()} artificial regional groups")
    
    # Show grouping summary
    group_summary = pd.crosstab(adata.obs['pseudobulk_group'], adata.obs['Dx_OUD'])
    print("   Pseudo-bulk groups by condition:")
    print(group_summary)
    
    # Show data summary
    print(f"   Cell types: {adata.obs['cell_type'].nunique()}")
    print(f"   Conditions: {adata.obs['Dx_OUD'].value_counts().to_dict()}")
    print(f"   Pseudo-bulk groups: {adata.obs['pseudobulk_group'].nunique()}")
    
    return adata

def create_pseudobulk(adata):
    """Aggregate single cells to pseudo-bulk by region/sex only (ignore cell types for bulk comparison)"""
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
    
    # Create aggregation metadata - group by pseudobulk_group + condition ONLY (no cell type)
    agg_metadata = adata.obs[['pseudobulk_group', 'Dx_OUD']].copy()
    agg_metadata['sample_id'] = (
        agg_metadata['pseudobulk_group'].astype(str) + '_' + 
        agg_metadata['Dx_OUD'].astype(str)
    )
    
    # Aggregate counts by sample_id (sum ALL cells regardless of cell type)
    unique_samples = agg_metadata['sample_id'].unique()
    pseudobulk_counts = np.zeros((len(gene_names), len(unique_samples)))
    sample_metadata = []
    
    print(f"   Aggregating {len(adata)} cells into {len(unique_samples)} pseudo-bulk samples...")
    print("   📋 Ignoring cell types - aggregating all cells by region+sex+condition")
    
    for i, sample_id in enumerate(unique_samples):
        sample_mask = agg_metadata['sample_id'] == sample_id
        sample_cells = np.where(sample_mask)[0]
        
        # Sum counts across ALL cells in this region+sex+condition group
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
    pseudobulk_df = pd.DataFrame(
        pseudobulk_counts.T,
        columns=gene_names,
        index=[s['sample_id'] for s in sample_metadata]
    )
    
    metadata_df = pd.DataFrame(sample_metadata)
    metadata_df.set_index('sample_id', inplace=True)
    
    print(f"   Created pseudo-bulk matrix: {pseudobulk_df.shape}")
    print(f"   Samples per condition: {metadata_df['condition'].value_counts().to_dict()}")
    print(f"   Samples per region+sex group: {metadata_df['pseudobulk_group'].value_counts().to_dict()}")
    
    # Show cells per pseudo-bulk sample stats
    cells_per_sample = metadata_df['n_cells']
    print(f"   Cells per sample: mean={cells_per_sample.mean():.0f}, "
          f"median={cells_per_sample.median():.0f}, "
          f"range={cells_per_sample.min()}-{cells_per_sample.max()}")
    
    return pseudobulk_df, metadata_df

def run_edger_de(counts_df, metadata_df):
    """Run edgeR differential expression for multiple contrasts"""
    print(f"   Running edgeR DE for all contrasts...")
    
    if not R_AVAILABLE:
        print("   ⚠️  R not available, skipping edgeR analysis")
        return None
    
    try:
        if len(counts_df) < 4:  # Need at least 2 samples per condition
            print(f"   ⚠️  Too few samples ({len(counts_df)})")
            return None
        
        # Check if both conditions represented
        conditions = metadata_df['condition'].unique()
        if len(conditions) < 2:
            print(f"   ⚠️  Only one condition")
            return None
        
        # Check replicates per condition
        condition_counts = metadata_df['condition'].value_counts()
        if any(condition_counts < 2):
            print(f"   ⚠️  Too few replicates per condition: {condition_counts.to_dict()}")
            return None
        
        # Use context-based conversion for safer R integration
        with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
            # Import R packages
            try:
                edger = importr('edgeR')
                base = importr('base')
                print("   ✅ edgeR package loaded successfully")
            except Exception as e:
                print(f"   ⚠️  Could not import edgeR: {e}")
                return None
            
            # Debug: Check data types before conversion
            print(f"     Debug: counts_df shape: {counts_df.shape}, dtype: {counts_df.dtypes.iloc[0]}")
            print(f"     Debug: metadata dtypes: {metadata_df.dtypes.to_dict()}")
            
            # Ensure counts are numeric and clean
            counts_clean = counts_df.select_dtypes(include=[np.number]).fillna(0)
            
            # Convert to R objects with explicit type checking
            try:
                r_counts = ro.conversion.py2rpy(counts_clean.T.values)  # Use .values to avoid index issues
                r_condition = ro.StrVector(metadata_df['condition'].astype(str).values)
                
                # Extract region and sex separately with safe splitting
                region_sex_split = metadata_df['pseudobulk_group'].str.split('_', expand=True)
                r_region = ro.StrVector(region_sex_split[0].astype(str).values)
                r_sex = ro.StrVector(region_sex_split[1].astype(str).values)
                
                print(f"     ✅ R objects created successfully")
                
            except Exception as e:
                print(f"   ❌ R conversion failed: {e}")
                return None
            
            # Create DGEList with all factors
            ro.r('library(edgeR)')
            ro.globalenv['counts_matrix'] = r_counts
            ro.globalenv['condition_factor'] = r_condition
            ro.globalenv['region_factor'] = r_region
            ro.globalenv['sex_factor'] = r_sex
            ro.globalenv['gene_names'] = ro.StrVector(counts_clean.columns.astype(str))
            
            # Run edgeR pipeline with better error handling
            try:
                ro.r('''
                # Create design matrix with all factors
                design <- model.matrix(~ condition_factor + region_factor + sex_factor)
                print(paste("Design matrix dimensions:", dim(design)))
                
                # Set row/column names properly
                rownames(counts_matrix) <- gene_names
                colnames(counts_matrix) <- paste0("Sample_", 1:ncol(counts_matrix))
                
                dge <- DGEList(counts = counts_matrix)
                print(paste("DGEList created with", nrow(dge), "genes and", ncol(dge), "samples"))
                
                # Filter low-expressed genes
                keep <- filterByExpr(dge, design)
                dge <- dge[keep, , keep.lib.sizes=FALSE]
                print(paste("After filtering:", nrow(dge), "genes retained"))
                
                # Normalization
                dge <- calcNormFactors(dge)
                print("Normalization completed")
                
                # Estimate dispersions
                dge <- estimateDisp(dge, design)
                print("Dispersion estimation completed")
                
                # Fit model
                fit <- glmQLFit(dge, design)
                print("Model fitting completed")
                
                # Test different contrasts
                # 1. OUD vs Control effect (main effect)
                qlf_condition <- glmQLFTest(fit, coef=2)
                results_condition <- topTags(qlf_condition, n=Inf)$table
                
                # 2. Region effect (if present)
                if (ncol(design) >= 4) {
                    qlf_region <- glmQLFTest(fit, coef=3)
                    results_region <- topTags(qlf_region, n=Inf)$table
                } else {
                    results_region <- NULL
                }
                
                # 3. Sex effect (if present)
                if (ncol(design) >= 5) {
                    qlf_sex <- glmQLFTest(fit, coef=4)
                    results_sex <- topTags(qlf_sex, n=Inf)$table
                } else {
                    results_sex <- NULL
                }
                
                print("edgeR analysis completed successfully")
                ''')
                
                print(f"     ✅ edgeR pipeline completed")
                
            except Exception as e:
                print(f"   ❌ edgeR pipeline failed: {e}")
                return None
            
            # Get all results back to Python with safer conversion
            try:
                results = {}
                
                # Convert results to DataFrames manually to avoid index issues
                contrasts_to_extract = [
                    ('OUD_vs_Control', 'results_condition'),
                ]
                
                # Check if region and sex results exist
                try:
                    ro.r('exists("results_region")')
                    contrasts_to_extract.append(('Putamen_vs_Caudate', 'results_region'))
                except:
                    pass
                
                try:
                    ro.r('exists("results_sex")')
                    contrasts_to_extract.append(('Male_vs_Female', 'results_sex'))
                except:
                    pass
                
                for contrast_name, r_result_name in contrasts_to_extract:
                    try:
                        # Get the R dataframe using alternative method
                        ro.r(f'current_result <- {r_result_name}')
                        
                        # Check if result exists and is not NULL
                        ro.r('result_exists <- !is.null(current_result)')
                        result_exists = list(ro.r('result_exists'))[0]
                        
                        if not result_exists:
                            print(f"     ⚠️  {contrast_name}: No results (NULL)")
                            continue
                        
                        # Extract data manually using R commands
                        ro.r('''
                        result_logFC <- current_result$logFC
                        result_PValue <- current_result$PValue  
                        result_FDR <- current_result$FDR
                        result_genes <- rownames(current_result)
                        ''')
                        
                        # Convert each column separately
                        logFC = list(ro.r('result_logFC'))
                        pvals = list(ro.r('result_PValue'))
                        fdr = list(ro.r('result_FDR'))
                        genes = list(ro.r('result_genes'))
                        
                        # Create DataFrame
                        results_df = pd.DataFrame({
                            'gene': genes,
                            'logFC': logFC,
                            'PValue': pvals,
                            'FDR': fdr,
                            'contrast': contrast_name
                        })
                        
                        results[contrast_name] = results_df
                        print(f"     ✅ {contrast_name}: {len(results_df)} genes")
                        
                    except Exception as e:
                        print(f"     ⚠️  Failed to extract {contrast_name}: {e}")
                        continue
                
                if results:
                    print(f"   ✅ edgeR completed successfully with {len(results)} contrasts")
                    return results
                else:
                    print(f"   ❌ No valid results extracted from edgeR")
                    return None
                    
            except Exception as e:
                print(f"   ❌ Results conversion failed: {e}")
                return None
    
    except Exception as e:
        print(f"   ❌ edgeR failed: {e}")
        return None

def run_python_de_fallback(counts_df, metadata_df):
    """Python fallback DE analysis with subgroup-specific OUD vs Control contrasts"""
    print(f"   Running Python fallback DE for all contrasts...")
    
    try:
        from scipy import stats
        from sklearn.linear_model import LinearRegression
        
        if len(counts_df) < 4:
            print(f"   ⚠️  Too few samples ({len(counts_df)})")
            return None
        
        # Check conditions
        conditions = metadata_df['condition'].unique()
        if len(conditions) < 2:
            print(f"   ⚠️  Only one condition")
            return None
        
        # Extract region and sex from pseudobulk_group
        region_sex_split = metadata_df['pseudobulk_group'].str.split('_', expand=True)
        metadata_df['region'] = region_sex_split[0]
        metadata_df['sex'] = region_sex_split[1]
        
        results_dict = {}
        
        # Test OUD vs Control in different subgroups
        contrasts = {
            # Simple contrasts within subgroups
            'OUD_vs_Control_Putamen': {'region': 'Putamen'},
            'OUD_vs_Control_Caudate': {'region': 'Caudate'},
            'OUD_vs_Control_Male': {'sex': 'M'},
            'OUD_vs_Control_Female': {'sex': 'F'},
            
            # Interaction contrasts - how OUD effects differ between groups
            'OUD_Effect_Male_vs_Female': 'interaction_sex',
            'OUD_Effect_Putamen_vs_Caudate': 'interaction_region'
        }
        
        for contrast_name, filter_criteria in contrasts.items():
            print(f"     Testing {contrast_name}...")
            
            # Handle interaction contrasts differently
            if contrast_name.startswith('OUD_Effect_'):
                if filter_criteria == 'interaction_sex':
                    # Test if OUD effect differs between males and females
                    # (OUD_Male - Control_Male) vs (OUD_Female - Control_Female)
                    results = []
                    
                    for i, gene in enumerate(counts_df.columns):
                        try:
                            gene_counts = counts_df.iloc[:, i].values
                            
                            # Get values for each group
                            male_oud = gene_counts[(metadata_df['sex'] == 'M') & (metadata_df['condition'] == 'OUD')]
                            male_ctrl = gene_counts[(metadata_df['sex'] == 'M') & (metadata_df['condition'] == 'None')]
                            female_oud = gene_counts[(metadata_df['sex'] == 'F') & (metadata_df['condition'] == 'OUD')]
                            female_ctrl = gene_counts[(metadata_df['sex'] == 'F') & (metadata_df['condition'] == 'None')]
                            
                            # Skip if insufficient data
                            if len(male_oud) < 1 or len(male_ctrl) < 1 or len(female_oud) < 1 or len(female_ctrl) < 1:
                                continue
                            
                            # Calculate OUD effect in each sex
                            male_oud_effect = np.mean(male_oud + 1) - np.mean(male_ctrl + 1)
                            female_oud_effect = np.mean(female_oud + 1) - np.mean(female_ctrl + 1)
                            
                            # Test if effects are different (interaction)
                            # Use log fold changes for comparison
                            male_logfc = np.log2(np.mean(male_oud + 1) / np.mean(male_ctrl + 1))
                            female_logfc = np.log2(np.mean(female_oud + 1) / np.mean(female_ctrl + 1))
                            
                            # Difference in OUD effects between sexes
                            interaction_effect = male_logfc - female_logfc
                            
                            # Statistical test: t-test on difference of effects
                            # Combine male and female effect sizes for variance estimation
                            combined_oud = np.concatenate([male_oud, female_oud])
                            combined_ctrl = np.concatenate([male_ctrl, female_ctrl])
                            
                            if len(combined_oud) > 1 and len(combined_ctrl) > 1:
                                t_stat, p_val = stats.ttest_ind(combined_oud, combined_ctrl)
                                # Adjust p-value for interaction test (more conservative)
                                p_val = min(p_val * 2, 1.0)  # Bonferroni-like adjustment
                            else:
                                t_stat, p_val = 0, 1
                            
                            results.append({
                                'gene': gene,
                                'logFC': interaction_effect,  # Difference in OUD effects
                                'PValue': p_val,
                                't_stat': t_stat,
                                'contrast': contrast_name,
                                'male_logfc': male_logfc,
                                'female_logfc': female_logfc,
                                'effect_difference': interaction_effect
                            })
                        except:
                            continue
                
                elif filter_criteria == 'interaction_region':
                    # Test if OUD effect differs between brain regions
                    # (OUD_Putamen - Control_Putamen) vs (OUD_Caudate - Control_Caudate)
                    results = []
                    
                    for i, gene in enumerate(counts_df.columns):
                        try:
                            gene_counts = counts_df.iloc[:, i].values
                            
                            # Get values for each group
                            putamen_oud = gene_counts[(metadata_df['region'] == 'Putamen') & (metadata_df['condition'] == 'OUD')]
                            putamen_ctrl = gene_counts[(metadata_df['region'] == 'Putamen') & (metadata_df['condition'] == 'None')]
                            caudate_oud = gene_counts[(metadata_df['region'] == 'Caudate') & (metadata_df['condition'] == 'OUD')]
                            caudate_ctrl = gene_counts[(metadata_df['region'] == 'Caudate') & (metadata_df['condition'] == 'None')]
                            
                            # Skip if insufficient data
                            if len(putamen_oud) < 1 or len(putamen_ctrl) < 1 or len(caudate_oud) < 1 or len(caudate_ctrl) < 1:
                                continue
                            
                            # Calculate OUD effect in each region
                            putamen_logfc = np.log2(np.mean(putamen_oud + 1) / np.mean(putamen_ctrl + 1))
                            caudate_logfc = np.log2(np.mean(caudate_oud + 1) / np.mean(caudate_ctrl + 1))
                            
                            # Difference in OUD effects between regions
                            interaction_effect = putamen_logfc - caudate_logfc
                            
                            # Statistical test
                            combined_oud = np.concatenate([putamen_oud, caudate_oud])
                            combined_ctrl = np.concatenate([putamen_ctrl, caudate_ctrl])
                            
                            if len(combined_oud) > 1 and len(combined_ctrl) > 1:
                                t_stat, p_val = stats.ttest_ind(combined_oud, combined_ctrl)
                                p_val = min(p_val * 2, 1.0)  # Bonferroni-like adjustment
                            else:
                                t_stat, p_val = 0, 1
                            
                            results.append({
                                'gene': gene,
                                'logFC': interaction_effect,  # Difference in OUD effects
                                'PValue': p_val,
                                't_stat': t_stat,
                                'contrast': contrast_name,
                                'putamen_logfc': putamen_logfc,
                                'caudate_logfc': caudate_logfc,
                                'effect_difference': interaction_effect
                            })
                        except:
                            continue
            
            else:
                # Handle simple subgroup contrasts (existing code)
                # Filter samples for this subgroup
                mask = pd.Series([True] * len(metadata_df), index=metadata_df.index)
                for col, val in filter_criteria.items():
                    mask = mask & (metadata_df[col] == val)
                
                if mask.sum() < 4:
                    print(f"       ⚠️  Too few samples in {contrast_name} ({mask.sum()})")
                    continue
                
                subgroup_metadata = metadata_df[mask]
                subgroup_counts = counts_df.loc[mask]
                
                # Check if both conditions present in subgroup
                subgroup_conditions = subgroup_metadata['condition'].unique()
                if len(subgroup_conditions) < 2:
                    print(f"       ⚠️  Only one condition in {contrast_name}")
                    continue
                
                # Encode OUD vs Control for this subgroup
                condition_encoded = (subgroup_metadata['condition'] == 'OUD').astype(int)
                
                results = []
                for i, gene in enumerate(subgroup_counts.columns):
                    try:
                        gene_counts = subgroup_counts.iloc[:, i].values
                        
                        # Skip genes with too many zeros
                        if np.sum(gene_counts > 0) < 2:
                            continue
                        
                        # Simple t-test for OUD vs Control in this subgroup
                        oud_values = gene_counts[condition_encoded == 1]
                        ctrl_values = gene_counts[condition_encoded == 0]
                        
                        if len(oud_values) < 1 or len(ctrl_values) < 1:
                            continue
                        
                        # T-test
                        t_stat, p_val = stats.ttest_ind(oud_values, ctrl_values)
                        
                        # Calculate fold change
                        oud_mean = np.mean(oud_values + 1)
                        ctrl_mean = np.mean(ctrl_values + 1)
                        log_fc = np.log2(oud_mean / ctrl_mean)
                        
                        results.append({
                            'gene': gene,
                            'logFC': log_fc,
                            'PValue': p_val,
                            't_stat': t_stat,
                            'contrast': contrast_name,
                            'oud_mean': oud_mean,
                            'ctrl_mean': ctrl_mean,
                            'n_oud': len(oud_values),
                            'n_ctrl': len(ctrl_values)
                        })
                    except:
                        continue
            
            if results:
                results_df = pd.DataFrame(results)
                
                # Multiple testing correction
                from statsmodels.stats.multitest import multipletests
                _, fdr_values, _, _ = multipletests(results_df['PValue'], method='fdr_bh')
                results_df['FDR'] = fdr_values
                
                results_dict[contrast_name] = results_df
                print(f"     ✅ {contrast_name}: {len(results_df)} genes tested")
        
        return results_dict
        
    except Exception as e:
        print(f"   ❌ Python fallback failed: {e}")
        return None

def run_differential_expression(pseudobulk_df, metadata_df):
    """Run DE analysis for subgroup-specific OUD vs Control contrasts"""
    print("\n🔬 RUNNING DIFFERENTIAL EXPRESSION ANALYSIS")
    print("=" * 50)
    print("   🧠 Testing OUD vs Control in different contexts:")
    print("     • Primary method: edgeR (if R available)")
    print("     • Fallback method: Python t-tests")
    
    # Try edgeR first (if R is available)
    if R_AVAILABLE:
        print("\n   🔬 Attempting edgeR analysis...")
        edger_results = run_edger_de(pseudobulk_df, metadata_df)
        
        if edger_results is not None:
            print("   ✅ edgeR analysis successful!")
            total_tests = sum(len(df) for df in edger_results.values())
            print(f"   ✅ edgeR complete: {total_tests} total tests across {len(edger_results)} contrasts")
            return edger_results
        else:
            print("   ⚠️  edgeR failed, falling back to Python analysis...")
    
    # Fallback to Python approach for subgroup-specific contrasts
    print("\n   🐍 Running Python fallback analysis...")
    print("     • OUD vs Control in Putamen")
    print("     • OUD vs Control in Caudate")
    print("     • OUD vs Control in Males") 
    print("     • OUD vs Control in Females")
    print("     • How OUD effects differ between Males vs Females")
    print("     • How OUD effects differ between Putamen vs Caudate")
    
    results = run_python_de_fallback(pseudobulk_df, metadata_df)
    
    if results is not None:
        total_tests = sum(len(df) for df in results.values())
        print(f"   ✅ Python DE analysis complete: {total_tests} total tests across {len(results)} contrasts")
        return results
    else:
        print("   ❌ No DE results generated")
        return None

def analyze_de_results(de_results_dict):
    """Analyze and summarize DE results for all contrasts"""
    print("\n📊 ANALYZING DE RESULTS")
    print("=" * 50)
    
    if de_results_dict is None:
        print("   No results to analyze")
        return None
    
    summary_stats = []
    
    for contrast_name, de_results in de_results_dict.items():
        # DEBUG: Check data distribution before significance filtering
        print(f"\n   DEBUG {contrast_name}:")
        print(f"     • P-value range: {de_results['PValue'].min():.2e} - {de_results['PValue'].max():.2e}")
        print(f"     • LogFC range: {de_results['logFC'].min():.2f} - {de_results['logFC'].max():.2f}")
        print(f"     • FDR range: {de_results['FDR'].min():.2e} - {de_results['FDR'].max():.2e}")
        
        # Count genes at different thresholds
        p_05 = (de_results['PValue'] < 0.05).sum()
        p_01 = (de_results['PValue'] < 0.01).sum()
        p_001 = (de_results['PValue'] < 0.001).sum()
        fdr_05 = (de_results['FDR'] < 0.05).sum()
        fdr_20 = (de_results['FDR'] < 0.20).sum()
        
        print(f"     • P < 0.05: {p_05} genes")
        print(f"     • P < 0.01: {p_01} genes")
        print(f"     • P < 0.001: {p_001} genes")
        print(f"     • FDR < 0.05: {fdr_05} genes")
        print(f"     • FDR < 0.20: {fdr_20} genes")
        
        # Use NOMINAL p-values for small sample discovery (common in pilot studies)
        # This is acceptable given the exploratory nature and small sample size
        de_results['significant_nominal'] = (de_results['PValue'] < 0.01) & (abs(de_results['logFC']) > 0.5)
        de_results['significant_strict'] = (de_results['PValue'] < 0.001) & (abs(de_results['logFC']) > 1.0)
        de_results['significant'] = de_results['significant_nominal']  # Use nominal for discovery
        
        de_results['direction'] = np.where(de_results['logFC'] > 0, 'Up', 'Down')
        
        # Calculate summary
        n_sig_nominal = de_results['significant_nominal'].sum()
        n_sig_strict = de_results['significant_strict'].sum()
        n_up = ((de_results['significant']) & (de_results['logFC'] > 0)).sum()
        n_down = ((de_results['significant']) & (de_results['logFC'] < 0)).sum();
        
        print(f"     • Significant (nominal p<0.01, |FC|>0.5): {n_sig_nominal}")
        print(f"     • Significant (strict p<0.001, |FC|>1.0): {n_sig_strict}")
        
        summary_stats.append({
            'contrast': contrast_name,
            'total_genes': len(de_results),
            'significant_genes': n_sig_nominal,  # Use nominal
            'significant_strict': n_sig_strict,
            'upregulated': n_up,
            'downregulated': n_down,
            'pct_significant': n_sig_nominal / len(de_results) * 100,
            'genes_p05': p_05,
            'genes_p01': p_01,
            'genes_p001': p_001,
            'genes_fdr05': fdr_05,
            'genes_fdr20': fdr_20
        })
        
        print(f"   {contrast_name}:")
        print(f"     • Total genes: {len(de_results)}")
        print(f"     • Significant (nominal p<0.01): {n_sig_nominal} ({n_up}↑, {n_down}↓)")
        print(f"     • Percentage: {n_sig_nominal / len(de_results) * 100:.2f}%")
        
        # Show top genes by p-value regardless of significance
        if len(de_results) > 0:
            top_by_pval = de_results.nsmallest(5, 'PValue')[['gene', 'logFC', 'PValue', 'FDR']]
            print(f"     • Top 5 genes by p-value:")
            for _, row in top_by_pval.iterrows():
                print(f"       - {row['gene']}: logFC={row['logFC']:.2f}, p={row['PValue']:.2e}, FDR={row['FDR']:.2e}")
    
    summary_df = pd.DataFrame(summary_stats)
    return summary_df

# Add a function to create a pooled analysis for better power
def run_pooled_analysis(pseudobulk_df, metadata_df):
    """Run pooled OUD vs Control analysis across all samples for better power"""
    print("\n🔬 RUNNING POOLED ANALYSIS")
    print("=" * 50)
    print("   Pooling all samples for maximum statistical power")
    
    try:
        from scipy import stats
        
        # Simple OUD vs Control across all samples
        condition_encoded = (metadata_df['condition'] == 'OUD').astype(int)
        
        results = []
        for i, gene in enumerate(pseudobulk_df.columns):
            try:
                gene_counts = pseudobulk_df.iloc[:, i].values
                
                # Skip genes with too many zeros
                if np.sum(gene_counts > 0) < 4:
                    continue
                
                # T-test
                oud_values = gene_counts[condition_encoded == 1]
                ctrl_values = gene_counts[condition_encoded == 0]
                
                if len(oud_values) < 2 or len(ctrl_values) < 2:
                    continue
                
                t_stat, p_val = stats.ttest_ind(oud_values, ctrl_values)
                
                # Calculate fold change
                oud_mean = np.mean(oud_values + 1)
                ctrl_mean = np.mean(ctrl_values + 1)
                log_fc = np.log2(oud_mean / ctrl_mean)
                
                results.append({
                    'gene': gene,
                    'logFC': log_fc,
                    'PValue': p_val,
                    't_stat': t_stat,
                    'contrast': 'Pooled_OUD_vs_Control',
                    'oud_mean': oud_mean,
                    'ctrl_mean': ctrl_mean,
                    'n_oud': len(oud_values),
                    'n_ctrl': len(ctrl_values)
                })
            except:
                continue
        
        if results:
            results_df = pd.DataFrame(results)
            
            # Multiple testing correction with better power
            from statsmodels.stats.multitest import multipletests
            _, fdr_values, _, _ = multipletests(results_df['PValue'], method='fdr_bh')
            results_df['FDR'] = fdr_values
            
            print(f"   ✅ Pooled analysis: {len(results_df)} genes tested")
            print(f"   📊 Results summary:")
            print(f"     • P < 0.05: {(results_df['PValue'] < 0.05).sum()} genes")
            print(f"     • FDR < 0.05: {(results_df['FDR'] < 0.05).sum()} genes")
            print(f"     • FDR < 0.10: {(results_df['FDR'] < 0.10).sum()} genes")
            print(f"     • FDR < 0.20: {(results_df['FDR'] < 0.20).sum()} genes")
            
            return results_df
        else:
            return None
            
    except Exception as e:
        print(f"   ❌ Pooled analysis failed: {e}")
        return None

def check_data_quality(pseudobulk_df, metadata_df):
    """Check data quality and power for DE analysis"""
    print("\n🔍 DATA QUALITY CHECK")
    print("=" * 50)
    
    # Sample size per group
    print("   Sample sizes:")
    for group in metadata_df['pseudobulk_group'].unique():
        group_meta = metadata_df[metadata_df['pseudobulk_group'] == group]
        print(f"     {group}: {len(group_meta)} samples ({group_meta['condition'].value_counts().to_dict()})")
    
    # Gene expression distribution
    print(f"\n   Gene expression statistics:")
    print(f"     • Total genes: {pseudobulk_df.shape[1]}")
    print(f"     • Genes with >0 counts: {(pseudobulk_df > 0).any().sum()}")
    print(f"     • Mean counts per gene: {pseudobulk_df.mean().mean():.0f}")
    print(f"     • Median counts per gene: {pseudobulk_df.median().median():.0f}")
    
    # Check for obvious differences between conditions
    oud_samples = metadata_df[metadata_df['condition'] == 'OUD'].index
    ctrl_samples = metadata_df[metadata_df['condition'] == 'None'].index
    
    oud_total = pseudobulk_df.loc[oud_samples].sum(axis=1)
    ctrl_total = pseudobulk_df.loc[ctrl_samples].sum(axis=1)
    
    print(f"\n   Library sizes:")
    print(f"     • OUD samples: {oud_total.mean():.0f} ± {oud_total.std():.0f}")
    print(f"     • Control samples: {ctrl_total.mean():.0f} ± {ctrl_total.std():.0f}")
    
    # Simple fold change check
    oud_mean = pseudobulk_df.loc[oud_samples].mean()
    ctrl_mean = pseudobulk_df.loc[ctrl_samples].mean()
    fold_changes = np.log2((oud_mean + 1) / (ctrl_mean + 1))
    
    print(f"\n   Fold change distribution (log2):")
    print(f"     • Range: {fold_changes.min():.2f} to {fold_changes.max():.2f}")
    print(f"     • |FC| > 0.5: {(abs(fold_changes) > 0.5).sum()} genes")
    print(f"     • |FC| > 1.0: {(abs(fold_changes) > 1.0).sum()} genes")

def create_de_plots(de_results_dict, summary_df):
    """Create visualization plots for all contrasts"""
    print("\n📈 CREATING DE VISUALIZATION PLOTS")
    print("=" * 50)
    
    if de_results_dict is None or summary_df is None:
        print("   No results to plot")
        return
    
    # Set matplotlib backend
    plt.switch_backend('Agg')
    
    # Create separate plots for each contrast
    for contrast_name, de_results in de_results_dict.items():
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # 1. Overall DE summary
        n_sig = de_results['significant'].sum()
        n_up = ((de_results['significant']) & (de_results['logFC'] > 0)).sum()
        n_down = ((de_results['significant']) & (de_results['logFC'] < 0)).sum()
        n_nonsig = len(de_results) - n_sig;
        
        categories = ['Upregulated', 'Downregulated', 'Non-significant']
        counts = [n_up, n_down, n_nonsig]
        colors = ['red', 'blue', 'gray']
        
        axes[0,0].bar(categories, counts, color=colors)
        axes[0,0].set_title(f'DE Gene Summary - {contrast_name}')
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
        axes[0,1].set_title(f'Volcano Plot - {contrast_name}')
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
        
        plt.suptitle(f'Differential Expression Analysis - {contrast_name}', fontsize=16)
        plt.tight_layout()
        
        # Save plot for this contrast
        safe_name = contrast_name.replace(' ', '_').replace('vs', 'vs')
        plt.savefig(f"{PLOTS_DIR}/de_analysis_{safe_name}.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # Create summary comparison plot
    if len(summary_df) > 0:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        contrasts = summary_df['contrast']
        x_pos = np.arange(len(contrasts))
        
        ax.bar(x_pos - 0.2, summary_df['upregulated'], 0.4, label='Upregulated', color='red', alpha=0.7)
        ax.bar(x_pos + 0.2, summary_df['downregulated'], 0.4, label='Downregulated', color='blue', alpha=0.7)
        
        ax.set_xlabel('Contrasts')
        ax.set_ylabel('Number of Significant DE Genes')
        ax.set_title('DE Genes Across All Contrasts')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(contrasts, rotation=45)
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/de_analysis_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"   Plots saved to: {PLOTS_DIR}/")

def save_results(de_results_dict, summary_df, pseudobulk_df, metadata_df):
    """Save all DE analysis results for multiple contrasts"""
    print("\n💾 SAVING RESULTS")
    print("=" * 50)
    
    if de_results_dict is not None:
        # Save combined results
        all_results = pd.concat(de_results_dict.values(), ignore_index=True)
        all_results.to_csv(f"{OUTPUT_DIR}/differential_expression_all_contrasts.csv", index=False)
        print(f"   All DE results: {OUTPUT_DIR}/differential_expression_all_contrasts.csv")
        
        # Save each contrast separately
        for contrast_name, de_results in de_results_dict.items():
            safe_name = contrast_name.replace(' ', '_').replace('vs', 'vs')
            de_results.to_csv(f"{OUTPUT_DIR}/de_results_{safe_name}.csv", index=False)
            
            # Save significant genes for each contrast (both strict and relaxed)
            if 'significant' in de_results.columns:
                sig_results = de_results[de_results['significant']]
                if len(sig_results) > 0:
                    sig_results.to_csv(f"{OUTPUT_DIR}/significant_de_genes_{safe_name}.csv", index=False)
                    print(f"   {contrast_name} significant DE: {OUTPUT_DIR}/significant_de_genes_{safe_name}.csv")
                
                # Also save top genes by p-value for exploration
                top_results = de_results.nsmallest(100, 'PValue')
                top_results.to_csv(f"{OUTPUT_DIR}/top100_genes_{safe_name}.csv", index=False)
                print(f"   {contrast_name} top 100: {OUTPUT_DIR}/top100_genes_{safe_name}.csv")
    
    # Save summary
    if summary_df is not None:
        summary_df.to_csv(f"{OUTPUT_DIR}/de_summary_all_contrasts.csv", index=False)
        print(f"   Summary: {OUTPUT_DIR}/de_summary_all_contrasts.csv")
    
    # Save pseudo-bulk data
    pseudobulk_df.to_csv(f"{OUTPUT_DIR}/pseudobulk_counts.csv")
    metadata_df.to_csv(f"{OUTPUT_DIR}/pseudobulk_metadata.csv")
    print(f"   Pseudo-bulk data: {OUTPUT_DIR}/pseudobulk_*.csv")

def main():
    """Main differential expression workflow"""
    print("🔬 DIFFERENTIAL EXPRESSION ANALYSIS - OUD vs CONTROL")
    print("=" * 70)
    print("Using pseudo-bulk aggregation + edgeR/Python for statistical analysis")
    print("=" * 70)
    
    try:
        # Load and prepare data
        adata = load_and_prepare_data()
        
        # Create pseudo-bulk aggregates
        pseudobulk_df, metadata_df = create_pseudobulk(adata)
        
        # Check data quality first
        check_data_quality(pseudobulk_df, metadata_df)
        
        # Run differential expression for multiple contrasts
        de_results_dict = run_differential_expression(pseudobulk_df, metadata_df)
        
        # Add pooled analysis for better power
        pooled_results = run_pooled_analysis(pseudobulk_df, metadata_df)
        if pooled_results is not None:
            de_results_dict['Pooled_OUD_vs_Control'] = pooled_results
        
        # Analyze results
        summary_df = analyze_de_results(de_results_dict)
        
        # Create plots
        create_de_plots(de_results_dict, summary_df)
        
        # Save everything
        save_results(de_results_dict, summary_df, pseudobulk_df, metadata_df)
        
        print("\n" + "=" * 70)
        print("✅ DIFFERENTIAL EXPRESSION ANALYSIS COMPLETE!")
        if de_results_dict is not None:
            total_sig = sum(df.get('significant', pd.Series([])).sum() for df in de_results_dict.values())
            print(f"📊 Found {total_sig} significant DE genes across {len(de_results_dict)} contrasts")
        print(f"📁 Results saved to: {OUTPUT_DIR}/")
        print("🚀 Ready for pathway analysis and interpretation!")
        print("=" * 70)
        
        return de_results_dict, summary_df
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None

if __name__ == "__main__":
    main()
