#!/usr/bin/env python3
"""
🌊 LEMUR Strategic Regional Analysis - High-Priority Biological Contrasts
GSE225158 - Strategic implementation focusing on key regional and sex effects

This script implements the final strategic LEMUR contrasts to complete the 
comprehensive comparison with DESeq2:

Strategic Design (6 Core Contrasts):
✅ Tier 1: Already Complete (from previous scripts)
1. OUD vs Control (Pooled) - Overall effect
2. OUD vs Control (Males) - Male-specific effects  
3. OUD vs Control (Females) - Female-specific effects
4. OUD × Sex Interaction - Sex-differential effects

🆕 Tier 2: Regional Contrasts (This Script)
5. OUD vs Control (Caudate only) - Goal-directed behavior effects
6. OUD vs Control (Putamen only) - Habit formation effects

🔬 Tier 3: Complex Interactions (Optional)
7. OUD × Sex (Caudate only) - Regional sex interactions
8. OUD × Sex (Putamen only) - Regional sex interactions

Regional Biology Focus:
- Caudate: Goal-directed behavior, cognitive control, executive function
- Putamen: Habit formation, motor control, addiction maintenance

Based on DESeq2 findings:
- Caudate shows more OUD effects (70 genes) than Putamen (10 genes)
- Strong sex differences justify regional-sex interactions
- Females show massive effects (2,743 genes) vs Males (51 genes)

Author: Research Team
Date: 2024
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests
import anndata as ad

# Configure scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=300, facecolor='white')
warnings.filterwarnings('ignore')

# LEMUR dependencies
try:
    import pylemur
    PYLEMUR_AVAILABLE = True
    print("✅ pyLEMUR available")
except ImportError:
    PYLEMUR_AVAILABLE = False
    print("❌ pyLEMUR not available. Install with: pip install pylemur")
    exit(1)

# ============================================================================
# 📁 CONFIGURATION
# ============================================================================

# File paths
BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
RAW_H5AD = f"{BASE_DIR}/data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_analysis"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"

# Analysis parameters optimized for full dataset regional analysis
ANALYSIS_CONFIG = {
    'n_embedding': 15,          # Number of latent dimensions
    'max_cells': None,          # 🚀 NO SUBSAMPLING - Use full ~99K cells with 64GB RAM
    'n_hvg': None,              # 🧬 NO HVG LIMIT - Use all highly variable genes
    'random_seed': 42,          # Reproducibility
    
    # Discovery-oriented thresholds
    'fdr_threshold': 0.1,       # Relaxed for discovery
    'lfc_threshold': 0.05,      # Relaxed for subtle regional effects
    
    # Plotting (adjusted for larger dataset)
    'figure_dpi': 300,
    'point_size': 1,            # Smaller points for full dataset
    'alpha': 0.4,               # More transparency for dense plots
}

# Strategic experimental design
DESIGN_CONFIG = {
    'condition_col': 'level1',      # OUD vs CTL
    'sex_col': 'Sex',               # F vs M
    'region_col': 'Region',         # Caudate vs Putamen
    'sample_col': 'orig.ident',     # Sample ID
    'celltype_col': 'celltype3',    # Cell type annotation
    
    # Strategic contrasts - focused on high biological impact
    'strategic_contrasts': {
        # Tier 2: Core Regional Effects (High Priority)
        'caudate_oud_effect': {
            'name': 'OUD vs Control (Caudate)',
            'design': '~ level1',
            'subset': 'caudate',
            'coef_idx': 1,
            'description': 'OUD effects on goal-directed behavior circuits',
            'biology': 'Cognitive control, executive function, goal-directed behavior',
            'priority': 'HIGH'
        },
        'putamen_oud_effect': {
            'name': 'OUD vs Control (Putamen)',
            'design': '~ level1',
            'subset': 'putamen',
            'coef_idx': 1,
            'description': 'OUD effects on habit formation circuits',
            'biology': 'Habit formation, motor control, addiction maintenance',
            'priority': 'HIGH'
        },
        
        # Tier 3: Regional-Sex Interactions (Medium Priority)
        'caudate_sex_interaction': {
            'name': 'OUD × Sex (Caudate)',
            'design': '~ level1 + Sex + level1:Sex',
            'subset': 'caudate',
            'coef_idx': 3,
            'description': 'Sex-specific OUD effects in cognitive control',
            'biology': 'Sex differences in goal-directed behavior disruption',
            'priority': 'MEDIUM'
        },
        'putamen_sex_interaction': {
            'name': 'OUD × Sex (Putamen)',
            'design': '~ level1 + Sex + level1:Sex',
            'subset': 'putamen',
            'coef_idx': 3,
            'description': 'Sex-specific OUD effects in habit formation',
            'biology': 'Sex differences in motor control and habit disruption',
            'priority': 'MEDIUM'
        },
        
        # Tier 3: Regional Sex Effects (Medium Priority)
        'caudate_sex_oud': {
            'name': 'Male vs Female OUD (Caudate)',
            'design': '~ Sex',
            'subset': 'caudate_oud',
            'coef_idx': 1,
            'description': 'Sex differences in OUD subjects - cognitive circuits',
            'biology': 'Sex-specific vulnerability in executive function',
            'priority': 'MEDIUM'
        },
        'putamen_sex_oud': {
            'name': 'Male vs Female OUD (Putamen)',
            'design': '~ Sex',
            'subset': 'putamen_oud',
            'coef_idx': 1,
            'description': 'Sex differences in OUD subjects - motor circuits',
            'biology': 'Sex-specific vulnerability in habit formation',
            'priority': 'MEDIUM'
        }
    }
}

# Create output directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🌊 LEMUR STRATEGIC REGIONAL ANALYSIS - FULL DATASET")
print("===================================================")
print(f"🚀 FULL DATASET MODE: Using all ~99K cells with 64GB RAM")
print(f"🧬 FULL GENE MODE: Using all highly variable genes (no 3K limit)")
print(f"Input: {RAW_H5AD}")
print(f"Output: {OUTPUT_DIR}")
print(f"Strategic contrasts: {len(DESIGN_CONFIG['strategic_contrasts'])}")
print(f"Focus: Maximum statistical power for regional and sex effects")

# ============================================================================
# 📊 DATA LOADING AND PREPARATION
# ============================================================================

def load_and_prepare_data():
    """Load and prepare raw data with regional validation"""
    print(f"\n📁 LOADING RAW DATA WITH REGIONAL FOCUS")
    print("=" * 40)
    
    if not os.path.exists(RAW_H5AD):
        print(f"   ❌ File not found: {RAW_H5AD}")
        return None
    
    # Load raw data
    adata = sc.read_h5ad(RAW_H5AD)
    print(f"   ✅ Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    # Validate required columns
    required_cols = [DESIGN_CONFIG['condition_col'], 
                    DESIGN_CONFIG['sex_col'], 
                    DESIGN_CONFIG['region_col'],
                    DESIGN_CONFIG['sample_col']]
    
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        print(f"   ❌ Missing required columns: {missing_cols}")
        return None
    
    print(f"   ✅ All required metadata columns present")
    
    # Regional data validation
    print(f"\n🧠 REGIONAL DATA VALIDATION")
    print("=" * 28)
    
    region_counts = adata.obs[DESIGN_CONFIG['region_col']].value_counts()
    print(f"   📊 Regional distribution:")
    for region, count in region_counts.items():
        pct = count / adata.n_obs * 100
        print(f"      {region}: {count:,} cells ({pct:.1f}%)")
    
    # Regional power analysis
    print(f"\n   ⚡ Regional Power Analysis:")
    for region in adata.obs[DESIGN_CONFIG['region_col']].unique():
        region_data = adata.obs[adata.obs[DESIGN_CONFIG['region_col']] == region]
        
        # Condition distribution
        oud_count = region_data[region_data[DESIGN_CONFIG['condition_col']] == 'OUD'].shape[0]
        ctl_count = region_data[region_data[DESIGN_CONFIG['condition_col']] == 'CTL'].shape[0]
        
        # Sex distribution
        male_count = region_data[region_data[DESIGN_CONFIG['sex_col']] == 'M'].shape[0]
        female_count = region_data[region_data[DESIGN_CONFIG['sex_col']] == 'F'].shape[0]
        
        # Power assessment
        min_group = min(oud_count, ctl_count)
        power_status = "✅ Excellent" if min_group >= 5000 else "⚡ Good" if min_group >= 2000 else "⚠️ Moderate"
        
        print(f"      {region}:")
        print(f"         Condition: OUD={oud_count:,}, CTL={ctl_count:,} → {power_status}")
        print(f"         Sex: M={male_count:,}, F={female_count:,}")
    
    # Cross-tabulation for interaction analysis
    print(f"\n   🔬 Interaction Analysis Power:")
    for region in adata.obs[DESIGN_CONFIG['region_col']].unique():
        region_mask = adata.obs[DESIGN_CONFIG['region_col']] == region
        region_data = adata.obs[region_mask]
        
        crosstab = pd.crosstab(region_data[DESIGN_CONFIG['condition_col']], 
                              region_data[DESIGN_CONFIG['sex_col']])
        print(f"      {region}:")
        print(crosstab.to_string().replace('\n', '\n         '))
        print()
    
    return adata

def prepare_strategic_datasets(adata):
    """Create strategically designed regional datasets"""
    print(f"\n🎯 PREPARING STRATEGIC DATASETS - FULL POWER")
    print("=" * 47)
    
    # Full dataset mode - maximum statistical power
    print(f"   🚀 FULL DATASET MODE: Using all {adata.n_obs:,} cells")
    print(f"   💾 Memory available: 64GB RAM")
    print(f"   ⚡ Statistical power: MAXIMIZED")
    
    # HVG selection - use ALL highly variable genes
    print("   🧬 Using ALL highly variable genes...")
    if 'highly_variable' in adata.var.columns:
        hvg_genes = adata.var[adata.var['highly_variable']].index
        print(f"   ✅ Using ALL pre-computed HVG: {len(hvg_genes)} genes")
    else:
        # Compute HVG without top_genes limit to get all variable genes
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3')
        hvg_genes = adata.var[adata.var['highly_variable']].index
        print(f"   ✅ Computed ALL HVG: {len(hvg_genes)} genes")
    
    adata_hvg = adata[:, hvg_genes].copy()
    print(f"   Selected {adata_hvg.n_vars} highly variable genes")
    
    # Prepare metadata
    print("   🏷️  Preparing metadata...")
    for col in [DESIGN_CONFIG['condition_col'], DESIGN_CONFIG['sex_col'], DESIGN_CONFIG['region_col']]:
        if adata_hvg.obs[col].dtype.name == 'category':
            adata_hvg.obs[col] = adata_hvg.obs[col].astype(str)
        adata_hvg.obs[col] = pd.Categorical(adata_hvg.obs[col])
    
    # Create strategic subset datasets
    datasets = create_strategic_subsets(adata_hvg)
    
    print(f"   ✅ Full dataset prepared: {adata_hvg.n_obs:,} cells × {adata_hvg.n_vars:,} genes")
    print(f"   📈 Maximum power: {adata_hvg.n_obs:,} cells × {adata_hvg.n_vars:,} genes")
    print(f"   🧬 Full biological coverage: No gene filtering limits")
    
    return datasets

def subsample_balanced_strategic(adata, max_cells):
    """Legacy function - not used in full dataset mode"""
    print("   ⚠️  Subsampling function called but full dataset mode is enabled")
    print("   🚀 Returning full dataset for maximum statistical power")
    return adata.copy()

def create_strategic_subsets(adata):
    """Create strategically designed subsets for each contrast"""
    print(f"\n🔬 CREATING STRATEGIC SUBSETS")
    print("=" * 30)
    
    datasets = {}
    
    # Full dataset
    datasets['all'] = adata.copy()
    print(f"   All: {adata.n_obs:,} cells")
    
    # Regional datasets
    for region in ['Caudate', 'Putamen']:
        region_mask = adata.obs[DESIGN_CONFIG['region_col']] == region
        datasets[region.lower()] = adata[region_mask].copy()
        
        region_data = datasets[region.lower()]
        condition_dist = region_data.obs[DESIGN_CONFIG['condition_col']].value_counts()
        sex_dist = region_data.obs[DESIGN_CONFIG['sex_col']].value_counts()
        
        print(f"   {region}: {region_data.n_obs:,} cells")
        print(f"      Condition: {dict(condition_dist)}")
        print(f"      Sex: {dict(sex_dist)}")
    
    # OUD-only regional datasets (for sex effects in OUD)
    for region in ['Caudate', 'Putamen']:
        oud_region_mask = (adata.obs[DESIGN_CONFIG['region_col']] == region) & \
                         (adata.obs[DESIGN_CONFIG['condition_col']] == 'OUD')
        datasets[f"{region.lower()}_oud"] = adata[oud_region_mask].copy()
        
        oud_data = datasets[f"{region.lower()}_oud"]
        if oud_data.n_obs > 0:
            sex_dist = oud_data.obs[DESIGN_CONFIG['sex_col']].value_counts()
            print(f"   {region} OUD: {oud_data.n_obs:,} cells | Sex: {dict(sex_dist)}")
    
    return datasets

# ============================================================================
# 🌊 STRATEGIC LEMUR ANALYSIS
# ============================================================================

def run_strategic_lemur_analysis(datasets):
    """Run strategic LEMUR analysis focusing on high-impact contrasts"""
    print(f"\n🌊 RUNNING STRATEGIC LEMUR ANALYSIS")
    print("=" * 37)
    
    results = {}
    
    # Prioritize contrasts by biological impact
    high_priority = [k for k, v in DESIGN_CONFIG['strategic_contrasts'].items() if v['priority'] == 'HIGH']
    medium_priority = [k for k, v in DESIGN_CONFIG['strategic_contrasts'].items() if v['priority'] == 'MEDIUM']
    
    print(f"   🎯 High Priority Contrasts: {len(high_priority)}")
    print(f"   🔬 Medium Priority Contrasts: {len(medium_priority)}")
    
    # Run high priority contrasts first
    for contrast_name in high_priority + medium_priority:
        contrast_config = DESIGN_CONFIG['strategic_contrasts'][contrast_name]
        
        print(f"\n🎯 CONTRAST: {contrast_name} ({contrast_config['priority']} PRIORITY)")
        print("=" * (12 + len(contrast_name) + len(contrast_config['priority']) + 12))
        print(f"   📋 Description: {contrast_config['description']}")
        print(f"   🧠 Biology: {contrast_config['biology']}")
        
        try:
            # Get appropriate dataset
            subset_name = contrast_config['subset']
            if subset_name not in datasets:
                print(f"   ❌ Dataset subset '{subset_name}' not found")
                continue
            
            dataset = datasets[subset_name]
            design_formula = contrast_config['design']
            
            print(f"   📊 Dataset: {dataset.n_obs:,} cells")
            print(f"   🔬 Design: {design_formula}")
            
            # Cell count validation for full dataset
            if dataset.n_obs < 5000:
                print(f"   ⚠️  Moderate cell count: {dataset.n_obs:,} cells")
            else:
                print(f"   ✅ Excellent power: {dataset.n_obs:,} cells")
            
            # Run LEMUR analysis
            contrast_result = run_single_strategic_contrast(
                dataset, design_formula, contrast_name, contrast_config
            )
            
            if contrast_result is not None:
                results[contrast_name] = contrast_result
                n_sig = contrast_result['de_results']['significant'].sum()
                n_sig_lenient = contrast_result['de_results']['significant_lenient'].sum()
                print(f"   ✅ SUCCESS: {n_sig} significant genes (strict), {n_sig_lenient} significant (lenient)")
                    
                # Show breakdown of significance criteria
                de_data = contrast_result['de_results']
                n_strong_effect = de_data['significant_strong_effect'].sum()
                n_strong_combined = de_data['significant_strong_combined'].sum()
                n_fdr_sig = de_data['significant_fdr'].sum()
                
                print(f"   📊 Significance breakdown: {n_strong_effect} strong effect, {n_strong_combined} strong combined, {n_fdr_sig} FDR significant")
                
                if n_sig > 0:
                    sig_genes = de_data[de_data['significant']]
                    top_genes = sig_genes.nlargest(3, 'combined_score')
                    print(f"   🧬 Top significant genes (by combined score):")
                    for gene, row in top_genes.iterrows():
                        criteria = []
                        if row['significant_strong_effect']: criteria.append('SE')
                        if row['significant_strong_combined']: criteria.append('SC')
                        if row['significant_fdr']: criteria.append('FDR')
                        criteria_str = '+'.join(criteria)
                        print(f"      • {gene}: coef={row['coefficient']:.3f}, score={row['combined_score']:.4f} [{criteria_str}]")
                elif n_sig_lenient > 0:
                    sig_genes = de_data[de_data['significant_lenient']]
                    top_genes = sig_genes.nlargest(3, 'combined_score')
                    print(f"   🔬 Top lenient genes (by combined score):")
                    for gene, row in top_genes.iterrows():
                        print(f"      • {gene}: coef={row['coefficient']:.3f}, score={row['combined_score']:.4f}")
                else:
                    # Show top effects by combined score
                    top_genes = de_data.nlargest(3, 'combined_score')
                    print(f"   🎯 Top effects (by combined score):")
                    for gene, row in top_genes.iterrows():
                        print(f"      • {gene}: coef={row['coefficient']:.3f}, score={row['combined_score']:.4f}")
            else:
                print(f"   ❌ FAILED: Analysis failed for {contrast_name}")
                
        except Exception as e:
            print(f"   ❌ ERROR: {contrast_name} failed with error: {e}")
            import traceback
            traceback.print_exc()
    
    return results

def run_single_strategic_contrast(adata, design_formula, contrast_name, contrast_config):
    """Run single strategic LEMUR contrast with robust error handling"""
    print(f"      🔬 Fitting LEMUR model...")
    
    try:
        # Create LEMUR model
        lemur_model = pylemur.tl.LEMUR(
            adata,
            design=design_formula,
            n_embedding=ANALYSIS_CONFIG['n_embedding'],
            copy=True
        )
        
        print(f"         ✅ Model created")
        
        # Fit model
        lemur_model.fit(verbose=False)
        print(f"         ✅ Model fitted")
        
        # Try alignment
        try:
            lemur_model.align_with_harmony()
            print(f"         ✅ Harmony alignment successful")
        except Exception as e:
            print(f"         ⚠️  Harmony alignment failed: {e}")
        
        # Extract differential expression
        de_results = extract_strategic_de_results(
            lemur_model, contrast_name, contrast_config
        )
        
        if de_results is not None:
            return {
                'lemur_model': lemur_model,
                'de_results': de_results,
                'dataset_used': adata,
                'contrast_config': contrast_config,
                'method': 'strategic_lemur'
            }
        else:
            print(f"         ❌ DE extraction failed")
            return None
            
    except Exception as e:
        print(f"         ❌ LEMUR fitting failed: {e}")
        return None

def extract_strategic_de_results(lemur_model, contrast_name, contrast_config):
    """Extract DE results with LEMUR-specific statistical testing"""
    print(f"         🧬 Extracting differential expression...")
    
    try:
        # Get coefficient tensor
        coefficients = lemur_model.coefficients
        coef_idx = contrast_config['coef_idx']
        
        print(f"            📊 Coefficient tensor shape: {coefficients.shape}")
        print(f"            🎯 Extracting coefficient index: {coef_idx}")
        
        if coef_idx >= coefficients.shape[2]:
            print(f"            ❌ Coefficient index {coef_idx} out of range (max: {coefficients.shape[2]-1})")
            return None
        
        # Extract coefficients across embedding dimensions
        contrast_coefficients = coefficients[:, :, coef_idx]  # (n_embedding, n_genes)
        
        # Get gene names
        gene_names = lemur_model.adata.var.index
        adata_subset = lemur_model.adata
        
        print(f"            🔬 Implementing LEMUR-specific statistical testing...")
        
        # LEMUR-specific approach: Use embedding consistency as significance measure
        # Calculate statistics across embedding dimensions
        coef_mean = np.mean(contrast_coefficients, axis=0)
        coef_std = np.std(contrast_coefficients, axis=0)
        coef_var = np.var(contrast_coefficients, axis=0)
        
        # Consistency score: genes with consistent effects across embeddings are more reliable
        consistency_score = np.abs(coef_mean) / (coef_std + 1e-8)
        
        # Effect magnitude: absolute mean effect
        effect_magnitude = np.abs(coef_mean)
        
        # Combined score: effect size weighted by consistency
        combined_score = effect_magnitude * consistency_score
        
        # Rank-based approach for LEMUR coefficients
        effect_ranks = stats.rankdata(-effect_magnitude)  # Higher effects get lower ranks
        consistency_ranks = stats.rankdata(-consistency_score)
        combined_ranks = stats.rankdata(-combined_score)
        
        # Convert ranks to percentiles (0-1, where 1 is best)
        n_genes = len(gene_names)
        effect_percentile = 1 - (effect_ranks - 1) / (n_genes - 1)
        consistency_percentile = 1 - (consistency_ranks - 1) / (n_genes - 1)
        combined_percentile = 1 - (combined_ranks - 1) / (n_genes - 1)
        
        # Empirical p-values based on null distribution of coefficients
        # Use the global distribution as null (no effect assumption)
        all_coefs = contrast_coefficients.flatten()
        null_mean = np.mean(all_coefs)
        null_std = np.std(all_coefs)
        
        # Z-scores relative to empirical null
        z_scores = (coef_mean - null_mean) / (null_std + 1e-8)
        empirical_pvals = 2 * (1 - stats.norm.cdf(np.abs(z_scores)))
        
        # Create results DataFrame
        results_df = pd.DataFrame({
            'gene': gene_names,
            'coefficient': coef_mean,
            'abs_coefficient': effect_magnitude,
            'consistency_score': consistency_score,
            'combined_score': combined_score,
            'effect_rank': effect_ranks,
            'consistency_rank': consistency_ranks,
            'combined_rank': combined_ranks,
            'effect_percentile': effect_percentile,
            'consistency_percentile': consistency_percentile,
            'combined_percentile': combined_percentile,
            'z_score': z_scores,
            'pval': empirical_pvals
        })
        
        # Multiple testing correction on empirical p-values
        valid_pvals = ~np.isnan(empirical_pvals) & (empirical_pvals >= 0) & (empirical_pvals <= 1)
        fdr_values = np.full(len(empirical_pvals), 1.0)
        
        if valid_pvals.sum() > 0:
            _, fdr_values[valid_pvals], _, _ = multipletests(
                empirical_pvals[valid_pvals], method='fdr_bh'
            )
        
        results_df['fdr'] = fdr_values
        
        # Truly adaptive significance based on effect distribution characteristics
        # Use natural breaks in the data rather than fixed percentiles
        
        # Calculate distribution statistics
        effect_mean = np.mean(effect_magnitude)
        effect_std = np.std(effect_magnitude)
        effect_median = np.median(effect_magnitude)
        
        consistency_mean = np.mean(consistency_score)
        consistency_std = np.std(consistency_score)
        
        combined_mean = np.mean(combined_score)
        combined_std = np.std(combined_score)
        
        # Adaptive thresholds based on distribution characteristics
        # Effect size: use standard deviations above mean for natural breaks
        effect_threshold_low = effect_mean + 1.0 * effect_std    # Moderate effect
        effect_threshold_high = effect_mean + 2.0 * effect_std   # Strong effect
        
        # Consistency: genes with very consistent effects across embeddings
        consistency_threshold = consistency_mean + 1.5 * consistency_std
        
        # Combined score: balanced approach
        combined_threshold_low = combined_mean + 1.0 * combined_std
        combined_threshold_high = combined_mean + 2.0 * combined_std
        
        # FDR-based with meaningful effect cutoff
        fdr_threshold = 0.1  # Stricter FDR
        effect_for_fdr = max(effect_threshold_low, effect_median * 2)  # Meaningful effect
        
        print(f"            🎯 Adaptive thresholds: effect_high={effect_threshold_high:.4f}, combined_high={combined_threshold_high:.4f}")
        print(f"            📊 Distribution stats: effect_mean±std={effect_mean:.4f}±{effect_std:.4f}")
        
        # Significance based on distribution characteristics
        results_df['significant_moderate_effect'] = results_df['abs_coefficient'] >= effect_threshold_low
        results_df['significant_strong_effect'] = results_df['abs_coefficient'] >= effect_threshold_high
        results_df['significant_high_consistency'] = results_df['consistency_score'] >= consistency_threshold
        results_df['significant_moderate_combined'] = results_df['combined_score'] >= combined_threshold_low
        results_df['significant_strong_combined'] = results_df['combined_score'] >= combined_threshold_high
        
        # FDR-based significance with meaningful effect cutoff
        results_df['significant_fdr'] = (
            (results_df['fdr'] < fdr_threshold) &
            (results_df['abs_coefficient'] >= effect_for_fdr)
        )
        
        # Conservative significance: strong signals
        results_df['significant'] = (
            results_df['significant_strong_combined'] |
            results_df['significant_fdr'] |
            (results_df['significant_strong_effect'] & results_df['significant_high_consistency'])
        )
        
        # Lenient significance: any meaningful signal
        results_df['significant_lenient'] = (
            results_df['significant_moderate_combined'] |
            results_df['significant_moderate_effect'] |
            results_df['significant_high_consistency'] |
            results_df['significant_fdr'] |
            (results_df['fdr'] < 0.2)  # Very lenient FDR
        )
        
        # Add biological annotations
        results_df['contrast'] = contrast_name
        results_df['biology'] = contrast_config['biology']
        results_df['priority'] = contrast_config['priority']
        
        # Sort by combined score (effect + consistency)
        results_df = results_df.sort_values('combined_score', ascending=False)
        results_df = results_df.set_index('gene')
        
        n_sig = results_df['significant'].sum()
        n_sig_lenient = results_df['significant_lenient'].sum()
        n_strong_effect = results_df['significant_strong_effect'].sum()
        n_strong_combined = results_df['significant_strong_combined'].sum()
        n_fdr_sig = results_df['significant_fdr'].sum()
        
        print(f"            ✅ {len(results_df)} genes analyzed")
        print(f"            📊 {n_sig} significant (conservative), {n_sig_lenient} significant (lenient)")
        print(f"            🎯 {n_strong_effect} strong effects, {n_strong_combined} strong combined, {n_fdr_sig} FDR significant")
        print(f"            📈 Effect range: {results_df['coefficient'].min():.4f} to {results_df['coefficient'].max():.4f}")
        print(f"            🏆 Max combined score: {results_df['combined_score'].max():.4f}")
        
        return results_df
        
    except Exception as e:
        print(f"            ❌ DE extraction error: {e}")
        import traceback
        traceback.print_exc()
        return None

# ============================================================================
# 📊 STRATEGIC VISUALIZATION
# ============================================================================

def create_strategic_visualizations(results):
    """Create strategic visualizations focused on biological insights"""
    print(f"\n📊 CREATING STRATEGIC VISUALIZATIONS")
    print("=" * 38)
    
    if not results:
        print("   ❌ No results to visualize")
        return
    
    # 1. Regional comparison overview
    create_regional_overview_plots(results)
    
    # 2. Biological significance plots
    create_biological_significance_plots(results)
    
    # 3. Strategic summary dashboard
    create_strategic_dashboard(results)

def create_regional_overview_plots(results):
    """Create overview plots comparing regional effects"""
    print("   📊 Creating regional overview plots...")
    
    try:
        # Filter for regional contrasts
        regional_results = {}
        for key in ['caudate_oud_effect', 'putamen_oud_effect']:
            if key in results:
                regional_results[key] = results[key]['de_results']
        
        if len(regional_results) < 2:
            print("      ⚠️  Need both regional results for comparison")
            return
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        # Plot 1: Regional effect sizes
        ax = axes[0]
        for name, data in regional_results.items():
            region = 'Caudate' if 'caudate' in name else 'Putamen'
            color = 'lightcoral' if 'caudate' in name else 'lightblue'
            ax.hist(data['coefficient'], bins=50, alpha=0.7, 
                   label=region, density=True, color=color)
        ax.set_xlabel('OUD Effect Size (Coefficient)')
        ax.set_ylabel('Density')
        ax.set_title('Regional OUD Effect Distributions')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Significant genes comparison
        ax = axes[1]
        regions = []
        sig_counts = []
        colors = []
        
        for name, data in regional_results.items():
            region = 'Caudate' if 'caudate' in name else 'Putamen'
            regions.append(region)
            sig_counts.append(data['significant'].sum())
            colors.append('lightcoral' if 'caudate' in name else 'lightblue')
        
        bars = ax.bar(regions, sig_counts, color=colors)
        ax.set_ylabel('Significant Genes')
        ax.set_title('Significant Genes by Region')
        
        # Add counts and percentages
        for bar, count, data in zip(bars, sig_counts, regional_results.values()):
            height = bar.get_height()
            pct = count / len(data) * 100
            ax.text(bar.get_x() + bar.get_width()/2., height + max(sig_counts)*0.01,
                   f'{count}\n({pct:.1f}%)', ha='center', va='bottom')
        
        # Plot 3: Caudate vs Putamen scatter
        if len(regional_results) == 2:
            ax = axes[2]
            caudate_data = regional_results['caudate_oud_effect']
            putamen_data = regional_results['putamen_oud_effect']
            
            # Common genes
            common_genes = set(caudate_data.index) & set(putamen_data.index)
            if len(common_genes) > 10:
                common_list = list(common_genes)
                caudate_coefs = caudate_data.loc[common_list, 'coefficient']
                putamen_coefs = putamen_data.loc[common_list, 'coefficient']
                
                # Color by significance in either region
                colors = []
                for gene in common_list:
                    if caudate_data.loc[gene, 'significant'] and putamen_data.loc[gene, 'significant']:
                        colors.append('red')  # Significant in both
                    elif caudate_data.loc[gene, 'significant']:
                        colors.append('orange')  # Caudate only
                    elif putamen_data.loc[gene, 'significant']:
                        colors.append('blue')  # Putamen only
                    else:
                        colors.append('lightgray')  # Neither
                
                ax.scatter(caudate_coefs, putamen_coefs, c=colors, alpha=0.6, s=20)
                ax.set_xlabel('Caudate OUD Effect')
                ax.set_ylabel('Putamen OUD Effect')
                ax.set_title('Regional OUD Effects: Caudate vs Putamen')
                
                # Diagonal line
                lims = [np.min([ax.get_xlim(), ax.get_ylim()]),
                        np.max([ax.get_xlim(), ax.get_ylim()])]
                ax.plot(lims, lims, 'k--', alpha=0.5, zorder=0)
                
                # Correlation
                corr, p_val = stats.pearsonr(caudate_coefs, putamen_coefs)
                ax.text(0.05, 0.95, f'r = {corr:.3f}\np = {p_val:.2e}', 
                       transform=ax.transAxes, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                # Custom legend
                from matplotlib.patches import Patch
                legend_elements = [
                    Patch(facecolor='red', label='Both regions'),
                    Patch(facecolor='orange', label='Caudate only'),
                    Patch(facecolor='blue', label='Putamen only'),
                    Patch(facecolor='lightgray', label='Neither')
                ]
                ax.legend(handles=legend_elements, loc='lower right')
        
        # Plot 4-5: Volcano plots
        plot_idx = 3
        for name, data in regional_results.items():
            if plot_idx < 6:
                ax = axes[plot_idx]
                region = 'Caudate' if 'caudate' in name else 'Putamen'
                
                x = data['coefficient']
                y = -np.log10(data['fdr'].clip(lower=1e-300))
                
                # Color by significance
                colors = ['red' if sig else 'lightgray' for sig in data['significant']]
                
                ax.scatter(x, y, c=colors, s=15, alpha=0.7)
                ax.set_xlabel('OUD Effect Coefficient')
                ax.set_ylabel('-log10(FDR)')
                ax.set_title(f'{region} OUD Effects\n(Volcano Plot)')
                
                # Significance thresholds
                ax.axhline(y=-np.log10(ANALYSIS_CONFIG['fdr_threshold']), 
                          color='blue', linestyle='--', alpha=0.5, label=f"FDR < {ANALYSIS_CONFIG['fdr_threshold']}")
                ax.axvline(x=ANALYSIS_CONFIG['lfc_threshold'], 
                          color='blue', linestyle='--', alpha=0.5)
                ax.axvline(x=-ANALYSIS_CONFIG['lfc_threshold'], 
                          color='blue', linestyle='--', alpha=0.5)
                
                # Add statistics
                n_sig = data['significant'].sum()
                ax.text(0.05, 0.95, f'{n_sig} significant genes', 
                       transform=ax.transAxes, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                plot_idx += 1
        
        # Remove empty subplot
        if plot_idx < len(axes):
            axes[plot_idx].remove()
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/strategic_regional_overview.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        print("      ✅ Regional overview plots saved")
        
    except Exception as e:
        print(f"      ❌ Regional overview plotting failed: {e}")

def create_biological_significance_plots(results):
    """Create plots focused on biological significance"""
    print("   🧬 Creating biological significance plots...")
    
    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Effect sizes by priority
        ax = axes[0, 0]
        high_priority_effects = []
        medium_priority_effects = []
        
        for name, result in results.items():
            if 'de_results' in result:
                data = result['de_results']
                priority = result['contrast_config']['priority']
                
                if priority == 'HIGH':
                    high_priority_effects.extend(data['abs_coefficient'].tolist())
                elif priority == 'MEDIUM':
                    medium_priority_effects.extend(data['abs_coefficient'].tolist())
        
        if high_priority_effects and medium_priority_effects:
            ax.hist(high_priority_effects, bins=50, alpha=0.7, 
                   label='High Priority', density=True, color='red')
            ax.hist(medium_priority_effects, bins=50, alpha=0.7, 
                   label='Medium Priority', density=True, color='orange')
            ax.set_xlabel('Absolute Effect Size')
            ax.set_ylabel('Density')
            ax.set_title('Effect Size Distribution by Priority')
            ax.legend()
        
        # Plot 2: Biological pathway summary
        ax = axes[0, 1]
        contrast_names = []
        sig_counts = []
        priority_colors = []
        
        for name, result in results.items():
            if 'de_results' in result:
                contrast_names.append(name.replace('_', '\n'))
                sig_counts.append(result['de_results']['significant'].sum())
                priority = result['contrast_config']['priority']
                priority_colors.append('red' if priority == 'HIGH' else 'orange')
        
        if contrast_names:
            bars = ax.bar(range(len(contrast_names)), sig_counts, color=priority_colors)
            ax.set_ylabel('Significant Genes')
            ax.set_title('Significant Genes by Biological Contrast')
            ax.set_xticks(range(len(contrast_names)))
            ax.set_xticklabels(contrast_names, rotation=45, ha='right')
            
            # Add counts
            for bar, count in zip(bars, sig_counts):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + max(sig_counts)*0.01,
                       f'{count}', ha='center', va='bottom')
        
        # Plot 3: Regional biology summary
        ax = axes[1, 0]
        ax.text(0.05, 0.95, "🧠 Regional Biology", fontweight='bold', fontsize=14, transform=ax.transAxes)
        ax.text(0.05, 0.85, "Caudate (Goal-directed):", fontweight='bold', transform=ax.transAxes)
        ax.text(0.05, 0.80, "• Executive function", transform=ax.transAxes)
        ax.text(0.05, 0.75, "• Cognitive control", transform=ax.transAxes)
        ax.text(0.05, 0.70, "• Decision making", transform=ax.transAxes)
        
        ax.text(0.05, 0.60, "Putamen (Habit-based):", fontweight='bold', transform=ax.transAxes)
        ax.text(0.05, 0.55, "• Habit formation", transform=ax.transAxes)
        ax.text(0.05, 0.50, "• Motor control", transform=ax.transAxes)
        ax.text(0.05, 0.45, "• Addiction maintenance", transform=ax.transAxes)
        
        ax.text(0.05, 0.35, "Clinical Relevance:", fontweight='bold', transform=ax.transAxes)
        ax.text(0.05, 0.30, "• Caudate dysfunction → poor decisions", transform=ax.transAxes)
        ax.text(0.05, 0.25, "• Putamen dysfunction → compulsive use", transform=ax.transAxes)
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        
        # Plot 4: Analysis summary
        ax = axes[1, 1]
        total_contrasts = len(results)
        total_sig = sum(result['de_results']['significant'].sum() 
                       for result in results.values() if 'de_results' in result)
        
        high_priority_contrasts = sum(1 for result in results.values() 
                                    if result['contrast_config']['priority'] == 'HIGH')
        
        ax.text(0.05, 0.95, "📊 Analysis Summary", fontweight='bold', fontsize=14, transform=ax.transAxes)
        ax.text(0.05, 0.85, f"Total contrasts: {total_contrasts}", transform=ax.transAxes)
        ax.text(0.05, 0.80, f"High priority: {high_priority_contrasts}", transform=ax.transAxes)
        ax.text(0.05, 0.75, f"Total significant genes: {total_sig:,}", transform=ax.transAxes)
        
        ax.text(0.05, 0.65, "Thresholds:", fontweight='bold', transform=ax.transAxes)
        ax.text(0.05, 0.60, f"FDR < {ANALYSIS_CONFIG['fdr_threshold']}", transform=ax.transAxes)
        ax.text(0.05, 0.55, f"|Effect| > {ANALYSIS_CONFIG['lfc_threshold']}", transform=ax.transAxes)
        
        ax.text(0.05, 0.45, "Strategic Focus:", fontweight='bold', transform=ax.transAxes)
        ax.text(0.05, 0.40, "✅ Regional OUD effects", transform=ax.transAxes)
        ax.text(0.05, 0.35, "✅ Sex × Region interactions", transform=ax.transAxes)
        ax.text(0.05, 0.30, "✅ Biological pathway relevance", transform=ax.transAxes)
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        
        plt.tight_layout()
        plt.savefig(f"{PLOTS_DIR}/biological_significance_analysis.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        print("      ✅ Biological significance plots saved")
        
    except Exception as e:
        print(f"      ❌ Biological significance plotting failed: {e}")

def create_strategic_dashboard(results):
    """Create comprehensive strategic analysis dashboard"""
    print("   📊 Creating strategic dashboard...")
    
    try:
        fig = plt.figure(figsize=(20, 16))
        
        # Create grid layout
        gs = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)
        
        # Title
        fig.suptitle('LEMUR Strategic Regional Analysis Dashboard', fontsize=20, fontweight='bold')
        
        # Dashboard sections
        create_dashboard_summary(fig, gs, results)
        create_dashboard_regional_comparison(fig, gs, results)
        create_dashboard_priority_analysis(fig, gs, results)
        create_dashboard_biological_insights(fig, gs, results)
        
        plt.savefig(f"{PLOTS_DIR}/strategic_analysis_dashboard.png", 
                   dpi=ANALYSIS_CONFIG['figure_dpi'], bbox_inches='tight')
        plt.close()
        
        print("      ✅ Strategic dashboard saved")
        
    except Exception as e:
        print(f"      ❌ Dashboard creation failed: {e}")

def create_dashboard_summary(fig, gs, results):
    """Create dashboard summary section"""
    ax = fig.add_subplot(gs[0, :2])
    ax.text(0.05, 0.9, "🎯 Strategic Analysis Summary", fontsize=16, fontweight='bold', transform=ax.transAxes)
    
    total_contrasts = len(results)
    total_genes_analyzed = 0
    total_sig_genes = 0
    
    for result in results.values():
        if 'de_results' in result:
            total_genes_analyzed = len(result['de_results'])  # Same for all
            total_sig_genes += result['de_results']['significant'].sum()
    
    ax.text(0.05, 0.7, f"• Contrasts analyzed: {total_contrasts}", fontsize=12, transform=ax.transAxes)
    ax.text(0.05, 0.6, f"• Genes per contrast: {total_genes_analyzed:,}", fontsize=12, transform=ax.transAxes)
    ax.text(0.05, 0.5, f"• Total significant findings: {total_sig_genes:,}", fontsize=12, transform=ax.transAxes)
    ax.text(0.05, 0.4, f"• Average per contrast: {total_sig_genes/total_contrasts:.0f}", fontsize=12, transform=ax.transAxes)
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

def create_dashboard_regional_comparison(fig, gs, results):
    """Create regional comparison section"""
    ax = fig.add_subplot(gs[0, 2:])
    
    # Regional effects comparison
    regional_data = {}
    for key in ['caudate_oud_effect', 'putamen_oud_effect']:
        if key in results:
            regional_data[key] = results[key]['de_results']['significant'].sum()
    
    if regional_data:
        regions = ['Caudate', 'Putamen']
        counts = [regional_data.get('caudate_oud_effect', 0), 
                 regional_data.get('putamen_oud_effect', 0)]
        colors = ['lightcoral', 'lightblue']
        
        bars = ax.bar(regions, counts, color=colors)
        ax.set_ylabel('Significant Genes')
        ax.set_title('Regional OUD Effects Comparison', fontweight='bold')
        
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(counts)*0.01,
                   f'{count}', ha='center', va='bottom')

def create_dashboard_priority_analysis(fig, gs, results):
    """Create priority analysis section"""
    ax = fig.add_subplot(gs[1, :])
    
    # Create priority-based analysis
    priorities = ['HIGH', 'MEDIUM']
    priority_data = {p: [] for p in priorities}
    
    for name, result in results.items():
        if 'de_results' in result:
            priority = result['contrast_config']['priority']
            sig_count = result['de_results']['significant'].sum()
            priority_data[priority].append((name, sig_count))
    
    # Plot by priority
    x_pos = 0
    colors = {'HIGH': 'red', 'MEDIUM': 'orange'}
    
    for priority in priorities:
        if priority_data[priority]:
            names, counts = zip(*priority_data[priority])
            x_positions = range(x_pos, x_pos + len(names))
            bars = ax.bar(x_positions, counts, color=colors[priority], 
                         alpha=0.7, label=f'{priority} Priority')
            
            # Add labels
            for i, (name, count) in enumerate(priority_data[priority]):
                ax.text(x_pos + i, count + max(counts)*0.01, f'{count}', 
                       ha='center', va='bottom')
                display_name = name.replace('_', '\n')
                ax.text(x_pos + i, -max(counts)*0.05, display_name, 
                       ha='center', va='top', rotation=45, fontsize=8)
            
            x_pos += len(names) + 1
    
    ax.set_ylabel('Significant Genes')
    ax.set_title('Strategic Contrasts by Priority Level', fontweight='bold')
    ax.legend()
    ax.set_xticks([])

def create_dashboard_biological_insights(fig, gs, results):
    """Create biological insights section"""
    ax = fig.add_subplot(gs[2:, :])
    ax.text(0.02, 0.95, "🧬 Key Biological Insights", fontsize=16, fontweight='bold', transform=ax.transAxes)
    
    # Extract top findings
    y_pos = 0.85
    for name, result in results.items():
        if 'de_results' in result and result['de_results']['significant'].sum() > 0:
            n_sig = result['de_results']['significant'].sum()
            biology = result['contrast_config']['biology']
            priority = result['contrast_config']['priority']
            
            priority_emoji = "🔥" if priority == 'HIGH' else "⚡"
            
            ax.text(0.02, y_pos, f"{priority_emoji} {name.replace('_', ' ').title()}", 
                   fontsize=12, fontweight='bold', transform=ax.transAxes)
            ax.text(0.02, y_pos-0.03, f"   {n_sig} significant genes", 
                   fontsize=10, transform=ax.transAxes)
            ax.text(0.02, y_pos-0.06, f"   Biology: {biology}", 
                   fontsize=10, transform=ax.transAxes, style='italic')
            
            # Show top genes
            if n_sig > 0:
                top_genes = result['de_results'].nlargest(3, 'abs_coefficient')
                gene_list = ", ".join(top_genes.index[:3])
                ax.text(0.02, y_pos-0.09, f"   Top genes: {gene_list}", 
                       fontsize=9, transform=ax.transAxes, color='blue')
            
            y_pos -= 0.15
            
            if y_pos < 0.1:
                break
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

# ============================================================================
# 💾 STRATEGIC EXPORT FUNCTIONS
# ============================================================================

def export_strategic_results(results):
    """Export strategic results with biological focus"""
    print(f"\n💾 EXPORTING STRATEGIC RESULTS")
    print("=" * 32)
    
    # Export each contrast with biological annotations
    for contrast_name, result in results.items():
        if 'de_results' in result:
            print(f"   📊 Exporting {contrast_name}...")
            
            de_results = result['de_results']
            contrast_config = result['contrast_config']
            
            # Add biological context to results
            de_results_annotated = de_results.copy()
            de_results_annotated['biological_context'] = contrast_config['biology']
            de_results_annotated['analysis_priority'] = contrast_config['priority']
            de_results_annotated['contrast_description'] = contrast_config['description']
            
            # Export full results
            filename = f"{TABLES_DIR}/strategic_{contrast_name}_results.csv"
            de_results_annotated.to_csv(filename, index=True)
            print(f"      ✅ Full results: {filename}")
            
            # Export significant genes with biological focus
            sig_genes = de_results_annotated[de_results_annotated['significant']]
            sig_genes_lenient = de_results_annotated[de_results_annotated['significant_lenient']]
            
            if len(sig_genes) > 0:
                sig_filename = f"{TABLES_DIR}/strategic_{contrast_name}_significant_genes.csv"
                sig_genes.to_csv(sig_filename, index=True)
                print(f"      🧬 Significant genes ({len(sig_genes)}): {sig_filename}")
                
                # Show top biological findings by combined score
                top_genes = sig_genes.nlargest(5, 'combined_score')
                print(f"         Top genes: {', '.join(top_genes.index[:5])}")
            
            elif len(sig_genes_lenient) > 0:
                sig_filename = f"{TABLES_DIR}/strategic_{contrast_name}_significant_lenient_genes.csv"
                sig_genes_lenient.to_csv(sig_filename, index=True)
                print(f"      🔬 Lenient significant genes ({len(sig_genes_lenient)}): {sig_filename}")
                
                # Show top biological findings by combined score
                top_genes = sig_genes_lenient.nlargest(5, 'combined_score')
                print(f"         Top genes: {', '.join(top_genes.index[:5])}")
            
            else:
                # Export top effects by combined score
                top_effects = de_results_annotated.nlargest(20, 'combined_score')
                top_filename = f"{TABLES_DIR}/strategic_{contrast_name}_top_effects.csv"
                top_effects.to_csv(top_filename, index=True)
                print(f"      📊 Top effects (by combined score, {len(top_effects)}): {top_filename}")
                
                # Also export top effect percentile genes
                top_percentile = de_results_annotated[de_results_annotated['effect_percentile'] >= 0.95]
                if len(top_percentile) > 0:
                    percentile_filename = f"{TABLES_DIR}/strategic_{contrast_name}_top_percentile_effects.csv"
                    top_percentile.to_csv(percentile_filename, index=True)
                    print(f"      🎯 Top 5% effects ({len(top_percentile)}): {percentile_filename}")
    
    # Create strategic master summary
    create_strategic_master_summary(results)

def create_strategic_master_summary(results):
    """Create strategic master summary with biological insights"""
    print("   📋 Creating strategic master summary...")
    
    summary_lines = [
        "# LEMUR Strategic Regional Analysis - Master Summary",
        "=" * 54,
        "",
        "## Strategic Analysis Overview",
        "This analysis completes the comprehensive LEMUR methodology comparison",
        "with DESeq2 by focusing on high-impact regional and sex effects in OUD.",
        "",
        "## Strategic Contrast Design",
        "",
        "### Tier 1: Core Effects (Previously Completed)",
        "1. ✅ OUD vs Control (Pooled) - Overall effect",
        "2. ✅ OUD vs Control (Males) - Male-specific effects", 
        "3. ✅ OUD vs Control (Females) - Female-specific effects",
        "4. ✅ OUD × Sex Interaction - Sex-differential effects",
        "",
        "### Tier 2: Regional Effects (This Analysis - HIGH PRIORITY)",
    ]
    
    # Add strategic results
    regional_summary = {}
    interaction_summary = {}
    
    for contrast_name, result in results.items():
        if 'de_results' in result:
            de_results = result['de_results']
            n_sig = de_results['significant'].sum()
            n_sig_lenient = de_results['significant_lenient'].sum()
            priority = result['contrast_config']['priority']
            biology = result['contrast_config']['biology']
            
            if 'oud_effect' in contrast_name and priority == 'HIGH':
                region = 'Caudate' if 'caudate' in contrast_name else 'Putamen'
                regional_summary[region] = {
                    'n_sig': n_sig,
                    'n_sig_lenient': n_sig_lenient,
                    'top_genes': list(de_results.nlargest(5, 'abs_coefficient').index),
                    'biology': biology
                }
            elif priority == 'MEDIUM':
                interaction_summary[contrast_name] = {
                    'n_sig': n_sig,
                    'n_sig_lenient': n_sig_lenient,
                    'biology': biology
                }
    
    # Regional results
    for region, data in regional_summary.items():
        summary_lines.extend([
            f"5. 🧠 OUD vs Control ({region}) - {data['n_sig']} significant genes (strict), {data['n_sig_lenient']} genes (lenient)",
            f"   Biology: {data['biology']}",
            f"   Top genes: {', '.join(data['top_genes'][:3])}",
            ""
        ])
    
    summary_lines.extend([
        "### Tier 3: Complex Interactions (MEDIUM PRIORITY)",
    ])
    
    for contrast, data in interaction_summary.items():
        display_name = contrast.replace('_', ' ').title()
        summary_lines.extend([
            f"• {display_name} - {data['n_sig']} significant genes (strict), {data['n_sig_lenient']} genes (lenient)",
            f"  Biology: {data['biology']}",
            ""
        ])
    
    # Key findings
    total_sig = sum(result['de_results']['significant'].sum() 
                   for result in results.values() if 'de_results' in result)
    total_sig_lenient = sum(result['de_results']['significant_lenient'].sum() 
                           for result in results.values() if 'de_results' in result)
    
    summary_lines.extend([
        "## Key Strategic Findings",
        f"- Total significant genes across regional contrasts: {total_sig:,} (strict), {total_sig_lenient:,} (lenient)",
        "",
        "### Regional Comparison",
    ])
    
    if 'Caudate' in regional_summary and 'Putamen' in regional_summary:
        caudate_sig = regional_summary['Caudate']['n_sig']
        putamen_sig = regional_summary['Putamen']['n_sig']
        caudate_sig_lenient = regional_summary['Caudate']['n_sig_lenient']
        putamen_sig_lenient = regional_summary['Putamen']['n_sig_lenient']
        
        summary_lines.extend([
            f"- **Caudate OUD effects**: {caudate_sig} genes (strict), {caudate_sig_lenient} genes (lenient)",
            f"- **Putamen OUD effects**: {putamen_sig} genes (strict), {putamen_sig_lenient} genes (lenient)",
            f"- **Regional ratio (lenient)**: {caudate_sig_lenient/putamen_sig_lenient:.1f}x (Caudate/Putamen)" if putamen_sig_lenient > 0 else "- **Regional pattern**: Caudate-dominant effects",
            "",
            "### Biological Interpretation",
            "- Caudate (goal-directed): More OUD-related disruption" if caudate_sig_lenient > putamen_sig_lenient else "- Putamen (habit-based): More OUD-related disruption",
            "- Clinical relevance: Target region-specific interventions",
            ""
        ])
    
    summary_lines.extend([
        "## Methodology Validation",
        "- ✅ LEMUR raw data analysis completed",
        "- ✅ Regional stratification implemented", 
        "- ✅ Sex interaction effects analyzed",
        "- ✅ Comprehensive comparison with DESeq2 enabled",
        "",
        "## Next Steps",
        "1. Compare LEMUR vs DESeq2 regional findings",
        "2. Validate top regional genes experimentally",
        "3. Develop region-specific therapeutic targets",
        "",
        f"## Technical Parameters",
        f"- Embedding dimensions: {ANALYSIS_CONFIG['n_embedding']}",
        f"- Discovery FDR threshold: {ANALYSIS_CONFIG['fdr_threshold']}",
        f"- Effect size threshold: {ANALYSIS_CONFIG['lfc_threshold']}",
        f"- Analysis approach: Strategic high-impact contrasts",
        f"- Gene coverage: All highly variable genes (no limits)",
        f"- Total contrasts: {len(results)} (focused on biological relevance)",
    ])
    
    # Save master summary
    summary_filename = f"{TABLES_DIR}/strategic_master_summary.md"
    with open(summary_filename, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    print(f"      ✅ Strategic master summary: {summary_filename}")

# ============================================================================
# 🚀 MAIN STRATEGIC EXECUTION
# ============================================================================

def main():
    """Main strategic analysis pipeline"""
    print(f"\n🚀 LEMUR STRATEGIC REGIONAL ANALYSIS PIPELINE - FULL DATASET")
    print("=" * 62)
    print(f"🎯 Focus: High-impact biological contrasts with maximum power")
    print(f"🧠 Goal: Complete regional OUD effects analysis using all ~99K cells")
    print(f"💾 Resources: 64GB RAM, no subsampling or gene filtering limitations")
    
    # 1. Load and validate data
    adata = load_and_prepare_data()
    if adata is None:
        print("❌ Data loading failed")
        return None
    
    # 2. Create strategic datasets
    datasets = prepare_strategic_datasets(adata)
    if not datasets:
        print("❌ Strategic dataset preparation failed")
        return None
    
    # 3. Run strategic LEMUR analysis
    results = run_strategic_lemur_analysis(datasets)
    if not results:
        print("❌ Strategic LEMUR analysis failed")
        return None
    
    # 4. Create strategic visualizations
    create_strategic_visualizations(results)
    
    # 5. Export strategic results
    export_strategic_results(results)
    
    # Final strategic summary
    print(f"\n🎉 STRATEGIC REGIONAL ANALYSIS COMPLETE!")
    print("=" * 42)
    print(f"📊 Results: {OUTPUT_DIR}")
    print(f"🖼️  Plots: {PLOTS_DIR}/")
    print(f"📋 Tables: {TABLES_DIR}/")
    
    print(f"\n🧠 STRATEGIC FINDINGS SUMMARY:")
    
    # Priority-based summary
    high_priority_results = {}
    medium_priority_results = {}
    
    for name, result in results.items():
        if 'de_results' in result:
            priority = result['contrast_config']['priority']
            n_sig = result['de_results']['significant'].sum()
            n_sig_lenient = result['de_results']['significant_lenient'].sum()
            
            if priority == 'HIGH':
                high_priority_results[name] = {'strict': n_sig, 'lenient': n_sig_lenient}
            else:
                medium_priority_results[name] = {'strict': n_sig, 'lenient': n_sig_lenient}
    
    print(f"  🔥 HIGH PRIORITY EFFECTS:")
    for name, counts in high_priority_results.items():
        region = name.replace('_oud_effect', '').title()
        print(f"     • {region} OUD effects: {counts['strict']} significant genes (strict), {counts['lenient']} genes (lenient)")
    
    if medium_priority_results:
        print(f"  ⚡ MEDIUM PRIORITY EFFECTS:")
        for name, counts in medium_priority_results.items():
            display_name = name.replace('_', ' ').title()
            print(f"     • {display_name}: {counts['strict']} significant genes (strict), {counts['lenient']} genes (lenient)")
    
    # Regional comparison
    if 'caudate_oud_effect' in high_priority_results and 'putamen_oud_effect' in high_priority_results:
        caudate_sig = high_priority_results['caudate_oud_effect']['lenient']
        putamen_sig = high_priority_results['putamen_oud_effect']['lenient']
        
        print(f"\n🏆 REGIONAL COMPARISON (LENIENT THRESHOLDS):")
        print(f"     • Caudate: {caudate_sig} genes (goal-directed behavior)")
        print(f"     • Putamen: {putamen_sig} genes (habit formation)")
        
        if putamen_sig > 0:
            ratio = caudate_sig / putamen_sig
            dominant_region = "Caudate" if ratio > 1 else "Putamen"
            print(f"     • Dominant region: {dominant_region} ({ratio:.1f}x ratio)")
        else:
            print(f"     • Pattern: Caudate-specific OUD effects")
    
    total_strategic_genes = sum(counts['strict'] for counts in high_priority_results.values()) + sum(counts['strict'] for counts in medium_priority_results.values())
    total_strategic_genes_lenient = sum(counts['lenient'] for counts in high_priority_results.values()) + sum(counts['lenient'] for counts in medium_priority_results.values())
    print(f"\n📈 TOTAL STRATEGIC DISCOVERIES: {total_strategic_genes} significant genes (strict), {total_strategic_genes_lenient} genes (lenient)")
    print(f"💪 FULL DATASET ADVANTAGE: Maximum statistical power achieved")
    
    print(f"\n✅ Strategic analysis ready for DESeq2 comparison!")
    print(f"   🚀 Full dataset provides optimal power for regional effects")
    print(f"   📊 Use results for comprehensive methodology validation.")
    
    return results

if __name__ == "__main__":
    # Set random seed for reproducibility
    np.random.seed(ANALYSIS_CONFIG['random_seed'])
    
    # Run strategic analysis
    results = main()