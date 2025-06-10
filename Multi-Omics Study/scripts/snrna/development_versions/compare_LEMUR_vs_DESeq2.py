#!/usr/bin/env python3
"""
🔬 LEMUR vs pyDESeq2 Comparison Analysis
Compare differential expression results between LEMUR and pyDESeq2 methods

This script compares:
1. Gene overlap between methods
2. Effect size correlations
3. Direction concordance
4. Statistical significance agreement
5. Top gene rankings

Author: Research Team
Date: 2024
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import warnings

warnings.filterwarnings('ignore')

# Set up plotting
plt.style.use('default')
sns.set_palette("husl")

# ============================================================================
# 📁 CONFIGURATION
# ============================================================================

BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"

# Input files
LEMUR_MAIN_EFFECT = f"{BASE_DIR}/results/snrna_scvi/lemur_final_working/tables/lemur_main_effect_results.csv"
LEMUR_SEX_EFFECT = f"{BASE_DIR}/results/snrna_scvi/lemur_final_working/tables/lemur_sex_effect_results.csv"
DESEQ2_MAIN_EFFECT = f"{BASE_DIR}/results/snrna_scvi/pydeseq2_analysis/tables/Full_Dataset/01_OUD_vs_Control_Pooled_results.csv"
DESEQ2_SEX_EFFECT = f"{BASE_DIR}/results/snrna_scvi/pydeseq2_analysis/tables/Full_Dataset/06_OUD_Effect_Male_vs_Female_results.csv"

# Output directory
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/method_comparison"
PLOTS_DIR = f"{OUTPUT_DIR}/plots"
TABLES_DIR = f"{OUTPUT_DIR}/tables"

# Create output directories
for directory in [OUTPUT_DIR, PLOTS_DIR, TABLES_DIR]:
    os.makedirs(directory, exist_ok=True)

print("🔬 LEMUR vs pyDESeq2 COMPARISON")
print("===============================")
print(f"Output: {OUTPUT_DIR}")

# ============================================================================
# 📊 DATA LOADING
# ============================================================================

def load_comparison_data():
    """Load both LEMUR and DESeq2 results"""
    print(f"\n📁 LOADING COMPARISON DATA")
    print("=" * 28)
    
    # Load LEMUR results
    print("   Loading LEMUR results...")
    try:
        lemur_main = pd.read_csv(LEMUR_MAIN_EFFECT, index_col=0)
        lemur_sex = pd.read_csv(LEMUR_SEX_EFFECT, index_col=0)
        print(f"      ✅ LEMUR main effect: {len(lemur_main)} genes")
        print(f"      ✅ LEMUR sex effect: {len(lemur_sex)} genes")
    except Exception as e:
        print(f"      ❌ LEMUR loading failed: {e}")
        return None, None, None, None
    
    # Load DESeq2 results
    print("   Loading DESeq2 results...")
    try:
        deseq2_main = pd.read_csv(DESEQ2_MAIN_EFFECT)
        deseq2_sex = pd.read_csv(DESEQ2_SEX_EFFECT)
        print(f"      ✅ DESeq2 main effect: {len(deseq2_main)} genes")
        print(f"      ✅ DESeq2 sex effect: {len(deseq2_sex)} genes")
    except Exception as e:
        print(f"      ❌ DESeq2 loading failed: {e}")
        return None, None, None, None
    
    return lemur_main, lemur_sex, deseq2_main, deseq2_sex

# ============================================================================
# 🔍 COMPARISON FUNCTIONS
# ============================================================================

def compare_main_effects(lemur_main, deseq2_main):
    """Compare main effect (OUD vs Control) between methods"""
    print(f"\n🎯 COMPARING MAIN EFFECTS (OUD vs Control)")
    print("=" * 42)
    
    # Merge datasets on gene names
    # DESeq2 has 'gene' column, LEMUR has gene names as index
    deseq2_main_clean = deseq2_main.copy()
    deseq2_main_clean = deseq2_main_clean.set_index('gene')
    
    # Find common genes
    common_genes = set(lemur_main.index) & set(deseq2_main_clean.index)
    common_genes_list = list(common_genes)
    print(f"   Common genes: {len(common_genes)}")
    print(f"   LEMUR only: {len(set(lemur_main.index) - common_genes)}")
    print(f"   DESeq2 only: {len(set(deseq2_main_clean.index) - common_genes)}")
    
    if len(common_genes) == 0:
        print("   ❌ No common genes found!")
        return None
    
    # Create merged dataset
    merged = pd.DataFrame(index=common_genes_list)
    merged['lemur_coef'] = lemur_main.loc[common_genes_list, 'coefficient']
    merged['lemur_fdr'] = lemur_main.loc[common_genes_list, 'fdr']
    merged['lemur_significant'] = lemur_main.loc[common_genes_list, 'significant']
    
    merged['deseq2_lfc'] = deseq2_main_clean.loc[common_genes_list, 'log2FoldChange']
    merged['deseq2_padj'] = deseq2_main_clean.loc[common_genes_list, 'padj']
    merged['deseq2_significant'] = deseq2_main_clean.loc[common_genes_list, 'significant']
    
    # Handle missing values
    merged = merged.dropna()
    print(f"   Final comparison dataset: {len(merged)} genes")
    
    # Calculate correlations
    print(f"\n   📊 EFFECT SIZE CORRELATIONS")
    print("   " + "=" * 28)
    
    pearson_r, pearson_p = pearsonr(merged['lemur_coef'], merged['deseq2_lfc'])
    spearman_r, spearman_p = spearmanr(merged['lemur_coef'], merged['deseq2_lfc'])
    
    print(f"   Pearson correlation: r = {pearson_r:.3f}, p = {pearson_p:.2e}")
    print(f"   Spearman correlation: ρ = {spearman_r:.3f}, p = {spearman_p:.2e}")
    
    # Direction concordance
    print(f"\n   📈 DIRECTION CONCORDANCE")
    print("   " + "=" * 24)
    
    lemur_direction = np.sign(merged['lemur_coef'])
    deseq2_direction = np.sign(merged['deseq2_lfc'])
    direction_agreement = (lemur_direction == deseq2_direction).mean()
    
    print(f"   Direction agreement: {direction_agreement:.1%}")
    
    # Significance overlap
    print(f"\n   🎯 SIGNIFICANCE OVERLAP")
    print("   " + "=" * 24)
    
    lemur_sig_count = merged['lemur_significant'].sum()
    deseq2_sig_count = merged['deseq2_significant'].sum()
    both_sig_count = (merged['lemur_significant'] & merged['deseq2_significant']).sum()
    
    print(f"   LEMUR significant: {lemur_sig_count}")
    print(f"   DESeq2 significant: {deseq2_sig_count}")
    print(f"   Both significant: {both_sig_count}")
    
    if deseq2_sig_count > 0:
        overlap_pct = both_sig_count / deseq2_sig_count * 100
        print(f"   Overlap (% of DESeq2): {overlap_pct:.1f}%")
    
    return merged

def compare_top_genes(lemur_main, deseq2_main, top_n=50):
    """Compare top genes between methods"""
    print(f"\n🔝 TOP {top_n} GENES COMPARISON")
    print("=" * 30)
    
    # Prepare DESeq2 data
    deseq2_main_clean = deseq2_main.copy()
    deseq2_main_clean = deseq2_main_clean.set_index('gene')
    
    # Get top genes by effect size magnitude
    lemur_top = lemur_main.nlargest(top_n, 'abs_coefficient')
    deseq2_top = deseq2_main_clean.nlargest(top_n, 'baseMean')  # or could use abs(log2FoldChange)
    
    # Alternative: top by significance
    deseq2_sig = deseq2_main_clean[deseq2_main_clean['significant'] == True]
    if len(deseq2_sig) > 0:
        deseq2_top_sig = deseq2_sig.nsmallest(min(top_n, len(deseq2_sig)), 'padj')
    else:
        deseq2_top_sig = deseq2_main_clean.nsmallest(top_n, 'pvalue')
    
    print(f"   LEMUR top {top_n} genes (by effect size):")
    for i, (gene, row) in enumerate(lemur_top.head(10).iterrows(), 1):
        print(f"      {i:2d}. {gene}: coef={row['coefficient']:.3f}, FDR={row['fdr']:.2e}")
    
    print(f"\n   DESeq2 top 10 significant genes:")
    for i, (gene, row) in enumerate(deseq2_top_sig.head(10).iterrows(), 1):
        print(f"      {i:2d}. {gene}: LFC={row['log2FoldChange']:.3f}, padj={row['padj']:.2e}")
    
    # Find overlaps
    lemur_top_genes = set(lemur_top.index)
    deseq2_top_genes = set(deseq2_top_sig.index)
    overlap_genes = lemur_top_genes & deseq2_top_genes
    
    print(f"\n   📊 TOP GENE OVERLAP")
    print("   " + "=" * 20)
    print(f"   LEMUR top {top_n}: {len(lemur_top_genes)} genes")
    print(f"   DESeq2 significant: {len(deseq2_top_genes)} genes")
    print(f"   Overlap: {len(overlap_genes)} genes")
    
    if len(overlap_genes) > 0:
        print(f"   Overlapping genes: {list(overlap_genes)[:10]}")  # Show first 10
    
    return lemur_top, deseq2_top_sig, overlap_genes

def analyze_specific_genes(lemur_main, deseq2_main):
    """Analyze specific genes of interest"""
    print(f"\n🧬 SPECIFIC GENE ANALYSIS")
    print("=" * 25)
    
    # Genes of interest from LEMUR analysis
    genes_of_interest = [
        'ST6GALNAC3', 'LRP1B', 'MTRNR2L12', 'PDE1A', 'DPP10',
        'LMCD1-AS1', 'LINC00882', 'ESRRG', 'MOBP', 'ANKRD44'
    ]
    
    # Prepare DESeq2 data
    deseq2_main_clean = deseq2_main.copy()
    deseq2_main_clean = deseq2_main_clean.set_index('gene')
    
    comparison_results = []
    
    for gene in genes_of_interest:
        if gene in lemur_main.index and gene in deseq2_main_clean.index:
            lemur_data = lemur_main.loc[gene]
            deseq2_data = deseq2_main_clean.loc[gene]
            
            comparison_results.append({
                'gene': gene,
                'lemur_coef': lemur_data['coefficient'],
                'lemur_fdr': lemur_data['fdr'],
                'lemur_sig': lemur_data['significant'],
                'deseq2_lfc': deseq2_data['log2FoldChange'],
                'deseq2_padj': deseq2_data['padj'],
                'deseq2_sig': deseq2_data['significant'],
                'direction_agree': np.sign(lemur_data['coefficient']) == np.sign(deseq2_data['log2FoldChange'])
            })
    
    if comparison_results:
        comparison_df = pd.DataFrame(comparison_results)
        print(f"   Found {len(comparison_df)} genes in both datasets:")
        
        for _, row in comparison_df.iterrows():
            direction_str = "✅" if row['direction_agree'] else "❌"
            lemur_sig_str = "***" if row['lemur_sig'] else ""
            deseq2_sig_str = "***" if row['deseq2_sig'] else ""
            
            print(f"   {row['gene']:12} | LEMUR: {row['lemur_coef']:6.3f} {lemur_sig_str:3} | "
                  f"DESeq2: {row['deseq2_lfc']:6.3f} {deseq2_sig_str:3} | {direction_str}")
        
        return comparison_df
    else:
        print("   ❌ No genes of interest found in both datasets")
        return None

# ============================================================================
# 📊 VISUALIZATION FUNCTIONS
# ============================================================================

def create_comparison_plots(merged_data, comparison_df=None):
    """Create comprehensive comparison visualizations"""
    print(f"\n📊 CREATING COMPARISON PLOTS")
    print("=" * 30)
    
    if merged_data is None:
        print("   ❌ No merged data available for plotting")
        return
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    
    # Plot 1: Effect size correlation
    ax1 = plt.subplot(2, 3, 1)
    scatter = ax1.scatter(merged_data['lemur_coef'], merged_data['deseq2_lfc'], 
                         c='blue', alpha=0.6, s=20)
    ax1.set_xlabel('LEMUR Coefficient')
    ax1.set_ylabel('DESeq2 Log2 Fold Change')
    ax1.set_title('Effect Size Correlation')
    
    # Add correlation line
    z = np.polyfit(merged_data['lemur_coef'], merged_data['deseq2_lfc'], 1)
    p = np.poly1d(z)
    ax1.plot(merged_data['lemur_coef'], p(merged_data['lemur_coef']), "r--", alpha=0.8)
    
    # Add correlation coefficient
    r, p_val = pearsonr(merged_data['lemur_coef'], merged_data['deseq2_lfc'])
    ax1.text(0.05, 0.95, f'r = {r:.3f}\np = {p_val:.2e}', 
             transform=ax1.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Plot 2: -log10(p-value) comparison
    ax2 = plt.subplot(2, 3, 2)
    lemur_log_p = -np.log10(merged_data['lemur_fdr'].clip(lower=1e-300))
    deseq2_log_p = -np.log10(merged_data['deseq2_padj'].clip(lower=1e-300))
    
    ax2.scatter(lemur_log_p, deseq2_log_p, c='green', alpha=0.6, s=20)
    ax2.set_xlabel('-log10(LEMUR FDR)')
    ax2.set_ylabel('-log10(DESeq2 padj)')
    ax2.set_title('Significance Correlation')
    
    # Plot 3: Significance agreement
    ax3 = plt.subplot(2, 3, 3)
    sig_categories = ['Both NS', 'LEMUR only', 'DESeq2 only', 'Both Sig']
    sig_counts = [
        (~merged_data['lemur_significant'] & ~merged_data['deseq2_significant']).sum(),
        (merged_data['lemur_significant'] & ~merged_data['deseq2_significant']).sum(),
        (~merged_data['lemur_significant'] & merged_data['deseq2_significant']).sum(),
        (merged_data['lemur_significant'] & merged_data['deseq2_significant']).sum()
    ]
    
    bars = ax3.bar(sig_categories, sig_counts, color=['gray', 'blue', 'red', 'purple'])
    ax3.set_ylabel('Number of Genes')
    ax3.set_title('Significance Agreement')
    plt.setp(ax3.get_xticklabels(), rotation=45, ha='right')
    
    # Add counts on bars
    for bar, count in zip(bars, sig_counts):
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height + max(sig_counts)*0.01,
                f'{count}', ha='center', va='bottom')
    
    # Plot 4: Effect size distribution comparison
    ax4 = plt.subplot(2, 3, 4)
    ax4.hist(merged_data['lemur_coef'], bins=50, alpha=0.5, label='LEMUR', density=True)
    ax4.hist(merged_data['deseq2_lfc'], bins=50, alpha=0.5, label='DESeq2', density=True)
    ax4.set_xlabel('Effect Size')
    ax4.set_ylabel('Density')
    ax4.set_title('Effect Size Distributions')
    ax4.legend()
    
    # Plot 5: Rank correlation
    ax5 = plt.subplot(2, 3, 5)
    lemur_ranks = merged_data['lemur_coef'].abs().rank(ascending=False)
    deseq2_ranks = merged_data['deseq2_lfc'].abs().rank(ascending=False)
    
    ax5.scatter(lemur_ranks, deseq2_ranks, c='orange', alpha=0.6, s=20)
    ax5.set_xlabel('LEMUR Rank (by |coefficient|)')
    ax5.set_ylabel('DESeq2 Rank (by |LFC|)')
    ax5.set_title('Gene Ranking Correlation')
    
    # Add Spearman correlation
    rho, p_val = spearmanr(lemur_ranks, deseq2_ranks)
    ax5.text(0.05, 0.95, f'ρ = {rho:.3f}\np = {p_val:.2e}', 
             transform=ax5.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Plot 6: Direction agreement
    ax6 = plt.subplot(2, 3, 6)
    lemur_sign = np.sign(merged_data['lemur_coef'])
    deseq2_sign = np.sign(merged_data['deseq2_lfc'])
    
    direction_agreement = (lemur_sign == deseq2_sign)
    agree_counts = [direction_agreement.sum(), (~direction_agreement).sum()]
    
    ax6.pie(agree_counts, labels=['Same Direction', 'Opposite Direction'], 
            autopct='%1.1f%%', colors=['lightgreen', 'salmon'])
    ax6.set_title('Direction Agreement')
    
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/lemur_vs_deseq2_comparison.png", 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print("   ✅ Comparison plots saved")

def create_gene_comparison_plot(comparison_df):
    """Create specific gene comparison plot"""
    if comparison_df is None:
        return
    
    print("   Creating gene-specific comparison...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Effect size comparison
    x_pos = np.arange(len(comparison_df))
    width = 0.35
    
    ax1.bar(x_pos - width/2, comparison_df['lemur_coef'], width, 
            label='LEMUR Coefficient', alpha=0.7)
    ax1.bar(x_pos + width/2, comparison_df['deseq2_lfc'], width, 
            label='DESeq2 Log2FC', alpha=0.7)
    
    ax1.set_xlabel('Genes')
    ax1.set_ylabel('Effect Size')
    ax1.set_title('Gene-Specific Effect Size Comparison')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(comparison_df['gene'], rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Significance comparison
    lemur_sig = comparison_df['lemur_sig'].astype(int)
    deseq2_sig = comparison_df['deseq2_sig'].astype(int)
    
    ax2.bar(x_pos - width/2, lemur_sig, width, 
            label='LEMUR Significant', alpha=0.7)
    ax2.bar(x_pos + width/2, deseq2_sig, width, 
            label='DESeq2 Significant', alpha=0.7)
    
    ax2.set_xlabel('Genes')
    ax2.set_ylabel('Significant (0/1)')
    ax2.set_title('Gene-Specific Significance Comparison')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(comparison_df['gene'], rotation=45, ha='right')
    ax2.legend()
    ax2.set_ylim(0, 1.2)
    
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/gene_specific_comparison.png", 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print("   ✅ Gene-specific plots saved")

# ============================================================================
# 💾 EXPORT FUNCTIONS
# ============================================================================

def export_comparison_results(merged_data, comparison_df, lemur_top, deseq2_top, overlap_genes):
    """Export all comparison results"""
    print(f"\n💾 EXPORTING RESULTS")
    print("=" * 20)
    
    # Export merged comparison data
    if merged_data is not None:
        merged_data.to_csv(f"{TABLES_DIR}/lemur_vs_deseq2_merged_comparison.csv")
        print("   ✅ Merged comparison data exported")
    
    # Export gene-specific comparison
    if comparison_df is not None:
        comparison_df.to_csv(f"{TABLES_DIR}/gene_specific_comparison.csv", index=False)
        print("   ✅ Gene-specific comparison exported")
    
    # Export top genes lists
    if lemur_top is not None:
        lemur_top.to_csv(f"{TABLES_DIR}/lemur_top_genes.csv")
        print("   ✅ LEMUR top genes exported")
    
    if deseq2_top is not None:
        deseq2_top.to_csv(f"{TABLES_DIR}/deseq2_top_genes.csv")
        print("   ✅ DESeq2 top genes exported")
    
    # Export overlap analysis
    if overlap_genes:
        overlap_df = pd.DataFrame({'overlapping_genes': list(overlap_genes)})
        overlap_df.to_csv(f"{TABLES_DIR}/overlapping_top_genes.csv", index=False)
        print("   ✅ Overlapping genes exported")
    
    # Create summary report
    create_summary_report(merged_data, comparison_df, overlap_genes)

def create_summary_report(merged_data, comparison_df, overlap_genes):
    """Create comprehensive summary report"""
    print("   Creating summary report...")
    
    report_lines = [
        "# LEMUR vs pyDESeq2 Comparison Report",
        "=" * 40,
        "",
        "## Analysis Overview",
        f"- LEMUR analysis: Single-cell latent space modeling",
        f"- DESeq2 analysis: Pseudobulk differential expression",
        f"- Comparison focus: Main effect (OUD vs Control)",
        "",
        "## Data Summary",
    ]
    
    if merged_data is not None:
        # Calculate statistics
        pearson_r, pearson_p = pearsonr(merged_data['lemur_coef'], merged_data['deseq2_lfc'])
        spearman_r, spearman_p = spearmanr(merged_data['lemur_coef'], merged_data['deseq2_lfc'])
        
        lemur_direction = np.sign(merged_data['lemur_coef'])
        deseq2_direction = np.sign(merged_data['deseq2_lfc'])
        direction_agreement = (lemur_direction == deseq2_direction).mean()
        
        lemur_sig_count = merged_data['lemur_significant'].sum()
        deseq2_sig_count = merged_data['deseq2_significant'].sum()
        both_sig_count = (merged_data['lemur_significant'] & merged_data['deseq2_significant']).sum()
        
        report_lines.extend([
            f"- Common genes analyzed: {len(merged_data)}",
            f"- LEMUR significant genes: {lemur_sig_count}",
            f"- DESeq2 significant genes: {deseq2_sig_count}",
            f"- Both methods significant: {both_sig_count}",
            "",
            "## Correlation Analysis",
            f"- Pearson correlation (effect sizes): r = {pearson_r:.3f}, p = {pearson_p:.2e}",
            f"- Spearman correlation (ranks): ρ = {spearman_r:.3f}, p = {spearman_p:.2e}",
            f"- Direction agreement: {direction_agreement:.1%}",
            "",
        ])
    
    if comparison_df is not None:
        direction_agree_count = comparison_df['direction_agree'].sum()
        report_lines.extend([
            "## Gene-Specific Analysis",
            f"- Genes of interest analyzed: {len(comparison_df)}",
            f"- Direction agreement in specific genes: {direction_agree_count}/{len(comparison_df)}",
            "",
        ])
    
    if overlap_genes:
        report_lines.extend([
            "## Top Gene Overlap",
            f"- Overlapping genes in top lists: {len(overlap_genes)}",
            f"- Overlapping genes: {', '.join(list(overlap_genes)[:10])}",
            "",
        ])
    
    report_lines.extend([
        "## Key Findings",
        "1. Method comparison reveals...",
        "2. Effect size correlations show...",
        "3. Significance agreement indicates...",
        "",
        "## Recommendations",
        "1. Consider both methods for comprehensive analysis",
        "2. Focus on genes significant in both methods",
        "3. Investigate discrepancies for biological insights",
    ])
    
    # Write report
    with open(f"{TABLES_DIR}/comparison_summary_report.md", 'w') as f:
        f.write('\n'.join(report_lines))
    
    print("   ✅ Summary report created")

# ============================================================================
# 🚀 MAIN EXECUTION
# ============================================================================

def main():
    """Main comparison workflow"""
    print(f"\n🚀 STARTING COMPARISON ANALYSIS")
    print("=" * 32)
    
    # 1. Load data
    lemur_main, lemur_sex, deseq2_main, deseq2_sex = load_comparison_data()
    if any(x is None for x in [lemur_main, lemur_sex, deseq2_main, deseq2_sex]):
        print("❌ Data loading failed")
        return
    
    # 2. Compare main effects
    merged_data = compare_main_effects(lemur_main, deseq2_main)
    
    # 3. Compare top genes
    lemur_top, deseq2_top, overlap_genes = compare_top_genes(lemur_main, deseq2_main)
    
    # 4. Analyze specific genes
    comparison_df = analyze_specific_genes(lemur_main, deseq2_main)
    
    # 5. Create visualizations
    create_comparison_plots(merged_data, comparison_df)
    if comparison_df is not None:
        create_gene_comparison_plot(comparison_df)
    
    # 6. Export results
    export_comparison_results(merged_data, comparison_df, lemur_top, deseq2_top, overlap_genes)
    
    print(f"\n✅ COMPARISON ANALYSIS COMPLETE")
    print("=" * 32)
    print(f"Results saved to: {OUTPUT_DIR}")
    print(f"• Plots: {PLOTS_DIR}/")
    print(f"• Tables: {TABLES_DIR}/")
    
    # Final summary
    if merged_data is not None:
        pearson_r, _ = pearsonr(merged_data['lemur_coef'], merged_data['deseq2_lfc'])
        direction_agreement = (np.sign(merged_data['lemur_coef']) == np.sign(merged_data['deseq2_lfc'])).mean()
        
        print(f"\n📊 KEY FINDINGS:")
        print(f"  • Effect size correlation: r = {pearson_r:.3f}")
        print(f"  • Direction agreement: {direction_agreement:.1%}")
        print(f"  • Common genes analyzed: {len(merged_data)}")
        
        if overlap_genes:
            print(f"  • Top gene overlap: {len(overlap_genes)} genes")

if __name__ == "__main__":
    main()