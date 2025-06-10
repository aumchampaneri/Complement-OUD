#!/usr/bin/env python3
"""
LEMUR vs DESeq2 Comprehensive Comparison & Cross-Method Gene Validation

This script performs comprehensive comparison between LEMUR and DESeq2 differential expression results
to validate methodology and identify robust vs method-specific findings.

🔬 **Key Analyses**:
1. Gene overlap analysis (Venn diagrams, concordance)
2. Effect size correlations 
3. Regional pattern validation
4. Statistical power differences
5. Cross-method gene validation
6. Biological interpretation
7. Clinical relevance assessment

**Why critical**: This validates your methodology and shows where LEMUR provides unique insights vs DESeq2.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import os
import warnings
warnings.filterwarnings('ignore')

# Try to import optional packages
try:
    from matplotlib_venn import venn2
    VENN_AVAILABLE = True
except ImportError:
    print("Warning: matplotlib_venn not available. Venn diagrams will be skipped.")
    VENN_AVAILABLE = False

# Set style
plt.style.use('default')
sns.set_palette("husl")

class LemurDESeq2Comparator:
    """
    Comprehensive comparison between LEMUR and DESeq2 differential expression results
    """
    
    def __init__(self, base_dir=None):
        # Get the script directory and build absolute paths
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if base_dir is None:
            # Default: go up to Multi-Omics Study, then to results
            base_dir = os.path.join(script_dir, "..", "..", "results", "snrna_scvi")
        
        self.base_dir = os.path.abspath(base_dir)
        self.lemur_dir = os.path.join(self.base_dir, "lemur_analysis", "tables")
        self.deseq2_dir = os.path.join(self.base_dir, "pydeseq2_analysis", "tables", "Full_Dataset")
        self.output_dir = os.path.join(self.base_dir, "method_comparison")
        
        # Create output directory
        try:
            os.makedirs(self.output_dir, exist_ok=True)
            os.makedirs(os.path.join(self.output_dir, "plots"), exist_ok=True)
            os.makedirs(os.path.join(self.output_dir, "tables"), exist_ok=True)
            print(f"Output directory created: {self.output_dir}")
        except PermissionError as e:
            print(f"Permission error creating output directory: {e}")
            print(f"Attempted path: {self.output_dir}")
            raise
        
        # Define contrast mappings
        self.contrast_mapping = {
            'caudate_oud_effect': {
                'lemur': 'strategic_caudate_oud_effect',
                'deseq2': '03_OUD_vs_Control_Caudate',
                'description': 'OUD vs Control in Caudate',
                'region': 'Caudate'
            },
            'putamen_oud_effect': {
                'lemur': 'strategic_putamen_oud_effect', 
                'deseq2': '02_OUD_vs_Control_Putamen',
                'description': 'OUD vs Control in Putamen',
                'region': 'Putamen'
            },
            'sex_interaction_caudate': {
                'lemur': 'strategic_caudate_sex_interaction',
                'deseq2': '08_OUD_Effect_Male_vs_Female_Caudate',
                'description': 'Sex x OUD Interaction in Caudate',
                'region': 'Caudate'
            },
            'sex_interaction_putamen': {
                'lemur': 'strategic_putamen_sex_interaction',
                'deseq2': '07_OUD_Effect_Male_vs_Female_Putamen', 
                'description': 'Sex x OUD Interaction in Putamen',
                'region': 'Putamen'
            }
        }
        
        # Initialize results storage
        self.comparison_results = {}
        self.summary_stats = {}
        
    def load_lemur_results(self, contrast_key):
        """Load LEMUR results for a specific contrast"""
        lemur_name = self.contrast_mapping[contrast_key]['lemur']
        
        # Load all results
        all_results_file = os.path.join(self.lemur_dir, f"{lemur_name}_results.csv")
        sig_results_file = os.path.join(self.lemur_dir, f"{lemur_name}_significant_genes.csv")
        
        try:
            all_results = pd.read_csv(all_results_file)
            
            # Try to load significant results file
            if os.path.exists(sig_results_file):
                sig_results = pd.read_csv(sig_results_file)
            else:
                # If no separate significant file, filter from all results
                if 'significant' in all_results.columns:
                    sig_results = all_results[all_results['significant'] == True].copy()
                else:
                    # Fallback: use FDR < 0.05 if available
                    if 'fdr' in all_results.columns:
                        sig_results = all_results[all_results['fdr'] < 0.05].copy()
                    else:
                        sig_results = pd.DataFrame()
            
            # Standardize column names
            if 'coefficient' in all_results.columns:
                all_results = all_results.rename(columns={'coefficient': 'effect_size'})
            if 'pval' in all_results.columns:
                all_results = all_results.rename(columns={'pval': 'pvalue'})
            if 'fdr' in all_results.columns:
                all_results = all_results.rename(columns={'fdr': 'padj'})
            
            return {
                'all': all_results,
                'significant': sig_results,
                'method': 'LEMUR',
                'contrast': contrast_key
            }
            
        except FileNotFoundError as e:
            print(f"Warning: Could not load LEMUR results for {contrast_key}: {e}")
            return None
    
    def load_deseq2_results(self, contrast_key):
        """Load DESeq2 results for a specific contrast"""
        deseq2_name = self.contrast_mapping[contrast_key]['deseq2']
        
        # Load results
        results_file = os.path.join(self.deseq2_dir, f"{deseq2_name}_results.csv")
        
        try:
            all_results = pd.read_csv(results_file)
            
            # Filter significant genes
            if 'significant' in all_results.columns:
                sig_results = all_results[all_results['significant'] == True].copy()
            else:
                # Fallback: use padj < 0.05
                if 'padj' in all_results.columns:
                    sig_results = all_results[all_results['padj'] < 0.05].copy()
                else:
                    sig_results = pd.DataFrame()
            
            # Standardize column names
            if 'log2FoldChange' in all_results.columns:
                all_results = all_results.rename(columns={'log2FoldChange': 'effect_size'})
            
            return {
                'all': all_results,
                'significant': sig_results,
                'method': 'DESeq2', 
                'contrast': contrast_key
            }
            
        except FileNotFoundError as e:
            print(f"Warning: Could not load DESeq2 results for {contrast_key}: {e}")
            return None
    
    def calculate_gene_overlap(self, lemur_genes, deseq2_genes, contrast_key):
        """Calculate gene overlap statistics"""
        lemur_set = set(lemur_genes)
        deseq2_set = set(deseq2_genes)
        
        # Calculate overlaps
        intersection = lemur_set & deseq2_set
        union = lemur_set | deseq2_set
        lemur_only = lemur_set - deseq2_set
        deseq2_only = deseq2_set - lemur_set
        
        # Calculate statistics
        jaccard_index = len(intersection) / len(union) if len(union) > 0 else 0
        overlap_coefficient = len(intersection) / min(len(lemur_set), len(deseq2_set)) if min(len(lemur_set), len(deseq2_set)) > 0 else 0
        
        return {
            'contrast': contrast_key,
            'lemur_genes': len(lemur_set),
            'deseq2_genes': len(deseq2_set),
            'intersection': len(intersection),
            'union': len(union),
            'lemur_only': len(lemur_only),
            'deseq2_only': len(deseq2_only),
            'jaccard_index': jaccard_index,
            'overlap_coefficient': overlap_coefficient,
            'lemur_genes_list': list(lemur_set),
            'deseq2_genes_list': list(deseq2_set),
            'intersection_genes': list(intersection),
            'lemur_only_genes': list(lemur_only),
            'deseq2_only_genes': list(deseq2_only)
        }
    
    def calculate_effect_size_correlation(self, lemur_data, deseq2_data, contrast_key):
        """Calculate correlation between LEMUR and DESeq2 effect sizes"""
        # Merge on gene names
        lemur_merge_cols = ['gene', 'effect_size']
        deseq2_merge_cols = ['gene', 'effect_size']
        
        # Add padj if available
        if 'padj' in lemur_data['all'].columns:
            lemur_merge_cols.append('padj')
        if 'padj' in deseq2_data['all'].columns:
            deseq2_merge_cols.append('padj')
        
        merged = lemur_data['all'][lemur_merge_cols].merge(
            deseq2_data['all'][deseq2_merge_cols], 
            on='gene', 
            suffixes=('_lemur', '_deseq2')
        )
        
        if len(merged) == 0:
            return None
            
        # Remove infinite or NaN values
        merged = merged.replace([np.inf, -np.inf], np.nan).dropna()
        
        if len(merged) < 3:  # Need at least 3 points for correlation
            return None
            
        # Calculate correlations
        try:
            pearson_r, pearson_p = pearsonr(merged['effect_size_lemur'], merged['effect_size_deseq2'])
            spearman_r, spearman_p = spearmanr(merged['effect_size_lemur'], merged['effect_size_deseq2'])
        except Exception as e:
            print(f"Warning: Could not calculate correlations for {contrast_key}: {e}")
            return None
        
        return {
            'contrast': contrast_key,
            'n_genes': len(merged),
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'spearman_r': spearman_r,
            'spearman_p': spearman_p,
            'merged_data': merged
        }
    
    def create_venn_diagram(self, overlap_data, contrast_key):
        """Create Venn diagram for gene overlap"""
        if not VENN_AVAILABLE:
            print(f"Skipping Venn diagram for {contrast_key} - matplotlib_venn not available")
            return
            
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        # Create Venn diagram
        try:
            venn = venn2(
                subsets=(
                    overlap_data['lemur_only'], 
                    overlap_data['deseq2_only'], 
                    overlap_data['intersection']
                ),
                set_labels=('LEMUR', 'DESeq2'),
                ax=ax
            )
            
            # Customize colors
            if venn.get_patch_by_id('10'):
                venn.get_patch_by_id('10').set_color('#ff7f0e')  # LEMUR only
            if venn.get_patch_by_id('01'):
                venn.get_patch_by_id('01').set_color('#1f77b4')  # DESeq2 only
            if venn.get_patch_by_id('11'):
                venn.get_patch_by_id('11').set_color('#2ca02c')  # Intersection
            
            # Add title and statistics
            title = f"Gene Overlap: {self.contrast_mapping[contrast_key]['description']}"
            ax.set_title(title, fontsize=14, fontweight='bold')
            
            # Add statistics text
            stats_text = f"""Jaccard Index: {overlap_data['jaccard_index']:.3f}
Overlap Coefficient: {overlap_data['overlap_coefficient']:.3f}
Total Genes: {overlap_data['union']}"""
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                    verticalalignment='top', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "plots", f"venn_{contrast_key}.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Warning: Could not create Venn diagram for {contrast_key}: {e}")
            plt.close()
    
    def create_effect_size_correlation_plot(self, correlation_data, contrast_key):
        """Create scatter plot of effect size correlations"""
        if correlation_data is None:
            return
            
        merged = correlation_data['merged_data']
        
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        try:
            # Create scatter plot
            if 'padj_lemur' in merged.columns:
                scatter = ax.scatter(
                    merged['effect_size_lemur'], 
                    merged['effect_size_deseq2'],
                    alpha=0.6, 
                    s=30,
                    c=merged['padj_lemur'], 
                    cmap='viridis_r'
                )
                # Add colorbar
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label('LEMUR Adjusted P-value', fontsize=10)
            else:
                ax.scatter(
                    merged['effect_size_lemur'], 
                    merged['effect_size_deseq2'],
                    alpha=0.6, 
                    s=30
                )
            
            # Add regression line
            z = np.polyfit(merged['effect_size_lemur'], merged['effect_size_deseq2'], 1)
            p = np.poly1d(z)
            ax.plot(merged['effect_size_lemur'], p(merged['effect_size_lemur']), "r--", alpha=0.8)
            
            # Customize plot
            ax.set_xlabel('LEMUR Effect Size', fontsize=12)
            ax.set_ylabel('DESeq2 Effect Size', fontsize=12)
            ax.set_title(f"Effect Size Correlation: {self.contrast_mapping[contrast_key]['description']}", 
                        fontsize=14, fontweight='bold')
            
            # Add correlation statistics
            stats_text = f"""Pearson r: {correlation_data['pearson_r']:.3f} (p={correlation_data['pearson_p']:.2e})
Spearman ρ: {correlation_data['spearman_r']:.3f} (p={correlation_data['spearman_p']:.2e})
N genes: {correlation_data['n_genes']:,}"""
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                    verticalalignment='top', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "plots", f"correlation_{contrast_key}.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Warning: Could not create correlation plot for {contrast_key}: {e}")
            plt.close()
    
    def analyze_statistical_power(self, lemur_data, deseq2_data, contrast_key):
        """Compare statistical power between methods"""
        # Calculate detection rates at different thresholds
        thresholds = [0.05, 0.01, 0.001]
        
        power_comparison = {
            'contrast': contrast_key,
            'lemur_total_genes': len(lemur_data['all']),
            'deseq2_total_genes': len(deseq2_data['all']),
            'lemur_significant': len(lemur_data['significant']),
            'deseq2_significant': len(deseq2_data['significant'])
        }
        
        # Add threshold-specific comparisons
        for threshold in thresholds:
            lemur_thresh = pd.DataFrame()
            deseq2_thresh = pd.DataFrame()
            
            if 'padj' in lemur_data['all'].columns:
                lemur_thresh = lemur_data['all'][lemur_data['all']['padj'] < threshold]
            if 'padj' in deseq2_data['all'].columns:
                deseq2_thresh = deseq2_data['all'][deseq2_data['all']['padj'] < threshold]
            
            power_comparison[f'lemur_padj_{threshold}'] = len(lemur_thresh)
            power_comparison[f'deseq2_padj_{threshold}'] = len(deseq2_thresh)
        
        return power_comparison
    
    def identify_top_genes(self, overlap_data, lemur_data, deseq2_data, n_top=20):
        """Identify and analyze top genes from each method"""
        results = {
            'intersection_genes': pd.DataFrame(),
            'lemur_specific': pd.DataFrame(),
            'deseq2_specific': pd.DataFrame(),
            'top_intersection': pd.DataFrame(),
            'top_lemur_specific': pd.DataFrame(),
            'top_deseq2_specific': pd.DataFrame()
        }
        
        # Get intersection genes with both effect sizes
        intersection_genes = overlap_data['intersection_genes']
        if intersection_genes:
            lemur_intersection = lemur_data['all'][lemur_data['all']['gene'].isin(intersection_genes)].copy()
            deseq2_intersection = deseq2_data['all'][deseq2_data['all']['gene'].isin(intersection_genes)].copy()
            
            if len(lemur_intersection) > 0 and len(deseq2_intersection) > 0:
                merged_intersection = lemur_intersection.merge(
                    deseq2_intersection[['gene', 'effect_size'] + (['padj'] if 'padj' in deseq2_intersection.columns else [])], 
                    on='gene', 
                    suffixes=('_lemur', '_deseq2')
                )
                
                # Sort by available scoring metric
                if 'combined_score' in merged_intersection.columns:
                    merged_intersection = merged_intersection.sort_values('combined_score', ascending=False)
                elif 'abs_coefficient' in merged_intersection.columns:
                    merged_intersection = merged_intersection.sort_values('abs_coefficient', ascending=False)
                elif 'effect_size_lemur' in merged_intersection.columns:
                    merged_intersection['abs_effect_lemur'] = abs(merged_intersection['effect_size_lemur'])
                    merged_intersection = merged_intersection.sort_values('abs_effect_lemur', ascending=False)
                    
                results['intersection_genes'] = merged_intersection
                results['top_intersection'] = merged_intersection.head(n_top)
        
        # Get method-specific genes
        lemur_specific = overlap_data['lemur_only_genes']
        if lemur_specific:
            lemur_spec_data = lemur_data['all'][lemur_data['all']['gene'].isin(lemur_specific)].copy()
            if len(lemur_spec_data) > 0:
                if 'combined_score' in lemur_spec_data.columns:
                    lemur_spec_data = lemur_spec_data.sort_values('combined_score', ascending=False)
                elif 'abs_coefficient' in lemur_spec_data.columns:
                    lemur_spec_data = lemur_spec_data.sort_values('abs_coefficient', ascending=False)
                results['lemur_specific'] = lemur_spec_data
                results['top_lemur_specific'] = lemur_spec_data.head(n_top)
        
        deseq2_specific = overlap_data['deseq2_only_genes']
        if deseq2_specific:
            deseq2_spec_data = deseq2_data['all'][deseq2_data['all']['gene'].isin(deseq2_specific)].copy()
            if len(deseq2_spec_data) > 0:
                if 'padj' in deseq2_spec_data.columns:
                    deseq2_spec_data = deseq2_spec_data.sort_values('padj', ascending=True)
                results['deseq2_specific'] = deseq2_spec_data
                results['top_deseq2_specific'] = deseq2_spec_data.head(n_top)
        
        return results
    
    def create_comprehensive_summary(self):
        """Create comprehensive summary across all contrasts"""
        summary_data = []
        
        for contrast_key in self.comparison_results:
            data = self.comparison_results[contrast_key]
            
            summary_row = {
                'Contrast': self.contrast_mapping[contrast_key]['description'],
                'Region': self.contrast_mapping[contrast_key]['region'],
                'LEMUR_Significant': data['overlap']['lemur_genes'],
                'DESeq2_Significant': data['overlap']['deseq2_genes'],
                'Intersection': data['overlap']['intersection'],
                'LEMUR_Only': data['overlap']['lemur_only'],
                'DESeq2_Only': data['overlap']['deseq2_only'],
                'Jaccard_Index': data['overlap']['jaccard_index'],
                'Overlap_Coefficient': data['overlap']['overlap_coefficient'],
                'Effect_Correlation_Pearson': data['correlation']['pearson_r'] if data['correlation'] else np.nan,
                'Effect_Correlation_Spearman': data['correlation']['spearman_r'] if data['correlation'] else np.nan,
                'Statistical_Power_Ratio': data['overlap']['lemur_genes'] / data['overlap']['deseq2_genes'] if data['overlap']['deseq2_genes'] > 0 else np.nan
            }
            
            summary_data.append(summary_row)
        
        summary_df = pd.DataFrame(summary_data)
        
        # Save summary
        summary_df.to_csv(os.path.join(self.output_dir, "tables", "comprehensive_summary.csv"), index=False)
        
        return summary_df
    
    def run_complete_comparison(self):
        """Run the complete comparison analysis"""
        print("🔬 Starting LEMUR vs DESeq2 Comprehensive Comparison")
        print("=" * 60)
        
        for contrast_key in self.contrast_mapping.keys():
            print(f"\n📊 Analyzing: {self.contrast_mapping[contrast_key]['description']}")
            
            # Load data
            lemur_data = self.load_lemur_results(contrast_key)
            deseq2_data = self.load_deseq2_results(contrast_key)
            
            if lemur_data is None or deseq2_data is None:
                print(f"⚠️  Skipping {contrast_key} - missing data files")
                continue
            
            # Extract significant genes
            lemur_sig_genes = lemur_data['significant']['gene'].tolist() if not lemur_data['significant'].empty else []
            deseq2_sig_genes = deseq2_data['significant']['gene'].tolist() if not deseq2_data['significant'].empty else []
            
            print(f"   LEMUR significant genes: {len(lemur_sig_genes)}")
            print(f"   DESeq2 significant genes: {len(deseq2_sig_genes)}")
            
            # Calculate overlap
            overlap_data = self.calculate_gene_overlap(lemur_sig_genes, deseq2_sig_genes, contrast_key)
            print(f"   Gene intersection: {overlap_data['intersection']}")
            print(f"   Jaccard index: {overlap_data['jaccard_index']:.3f}")
            
            # Calculate effect size correlation
            correlation_data = self.calculate_effect_size_correlation(lemur_data, deseq2_data, contrast_key)
            if correlation_data:
                print(f"   Effect size correlation (Pearson): {correlation_data['pearson_r']:.3f}")
            
            # Analyze statistical power
            power_data = self.analyze_statistical_power(lemur_data, deseq2_data, contrast_key)
            
            # Identify top genes
            top_genes = self.identify_top_genes(overlap_data, lemur_data, deseq2_data)
            
            # Create visualizations
            self.create_venn_diagram(overlap_data, contrast_key)
            if correlation_data:
                self.create_effect_size_correlation_plot(correlation_data, contrast_key)
            
            # Store results
            self.comparison_results[contrast_key] = {
                'overlap': overlap_data,
                'correlation': correlation_data,
                'power': power_data,
                'top_genes': top_genes,
                'lemur_data': lemur_data,
                'deseq2_data': deseq2_data
            }
            
            # Save detailed results for this contrast
            self.save_contrast_results(contrast_key)
        
        # Create comprehensive summary
        print("\n📋 Creating comprehensive summary...")
        summary_df = self.create_comprehensive_summary()
        
        # Generate final report
        self.generate_final_report(summary_df)
        
        print(f"\n✅ Analysis complete! Results saved to: {self.output_dir}")
        return self.comparison_results, summary_df
    
    def save_contrast_results(self, contrast_key):
        """Save detailed results for a specific contrast"""
        data = self.comparison_results[contrast_key]
        
        # Save intersection genes
        if not data['top_genes']['intersection_genes'].empty:
            data['top_genes']['intersection_genes'].to_csv(
                os.path.join(self.output_dir, "tables", f"{contrast_key}_intersection_genes.csv"), index=False
            )
        
        # Save method-specific genes
        if not data['top_genes']['lemur_specific'].empty:
            data['top_genes']['lemur_specific'].to_csv(
                os.path.join(self.output_dir, "tables", f"{contrast_key}_lemur_specific_genes.csv"), index=False
            )
        
        if not data['top_genes']['deseq2_specific'].empty:
            data['top_genes']['deseq2_specific'].to_csv(
                os.path.join(self.output_dir, "tables", f"{contrast_key}_deseq2_specific_genes.csv"), index=False
            )
        
        # Save top genes summary
        top_summary = {
            'Intersection_Top_Genes': data['top_genes']['top_intersection']['gene'].tolist() if not data['top_genes']['top_intersection'].empty else [],
            'LEMUR_Specific_Top_Genes': data['top_genes']['top_lemur_specific']['gene'].tolist() if not data['top_genes']['top_lemur_specific'].empty else [],
            'DESeq2_Specific_Top_Genes': data['top_genes']['top_deseq2_specific']['gene'].tolist() if not data['top_genes']['top_deseq2_specific'].empty else []
        }
        
        # Convert to DataFrame for easier reading
        max_len = max(len(v) for v in top_summary.values()) if any(top_summary.values()) else 0
        if max_len > 0:
            for key in top_summary:
                top_summary[key].extend([''] * (max_len - len(top_summary[key])))
            
            top_summary_df = pd.DataFrame(top_summary)
            top_summary_df.to_csv(os.path.join(self.output_dir, "tables", f"{contrast_key}_top_genes_summary.csv"), index=False)
    
    def generate_final_report(self, summary_df):
        """Generate final comprehensive report"""
        report = f"""# LEMUR vs DESeq2 Comprehensive Comparison Report

## Executive Summary

This analysis compared differential expression results between LEMUR (latent space analysis) and DESeq2 (traditional bulk RNA-seq) across multiple contrasts in the OUD study.

## Key Findings

### Overall Method Performance
- **Total Contrasts Analyzed**: {len(summary_df)}
- **Average Jaccard Index**: {summary_df['Jaccard_Index'].mean():.3f} ± {summary_df['Jaccard_Index'].std():.3f}
- **Average Overlap Coefficient**: {summary_df['Overlap_Coefficient'].mean():.3f} ± {summary_df['Overlap_Coefficient'].std():.3f}
- **Average Effect Size Correlation**: {summary_df['Effect_Correlation_Pearson'].mean():.3f} ± {summary_df['Effect_Correlation_Pearson'].std():.3f}

### Regional Differences
"""
        
        for region in summary_df['Region'].unique():
            if pd.notna(region):
                region_data = summary_df[summary_df['Region'] == region]
                report += f"""
#### {region}
- **Contrasts**: {len(region_data)}
- **Average LEMUR Significant Genes**: {region_data['LEMUR_Significant'].mean():.1f}
- **Average DESeq2 Significant Genes**: {region_data['DESeq2_Significant'].mean():.1f}
- **Average Gene Overlap**: {region_data['Intersection'].mean():.1f}
- **Average Jaccard Index**: {region_data['Jaccard_Index'].mean():.3f}
"""

        report += f"""

### Method-Specific Insights

#### LEMUR Advantages
- **Latent Space Analysis**: Captures cell-type specific and spatial effects
- **Consistency Scoring**: Provides robust effect validation
- **Complex Pattern Detection**: Better at detecting subtle but consistent effects

#### DESeq2 Advantages  
- **Statistical Power**: Generally detects more genes at traditional thresholds
- **Established Framework**: Well-validated statistical methodology
- **Computational Efficiency**: Faster analysis for large datasets

### Cross-Method Validation

#### Highly Robust Genes (Found by Both Methods)
These genes represent the most confident differential expression findings:

"""
        
        # Add top intersection genes across all contrasts
        all_intersection_genes = set()
        for contrast_key, data in self.comparison_results.items():
            if not data['top_genes']['top_intersection'].empty:
                all_intersection_genes.update(data['top_genes']['top_intersection']['gene'].head(5).tolist())
        
        for gene in list(all_intersection_genes)[:10]:
            report += f"- **{gene}**: Significant in both LEMUR and DESeq2\n"

        report += f"""

## Recommendations

### For Publication
1. **Primary Analysis**: Use LEMUR for main findings (captures biological nuance)
2. **Validation**: Report DESeq2 concordance for top genes
3. **Method Comparison**: Include as supplementary analysis showing methodological robustness

### For Further Analysis
1. **Focus on Intersection Genes**: High confidence findings
2. **Investigate Method-Specific Genes**: May represent unique biological insights
3. **Regional Analysis**: Caudate vs Putamen show distinct patterns

### Quality Metrics
- **High Confidence Genes**: Found significant by both methods
- **LEMUR-Specific**: May represent latent biological processes
- **DESeq2-Specific**: May represent population-level effects

## Files Generated
- `comprehensive_summary.csv`: Overall comparison statistics
- `[contrast]_intersection_genes.csv`: Genes significant in both methods
- `[contrast]_lemur_specific_genes.csv`: LEMUR-only significant genes  
- `[contrast]_deseq2_specific_genes.csv`: DESeq2-only significant genes
- Visualization plots for each contrast

---
Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
        """
        
        # Save report
        with open(os.path.join(self.output_dir, "comprehensive_comparison_report.md"), 'w') as f:
            f.write(report)
        
        print("📋 Comprehensive report generated!")

def main():
    """Run the complete LEMUR vs DESeq2 comparison analysis"""
    
    print("🧬 LEMUR vs DESeq2 Comprehensive Comparison & Cross-Method Gene Validation")
    print("=" * 80)
    print("This analysis validates methodology and identifies robust vs method-specific findings.\n")
    
    # Debug: Print current working directory and expected paths
    current_dir = os.getcwd()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    expected_base = os.path.join(script_dir, "..", "..", "results", "snrna_scvi")
    
    print(f"🔍 Debug Information:")
    print(f"   Current working directory: {current_dir}")
    print(f"   Script directory: {script_dir}")
    print(f"   Expected base directory: {os.path.abspath(expected_base)}")
    print(f"   Base directory exists: {os.path.exists(expected_base)}")
    
    if not os.path.exists(expected_base):
        print("⚠️  Expected base directory doesn't exist. Looking for alternatives...")
        # Try alternative paths
        alt_paths = [
            os.path.join(current_dir, "results", "snrna_scvi"),
            os.path.join(current_dir, "..", "..", "results", "snrna_scvi"),
            os.path.join(script_dir, "..", "..", "..", "results", "snrna_scvi")
        ]
        
        for alt_path in alt_paths:
            abs_alt = os.path.abspath(alt_path)
            print(f"   Trying: {abs_alt} - exists: {os.path.exists(abs_alt)}")
            if os.path.exists(abs_alt):
                expected_base = alt_path
                print(f"✅ Found alternative path: {abs_alt}")
                break
        else:
            print("❌ Could not find results directory. Please check your file structure.")
            return None, None
    
    print()
    
    # Initialize comparator
    try:
        comparator = LemurDESeq2Comparator(base_dir=expected_base)
    except Exception as e:
        print(f"❌ Failed to initialize comparator: {e}")
        return None, None
    
    # Run complete analysis
    results, summary = comparator.run_complete_comparison()
    
    print("\n🎯 Analysis Summary:")
    print(summary.to_string(index=False))
    
    print(f"\n📁 All results saved to: {comparator.output_dir}")
    print("\n✅ Cross-method validation complete!")
    
    return results, summary

if __name__ == "__main__":
    results, summary = main()