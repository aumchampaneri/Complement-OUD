#!/usr/bin/env python3
"""
Quick Regional Distribution Check
Validate regional breakdown for LEMUR contrast design
"""

import scanpy as sc
import pandas as pd
import numpy as np

# Load data
print("🔍 Loading data...")
adata = sc.read_h5ad('/Users/aumchampaneri/Complement-OUD/Multi-Omics Study/data/processed/snrna_scvi/GSE225158_annotated_scvi.h5ad')

print(f"📊 Dataset Overview:")
print(f"   Total cells: {adata.n_obs:,}")
print(f"   Total genes: {adata.n_vars:,}")

# Check regional distribution
if 'Region' in adata.obs.columns:
    print(f"\n🧠 Regional Distribution:")
    region_counts = adata.obs['Region'].value_counts()
    for region, count in region_counts.items():
        print(f"   {region}: {count:,} cells ({count/adata.n_obs*100:.1f}%)")
    
    # Regional breakdown by condition
    print(f"\n📋 Region × Condition:")
    region_condition = pd.crosstab(adata.obs['Region'], adata.obs['Dx_OUD'])
    print(region_condition)
    
    # Regional breakdown by sex
    if 'Sex' in adata.obs.columns:
        print(f"\n⚧️ Region × Sex:")
        region_sex = pd.crosstab(adata.obs['Region'], adata.obs['Sex'])
        print(region_sex)
    
    # Full breakdown: Region × Condition × Sex
    print(f"\n🔬 Complete Stratification:")
    for region in adata.obs['Region'].unique():
        print(f"\n   {region}:")
        region_data = adata.obs[adata.obs['Region'] == region]
        
        for condition in ['OUD', 'None']:
            cond_data = region_data[region_data['Dx_OUD'] == condition]
            condition_label = 'OUD' if condition == 'OUD' else 'Control'
            
            if 'Sex' in adata.obs.columns:
                for sex in ['M', 'F']:
                    sex_count = cond_data[cond_data['Sex'] == sex].shape[0]
                    print(f"     {condition_label} {sex}: {sex_count:,} cells")
            else:
                print(f"     {condition_label}: {cond_data.shape[0]:,} cells")

    # Power analysis for regional contrasts
    print(f"\n⚡ Power Analysis for Regional Contrasts:")
    print(f"   Minimum recommended cells per group: 100")
    
    for region in adata.obs['Region'].unique():
        region_mask = adata.obs['Region'] == region
        oud_count = adata.obs[(region_mask) & (adata.obs['Dx_OUD'] == 'OUD')].shape[0]
        ctrl_count = adata.obs[(region_mask) & (adata.obs['Dx_OUD'] == 'None')].shape[0]
        
        power_status = "✅ Sufficient" if min(oud_count, ctrl_count) >= 100 else "⚠️ Low power"
        print(f"   {region}: OUD={oud_count:,}, Control={ctrl_count:,} → {power_status}")

else:
    print("❌ No 'Region' column found!")

print(f"\n✅ Regional validation complete!")
