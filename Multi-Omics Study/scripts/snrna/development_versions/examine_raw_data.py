#!/usr/bin/env python3
"""
🔍 Examine Raw Data Structure for LEMUR Analysis
GSE225158 - Explore raw unintegrated data structure and metadata

This script examines the raw data to understand:
1. Available metadata columns
2. Sample distribution across conditions
3. Data structure and format
4. Prepare for proper LEMUR experimental design

Author: Research Team
Date: 2024
"""

import os
import pandas as pd
import scanpy as sc
import numpy as np

# Configure scanpy
sc.settings.verbosity = 1

# File paths
BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
RAW_H5AD = f"{BASE_DIR}/data/raw/snrna/GSE225158_BU_OUD_Striatum_refined_all_SeuratObj_N22.h5ad"

print("🔍 EXAMINING RAW DATA STRUCTURE")
print("===============================")
print(f"File: {RAW_H5AD}")

# Check if file exists
if not os.path.exists(RAW_H5AD):
    print(f"❌ File not found: {RAW_H5AD}")
    exit(1)

# Load the raw data
print("\n📁 Loading raw data...")
adata = sc.read_h5ad(RAW_H5AD)

print(f"✅ Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Examine data structure
print("\n📊 DATA STRUCTURE")
print("=" * 20)

print(f"Data type: {type(adata.X)}")
print(f"Data shape: {adata.X.shape}")
print(f"Data dtype: {adata.X.dtype}")

# Check if sparse
if hasattr(adata.X, 'nnz'):
    print(f"Sparse matrix - Non-zero elements: {adata.X.nnz:,}")
    print(f"Sparsity: {(1 - adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1])) * 100:.2f}%")

# Sample some expression values
print(f"\nSample expression values:")
if hasattr(adata.X, 'todense'):
    sample_vals = adata.X[:5, :5].todense()
else:
    sample_vals = adata.X[:5, :5]
print(sample_vals)

print("\n🏷️  METADATA COLUMNS")
print("=" * 20)
print("Available obs columns:")
for i, col in enumerate(adata.obs.columns, 1):
    print(f"  {i:2d}. {col}")

print(f"\nTotal metadata columns: {len(adata.obs.columns)}")

# Examine key columns for LEMUR design
print("\n🧬 KEY METADATA ANALYSIS")
print("=" * 25)

# Check for condition/diagnosis column
condition_cols = [col for col in adata.obs.columns if any(term in col.lower() for term in ['condition', 'diagnosis', 'group', 'oud', 'control'])]
print(f"Potential condition columns: {condition_cols}")

# Check for sex column
sex_cols = [col for col in adata.obs.columns if any(term in col.lower() for term in ['sex', 'gender', 'male', 'female'])]
print(f"Potential sex columns: {sex_cols}")

# Check for sample/patient ID columns
id_cols = [col for col in adata.obs.columns if any(term in col.lower() for term in ['sample', 'patient', 'subject', 'id', 'donor'])]
print(f"Potential ID columns: {id_cols}")

# Check for cell type columns
celltype_cols = [col for col in adata.obs.columns if any(term in col.lower() for term in ['cell', 'type', 'cluster', 'annotation'])]
print(f"Potential cell type columns: {celltype_cols}")

print("\n📋 DETAILED COLUMN ANALYSIS")
print("=" * 28)

# Analyze each column
for col in adata.obs.columns:
    unique_vals = adata.obs[col].nunique()
    print(f"\n{col}:")
    print(f"  Type: {adata.obs[col].dtype}")
    print(f"  Unique values: {unique_vals}")
    
    if unique_vals <= 20:  # Show unique values if reasonable number
        value_counts = adata.obs[col].value_counts()
        print(f"  Values:")
        for val, count in value_counts.items():
            print(f"    {val}: {count:,} cells")
    else:
        print(f"  Sample values: {list(adata.obs[col].unique()[:5])}")

# Check layers
print(f"\n🧭 DATA LAYERS")
print("=" * 15)
print(f"Available layers: {list(adata.layers.keys())}")

for layer_name in adata.layers.keys():
    layer_data = adata.layers[layer_name]
    print(f"\n{layer_name}:")
    print(f"  Type: {type(layer_data)}")
    print(f"  Dtype: {layer_data.dtype}")
    if hasattr(layer_data, 'nnz'):
        print(f"  Sparsity: {(1 - layer_data.nnz / (layer_data.shape[0] * layer_data.shape[1])) * 100:.2f}%")

# Check var (gene) information
print(f"\n🧬 GENE INFORMATION")
print("=" * 19)
print(f"Gene columns: {list(adata.var.columns)}")

for col in adata.var.columns:
    unique_vals = adata.var[col].nunique()
    print(f"\n{col}:")
    print(f"  Type: {adata.var[col].dtype}")
    print(f"  Unique values: {unique_vals}")
    if unique_vals <= 10:
        print(f"  Values: {list(adata.var[col].unique())}")

# Check uns (unstructured) information
print(f"\n📊 UNSTRUCTURED DATA")
print("=" * 20)
print(f"Available uns keys: {list(adata.uns.keys())}")

# Generate summary for LEMUR design
print(f"\n🌊 LEMUR DESIGN RECOMMENDATIONS")
print("=" * 32)

# Try to identify the best columns for LEMUR design
recommended_design = []

# Condition column
if condition_cols:
    main_condition = condition_cols[0]
    print(f"✅ Condition column: '{main_condition}'")
    values = adata.obs[main_condition].value_counts()
    print(f"   Values: {dict(values)}")
    recommended_design.append(main_condition)
else:
    print("❌ No clear condition column found")

# Sex column
if sex_cols:
    main_sex = sex_cols[0]
    print(f"✅ Sex column: '{main_sex}'")
    values = adata.obs[main_sex].value_counts()
    print(f"   Values: {dict(values)}")
    recommended_design.append(main_sex)
else:
    print("❌ No clear sex column found")

# Sample ID column
if id_cols:
    main_id = id_cols[0]
    print(f"✅ Sample ID column: '{main_id}'")
    n_samples = adata.obs[main_id].nunique()
    print(f"   Number of samples: {n_samples}")
    recommended_design.append(main_id)
else:
    print("❌ No clear sample ID column found")

# Suggest LEMUR design formula
if len(recommended_design) >= 2:
    if len(recommended_design) == 3:  # condition, sex, sample_id
        design_formula = f"~ {recommended_design[2]} + {recommended_design[0]} + {recommended_design[1]} + {recommended_design[0]}:{recommended_design[1]}"
        print(f"\n🔬 Suggested LEMUR design formula:")
        print(f"   {design_formula}")
        print(f"   This includes:")
        print(f"   - Random effect: {recommended_design[2]} (sample ID)")
        print(f"   - Main effect: {recommended_design[0]} (condition)")
        print(f"   - Main effect: {recommended_design[1]} (sex)")
        print(f"   - Interaction: {recommended_design[0]}:{recommended_design[1]} (condition × sex)")
    else:
        print(f"\n⚠️  Limited design options with available columns")

# Check data quality for LEMUR
print(f"\n✅ DATA QUALITY CHECKS")
print("=" * 23)

# Check for missing values
missing_obs = adata.obs.isnull().sum().sum()
print(f"Missing values in obs: {missing_obs}")

# Check expression data
expr_zeros = (adata.X == 0).sum()
if hasattr(expr_zeros, 'item'):
    expr_zeros = expr_zeros.item()
total_expr = adata.X.shape[0] * adata.X.shape[1]
zero_fraction = expr_zeros / total_expr
print(f"Zero expression fraction: {zero_fraction:.3f}")

# Gene expression summary
if hasattr(adata.X, 'todense'):
    sample_expr = adata.X[:1000, :].todense()
else:
    sample_expr = adata.X[:1000, :]

print(f"Expression range: {sample_expr.min():.3f} to {sample_expr.max():.3f}")
print(f"Mean expression: {sample_expr.mean():.3f}")

print(f"\n🎯 NEXT STEPS FOR LEMUR")
print("=" * 24)
print("1. ✅ Raw data is available and loaded successfully")
print("2. 🔍 Identify correct column names from analysis above")
print("3. 🌊 Create proper LEMUR script with:")
print("   - Raw unintegrated data")
print("   - Proper experimental design formula")
print("   - test_fraction parameter")
print("   - Multiple contrasts")
print("   - Neighborhood-based DE analysis")
print("4. 📊 Run comprehensive analysis with sex interactions")

print(f"\n✅ EXAMINATION COMPLETE")
print("=" * 23)
print("Review the column analysis above to determine the correct")
print("column names for condition, sex, and sample ID.")