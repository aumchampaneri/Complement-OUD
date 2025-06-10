#!/usr/bin/env python3
"""
🛤️ Pathway Analysis using Python decoupleR
Comprehensive functional analysis integrating pyDESeq2 and LEMUR results

Strategy:
1. Load differential expression results from both pyDESeq2 and LEMUR
2. Apply decoupleR pathway analysis to both datasets
3. Compare pathway activities between discrete and continuous approaches
4. Generate integrated functional insights for OUD research
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from pathlib import Path
import pickle
from collections import defaultdict

# decoupleR for pathway analysis
try:
    import decoupler as dc
    DECOUPLER_AVAILABLE = True
    print("✅ decoupleR available")
except ImportError:
    DECOUPLER_AVAILABLE = False
    print("❌ decoupleR not available. Install with: pip install decoupler")

warnings.filterwarnings('ignore')

# ============================================================================
# 📁 CONFIGURATION
# ============================================================================

BASE_DIR = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
PYDESEQ2_DIR = f"{BASE_DIR}/results/snrna_scvi/pydeseq2_analysis"
LEMUR_DIR = f"{BASE_DIR}/results/snrna_scvi/lemur_analysis"
OUTPUT_DIR = f"{BASE_DIR}/results/snrna_scvi/pathway_analysis"
PLOTS_DIR = f"{
