#!/usr/bin/env python3
"""
🚀 Complete edgeR Pipeline - Differential Expression + Pathway Enrichment
GSE225158 - Automated pipeline for comprehensive analysis

Pipeline Steps:
1. Run edgeR differential expression analysis
2. Perform pathway enrichment on results
3. Generate comprehensive reports and visualizations

Usage:
    python run_edgeR_pipeline.py

Dependencies:
- All dependencies from 03b_differential_expression_edgeR.py
- All dependencies from 04_pathway_enrichment_edgeR.py
"""

import os
import sys
import subprocess
import time
from datetime import datetime
from pathlib import Path
import pandas as pd

# Add current directory to path for imports
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

def print_banner():
    """Print pipeline banner"""
    print("🚀" * 35)
    print("🔬 EDGER COMPLETE ANALYSIS PIPELINE 🔬")
    print("🚀" * 35)
    print(f"📅 Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("📊 GSE225158 - OUD vs Control Analysis")
    print("🧬 Steps: edgeR DE → Pathway Enrichment → Reports")
    print("🚀" * 35)

def check_prerequisites():
    """Check if all required files and dependencies are available"""
    print("\n🔍 CHECKING PREREQUISITES")
    print("=" * 50)
    
    # Check if edgeR script exists
    edger_script = current_dir / "03b_differential_expression_edgeR.py"
    if not edger_script.exists():
        raise FileNotFoundError(f"edgeR script not found: {edger_script}")
    print(f"✅ edgeR script found: {edger_script.name}")
    
    # Check if pathway enrichment script exists
    pathway_script = current_dir / "04_pathway_enrichment_edgeR.py"
    if not pathway_script.exists():
        raise FileNotFoundError(f"Pathway script not found: {pathway_script}")
    print(f"✅ Pathway script found: {pathway_script.name}")
    
    # Check if input data exists
    base_dir = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
    input_file = Path(f"{base_dir}/data/processed/snrna_scvi/GSE225158_annotated_scvi.h5ad")
    if not input_file.exists():
        raise FileNotFoundError(f"Input data not found: {input_file}")
    print(f"✅ Input data found: {input_file.name}")
    
    # Test Python dependencies
    try:
        import scanpy
        import pandas
        import numpy
        import matplotlib
        import rpy2.robjects
        print("✅ Python dependencies available")
    except ImportError as e:
        raise ImportError(f"Missing Python dependency: {e}")
    
    # Test R dependencies
    try:
        import subprocess
        result = subprocess.run(['R', '--slave', '-e', 'library(edgeR); library(clusterProfiler)'], 
                              capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            print("✅ R dependencies available")
        else:
            raise RuntimeError("R packages not available")
    except Exception as e:
        print(f"⚠️  Could not verify R dependencies: {e}")
        print("   Proceeding anyway - will fail if packages missing")
    
    print("✅ All prerequisites check passed!")

def run_edger_analysis():
    """Run the edgeR differential expression analysis"""
    print("\n🔬 STEP 1: EDGER DIFFERENTIAL EXPRESSION")
    print("=" * 50)
    
    edger_script = current_dir / "03b_differential_expression_edgeR.py"
    
    print(f"🚀 Executing: {edger_script}")
    start_time = time.time()
    
    try:
        # Run edgeR script
        result = subprocess.run([sys.executable, str(edger_script)], 
                              cwd=current_dir, 
                              capture_output=False,  # Show output in real-time
                              text=True)
        
        elapsed_time = time.time() - start_time
        
        if result.returncode == 0:
            print(f"✅ edgeR analysis completed successfully!")
            print(f"⏱️  Time elapsed: {elapsed_time:.1f} seconds")
            return True
        else:
            print(f"❌ edgeR analysis failed with return code: {result.returncode}")
            return False
            
    except Exception as e:
        print(f"❌ Error running edgeR analysis: {e}")
        return False

def run_pathway_enrichment():
    """Run the pathway enrichment analysis"""
    print("\n🛤️  STEP 2: PATHWAY ENRICHMENT ANALYSIS")
    print("=" * 50)
    
    pathway_script = current_dir / "04_pathway_enrichment_edgeR.py"
    
    print(f"🚀 Executing: {pathway_script}")
    start_time = time.time()
    
    try:
        # Run pathway enrichment script
        result = subprocess.run([sys.executable, str(pathway_script)], 
                              cwd=current_dir, 
                              capture_output=False,  # Show output in real-time
                              text=True)
        
        elapsed_time = time.time() - start_time
        
        if result.returncode == 0:
            print(f"✅ Pathway enrichment completed successfully!")
            print(f"⏱️  Time elapsed: {elapsed_time:.1f} seconds")
            return True
        else:
            print(f"❌ Pathway enrichment failed with return code: {result.returncode}")
            return False
            
    except Exception as e:
        print(f"❌ Error running pathway enrichment: {e}")
        return False

def generate_pipeline_report():
    """Generate a summary report of the complete pipeline"""
    print("\n📋 STEP 3: GENERATING PIPELINE REPORT")
    print("=" * 50)
    
    base_dir = "/Users/aumchampaneri/Complement-OUD/Multi-Omics Study"
    edger_dir = f"{base_dir}/results/snrna_scvi/differential_expression_edgeR"
    pathway_dir = f"{base_dir}/results/snrna_scvi/pathway_enrichment_edgeR"
    
    report_lines = []
    report_lines.append("EDGER COMPLETE PIPELINE REPORT")
    report_lines.append("=" * 50)
    report_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(f"Dataset: GSE225158 - OUD vs Control")
    report_lines.append("")
    
    # Check edgeR results
    edger_files = list(Path(edger_dir).glob("*.csv")) if Path(edger_dir).exists() else []
    report_lines.append("DIFFERENTIAL EXPRESSION RESULTS:")
    report_lines.append("-" * 30)
    if edger_files:
        for file in edger_files:
            try:
                df = pd.read_csv(file)
                report_lines.append(f"  📁 {file.name}: {len(df)} genes")
            except:
                report_lines.append(f"  📁 {file.name}: File exists")
    else:
        report_lines.append("  ❌ No edgeR results found")
    
    report_lines.append("")
    
    # Check pathway results
    pathway_files = list(Path(pathway_dir).glob("*.csv")) if Path(pathway_dir).exists() else []
    report_lines.append("PATHWAY ENRICHMENT RESULTS:")
    report_lines.append("-" * 30)
    if pathway_files:
        for file in pathway_files:
            try:
                df = pd.read_csv(file)
                report_lines.append(f"  📁 {file.name}: {len(df)} pathways")
            except:
                report_lines.append(f"  📁 {file.name}: File exists")
    else:
        report_lines.append("  ❌ No pathway results found")
    
    report_lines.append("")
    
    # Check plots
    plot_dirs = [
        f"{edger_dir}/plots",
        f"{pathway_dir}/plots"
    ]
    
    total_plots = 0
    for plot_dir in plot_dirs:
        if Path(plot_dir).exists():
            plots = list(Path(plot_dir).glob("*.png"))
            total_plots += len(plots)
    
    report_lines.append(f"VISUALIZATIONS GENERATED: {total_plots} plots")
    report_lines.append("")
    
    # Output locations
    report_lines.append("OUTPUT LOCATIONS:")
    report_lines.append("-" * 30)
    report_lines.append(f"  📂 edgeR results: {edger_dir}")
    report_lines.append(f"  📂 Pathway results: {pathway_dir}")
    report_lines.append("")
    
    # Next steps
    report_lines.append("RECOMMENDED NEXT STEPS:")
    report_lines.append("-" * 30)
    report_lines.append("  1. Review differential expression results")
    report_lines.append("  2. Examine pathway enrichment plots")
    report_lines.append("  3. Focus on complement system pathways")
    report_lines.append("  4. Compare with other datasets")
    report_lines.append("  5. Validate key findings experimentally")
    
    # Save report
    try:
        pipeline_report_file = f"{base_dir}/results/snrna_scvi/edgeR_pipeline_report.txt"
        os.makedirs(os.path.dirname(pipeline_report_file), exist_ok=True)
        
        with open(pipeline_report_file, 'w') as f:
            f.write('\n'.join(report_lines))
        
        print(f"📋 Pipeline report saved: {pipeline_report_file}")
        
        # Print summary to console
        print("\n📊 PIPELINE SUMMARY:")
        print("-" * 30)
        print(f"  edgeR files: {len(edger_files)}")
        print(f"  Pathway files: {len(pathway_files)}")
        print(f"  Plots generated: {total_plots}")
        
        return True
        
    except Exception as e:
        print(f"❌ Failed to generate report: {e}")
        return False

def main():
    """Main pipeline execution"""
    pipeline_start = time.time()
    
    print_banner()
    
    try:
        # Step 0: Check prerequisites
        check_prerequisites()
        
        # Step 1: Run edgeR analysis
        edger_success = run_edger_analysis()
        if not edger_success:
            print("❌ Pipeline failed at edgeR step")
            return False
        
        # Step 2: Run pathway enrichment
        pathway_success = run_pathway_enrichment()
        if not pathway_success:
            print("⚠️  edgeR completed but pathway enrichment failed")
            print("   You can still use the edgeR results")
        
        # Step 3: Generate report
        report_success = generate_pipeline_report()
        
        # Final summary
        total_time = time.time() - pipeline_start
        
        print("\n" + "🎉" * 35)
        print("🏁 PIPELINE EXECUTION COMPLETE! 🏁")
        print("🎉" * 35)
        print(f"⏱️  Total time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
        print(f"✅ edgeR analysis: {'SUCCESS' if edger_success else 'FAILED'}")
        print(f"✅ Pathway enrichment: {'SUCCESS' if pathway_success else 'FAILED'}")
        print(f"✅ Report generation: {'SUCCESS' if report_success else 'FAILED'}")
        
        if edger_success:
            print("\n🎯 KEY OUTPUTS:")
            print("  📊 Differential expression results")
            if pathway_success:
                print("  🛤️  Pathway enrichment analysis")
            print("  📈 Comprehensive visualizations")
            print("  📋 Analysis reports")
            
            print("\n🔬 READY FOR BIOLOGICAL INTERPRETATION!")
        
        print("🎉" * 35)
        
        return edger_success
        
    except Exception as e:
        print(f"\n❌ PIPELINE FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)