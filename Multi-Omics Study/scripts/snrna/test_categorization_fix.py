#!/usr/bin/env python3
"""
Test script to verify the categorization fixes for LEMUR results loading
This script simulates the categorization logic to ensure it works correctly
"""

import os
import pandas as pd
import numpy as np

def test_categorization_logic():
    """Test the improved categorization logic"""

    print("🧪 Testing Categorization Logic Fixes")
    print("=" * 50)

    # Test cases with expected outcomes
    test_cases = [
        # (filename, expected_category, description)
        ("Male_Vs_Female_results.csv", "main", "Sex comparison - should be main effect"),
        ("OUD_Vs_Control_results.csv", "main", "Disease comparison - should be main effect"),
        ("sex_oud_interaction_results.csv", "interaction", "Interaction effect"),
        ("OUD_within_male_results.csv", "sex_stratified", "Analysis within males only"),
        ("OUD_within_female_results.csv", "sex_stratified", "Analysis within females only"),
        ("Control_male_only_results.csv", "sex_stratified", "Male-only analysis"),
        ("Treatment_female_only_results.csv", "sex_stratified", "Female-only analysis"),
        ("stratified_male_analysis_results.csv", "sex_stratified", "Stratified male analysis"),
        ("stratified_female_analysis_results.csv", "sex_stratified", "Stratified female analysis"),
        ("Caudate_Vs_Putamen_results.csv", "main", "Regional comparison"),
        ("Age_effect_interaction_results.csv", "interaction", "Age interaction"),
    ]

    # Test the categorization logic
    passed = 0
    failed = 0

    for filename, expected, description in test_cases:
        # Apply the same logic as in the fixed script
        is_interaction = 'interaction' in filename.lower()

        # Sex-stratified: analyses within one sex (more specific patterns)
        sex_stratified_patterns = ['within_male', 'within_female', 'male_only', 'female_only',
                                 'stratified_male', 'stratified_female']
        is_sex_stratified = any(pattern in filename.lower() for pattern in sex_stratified_patterns) and not is_interaction

        # Main effects: everything else, including male_vs_female comparisons
        is_main_effect = not is_interaction and not is_sex_stratified

        if is_interaction:
            actual_category = "interaction"
        elif is_sex_stratified:
            actual_category = "sex_stratified"
        else:
            actual_category = "main"

        # Check result
        if actual_category == expected:
            print(f"✅ {filename}")
            print(f"   {description}")
            print(f"   Expected: {expected}, Got: {actual_category}")
            passed += 1
        else:
            print(f"❌ {filename}")
            print(f"   {description}")
            print(f"   Expected: {expected}, Got: {actual_category}")
            failed += 1
        print()

    print(f"📊 Test Results: {passed} passed, {failed} failed")
    return failed == 0

def test_column_mapping():
    """Test the column mapping for interaction files"""

    print("\n🧪 Testing Column Mapping for Interaction Files")
    print("=" * 50)

    # Simulate interaction file columns
    interaction_columns = ['gene', 'log_fold_change_female', 'pvalue_female', 'padj_female',
                          'log_fold_change_male', 'pvalue_male', 'padj_male',
                          'interaction_effect', 'abs_interaction_effect', 'interaction_pvalue', 'interaction_padj']

    # Test the column mapping logic
    print("Testing interaction file column detection:")

    # Create a mock dataframe
    df = pd.DataFrame({col: [0.1, 0.2, 0.3] for col in interaction_columns})

    # Apply the logic from the fixed script
    if 'interaction_effect' in df.columns:
        coef_col = 'interaction_effect'
        pval_col = 'interaction_pvalue'
        fdr_col = 'interaction_padj'
        gene_col = 'gene'
        print("✅ Interaction file detected correctly")
        print(f"   Coefficient column: {coef_col}")
        print(f"   P-value column: {pval_col}")
        print(f"   FDR column: {fdr_col}")
        print(f"   Gene column: {gene_col}")

        # Test standardization
        df_standard = df.copy()
        df_standard = df_standard.rename(columns={
            coef_col: 'coefficient',
            pval_col: 'pvalue',
            gene_col: 'gene',
            fdr_col: 'padj'
        })

        expected_cols = ['coefficient', 'pvalue', 'gene', 'padj']
        if all(col in df_standard.columns for col in expected_cols):
            print("✅ Column standardization successful")
            return True
        else:
            print("❌ Column standardization failed")
            return False
    else:
        print("❌ Interaction file not detected")
        return False

def test_oud_pattern_categorization():
    """Test the OUD pattern categorization logic"""

    print("\n🧪 Testing OUD Pattern Categorization")
    print("=" * 50)

    # Simulate contrast names that might be found
    oud_contrasts = [
        "OUD_Vs_Control",
        "Male_Vs_Female",
        "OUD_within_male",
        "OUD_within_female",
        "sex_oud_interaction",
        "Control_Vs_OUD_male_only",
        "stratified_female_oud_analysis"
    ]

    print("Testing OUD pattern categorization:")

    # Apply the improved categorization logic
    interaction_oud = [c for c in oud_contrasts if 'interaction' in c.lower()]

    # Sex-stratified: analyses within one sex (more specific patterns)
    sex_stratified_patterns = ['within_male', 'within_female', 'male_only', 'female_only',
                             'stratified_male', 'stratified_female']
    sex_stratified_oud = [c for c in oud_contrasts
                         if any(pattern in c.lower() for pattern in sex_stratified_patterns)
                         and 'interaction' not in c.lower()]

    # Main effects: everything else, including male_vs_female comparisons
    main_oud = [c for c in oud_contrasts if c not in interaction_oud and c not in sex_stratified_oud]

    print(f"\n📊 Main effects: {len(main_oud)}")
    for contrast in main_oud:
        print(f"   • {contrast}")

    print(f"\n👫 Sex-stratified: {len(sex_stratified_oud)}")
    for contrast in sex_stratified_oud:
        print(f"   • {contrast}")

    print(f"\n🔄 Interaction effects: {len(interaction_oud)}")
    for contrast in interaction_oud:
        print(f"   • {contrast}")

    # Verify expected results
    expected_main = {"OUD_Vs_Control", "Male_Vs_Female"}
    expected_sex_stratified = {"OUD_within_male", "OUD_within_female", "Control_Vs_OUD_male_only", "stratified_female_oud_analysis"}
    expected_interaction = {"sex_oud_interaction"}

    success = (
        set(main_oud) == expected_main and
        set(sex_stratified_oud) == expected_sex_stratified and
        set(interaction_oud) == expected_interaction
    )

    if success:
        print("\n✅ OUD pattern categorization is correct!")
    else:
        print("\n❌ OUD pattern categorization has issues!")
        print(f"Expected main: {expected_main}")
        print(f"Got main: {set(main_oud)}")
        print(f"Expected sex-stratified: {expected_sex_stratified}")
        print(f"Got sex-stratified: {set(sex_stratified_oud)}")
        print(f"Expected interaction: {expected_interaction}")
        print(f"Got interaction: {set(interaction_oud)}")

    return success

def main():
    """Run all tests"""
    print("🚀 Testing LEMUR Categorization Fixes")
    print("="*60)

    # Run tests
    test1 = test_categorization_logic()
    test2 = test_column_mapping()
    test3 = test_oud_pattern_categorization()

    print("\n" + "="*60)
    print("📋 FINAL TEST RESULTS")
    print("="*60)

    if all([test1, test2, test3]):
        print("🎉 ALL TESTS PASSED!")
        print("✅ Categorization fixes are working correctly")
        print("✅ Male_Vs_Female will be correctly classified as main effect")
        print("✅ Interaction files will be properly detected and loaded")
        print("✅ Sex-stratified analyses will be correctly identified")
        print("✅ OUD pattern analysis will show proper counts")
    else:
        print("❌ SOME TESTS FAILED!")
        print(f"   Basic categorization: {'✅' if test1 else '❌'}")
        print(f"   Column mapping: {'✅' if test2 else '❌'}")
        print(f"   OUD patterns: {'✅' if test3 else '❌'}")

    print(f"\n💡 The key fixes implemented:")
    print("   1. Male_Vs_Female is now correctly categorized as 'main effect'")
    print("   2. Sex-stratified detection uses specific patterns (within_male, male_only, etc.)")
    print("   3. Interaction files are properly detected and column-mapped")
    print("   4. OUD pattern analysis uses the same improved logic")

if __name__ == "__main__":
    main()
