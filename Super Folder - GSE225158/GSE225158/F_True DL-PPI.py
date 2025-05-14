import pandas as pd

# Process male results
male_df = pd.read_csv("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/True M_DL-PPI outputs/dl_ppi_predictions.csv")
male_viz = pd.DataFrame({
    'gene1': male_df['acc1'],
    'gene2': male_df['acc2'],
    'prediction': male_df['prediction'],
    'method': male_df['method']
})
male_filtered = male_viz[male_viz['prediction'] >= 0.4]
print(f"Male interactions at 0.4 threshold: {len(male_filtered)}")
male_filtered.to_csv("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/dl_ppi_viz_data_male_0.4.csv", index=False)

# Process female results
female_df = pd.read_csv("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/True F_DL-PPI outputs - 0.5/dl_ppi_predictions.csv")
female_viz = pd.DataFrame({
    'gene1': female_df['acc1'],
    'gene2': female_df['acc2'],
    'prediction': female_df['prediction'],
    'method': female_df['method']
})
female_filtered = female_viz[female_viz['prediction'] >= 0.4]
print(f"Female interactions at 0.4 threshold: {len(female_filtered)}")
female_filtered.to_csv("/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE225158/DL-PPI outputs/dl_ppi_viz_data_female_0.4.csv", index=False)