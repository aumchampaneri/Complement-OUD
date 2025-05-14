import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

results_dir = '/Users/aumchampaneri/PycharmProjects/Complement-OUD/GSE289002/DE_results/sex_treatment'
plots_dir = os.path.join(results_dir, 'plots')
os.makedirs(plots_dir, exist_ok=True)

for fname in os.listdir(results_dir):
    if fname.endswith('.csv'):
        contrast = os.path.splitext(fname)[0]
        csv_path = os.path.join(results_dir, fname)
        df = pd.read_csv(csv_path, index_col=0)

        # Volcano plot
        plt.figure(figsize=(8,6))
        sns.scatterplot(
            x='logFC', y=-np.log10(df['FDR']),
            data=df, hue=(df['FDR'] < 0.05) & (abs(df['logFC']) > 1),
            palette={True: 'red', False: 'grey'}, legend=False, alpha=0.7
        )
        plt.axhline(-np.log10(0.05), color='blue', linestyle='--')
        plt.axvline(1, color='green', linestyle='--')
        plt.axvline(-1, color='green', linestyle='--')
        plt.xlabel('log2 Fold Change')
        plt.ylabel('-log10(FDR)')
        plt.title(f'Volcano Plot: {contrast}')
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, f'{contrast}_volcano.png'))
        plt.close()

        # MA plot
        plt.figure(figsize=(8,6))
        sns.scatterplot(
            x='logCPM', y='logFC',
            data=df, hue=(df['FDR'] < 0.05),
            palette={True: 'red', False: 'grey'}, legend=False, alpha=0.7
        )
        plt.axhline(0, color='black', linestyle='--')
        plt.xlabel('logCPM')
        plt.ylabel('log2 Fold Change')
        plt.title(f'MA Plot: {contrast}')
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, f'{contrast}_MA.png'))
        plt.close()