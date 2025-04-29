import pandas as pd
'''
This script subsets the `gene_library` DataFrame to make a new `.csv` file
    - The new file only contains the unique gene names
'''

# Load the gene library
gene_library = pd.read_csv('gene_library.csv')

# Subset the gene library to only include unique gene names
gene_library_subset = gene_library[['Gene Name']].drop_duplicates()

# Save the subsetted gene library to a new file
gene_library_subset.to_csv('gene_names.csv', index=False)