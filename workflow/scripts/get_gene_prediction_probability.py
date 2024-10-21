# Copyright 2024 Josef Hackl and Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import pandas as pd


# Read the windows with probabilities
windows_df = pd.read_csv(snakemake.input.predictions)

# Read the genes with their start and end positions
genes_df = pd.read_csv(snakemake.input.goi, names=['gene', 'start', 'end'])

def get_gene_stats(windows_df, genes_df):
    additional_df = []

    # Iterate over each gene
    for _, row in genes_df.iterrows():
        gene_name = row['gene']
        gene_start = row['start']
        gene_end = row['end']

        # Find overlapping windows
        overlapping_windows = windows_df[(windows_df['start'] <= gene_end) & (windows_df['end'] >= gene_start)]

        # Calculate statistics for prob1
        if not overlapping_windows.empty:
            mean_prob = overlapping_windows['prob1'].mean()
            min_prob = overlapping_windows['prob1'].min()
            max_prob = overlapping_windows['prob1'].max()
            std_prob = overlapping_windows['prob1'].std()

            additional_df.append([gene_name, mean_prob, min_prob, max_prob, std_prob])
        else:
            # If no overlapping windows, append NaNs
            additional_df.append([gene_name, None, None, None, None])

    # Create a DataFrame to store the results
    result_df = pd.DataFrame(additional_df, columns=['gene', 'mean', 'min', 'max', 'std'])

    return result_df

# Call the function
gene_stats_df = get_gene_stats(windows_df, genes_df)

# Save the results to a CSV file
gene_stats_df.to_csv(snakemake.output.gene_stats, index=False)
