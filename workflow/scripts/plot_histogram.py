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
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import pandas as pd
import matplotlib.pyplot as plt


def return_percentage(df, given_probability=0.8):
    higher_count = df[df['prob1'] > given_probability].shape[0]
    total_rows = df.shape[0]
    percentage_higher = (higher_count / total_rows) * 100
    return percentage_higher


def plot_stat_distribution(ax, stats, df, stat='max', colors=None, default_color='grey', bins=30, title=None):
    # Extract the selected statistic and gene names from the stats DataFrame
    values = stats[stat].tolist()
    names = stats["gene"].tolist()
    percentages = [return_percentage(df, val) for val in values]

    # Set colors
    if colors is None:
        colors = ['red', 'green', 'blue', 'yellow', 'orange', 'magenta']

    # Create histogram of the 'prob1' column
    n, bins, _ = ax.hist(df["prob1"], bins=bins, edgecolor='black', alpha=0)
    ax.hist(df["prob1"], bins=bins, edgecolor='black', color=default_color, alpha=0.7, label='predicted probabilities of windows')

    # Add vertical lines for each statistic value
    for value, perc, name, color in zip(values, percentages, names, colors):
        ax.axvline(x=value, color=color, linestyle='--', linewidth=2, label=f'{name}, {stat}: {value:.3f}, top: {perc:.3f} %')

    # Set plot labels and title
    ax.set_xlabel('Probability')
    ax.set_ylabel('Count')
    ax.set_xlim(0)
    
    # Set individual subplot title
    if title:
        ax.set_title(title)
    
    ax.legend()


df = pd.read_csv(snakemake.input.df)
stats = pd.read_csv(snakemake.input.stats)

# Create subplots
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# Plot 'max' in the first subplot with a title
plot_stat_distribution(axs[0], stats, df, stat='max', title='Max Probabilities of HLA Genes')

# Plot 'mean' in the second subplot with a title
plot_stat_distribution(axs[1], stats, df, stat='mean', title='Mean Probabilities of HLA Genes')

# Adjust layout and save the figure
plt.suptitle(snakemake.params.title)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(snakemake.output.plot)
plt.close()
