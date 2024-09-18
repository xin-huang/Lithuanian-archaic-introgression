# GNU General Public License v3.0
# Copyright 2024 Josef Hackl and Xin Huang
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


rule plot_betascan_results:
    input:
        scores = rules.estimate_b1.output.scores,
    output:
        plot = "results/plots/lit.chr6.imputed.biallelic.snps.b1.scores.pdf",
    params:
        cutoff = 71.572528,
    shell:
        """
        Rscript workflow/scripts/manhattan.R {input.scores} {output.plot} {params.cutoff}
        """


rule process_predictions:
    input:
        predictions="eurref/Nea_to_CEU_af-0.05_2250018620/predictions.txt"
    output:
        df="processed_data.csv"
    run:
        import pandas as pd
        df = pd.read_csv(input.predictions, sep="\t")
        df.rename(columns={'Pr{AI}': 'prob1'}, inplace=True)
        df = df[["start", "end", "prob1"]]
        df.to_csv(output.df, index=False)

rule read_genes_of_interest:
    input:
        "hla_genes_paper.txt"
    output:
        goi="processed_goi.csv"
    run:
        import pandas as pd
        goi = pd.read_csv(input[0], header=None)
        goi.columns = ["gene", "start", "end"]
        goi.to_csv(output.goi, index=False)

rule calculate_overlap_intervals:
    input:
        df="processed_data.csv"
    output:
        overlap_df="overlap_intervals.csv"
    run:
        import pandas as pd
        df = pd.read_csv(input.df)

        window_size = 100000
        stride = 20000
        overall_start = df['start'].min()
        overall_end = df['end'].max()

        overlap_intervals = []
        current_start = overall_start

        while current_start <= overall_end:
            current_end = current_start + stride

            overlapping_windows = df[(df['start'] <= current_end) & (df['end'] >= current_start)]
            if not overlapping_windows.empty:
                average_prob = overlapping_windows['prob1'].mean()
                overlap_intervals.append({'chrom': 6, 'start': current_start, 'end': current_end, 'average_prob1': average_prob})

            current_start = current_start + stride

        overlap_df = pd.DataFrame(overlap_intervals)
        overlap_df.to_csv(output.overlap_df, index=False)

rule calculate_gene_statistics:
    input:
        overlap_df="overlap_intervals.csv",
        goi="processed_goi.csv"
    output:
        gene_stats="gene_stats.csv"
    run:
        import pandas as pd

        def get_gene_df(odf, goi, prob_column="average_prob1", gene_column="gene", only_full=False):
            additional_df = []
            odf = odf.copy()

            for _, row in goi.iterrows():
                overlapping_windows = odf[(row["start"] <= odf['end'].astype(int)) & (row["end"] >= odf['start'].astype(int))]

                average_prob = overlapping_windows[prob_column].mean()
                min_prob = overlapping_windows[prob_column].min()
                max_prob = overlapping_windows[prob_column].max()
                std_prob = overlapping_windows[prob_column].std()

                additional_df.append([overlapping_windows, average_prob, min_prob, max_prob, std_prob, row["gene"]])

            return additional_df

        odf = pd.read_csv(input.overlap_df)
        goi = pd.read_csv(input.goi)
        adf = get_gene_df(odf, goi)

        goi["length"] = goi["end"] - goi["start"]
        adf_ranges = [x[1:-1] for x in adf]

        names = [x[-1] for x in adf]
        means = [x[1] for x in adf]
        mins = [x[2] for x in adf]
        maxs = [x[3] for x in adf]
        stds = [x[4] for x in adf]

        gene_stats = pd.DataFrame({
            "gene": names,
            "mean": means,
            "min": mins,
            "max": maxs,
            "std": stds
        })

        gene_stats.to_csv(output.gene_stats, index=False)

rule plot_max_values:
    input:
        df="processed_data.csv",
        stats="gene_stats.csv"
    output:
        "max_genomatnn_005.png"
    run:
        import pandas as pd
        import matplotlib.pyplot as plt

        df = pd.read_csv(input.df)
        stats = pd.read_csv(input.stats)

        def return_percentage(df, given_probability=0.8):
            higher_count = df[df['prob1'] > given_probability].shape[0]
            total_rows = df.shape[0]
            percentage_higher = (higher_count / total_rows) * 100
            return percentage_higher

        maxs = stats["max"].tolist()
        names = stats["gene"].tolist()
        max_percs = [return_percentage(df, maxv) for maxv in maxs]

        colors = ['red', 'green', 'blue', 'yellow', 'orange', 'magenta']
        default_color = 'grey'

        n, bins, _ = plt.hist(df["prob1"], bins=30, edgecolor='black', alpha=0)
        plt.hist(df["prob1"], bins=bins, edgecolor='black', color=default_color, alpha=0.7, label='predicted values of windows')

        for max_value, maxp, name, color in zip(maxs, max_percs, names, colors):
            plt.axvline(x=max_value, color=color, linestyle='--', linewidth=2, label=f'{name}, max: {max_value:.3f}, top: {maxp:.3f} %')

        plt.xlabel('Value')
        plt.ylabel('Frequency')
        plt.xlim(0)
        plt.title('Genomatnn, chr 6, 0.05 model, beagle with reference \n Lines indicate max values')
        plt.legend()
        plt.savefig(output[0], format='png', dpi=300)
        plt.show()

rule plot_mean_values:
    input:
        overlap_df="overlap_intervals.csv",
        stats="gene_stats.csv"
    output:
        "mean_genomatnn_005.png"
    run:
        import pandas as pd
        import matplotlib.pyplot as plt

        overlap_df = pd.read_csv(input.overlap_df)
        stats = pd.read_csv(input.stats)

        def return_percentage(df, given_probability=0.8):
            higher_count = df[df['average_prob1'] > given_probability].shape[0]
            total_rows = df.shape[0]
            percentage_higher = (higher_count / total_rows) * 100
            return percentage_higher

        means = stats["mean"].tolist()
        names = stats["gene"].tolist()
        mean_percs = [return_percentage(overlap_df, meanv) for meanv in means]

        colors = ['red', 'green', 'blue', 'yellow', 'orange', 'magenta']
        default_color = 'grey'

        n, bins, _ = plt.hist(overlap_df["average_prob1"], bins=30, edgecolor='black', alpha=0)
        plt.hist(overlap_df["average_prob1"], bins=bins, edgecolor='black', color=default_color, alpha=0.7, label='averaged predicted values of windows')

        for mean_value, meanp, name, color in zip(means, mean_percs, names, colors):
            plt.axvline(x=mean_value, color=color, linestyle='--', linewidth=2, label=f'{name}, mean: {mean_value:.3f}, top: {meanp:.3f} %')

        plt.xlabel('Value')
        plt.ylabel('Frequency')
        plt.xlim(0)
        plt.title('Genomatnn, chr 6, 0.05 model, beagle with reference \n Histogram shows averaged fragments of length 20k \n Lines indicate mean values')
        plt.legend()
        plt.savefig(output[0], format='png', dpi=300)
        plt.show()
