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


import matplotlib
matplotlib.use('Agg')


rule process_maladapt_predictions:
    input:
        predictions = rules.run_maladapt_prediction.output.predictions,
    output:
        predictions = "results/maladapt/6window_nea_LIT.processed.predictions",
    run:
        import pandas as pd
        df = pd.read_csv(input.predictions, sep=",")
        df.rename(columns={'Pr{AI}': 'prob1'}, inplace=True)
        df = df[["start", "end", "prob1"]]
        df.to_csv(output.predictions, index=False)


rule calculate_gene_statistics:
    input:
        predictions = rules.process_maladapt_predictions.output.predictions,
        goi = "config/hla_genes_paper.txt",
    output:
        gene_stats = "results/plots/maladapt_gene_prediction_probability_stats.csv",
    script:
        "../scripts/get_gene_prediction_probability.py"


rule plot_maladapt_histogram:
    input:
        df = rules.process_maladapt_predictions.output.predictions,
        stats = rules.calculate_gene_statistics.output.gene_stats,
    output:
        plot = "results/plots/maladapt_histogram.pdf",
    params:
        title = "Distributions of Adaptive Introgression Probabilities on Chromosome 6 in Lithuanian Genomes \n Using the Pretrained MaLAdapt_25_-sweep-all_model Model",
    script:
        "../scripts/plot_histogram.py"
