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


rule plot_betascan_results:
    input:
        scores = rules.estimate_b1.output.scores,
    output:
        plot = "results/plots/b1.scores.png",
    params:
        cutoff1 = 43.14256,
        cutoff2 = 18.62207,
    shell:
        """
        Rscript workflow/scripts/manhattan.R {input.scores} {output.plot} {params.cutoff1} {params.cutoff2}
        """
