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


rule test_hwe:
    input:
        vcf = rules.extract_biallelic_snps.output.vcf,
    output:
        hwe_outliers = "results/processed_data/Lithuanians/lit.chr6.imputed.biallelic.snps.hwe.outliers",
    shell:
        """
        plink --vcf {input.vcf} --hardy --out results/processed_data/Lithuanians/lit.chr6.imputed.biallelic.snps --set-missing-var-ids @:# --const-fid
        awk '$7>$8' results/processed_data/Lithuanians/lit.chr6.imputed.biallelic.snps.hwe | \
        sed '1d' | \
        awk '$9<0.001{{print $2}}' | \
        awk -F ":" '{{print $1"\\t"$2-1"\\t"$2}}' > {output.hwe_outliers}
        """


rule convert_repeat_files:
    input:
        rmsk = rules.download_repeats.output.rmsk,
        segdup = rules.download_repeats.output.segdup,
        simrep = rules.download_repeats.output.simrep,
    output:
        rmsk = "results/processed_data/repeats/hg19.rmsk.autosomes.bed",
        segdup = "results/processed_data/repeats/hg19.seg.dups.autosomes.bed",
        simrep = "results/processed_data/repeats/hg19.simple.repeats.autosomes.bed",
    shell:
        """
        zcat {input.rmsk} | awk 'BEGIN{{OFS="\\t"}}$6!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $6,$7,$8,$11,$2,$10}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.rmsk}
        zcat {input.segdup} | awk 'BEGIN{{OFS="\\t"}}$2!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $2,$3,$4,$5,$6,$7}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.segdup}
        zcat {input.simrep} | awk 'BEGIN{{OFS="\\t"}}$2!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $2,$3,$4,$5,$11}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.simrep}
        """


rule get_allele_counts:
    input:
        vcf = rules.extract_biallelic_snps.output.vcf,
        hwe_outliers = rules.test_hwe.output.hwe_outliers,
        rmsk = rules.convert_repeat_files.output.rmsk,
        segdup = rules.convert_repeat_files.output.segdup,
        simrep = rules.convert_repeat_files.output.simrep,
    output:
        ac = "results/betascan/lit.chr6.imputed.biallelic.snps.ac",
    params:
        maf1 = 100*0.05,
        maf2 = 100*0.95,
    shell:
        """
        if [ -s {input.hwe_outliers} ]; then \
            bcftools view -T ^{input.rmsk} {input.vcf} | \
            bcftools view -T ^{input.simrep} | \
            bcftools view -T ^{input.segdup} | \
            bcftools view -T ^{input.hwe_outliers} | \
            bcftools query -f "%POS\\t%INFO/AC\\t%INFO/AN\\n" | \
            awk '($2>{params.maf1})&&($2<{params.maf2})' > {output.ac}; \
        else \
            bcftools view -T ^{input.rmsk} {input.vcf} | \
            bcftools view -T ^{input.simrep} | \
            bcftools view -T ^{input.segdup} | \
            bcftools query -f "%POS\\t%INFO/AC\\t%INFO/AN\\n" | \
            awk '($2>{params.maf1})&&($2<{params.maf2})' > {output.ac}; \
        fi
        """


rule estimate_b1:
    input:
        ac = rules.get_allele_counts.output.ac,
        flag = rules.download_betascan.output.download_flag,
    output:
        scores = "results/betascan/lit.chr6.imputed.biallelic.snps.b1.scores",
    shell:
        """
        python resources/tools/BetaScan/BetaScan.py -i {input.ac} -m .15 -fold | grep -v Position | awk -v chr=6 '{{print chr"\\t"$0}}' | awk '{{print $1":"$2"\\t"$1"\\t"$2"\\t"$3}}' | sed '1iSNP\\tCHR\\tBP\\tB1'> {output.scores}
        """


rule b1_summary:
    input:
        scores = rules.estimate_b1.output.scores,
    output:
        summary = "results/betascan/lit.chr6.imputed.biallelic.snps.b1.scores.summary",
    shell:
        """
        set +o pipefail
        echo "#top 0.05% cutoff" > {output.summary}
        sed '1d' {input.scores} | sort -rnk 4,4 | head -55 | tail -1 >> {output.summary}
        echo "#top 1% cutoff" >> {output.summary}
        sed '1d' {input.scores} | sort -rnk 4,4 | head -1111 | tail -1 >> {output.summary}
        echo "#HLA-B; HLA-C highest score" >> {output.summary}
        awk '$3>=31198205&&$3<=31348022' {input.scores} | sort -rnk 4,4 | head -1 >> {output.summary}
        echo "#HLA-DQA1; HLA-DQB1 highest score" >> {output.summary}
        awk '$3>=32598042&&$3<=32644388' {input.scores} | sort -rnk 4,4 | head -1 >> {output.summary}
        echo "#HLA-DQA2; HLA-DQB2 highest score" >> {output.summary}
        awk '$3>=32698044&&$3<=32748039' {input.scores} | sort -rnk 4,4 | head -1 >> {output.summary}
        """
