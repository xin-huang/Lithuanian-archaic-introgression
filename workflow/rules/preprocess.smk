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


rule plink_to_vcf:
    input:
        ped = rules.download_lithuanian_genomes.output.ped,
        map = rules.download_lithuanian_genomes.output.map,
    output:
        vcf = "results/processed_data/Lithuanians/lit.vcf.gz",
    shell:
        """
        plink --ped {input.ped} --map {input.map} --recode vcf --out results/processed_data/Lithuanians/lit
        bgzip -c results/processed_data/Lithuanians/lit.vcf > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule extract_chr6:
    input:
        vcf = rules.plink_to_vcf.output.vcf,
    output:
        vcf = "results/processed_data/Lithuanians/lit.chr6.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -r 6 | bgzip -c > {output.vcf}
        """


rule align_vcf_to_reference:
    input:
        vcf = rules.extract_chr6.output.vcf,
        reference = rules.download_1KG_genomes.output.ref,
    output:
        vcf = "results/processed_data/Lithuanians/lit.chr6.aligned.vcf.gz",
    script:
        "../scripts/align_reference.py"


rule tabix_index:
    input:
        vcf = rules.align_vcf_to_reference.output.vcf,
    output:
        vcf = "results/processed_data/Lithuanians/lit.chr6.aligned.vcf.gz.tbi",
    shell:
        """
        tabix -p vcf {input.vcf}
        """


rule extract_yoruba_eur_samples:
    input:
        samples = rules.download_1KG_genomes.output.samples,
    output:
        yri_samples = "results/processed_data/1KG/yri_samples.txt",
        eur_samples = "results/processed_data/1KG/european_samples.txt",
    shell:
        """
        grep "YRI" {input.samples} | awk '{{print $1}}' > {output.yri_samples}
        grep -E "CEU|TSI|FIN|GBR|IBS" {input.samples} | awk '{{print $1}}' > {output.eur_samples}
        """


rule filter_for_eur:
    input:
        vcf = rules.download_1KG_genomes.output.vcf,
        samples = rules.extract_yoruba_eur_samples.output.eur_samples,
    output:
        vcf = "results/processed_data/1KG/chr6.EUR.vcf.gz",
    shell:
        """
        bcftools view -S {input.samples} -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """


rule filter_for_yri:
    input:
        vcf = rules.download_1KG_genomes.output.vcf,
        samples = rules.extract_yoruba_eur_samples.output.yri_samples,
    output:
        vcf = "results/processed_data/1KG/chr6.YRI.vcf.gz",
    shell:
        """
        bcftools view -S {input.samples} -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """


rule beagle_imputation_with_ref:
    input:
        vcf = rules.align_vcf_to_reference.output.vcf,
        ref = rules.filter_for_eur.output.vcf,
        flag = rules.download_beagle.output.download_flag,
    output:
        vcf = "results/processed_data/Lithuanians/lit.chr6.imputed.vcf.gz",
    resources:
        cpus = 2,
    shell:
        """
        java -Xmx4g -jar resources/tools/beagle/beagle.06Aug24.a91.jar gt={input.vcf} ref={input.ref} out=results/processed_data/Lithuanians/lit.chr6.imputed nthreads={resources.cpus} ne=10000
        tabix -p vcf {output.vcf}
        """


rule extract_biallelic_snps:
    input:
        vcf = rules.beagle_imputation_with_ref.output.vcf,
    output:
        vcf = "results/processed_data/Lithuanians/lit.chr6.imputed.biallelic.snps.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -v snps -m 2 -M 2 | bcftools view -c 1:minor | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
