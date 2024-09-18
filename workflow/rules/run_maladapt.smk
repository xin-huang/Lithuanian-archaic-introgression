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


rule merge_lit_yri:
    input:
        vcf1 = rules.beagle_imputation_with_ref.output.vcf,
        vcf2 = rules.filter_for_yri.output.vcf,
    output:
        vcf = "results/maladapt/vcf_modern/merged_6_output.vcf.gz",
    shell:
        """
        bcftools isec -n=2 -p dir {input.vcf1} {input.vcf2}
        bgzip -c dir/0000.vcf > dir/0000.vcf.gz
        tabix -p vcf dir/0000.vcf.gz
        bgzip -c dir/0001.vcf > dir/0001.vcf.gz
        tabix -p vcf dir/0001.vcf.gz
        bcftools merge dir/0000.vcf.gz dir/0001.vcf.gz | bcftools view -v snps -m 2 -M 2 -c 1:minor | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        rm -r dir
        """


rule merge_lit_yri_nea:
    input:
        vcf1 = rules.merge_lit_yri.output.vcf,
        vcf2 = rules.download_nea_genome.output.vcf,
    output:
        vcf = "results/maladapt/merged_6_output_with_nea.vcf.gz",
    shell:
        """
        bcftools isec -n=2 -p dir {input.vcf1} {input.vcf2}
        bgzip -c dir/0000.vcf > dir/0000.vcf.gz
        tabix -p vcf dir/0000.vcf.gz
        bgzip -c dir/0001.vcf > dir/0001.vcf.gz
        tabix -p vcf dir/0001.vcf.gz
        bcftools merge -Oz -o {output.vcf} dir/0000.vcf.gz dir/0001.vcf.gz
        tabix -p vcf {output.vcf}
        rm -r dir
        """


rule create_stepsize_file:
    input:
        vcf = rules.merge_lit_yri.output.vcf,
    output:
        stepsize_file = "results/maladapt/window_range_custom/CHR6_windows_overlap.txt",
    shell:
        """
        python workflow/scripts/window_range_creator_stepsize.py {input.vcf} {output.stepsize_file} 50000
        """


rule get_sample_names:
    input:
        yri_samples = rules.extract_yoruba_eur_samples.output.yri_samples,
        lit_vcf = rules.beagle_imputation_with_ref.output.vcf,
    output:
        yri_samples = "results/maladapt/dir_pop/YRI.txt",
        lit_samples = "results/maladapt/dir_pop/LIT.txt",
    shell:
        """
        cp {input.yri_samples} {output.yri_samples}
        bcftools query -l {input.lit_vcf} > {output.lit_samples}
        """


rule run_maladapt_feature_extraction:
    input:
        merged_vcf = rules.merge_lit_yri.output.vcf,
        archaic_vcf = rules.download_nea_genome.output.vcf,
        yri_samples = rules.get_sample_names.output.yri_samples,
        lit_samples = rules.get_sample_names.output.lit_samples,
        overlap = rules.create_stepsize_file.output.stepsize_file,
    output:
        feature = "results/maladapt/dir_stat_out/6window_nea_LIT.txt",
    resources:
        time = 43200, mem_gb = 100, partition = 'basic',    
    conda: "../envs/maladapt-env.yaml",
    shell:
        """
        mkdir -p results/maladapt/dir_stat_out
        mkdir -p results/maladapt/dir_out
        mkdir -p results/maladapt/vcf_modern_out
        mkdir -p results/maladapt/vcf_archaic_out
        mkdir -p results/maladapt/vcf_archaic
        cp {input.archaic_vcf} results/maladapt/vcf_archaic/
        python workflow/scripts/1empirical_compute-features_custom_chr6_final.py
        """
