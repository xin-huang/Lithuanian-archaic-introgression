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


rule create_stepsize_file:
    input:
        vcf = rules.merge_lit_yri.output.vcf,
    output:
        stepsize_file = "results/maladapt/window_range_custom/CHR6_windows_overlap.txt",
    shell:
        """
        python workflow/scripts/window_range_creator_stepsize.py {input.vcf} {output.stepsize_file} 50000
        """


rule split_stepsize_file:
    input:
        stepsize_file = rules.create_stepsize_file.output.stepsize_file,
    output:
        expand("results/maladapt/window_range_custom/CHR6_windows_overlap_{n}.txt", n = range(1, 1001))
    params:
        output_dir = "results/maladapt/window_range_custom/",
        output_prefix = "CHR6_windows_overlap",
        nfiles = 1000,
    script:
        "../scripts/split_window_ranges.py"


rule get_sample_names:
    input:
        yri_samples = rules.extract_yoruba_eur_samples.output.yri_samples,
        lit_vcf = rules.beagle_imputation_with_ref.output.vcf,
        archaic_vcf = rules.download_nea_genome.output.vcf,
    output:
        yri_samples = "results/maladapt/dir_pop/YRI.txt",
        lit_samples = "results/maladapt/dir_pop/LIT.txt",
    params:
        vcf_archaic = "results/maladapt/vcf_archaic/",
    shell:
        """
        cp {input.yri_samples} {output.yri_samples}
        bcftools query -l {input.lit_vcf} > {output.lit_samples}
        mkdir -p {params.vcf_archaic}
        cp {input.archaic_vcf} {params.vcf_archaic}
        """


rule run_maladapt_feature_extraction:
    input:
        merged_vcf = rules.merge_lit_yri.output.vcf,
        yri_samples = rules.get_sample_names.output.yri_samples,
        lit_samples = rules.get_sample_names.output.lit_samples,
        overlap = "results/maladapt/window_range_custom/CHR6_windows_overlap_{n}.txt",
    output:
        feature = "results/maladapt/dir_stat_out_{n}/6window_nea_LIT.txt",
    params:
        dir_stat_out = "results/maladapt/dir_stat_out_{n}/",
        dir_out = "results/maladapt/dir_out_{n}/",
        vcf_modern_out = "results/maladapt/vcf_modern_out_{n}/",
        vcf_archaic_out = "results/maladapt/vcf_archaic_out_{n}/",
        script = "resources/tools/maladapt/empirical/1empirical_compute-features.py",
    resources:
        mem_gb = lambda wildcards: 32 if (int(wildcards.n) < 175) or (int(wildcards.n) > 195) else 1000, 
        partition = lambda wildcards: "short" if (int(wildcards.n) < 175) or (int(wildcards.n) > 195) else "himem",
        time = lambda wildcards: 360 if (int(wildcards.n) < 175) or (int(wildcards.n) > 195) else 1440,
    conda: "../envs/maladapt-env.yaml",
    shell:
        """
        python {params.script} --vcfm {params.vcf_modern_out} --vcfa {params.vcf_archaic_out} --out {params.dir_out} --stat {params.dir_stat_out} --range {input.overlap}
        """


rule merge_maladapt_features:
    input:
        expand("results/maladapt/dir_stat_out_{n}/6window_nea_LIT.txt", n=range(1,1001)),
    output:
        features = "results/maladapt/6window_nea_LIT.txt",
    shell:
        """
        head -1 results/maladapt/dir_stat_out_1/6window_nea_LIT.txt > {output.features}
        cat {input} | grep -v chr >> {output.features}
        rm -r results/maladapt/dir_stat_out_*
        rm -r results/maladapt/dir_out_*
        rm -r results/maladapt/vcf_modern_out_*
        rm -r results/maladapt/vcf_archaic_out_*
        """


rule run_maladapt_prediction:
    input:
        features = rules.merge_maladapt_features.output.features,
        flag = rules.download_maladapt_pretrained_model.output.download_flag,
    output:
        predictions = "results/maladapt/6window_nea_LIT.predictions",
    params:
        model = "resources/tools/maladapt/pretrained_models/MaLAdapt_25_-sweep-all_model-001.sav",
        feature_config = "resources/tools/maladapt/feature/set6-all.txt",
        script = "resources/tools/maladapt/empirical/3empirical_prediction.py",
    resources:
        mem_gb = 128,
    conda: "../envs/maladapt-env.yaml",
    shell:
        """
        python {params.script} --input {input.features} --output {output.predictions} --model {params.model} --feature {params.feature_config}
        """
