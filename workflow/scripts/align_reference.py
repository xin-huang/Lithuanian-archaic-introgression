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


import pysam
from copy import deepcopy


reference = pysam.FastaFile(snakemake.input.reference)
vcf_in = pysam.VariantFile(snakemake.input.vcf)
vcf_out = pysam.VariantFile(snakemake.output.vcf, 'w', header=vcf_in.header)

not_present_var = 0
for record in vcf_in:
    ref_allele = reference.fetch(record.chrom, record.pos-1, record.pos)

    if record.ref != ref_allele and ref_allele not in record.alts:
        print("WARNING! It seems that the reference allele from the fasta file is not one of the variants in the vcf-file!")
        not_present_var += 1

    if record.ref != ref_allele:
        print(f"Mismatch at {record.chrom}:{record.pos} - VCF: {record.ref}, REF: {ref_allele}")
        if len(record.alts) > 1:
            print("WARNING! There are multiallelic sites! The code should be modified accordingly!")
        old_ref = deepcopy(record.ref)
        record.ref = ref_allele

        for sample in record.samples:
            gt = record.samples[sample]['GT']
            gtnew = list(deepcopy(gt))
            for ie, elem in enumerate(gt):
                if elem == 0:
                    gtnew[ie] = 1
                if elem == 1:
                    gtnew[ie] = 0
            record.samples[sample]['GT'] = tuple(gtnew)

        record.alts = [old_ref]

    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
reference.close()
