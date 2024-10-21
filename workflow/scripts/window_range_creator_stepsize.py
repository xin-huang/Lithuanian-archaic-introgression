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
import argparse


def get_max_position(vcf_file, chromosome):

    vcf = pysam.VariantFile(vcf_file)

    max_position = 0

    # Iterate
    for record in vcf.fetch(chromosome):
        if record.pos > max_position:
            max_position = record.pos

    return max_position



def scan_vcf(vcf_file, window_size=50000, stepsize=10000, start_snp=0):
    vcf = pysam.VariantFile(vcf_file)
    window_starts = []
    window_ends = []
    window_lengths = []
    snps_counts = []

    for chrom in vcf.header.contigs:
        chrom_length = vcf.header.contigs[chrom].length

        if not chrom_length:
            chrom_length = get_max_position(vcf_file, chrom)

        start = start_snp

        while start < chrom_length:
            end = min(start + window_size, chrom_length)
            snps_in_window = sum(1 for _ in vcf.fetch(chrom, start, end))

            length = end - start

            if snps_in_window > 0:

                window_starts.append(start + 1)
                window_ends.append(end)
                window_lengths.append(length)
                snps_counts.append(snps_in_window)

            start = start + stepsize

    return window_starts, window_ends, window_lengths, snps_counts


def write_output(window_starts, window_ends, window_lengths, snps_counts, output_file):
    with open(output_file, 'w') as f:
        f.write(' '.join(map(str, window_starts)) + '\n')
        f.write(' '.join(map(str, window_ends)) + '\n')
        f.write(' '.join(map(str, window_lengths)) + '\n')
        f.write(' '.join(map(str, snps_counts)) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf_file', type=str)
    parser.add_argument('output_file', type=str)
    parser.add_argument('window_size', type=int, default=50000)
    parser.add_argument('stepsize', nargs='?', type=int, default=10000)

    args = parser.parse_args()

    window_starts, window_ends, window_lengths, snps_counts = scan_vcf(args.vcf_file, args.window_size, args.stepsize)
    write_output(window_starts, window_ends, window_lengths, snps_counts, args.output_file)
    print(f"Output written to {args.output_file}")


if __name__ == "__main__":
    main()
