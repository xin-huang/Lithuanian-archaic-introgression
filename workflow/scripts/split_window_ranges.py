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

import os
import warnings


def split_rows_into_files(input_filename, output_dir, output_prefix, n):
    if n <= 0:
        raise ValueError(f"Invalid value for number of files: {n}. Please provide a positive integer greater than 0.")

    os.makedirs(output_dir, exist_ok=True)

    with open(input_filename, 'r') as infile:
        input_data = [list(map(int, line.strip().split())) for line in infile if line.strip()]

    elements_per_row = len(input_data[0])

    # If n is greater than the number of elements in a row, adjust n and raise a warning
    if n > elements_per_row:
        warnings.warn(f"Requested {n} files, but each row only has {elements_per_row} elements. Creating {elements_per_row} files instead.")
        n = elements_per_row

    elements_per_file = elements_per_row // n
    remainder = elements_per_row % n

    for file_counter in range(1, n + 1):
        output_file_path = os.path.join(output_dir, f'{output_prefix}_{file_counter}.txt')
        with open(output_file_path, 'w') as outfile:
            for row in input_data:
                start_index = (file_counter - 1) * elements_per_file + min(file_counter - 1, remainder)
                end_index = start_index + elements_per_file + (1 if file_counter <= remainder else 0)

                if start_index < len(row):  # Only write if there are elements in the slice
                    outfile.write(' '.join(map(str, row[start_index:end_index])) + '\n')


split_rows_into_files(
    input_filename=snakemake.input.stepsize_file, 
    output_dir=snakemake.params.output_dir,
    output_prefix=snakemake.params.output_prefix,
    n=snakemake.params.nfiles,
)
