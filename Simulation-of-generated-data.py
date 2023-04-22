# coding:utf-8
"""
Author  : Tian
Time    : 2023-04-22 13:40
Desc:
"""
import os
import subprocess

def generate_data(input_file_path, output_folder, ss, read_length, coverage, mean_fragment_length, sd):
    file_prefix = os.path.splitext(os.path.basename(input_file_path))[0]
    output_file_prefix = os.path.join(output_folder, f'{file_prefix}_{ss}_{read_length}_{coverage}_{mean_fragment_length}_{sd}')

    art_command = f'~/software/ART/art_bin_MountRainier/art_illumina -ss {ss} -i {input_file_path} -p -l {read_length} -f {coverage} -m {mean_fragment_length} -s {sd} -o {output_file_prefix}'
    subprocess.run(art_command, shell=True, check=True)

    print(f'Successfully generated simulated data for {file_prefix} with parameters {ss}, {read_length}, {coverage}, {mean_fragment_length}, {sd}.')

input_folder = os.path.expanduser('~/software/ART/datasets/fasta_sequences')
output_folder = os.path.expanduser('~/software/ART/datasets/simulated_data')

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

sequencing_systems = ['HS25', 'HS20']
read_lengths = { 'HS25': [100, 150], 'HS20': [100] }
coverages = [10, 20, 30]
mean_fragment_lengths = [200, 300, 400]
standard_deviations = [10, 25, 50]

for fasta_file in os.listdir(input_folder):
    if fasta_file.endswith('.fasta'):
        input_file_path = os.path.join(input_folder, fasta_file)

        for ss in sequencing_systems:
            for read_length in read_lengths[ss]:
                for coverage in coverages:
                    for mean_fragment_length in mean_fragment_lengths:
                        for sd in standard_deviations:
                            generate_data(input_file_path, output_folder, ss, read_length, coverage, mean_fragment_length, sd)