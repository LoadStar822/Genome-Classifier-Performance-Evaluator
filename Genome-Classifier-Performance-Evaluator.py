# coding:utf-8
"""
Author  : Tian
Time    : 2023-04-21 12:01
Desc:
"""
import os
import subprocess
import time
import psutil
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
from colorama import Fore, Style
import logging
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Run experiments with specified classifiers and thread number.")
    parser.add_argument('-c', '--classifiers', type=str, required=True, nargs='+',
                        help="Space separated classifiers to run. e.g. clark-s krakenuniq")
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help="Number of threads to use.")
    return parser.parse_args()

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)
divider = '-' * 80


# Define a helper function to print in color
def print_color(text, color):
    logger.info(color + text + Style.RESET_ALL)

def get_memory_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss / (1024 ** 2)  # Returns the memory usage in MB


def run_experiment(command, conda_env=None, working_directory=None):
    if conda_env:
        command = f'conda run -n {conda_env} {command}'

    # Use /usr/bin/time to get time and memory usage
    command = f'/usr/bin/time -v {command}'

    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               cwd=working_directory)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"Command execution failed: {stderr.decode('utf-8')}")
        return None

    # Parse /usr/bin/time output
    time_output = stderr.decode('utf-8')

    # Using regex to extract time and memory
    time_regex = r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ((\d+:)?\d+:\d+\.\d+)'
    mem_regex = r'Maximum resident set size \(kbytes\): (\d+)'

    time_match = re.search(time_regex, time_output)
    mem_match = re.search(mem_regex, time_output)

    if time_match and mem_match:
        elapsed_time_str = time_match.group(1)
        memory_usage_kbytes = int(mem_match.group(1))

        # Convert elapsed time to seconds
        elapsed_time_parts = elapsed_time_str.split(':')
        if len(elapsed_time_parts) == 3:
            hours, minutes, seconds = elapsed_time_parts
            elapsed_time = float(hours) * 3600 + float(minutes) * 60 + float(seconds)
        else:
            minutes, seconds = elapsed_time_parts
            elapsed_time = float(minutes) * 60 + float(seconds)

        # Convert memory usage to megabytes
        memory_usage = memory_usage_kbytes / 1024

        return elapsed_time, memory_usage




conda_envs = {
    'kraken2': 'base',
    'clark': 'base',
    'megablast': 'base',
    'krakenuniq': 'base',
    'k-SLAM': 'kslam',
    'taxmaps': 'taxmaps',
    'pathseq': 'bioinformatics',
    'clark-s': 'base',
    'centrifuge': 'base',
    'diamond': 'bioinformatics',
}


def run_experiments(classifiers, input_folder, output_folders, databases, num_threads=1, num_repeats=1, executor=None):
    results = {classifier: {'total_time': 0, 'total_memory': 0} for classifier in classifiers}
    commands = []

    for file_prefix in os.listdir(input_folder):
        if file_prefix.endswith('1.fq'):
            file_base = file_prefix[:-4]
            forward_reads = os.path.join(input_folder, f'{file_base}1.fq')
            reverse_reads = os.path.join(input_folder, f'{file_base}2.fq')

            for classifier in classifiers:
                output_folder = output_folders[classifier]
                output_file = os.path.join(output_folder, f'{file_base}_{classifier}.out')

                if not os.path.exists(output_folder):
                    os.makedirs(output_folder)

                has_existing_output = any(
                    name.startswith(f"{file_base}_{classifier}") for name in os.listdir(output_folder)
                )
                if has_existing_output:
                    print(f"Skipping {file_base} for {classifier} because output already exists")
                    continue

                working_directory = os.path.expanduser('~')
                if classifier == 'kraken2':
                    command = f'~/software/kraken2/kraken2 --db {databases[classifier]} --paired --output {output_file} {forward_reads} {reverse_reads} --threads {num_threads}'
                elif classifier == 'clark':
                    command = f'~/software/CLARKSCV1.2.6.1/classify_metagenome.sh -n {num_threads} -P {forward_reads} {reverse_reads} -R {output_file}'
                    working_directory = os.path.expanduser('~/software/CLARKSCV1.2.6.1')
                elif classifier == 'clark-s':
                    command = f'~/software/CLARKSCV1.2.6.1/classify_metagenome.sh -n {num_threads} -P {forward_reads} {reverse_reads} -R {output_file} --spaced'
                    working_directory = os.path.expanduser('~/software/CLARKSCV1.2.6.1')
                elif classifier == 'krakenuniq':
                    command = f'krakenuniq --db {databases[classifier]} --paired --output {output_file} {forward_reads} {reverse_reads} --threads {num_threads}'
                elif classifier == 'megablast':
                    forward_reads_fasta = os.path.join(input_folder, f'{file_base}1.fasta')
                    command = f'blastn -task megablast -query {forward_reads_fasta} -db {databases[classifier]} -out {output_file} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames stitle"'
                elif classifier == 'k-SLAM':
                    command = f'~/software/k-SLAM/SLAM --db={databases[classifier]} --output-file={output_file} {forward_reads} {reverse_reads}'
                elif classifier == 'taxmaps':
                    command = f'~/software/taxmaps/taxMaps -1 {forward_reads} -2 {reverse_reads}  -d {databases[classifier] + "/refseq_microbial.lcak300.gem"} -t {databases[classifier] + "/taxonomy.tbl"} -p {file_base} -o {output_folder}/{file_base}_{classifier}'
                elif classifier == 'pathseq':
                    bam_file = os.path.join(bam_input_folder, f'{file_base}.bam')
                    command = f'gatk PathSeqPipelineSpark --input {bam_file} --microbe-bwa-image {databases[classifier]["microbe_bwa_image"]} --microbe-dict {databases[classifier]["microbe_dict"]} --taxonomy-file {databases[classifier]["taxonomy_file"]} --scores-output {output_file} --output {output_file.replace(".out", ".bam")}'
                elif classifier == 'centrifuge':
                    command = f'centrifuge -x {databases[classifier]} -1 {forward_reads} -2 {reverse_reads} -S {output_file} --threads {num_threads}'
                elif classifier == 'diamond':
                    command = f'diamond blastp -d {databases[classifier]} -q {protein_folder} -o {output_file} --threads {num_threads}'

                for _ in range(num_repeats):
                    conda_env = conda_envs.get(classifier)
                    command_info = {
                        'command': command,
                        'conda_env': conda_env,
                        'working_directory': working_directory,
                        'classifier': classifier
                    }
                    commands.append(command_info)

    futures = []
    future_to_classifier = {}
    total_commands = len(commands)
    completed_commands = 0
    for command_info in commands:
        future = executor.submit(run_experiment, command_info['command'],
                                 conda_env=command_info['conda_env'],
                                 working_directory=command_info['working_directory'])
        future_to_classifier[future] = command_info['classifier']
        futures.append(future)

    for future in as_completed(futures):
        classifier = future_to_classifier[future]  # Get the classifier for this future
        try:
            result = future.result()
            if result is not None:
                elapsed_time, memory_usage = result
                results[classifier]['total_time'] += elapsed_time
                results[classifier]['total_memory'] += memory_usage
                print_color(
                    f"Task for {classifier} finished. Elapsed time: {elapsed_time:.2f}s, Memory usage: {memory_usage:.2f}MB",
                    Fore.BLUE)
        except Exception as e:
            print_color(f"Error while running experiment for {classifier}: {e}", Fore.RED)
        completed_commands += 1
        remaining_commands = total_commands - completed_commands
        print(f"Completed {completed_commands} of {total_commands} commands. Remaining: {remaining_commands} commands.")

    return results


def save_results_to_file(results, output_file):
    with open(output_file, 'w') as f:
        f.write("Classifier\tTotal running time (s)\tTotal Memory Usage (MB)\n")
        for classifier in classifiers:
            total_time = results[classifier]['total_time']
            total_memory = results[classifier]['total_memory']
            f.write(f"{classifier}\t{total_time:.2f}\t{total_memory:.2f}\n")
    print(f"Results saved to {output_file}")


input_folder = os.path.expanduser('~/software/ART/datasets/simulated_data_new')
bam_input_folder = os.path.expanduser('~/software/ART/datasets/bam_files')
protein_folder = os.path.expanduser('~/software/ART/datasets/simulated_protein')

classifiers = ['clark-s', 'clark', 'krakenuniq', 'pathseq', 'kraken2', 'taxmaps', 'k-SLAM', 'megablast', 'centrifuge', 'diamond']
output_folders = {
    'kraken2': os.path.expanduser('~/software/ART/datasets/kraken2_results'),
    'clark': os.path.expanduser('~/software/ART/datasets/clark_results'),
    'megablast': os.path.expanduser('~/software/ART/datasets/megablast_results'),
    'krakenuniq': os.path.expanduser('~/software/ART/datasets/krakenuniq_results'),
    'k-SLAM': os.path.expanduser('~/software/ART/datasets/k-SLAM_results'),
    'taxmaps': os.path.expanduser('~/software/ART/datasets/taxmaps_results'),
    'pathseq': os.path.expanduser('~/software/ART/datasets/pathseq_results'),
    'clark-s': os.path.expanduser('~/software/ART/datasets/clark-s_results'),
    'centrifuge': os.path.expanduser('~/software/ART/datasets/centrifuge_results'),
    'diamond': os.path.expanduser('~/software/ART/datasets/diamond_results')
}
databases = {
    'kraken2': os.path.expanduser('~/software/kraken2/standard'),
    'clark': os.path.expanduser('~/software/CLARKSCV1.2.6.1/DIR_DB'),
    'clark-s': os.path.expanduser('~/software/CLARKSCV1.2.6.1/DIR_DB'),
    'megablast': os.path.expanduser('~/software/MegaBLAST/refseq_rna'),
    'krakenuniq': os.path.expanduser('~/software/krakenuniq'),
    'k-SLAM': os.path.expanduser('~/software/k-SLAM/refseq'),
    'taxmaps': os.path.expanduser('~/software/taxmaps'),
    'pathseq': {
        'microbe_bwa_image': os.path.expanduser('~/software/pathseq/pathseq_microbe.fa.img'),
        'microbe_dict': os.path.expanduser('~/software/pathseq/pathseq_microbe.dict'),
        'taxonomy_file': os.path.expanduser('~/software/pathseq/pathseq_taxonomy.db')
    },
    'centrifuge': os.path.expanduser('~/software/centrifuge-master/refseq/p+h+v'),
    'diamond': os.path.expanduser('~/software/diamond/nr.dmnd')
}
if __name__ == '__main__':
    args = parse_args()
    classifiers = args.classifiers
    num_threads = args.threads

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        future_to_classifier = {}
        results = run_experiments(classifiers, input_folder, output_folders, databases, executor=executor)
    current_time = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
    result_output_file = f"results_{current_time}.txt"
    save_results_to_file(results, result_output_file)

    for classifier in classifiers:
        total_time = results[classifier]['total_time']
        total_memory = results[classifier]['total_memory']
        print(f'{classifier}Total running time: {total_time:.2f}s，Total Memory Usage: {total_memory:.2f}MB')
