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


def get_memory_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss / (1024 ** 2)  # Returns the memory usage in MB


def run_experiment(command):
    start_time = time.time()
    memory_before = get_memory_usage()
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Command execution failed: {e}")
        return None
    end_time = time.time()
    memory_after = get_memory_usage()
    elapsed_time = end_time - start_time
    memory_usage = memory_after - memory_before
    return elapsed_time, memory_usage


def run_experiments(classifiers, input_folder, output_folders, databases, num_threads=1, num_repeats=5):
    results = {classifier: {'total_time': 0, 'total_memory': 0} for classifier in classifiers}

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

                if classifier in ['kraken2']:
                    command = f'~/software/kraken2/kraken2 --db {databases[classifier]} --paired --output {output_file} {forward_reads} {reverse_reads} --threads {num_threads}'
                elif classifier == 'clark':
                    command = f'~/software/CLARKSCV1.2.6.1/classify_metagenome.sh -n {num_threads} -P {forward_reads} {reverse_reads} -R {output_file}'
                elif classifier == 'krakenuniq':
                    command = f'krakenuniq --db {databases[classifier]} --paired --output {output_file} {forward_reads} {reverse_reads} --threads {num_threads}'
                elif classifier == 'megablast':
                    forward_reads = os.path.join(input_folder, f'{file_base}1.fasta')
                    reverse_reads = os.path.join(input_folder, f'{file_base}2.fasta')
                    command = f'blastn -task megablast -query {forward_reads} -db {databases[classifier]} -out {output_file} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames stitle"'

                for _ in range(num_repeats):
                    result = run_experiment(command)
                    if result is not None:
                        elapsed_time, memory_usage = result
                        results[classifier]['total_time'] += elapsed_time
                        results[classifier]['total_memory'] += memory_usage

    return results


input_folder = os.path.expanduser('~/software/ART/datasets/simulated_data')

classifiers = ['k-SLAM' , 'megablast', 'clark', 'krakenuniq', 'kraken2']
output_folders = {
    'kraken2': os.path.expanduser('~/software/ART/datasets/kraken2_results'),
    'clark': os.path.expanduser('~/software/ART/datasets/clark_results'),
    'megablast': os.path.expanduser('~/software/ART/datasets/megablast_results'),
    'krakenuniq': os.path.expanduser('~/software/ART/datasets/krakenuniq_results'),
    'k-SLAM': os.path.expanduser('~/software/ART/datasets/k-SLAM_results')
}
databases = {
    'kraken2': os.path.expanduser('~/software/kraken2/standard'),
    'clark': os.path.expanduser('~/software/CLARKSCV1.2.6.1/DIR_DB'),
    'megablast': os.path.expanduser('~/software/MegaBLAST/refseq_rna'),
    'krakenuniq': os.path.expanduser('~/software/krakenuniq')
}
results = run_experiments(classifiers, input_folder, output_folders, databases)

for classifier in classifiers:
    total_time = results[classifier]['total_time']
    total_memory = results[classifier]['total_memory']
    print(f'{classifier}Total running time: {total_time:.2f}秒，Total Memory Usage: {total_memory:.2f}MB')
