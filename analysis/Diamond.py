# coding:utf-8
"""
Author  : Tian
Time    : 2023-06-23 15:14
Desc:
"""
from Bio import Entrez
import os
import re
import csv
from collections import defaultdict

Entrez.email = ""


def get_taxid(species_name):
    handle = Entrez.esearch(db="taxonomy", term=species_name)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"][0]


def analyze_file(filename, species_id, results):
    base_name = os.path.basename(filename).replace('_diamond.out', '1.fasta')
    fasta_file = os.path.join('/home/zqtianqinzhong/software/ART/datasets/simulated_data_new', base_name)

    groups = defaultdict(list)
    processed_keys = set()
    with open(filename, 'r') as f:
        for line in f:
            line_parts = line.strip().split('\t')
            if line_parts[0] not in processed_keys:  # Skip if key already processed
                groups[line_parts[0]].append(line_parts)
                processed_keys.add(line_parts[0])

    for key, lines in groups.items():
        line_parts = lines[0]  # Use the first line instead of maximum
        if len(line_parts) < 2:
            results['total_classified'] += 1
            results['incorrect_classifications'] += 1
            continue
        results['total_classified'] += 1
        classification_id = line_parts[1].split(';')[0]
        if classification_id == species_id:
            results['correct_classifications'] += 1
        else:
            results['incorrect_classifications'] += 1

    with open(fasta_file, 'r') as f:
        fasta_lines = sum(1 for _ in f)
    results['total_unclassified'] += fasta_lines/2 - results['total_classified']

    return results



folder_path = '/home/zqtianqinzhong/software/ART/datasets/diamond_results'

file_results_list = []

global_counter = {
    'total_classified': 0,
    'total_unclassified': 0,
    'correct_classifications': 0,
    'incorrect_classifications': 0
}

for filename in os.listdir(folder_path):
    species_name = re.match(r'(.+?)_HS', filename).group(1)

    species_id = get_taxid(species_name)

    file_results = {
        'filename': filename,
        'total_classified': 0,
        'total_unclassified': 0,
        'correct_classifications': 0,
        'incorrect_classifications': 0
    }

    file_results = analyze_file(os.path.join(folder_path, filename), species_id, file_results)

    TP = file_results['correct_classifications']
    FP = file_results['incorrect_classifications']
    FN = file_results['total_unclassified']

    precision = TP / (TP + FP) if TP + FP > 0 else 0
    recall = TP / (TP + FN) if TP + FN > 0 else 0
    f1_score = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0
    accuracy = TP / (TP + FP + FN) if TP + FP + FN > 0 else 0

    file_results.update({
        'precision': precision,
        'recall': recall,
        'f1_score': f1_score,
        'accuracy': accuracy
    })

    for key in global_counter:
        global_counter[key] += file_results[key]

    file_results_list.append(file_results)

TP = global_counter['correct_classifications']
FP = global_counter['incorrect_classifications']
FN = global_counter['total_unclassified']

precision = TP / (TP + FP) if TP + FP > 0 else 0
recall = TP / (TP + FN) if TP + FN > 0 else 0
f1_score = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0
accuracy = TP / (TP + FP + FN) if TP + FP + FN > 0 else 0

summary_row = {
    'filename': 'Overall',
    'total_classified': global_counter['total_classified'],
    'total_unclassified': global_counter['total_unclassified'],
    'correct_classifications': TP,
    'incorrect_classifications': FP,
    'precision': precision,
    'recall': recall,
    'f1_score': f1_score,
    'accuracy': accuracy
}

file_results_list.append(summary_row)

with open('diamond_results.csv', 'w', newline='') as f:
    fieldnames = ['filename', 'total_classified', 'total_unclassified', 'correct_classifications',
                  'incorrect_classifications', 'precision', 'recall', 'f1_score', 'accuracy']
    writer = csv.DictWriter(f, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(file_results_list)
