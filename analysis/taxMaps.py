# coding:utf-8
"""
Author  : Tian
Time    : 2023-06-30 15:20
Desc:
"""
from Bio import Entrez
import os
import re
import csv
import time

Entrez.email = ""


def safe_entrez_request(func, *args, max_retries=10, wait_time=5, **kwargs):
    retries = 0
    while retries < max_retries:
        try:
            return func(*args, **kwargs)
        except Exception as e:
            print(f"Request failed with error {str(e)}, retrying ({retries + 1} of {max_retries})")
            time.sleep(wait_time)
            retries += 1
    raise Exception("Max retries exceeded")
def get_taxid(species_name):
    handle = safe_entrez_request(Entrez.esearch, db="taxonomy", term=species_name)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"][0]


def analyze_file(filename, species_id, results):
    with open(filename, 'r') as f:
        for line in f:
            line_parts = line.strip().split('\t')
            classification_status = line_parts[0]
            if classification_status == '0' or classification_status == '-':
                results['total_unclassified'] += int(line_parts[4])
            else:
                results['total_classified'] += int(line_parts[4])
                classification_id = line_parts[0]
                if classification_id == species_id:
                    results['correct_classifications'] += int(line_parts[4])
                else:
                    results['incorrect_classifications'] += int(line_parts[4])

    return results


folder_path = '/home/zqtianqinzhong/software/ART/datasets/taxmaps_results'

file_results_list = []

global_counter = {
    'total_classified': 0,
    'total_unclassified': 0,
    'correct_classifications': 0,
    'incorrect_classifications': 0
}

for root, dirs, files in os.walk(folder_path):
    for filename in files:
        if filename.endswith('.summary'):  # This ensures we're only dealing with the .summary files
            full_file_path = os.path.join(root, filename)
            species_name = re.match(r'(.+?)_HS', filename).group(1)

            species_id = get_taxid(species_name)

            file_results = {
                'filename': filename,
                'total_classified': 0,
                'total_unclassified': 0,
                'correct_classifications': 0,
                'incorrect_classifications': 0
            }

            file_results = analyze_file(full_file_path, species_id, file_results)

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

with open('taxmaps_results.csv', 'w', newline='') as f:
    fieldnames = ['filename', 'total_classified', 'total_unclassified', 'correct_classifications',
                  'incorrect_classifications', 'precision', 'recall', 'f1_score', 'accuracy']
    writer = csv.DictWriter(f, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(file_results_list)
