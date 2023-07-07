# coding:utf-8
"""
Author  : Tian
Time    : 2023-06-22 10:46
Desc:
"""
import time

from Bio import Entrez
import os
import re
import csv

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
    handle = Entrez.esearch(db="taxonomy", term=species_name)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"][0]


def analyze_file(filename, genus_id, results):
    base_name = os.path.basename(filename).replace('_pathseq.out', '1.fasta')
    fasta_file = os.path.join('/home/zqtianqinzhong/software/ART/datasets/simulated_data_new', base_name)

    with open(filename, 'r') as f:
        for line in f:
            line_parts = line.strip().split('\t')
            if line_parts[2] == 'genus':
                classification_id = line_parts[0]
                if classification_id == genus_id:
                    results['correct_classifications'] += int(line_parts[7])
                    break
            elif line_parts[2] == 'root':
                results['total_classified'] = int(line_parts[7])

        results['incorrect_classifications'] = results['total_classified'] - results['correct_classifications']

    with open(fasta_file, 'r') as f:
        fasta_lines = sum(1 for _ in f)
    results['total_unclassified'] = fasta_lines - results['total_classified']

    return results


folder_path = '/home/zqtianqinzhong/software/ART/datasets/pathseq_results'

file_results_list = []

global_counter = {
    'total_classified': 0,
    'total_unclassified': 0,
    'correct_classifications': 0,
    'incorrect_classifications': 0
}


def get_genus_taxids(species_ids):
    handle = safe_entrez_request(Entrez.efetch, db="taxonomy", id=species_ids)
    records = Entrez.read(handle)
    handle.close()

    taxid_to_genus = {}
    for record in records:
        if record["Rank"] == "genus":
            taxid_to_genus[record["TaxId"]] = record["TaxId"]
        else:
            genus_taxid = None
            if "LineageEx" in record:
                for lineage in record["LineageEx"]:
                    if lineage["Rank"] == "genus":
                        genus_taxid = lineage["TaxId"]
                        break
            taxid_to_genus[record["TaxId"]] = genus_taxid

    return taxid_to_genus


for filename in os.listdir(folder_path):
    if filename.endswith('.out'):  # Add this line to filter for .out_PerRead files
        species_name = re.match(r'(.+?)_HS', filename).group(1)

        species_id = get_taxid(species_name)
        genus_id = get_genus_taxids(species_id)[species_id]

        file_results = {
            'filename': filename,
            'total_classified': 0,
            'total_unclassified': 0,
            'correct_classifications': 0,
            'incorrect_classifications': 0
        }

        file_results = analyze_file(os.path.join(folder_path, filename), genus_id, file_results)

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

with open('pathseq_results_genus.csv', 'w', newline='') as f:
    fieldnames = ['filename', 'total_classified', 'total_unclassified', 'correct_classifications',
                  'incorrect_classifications', 'precision', 'recall', 'f1_score', 'accuracy']
    writer = csv.DictWriter(f, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(file_results_list)
