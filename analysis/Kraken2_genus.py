# coding:utf-8
"""
Author  : Tian
Time    : 2023-06-18 15:25
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




def collect_taxids(folder_path):
    taxids = set()

    for filename in os.listdir(folder_path):
        with open(os.path.join(folder_path, filename), 'r') as f:
            for line in f:
                line_parts = line.strip().split('\t')
                classification_status = line_parts[0]

                if classification_status != 'U':
                    taxids.add(line_parts[2])

    return list(taxids)


def analyze_file(filename, taxid_to_genus, results):
    with open(filename, 'r') as f:
        for line in f:
            line_parts = line.strip().split('\t')
            classification_status = line_parts[0]

            if classification_status == 'U':
                results['total_unclassified'] += 1
            else:
                results['total_classified'] += 1
                classification_id = taxid_to_genus.get(line_parts[2])
                if classification_id:
                    results['correct_classifications'] += 1
                else:
                    results['incorrect_classifications'] += 1

    return results


folder_path = '/home/zqtianqinzhong/software/ART/datasets/krakenuniq_results'

file_results_list = []

global_counter = {
    'total_classified': 0,
    'total_unclassified': 0,
    'correct_classifications': 0,
    'incorrect_classifications': 0
}

taxids = collect_taxids(folder_path)
taxid_to_genus = get_genus_taxids(taxids)

for filename in os.listdir(folder_path):
    species_name = re.match(r'(.+?)_HS', filename).group(1)

    species_id = get_taxid(species_name)
    genus_id = taxid_to_genus.get(species_id)

    file_results = {
        'filename': filename,
        'total_classified': 0,
        'total_unclassified': 0,
        'correct_classifications': 0,
        'incorrect_classifications': 0
    }

    file_results = analyze_file(os.path.join(folder_path, filename), taxid_to_genus, file_results)

    print(file_results)

    if file_results['correct_classifications'] > 0:
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

        file_results_list.append(file_results)

        for key in global_counter.keys():
            global_counter[key] += file_results[key]

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

with open('kraken2_results_genus.csv', 'w', newline='') as f:
    fieldnames = ['filename', 'total_classified', 'total_unclassified', 'correct_classifications',
                  'incorrect_classifications', 'precision', 'recall', 'f1_score', 'accuracy']
    writer = csv.DictWriter(f, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(file_results_list)
