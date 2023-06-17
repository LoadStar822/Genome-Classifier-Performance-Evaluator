# coding:utf-8
"""
Author  : Tian
Time    : 2023-06-17 14:14
Desc:
"""
from Bio import Entrez
import os
import re

Entrez.email = "@gmail.com"

def get_taxid(species_name):
    handle = Entrez.esearch(db="taxonomy", term=species_name)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"][0]

def analyze_file(filename, species_id, results):
    with open(filename, 'r') as f:
        for line in f:
            line_parts = line.strip().split('\t')
            classification_status = line_parts[0]

            if classification_status == 'U':
                results['total_unclassified'] += 1
            else:
                results['total_classified'] += 1
                classification_id = line_parts[2]
                if classification_id == species_id:
                    results['correct_classifications'] += 1
                else:
                    results['incorrect_classifications'] += 1

    return results

folder_path = '/home/zqtianqinzhong/software/ART/datasets/krakenuniq_results'

global_results = {
    'total_classified': 0,
    'total_unclassified': 0,
    'correct_classifications': 0,
    'incorrect_classifications': 0
}

for filename in os.listdir(folder_path):
    species_name = re.match(r'(.+?)_HS', filename).group(1)

    species_id = get_taxid(species_name)

    global_results = analyze_file(os.path.join(folder_path, filename), species_id, global_results)

TP = global_results['correct_classifications']
FP = global_results['incorrect_classifications']
FN = global_results['total_unclassified']

precision = TP / (TP + FP)
recall = TP / (TP + FN)
f1_score = 2 * precision * recall / (precision + recall)
accuracy = TP / (TP + FP + FN)

print(f"Global Results:\nPrecision: {precision}\nRecall: {recall}\nF1-Score: {f1_score}\nAccuracy: {accuracy}\n")

with open('krakenuniq_results.txt', 'w') as f:
    f.write(f"TP: {TP}\nFP: {FP}\nFN: {FN}\n")
    f.write(f"Global Results:\nPrecision: {precision}\nRecall: {recall}\nF1-Score: {f1_score}\nAccuracy: {accuracy}\n")

