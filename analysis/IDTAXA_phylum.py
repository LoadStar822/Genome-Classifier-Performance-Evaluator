"""
@Name: IDTAXA_phylum.py
@Auth: Tian
@Date: 2023/7/14-18:45
@Desc: 
@Ver : 0.0.0
"""
import time

from Bio import Entrez
import os
import re
import csv

Entrez.email = ""
species_to_phylum = {
    'Lactobacillus_acidophilus': 'Firmicutes',
    'Escherichia_coli': 'Proteobacteria',
    'Bacillus_subtilis': 'Firmicutes',
    'Methanosarcina_mazei': 'Halobacterota',
    'Methanocaldococcus_jannaschii': 'Euryarchaeota',
    'Haemophilus_influenzae': 'Proteobacteria',
    'Thermococcus_kodakarensis': 'Euryarchaeota',
    'Neisseria_meningitidis': 'Proteobacteria',
    'Nanoarchaeum_equitans': 'Nanoarchaeota',
    'Pseudomonas_aeruginosa': 'Proteobacteria',
    'Aeropyrum_pernix': 'Crenarchaeota',
    'Helicobacter_pylori': 'Campylobacterota',
    'Listeria_monocytogenes': 'Firmicutes',
    'Halobacterium_salinarum': 'Halobacterota',
    'Methanobrevibacter_smithii': 'Thermoplasmatota',
    'Chlamydia_trachomatis': 'Verrucomicrobiota',
    'Staphylococcus_aureus': 'Firmicutes',
    'Pyrococcus_furiosus': 'Euryarchaeota',
    'Vibrio_cholerae': 'Proteobacteria',
    'Mycobacterium_tuberculosis': 'Actinobacteriota',
    'Sulfolobus_solfataricus': 'Crenarchaeota',
    'Legionella_pneumophila': 'Proteobacteria',
    'Enterococcus_faecalis': 'Firmicutes',
    'Methanococcus_maripaludis': 'Euryarchaeota',
    'Salmonella_enterica': 'Proteobacteria',
    'Streptococcus_pyogenes': 'Firmicutes',
    'Campylobacter_jejuni': 'Campylobacterota',
    'Clostridium_difficile': 'Firmicutes',
    'Corynebacterium_diphtheriae': 'Actinobacteriota',
    'Shigella_flexneri': 'Proteobacteria',
}

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


def analyze_file(filename, taxonomy_hierarchy, results):
    with open(filename, 'r') as f:
        next(f)
        for line in f:
            line = line.replace("\"", "")
            line_parts = line.strip().split(',')
            classification_status = line_parts[14]

            if classification_status == 'NA':
                results['total_unclassified'] += 1
            else:
                results['total_classified'] += 1
                classification_id = line_parts[6]
                if classification_id == taxonomy_hierarchy:
                    results['correct_classifications'] += 1
                else:
                    results['incorrect_classifications'] += 1

    return results

def get_taxonomy_hierarchy(species_name):
    if species_name in species_to_phylum:
        return species_to_phylum[species_name]

# def get_taxonomy_hierarchy(species_name):
#     handle = safe_entrez_request(Entrez.esearch, db="taxonomy", term=species_name)
#     record = Entrez.read(handle)
#     handle.close()
#     if record["IdList"]:
#         taxid = record["IdList"][0]
#         handle = safe_entrez_request(Entrez.efetch, db="taxonomy", id=taxid, retmode="xml")
#         records = Entrez.read(handle)
#         handle.close()
#         if records:
#             lineage = records[0]["Lineage"].split("; ")
#             return lineage
#     return []


folder_path = '/dev/disk5/zqtqz/project/zongshu/result/idtaxa_results'

file_results_list = []

global_counter = {
    'total_classified': 0,
    'total_unclassified': 0,
    'correct_classifications': 0,
    'incorrect_classifications': 0
}

for filename in os.listdir(folder_path):
    if filename.endswith('.csv'):
        species_name = re.match(r'(.+?)_HS', filename).group(1)
        taxonomy_hierarchy = get_taxonomy_hierarchy(species_name)

        file_results = {
            'filename': filename,
            'total_classified': 0,
            'total_unclassified': 0,
            'correct_classifications': 0,
            'incorrect_classifications': 0
        }

        file_results = analyze_file(os.path.join(folder_path, filename), taxonomy_hierarchy, file_results)

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

with open('idtaxa_results_phylum.csv', 'w', newline='') as f:
    fieldnames = ['filename', 'total_classified', 'total_unclassified', 'correct_classifications',
                  'incorrect_classifications', 'precision', 'recall', 'f1_score', 'accuracy']
    writer = csv.DictWriter(f, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(file_results_list)
