# coding:utf-8
"""
Author  : Tian
Time    : 2023-07-05 13:20
Desc:
"""
import time
from collections import defaultdict

from Bio import Entrez
import os
import re
import csv

Entrez.email = ""

centrifuge_folder_path = '/home/zqtianqinzhong/software/ART/datasets/centrifuge_results'
kraken2_folder_path = '/home/zqtianqinzhong/software/ART/datasets/kraken2_results'
krakenuniq_folder_path = '/home/zqtianqinzhong/software/ART/datasets/krakenuniq_results'
clark_folder_path = '/home/zqtianqinzhong/software/ART/datasets/clark_results'
clarks_folder_path = '/home/zqtianqinzhong/software/ART/datasets/clark-s_results'
kslam_folder_path = '/home/zqtianqinzhong/software/ART/datasets/k-SLAM_results'
megablast_folder_path = '/home/zqtianqinzhong/software/ART/datasets/megablast_results'
kaiju_folder_path = '/home/zqtianqinzhong/software/ART/datasets/kaiju_results'
group_results = {}
global_genus = {}
taxid_to_genus = {}

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
    if isinstance(species_ids, list) and len(species_ids) == 1 and species_ids[0] == '0':
        return {species_ids[0]: '0'}
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
                    if "Rank" in lineage and lineage["Rank"] == "genus":
                        genus_taxid = lineage["TaxId"]
                        break

            taxid_to_genus[record["TaxId"]] = genus_taxid

    return taxid_to_genus



def analyze_file_centrifuge(filename):
    groups = {}
    with open(filename, 'r') as f:
        next(f)
        for line in f:
            line_parts = line.strip().split('\t')
            group_key = line_parts[0]

            if group_key not in groups:
                groups[group_key] = line_parts
            else:
                continue

            classification_status = line_parts[1]
            group_results[group_key] = {}
            if classification_status == 'unclassified':
                if 'unclassified' not in group_results[group_key]:
                    group_results[group_key]['unclassified'] = 0
                group_results[group_key]['unclassified'] += 1
            else:
                classification_id = line_parts[2]
                classification_id = taxid_to_genus.get(classification_id)
                if classification_id not in group_results[group_key]:
                    group_results[group_key][classification_id] = 0
                group_results[group_key][classification_id] += 2


def analyze_file_kraken2(filename):
    with open(filename, 'r') as f:
        for line in f:
            line_parts = line.strip().split('\t')
            classification_status = line_parts[0]
            group_key = line_parts[1]
            if classification_status == 'U':
                if 'unclassified' not in group_results[group_key]:
                    group_results[group_key]['unclassified'] = 0
                group_results[group_key]['unclassified'] += 1
            else:
                classification_id = line_parts[2]
                classification_id = taxid_to_genus.get(classification_id)
                if classification_id not in group_results[group_key]:
                    group_results[group_key][classification_id] = 0
                group_results[group_key][classification_id] += 2


def analyze_file_krakenuniq(filename):
    with open(filename, 'r') as f:
        for line in f:
            line_parts = line.strip().split('\t')
            classification_status = line_parts[0]
            group_key = line_parts[1]
            if classification_status == 'U':
                if 'unclassified' not in group_results[group_key]:
                    group_results[group_key]['unclassified'] = 0
                group_results[group_key]['unclassified'] += 1
            else:
                classification_id = line_parts[2]
                classification_id = taxid_to_genus.get(classification_id)
                if classification_id not in group_results[group_key]:
                    group_results[group_key][classification_id] = 0
                group_results[group_key][classification_id] += 2


def analyze_file_clark(filename):
    with open(filename, 'r') as f:
        next(f)
        for line in f:
            line_parts = line.strip().split(',')
            classification_status = line_parts[2]
            group_key = line_parts[0]
            if classification_status == 'NA':
                if 'unclassified' not in group_results[group_key]:
                    group_results[group_key]['unclassified'] = 0
                group_results[group_key]['unclassified'] += 1
            else:
                classification_id = line_parts[2]
                classification_id = taxid_to_genus.get(classification_id)
                if classification_id not in group_results[group_key]:
                    group_results[group_key][classification_id] = 0
                group_results[group_key][classification_id] += 60


def analyze_file_clarks(filename):
    with open(filename, 'r') as f:
        next(f)
        for line in f:
            line_parts = line.strip().split(',')
            classification_status = line_parts[3]
            group_key = line_parts[0]
            if classification_status == 'NA':
                if 'unclassified' not in group_results[group_key]:
                    group_results[group_key]['unclassified'] = 0
                group_results[group_key]['unclassified'] += 1
            else:
                classification_id = line_parts[3]
                classification_id = taxid_to_genus.get(classification_id)
                if classification_id not in group_results[group_key]:
                    group_results[group_key][classification_id] = 0
                group_results[group_key][classification_id] += 20


def analyze_file_kslam(filename):
    with open(filename, 'r') as f:
        for line in f:
            line_parts = line.strip().split('\t')
            classification_status = line_parts[1]
            group_key = line_parts[0]
            if classification_status == '0':
                if 'unclassified' not in group_results[group_key]:
                    group_results[group_key]['unclassified'] = 0
                group_results[group_key]['unclassified'] += 1
            else:
                classification_id = line_parts[1]
                classification_id = taxid_to_genus.get(classification_id)
                if classification_id not in group_results[group_key]:
                    group_results[group_key][classification_id] = 0
                group_results[group_key][classification_id] += 1


def analyze_file_megablast(filename):
    groups = defaultdict(list)
    with open(filename, 'r') as f:
        for line in f:
            line_parts = line.strip().split('\t')
            groups[line_parts[0]].append(line_parts)

    max_lines = [max(lines, key=lambda x: float(x[11])) for lines in groups.values()]
    group_key = line_parts[0].split('/')[0]
    for line_parts in max_lines:
        classification_id = line_parts[12].split(';')[0]
        classification_id = taxid_to_genus.get(classification_id)
        if classification_id not in group_results[group_key]:
            group_results[group_key][classification_id] = 0
        group_results[group_key][classification_id] += 1


def analyze_file_kaiju(filename):
    with open(filename, 'r') as f:
        for line in f:
            line_parts = line.strip().split('\t')
            classification_status = line_parts[0]
            group_key = line_parts[1]
            if classification_status == 'U':
                if 'unclassified' not in group_results[group_key]:
                    group_results[group_key]['unclassified'] = 0
                group_results[group_key]['unclassified'] += 1
            else:
                classification_id = line_parts[2]
                classification_id = taxid_to_genus.get(classification_id)
                if classification_id not in group_results[group_key]:
                    group_results[group_key][classification_id] = 0
                group_results[group_key][classification_id] += 1


def analyze_file(filename, genus_id, results):
    centrifuge_name = centrifuge_folder_path + '/' + filename + '_centrifuge.out'
    kraken2_name = kraken2_folder_path + '/' + filename + '_kraken2.out'
    krakenuniq_name = krakenuniq_folder_path + '/' + filename + '_krakenuniq.out'
    clark_name = clark_folder_path + '/' + filename + '_clark.out.csv'
    clarks_name = clarks_folder_path + '/' + filename + '_clark-s.out.csv'
    kslam_name = kslam_folder_path + '/' + filename + '_k-SLAM.out_PerRead'
    megablast_name = megablast_folder_path + '/' + filename + '_megablast.out'
    kaiju_name = kaiju_folder_path + '/' + filename + '_kaiju.out'
    # 合并为一个列表
    analyze_file_centrifuge(centrifuge_name)
    analyze_file_kraken2(kraken2_name)
    analyze_file_krakenuniq(krakenuniq_name)
    analyze_file_clark(clark_name)
    analyze_file_clarks(clarks_name)
    analyze_file_kslam(kslam_name)
    analyze_file_megablast(megablast_name)
    analyze_file_kaiju(kaiju_name)

    for group_key in group_results:
        # If the group only contains 'unclassified', count it as unclassified
        if set(group_results[group_key].keys()) == {'unclassified'}:
            results['total_unclassified'] += 1
        else:
            results['total_classified'] += 1
            # Exclude 'unclassified' and find the classification_status with the highest count
            max_classification_status = max((k for k in group_results[group_key] if k != 'unclassified'),
                                            key=group_results[group_key].get)
            if max_classification_status == genus_id:
                results['correct_classifications'] += 1
            else:
                results['incorrect_classifications'] += 1

    return results

def collect_taxids_centrifuge(folder_path):
    taxids = set()

    for filename in os.listdir(folder_path):
        with open(os.path.join(folder_path, filename), 'r') as f:
            for line in f:
                line_parts = line.strip().split('\t')
                classification_status = line_parts[0]

                if classification_status != 'U':
                    taxids.add(line_parts[2])

    return list(taxids)

def collect_taxids_clark(folder_path):
    taxids = set()

    for filename in os.listdir(folder_path):
        with open(os.path.join(folder_path, filename), 'r') as f:
            next(f)
            for line in f:
                line_parts = line.strip().split(',')
                classification_status = line_parts[2]

                if classification_status != 'NA':
                    taxids.add(line_parts[2])

    return list(taxids)

def collect_taxids_clarks(folder_path):
    taxids = set()

    for filename in os.listdir(folder_path):
        with open(os.path.join(folder_path, filename), 'r') as f:
            next(f)
            for line in f:
                line_parts = line.strip().split(',')
                classification_status = line_parts[3]

                if classification_status != 'NA':
                    taxids.add(line_parts[3])

    return list(taxids)

def collect_taxids_kslam(folder_path):
    taxids = set()

    for filename in os.listdir(folder_path):
        if filename.endswith('.out_PerRead'):  # Add this line to filter for .out_PerRead files

            with open(os.path.join(folder_path, filename), 'r') as f:
                for line in f:
                    line_parts = line.strip().split('\t')
                    classification_status = line_parts[1]

                    if classification_status != 'U':
                        taxids.add(line_parts[1])

    return list(taxids)

def collect_taxids_kaiju(folder_path):
    taxids = set()

    for filename in os.listdir(folder_path):
        with open(os.path.join(folder_path, filename), 'r') as f:
            for line in f:
                line_parts = line.strip().split('\t')
                classification_status = line_parts[0]

                if classification_status != 'U':
                    taxids.add(line_parts[2])

    return list(taxids)

def collect_taxids_kraken2(folder_path):
    taxids = set()

    for filename in os.listdir(folder_path):
        with open(os.path.join(folder_path, filename), 'r') as f:
            for line in f:
                line_parts = line.strip().split('\t')
                classification_status = line_parts[0]

                if classification_status != 'U':
                    taxids.add(line_parts[2])

    return list(taxids)

def collect_taxids_krakenuniq(folder_path):
    taxids = set()

    for filename in os.listdir(folder_path):
        with open(os.path.join(folder_path, filename), 'r') as f:
            for line in f:
                line_parts = line.strip().split('\t')
                classification_status = line_parts[0]

                if classification_status != 'U':
                    taxids.add(line_parts[2])

    return list(taxids)

def collect_taxids_megablast(folder_path):
    taxids = set()

    for filename in os.listdir(folder_path):
        with open(os.path.join(folder_path, filename), 'r') as f:
            for line in f:
                line_parts = line.strip().split('\t')
                taxids.add(line_parts[12].split(';')[0])

    return list(taxids)

file_results_list = []

global_counter = {
    'total_classified': 0,
    'total_unclassified': 0,
    'correct_classifications': 0,
    'incorrect_classifications': 0
}
taxids = collect_taxids_centrifuge(centrifuge_folder_path)
taxid_to_genus.update(get_genus_taxids(taxids))
print('centrifuge done' + str(len(taxid_to_genus)))
taxids = collect_taxids_clark(clark_folder_path)
taxid_to_genus.update(get_genus_taxids(taxids))
print('clark done' + str(len(taxid_to_genus)))
taxids = collect_taxids_clarks(clarks_folder_path)
taxid_to_genus.update(get_genus_taxids(taxids))
print('clarks done' + str(len(taxid_to_genus)))
taxids = collect_taxids_kslam(kslam_folder_path)
taxid_to_genus.update(get_genus_taxids(taxids))
print('kslam done' + str(len(taxid_to_genus)))
taxids = collect_taxids_kaiju(kaiju_folder_path)
taxid_to_genus.update(get_genus_taxids(taxids))
print('kaiju done' + str(len(taxid_to_genus)))
taxids = collect_taxids_kraken2(kraken2_folder_path)
taxid_to_genus.update(get_genus_taxids(taxids))
print('kraken2 done' + str(len(taxid_to_genus)))
taxids = collect_taxids_krakenuniq(krakenuniq_folder_path)
taxid_to_genus.update(get_genus_taxids(taxids))
print('krakenuniq done' + str(len(taxid_to_genus)))
taxids = collect_taxids_megablast(megablast_folder_path)
taxid_to_genus.update(get_genus_taxids(taxids))
print('megablast done' + str(len(taxid_to_genus)))

for filename in os.listdir(centrifuge_folder_path):
    print(filename + ' started')
    group_results = {}

    species_name = re.match(r'(.+?)_HS', filename).group(1)

    species_id = get_taxid(species_name)
    genus_id = taxid_to_genus.get(species_id)

    name = filename.rsplit("_centrifuge.out", 1)[0]

    file_results = {
        'basename': filename,
        'total_classified': 0,
        'total_unclassified': 0,
        'correct_classifications': 0,
        'incorrect_classifications': 0
    }

    file_results = analyze_file(name, genus_id, file_results)

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
    'basename': 'Overall',
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

with open('merge_results_genus.csv', 'w', newline='') as f:
    fieldnames = ['basename', 'total_classified', 'total_unclassified', 'correct_classifications',
                  'incorrect_classifications', 'precision', 'recall', 'f1_score', 'accuracy']
    writer = csv.DictWriter(f, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(file_results_list)
