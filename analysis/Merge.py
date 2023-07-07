# coding:utf-8
"""
Author  : Tian
Time    : 2023-07-05 13:20
Desc:
"""
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


def get_taxid(species_name):
    handle = Entrez.esearch(db="taxonomy", term=species_name)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"][0]


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
                if classification_id not in group_results[group_key]:
                    group_results[group_key][classification_id] = 0
                group_results[group_key][classification_id] += 50


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
                if classification_id not in group_results[group_key]:
                    group_results[group_key][classification_id] = 0
                group_results[group_key][classification_id] += 5


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
                if classification_id not in group_results[group_key]:
                    group_results[group_key][classification_id] = 0
                group_results[group_key][classification_id] += 1


def analyze_file(filename, species_id, results):
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
            max_classification_status = max((k for k in group_results[group_key] if k != 'unclassified'), key=group_results[group_key].get)
            if max_classification_status == species_id:
                results['correct_classifications'] += 1
            else:
                results['incorrect_classifications'] += 1

    return results


file_results_list = []

global_counter = {
    'total_classified': 0,
    'total_unclassified': 0,
    'correct_classifications': 0,
    'incorrect_classifications': 0
}

for filename in os.listdir(centrifuge_folder_path):
    group_results = {}

    species_name = re.match(r'(.+?)_HS', filename).group(1)

    species_id = get_taxid(species_name)

    name = filename.rsplit("_centrifuge.out", 1)[0]

    file_results = {
        'basename': filename,
        'total_classified': 0,
        'total_unclassified': 0,
        'correct_classifications': 0,
        'incorrect_classifications': 0
    }

    file_results = analyze_file(name, species_id, file_results)

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

with open('merge_results.csv', 'w', newline='') as f:
    fieldnames = ['basename', 'total_classified', 'total_unclassified', 'correct_classifications',
                  'incorrect_classifications', 'precision', 'recall', 'f1_score', 'accuracy']
    writer = csv.DictWriter(f, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(file_results_list)
