# coding:utf-8
"""
Author  : Tian
Time    : 2023-06-17 16:47
Desc:
"""
import os
import subprocess

directory = "/home/zqtianqinzhong/software/ART/datasets/simulated_protein"

for filename in os.listdir(directory):
    if filename.endswith(".fasta"):
        filepath = os.path.join(directory, filename)
        db_name = filepath.replace(".fasta", "DB")
        cmd = f"mmseqs createdb {filepath} {db_name}"
        subprocess.run(cmd, shell=True)
