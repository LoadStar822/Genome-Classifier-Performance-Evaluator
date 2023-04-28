#!/bin/bash

input_directory="/home/zqtianqinzhong/software/ART/datasets/simulated_data"
output_directory="/home/zqtianqinzhong/software/ART/datasets/bam_files"

mkdir -p "${output_directory}"

for read1_file in "${input_directory}"/*1.fq; do
    base_name=$(basename "${read1_file}" 1.fq)
    read2_file="${input_directory}/${base_name}2.fq"

    if [ -e "${read2_file}" ]; then
        echo "Converting ${base_name} to BAM format..."

        picard FastqToSam \
            F1="${read1_file}" \
            F2="${read2_file}" \
            O="${output_directory}/${base_name}.bam" \
            SM="${base_name}"
    else
        echo "Error: Paired read file for ${base_name} not found."
    fi
done