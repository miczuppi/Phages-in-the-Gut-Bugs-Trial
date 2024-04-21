#!/bin/bash

# Input file
file=$1
# Extract sample ID from the input file name
sample_id=$(echo "$file" | cut -d_ -f1-2)

# Paths for various data directories
base_dir="path_to_dir"
rel_abund_dir="${base_dir}/rel_abund"
cluster_dir="${base_dir}/cluster"
python_scripts="${base_dir}/python_scripts"

# Paths for the sample's data files
input_file1="${sample_id}_sortedR1.fq.gz"
input_file2="${sample_id}_sortedR2.fq.gz"
output_bam_file="${rel_abund_dir}/${sample_id}.sort.bam"
filtered_bam_file="${rel_abund_dir}/${sample_id}.sort.filtered.ID95.bam"
flagstat_file="${rel_abund_dir}/${sample_id}.ID95.flagstat.txt"
abundance_file="${rel_abund_dir}/${sample_id}.abundance.txt"
count_file="${rel_abund_dir}/${sample_id}.count"
learn_file="${rel_abund_dir}/${sample_id}.learn"

# Paths for python scripts
#These require python v2.7.17
bamfilter="${python_scripts}/bamfilter_mate_saver.py"
bamcount="${python_scripts}/bam2counts.py"
#Python files can be found on the GitHub repository

# Reference catalog
ref_catalog="${cluster_dir}/gbt_phages_gene_catalog"
ref_catalog_fai="${cluster_dir}/gbt_phages_gene_catalog.fai"

# Number of threads for parallel processing
threads=2

# Memory per thread for sorting
memory_per_thread="2G"

# Align reads using BWA-MEM and filter the alignment with samtools
bwa mem -t "$threads" -M "$ref_catalog" "$input_file1" "$input_file2" | \
samtools view - -h -Su -F 2308 -q 0 |
samtools sort -n -@ "$threads" -m "$memory_per_thread" -O bam -o "$output_bam_file"
"$bamfilter" -i "$output_bam_file" -o "$filtered_bam_file" -f 0.95
samtools flagstat -@ "$threads" "$filtered_bam_file" > "$flagstat_file"
"$bamcount" -t "$threads"  -y 1 -i "$filtered_bam_file" -s "${rel_abund_dir}/$sample_id" -v "$ref_catalog_fai" -a "$abundance_file" -c "$count_file" -l "$learn_file"


