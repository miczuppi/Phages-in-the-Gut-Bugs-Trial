#!/bin/bash

#Set folders
base_folder="path_to_dir"
cluster_dir="${base_dir}/cluster"

#set bwamem file location
bwamem_file="${base_folder}/bwamem.sh" 
#bwamem.sh file can be found on the GitHub repository
# Reference catalog
ref_catalog="${cluster_dir}/gbt_phages_gene_catalog"

#Index the reference catalog
bwa index "$ref_catalog"
samtools faidx "$ref_catalog"


# Move to the directory containing QC reads and run
cd /mnt/projects/gutbugs/mgx/qc_reads_sorted
ls -1 *sortedR1* | xargs -P 5 -I {}  "$bwamem_file" {} 

