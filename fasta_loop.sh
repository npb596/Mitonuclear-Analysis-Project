#!/bin/bash

# Replace bed file as needed

filename=random_nuclear_170_exons.bed

# Load Anaconda module for the python script in loop

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/2-5.2.0

# Define arrays for each column of the bed file

chrom=(`awk '{print $1}' $filename`)
start_points=(`awk '{print $2}' $filename`)
end_points=(`awk '{print $3}' $filename`)
gene_name=(`awk '{print $4}' $filename`)

# Define a sequence of the number of exons/genes in the bed file

gene_number=`grep -c "$chrom\s" $filename`
gene_number_minus_one=`expr $gene_number - 1`
gene_sequence=`seq 0 $gene_number_minus_one`

# Create a directory to deposit all the fasta output generated from this script

mkdir fasta_output

# Obtain fasta files for all exons in the bed file
# The python script originally comes from https://github.com/Bahler-Lab/alignment-from-vcf and the file location should be changed accordingly.
# The version of the python script used here is modified by Dr. Nathan Hall at Auburn University to print fasta sequences with no variants
# Usage of the script can be found on the Bahler-Lab so input files should be changed accordingly 

for d in $gene_sequence; do
 python ../../modified-alignment-from-vcf.py rheMac8.fa ../../unzip_chr_macaque_files/combined_filtered_by_chr/Macaque_merged.${chrom[$d]} ${chrom[$d]} ${start_points[$d]} ${end_points[$d]} 2 ${gene_name[$d]} Macaque_names.txt
 mv ${gene_name[$d]} fasta_output/
done 
