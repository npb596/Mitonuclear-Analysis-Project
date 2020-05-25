#!/bin/bash

# Set first argument as input (must be a bed file with gene name at the end)

filename=$1

# Load Anaconda module for the python script in loop

#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load anaconda/2-5.2.0

# Define arrays for each column of the bed file

gene_name=(`awk '{print $4}' $filename`)

for x in ${gene_name[@]}; do

chrom=(`grep "\s${x}$" $filename | awk '{print $1}'`)
start_point=(`grep "\s${x}$" $filename | awk '{print $2}'`)
end_point=(`grep "\s${x}$" $filename | awk '{print $3}'`)

# Obtain fasta files for all exons in the bed file
# The python script originally comes from https://github.com/Bahler-Lab/alignment-from-vcf and the file location should be changed accordingly.
# The version of the python script used here is modified by Dr. Nathan Hall at Auburn University to print fasta sequences with no variants
# Usage of the script can be found on the Bahler-Lab so input files should be changed accordingly 

python ~/modified-alignment-from-vcf.py ~/rheMac8.fa ~/unzip_chr_macaque_files/combined_filtered_by_chr/Macaque_merged.${chrom}.gz ${chrom} ${start_point} ${end_point} 2 ${x}.fa

done 
