#!/bin/bash

# Load python module on ASC necessary for GC_script.py

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load python/2.7.9

query_genes=$1

# Define arrays of gene lengths, gene names, and chromosomes

gene_array=(`awk '{print $4}' ${query_genes}`)

for x in ${gene_array[@]}; do

chrom=`grep "\s${x}$" ${query_genes} | awk '{print $1}'`
start_point=`grep "\s${x}$" ${query_genes} | awk '{print $2}'`
end_point=`grep "\s${x}$" ${query_genes} | awk '{print $3}'`

# Run a loop where GC content is found for specific positions in a reference genome file using GC_script.py and 
# then output with the bed columns defined earlier and with a 5th column for GC content 
# NOTE: The script takes a few minutes (3-5) per gene. Therefore, it takes a long time to use for many genes.
# Hopefully I will find a way to speed it up in the future.

for x in ${gene_sequence[@]}; do
GC_content=`python GC_script.py ../fasta_experimentation/rheMac8.fa ${chrom} ${start_point} ${end_point}`
echo "${chrom}	${start_point}	${end_point}	${x}	$GC_content" >> output.bed
done
