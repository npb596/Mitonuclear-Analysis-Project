#!/bin/bash

# Load python module on ASC necessary for GC_script.py

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load python/2.7.9

# Define input bed file to determine GC content of genes within and define arrays of each column

filename=truest_nuclear_CDS.bed

chrom=(`awk '{print $1}' $filename`)
start_points=(`awk '{print $2}' $filename`)
end_points=(`awk '{print $3}' $filename`)
gene_names=(`awk '{print $4}' $filename`)

# Define array of 0 to the number of genes inside the bed file subtracted by one, which ultimately is the number of genes in the bed file

gene_number=`grep -c "chr" $filename`
gene_number_minus_one=`expr $gene_number - 1`
gene_sequence=(`seq 0 $gene_number_minus_one`)

# Run a loop where GC content is found for specific positions in a reference genome file using GC_script.py and 
# then output with the bed columns defined earlier and with a 5th column for GC content 
# NOTE: The script takes a few minutes (3-5) per gene. Therefore, it takes a long time to use for many genes.
# Hopefully I will find a way to speed it up in the future.

for x in ${gene_sequence[@]}; do
GC_content=`python GC_script.py ../fasta_experimentation/rheMac8.fa ${chrom[${x}]} ${start_points[${x}]} ${end_points[${x}]}`
echo "${chrom[${x}]}	${start_points[${x}]}	${end_points[${x}]}	${gene_names[${x}]}	$GC_content" >> nuclear_GC.bed
done
