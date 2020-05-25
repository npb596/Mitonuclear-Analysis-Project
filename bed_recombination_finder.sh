#!/bin/bash

# Load bedtools on the Alabama Supercomputer
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bedtools/2.26.0

# Create a loop to go through the numbers 1-20 for chromosomes 1-20

for i in {1..20}; do 

# Use the intersectBed command to intersect a bed file of mitonuclear genes and recombination maps for 20 chromosomes (created by Dr. Zachary Szpiech), 
# then extract only the gene names, positions and recombination values from this.

intersectBed -a NMT_genes.bed -b low_chr${i}_rmap.bed -wb | awk '{print $1,$2,$3,$4,$8}' > intermediate_file.txt

# Create an array of gene names only if they are present in the chromosome being looped through at a given time

name_array=(`grep "chr${i}" NMT_genes.bed | awk '{print $4}'`)

# Loop through the array extracting the recombination values corresponding to specific genes and summing them together to give a single value for the gene

for x in ${name_array[@]}; do
grep ${x} intermediate_file.txt | awk '{sum+=$5} END {OFS="\t" ; print $1,$2,$3,$4,sum}' >> intermediate_rates.txt
done

# Removed duplicated genes before adding to a recombination rate file for the chromosome and remove the intermediate files up to this point

sort -u intermediate_rates.txt > chr${i}_recombination_rates.txt
rm intermediate_*

done

#combine all chromosomes into a single file with all mitonuclear recombination rates

cat *recombination_rates.txt >> NMT_recombination_rates_inter.txt

# Use sort_by_chr.sh script written by Dr. Laurie Stevison to sort the file generated from the last step and then remove that file 

~/random-scripts/sort_by_chr.sh ~/mitonuclear_rhesus_files/recode_experimentation/list_chr.txt NMT_recombination_rates_inter.txt NMT_recombination_rates.txt
rm NMT_recombination_rates_inter.txt
