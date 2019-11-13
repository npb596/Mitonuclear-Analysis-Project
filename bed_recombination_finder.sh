#!/bin/bash

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bedtools/2.26.0

for i in {1..20}; do 
intersectBed -a NMT_genes.bed -b low_chr${i}_rmap.bed -wb | awk '{print $4,$8}' > intermediate_file.txt
name_array=(`grep "chr${i}" NMT_genes.bed | awk '{print $4}'`)
for x in ${name_array[@]}; do
grep ${x} intermediate_file.txt | awk '{sum+=$2} END {print $1,sum}' >> intermediate_rates.txt
done
sort -u intermediate_rates.txt > chr${i}_recombination_rates.txt
rm intermediate_*
done
