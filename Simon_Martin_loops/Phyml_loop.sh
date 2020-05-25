#!/bin/bash

# Load Anaconda and VCFtools modules 

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/2-5.2.0
module load vcftools/0.1.14

# Define bed file with gene names at end as query_genes, change bed file as needed.

query_genes=$1

# Define arrays of gene lengths, gene names, and chromosomes

gene_array=(`awk '{print $4}' ${query_genes}`)

for x in ${gene_array[@]}; do

length=`grep "\s${x}$" ${query_genes} | awk '{diff = $3-$2 ; print diff}'`
chrom=`grep "\s${x}$" ${query_genes} | awk '{print $1}'`

# Create trees of genes individually by recoding each gene as a temporary vcf and temporary geno file,
# window sizes are the same as the gene length and minimums are defined based on gene length as well as
# The exons.bed file in the min_length variable and the vcf exon recode files should be changed to match the query gene file
# The python scripts were obtained from https://github.com/simonhmartin/genomics_general and are a necessary part of this script's function
# Location of python scripts and population/individual names should be changed accordingly

    min_length=`grep "[0-9]_${x}$" exons.bed | awk '{diff = $3 - $2; summed += diff} END {print summed}'`
    echo $min_length > tree_minimum.txt
    minimum_variable=`awk '{ minimum = $1/1000 * 5; print minimum }' tree_minimum.txt | sed 's/\.[0-9]*//'`
    cat header.txt > tree_temp_file.bed
    grep "\s${x}$" $query_genes >> tree_temp_file.bed
    grep "]s${x}$" $query_genes | awk '{print $1,$2,$3}' > tree_temp_coordinates.bed
    vcftools --vcf ~/mitonuclear_rhesus_files/recode_experimentation/recodes_without_suffix/NMT_exons.${chrom} --bed tree_temp_file.bed --recode --out tree_temp_file
    python ~/pop_gen_scripts/genomics_general/parseVCF.py -i tree_temp_file.recode.vcf -o tree_temp_file.geno --skipIndels --skipMono
    python phyml_sliding_windows.py -T 4 -g tree_temp_file.geno -p ${x} -w ${length} -M ${minimum_variable} --windType predefined --windCoords tree_temp_coordinates.bed --model GTR

done

# Remove extraneous files
rm tree_temp*
