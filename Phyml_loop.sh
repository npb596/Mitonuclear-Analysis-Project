#!/bin/bash

# Load genomics_general module for python scripts in loop

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load genomics_general/24may2019

# Load anaconda module necessary for python scripts in loop

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/2-5.2.0

# Load VCFtools for loop VCF recode

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load vcftools/0.1.14

# Define bed file as query_genes, change bed file as needed.

query_genes=~/mitonuclear_rhesus_files/length_experimentation/NMT_CDS.bed

# Create sequence of numbers equal to the number of genes

gene_number=`grep -c "chr" $query_genes`
gene_number_minus_one=`expr $gene_number - 1`
gene_sequence=(`seq 0 $gene_number_minus_one`)

# Define arrays of gene lengths, gene names, and chromosomes

length_array=(`awk '{ diff = $3 - $2 ; print diff }' $query_genes`)
gene_array=(`awk '{ print $4 }' $query_genes`)
chr_array=(`awk '{ print $1 }' $query_genes`)

# Create trees of genes individually by recoding each gene as a temporary vcf and temporary geno file,
# window sizes are the same as the gene length and minimums are defined based on gene length as well as
# The exon.bed file in the min_length variable and the vcf exon recode files should be changed to match the query gene file
# The python scripts were obtained from https://github.com/simonhmartin/genomics_general and are a necessary part of this script's function
# Location of python scripts and population/individual names should be changed accordingly

for d in ${gene_sequence[@]}; do
    min_length=`grep "[0-9]_${gene_array[$d]}$" ~/mitonuclear_rhesus_files/exon_experimentation/NMT_exons.bed | awk '{diff = $3 - $2; summed += diff} END {print summed}'`
    echo $min_length > tree_minimum.txt
    minimum_variable=`awk '{ minimum = $1/1000 * 5; print minimum }' tree_minimum.txt | sed 's/\.[0-9]*//'`
    cat header.txt > tree_temp_file.bed
    grep "${gene_array[$d]}$" $query_genes >> tree_temp_file.bed
    grep "${gene_array[$d]}$" $query_genes | awk '{print $1,$2,$3}' > tree_temp_coordinates.bed
    vcftools --vcf ~/mitonuclear_rhesus_files/recode_experimentation/recodes_without_suffix/NMT_exons.${chr_array[${d}]} --bed tree_temp_file.bed --recode --out tree_temp_file
    python ~/pop_gen_scripts/genomics_general/parseVCF.py -i tree_temp_file.recode.vcf -o tree_temp_file.geno --skipIndels --skipMono
    python phyml_sliding_windows.py -T 4 -g tree_temp_file.geno -p tree_${gene_array[$d]} -w ${length_array[$d]} -M $minimum_variable --windType predefined --windCoords tree_temp_coordinates.bed --model GTR
done

# Remove extraneous files
rm tree_temp*

# Move files to output folder

mv tree* output/
