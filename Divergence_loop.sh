#!/bin/bash

# Load Anaconda module for python script in loop

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/2-5.2.0

# Load VCFtools for loop recode

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load vcftools/0.1.14

# Define bed file as query_genes, change bed file as needed.

query_genes=random_nuc_CDS_tabbed.bed

# Create sequence of numbers equal to the number of genes

gene_number=`grep -c "chr" $query_genes`
gene_number_minus_one=`expr $gene_number - 1`
gene_sequence=`seq 0 $gene_number_minus_one`

# Define arrays of gene lengths, gene names, and chromosomes

length_array=`awk '{ diff = $3 - $2 ; print diff }' $query_genes`
gene_array=`awk '{ print $4 }' $query_genes`
chr_array=`awk '{ print $1 }' $query_genes`

# Run divergence analyses on genes individually by recoding each gene as a temporary vcf and temporary geno file,
# window sizes are the same as the gene length and minimums are defined based on gene length as well
# The exon.bed file in the min_length variable and the vcf exon recode files should be changed to match the query gene file
# The python scripts were obtained from https://github.com/simonhmartin/genomics_general and are a necessary part of this script's function
# Location of python scripts and population/individual names should be changed accordingly 

for d in $gene_sequence; do
    min_length=`grep "[0-9]_${gene_array[$d]}$" ../exon_experimentation/random_nuclear_exons.bed | awk '{diff = $3 - $2; summed += diff} END {print summed}'`
    echo $min_length > divergence_minimum.txt
    minimum_variable=`awk '{ minimum = $1/1000 * 5; print minimum }' divergence_minimum.txt | sed 's/\.[0-9]*//'` 
    cat header.txt > temp_file.bed
    grep "${gene_array[$d]}$" $query_genes >> temp_file.bed
    vcftools --vcf ../recode_experimentation/nuclear_exons.${chr_array[$d]} --bed temp_file.bed --recode --out temp_file
    python ~/pop_gen_scripts/genomics_general/parseVCF.py -i temp_file.recode.vcf -o temp_file.geno --skipIndels --skipMono 
    python ~/pop_gen_scripts/genomics_general/popgenWindows.py -g temp_file.geno -o divergence_${gene_array[$d]} -f phased -w ${length_array[$d]} -m ${minimum_variable} -p Sinica Tibetan-macaque-NO.3,XH1.Assamensis,A20 -p Fascicularis BGI-96346.mulatta,BGI-CE-4.fascicularis,BGI.Mulatta-1,CR-5.Mulatta-2 -p Arctoides Malaya,SM1.Arctoides-1,SM2.Arctoides-2 -p Baboon papAnu4
done

# Create a .csv file that includes all measurements obtained
# The name of introgression_GENENAME should be changed to make sure it's a gene actually present in the measurements
# and final output name should accordingly be changed to match the query gene file

grep "scaffold" divergence_GTSE1 > Divergence_nuclear_measurements.csv
cat divergence* > divergence_measurements_temp
grep -v "scaffold" divergence_measurements_temp >> Divergence_nuclear_measurements.csv

# Remove extraneous files

rm divergence* temp_file*
