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

# Run divergence analyses on genes individually by recoding each gene as a temporary vcf and temporary geno file,
# window sizes are the same as the gene length and minimums are defined based on gene length as well
# The python scripts were obtained from https://github.com/simonhmartin/genomics_general and are a necessary part of this script's function
# Location of python scripts and population/individual names should be changed accordingly 

    echo ${length} > temp_file_minimum.txt
    minimum=`awk '{ minimum = $1/1000 * 3; print minimum }' temp_file_minimum.txt | sed 's/\.[0-9]*//'` 
    cat header.txt > temp_file.bed
    grep "\s${x}$" $query_genes >> temp_file.bed
    grep "\s${x}$" $query_genes | awk '{OFS="\t" ; print $1,$2,$3}' > temp_file_windcoord.bed
    vcftools --gzvcf ../recode_experimentation/mitonuclear_CDS_recodes/nmt_cds.${chrom} --bed temp_file.bed --recode --out temp_file
    python ~/genomics_general/VCF_processing/parseVCF.py -i temp_file.recode.vcf -o temp_file.geno --skipIndels --skipMono 
    python ~/genomics_general/popgenWindows.py -g temp_file.geno -o divergence_${x} -f phased --windType predefined --windCoords temp_file_windcoord.bed -w ${length} -m ${minimum} -p Sinica Tibetan-macaque-NO.3,XH1.Assamensis,A20 -p Fascicularis BGI-96346.mulatta,BGI-CE-4.fascicularis,BGI.Mulatta-1,CR-5.Mulatta-2 -p Arctoides Malaya,SM1.Arctoides-1,SM2.Arctoides-2 -p Silenus PM664.nemestrina,PF660.nigra,PM592.tonkeana -T 4
done

# Create a .csv file that includes all measurements obtained
# The name of introgression_GENENAME should be changed to make sure it's a gene actually present in the measurements
# and final output name should accordingly be changed to match the query gene file

grep "scaffold" divergence_GENENAME > output.csv
cat divergence* > divergence_measurements_temp
grep -v "scaffold" divergence_measurements_temp >> output.csv

# Remove extraneous files

rm divergence* temp_file*
