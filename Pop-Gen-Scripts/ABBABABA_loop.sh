#!/bin/bash

# Load Anaconda and VCFtools modules 

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/2-5.2.0
module load vcftools/0.1.14

# Define bed file with gene names at end as query_genes, change bed file as needed.

query_genes=$1
vcf_file=$2
output_file=$3

# Define arrays of gene lengths, gene names, and chromosomes

gene_array=(`awk '{print $4}' ${query_genes}`)

for x in ${gene_array[@]}; do

length=`grep "\s${x}$" ${query_genes} | awk '{diff = $3-$2 ; print diff}'`
chrom=`grep "\s${x}$" ${query_genes} | awk '{print $1}'`

# Run introgression analyses on genes individually by recoding each gene as a temporary vcf and temporary geno file,
# window sizes are the same as the gene length and minimums are defined based on gene length as well
# The python scripts were obtained from https://github.com/simonhmartin/genomics_general and are a necessary part of this script's function
# Location of python scripts and population/individual names should be changed accordingly

    echo ${length} > FDM_minimum.txt
    fdM_minimum_variable=`awk '{ minimum = $1/1000 * 3 ; print minimum }' FDM_minimum.txt | sed 's/\.[0-9]*//'`
    cat header.txt > FDM_temp_file.bed
    grep "${x}$" $query_genes >> FDM_temp_file.bed
    grep "${x}$" $query_genes | awk '{OFS="\t" ; print $1,$2,$3}' > FDM_windcoords.bed
    vcftools --gzvcf ${vcf_file} --bed FDM_temp_file.bed --recode --out FDM_temp_file
    python ~/genomics_general/VCF_processing/parseVCF.py -i FDM_temp_file.recode.vcf -o FDM_temp_file.geno --skipIndels --skipMono
    python ~/genomics_general/ABBABABAwindows.py -g FDM_temp_file.geno -o FDM_${gene_array[$d]} -f phased --windType predefined --windCoords FDM_windcoords.bed -w ${length_array[$d]} -m ${fdM_minimum_variable} -P1 Sinica Tibetan-macaque-NO.3,XH1.Assamensis,A20 -P2 Fascicularis BGI-96346.mulatta,BGI-CE-4.fascicularis,BGI.Mulatta-1,CR-5.Mulatta-2 -P3 Arctoides Malaya,SM1.Arctoides-1,SM2.Arctoides-2 -O Silenus PM664.nemestrina,PF660.nigra,PM592.tonkeana -T 4
done

# Create a .csv file that includes all measurements obtained
# The name of introgression_GENENAME should be changed to make sure it's a gene actually present in the measurements
# and final output name should accordingly be changed to match the query gene file

grep "scaffold" FDM_GENENAME > ${output_file}
cat introgression* > FDM_measurements_temp
grep -v "scaffold" FDM_measurements_temp >> ${output_file}

# Remove extraneous files

rm FDM_*
