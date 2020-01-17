#!/bin/bash

# INTRODUCTION: This script is intended for use with Ensembl ensGene CDS and exon annotations. So "CDS" and "exon" bed files referred to from this point on
# are bed files made from those annotation positions. Generally, the actual protein-coding portion of a gene is a combination of the positions
# given for both CDS and exons so this script is designed to combine the positions into one file. It is only effective up to a certain point.
# When it operates on an exon bed file you will notice that the file should have extra commas left over, these will indicate
# where this script has operated but also shows the parts that must be manually edited to yield accurate CDS/exon combinations.
# As the script currently operates, it is normal for it to generate awk and sed errors.

# Define array of gene names. Replace bed file as needed.
gene_array=(`awk '{print $4}' matching_nuclear_CDS.bed`)

# For every gene in the array do...
for s in ${gene_array[@]}; do

# Define CDS and first exon start and end positions.
CDS_start=`grep "\s${s}$" matching_nuclear_CDS.bed | awk '{print $2}'`
CDS_end=`grep "\s${s}$" matching_nuclear_CDS.bed | awk '{print $3}'`
exon_start=`grep "\s${s}$" matching_nuclear_Exons_inter.bed | awk '{print $2}' | awk 'BEGIN {FS=","} {print $1}'`
exon_end=`grep "\s${s}$" matching_nuclear_Exons_inter.bed | awk '{print $3}' | awk 'BEGIN {FS=","} {print $(NF-1)}'`

# Replace the first exon start and end positions with the respective CDS start and end positions.
sed -i "s/${exon_start}/${CDS_start}/" matching_nuclear_Exons_inter.bed
sed -i "s/${exon_end}/${CDS_end}/" matching_nuclear_Exons_inter.bed

# Define various positions. In order: start of second exon, start of third exon, start of fourth exon, start of fifth exon, end of second-to-last exon,
# end of third-to-last exon, end of first exon, end of second exon.
exon_second_start=`grep "\s${s}$" matching_nuclear_Exons_inter.bed | awk '{print $2}' | awk 'BEGIN {FS=","} {print $2}'`
exon_third_start=`grep "\s${s}$" matching_nuclear_Exons_inter.bed | awk '{print $2}' | awk 'BEGIN {FS=","} {print $3}'`
exon_fourth_start=`grep "\s${s}$" matching_nuclear_Exons_inter.bed | awk '{print $2}' | awk 'BEGIN {FS=","} {print $4}'`
exon_fifth_start=`grep "\s${s}$" matching_nuclear_Exons_inter.bed | awk '{print $2}' | awk 'BEGIN {FS=","} {print $5}'`
exon_second_from_end=`grep "\s${s}$" matching_nuclear_Exons_inter.bed | awk '{print $3}' | awk 'BEGIN {FS=","} {print $(NF-2)}'`
exon_third_from_end=`grep "\s${s}$" matching_nuclear_Exons_inter.bed | awk '{print $3}' | awk 'BEGIN {FS=","} {print $(NF-3)}'`
exon_first_end=`grep "\s${s}$" matching_nuclear_Exons_inter.bed | awk '{print $3}' | awk 'BEGIN {FS=","} {print $1}'`
exon_second_end=`grep "\s${s}$" matching_nuclear_Exons_inter.bed | awk '{print $3}' | awk 'BEGIN {FS=","} {print $2}'`

# Run various conditional statements determining if the CDS positions contradict the exon positions defined above. If so the exon positions are removed.
if [[ $CDS_start > $exon_second_start ]]; then
	sed -i "s/${exon_second_start}//" matching_nuclear_Exons_inter.bed
fi

if [[ $CDS_start > $exon_third_start ]]; then
	sed -i "s/${exon_third_start}//" matching_nuclear_Exons_inter.bed
fi

if [[ $CDS_start > $exon_fourth_start ]]; then
        sed -i "s/${exon_fourth_start}//" matching_nuclear_Exons_inter.bed
fi

if [[ $CDS_start > $exon_fifth_start ]]; then
        sed -i "s/${exon_fifth_start}//" matching_nuclear_Exons_inter.bed
fi

if [[ $CDS_start > $exon_first_end ]]; then
        sed -i "s/${exon_first_end}//" matching_nuclear_Exons_inter.bed
fi

if [[ $CDS_start > $exon_second_end ]]; then
        sed -i "s/${exon_second_end}//" matching_nuclear_Exons_inter.bed
fi

if [[ $CDS_end < $exon_second_from_end && $exon_second_from_end != *","* ]]; then
	sed -i "s/${exon_second_from_end}//" matching_nuclear_Exons_inter.bed
fi

if [[ $CDS_end < $exon_third_from_end && $exon_second_from_end != *","* ]]; then
        sed -i "s/${exon_third_from_end}//" matching_nuclear_Exons_inter.bed
fi

done

# When this script is run on the Alabama Supercomputer, multiple "sed" files are generated and error messages are given. I do not know why this is.
# The below commands will give the user permission to delete all these files and then delete them.
chmod +rwx sed*
rm sed*
