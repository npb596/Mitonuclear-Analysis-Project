#!/bin/bash

# Define arrays based on Macaque individual names and the names of N-mt genes.

monkey_array=(`cat Macaque_names.txt`)
gene_array=(`awk '{print $4}' ../NMT_CDS.bed`)

# Start loops for every gene name and for every fasta file in the folder (ordered numerically), which are named after the genes

for s in ${gene_array[@]}; do

exon_array=(`ls *$s | sort -n`) 

	for filename in ${exon_array[@]}
	do

# Take the number of lines for each fasta file and subtract this number by 26 then divide it by 26.
# 26 is the number of individuals per file (13 individuals but split into two fasta sequences because they are diploid)
# So all this math gives us the amount of sequence lines per exon

	file_count=`wc -l $filename | awk '{print $1}'`
	minus_twenty=`expr $file_count - 26`
	final_number=`expr $minus_twenty / 26`

# Start a loop for each macaque individual name and create files with all exons of a particular gene for each individual macaque

		for x in ${monkey_array[@]}; do 
                grep -A $final_number "${x}-0" $filename | grep -v "${x}-0" >> ${x}-0_${s}.temp
		echo ">${s}_${x}-0" > ${x}-0_${s}_final.temp
                cat ${x}-0_${s}.temp >> ${x}-0_${s}_final.temp
                done
	done
cat *${s}_final.temp >> master0_${s}.fasta
done
rm *temp
