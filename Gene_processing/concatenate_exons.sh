#!/bin/bash

# Define arrays based on Macaque individual names and the names of N-mt genes.

indv_file=$1
gene_file=$2

monkey_array=(`cat ${indv_file}`)
gene_array=(`awk '{print $4}' ${gene_file}`)

# Start loops for every gene name and for every fasta file in the folder (ordered numerically), which are named after the genes

for s in ${gene_array[@]}; do

exon_array=(`ls *${s}.fa | sort -n`) 

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
		grep -A $final_number "${x}-1" $filename | grep -v "${x}-1" >> ${x}-1_${s}.temp
		echo ">${x}_A" > ${x}-0_${s}_final.temp
		echo ">${x}_B" > ${x}-1_${s}_final.temp
                cat ${x}-0_${s}.temp >> ${x}-0_${s}_final.temp
		cat ${x}-1_${s}.temp >> ${x}-1_${s}_final.temp
                done
	done
cat *${s}_final.temp >> master_${s}.fasta
done
rm *temp
