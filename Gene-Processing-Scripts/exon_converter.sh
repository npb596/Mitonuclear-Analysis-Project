#!/bin/bash

# Define filename as first argument (which should be a bed file with tab separated fields without the header), create a text file and subsequent array of just the gene names.

filename=$1
output_file=$2

gene_name=(`awk '{print $4}' $filename`)

# Define an array that is the amount of exons per gene

for p in ${gene_name[@]}; do

chrom=`grep "\s${p}$" $filename | awk '{print $1}'`
exon_start_array=(`grep "\s${p}$" $filename | awk '{print $2}' | sed "s/,/\t/g"`)
exon_end_array=(`grep "\s${p}$" $filename | awk '{print $3}' | sed "s/,/\t/g"`)
exon_number=`grep "\s${p}$" $filename | awk '{print $3}' | sed "s/,/\t/g" | awk '{print NF}'`  


# Take the amount of exons and subtract it by one and create a sequence of numbers from zero to this number

 exon_number_array_subtracted=`expr $exon_number - 1`
 exon_number_seq=`seq 0 ${exon_number_array_subtracted[@]}` 
 
 # Print every exon start and end position to another file

   for d in $exon_number_seq; do
	echo -e "${chrom}\t${exon_start_array[$d]}\t${exon_end_array[$d]}\t${d}_${p}" >> $output_file
   done
done

