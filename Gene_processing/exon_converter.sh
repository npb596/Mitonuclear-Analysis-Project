#!/bin/bash

# Define filename as first argument (which should be a bed file without the header), create a text file and subsequent array of just the gene names.

filename=$1
output_file=$2

# Create a variable that is the sequence of numbers from zero to the amount of genes subtracted by one

gene_count=`grep -c "chr" $filename`
gene_count_subtracted=`expr $gene_count - 1`
gene_count_seq=`seq 0 $gene_count_subtracted`

# Define arrays for each column of the bed file 

chrom_number_array=(`awk '{print $1}' $filename`)
exon_start_array=(`awk '{print $2}' $filename`)
exon_end_array=(`awk '{print $3}' $filename`)
locus_name_array=(`awk '{print $4}' $filename`)

# Define an array that is the amount of exons per gene

exon_number_array=(`awk '{print $3}' $filename | sed 's/,/ /g' | awk '{print NF}'`) 

# A loop that will act on every number in the above sequence

for p in $gene_count_seq; do 
 
 # Take the amount of exons and subtract it by one and create a sequence of numbers from zero to this number

 for s in ${exon_number_array[$p]}; do exon_number_array_subtracted=`expr $s - 1`; done
 exon_number_seq=`seq 0 ${exon_number_array_subtracted[@]}` 
 
 # Create files for all exon start and end positions of a gene and then make an array of them

 echo ${exon_start_array[$p]} > start_array_file.txt
 echo ${exon_end_array[$p]} > end_array_file.txt
 final_exon_start_array=(`sed 's/,/ /g' start_array_file.txt`)
 final_exon_end_array=(`sed 's/,/ /g' end_array_file.txt`)
   
 # Print every exon start and end position to another file

   for d in $exon_number_seq; do
    echo -e "${chrom_number_array[$p]}\t${final_exon_start_array[$d]}\t${final_exon_end_array[$d]}\t${d}_${locus_name_array[$p]}" >> $output_file
   done
done

# Remove extraneous files

rm start_array_file.txt end_array_file.txt
