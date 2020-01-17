#!/bin/bash

# Introduction: This script will take two bed files containing recombination rates, GC content, and gene length values on top of standard bed information.
# The gene lengths can be found by simply subtracting bed start from end values. Finding recom rate and GC content will require additional scripts and data.
# This script is fairly repetetive, the mechanism of matching recom rate to recom rate and GC content to GC content is pretty similar. Thus, there are
# three seperate but very similar python scripts involved (named match_math[1-3].py). More specific detail can be found in each of these scripts
# but in essence they define a range of values for a matching gene to fit within. Recom rates, GC content, and gene lengths are orders of magnitude apart
# so different stringencies are required, hence marginally different python scripts.
# Additionally, this script will generate ordered "contender" and "unique" files, where the former define increasingly similar matching genes and the latter 
# simply remove duplicates from the former. A good way to be sure the script is functioning correctly is to NOT remove the files (as the script does by default)
# and then line count each contender file. Since they become increasingly stringent the lines should decrease (line count of: fourth_unique.txt < third_unique.txt etc.)

# Input bed files of genes to query and genes you want to match against query genes. These files must already contain 
# recombination rate, GC and gene (or CDS) length in the columns defined by awk below
query_genes=chrX_NMT_data.bed
matching_genes=chrX_nuclear_data.bed

# While reading the query genes do...
while read p; do

# Print each line into an individualized bed file for only one gene. Then define arrays for the chromosome, recombination rate, GC content and gene (or CDS) length.
echo ${p} > individual_gene.bed

query_chr=`awk '{print $1}' individual_gene.bed`
query_recom=`awk '{print $5}' individual_gene.bed`
query_GC=`awk '{print $6}' individual_gene.bed`
query_length=`awk '{print $7}' individual_gene.bed`

# Search the matching genes file for all genes on the same chromosome as the current query gene and output them to a temporary file.
grep "^${query_chr}\s" ${matching_genes} > first_contender.txt
	
# Define an array for recombination rates of the file generated above and determine if each recombination rate is sufficiently
# close to the query recombination rate using the associated "match_math_1.py" script. Qualifying genes will be output to another temporary file.
match_recom=(`awk '{print $5}' first_contender.txt`)
	for x in ${match_recom[@]}; do
		python match_math_1.py ${query_recom} ${x} > temp_file.txt
		grep -f temp_file.txt first_contender.txt >> second_contender.txt
	done
	sort -u second_contender.txt > second_unique.txt

# Define an array for GC content of the file generated above and determine if each GC content value is sufficiently
# close to the query GC content value  using the associated "match_math_2.py" script. Qualifying genes will be output to another temporary file.
match_GC=(`awk '{print $6}' second_unique.txt`)
	for y in ${match_GC[@]}; do
		python match_math_2.py ${query_GC} ${y} > temp_file.txt
		grep -f temp_file.txt second_unique.txt >> third_contender.txt
	done
	sort -u third_contender.txt > third_unique.txt	

# Define an array for gene length of the file generated above and determine if each gene length is sufficiently
# close to the query gene length using the associated "match_math_3.py" script. Qualifying genes will be output to another temporary file.
# Note: The value contained in temp_file.txt is not directly sought after because it is generally a four/five digit number that could randomly
# be found inside the several-digit numbers preceding it in a bed file. Thus this loop requires the number be found at the end of the line 
# and that it is not preceded by any other numbers.
match_length=(`awk '{print $7}' third_unique.txt`)
	for z in ${match_length[@]}; do
		python match_math_3.py ${query_length} ${z} > temp_file.txt
		final_length=(`cat temp_file.txt`)
		grep "\s${final_length}$" third_unique.txt >> fourth_contender.txt
	done
	sort -u fourth_contender.txt > fourth_unique.txt

# Take one random value from the file generated above and output it into the final output file. Now you have one gene matching your query gene.
shuf -n 1 fourth_unique.txt >> final_nuclear_gene_matching.bed

# Remove extraneous files within and outside the loop.
# IMPORTANT NOTE: Check that these names will not include files already in your directory.
rm *contender.txt *unique.txt 
done < ${query_genes}
rm temp_file.txt individual_gene.bed
