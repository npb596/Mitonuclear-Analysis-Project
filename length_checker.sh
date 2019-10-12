#!/bin/bash

# The following conditional statement will print a help message if -h is entered as an argument.

if [[ $1 == "-h" ]]; then
 echo -e "\nHelp page for length_checker.sh"
 echo -e "\nEnter arguments in the following order:"
 echo -e "\t1. The name of the file with all available gene lengths to query."
 echo -e "\t2. The name of the file with genes you want to match."
 echo -e "\t3. The name of the output file."
 echo -e "\nor just type "-h" to pull up this help page \n"
 exit 0

else

# The following defines variables based on the arguments provided by the user when using this script.

query_genes=$1
genes_to_match=$2
output_file=$3

# The following commands will create a list of numbers 10% less than lengths of the genes to match, a list of numbers 10% greater than lengths of the genes to match
# and a list of the lengths of the query genes.

lower_bound=(`awk '{OFS="\t"} { diff = $3 - $2 ; lower = diff - diff*0.1 ; print lower }' $genes_to_match | sed 's/\.[0-9]//'`)
upper_bound=(`awk '{OFS="\t"} { diff = $3 - $2 ; upper = diff + diff*0.1 ; print upper}' $genes_to_match | sed 's/\.[0-9]//'`)
awk '{OFS="\t"} { diff = $3 - $2 ; print diff }' $query_genes > length.txt

# The following defines a variable as the number of lines in lower lengths to obtain the number of genes the user will want 
# outputted and another variable that subtracts this by 1.

gene_intermediate=`grep -c "chr" $genes_to_match`
gene_number=`expr $gene_intermediate - 1`

# The following command creates a variable that is just a sequence of numbers from 0 to the gene_number variable above. Bash arrays start at 0
# so that is why gene_number subtracts the amount of desired genes by 1.
seq_num=`seq 0 $gene_number`

# The following loop will look through the values in the arrays above, one by one, and find the lengths of query genes
# that fit within the bounds. These query gene lengths will then be placed into separate files corresponding to 
# each gene and then duplicate lengths will be removed for each file.
# Then a single desired gene length will be obtained from each of the query gene files and they will be concatenated into a single
# file, which results in one file with genes roughly corresponding in length to each of the gene lengths.
# Intermediate files are all removed at the end. 
for d in $seq_num; do 
 filename=length.txt
 while read p; do
  if [[ $p -gt ${lower_bound[$d]} && $p -lt ${upper_bound[$d]} ]]
  then
      echo "$p" >> numeros_in_ranges"$d".txt
      sort numeros_in_ranges"$d".txt | uniq > numeros_unique"$d".txt
      shuf -n 1 numeros_unique"$d".txt -o randomized_gene"$d".txt 
  else
      echo ''
  fi
 done < $filename 
done
cat randomized* >> $output_file
rm randomized* numero* length.txt

# The following commands will obtain the full bed information of the genes corresponding to the lengths obtained above. 

awk '{OFS="\t"} { diff = $3 - $2 ; print $1,$2,$3,$4,diff }' $query_genes > length_and_names.txt
readarray output_array < $output_file
for p in ${output_array[@]}; do grep "\s$p$" length_and_names.txt >> "$output_file".inter; done
awk '{print $1,$2,$3,$4}' "$output_file".inter > "$output_file".bed
rm $output_file "$output_file".inter length_and_names.txt
fi
