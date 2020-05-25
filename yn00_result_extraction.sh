#!/bin/bash

# This script has been optimized to my dataset but if altered can be used to parse yn00 output in general

# Define an array of the names of all genes in a bed file

filename=$1

gene_names=(`awk '{print $4}' $filename`)

# For every name in this array, grep only the yang-nielsen output to a file and define arrays of this file
# The arrays are the individual macaques being compared and the omega (dN/dS) values

for x in ${gene_names[@]}; do
grep -A 326 "omega" ${x}.yn00 > ${x}_output.txt
individual_one=(`awk 'NR>2 {print $1}' ${x}_output.txt`)
individual_two=(`awk 'NR>2 {print $2}' ${x}_output.txt`)
omega_values=(`awk 'NR>2 {print $7}' ${x}_output.txt`)

# Define a sequence of numbers equal to the amount of comparisons in the file produce above

individual_sequence=`seq 0 325`

# For every value in this sequence find out if the comparisons involve Arctoides individuals (in paml individuals are represented as numbers)
# so Arctoides is 6,7, and 8. If the comparisons are between Arctoides and something not Arctoides, output it to a file.
# If the comparison is between two Arctoides individuals or no Arctoides individuals than don't output files but output "error" messages

for d in $individual_sequence; do
if [[ ${individual_one[$d]} -eq 1 || ${individual_one[$d]} -eq 12 || ${individual_one[$d]} -eq 13 || ${individual_one[$d]} -eq 14 || ${individual_one[$d]} -eq 25 || ${individual_one[$d]} -eq 26 ]]; then
 if [[ ${individual_two[$d]} -eq 7 || ${individual_two[$d]} -eq 8 || ${individual_two[$d]} -eq 9 || ${individual_two[$d]} -eq 20 || ${individual_two[$d]} -eq 21 || ${individual_two[$d]} -eq 22 ]]; then 
  echo -e "${individual_one[$d]}\t${individual_two[$d]}\t${omega_values[$d]}" | grep -v "99.0000" >> ${x}_Sin-Sil_comparison.txt
 else 
  echo -e "${individual_one[$d]}\t${individual_two[$d]}\t${omega_values[$d]}" >> Other_comparisons.txt
 fi
elif [[ ${individual_one[$d]} -eq 7 || ${individual_one[$d]} -eq 8 || ${individual_one[$d]} -eq 9 || ${individual_one[$d]} -eq 20 || ${individual_one[$d]} -eq 21 || ${individual_one[$d]} -eq 22 ]]; then 
 if [[ ${individual_two[$d]} -eq 1 || ${individual_two[$d]} -eq 12 || ${individual_two[$d]} -eq 13 || ${individual_two[$d]} -eq 14 || ${individual_two[$d]} -eq 25 || ${individual_two[$d]} -eq 26 ]]; then 
  echo -e "${individual_one[$d]}\t${individual_two[$d]}\t${omega_values[$d]}" | grep -v "99.0000" >> ${x}_Sin-Sil_comparison.txt
 else 
  echo -e "${individual_one[$d]}\t${individual_two[$d]}\t${omega_values[$d]}" >> Other_comparisons.txt
 fi
fi
done
echo ${x} > ${x}_omega.txt
awk '{sum+=$3} END {print sum/NR}' ${x}_Sin-Sil_comparison.txt >> ${x}_omega.txt
done
cat *_omega.txt > Master_Sin-Sil.dnds.txt
rm *_Sin-Sil_comparison.txt Other_comparisons.txt *_omega.txt *output.txt
