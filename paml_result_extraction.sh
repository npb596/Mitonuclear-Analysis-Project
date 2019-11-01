#!/bin/bash

# This script has been highly optimized to my dataset but if altered can be used to parse yn00 output in general

# Define an array of the names of all mitonuclear genes

gene_names=(`cat ../NMT_gene_names.txt`)

# For every name in this array, grep only the yang-nielsen output to a file and define arrays of this file
# The arrays are the individual macaques being compared and the omega (dN/dS) values

for x in ${gene_names[@]}; do
grep -A 46 "omega" output_${x} > grepped_output.txt
individual_one=(`awk 'NR>1 {print $1}' grepped_output.txt`)
individual_two=(`awk 'NR>1 {print $2}' grepped_output.txt`)
omega_values=(`awk 'NR>1 {print $7}' grepped_output.txt`)

# Define a sequence of numbers equal to the amount of comparisons in the file produce above

individual_sequence=`seq 0 44`

# For every value in this sequence find out if the comparisons involve Arctoides individuals (in paml individuals are represented as numbers)
# so Arctoides is 6,7, and 8. If the comparisons are between Arctoides and something not Arctoides, output it to a file.
# If the comparison is between two Arctoides individuals or no Arctoides individuals than don't output files but output "error" messages

for d in $individual_sequence; do
if [[ ${individual_one[$d]} -eq 6 || ${individual_one[$d]} -eq 7 || ${individual_one[$d]} -eq 8 ]]; then
 if [[ ${individual_two[$d]} -eq 6 || ${individual_two[$d]} -eq 7 || ${individual_two[$d]} -eq 8 ]]; then
  echo "Arctoides comparison" 
  echo -e "${individual_one[$d]}\t${individual_two[$d]}\t${omega_values[$d]}" >> Arc-Arc_comparison.txt
 else 
  echo -e "${individual_one[$d]}\t${individual_two[$d]}\t${omega_values[$d]}" >> Arctoides_comparisons.txt
  echo "solid"
 fi
elif [[ ${individual_two[$d]} -eq 6 || ${individual_two[$d]} -eq 7 || ${individual_two[$d]} -eq 8 ]]; then
 if [[ ${individual_one[$d]} -eq 6 || ${individual_one[$d]} -eq 7 || ${individual_one[$d]} -eq 8 ]]; then
  echo "nonsense comparison"
  echo -e "${individual_one[$d]}\t${individual_two[$d]}\t${omega_values[$d]}" >> Arc-Arc_comparison.txt
 else 
  echo -e "${individual_one[$d]}\t${individual_two[$d]}\t${omega_values[$d]}" >> Arctoides_comparisons.txt
  echo "solid"
 fi
else 
 echo -e "${individual_one[$d]}\t${individual_two[$d]}\t${omega_values[$d]}" >> Sin_Fas_comparison_temp.txt
fi
done
done

individual_one_again=(`awk '{print $1}' Arctoides_comparisons.txt`)
individual_two_again=(`awk '{print $2}' Arctoides_comparisons.txt`)
omega_values_again=(`awk '{print $3}' Arctoides_comparisons.txt`)

individual_number=`grep -c "[0-9]" Arctoides_comparisons.txt` 
individual_sequence_again=`seq 0 ${individual_number}`

for d in $individual_sequence_again; do
if [[ ${omega_values_again[$d]} != 99.0000 ]]; then
 if [[ ${individual_one_again[$d]} -eq 6 || ${individual_one_again[$d]} -eq 7 || ${individual_one_again[$d]} -eq 8 || ${individual_one_again[$d]} -eq 9 || ${individual_one_again[$d]} -eq 10 ]]; then
  if [[ ${individual_two_again[$d]} -eq 1 || ${individual_two_again[$d]} -eq 6 || ${individual_two_again[$d]} -eq 7 || ${individual_two_again[$d]} -eq 8 || ${individual_two_again[$d]} -eq 9 || ${individual_two[$d]} -eq 10 ]]; then
   echo "Sinica comparison"
   echo -e "${omega_values_again[$d]}" >> Sin_Arc_comparisons.txt
  else 
   echo "Fascicularis comparison"
   echo -e "${omega_values_again[$d]}" >> Fas_Arc_comparisons.txt
  fi
 else
  echo "Nothing"
 fi
else
 echo "Nothing"
fi
done

omega_values_sinfas=(`awk '{print $3}' Sin_Fas_comparison_temp.txt`)
individual_number_sinfas=`grep -c "[0-9]" Sin_Fas_comparison_temp.txt`
individual_sequence_sinfas=`seq 0 ${individual_number}`

for d in $individual_sequence_again; do
if [[ ${omega_values_sinfas[$d]} != 99.0000 ]]; then
 echo -e "${omega_values_sinfas[$d]}" >> Sin_Fas_comparisons.txt  
else
 echo "nothing"
fi
done
