#!/bin/bash

# Load Anaconda module that contains paml4.8

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/2-5.2.0

# Define a sequence of numbers from 0 to 9, this corresponds to how many individual macaques there are in my analysis
# The above numbers should be changed according to the actual number of individuals in the analysis.
# Define an array of all gene names, file should be changed accordingly

numbers=(`seq 0 9`)
gene_array=(`cat ../NMT_gene_names.txt`) 

# For every gene in the array, go through the corresponding files and replace the stop codons at ends of sequences with 3 dashes

for x in ${gene_array[@]}; do
sed "s/TAG$/---/g" remaster_${x}.fasta.phylip | sed "s/TGA$/---/g" | sed "s/TAA$/---/g" | sed "s/taa$/---/g" | sed "s/tag$/---/g" | sed "s/tga$/---/g" > temp_${x}

# Then take the first line of the resulting file and place it as a header into the expected final phylip file
# Define the second column of this header as the length of the sequences in the file (since this value in a phylip file actually is that)

head -n 1 temp_${x} > phylip_${x}

# Define arrays of the individual organism names and respective gene sequences

individual_name=(`awk 'NR>1 {print $1}' temp_${x}`)
sequence=(`awk 'NR>1 {print $2}' temp_${x}`)

# Now for every number in the first number sequence, echo the monkey names and then sequence names into the final phylip file

for s in ${numbers[@]}; do
echo ${individual_name[${s}]} >> phylip_${x}
echo ${sequence[${s}]} >> phylip_${x}
done

# Run a loop where a yn00 control file is defined for each gene and then run yn00 on the control file, this finishes the original loop defined at the beginning 
# yn00 is a program available as part of paml4.8, developed by Ziheng Yang, and based on a dN/dS model developed by Ziheng Yang and Rasmus Nielsen in 2000.

sed -e "s/JEAN/${x}/g" yn00_permanent.ctl > yn00.ctl
yn00
done

# Remove extraneous files
rm phylip_* temp_*
