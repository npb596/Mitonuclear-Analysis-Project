#!/bin/bash

# Load Anaconda and VCFtools modules 

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/3-2019.03 

# Define bed file with gene names at end as query_genes, change bed file as needed.

query_genes=$1

# Define arrays of gene lengths, gene names, and chromosomes

gene_array=(`awk '{print $4}' ${query_genes}`)

for x in ${gene_array[@]}; do

# Run divergence analyses on genes individually. Change output and twisst parameters accordingly.

    python ~/twisst/twisst.py --method complete -t ${x}.trees -w weights_${x} -g Sinica -g Fascicularis -g Arctoides -g Silenus --groupsFile Macaque_groups.tsv --iterations 100
    echo "${x}" >> output.tsv
    cat weights_${x} | grep "[0-9][0-9]" >> output.tsv
done

# Remove extraneous files

rm weights*
