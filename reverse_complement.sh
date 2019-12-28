#!/bin/bash

# Define array of gene names of interest from a bed file
gene_array=(`awk '{print $4}'  ../../../NMT_CDS.bed`)

# For each gene in the array make a variable that states the orientation (positive or negative) of the sequence using
# an annotation downloaded from ENSEMBL. 
for x in ${gene_array[@]}; do
orientation=`grep "${x}" ~/ensGene.combined.txt | awk '{print $4}' | head -n 1`

# Then, if the gene is in negative orientation, use EMBOSS (compiled by user according to instructions in http://emboss.sourceforge.net/download/) 
# to reverse compliment the fasta sequence of that gene.
if [[ $orientation == - ]]; then
~/EMBOSS-6.6.0/emboss/seqret -srev master_${x}.fasta -outseq re_master_${x}.fasta
fi
done
