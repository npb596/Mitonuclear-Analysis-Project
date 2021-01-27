# Mitonuclear Analysis Project
This repository contains data files and scripts to reproduce the analysis for the manuscript titled: **"Mitonuclear conflict in a macaque species exhibiting phylogenomic discordance"**

This repository will continue to be updated with better-optimized scripts accomplishing the purposes of this experiment. In particular, many of the shell scripts may be better optimized as python scripts. The version with all scripts originally used in writing this manuscript are available as the version 1.0 release on the github page.

# Abstract 

Speciation and hybridization are intertwined processes in the study of evolution. Hybridization between sufficiently diverged populations can result in genomic conflict within offspring, causing reduced viability and fertility, thus increasing divergence between populations. Conflicts between mitochondrial and nuclear genes are increasingly found to play a role in this process in various systems. We examine the possibility of this conflict in the bear macaque, Macaca arctoides, a primate species exhibiting mitonuclear discordance due to extensive hybridization with species in the sinica and fascicularis groups. Here, divergence, introgression, and natural selection of mitonuclear genes (N = 143) relative to nuclear control genes (N = 127) were analyzed to determine if there are evolutionary processes involved in resolving the potential conflict caused by mitonuclear discordance. Nucleotide divergence of mitonuclear genes is increased relative to nuclear genes between M. arctoides and the species sharing its nuclear ancestry (P = 0.028), consistent with genetic conflict. However, measures of introgression and selection do not identify large-scale co-introgression or co-evolution as means to resolve mitonuclear conflict. The methodology implemented provides a framework that can be used to examine the effects of mitonuclear co-introgression and co-evolution on a genomic scale in a variety of systems.

# Data Analysis

The directory ["Gene-Processing-Scripts"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/tree/master/Gene-Processing-Scripts)
contains scripts used to convert gene data between different formats and to obtain "matching" data for N-mt and nuclear genes. Dependendencies and authorship of source scripts are stated within scripts.

The directory ["Pop-Gen-Scripts"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/tree/master/Pop-Gen-Scripts) contains scripts used to run divergence, introgression, and selection analyses on N-mt and nuclear genes. These are all bash scripts designed to loop through programs computing these statistics.

A single data file ["Gene_Data.csv"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/blob/master/Gene_Data.csv) is provided that contains all information gathered from scripts in the above directories. A single R script ["Mitonuclear\_Stats\_Plots.R"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/blob/master/Mitonuclear_Stats_Plots.R) is provided that will analyze this data file and produce Figures and text files with statistical analyses and output these respectively to ["Figures"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/tree/master/Figures) and ["Stats"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/tree/master/Stats).

# Usage of Gene Processing Scripts

## Matching N-mt genes to control nuclear genes

N-mt genes were matched to control nuclear genes based if they were present on the same chromosome and had similar lengths, recombination rates, and GC content. Therefore, this should be obtained for the N-mt genes of interest and a large database of other nuclear genes (in this case, all Ensembl annotated genes not included in this study).

Chromosome information is given in the Ensembl annotation. Gene start and end positions are also given in the Ensembl annotation so the length of a gene can be given with an awk command such as (assuming a format such as bed where the gene start is column 2 and gene end is column 3):

```
awk '{diff=$3-$2 ; print diff}' filename.bed
```

The scripts GC\_loop.sh and GC\_script.py have been written to obtain GC content per geen for this study.

GC\_loop.sh requires a bed file as an argument and will then use the script GC\_script.py in a loop to find the GC content of those genes in a reference fasta file. The path to GC\_script.py and the reference fasta must be correctly defined within GC\_loop.sh

Recombination rates can be obtained from a bigwig file such as: ftp://ftp.hgsc.bcm.edu/ucscHub/rhesusSNVs/rheMac8/all.rate.bw

Then this can be lifted over to a bed file using UCSC liftOver. Given a bed file with gene coordinates and a bed file with recombination rates a new output file can be created using bedtools map with the gene coordinates and averaged recombination rate.

```
bedtools map -a NMT.bed -b RecomRate_CM.bed -c 5 -o mean > output.bed
```
Once all these parameters are found one can create a bed file with with the following tab-separated values by combining all the data:

```
chromosome	start	end	gene	recom_rate	GC_content	length
chr1    15716094        15749833        SDHB    0.292168        44.257387593    33739
```

The header is included here so you can know what the parameters are and the order they are in. gene\_matcher.py assumes there is no header so it would be best to make the file without the header as well. Once files of this sort are made for both the query genes (N-mt genes in our case) and a large set of potentially matching genes (subject or other nuclear genes in our case) then one can use the query file as first argument, subject file as second argument, and desired final output as third argument.

```
python3 gene_matcher.py Query.bed Subject.bed Output.bed
```

The gene\_matcher.py script uses the python random module to pick a random matching gene if multiple match all the criteria of the query gene. This necessarily adds randomness to the study. Additionally, it is specified in the manuscript that recombination rate is not used to match chromosome X genes, so this script can be modified by simply commenting out the relevant portions.

## Converting Ensembl exon coordinates to bed coordinates

The script exon\_converter.py is provided to take exon coordinates in Ensembl format and convert them to a bed file. For example, given a transcript with four exons where the exon positions are 100-150, 200-250, 300-350, and 400-450, Ensembl will give the following two tab-separated columns:

```
100,200,300,400,	150,250,350,450,
```
We want this to be represented in bed format such as:

```
100	150
200	250
300	350
400	450
```
Therefore, we can convert the Ensembl data to bed data with a command like the following:

```
./exon_converter.py ensGene.combined.txt Exons.bed
```

The second argument is the desired output file and the first argument is an Ensembl annotation file where HUGO gene names have been added using an awk command:

```
awk 'NR==FNR {file1[$2]=$0; next} {$1=file1[$1]; print}' ensGene.txt ensemblToGeneName.txt > ensGene.combined.txt 
```

Both the ensGene.txt and ensemblToGeneName.txt files can be found from the UCSC Genome Browser and the Ensembl site for a given reference genome (in this study rheMac8). The output of the exon\_converter.py will look something like this:

```
chr1    8230    10763   ENSMMUT00000032773.4    PGBD2
chr1    13490   13554   ENSMMUT00000032773.4    PGBD2
```

Where the positions are exon coordinates and the fourth and fifth columns are respectively transcript and gene names.

The Ensembl annotation file also provides CDS coordinates for each transcript. So for someone interested in obtaining the true protein-coding positions of a gene, these must be intersected with the given exon coordinates, either before or after the above conversion process. 
