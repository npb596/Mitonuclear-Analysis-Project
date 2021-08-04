# Mitonuclear Analysis Project
This repository contains data files and scripts to reproduce the analysis for the manuscript titled: **"Mitonuclear conflict in a macaque species exhibiting phylogenomic discordance"**

This repository will continue to be updated with better-optimized scripts accomplishing the purposes of this experiment. In particular, many of the shell scripts may be better optimized as python scripts. The version with all scripts originally used in writing this manuscript are available as the version 1.0 release on the github page.

# Abstract 

Speciation and hybridization are intertwined processes in the study of evolution. Hybridization between sufficiently diverged populations can result in genomic conflict within offspring, causing reduced viability and fertility, thus increasing divergence between populations. Conflicts between mitochondrial and nuclear genes are increasingly found to play a role in this process in various systems. We examine the possibility of this conflict in the bear macaque, *Macaca arctoides* (Primates: Cercopithecidae), a primate species exhibiting mitonuclear discordance due to extensive hybridization with species in the sinica and fascicularis groups. Here, divergence, introgression, and natural selection of mitonuclear genes (N = 160) relative to nuclear control genes (N = 144) were analyzed to determine if there are evolutionary processes involved in resolving the potential conflict caused by mitonuclear discordance. Nucleotide divergence of mitonuclear genes is increased relative to control nuclear genes between *M. arctoides* and the species sharing its nuclear ancestry (P = 0.007), consistent with genetic conflict. However, measures of introgression and selection do not identify large-scale co-introgression or co-evolution as means to resolve mitonuclear conflict. Nonetheless, mitochondrial tRNA-synthetases stand out in analyses using dN/dS and extended branch lengths as potential targets of selection. The methodology implemented provides a framework that can be used to examine the effects of mitonuclear co-introgression and co-evolution on a genomic scale in a variety of systems.

# Data Analysis

The directory ["Gene-Processing-Scripts"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/tree/master/Gene-Processing-Scripts)
contains scripts used to extract coding transcripts form genome sequences and to obtain "matching" data for N-mt and nuclear genes. Dependendencies and authorship of source scripts are stated within scripts.

A single data file ["Gene_Data.csv"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/blob/master/Gene_Data.csv) is provided that contains all information gathered from scripts in the above directories. A single R script ["Mitonuclear\_Stats\_Plots.R"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/blob/master/Mitonuclear_Stats_Plots.R) is provided that will analyze this data file and produce Figures and text files with statistical analyses and output these respectively to ["Figures"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/tree/master/Figures) and ["Stats"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/tree/master/Stats). When this git repo is cloned it should be possible to simply run "Rscript Mitonuclear_Stats_Plots.R" and figures and stat files should be output to their respective directories.

# Usage of Gene Processing Scripts

## Matching N-mt genes to control nuclear genes

N-mt genes were matched to control nuclear genes if they were (1) present on the same chromosome and (2) had similar lengths, (3) recombination rates, and (4) GC content. Therefore, this info should be obtained for the N-mt genes of interest and a large database of other nuclear genes (in this case, all Ensembl annotated genes not included in this study). The annotation used here is the rheMac8 ensGene gtf: https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/genes/rheMac8.ensGene.gtf.gz

Chromosome information and transcript coordinates are given in the Ensembl annotation. Therefore, a set of coordinates for each gene, inclusive of all transcripts for the gene, can be given by using the following script: https://github.com/StevisonLab/Arctoides-Hybridization/blob/main/scripts/gene_list_filter.pl

gene_matcher.py will use the start and end positions of each gene in the resulting bed file to calculate the length of each gene, though they could also be obtained from the resulting bed file with an awk command like so:

```
awk '{diff=$3-$2 ; print diff}' genes.bed
```

To obtain the GC content per gene using the masked rheMac8 FASTA (https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/rheMac8.fa.masked.gz) as a reference, bedtools 2.30.0 was used to obtain GC content using a command like so:

```
bedtools nuc -fi rheMac8.masked.fa -bed genes.bed
```

Recombination rates can be obtained from a bigwig file such as: ftp://ftp.hgsc.bcm.edu/ucscHub/rhesusSNVs/rheMac8/all.rate.bw

Then this can be lifted over to a bed file using UCSC liftOver. Given a bed file with gene coordinates and a bed file with recombination rates a new output file can be created using bedtools map with the gene coordinates and averaged recombination rate.

```
bedtools map -a genes.bed -b RecomRate_CM.bed -c 5 -o mean > output.bed
```
Once all these parameters are found one can create a bed file with with the following tab-separated values by combining all the data:

```
chromosome	start	end	gene	recom_rate	GC_content
chr1    15716094        15749833        SDHB    0.292168        44.257387593
...
```

The header is included here so you can know what the parameters are and the order they are in. gene_matcher.py assumes there is no header so it would be best to make the file without the header as well. Once files of this sort are made for both the query genes (N-mt genes in our case) and a large set of potentially matching genes (subject or other nuclear genes in our case) then one can use the query file as first argument, subject file as second argument, and desired final output as third argument.

```
python3 gene_matcher.py Query.bed Subject.bed Output.bed
```

The gene_matcher.py script uses the python random module to pick a random matching gene if multiple match all the criteria of the query gene. This necessarily adds randomness to the study.

## Extracting coding sequences for dN/dS analysis

The shell script ExtractCodingSequences.sh does most of the heavy lifting here. Details are explained in comments in the script but briefly it requires a VCF file with all samples of interest and a reference FASTA genome. It also requires several third-party programs to run (module loads in shell script) and requires the ReplaceFastaNames R and Python scripts to be present in the same directory. Other than what is in the script, it is worth noting that the VCF file should contain both invariant sites and SNP's, but not indels. We used bcftools view to retain only "ref" and "snp" types (http://samtools.github.io/bcftools/bcftools.html#expressions). The invariant sites and SNP's are both needed to create an accurate variant consensus FASTA sequence. The shell script will mask missing sites but not invariant, so the former should be explicitly given in the VCF file. Indels should be removed because they will throw off alignment given by gffread in the shell script.

# Population Genetic Analysis

Then methodology used to obtained topology weights, fdM values and Dxy values are described in the following repository: https://github.com/StevisonLab/Arctoides-Hybridization and largely implements scripts found in https://github.com/simonhmartin/genomics_general and https://github.com/simonhmartin/twisst

## dN/dS Analyses

Scripts used to implement PAML for dN/dS analyses and to properly weight results of these analyses for a per-gene estimate are provided in ["dNdS-Analysis-Scripts"](https://github.com/StevisonLab/Mitonuclear-Analysis-Project/tree/master/dNdS-Analysis-Scripts).

Usage of these scripts is largely explained in the scripts themselves. For a brief explanation, codeml.sh should be run on gene seqeunce files coming from ExtractCodingSequences.sh and dNdS_extraction.sh should be used after codeml is completed. File locations are hard-coded into the scripts so may need to be modified for future usage.

The Macaque_Trees.nwk file is the file used for the codeml.ctl file for all transcripts. The tree includes all possible topologies of the relationship between Arc, Sin and Fas with labels for each of these groups (described in more detail in the manuscript methods section). 
