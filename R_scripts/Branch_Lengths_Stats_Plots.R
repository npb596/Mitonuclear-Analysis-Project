# Load module necessary for generating contingency table image
library(gridExtra)

# Set working directory to Selection Data
setwd("~/Desktop/Macaque_files/mitonuclear_project/pop_gen_analysis/Selection_Data")

# Replace file as needed to enter contingency data
selection=read.table("functional_category_contingency.tsv",header=T,stringsAsFactors = T)

# Create contingency table
selection_table=table(selection$Category, selection$Extension)

# Crate nice simple image of contingency table
grid.table(selection_table)

# Conduct Fisher's exact test on table, change alternative as needed
fisher.test(selection_table, alternative = "greater")
