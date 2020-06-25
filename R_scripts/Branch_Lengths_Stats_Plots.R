library(gridExtra)

setwd("~/Desktop/Macaque_files/Altitude_Project")


selection=read.table("contingency_table_mitonuclear.tsv",header=T,stringsAsFactors = T)

selection_table=table(selection$Source, selection$Extension)

table(selection$Altitude, selection$Base)

#pdf("selection_contingency.pdf")
grid.table(selection_table)
#dev.off()

#sink("selection_fisher_test.txt")
fisher.test(selection_table, alternative = "less")
#sink()
