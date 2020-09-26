library(ggplot2)
library(plyr)

setwd("~/Desktop/Macaque_files/mitonuclear_project/pop_gen_analysis")

selection=read.table("Selection_Data/Master_dNdS.txt",header=F,stringsAsFactors = T)

colnames(selection)<-c("dnds","Pairing")

res.aov <- aov(dnds ~ Pairing, data = selection)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)

#pdf("Selection_results/dNdS_boxplot.pdf")
boxplot(dnds~Pairing, data=selection, 
        col=c("mediumpurple1","lightcoral","deepskyblue3"),
        ylab="dN/dS", xlab="",
        #names=c('Fas-Arc', 'Fas-Arc','Sin-Arc','Sin-Arc','Sin-Fas','Sin-Fas'),
        cex.main=1.8, cex.axis=1.15, cex.lab=1.25)
#dev.off()
