library(ggplot2)
library(ggpubr)
library(plyr)
library(gridExtra)

setwd("~/Desktop/Macaque_files/mitonuclear_project/pop_gen_analysis")

ALL.introgression=read.table("Introgression_Data/fdM_CDS.txt",header=F,stringsAsFactors = T)
ETS.introgression=read.table("Introgression_Data/ETS_fdM.txt",header=F,stringsAsFactors = T)
MRP.introgression=read.table("Introgression_Data/ribosomal_fdM.txt",header=F,stringsAsFactors = T)
MTS.introgression=read.table("Introgression_Data/synthetase_fdM.txt",header=F,stringsAsFactors = T)
MISC.introgression=read.table("Introgression_Data/miscellaneous_fdM.txt",header=F,stringsAsFactors = T)

colnames(ALL.introgression)<-c("ALL.fdm","ALL.DNAType")
colnames(ETS.introgression)<-c("ETS.fdm","Gene","ETS.DNAType")
colnames(MRP.introgression)<-c("MRP.fdm","Gene","MRP.DNAType")
colnames(MTS.introgression)<-c("MTS.fdm","Gene","MTS.DNAType")
colnames(MISC.introgression)<-c("MISC.fdm","Gene","MISC.DNAType")

mod_ALL.fdm <- ddply(ALL.introgression, "ALL.DNAType", summarise, grp.mean=mean(ALL.fdm))
mod_ETS.fdm <- ddply(ETS.introgression, "ETS.DNAType", summarise, grp.mean=mean(ETS.fdm))
mod_MRP.fdm <- ddply(MRP.introgression, "MRP.DNAType", summarise, grp.mean=mean(MRP.fdm))
mod_MTS.fdm <- ddply(MTS.introgression, "MTS.DNAType", summarise, grp.mean=mean(MTS.fdm))
mod_MISC.fdm <- ddply(MISC.introgression, "MISC.DNAType", summarise, grp.mean=mean(MISC.fdm))

ALL_fdm_plot <- ggplot(ALL.introgression, aes(x=ALL.fdm, color=ALL.DNAType)) +
  geom_density() + labs(x="fdM", y = "Density") + theme(legend.position = "none") +
  theme(legend.title=element_text(size=50)) + theme(legend.text=element_text(size=45)) +
  scale_color_manual(name = "Gene Type", labels = c("N-mt", "Nuclear"), values=c("darkorchid3","dodgerblue1")) + 
  geom_vline(data=mod_ALL.fdm, aes(xintercept=grp.mean, color=ALL.DNAType), linetype="dashed") + 
  theme(text = element_text(size=29))
ETS_fdm_plot <- ggplot(ETS.introgression, aes(x=ETS.fdm, color=ETS.DNAType)) +
  geom_density() + labs(x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(name = "Gene Type", labels = c("N-mt", "Nuclear"), values=c("darkorchid3","dodgerblue1")) + 
  geom_vline(data=mod_ETS.fdm, aes(xintercept=grp.mean, color=ETS.DNAType), linetype="dashed") +
  theme(text = element_text(size=19))
MRP_fdm_plot <- ggplot(MRP.introgression, aes(x=MRP.fdm, color=MRP.DNAType)) +
  geom_density() + labs(x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(name = "Gene Type", labels = c("N-mt", "Nuclear"), values=c("darkorchid3","dodgerblue1")) + 
  geom_vline(data=mod_MRP.fdm, aes(xintercept=grp.mean, color=MRP.DNAType), linetype="dashed") +
  theme(text = element_text(size=19))
MTS_fdm_plot <- ggplot(MTS.introgression, aes(x=MTS.fdm, color=MTS.DNAType)) +
  geom_density() + labs(x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(name = "Gene Type", labels = c("N-mt", "Nuclear"), values=c("darkorchid3","dodgerblue1")) + 
  geom_vline(data=mod_MTS.fdm, aes(xintercept=grp.mean, color=MTS.DNAType), linetype="dashed") + 
  theme(text = element_text(size=19))
MISC_fdm_plot <- ggplot(MISC.introgression, aes(x=MISC.fdm, color=MISC.DNAType)) +
  geom_density() + labs(x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(name = "Gene Type", labels = c("N-mt", "Nuclear"), values=c("darkorchid3","dodgerblue1")) + 
  geom_vline(data=mod_MISC.fdm, aes(xintercept=grp.mean, color=MISC.DNAType), linetype="dashed") + 
  theme(text = element_text(size=19))

ALL_fdm_legend <- get_legend(ALL_fdm_plot, position = "bottom")

ggarrange(ALL_fdm_plot,
          ggarrange(ETS_fdm_plot,MRP_fdm_plot,MTS_fdm_plot,MISC_fdm_plot,
                    labels = c("B", "C", "D","E"), label.x = 0.13, font.label = list(size = 27),
                    ncol = 2, nrow = 2), labels="A", label.x = 0.13, font.label = list(size = 37), 
          legend.grob = ALL_fdm_legend)

#sink("miscellaneous_fdM_wilcoxon.txt")
wilcox.test(ALL.fdm ~ ALL.DNAType, data=ALL.introgression, alternative="greater")
#sink()