library(ggplot2)
library(ggpubr)
library(plyr)
library(gridExtra)

setwd("~/Desktop/Macaque_files/mitonuclear_project/pop_gen_analysis")

ALL.dxy=read.table("Divergence_Data/Dxy_Sin_Arc_CDS.txt",header=F,stringsAsFactors = T)
ETS.dxy=read.table("Divergence_Data/ETS_Sin_Arc_Dxy.txt",header=F,stringsAsFactors = T)
MRP.dxy=read.table("Divergence_Data/ribosomal_Sin_Arc_Dxy.txt",header=F,stringsAsFactors = T)
MTS.dxy=read.table("Divergence_Data/synthetase_Sin_Arc_Dxy.txt",header=F,stringsAsFactors = T)
MISC.dxy=read.table("Divergence_Data/miscellaneous_Sin_Arc_Dxy.txt",header=F,stringsAsFactors = T)
colnames(ALL.dxy)<-c("ALL.dxy","ALL.DNAType")
colnames(ETS.dxy)<-c("ETS.dxy","Gene","ETS.DNAType")
colnames(MRP.dxy)<-c("MRP.dxy","Gene","MRP.DNAType")
colnames(MTS.dxy)<-c("MTS.dxy","Gene","MTS.DNAType")
colnames(MISC.dxy)<-c("MISC.dxy","Gene","MISC.DNAType")

mod_ALL.dxy <- ddply(ALL.dxy, "ALL.DNAType", summarise, grp.mean=mean(ALL.dxy))
mod_ETS.dxy <- ddply(ETS.dxy, "ETS.DNAType", summarise, grp.mean=mean(ETS.dxy))
mod_MRP.dxy <- ddply(MRP.dxy, "MRP.DNAType", summarise, grp.mean=mean(MRP.dxy))
mod_MTS.dxy <- ddply(MTS.dxy, "MTS.DNAType", summarise, grp.mean=mean(MTS.dxy))
mod_MISC.dxy <- ddply(MISC.dxy, "MISC.DNAType", summarise, grp.mean=mean(MISC.dxy))

ALL_plot <- ggplot(ALL.dxy, aes(x=ALL.dxy, color=ALL.DNAType)) +
  geom_density() + labs(x="Dxy", y = "Density") + theme(legend.position = "none") +
  theme(legend.title=element_text(size=50)) + theme(legend.text=element_text(size=45)) +
  scale_color_manual(name = "Gene Type", labels = c("N-mt", "Nuclear"), values=c("darkorchid3","dodgerblue1")) + 
  geom_vline(data=mod_ALL.dxy, aes(xintercept=grp.mean, color=ALL.DNAType), linetype="dashed") + 
  theme(text = element_text(size=29))
ETS_plot <- ggplot(ETS.dxy, aes(x=ETS.dxy, color=ETS.DNAType)) +
  geom_density() + labs(x="Dxy", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(name = "Gene Type", labels = c("N-mt", "Nuclear"), values=c("darkorchid3","dodgerblue1")) + 
  geom_vline(data=mod_ETS.dxy, aes(xintercept=grp.mean, color=ETS.DNAType), linetype="dashed") +
  theme(text = element_text(size=19))
MRP_plot <- ggplot(MRP.dxy, aes(x=MRP.dxy, color=MRP.DNAType)) +
  geom_density() + labs(x="Dxy", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(name = "Gene Type", labels = c("N-mt", "Nuclear"), values=c("darkorchid3","dodgerblue1")) + 
  geom_vline(data=mod_MRP.dxy, aes(xintercept=grp.mean, color=MRP.DNAType), linetype="dashed") +
  theme(text = element_text(size=19))
MTS_plot <- ggplot(MTS.dxy, aes(x=MTS.dxy, color=MTS.DNAType)) +
  geom_density() + labs(x="Dxy", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(name = "Gene Type", labels = c("N-mt", "Nuclear"), values=c("darkorchid3","dodgerblue1")) + 
  geom_vline(data=mod_MTS.dxy, aes(xintercept=grp.mean, color=MTS.DNAType), linetype="dashed") + 
  theme(text = element_text(size=19))
MISC_plot <- ggplot(MISC.dxy, aes(x=MISC.dxy, color=MISC.DNAType)) +
  geom_density() + labs(x="Dxy", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(name = "Gene Type", labels = c("N-mt", "Nuclear"), values=c("darkorchid3","dodgerblue1")) + 
  geom_vline(data=mod_MISC.dxy, aes(xintercept=grp.mean, color=MISC.DNAType), linetype="dashed") + 
  theme(text = element_text(size=19))

ALL_legend <- get_legend(ALL_plot, position = "bottom")

ggarrange(ALL_plot,
          ggarrange(ETS_plot,MRP_plot,MTS_plot,MISC_plot,
                    labels = c("B", "C", "D","E"), label.x = 0.12, font.label = list(size = 27),
                    ncol = 2, nrow = 2), labels="A", label.x = 0.12, font.label = list(size = 37), 
          legend.grob = ALL_legend)
