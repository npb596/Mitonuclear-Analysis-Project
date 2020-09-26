# Load modules for plots and summary statistics
library(ggplot2)
library(ggpubr)
library(plyr)
library(gridExtra)

# Set working directory to Introgression Data
setwd("~/Desktop/Macaque_files/mitonuclear_project/pop_gen_analysis/Introgression_Data")

# Define tables from introgression data files
ALL.introgression=read.table("fdM_CDS.txt",header=T,stringsAsFactors = T)
ETS.introgression=read.table("ETS_fdM.txt",header=T,stringsAsFactors = T)
MRP.introgression=read.table("ribosomal_fdM.txt",header=T,stringsAsFactors = T)
MTS.introgression=read.table("synthetase_fdM.txt",header=T,stringsAsFactors = T)
MISC.introgression=read.table("miscellaneous_fdM.txt",header=T,stringsAsFactors = T)

# Define variables containing means of data tables
mod_ALL.fdm <- ddply(ALL.introgression, "ALL.DNAType", summarise, grp.mean=mean(ALL.fdm))
mod_ETS.fdm <- ddply(ETS.introgression, "ETS.DNAType", summarise, grp.mean=mean(ETS.fdm))
mod_MRP.fdm <- ddply(MRP.introgression, "MRP.DNAType", summarise, grp.mean=mean(MRP.fdm))
mod_MTS.fdm <- ddply(MTS.introgression, "MTS.DNAType", summarise, grp.mean=mean(MTS.fdm))
mod_MISC.fdm <- ddply(MISC.introgression, "MISC.DNAType", summarise, grp.mean=mean(MISC.fdm))

# Plot all data with titles, x and y labels, coloring lines in density plots
ALL_fdm_plot <- ggplot(ALL.introgression, aes(x=ALL.fdm, color=ALL.DNAType)) +
  geom_text(x=0.15, y=2.2, label="N = 100, P = 0.35", size=10, show.legend = FALSE, colour="black") +
  geom_density() + labs(title="All Functional Categories", x="fdM", y = "Density") + theme(legend.position = "none") +
  theme(legend.title=element_text(size=50)) + theme(legend.text=element_text(size=45)) +
  scale_color_manual(name = "Gene Type", values=c("#af8dc3","#7fbf7b")) + 
  geom_vline(data=mod_ALL.fdm, aes(xintercept=grp.mean, color=ALL.DNAType), linetype="dashed") + 
  theme(text = element_text(size=29))
ETS_fdm_plot <- ggplot(ETS.introgression, aes(x=ETS.fdm, color=ETS.DNAType)) +
  geom_density() + labs(title="Electron Transport System", x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(values=c("#af8dc3","#7fbf7b")) + 
  geom_vline(data=mod_ETS.fdm, aes(xintercept=grp.mean, color=ETS.DNAType), linetype="dashed") +
  theme(text = element_text(size=19)) + 
  geom_text(x=-0.5, y=0.0, label="N = 48, P = 0.382", size=5, colour="black")
MRP_fdm_plot <- ggplot(MRP.introgression, aes(x=MRP.fdm, color=MRP.DNAType)) +
  geom_density() + labs(title="Mito Ribosomal Protein", x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(values=c("#af8dc3","#7fbf7b")) + 
  geom_vline(data=mod_MRP.fdm, aes(xintercept=grp.mean, color=MRP.DNAType), linetype="dashed") +
  theme(text = element_text(size=19)) +
  geom_text(x=0.15, y=2.2, label="N = 29, P = 0.354", size=5, colour="black")
MTS_fdm_plot <- ggplot(MTS.introgression, aes(x=MTS.fdm, color=MTS.DNAType)) +
  geom_density() + labs(title="Mito tRNA-synthetase", x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(values=c("#af8dc3","#7fbf7b")) + 
  geom_vline(data=mod_MTS.fdm, aes(xintercept=grp.mean, color=MTS.DNAType), linetype="dashed") + 
  theme(text = element_text(size=19)) +
  geom_text(x=0.025, y=3.27, label="N = 10, P = 0.685", size=5, colour="black")
MISC_fdm_plot <- ggplot(MISC.introgression, aes(x=MISC.fdm, color=MISC.DNAType)) +
  geom_density() + labs(title="Transcription/translation", x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(values=c("#af8dc3","#7fbf7b")) + 
  geom_vline(data=mod_MISC.fdm, aes(xintercept=grp.mean, color=MISC.DNAType), linetype="dashed") + 
  theme(text = element_text(size=19)) +
  geom_text(x=-0.5, y=0.0, label="N = 13, P = 0.35", size=5, colour="black")

# Create legend using info from ALL_fdm_plot above
ALL_fdm_legend <- get_legend(ALL_fdm_plot, position = "bottom")

# Combine all plots defined above and add labels
ggarrange(ALL_fdm_plot,
          ggarrange(ETS_fdm_plot,MRP_fdm_plot,MTS_fdm_plot,MISC_fdm_plot,
                    labels = c("B", "C", "D","E"), label.x = 0.13, label.y = 0.91, font.label = list(size = 27),
                    ncol = 2, nrow = 2), labels="A", label.x = 0.13, label.y = 0.93, font.label = list(size = 37), 
          legend.grob = ALL_fdm_legend)

#Change variables according to input file and conduct wilcoxon rank-sum test
wilcox.test(ALL.fdm ~ ALL.DNAType, data=ALL.introgression, alternative="greater")
