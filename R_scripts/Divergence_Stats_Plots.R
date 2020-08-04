# Load modules for plots
library(ggplot2)
library(ggpubr)
library(plyr)

# Set working directory to Divergence Data
setwd("~/Desktop/Macaque_files/mitonuclear_project/pop_gen_analysis/Divergence_Data")

# Define tables from divergence data files
ALL.dxy=read.table("ALL_Combined_Dxy.txt",header=T,stringsAsFactors = T)
ETS.dxy=read.table("ETS_Combined_Dxy.txt",header=T,stringsAsFactors = T)
MRP.dxy=read.table("ribosomal_Combined_Dxy.txt",header=T,stringsAsFactors = T)
MTS.dxy=read.table("synthetase_Combined_Dxy.txt",header=T,stringsAsFactors = T)
MISC.dxy=read.table("miscellaneous_Combined_Dxy.txt",header=T,stringsAsFactors = T)

# Plot all data with titles, x and y labels, colors filling violin plots
ALL_plot <- ggplot(ALL.dxy, aes(x=ALL.Comparison, y=ALL.dxy, fill = ALL.DNAType)) +
  labs(title="All Functional Categories", x ="Species", y="Dxy") +
  geom_violin(trim=FALSE) + theme(legend.position = "none") +
  theme(legend.title=element_text(size=50)) + theme(legend.text=element_text(size=45)) +
  scale_fill_manual(name = "Gene Type", values=c("#af8dc3","#7fbf7b")) + 
  theme(text = element_text(size=29)) + geom_text(x=2.1, y=0.65, label="N = 101, P = 0.009", size=10) + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9))
ETS_plot <- ggplot(ETS.dxy, aes(x=ETS.Comparison, y=ETS.dxy, fill = ETS.DNAType)) +
  labs(title="Electron Transport System", x ="Species", y="Dxy") + 
  geom_violin(trim=FALSE) + theme(legend.position = "none") +
  scale_fill_manual(values=c("#af8dc3","#7fbf7b")) +
  theme(text = element_text(size=19)) + geom_text(x=2.1, y=0.65, label="N = 49, P = 0.047", size=5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9))
MRP_plot <- ggplot(MRP.dxy, aes(x=MRP.Comparison, y=MRP.dxy, fill = MRP.DNAType)) +
  labs(title="Mito Ribosomal Protein", x ="Species", y="Dxy") + 
  geom_violin(trim=FALSE) + theme(legend.position = "none") +
  scale_fill_manual(name = "Gene Type", values=c("#af8dc3","#7fbf7b")) +
  theme(text = element_text(size=19)) + geom_text(x=2.1, y=0.6, label="N = 29, P = 0.036", size=5) +
  geom_boxplot(width = 0.15, position = position_dodge(0.9))
MTS_plot <- ggplot(MTS.dxy, aes(x=MTS.Comparison, y=MTS.dxy, fill = MTS.DNAType)) +
  labs(title="Mito tRNA-synthetase", x ="Species", y="Dxy") + 
  geom_violin(trim=FALSE) + theme(legend.position = "none") +
  scale_fill_manual(name = "Gene Type", values=c("#af8dc3","#7fbf7b")) +
  theme(text = element_text(size=19)) + geom_text(x=2.1, y=0.4, label="N = 10, P = 0.804", size=5) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9))
MISC_plot <- ggplot(MISC.dxy, aes(x=MISC.Comparison, y=MISC.dxy, fill = MISC.DNAType)) +
  labs(title="Transcription/translation", x ="Species", y="Dxy") + 
  geom_violin(trim=FALSE) + theme(legend.position = "none") +
  scale_fill_manual(name = "Gene Type", values=c("#af8dc3","#7fbf7b")) +
  theme(text = element_text(size=19)) + geom_text(x=2.1, y=0.45, label="N = 13, P = 0.051", size=5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9))

# Create legend using info from ALL_plot above
ALL_legend <- get_legend(ALL_plot, position = "bottom")

# Combine all plots defined above and add labels
ggarrange(ALL_plot,
          ggarrange(ETS_plot,MRP_plot,MTS_plot,MISC_plot,
                    labels = c("B", "C", "D","E"), label.x = 0.17, label.y = 0.91, font.label = list(size = 27),
                    ncol = 2, nrow = 2), labels="A", label.x = 0.12, label.y = 0.93, font.label = list(size = 37), 
          legend.grob = ALL_legend)

# Input file must be redefined to one including only Arc-Sin or Arc-Fas data
# Then change variables and alternative accordingly, then conduct wilcoxon rank-sum test
wilcox.test(ALL.dxy ~ ALL.DNAType, data=ALL.dxy, alternative="greater")
