# This script includes analysis for the manuscript titled:
# Examination of Mitonuclear Conflict in a Macaque Species Exhibiting Mitonuclear Discordance

# It is split into 4 code sections which can be expanded/collapsed in RStudio 
# Expand: either click on the arrow in the gutter or on the icon that overlays the folded code 
# Collapse: click on the arrow in the gutter

# This script can be run on a command line by simply running "Rscript Mitonuclear_Stats_Plots.R"
# in the directory containing the script. 
# If running on RStudio each command can be run individually
# Once all commands are run the script will output figures and statistical results
# to the Figures and Stats directories respectively
# Currently, sourcing the script in RStudio will not produce some figures and results

# The following libraries are necessary for the below analyses
# If running on RStudio then running each line will prompt the user to install the library
library(ggplot2)
library(ggpubr)
library(plyr)
library(gridExtra)
library(dplyr)

# Set working directory according to location of data file if using through RStudio
# Working directory should be set to directory containing the Git Repo
# As of now it assumes the home directory contains the git repo
#setwd("~/Mitonuclear-Analysis-Project")

# Read in complete data file for this study
# Also read in a data file containing only paired N-mt and nuclear genes
Gene_Data=read.csv("Gene_Data.csv",header=T,stringsAsFactors = T)
Paired_Gene_Data <- subset(Gene_Data, (!is.na(Gene_Data[,19])))

# Divergence ---------------------------------------------------------------

# Define a data frame with Dxy statistics between Arctoides and Sinica from paired data
# Then assign names to columns to allow for consistency with Arc-Fas statistics
# Then exclude genes with no Dxy value and genes paired to these
Sin_Dxy_one <- data.frame(Paired_Gene_Data$Sin.Dxy,Paired_Gene_Data$Type,
                      Paired_Gene_Data$FuncCategory,
                      Paired_Gene_Data$Pairing,Paired_Gene_Data$Gene,"Sinica")
names(Sin_Dxy_one)<-c("dxy","Type","Category","Pairing","Gene","Comparison")
Sin_Dxy_two <- Sin_Dxy_one[complete.cases(Sin_Dxy_one), ]
Sin_Dxy <- subset(Sin_Dxy_two,duplicated(Sin_Dxy_two$Pairing) | 
                    duplicated(Sin_Dxy_two$Pairing, fromLast=TRUE))

# Perform the steps described above but with Arctoides and Fascicularis data
Fas_Dxy_one <- data.frame(Paired_Gene_Data$Fas.Dxy,Paired_Gene_Data$Type,
                      Paired_Gene_Data$FuncCategory,
                      Paired_Gene_Data$Pairing,Paired_Gene_Data$Gene,"Fascicularis")
names(Fas_Dxy_one)<-c("dxy","Type","Category","Pairing","Gene","Comparison")
Fas_Dxy_two <- Fas_Dxy_one[complete.cases(Fas_Dxy_one), ]
Fas_Dxy <- subset(Fas_Dxy_two,duplicated(Fas_Dxy_two$Pairing) | 
                    duplicated(Fas_Dxy_two$Pairing, fromLast=TRUE))

# Combine the above datasets and subset all functional categories
Dxy <- rbind(Sin_Dxy, Fas_Dxy)
ETS=Dxy[Dxy$Category == "ETS", ]
MRP=Dxy[Dxy$Category == "MRP", ]
MTS=Dxy[Dxy$Category == "MTS", ]
MISC=Dxy[Dxy$Category == "MISC", ]
Sin_ETS=Sin_Dxy[Sin_Dxy$Category == "ETS", ]
Sin_MRP=Sin_Dxy[Sin_Dxy$Category == "MRP", ]
Sin_MTS=Sin_Dxy[Sin_Dxy$Category == "MTS", ]
Sin_MISC=Sin_Dxy[Sin_Dxy$Category == "MISC", ]
Fas_ETS=Fas_Dxy[Fas_Dxy$Category == "ETS", ]
Fas_MRP=Fas_Dxy[Fas_Dxy$Category == "MRP", ]
Fas_MTS=Fas_Dxy[Fas_Dxy$Category == "MTS", ]
Fas_MISC=Fas_Dxy[Fas_Dxy$Category == "MISC", ]

# Plot all data with titles, x (Sinica or Fascicularis) and y (Dxy values) labels, 
# and colors (N-mt and Nuclear) filling in violin plots 
ALL_plot <- ggplot(Dxy, aes(x=Comparison, y=dxy, fill = Type)) +
  labs(title="All Functional Categories", x ="Species", y="Dxy") +
  geom_violin(trim=FALSE) + theme(legend.position = "none") +
  theme(legend.title=element_text(size=25)) + theme(legend.text=element_text(size=23)) +
  scale_fill_manual(name = "Gene Type", values=c("#af8dc3","#7fbf7b")) + 
  theme(text = element_text(size=15)) + geom_text(x=2.1, y=0.65, label="N = 101, P = 0.009", size=5) + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9))
ETS_plot <- ggplot(ETS, aes(x=Comparison, y=dxy, fill = Type)) +
  labs(title="Electron Transport System", x ="Species", y="Dxy") + 
  geom_violin(trim=FALSE) + theme(legend.position = "none") +
  scale_fill_manual(values=c("#af8dc3","#7fbf7b")) +
  theme(text = element_text(size=10)) + geom_text(x=2.1, y=0.65, label="N = 49, P = 0.047", size=2) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9))
MRP_plot <- ggplot(MRP, aes(x=Comparison, y=dxy, fill = Type)) +
  labs(title="Mito Ribosomal Protein", x ="Species", y="Dxy") + 
  geom_violin(trim=FALSE) + theme(legend.position = "none") +
  scale_fill_manual(name = "Gene Type", values=c("#af8dc3","#7fbf7b")) +
  theme(text = element_text(size=10)) + geom_text(x=2.1, y=0.6, label="N = 29, P = 0.036", size=2) +
  geom_boxplot(width = 0.15, position = position_dodge(0.9))
MTS_plot <- ggplot(MTS, aes(x=Comparison, y=dxy, fill = Type)) +
  labs(title="Mito tRNA-synthetase", x ="Species", y="Dxy") + 
  geom_violin(trim=FALSE) + theme(legend.position = "none") +
  scale_fill_manual(name = "Gene Type", values=c("#af8dc3","#7fbf7b")) +
  theme(text = element_text(size=10)) + geom_text(x=2.1, y=0.4, label="N = 10, P = 0.804", size=2) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9))
MISC_plot <- ggplot(MISC, aes(x=Comparison, y=dxy, fill = Type)) +
  labs(title="Transcription/translation", x ="Species", y="Dxy") + 
  geom_violin(trim=FALSE) + theme(legend.position = "none") +
  scale_fill_manual(name = "Gene Type", values=c("#af8dc3","#7fbf7b")) +
  theme(text = element_text(size=10)) + geom_text(x=2.1, y=0.45, label="N = 13, P = 0.051", size=2) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9))

# Create legend using info from ALL_plot above
ALL_legend <- get_legend(ALL_plot, position = "bottom")

# Combine all plots defined above and add labels
png("Figures/Figure2.png",res=300,height = 5,width = 10,units="in", pointsize=1)
ggarrange(ALL_plot,
          ggarrange(ETS_plot,MRP_plot,MTS_plot,MISC_plot,
                    labels = c("B", "C", "D","E"), font.label = list(size = 14),
                    ncol = 2, nrow = 2), labels="A", font.label = list(size = 19), 
          legend.grob = ALL_legend)
dev.off()

# Conduct Wilcoxon Rank-Sum tests on all genes as well as category subsets
# and print these out as a single text file
# Expectation for Sinica Dxy is N-mt greater than nuclear
# Expectation for Fascicularis Dxy is the opposite
sink("Stats/Dxy_Wilcoxon_results.txt")
print("Dxy between Sinica and Arctoides for all paired genes")
wilcox.test(dxy ~ Type, data=Sin_Dxy, alternative="greater")
print("Dxy between Sinica and Arctoides for paired ETS genes")
wilcox.test(dxy ~ Type, data=Sin_ETS, alternative="greater")
print("Dxy between Sinica and Arctoides for paired MRP genes")
wilcox.test(dxy ~ Type, data=Sin_MRP, alternative="greater")
print("Dxy between Sinica and Arctoides for paired MTS genes")
wilcox.test(dxy ~ Type, data=Sin_MTS, alternative="greater")
print("Dxy between Sinica and Arctoides for paired MISC genes")
wilcox.test(dxy ~ Type, data=Sin_MISC, alternative="greater")
print("Dxy between Fascicularis and Arctoides for all paired genes")
wilcox.test(dxy ~ Type, data=Fas_Dxy, alternative="less")
print("Dxy between Fascicularis and Arctoides for paired ETS genes")
wilcox.test(dxy ~ Type, data=Fas_ETS, alternative="less")
print("Dxy between Fascicularis and Arctoides for paired MRP genes")
wilcox.test(dxy ~ Type, data=Fas_MRP, alternative="less")
print("Dxy between Fascicularis and Arctoides for paired MTS genes")
wilcox.test(dxy ~ Type, data=Fas_MTS, alternative="less")
print("Dxy between Fascicularis and Arctoides for paired MISC genes")
wilcox.test(dxy ~ Type, data=Fas_MISC, alternative="less")
sink()

# Introgression -----------------------------------------------------------

# Define a data frame with fdM statistics from paired data
# Then exclude genes with no fdM value and genes paired to these
Fdm_one <- data.frame(Paired_Gene_Data$fdM,Paired_Gene_Data$Type,
                          Paired_Gene_Data$FuncCategory,
                          Paired_Gene_Data$Pairing,Paired_Gene_Data$Gene)
Fdm_two <- Fdm_one[complete.cases(Fdm_one), ]
Fdm <- subset(Fdm_two,duplicated(Fdm_two$Paired_Gene_Data.Pairing) | 
                    duplicated(Fdm_two$Paired_Gene_Data.Pairing, fromLast=TRUE))

# Subset all functional categories and obtain means for N-mt and nuclear genes
Fdm_ETS=Fdm[Fdm$Paired_Gene_Data.FuncCategory == "ETS", ]
Fdm_MRP=Fdm[Fdm$Paired_Gene_Data.FuncCategory == "MRP", ]
Fdm_MTS=Fdm[Fdm$Paired_Gene_Data.FuncCategory == "MTS", ]
Fdm_MISC=Fdm[Fdm$Paired_Gene_Data.FuncCategory == "MISC", ]

mod_ALL.fdm <- ddply(Fdm, "Paired_Gene_Data.Type", summarise, grp.mean=mean(Paired_Gene_Data.fdM))
mod_ETS.fdm <- ddply(Fdm_ETS, "Paired_Gene_Data.Type", summarise, grp.mean=mean(Paired_Gene_Data.fdM))
mod_MRP.fdm <- ddply(Fdm_MRP, "Paired_Gene_Data.Type", summarise, grp.mean=mean(Paired_Gene_Data.fdM))
mod_MTS.fdm <- ddply(Fdm_MTS, "Paired_Gene_Data.Type", summarise, grp.mean=mean(Paired_Gene_Data.fdM))
mod_MISC.fdm <- ddply(Fdm_MISC, "Paired_Gene_Data.Type", summarise, grp.mean=mean(Paired_Gene_Data.fdM))

# Plot all data with titles, x (fdM) and y (density) labels, 
# and coloring lines (N-mt and nuclear) in density plots

ALL_fdm_plot <- ggplot(Fdm, aes(x=Paired_Gene_Data.fdM, color=Paired_Gene_Data.Type)) +
  geom_text(x=0.15, y=2.2, label="N = 100, P = 0.35", size=10, show.legend = FALSE, colour="black") +
  geom_density() + labs(title="All Functional Categories", x="fdM", y = "Density") + theme(legend.position = "none") +
  theme(legend.title=element_text(size=25)) + theme(legend.text=element_text(size=23)) +
  scale_color_manual(name = "Gene Type", values=c("#af8dc3","#7fbf7b")) + 
  geom_vline(data=mod_ALL.fdm, aes(xintercept=grp.mean, color=Paired_Gene_Data.Type), linetype="dashed") + 
  theme(text = element_text(size=15))
ETS_fdm_plot <- ggplot(Fdm_ETS, aes(x=Paired_Gene_Data.fdM, color=Paired_Gene_Data.Type)) +
  geom_density() + labs(title="Electron Transport System", x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(values=c("#af8dc3","#7fbf7b")) + 
  geom_vline(data=mod_ETS.fdm, aes(xintercept=grp.mean, color=Paired_Gene_Data.Type), linetype="dashed") +
  theme(text = element_text(size=10)) + 
  geom_text(x=-0.5, y=0.0, label="N = 48, P = 0.382", size=5, colour="black")
MRP_fdm_plot <- ggplot(Fdm_MRP, aes(x=Paired_Gene_Data.fdM, color=Paired_Gene_Data.Type)) +
  geom_density() + labs(title="Mito Ribosomal Protein", x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(values=c("#af8dc3","#7fbf7b")) + 
  geom_vline(data=mod_MRP.fdm, aes(xintercept=grp.mean, color=Paired_Gene_Data.Type), linetype="dashed") +
  theme(text = element_text(size=10)) +
  geom_text(x=0.15, y=2.2, label="N = 29, P = 0.354", size=5, colour="black")
MTS_fdm_plot <- ggplot(Fdm_MTS, aes(x=Paired_Gene_Data.fdM, color=Paired_Gene_Data.Type)) +
  geom_density() + labs(title="Mito tRNA-synthetase", x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(values=c("#af8dc3","#7fbf7b")) + 
  geom_vline(data=mod_MTS.fdm, aes(xintercept=grp.mean, color=Paired_Gene_Data.Type), linetype="dashed") + 
  theme(text = element_text(size=10)) +
  geom_text(x=0.025, y=3.27, label="N = 10, P = 0.685", size=5, colour="black")
MISC_fdm_plot <- ggplot(Fdm_MISC, aes(x=Paired_Gene_Data.fdM, color=Paired_Gene_Data.Type)) +
  geom_density() + labs(title="Transcription/translation", x="fdM", y = "Density") + theme(legend.position = "none") +
  scale_color_manual(values=c("#af8dc3","#7fbf7b")) + 
  geom_vline(data=mod_MISC.fdm, aes(xintercept=grp.mean, color=Paired_Gene_Data.Type), linetype="dashed") + 
  theme(text = element_text(size=10)) +
  geom_text(x=-0.5, y=0.0, label="N = 13, P = 0.35", size=5, colour="black")

# Create legend using info from ALL_fdm_plot above
ALL_fdm_legend <- get_legend(ALL_fdm_plot, position = "bottom")

# Combine all plots defined above and add labels
png("Figures/Figure3.png",res=300,height = 5,width = 10,units="in",pointsize = 1)
ggarrange(ALL_fdm_plot,
          ggarrange(ETS_fdm_plot,MRP_fdm_plot,MTS_fdm_plot,MISC_fdm_plot,
                    labels = c("B", "C", "D","E"), font.label = list(size = 14),
                    ncol = 2, nrow = 2), labels="A", font.label = list(size = 24), 
          legend.grob = ALL_fdm_legend)
dev.off()

# Conduct Wilcoxon rank-sum test on combined fdM data and all functional categories
# Expectation is N-mt introgression is more positive than nuclear
sink("Stats/fdM_Wilcoxon_results.txt")
print("fdM for all paired genes")
wilcox.test(Paired_Gene_Data.fdM ~ Paired_Gene_Data.Type, data=Fdm, alternative="greater")
print("fdM for ETS paired genes")
wilcox.test(Paired_Gene_Data.fdM ~ Paired_Gene_Data.Type, data=Fdm_ETS, alternative="greater")
print("fdM for MRP paired genes")
wilcox.test(Paired_Gene_Data.fdM ~ Paired_Gene_Data.Type, data=Fdm_MRP, alternative="greater")
print("fdM for MTS paired genes")
wilcox.test(Paired_Gene_Data.fdM ~ Paired_Gene_Data.Type, data=Fdm_MTS, alternative="greater")
print("fdM for MISC paired genes")
wilcox.test(Paired_Gene_Data.fdM ~ Paired_Gene_Data.Type, data=Fdm_MISC, alternative="greater")
sink()

# Selection ---------------------------------------------------------------

# Define a data frame with dNdS statistics between Arctoides and Silenus
# Then assign names to columns to allow for consistency with other dNdS statistics
# Then exclude genes with no fdM value and nuclear genes
Arc_dNdS_one <- data.frame(Gene_Data$Arc.dN.dS,Gene_Data$Type,
                      Gene_Data$FuncCategory,Gene_Data$Gene,"Arctoides")
names(Arc_dNdS_one)<-c("dNdS","Type","Category","Gene","Comparison")
Arc_dNdS_two <- Arc_dNdS_one[complete.cases(Arc_dNdS_one), ]
Arc_dNdS=Arc_dNdS_two[Arc_dNdS_two$Type == "N-mt", ]

# Perform the steps described above but with Sinica and Silenus data
Sin_dNdS_one <- data.frame(Gene_Data$Sin.dN.dS,Gene_Data$Type,
                       Gene_Data$FuncCategory,Gene_Data$Gene,"Sinica")
names(Sin_dNdS_one)<-c("dNdS","Type","Category","Gene","Comparison")
Sin_dNdS_two <- Sin_dNdS_one[complete.cases(Sin_dNdS_one), ]
Sin_dNdS=Sin_dNdS_two[Sin_dNdS_two$Type == "N-mt", ]

# Perform the steps described above but with Fascicularis and Silenus data
Fas_dNdS_one <- data.frame(Gene_Data$Fas.dN.dS,Gene_Data$Type,
                       Gene_Data$FuncCategory,Gene_Data$Gene,"Fascicularis")
names(Fas_dNdS_one)<-c("dNdS","Type","Category","Gene","Comparison")
Fas_dNdS_two <- Fas_dNdS_one[complete.cases(Fas_dNdS_one), ]
Fas_dNdS=Fas_dNdS_two[Fas_dNdS_two$Type == "N-mt", ]

# Combine the data above
Dnds <- rbind(Arc_dNdS,Sin_dNdS, Fas_dNdS)

# Create a boxplot with three separate boxes for each species pair
png("Figures/Figure5.png",res=300,height = 3.5,width = 4.5,units="in",pointsize = 10)
boxplot(dNdS~Comparison, data=Dnds, 
        col=c("#fdae61","#d7191c","#2c7bb6"),
        ylab="dN/dS", xlab="",
        cex.main=1.8, cex.axis=1.15, cex.lab=1.25)
dev.off()
# Conduct ANOVA and summarize results
res.aov <- aov(dNdS ~ Comparison, data = Dnds)
sink("Stats/dNdS_ANOVA_results.txt")
summary(res.aov)
sink()

# Define a data frame with Arctoides branch extension
# Then assign names to columns to allow for consistency with other branch extension data
# Then exclude genes with no branch extension data and nuclear genes
# Then include only genes with a majority Arc-Sin topology
Sin_Branches_one=data.frame(Gene_Data$Arc.Sin,Gene_Data$Extension,Gene_Data$Type,
                        Gene_Data$FuncCategory,Gene_Data$Gene,"Sinica")
names(Sin_Branches_one)<-c("Topo","Extension","Type","Category","Gene","Comparison")
Sin_Branches_two <- Sin_Branches_one[complete.cases(Sin_Branches_one), ]
Sin_Branches_three=Sin_Branches_two[Sin_Branches_two$Type == "N-mt", ]
Sin_Branches=Sin_Branches_three[which(Sin_Branches_three[,1]>864),]

# Perform the steps described above but with Arc-Fas topology
Fas_Branches_one=data.frame(Gene_Data$Arc.Fas,Gene_Data$Extension,Gene_Data$Type,
                            Gene_Data$FuncCategory,Gene_Data$Gene,"Fascicularis")
names(Fas_Branches_one)<-c("Topo","Extension","Type","Category","Gene","Comparison")
Fas_Branches_two <- Fas_Branches_one[complete.cases(Fas_Branches_one), ]
Fas_Branches_three=Fas_Branches_two[Fas_Branches_two$Type == "N-mt", ]
Fas_Branches=Fas_Branches_three[which(Fas_Branches_three[,1]>864),]

# Combine the data above then create a table of based on majority topology and
# branch extension. Then make this table look clean and conduct a Fisher's
# exact test where the hypothesis is that there is a greater amount of 
# extended branch lengths in Arc-Sin topology data
Branches=rbind(Sin_Branches, Fas_Branches)
Branches_table=table(Branches$Comparison, Branches$Extension)
png("Figures/Branch_Table.png",res=300,height = 3.5,width = 4.5,units="in",pointsize = 10)
grid.table(Branches_table)
dev.off()

sink("Stats/Branch_Exact_results.txt")
fisher.test(Branches_table, alternative = "greater")
sink()

# Recombination -----------------------------------------------------------

# Define a data frame with recombination and fdM statistics
# Then assign names to columns to allow for consistency with other statistics
# Then exclude genes with no fdM or recombination values
# Then only include genes with majority Arc-Sin topology
Sin_Recom_one=data.frame(Gene_Data$Arc.Sin,Gene_Data$RecomRate,Gene_Data$fdM,
                         Gene_Data$Type,Gene_Data$FuncCategory,Gene_Data$Gene,"Sinica")
names(Sin_Recom_one)<-c("Topo","Recombination","fdM","Type","Category","Gene","Comparison")
Sin_Recom_two <- Sin_Recom_one[complete.cases(Sin_Recom_one), ]
Sin_Recom_three=Sin_Recom_two[Sin_Recom_two$Type == "N-mt", ]
Sin_Recom=Sin_Recom_three[which(Sin_Recom_three[,1]>864),]

# Create table from above data and compute linear regression
Sin_table = table(Sin_Recom$Recombination, Sin_Recom$fdM)
SinMod <- lm(formula = Sin_Recom$Recombination ~ Sin_Recom$fdM)

# Perform the steps described above but with Arc-Fas topology genes
Fas_Recom_one=data.frame(Gene_Data$Arc.Fas,Gene_Data$RecomRate,Gene_Data$fdM,
                         Gene_Data$Type,Gene_Data$FuncCategory,Gene_Data$Gene,"Fascicularis")
names(Fas_Recom_one)<-c("Topo","Recombination","fdM","Type","Category","Gene","Comparison")
Fas_Recom_two <- Fas_Recom_one[complete.cases(Fas_Recom_one), ]
Fas_Recom_three=Fas_Recom_two[Fas_Recom_two$Type == "N-mt", ]
Fas_Recom=Fas_Recom_three[which(Fas_Recom_three[,1]>864),]

Fas_table = table(Fas_Recom$Recombination, Fas_Recom$fdM)
FasMod <- lm(formula = Fas_Recom$Recombination ~ Fas_Recom$fdM)

sink("Stats/Recombination_fdM_Regression_results.txt")
print("Arc-Sin Topology CM/Mb-fdM Regression")
summary(SinMod)
print("Arc-Fas Topology CM/Mb-fdM Regression")
summary(FasMod)
sink()

# Create expressions to paste over coming scatter plots
Sin_expression <- expression(paste("N = 56, ","P = 0.204, ",R^2," = 0.009"))
Fas_expression <- expression(paste("N = 23, ","P = 0.134, ",R^2," = 0.076"))

# Create scatter plots with labeled x (recombination) and y (fdM) axes
# Paste text over the scatter plots including the above expressions
png("Figures/Figure4.png",res=300,height = 5,width = 10,units="in",pointsize = 10)
par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow=c(1,2))
scatter.smooth(x=Sin_Recom$Recombination, y=Sin_Recom$fdM,
               pch = 19,
               col=c("#af8dc3"),
               xlab='Recombination Rate (cM/Mb)',
               ylab="fdM",
               cex.lab=2.3, cex.axis=2)
text(4, 0.05, labels = "A", cex=3.5)
text(3, -0.55, labels = Sin_expression, cex=2)
scatter.smooth(x=Fas_Recom$Recombination, y=Fas_Recom$fdM,
               pch = 19,
               col=c("#af8dc3"),
               xlab='Recombination Rate (cM/Mb)',
               ylab="",
               cex.lab=2.3, cex.axis=2)
text(0.35, 0.3, labels = "B",cex=3.5)
text(0.1, -0.45, labels = Fas_expression, cex=2)
dev.off()
