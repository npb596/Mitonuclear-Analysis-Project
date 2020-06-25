setwd("~/Desktop/Macaque_files/mitonuclear_project/pop_gen_analysis")

Sin_table = read.table("Introgression_Data/Arc_Sin_recom_data.txt")

#SinRecom <- log(Sin_table$V2)

SinMod <- lm(formula = Sin_table$V2 ~ Sin_table$V3)

summary(SinMod)

cor(Sin_table$V2, Sin_table$V3, use = "complete.obs")

#pdf("Arc_Sin_recom_plots.pdf")
plot(SinMod)

#dev.off()

#pdf("Introgression_results/Sinica_recombination_CDS.pdf")

#dev.off()

Fas_table = read.table("Introgression_Data/Arc_Fas_recom_data.txt")

#FasRecom <- log(Fas_table$V2)

FasMod <- lm(formula = Fas_table$V2 ~ Fas_table$V3)

summary(FasMod)

cor(Fas_table$V2, Fas_table$V3)

#pdf("Arc_Fas_recom_plots.pdf")
plot(FasMod)
#dev.off()

#pdf("Introgression_results/Fascicularis_recombination_CDS.pdf")
par(mfrow=c(1,2))
scatter.smooth(x=Sin_table$V2, y=Sin_table$V3,
               pch = 19,
               col=c("purple"),
               xlab='Recombination Rate (cM/Mb)',
               ylab="fdM",
               cex.lab=2.5, cex.axis=2)
               text(4, 0.05, labels = "A", cex=3)
scatter.smooth(x=Fas_table$V2, y=Fas_table$V3,
               pch = 19,
               col=c("purple"),
               xlab='Recombination Rate (cM/Mb)',
               ylab="",
               cex.lab=2.5, cex.axis=2)
               text(0.35, 0.3, labels = "B",cex=3)
#dev.off()
