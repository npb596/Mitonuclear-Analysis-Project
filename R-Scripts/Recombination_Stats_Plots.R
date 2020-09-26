# Set working directory to introgression data
setwd("~/Desktop/Macaque_files/mitonuclear_project/pop_gen_analysis/Introgression_Data")

# Define table from Arc-Sin recombinationa data file
# Calculate linear regression model and then summarize the ouput on console
Sin_table = read.table("Arc_Sin_recom_data.txt")
SinMod <- lm(formula = Sin_table$V2 ~ Sin_table$V3)
summary(SinMod)

# Define table from Arc-Fas recombinationa data file
# Calculate linear regression model and then summarize the ouput on console
Fas_table = read.table("Arc_Fas_recom_data.txt")
FasMod <- lm(formula = Fas_table$V2 ~ Fas_table$V3)
summary(FasMod)

# Create text containing summary stats of Arc-Sin and Arc-Fas files
Sin_expression <- expression(paste("N = 56, ","P = 0.204, ",R^2," = 0.009"))
Fas_expression <- expression(paste("N = 23, ","P = 0.134, ",R^2," = 0.076"))

# Plot Arc-Sin and Arc-Fas data together including text generated above
par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow=c(1,2))
scatter.smooth(x=Sin_table$V2, y=Sin_table$V3,
               pch = 19,
               col=c("#af8dc3"),
               xlab='Recombination Rate (cM/Mb)',
               ylab="fdM",
               cex.lab=2.5, cex.axis=2)
               text(4, 0.05, labels = "A", cex=3.5)
               text(3, -0.55, labels = Sin_expression, cex=2)
scatter.smooth(x=Fas_table$V2, y=Fas_table$V3,
               pch = 19,
               col=c("#af8dc3"),
               xlab='Recombination Rate (cM/Mb)',
               ylab="",
               cex.lab=2.5, cex.axis=2)
               text(0.35, 0.3, labels = "B",cex=3.5)
               text(0.1, -0.45, labels = Fas_expression, cex=2)
