#Author Akram Mohammed
#November 2018
#Kamaleswaran Lab

#Set working directory for download
setwd("C:/Users/amoham18/Documents/SepsisGenomics")

library(affy)
library(limma)

data1 <- read.csv("groupbyGene_SepticShock_2.csv", comment.char="#")

table <-read.csv("groupbyGene_SepticShock_2_Transpose.csv",header=T)
tableSub <- subset(table, select = c(X0:X180))

genes <- table$index

design <- model.matrix(~ 0+factor(c(data1$Class)))
colnames(design) <- c("S", "NS")
contrast.matrix <- makeContrasts(diff=NS-S, levels=design)

fit <- lmFit(tableSub, design = design) # What should the design be?
fit2 <- contrasts.fit(fit, contrast.matrix)
options(digits=3)
fit3 <- eBayes(fit2)

writefile = topTable(fit3, number=Inf, genelist=genes, adjust.method = "FDR", sort.by="logFC")

write.csv(writefile, file="DEG_Septic_Shock_2_FDR.csv")
