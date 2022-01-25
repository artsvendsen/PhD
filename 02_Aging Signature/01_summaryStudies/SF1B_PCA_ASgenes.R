## PCA plot for ReAnl aging signature
#Source file generated from Make_AS_table.ipynb
#File compiles logFCs of AS genes (0 values mean not found, so should be kept)
#A. Svendsen Feb 2020

setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
library(tidyverse)
library(patchwork)
library(ggrepel)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

#PCA calculation
mydata <- read.csv("00_sourceFiles/06_Aging_Signature_REAL_transposed_PCA.csv", header = T, stringsAsFactors = F)[,c(1:13)]
matrix <- data.matrix(mydata[,-1])
rownames(matrix) <- mydata$Genes
head(matrix)

#PCA
pca <- prcomp(t(matrix), scale. = TRUE, center = TRUE)

#PCA variation contribution
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 2)

#Labels
pca.data <- data.frame(Sample=rownames(pca$x), 
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       platform = c("microarray", "microarray", "microarray", "microarray",
                                    "microarray", "scRNA-seq", "scRNA-seq", "scRNA-seq",
                                    "scRNA-seq", "bRNA-seq","bRNA-seq", "bRNA-seq"))
pca.data$platform <- factor(pca.data$platform, levels = c("microarray", "bRNA-seq", "scRNA-seq"))

#PCA plot
ggplot(pca.data, aes(x=X, y=Y, color = platform)) +
  geom_point(shape = 1, size = 7, width = 0.1) +
  geom_text_repel(aes(label = Sample), show.legend = FALSE) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  scale_y_continuous(limits = c(-20,20))+
  scale_x_continuous(limits = c(-20,20))+
  ggtitle("PCA - ReAnl aging sig") +
  theme_pb() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)
