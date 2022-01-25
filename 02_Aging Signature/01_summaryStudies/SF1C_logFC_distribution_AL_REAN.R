## logFC of all genes per REAN study
#Source file generated from Make_AS_table.ipynb
#A. Svendsen Feb 2020

setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
library(tidyverse)
library(patchwork)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

mydata <- read.csv("00_sourceFiles/06_Aging_List_REAL_transposed_PCA.csv", header = T, stringsAsFactors = F)[,c(1:13)]
mydata.long <- gather(mydata[,-1], key = study, value = logFC, factor_key = TRUE)
mydata.long$platform <- c(rep("microarray", times = 1096 * 5), rep("scRNA-seq", times = 1096 *4), rep("bRNA-seq", times = 1096 * 3))
mydata.long$platform <- factor(mydata.long$platform, levels = c("microarray", "bRNA-seq", "scRNA-seq"))
mydata.long$study <- factor(mydata.long$study, levels = unique(mydata.long$study[order(mydata.long$platform)]))

ggplot(mydata.long, aes(x = study, y = logFC, fill = platform)) + 
  geom_boxplot(outlier.shape = 1,outlier.size = 7, outlier.alpha = 0.5) + 
  theme_pb() +
  theme(axis.text.x = element_text(angle = 90),
        aspect.ratio = 0.2)
