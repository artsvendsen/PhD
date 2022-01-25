## Distribution of number of genes and their consistency - pyramid plot
#Data extracted from gene_in_AS_pyramid ("generate freqs for pyramid")
#A. Svendsen Feb 2020

setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
library(tidyverse)
library(patchwork)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

mydata <- read.csv("00_sourceFiles/07_summary_num_genes_studies.csv", header = T, stringsAsFactors = F, na.strings = "NA")

#Top 5% stacked plot
mydata$analysis.top5 <- factor(mydata$analysis.top5, levels = c("meta-analysis", "re-analysis", "NS"))
ggplot(mydata, aes(x = analysis, y = num.genes, fill = analysis.top5)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#984ea3", "#4daf4a", "#bdbdbd")) +
  labs(x = " ", y = "# of DE genes") +
  theme_pb() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)