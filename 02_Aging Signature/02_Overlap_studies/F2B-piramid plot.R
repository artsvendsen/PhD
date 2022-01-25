## Distribution of number of genes and their consistency - pyramid plot
#Data extracted from gene_in_AS_pyramid ("generate freqs for pyramid")
#A. Svendsen Feb 2020

setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
library(tidyverse)
library(patchwork)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

mydata <- read.csv("00_sourceFiles/07_summary_num_genes_studies.csv", header = T, stringsAsFactors = F, na.strings = "NA")

#Pyramid Plot
ggplot(data = mydata, aes(x = as.factor(studies),
                          fill = analysis)) +
  geom_bar(data=subset(mydata, analysis =="meta"), aes(y = log2), stat = 'identity') + 
  geom_text(data=subset(mydata, analysis =="meta" & studies < 12), aes(y = log2, label=num.genes), hjust=-0.5, color = "white") +
  geom_bar(data=subset(mydata, analysis =="reanl"), aes(y = log2),stat = 'identity') + 
  geom_text(data=subset(mydata, analysis =="reanl" & studies <= 12), aes(y = log2, label=num.genes), hjust=1.5, color = "black") +
  scale_fill_manual(values = c("#984ea3", "#4daf4a")) +
  scale_y_continuous(breaks=seq(-10,10,5),labels=abs(seq(-10,10,5))) +
  labs(x = "number of studies", y = "log2(number of reported genes)") +
  coord_flip() +
  theme_pb() +
  theme(
    legend.text = element_text(size = 15, color = 'black'),
    legend.title = element_text(size = 15, color = 'black'),
    aspect.ratio = 1,
    panel.grid = element_blank()
  )