## Volcano plot of AS genes
#A. Svendsen Feb 2020

setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
library(tidyverse)
library(patchwork)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

logFC <- read.csv("00_sourceFiles/04_Aging_Signature_REAN_FC_NA.csv", header = T, stringsAsFactors = F, na.strings = "NA")[,c(1,14,15)]
pval <- read.csv("00_sourceFiles/04_Aging_Signature_REAN_pval_NA.csv", header = T, stringsAsFactors = F, na.strings = "NA")[,c(1,14)]


df <- merge(x = logFC, y = pval, by = "Genes")
colnames(df) <- c("Genes", "logFC", "consistency", "adjpval")

df$consistency <- factor(df$consistency, levels = levels(as.factor(df$consistency)))

ggplot(df, aes(x = logFC, y = -log(adjpval))) +
  #scale_color_viridis(discrete = TRUE, option = "D") +
  geom_vline(xintercept = c(-1,1),  color = "grey") +
  geom_hline(yintercept = 1, color = "grey") +
  scale_x_continuous(limits = c(-4,4)) +
  geom_point(alpha = 0.5, size = 7,shape = 21, aes(color = ifelse(logFC > 0, 'red', 'blue'))) +
  scale_color_identity() +
  theme_pb() + 
  theme(aspect.ratio = 1)

same.GO <- read.csv("03_Aging_Signature/merged_UP_DO_filtered.csv", stringsAsFactors = F, header = T)
same.GO.cur <- read.csv("03_Aging_Signature/merged_UP_DO_filtered_curated.csv", stringsAsFactors = F, header = T)
library(ggrepel)
ggplot(same.GO.cur, aes(x = -log(FDR.do,2), y = -log(FDR.up,2), label = term.description)) +
  geom_point(shape = 1, size = 7, color = "red") +
  geom_text_repel(nudge_y = 2,nudge_x = 5) +
  geom_abline() +
  scale_x_continuous(limits = c(0,30)) +
  scale_y_continuous(limits = c(0,30)) +
  theme_pb() +
  theme(aspect.ratio = 1)
