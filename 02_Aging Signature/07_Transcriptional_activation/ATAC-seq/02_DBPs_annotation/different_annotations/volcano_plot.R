## Broad peaks gene annotations
#A. Svendsen Jan 2021
library(tidyverse)
library(ggrepel)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R") 
setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/07_Transcriptional_activation/ATAC-seq/")

exon <- read.table("02_DBPs_annotation/different_annotations/04_overlap_broad_exons_wo.bed", stringsAsFactors = F, sep = "\t", header = F) %>%
  select(V1,V2,V3, V9, V11, V15, V17)
# intra.gen <- read.table("02_DBPs_annotation/different_annotations/04_overlap_broad_intragenic_wo.bed", stringsAsFactors = F, sep = "\t", header = F) %>%
#   select(V1,V2,V3, V9, V11, V15, V17)
inter.gen <-read.table("02_DBPs_annotation/different_annotations/04_overlap_broad_intrergenic_wo.bed", stringsAsFactors = F, sep = "\t", header = F) %>%
  select(V1,V2,V3, V9, V11)
inter.gen$V15 <- "intergenic"
inter.gen$V17 <- "."

rbind.data.frame(exon, inter.gen) %>%
  ggplot(aes(x = V9, y = -log10(V11), color = V15)) +
  geom_vline(xintercept = c(-1,1),  color = "grey") +
  geom_hline(yintercept = 1, color = "grey") +
  geom_point(size = 7, alpha = 0.5, shape = 21) +
  facet_grid(. ~ V15) +
  labs(y = "-log10(FDR)")+
  xlab(expression(atop("accessible sites",young %<->% aged))) +
  scale_y_continuous(limits = c(0,15)) +
  scale_x_continuous(limits = c(-5,5)) +
  theme_pb() +
  theme(aspect.ratio = 1,
        legend.position = "none")

