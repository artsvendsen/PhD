library(tidyverse)
library(ggrepel)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R") 
setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/07_Transcriptional_activation/ATAC-seq/")

broad <- read.csv("01_diffbind/01_ATACseq_broad_p0_01_DBP_diffbind.csv", header = T, stringsAsFactors = F)
broad$peak <- "broad"
narrow <- read.csv("01_diffbind/02_ATACseq_narrow_p0_01_DBP_diffbind.csv", header = T, stringsAsFactors = F)
narrow$peak <- "narrow"


df <- rbind.data.frame(broad, narrow) %>%
  ggplot(aes(x = Fold, y = -log10(FDR), color = peak)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(. ~ peak) +
  labs(y = "-log10(FDR)")+
  xlab(expression(atop("accessible sites",young %<->% aged))) +
  scale_x_continuous(limits = c(-5,5)) +
  theme_bw() +
  theme(aspect.ratio = 1)

df
