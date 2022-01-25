#Age distribution of individual studies
#A.Svendsen Nov 2019

library(tidyverse)
library(patchwork)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/01_summaryStudies/")

df <- read.csv("SF1A_age_distribution.csv", header = T, stringsAsFactors = F)
df$age <- factor(df$age, levels = c("young", "aged"))
head(df)
class(df$months)

p <- ggplot(df, aes(x = age, y = months, fill = age)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") + 
  geom_jitter(aes(color = age),shape = 1, size = 7, width = 0.1) +
  scale_y_continuous(breaks = c(1,5,10,20,25,30)) +
  scale_color_manual(values = c("grey", "red")) +
  scale_fill_manual(values = c("grey", "red")) +
  labs(x = " ", y = "lifespan (months)") +
  theme_pb() +
  theme(aspect.ratio = 1,
        legend.position = "none")

p
#(p | p) / (p | p)
