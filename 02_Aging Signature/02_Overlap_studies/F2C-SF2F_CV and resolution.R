## Saturation plots
#Data extracted from source file subsets_comp_devs.csv and 05_resolution_prediction.csv
#A. Svendsen Feb 2020

setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
library(tidyverse)
library(patchwork)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")


#CV
mydata1 <- read.csv("00_sourceFiles/08_resolution_CV.csv", header = T, stringsAsFactors = F) %>%
  gather(-parameter, -studies, key = "gene", value = "value") %>%
  spread(parameter, value = value)

p1 <- ggplot(mydata1, aes(x = as.numeric(mean), y = as.numeric(cv), color = as.factor(studies))) + geom_point(shape = 1, size = 7, alpha = 0.2) +
  scale_color_manual(values = c('#313695','#4575b4', '#74add1', '#abd9e9', '#e0f3f8', '#ffffbf', '#fee090', '#fdae61', '#f46d43', '#d73027', '#a50026')) +
  scale_y_continuous(limits = c(0,3)) +
  geom_smooth(method = "glm", se = FALSE) +
  labs(x = "mean score", y = "coeficient of variation (CV)", color = "subset of studies") +
  theme_pb() +
  theme(aspect.ratio = 1)
p1

#Resolution
mydata2 <- read.csv("00_sourceFiles/08_resolution_prediction.csv", header = T, stringsAsFactors = F)[1:19,] %>%
  gather("parameter", "fraction", -studies) %>%
  mutate(fraction = ifelse(parameter == "gain", fraction * 2, fraction))


p2 <- ggplot(filter(mydata2, parameter == "resolution"), aes(x = studies, y = fraction)) +
  geom_point(aes(), shape = 1, size = 7) +
  geom_smooth(se=F)+
  labs(x = "number of studies", y = "resolution") +
  theme_pb() +
  theme(legend.position = "none",
        aspect.ratio = 1)
p2
