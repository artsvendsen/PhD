## GO terms of Aging signature
#Data was extracted from downloaded STRING enrichment.function.tsv files for ALL, UP and DO. Further filtering was done using REVIGO
#A. Svendsen Feb 2020

setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
library(tidyverse)
library(patchwork)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

all <- read.csv("03_Aging_Signature/F3A_GO_function_ALL_filtered.csv", header = T, stringsAsFactors = F)[,c(1:5)] %>%
  mutate(log = -log(FDR,2)) %>%
  arrange(log) %>%
  mutate(term.description = factor(term.description, levels = .$term.description)) %>%
  ggplot(aes(x = term.description, y = log)) + 
  coord_flip() +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(observed,"/",background), hjust = 1), color = "white") +
  labs(x = "", y = "-log2FDR") +
  theme_pb() +
  theme(aspect.ratio = 2)

all

UP <- read.csv("03_Aging_Signature/Enrichment_analysis/zz_unique_Process_UP_filtered_original_merged.csv", header = T, stringsAsFactors = F)[,c(2:6)] %>%
  filter(observed.gene.count / background.gene.count > 0.01) %>%
  mutate(log = -log(false.discovery.rate,2)) %>%
  arrange(log) %>%
  mutate(term.description = factor(term.description, levels = .$term.description)) %>%
  tail(10) %>%
  ggplot(aes(x = term.description, y = log)) + 
  coord_flip() +
  geom_bar(stat = "identity", fill = "red") +
  geom_text(aes(label = paste0(observed.gene.count,"/",background.gene.count), hjust = 1), color = "white") +
  labs(x = "", y = "-log2FDR") +
  theme_pb() +
  theme(aspect.ratio = 2)

UP

DO <- read.csv("03_Aging_Signature/Enrichment_analysis/zz_unique_Process_DO_filtered_original_merged.csv", header = T, stringsAsFactors = F)[,c(2:6)] %>%
  filter(observed.gene.count / background.gene.count > 0.01) %>%
  mutate(log = -log(false.discovery.rate,2)) %>%
  arrange(log) %>%
  mutate(term.description = factor(term.description, levels = .$term.description)) %>%
  tail(10) %>%
  ggplot(aes(x = term.description, y = log)) + 
  coord_flip() +
  geom_bar(stat = "identity", fill = "blue") +
  geom_text(aes(label = paste0(observed.gene.count,"/",background.gene.count), hjust = 1), color = "white") +
  labs(x = "", y = "-log2FDR") +
  theme_pb() +
  theme(aspect.ratio = 2)

DO

path.UP <- read.csv("03_Aging_Signature/F3D_pathway_UP.csv", header = T, stringsAsFactors = F)[,c(1:5)] %>%
  mutate(log = -log(FDR,2)) %>%
  arrange(log) %>%
  mutate(term.description = factor(term.description, levels = .$term.description)) %>%
  ggplot(aes(x = term.description, y = log)) + 
  coord_flip() +
  geom_bar(stat = "identity", fill = "red") +
  geom_text(aes(label = paste0(observed,"/",background), hjust = 1), color = "white") +
  labs(x = "", y = "-log2FDR") +
  theme_pb() +
  theme(aspect.ratio = 2)

path.UP

path.DO <- read.csv("03_Aging_Signature/F3D_pathway_DO.csv", header = T, stringsAsFactors = F)[,c(1:5)] %>%
  mutate(log = -log(FDR,2)) %>%
  arrange(log) %>%
  mutate(term.description = factor(term.description, levels = .$term.description)) %>%
  ggplot(aes(x = term.description, y = log)) + 
  coord_flip() +
  geom_bar(stat = "identity", fill = "blue") +
  geom_text(aes(label = paste0(observed,"/",background), hjust = 1), color = "white") +
  labs(x = "", y = "-log2FDR") +
  theme_pb() +
  theme(aspect.ratio = 2)

path.DO