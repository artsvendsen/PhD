#Stacked bar plot with number of DE genes from meta and renalysis
#Data taken from summary table (Pub_list.xlsx)
#A.Svendsen Jul 2019 (reviwed Jan 2020)
library(tidyverse)
library(patchwork)
setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

#Number of DE genes in meta and reanalys
mydata <- read.csv("01_summaryStudies/F1B:C_Summary_number_DEG.csv.csv", header = T, stringsAsFactors = F)
mydata$total.genes <- mydata$up +mydata$down

#Metanalysis
meta <- mydata %>%
  filter(analysis == "meta") %>%
  arrange(desc(total.genes)) %>%
  mutate(study = factor(study, levels = .$study)) %>%
  gather(key =  "direction", value = "genes", up, down) %>%
  mutate(direction = factor(direction, levels = c("up", "down"))) %>%
  ggplot(aes(x = study, y = genes, fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('red', 'blue')) +
  labs(y = "# of reported DE genes",
       x = " ",
       title = "metanalysis") +
  theme_pb() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = 'black'),
    #legend.position = 'none',
    aspect.ratio = 1)

#Reanalysis
re <- mydata %>%
  filter(analysis == "reanalysis") %>%
  arrange(desc(total.genes)) %>%
  mutate(study = factor(study, levels = .$study)) %>%
  gather(key =  "direction", value = "genes", up, down) %>%
  mutate(direction = factor(direction, levels = c("up", "down"))) %>%
  ggplot(aes(x = study, y = genes, fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('red', 'blue')) +
  labs(y = " ",
       x = " ",
       title = "reanalysis") +
  theme_pb() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = 'black'),
    legend.text = element_text(size = 15, color = 'black'),
    legend.title = element_text(size = 15, color = 'black'),
    aspect.ratio = 1
  )


#Summary of Age of all studies
mydata2 <- read.csv("01_summaryStudies/F1B:C_Studies_ages_summary.csv", header = T, stringsAsFactors = F)
cronology <- c("Rossi", "Chambers", "Noda", "Norddahl", "Bersenev", "Wahlestedt", 
               "Beerman", "Sun", "Quere", "Flach", "Kowalczyk", "Grover", "Kirshner", "Maryanovich", "Mann", "Lazare")
summary.age <- mydata2 %>%
  mutate(study = factor(study),
         study = factor(study, levels = rev(cronology))) %>%
  ggplot(aes(x = age.start, xend = age.end, y = study, yend = study, color = platform)) +
  geom_segment(size = 3) +
  labs(y = " ",
       x = "Age (in months)") +
  theme_pb() +
  theme(
    legend.text = element_text(size = 15, color = 'black'),
    legend.title = element_text(size = 15, color = 'black')
  ) 
meta
re
summary.age
