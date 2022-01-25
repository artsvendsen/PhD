#Correlation degree between meta and reanalyzed sets:
#A.Svendsen Jul 2019 (Reviewed Feb 2020)

#Sets that have both meta and reAnl are:
#Chambers/Wahlestedt/Sun/Grover

setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
library(tidyverse)
library(patchwork)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

#Sun
Sun.meta <- read.csv("00_sourceFiles/meta/Sun_2014.txt", header = T, stringsAsFactors = F, sep = "\t", skip = 1)[,c(2:3)]
Sun.reanl <- read.csv("00_sourceFiles/rean/Sun_GSE47817_DE_YO.csv", header = T, stringsAsFactors = F, sep = "\t", skip = 1)[,c(2:3)]
colnames(Sun.reanl) <- c("geneSymbol", "log2FC")

Sun <- inner_join(Sun.meta, Sun.reanl, by = "geneSymbol")
#plot(Sun$log2FC.x ~ Sun$log2FC.y)
fit1 <- lm(Sun$log2FC.x ~ Sun$log2FC.y)
summary(fit1)

p1 <- ggplot(Sun, aes(x = log2FC.x, y = log2FC.y)) +
  geom_hline(yintercept=0, color = "grey") +
  geom_vline(xintercept=0, color = "grey") +
  geom_point(shape = 1, size = 7, width = 0.1, alpha = 0.4) +
  geom_smooth(method = "lm", se = F) +
  scale_x_continuous(limits = c(-4,6)) +
  scale_y_continuous(limits = c(-4,6)) +
  labs(x = "log2FC-meta",
       y=  "log2FC-ReAn",
       title = "Sun") +
  geom_label(
    aes(x = 1, y = -2), hjust = 0,
    #color = "purple",
    label = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5)),
    size = 4) +
  theme_pb() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())
p1

#Wahlestedt
Wahlestedt.meta <- read.csv("00_sourceFiles/meta/Wahlestedt_2013.txt", header = T, stringsAsFactors = F, sep = "\t", skip = 1)[,c(2:3)]
Wahlestedt.reanl <- read.csv("00_sourceFiles/rean/Wahlestedt_GSE44923_DE_YOn.csv", header = T, stringsAsFactors = F, sep = "\t", skip = 1)[,c(2:3)]
colnames(Wahlestedt.reanl) <- c("gene", "log2FC")

Wahlestedt <- inner_join(Wahlestedt.meta, Wahlestedt.reanl, by = "gene")
#plot(Wahlestedt$log2FC.x ~ Wahlestedt$log2FC.y)
fit1 <- lm(Wahlestedt$log2FC.x ~ Wahlestedt$log2FC.y)
summary(fit1)

p2 <- ggplot(Wahlestedt, aes(x = log2FC.x, y = log2FC.y)) +
  geom_hline(yintercept=0, color = "grey") +
  geom_vline(xintercept=0, color = "grey") +
  geom_point(shape = 1, size = 7, width = 0.1, alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  scale_y_continuous(limits = c(-3,6)) +
  scale_x_continuous(limits = c(-3,6)) +
  labs(x = "log2FC-meta",
       y=  "log2FC-ReAnl",
       title = "Wahlestedt") +
  geom_label(
    aes(x = 0.25, y = -2), hjust = 0,
    #color = "purple",
    label = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5)),
    size = 4) +
  theme_pb() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())
p2

#Chambers
Chambers.meta <- read.csv("00_sourceFiles/meta/Chambers_2007.txt", header = T, stringsAsFactors = F, sep = "\t", skip = 1)[,c(2:3)]
Chambers.reanl <- read.csv("00_sourceFiles/rean/Chambers_GSE6503_YOn.csv", header = T, stringsAsFactors = F, sep = "\t", skip = 1)[,c(2:3)]
colnames(Chambers.reanl) <- c("Genes", "log2FC")

Chambers <- inner_join(Chambers.meta, Chambers.reanl, by = "Genes")
#plot(Chambers$log2FC.x ~ Chambers$log2FC.y)
fit1 <- lm(Chambers$log2FC.x ~ Chambers$log2FC.y)
summary(fit1)

p3 <- ggplot(Chambers, aes(x = log2FC.x, y = log2FC.y)) +
  geom_hline(yintercept=0, color = "grey") +
  geom_vline(xintercept=0, color = "grey") +
  geom_point(shape = 1, size = 7, width = 1, alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  scale_x_continuous(limits = c(-5,7)) +
  scale_y_continuous(limits = c(-5,7)) +
  labs(x = "log2FC-meta",
       y=  "log2FC-ReAnl",
       title = "Chambers") +
  geom_label(
    aes(x = 0.25, y = -4), hjust = 0,
    #color = "purple",
    label = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5)),
    size = 4) +
  theme_pb() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())
p3

#Grover
Grover.meta <- read.csv("00_sourceFiles/meta/Grover_2016.txt", header = T, stringsAsFactors = F, sep = "\t", skip = 1)[,c(2:3)]
colnames(Grover.meta) <- c("Gene", "log2FC")
Grover.reanl <- read.csv("00_sourceFiles/rean/Grover_GSE70657_DE_OY.csv", header = T, stringsAsFactors = F, sep = "\t", skip = 1)[,c(2:3)]
colnames(Grover.reanl) <- c("Gene", "log2FC")

Grover <- inner_join(Grover.meta, Grover.reanl, by = "Gene")
#plot(Grover$log2FC.x ~ Grover$log2FC.y)
fit1 <- lm(Grover$log2FC.x ~ Grover$log2FC.y)
summary(fit1)

p4 <- ggplot(Grover, aes(x = log2FC.x, y = log2FC.y)) +
  geom_hline(yintercept=0, color = "grey") +
  geom_vline(xintercept=0, color = "grey") +
  geom_point(shape = 1, size = 7, witdh = 1, alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  scale_x_continuous(limits = c(-0,6)) +
  scale_y_continuous(limits = c(-0,6)) +
  labs(x = "log2FC-meta",
       y=  "log2FC-ReAnl",
       title = "Grover") +
  geom_label(
    aes(x = 0.25, y = 4), hjust = 0,
    #color = "purple",
    label = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5)),
    size = 4) +
  theme_pb() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())
p4
