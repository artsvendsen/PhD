## Number of genes per cecullar compartment
#Aging signature genes were inquired in Uniprot for their subcecullar localization (uniprot folder)
#A. Svendsen Feb 2020

setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
library(tidyverse)
library(patchwork)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

ReAnl <- read.csv("00_sourceFiles/05_Aging_Signature_REAN_logFC_piramid.csv",
                  header = T, stringsAsFactors = F)[,c(1,5)]
head(ReAnl)
colnames(ReAnl) <- c("Gene", "direction")

Compartmetns <- read.csv("00_sourceFiles/Uniprot/Uniprot_query_location_CURATED_final2.csv", header = T, stringsAsFactors = F)[,-c(1,3)]

comp.ReAnl <- merge(x = ReAnl, y = Compartmetns, by = "Gene")

#Keyword count
keywords <- comp.ReAnl %>%
  gather(level, Keyword, -c(Gene, direction), factor_key = TRUE) %>%
  filter(Keyword != "") %>%
  arrange(Gene) %>%
  select(-3)

keywords$Keyword <- as.factor(keywords$Keyword)
levels(keywords$Keyword)

dat <- keywords%>%
  group_by(Keyword) %>%
  summarise(count.up = length(Keyword[direction == "up"]),
            count.do = length(Keyword[direction == "do"]),
            count.tot = length(Keyword)) %>%
  mutate(per.up = count.up / count.tot,
         per.do = count.do / count.tot) %>%
  arrange(desc(count.tot)) %>%
  mutate(Keyword = factor(Keyword, rev(Keyword)))

dat

p1 <- ggplot(data = dat, aes(x = Keyword)) +
  geom_bar(aes(y = count.up), stat = "identity", fill = "red") +
  geom_bar(aes(y = -count.do), stat = "identity", fill = "blue") +
  geom_text(aes(y = count.up, label = count.up), hjust = 1) +
  geom_text(aes(y = -count.do, label = count.do), hjust = -0.7) +
  labs(x = " ", y = "number of genes") +
  scale_y_continuous(limits = c(-20, 100), breaks=seq(-20,100,20), labels=abs(seq(-20,100,20))) +
  coord_flip() +
  theme_pb() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)
p1
