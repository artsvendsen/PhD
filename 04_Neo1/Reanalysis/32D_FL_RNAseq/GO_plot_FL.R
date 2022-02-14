library(tidyverse)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
setwd("/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_FL_RNAseq/")

#Import in one go
files <- list.files(pattern = "*table.txt")
GOs <- c()
for (file in files){
  print(file)
  temp <- c()
  cat <- unlist(strsplit(x = file, fixed = TRUE, "_"))[2]
  temp <- read.table(file, stringsAsFactors = F, header = T, sep = "\t") %>%
    #filter(Adjusted.P.value <= 0.05) %>%
    dplyr::select(Term, Overlap, P.value, Adjusted.P.value, Combined.Score) %>%
    mutate(type = rep(cat, times = length(.$Term)))
  
  GOs <- rbind.data.frame(GOs, temp)
}

#clean up
GO.cur <- GOs %>%
  as_tibble() %>%
  mutate(type = factor(type, levels = c("Biological", "Molecular", "Cellular"))) %>%
  separate(Overlap, sep = "/", into = c("genes", "total"),remove = F) %>%
  mutate(log2 = -log(Adjusted.P.value, base = 2)) %>%
  filter(genes > 2) %>%
  arrange(-Combined.Score) %>%
  slice(1:10) %>%
  arrange(log2) %>%
  mutate(Term = factor(Term, levels = .$Term))

#plot
GO.cur %>%
  ggplot(aes(x = Term, y = log2, fill = type)) + 
  coord_flip() +
  geom_bar(stat = "identity") + 
  #facet_wrap(vars(type)) +
  geom_text(aes(label = Overlap, hjust = 1), color = "white") +
  labs(x = "", y = "-log2(p-val)") +
  theme_pb() +
  theme(aspect.ratio = 2)
