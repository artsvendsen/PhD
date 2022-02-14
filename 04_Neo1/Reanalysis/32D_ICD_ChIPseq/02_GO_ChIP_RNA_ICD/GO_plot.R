library(tidyverse)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
setwd("/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/02_GO_ChIP_RNA_ICD/")

#Import in one go
files <- list.files()
GOs <- c()
for (file in files){
  print(file)
  temp <- c()
  cat <- unlist(strsplit(x = file, fixed = TRUE, "_"))[2]
  temp <- read.table(file, stringsAsFactors = F, header = T, sep = "\t") %>%
    dplyr::select(Term, Overlap, P.value) %>%
    mutate(type = rep(cat, times = length(.$Term)))
  
  GOs <- rbind.data.frame(GOs, temp)
}

#clean up
GO.cur <- GOs %>%
  mutate(type = factor(type, levels = c("Biological", "Molecular", "Cellular"))) %>%
  filter(P.value <= 0.05) %>%
  mutate(log2 = -log(P.value, base = 2)) %>%
  arrange(log2) %>%
  mutate(Term = factor(Term, levels = .$Term))

#plot
GO.cur %>%
  top_n(10) %>%
ggplot(aes(x = Term, y = log2, fill = type)) + 
  coord_flip() +
  geom_bar(stat = "identity") +
  #geom_text(aes(label = Overlap, hjust = 1), color = "white") +
  labs(x = "", y = "-log2(p-val)") +
  theme_pb() +
  theme(aspect.ratio = 2)
