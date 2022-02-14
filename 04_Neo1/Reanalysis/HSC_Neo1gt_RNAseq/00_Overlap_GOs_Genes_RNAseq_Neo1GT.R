library(tidyverse)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
setwd("/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/HSC_Neo1gt_RNAseq/02_GOs/")


#Young data
y <- list.files(pattern = "y_G*")
y.GOs <- c()

# import
for (file in y){
  print(file)
  temp <- c()
  cat <- unlist(strsplit(x = file, fixed = TRUE, "_"))[3]
  temp <- read.csv(file, stringsAsFactors = F, header = T, sep = "\t") %>%
    filter(Adjusted.P.value <= 0.05) %>%
    dplyr::select(Term, Overlap, Adjusted.P.value, Combined.Score) %>%
    mutate(type = rep(cat, times = length(.$Term)))
  
  y.GOs <- rbind.data.frame(y.GOs, temp)
}

#clean up
y.GO.cur <- y.GOs %>%
  mutate(type = factor(type, levels = c("Biological", "Molecular", "Cellular"))) %>%
  mutate(log2 = -log(Adjusted.P.value, base = 2)) %>%
  arrange(-Combined.Score) %>%
  slice(1:10) %>%
  arrange(log2) %>%
  mutate(Term = factor(Term, levels = .$Term))

#plot
y.GO.cur %>%
  ggplot(aes(x = Term, y = log2, fill = type)) + 
  coord_flip() +
  geom_bar(stat = "identity") + 
  #facet_wrap(vars(type), scales = "free") +
  geom_text(aes(label = Overlap, hjust = 1), color = "white") +
  labs(x = "", y = "-log2(p-val)", title = "young RNA-seq Neo1GT") +
  theme_pb() +
  theme(aspect.ratio = 2)


# Aged  
o <- list.files(pattern = "o_G*")
o.GOs <- c()
for (file in o){
  print(file)
  temp <- c()
  cat <- unlist(strsplit(x = file, fixed = TRUE, "_"))[3]
  temp <- read.csv(file, stringsAsFactors = F, header = T, sep = "\t") %>%
    filter(Adjusted.P.value <= 0.05) %>%
    dplyr::select(Term, Overlap, Adjusted.P.value, Combined.Score) %>%
    mutate(type = rep(cat, times = length(.$Term)))
  
  o.GOs <- rbind.data.frame(o.GOs, temp)
}

#clean up
o.GO.cur <- o.GOs %>%
  mutate(type = factor(type, levels = c("Biological", "Molecular", "Cellular"))) %>%
  mutate(log2 = -log(Adjusted.P.value, base = 2)) %>%
  arrange(-Combined.Score) %>%
  slice(1:10) %>%
  arrange(log2) %>%
  mutate(Term = factor(Term, levels = .$Term))

#plot
o.GO.cur %>%
  ggplot(aes(x = Term, y = log2, fill = type)) + 
  coord_flip() +
  geom_bar(stat = "identity") + 
  #facet_wrap(vars(type), scales = "free") +
  geom_text(aes(label = Overlap, hjust = 1), color = "white") +
  labs(x = "", y = "-log2(p-val)", title = "aged RNA-seq Neo1GT") +
  theme_pb() +
  theme(aspect.ratio = 2)

#common data
c <- list.files(pattern = "common_G*")
c.GOs <- c()

# import
for (file in c){
  print(file)
  temp <- c()
  cat <- unlist(strsplit(x = file, fixed = TRUE, "_"))[3]
  temp <- read.csv(file, stringsAsFactors = F, header = T, sep = "\t") %>%
    #filter(Adjusted.P.value <= 0.05) %>%
    dplyr::select(Term, Overlap, Adjusted.P.value, Combined.Score) %>%
    mutate(type = rep(cat, times = length(.$Term)))
  
  c.GOs <- rbind.data.frame(c.GOs, temp)
}

#clean up
c.GO.cur <- c.GOs %>%
  mutate(type = factor(type, levels = c("Biological", "Molecular", "Cellular"))) %>%
  mutate(log2 = -log(Adjusted.P.value, base = 2)) %>%
  arrange(-Combined.Score) %>%
  slice(1:10) %>%
  arrange(log2) %>%
  mutate(Term = factor(Term, levels = .$Term))

#plot
c.GO.cur %>%
  ggplot(aes(x = Term, y = log2, fill = type)) + 
  coord_flip() +
  geom_bar(stat = "identity") + 
  #facet_wrap(vars(type), scales = "free") +
  geom_text(aes(label = Overlap, hjust = 1), color = "white") +
  labs(x = "", y = "-log2(p-val)", title = "common RNA-seq Neo1GT") +
  theme_pb() +
  theme(aspect.ratio = 2)


#Gene overlap
setwd("/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/HSC_Neo1gt_RNAseq/")
y <- read.csv("01_RNAseq_DEG_20211210_HSC_GT.csv", stringsAsFactors = F, header = T)[,2]
o <- read.csv("01_RNAseq_DEG_20211210_HSC_GT_old.csv", stringsAsFactors = F, header = T)[,3]

myV2 <- nVennR::plotVenn(list(
  Young = y,
  Aged = o))

myV2 <- nVennR::plotVenn(nVennObj = myV2)
nVennR::showSVG(nVennObj = myV2,
        opacity = 0.1,
        borderWidth = 3,
        labelRegions = F,
        setColors = c("#797979",
                      "#FB0106")
)
