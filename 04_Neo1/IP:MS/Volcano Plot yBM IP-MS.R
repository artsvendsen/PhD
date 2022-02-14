# Volcano plot for IP/MS data derived from young and old BM samples + GO plots
# see Exp 18
# input was dirved form analysis done using Provision (https://provision.shinyapps.io/provision/_w_17cb6b81/#)
# A. Svendsen dez 2021
library(biomaRt)
library(tidyverse)
library(ggrepel)
setwd("/Users/Art/Drive/PhD/Experiments/Neogenin1/Exp 18 - Proteomics (coIP:WES)/IP:MS/")
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")


#### Volcano plot (young0)
#Biomart info
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)

#Provision input
data <- read.csv("ProVisionProVision_YNeo1vsIgG.csv", header = T, stringsAsFactors = F) %>%
  as_tibble() %>%
  dplyr::select(GeneName,pValue, EffectSize)
colnames(data) <- c("uniprot_gn_id", "pval", "EffectSize")

anno <- getBM(
  attributes = c('external_gene_name', 'uniprot_gn_id'),
  filters    = 'uniprot_gn_id',
  values     = data$uniprot_gn_id,
  mart       = ensembl)

merged_y <- merge(anno, data, by = "uniprot_gn_id")
ggplot() +
  geom_vline(xintercept = c(-1,1),  color = "grey") +
  geom_hline(yintercept = 10^-0.05, color = "grey") +
  geom_point(data = filter(merged_y, pval > 0.1), size = 7,shape = 21,aes(x = EffectSize, y = -log10(pval)),color = "grey", alpha = 0.5) +
  geom_point(data = filter(merged_y, pval <= 0.05 & EffectSize > 0), size = 7,shape = 21, aes(x = EffectSize, y = -log10(pval)),color = "red", alpha = 0.5) +
  geom_text_repel(data = filter(merged_y, pval <= 0.05 & EffectSize > 0),
                  aes(x = EffectSize, y = -log10(pval), label = external_gene_name)) +
  labs(title ="IP/MS yBM vs. IgG") +
  scale_y_continuous(limits = c(0,6)) +
  scale_x_continuous(limits = c(-3,7.5)) +
  theme_pb() +
  theme(aspect.ratio = 1)


#### Volcano plot (aged)
#Provision input
data <- read.csv("ProVisionProVision_ONeo1vsIgG.csv", header = T, stringsAsFactors = F) %>%
  as_tibble() %>%
  dplyr::select(GeneName,pValue, EffectSize)
colnames(data) <- c("uniprot_gn_id", "pval", "EffectSize")

anno <- getBM(
  attributes = c('external_gene_name', 'uniprot_gn_id'),
  filters    = 'uniprot_gn_id',
  values     = data$uniprot_gn_id,
  mart       = ensembl)

merged_o <- merge(anno, data, by = "uniprot_gn_id")
ggplot() +
  geom_vline(xintercept = c(-1,1),  color = "grey") +
  geom_hline(yintercept = 10^-0.05, color = "grey") +
  geom_point(data = filter(merged_o, pval > 0.1), size = 7,shape = 21,aes(x = EffectSize, y = -log10(pval)),color = "grey", alpha = 0.5) +
  geom_point(data = filter(merged_o, pval <= 0.1 & EffectSize > 0), size = 7,shape = 21, aes(x = EffectSize, y = -log10(pval)),color = "red", alpha = 0.5) +
  geom_text_repel(data = filter(merged_o, pval <= 0.1 & EffectSize > 0),
                  aes(x = EffectSize, y = -log10(pval), label = external_gene_name)) +
  labs(title ="IP/MS agedBM vs. IgG") +
  #scale_y_continuous(limits = c(0,6)) +
  #scale_x_continuous(limits = c(-3,7.5)) +
  theme_pb() +
  theme(aspect.ratio = 1)


#Overlap young and aged (FDR 0.1)
overlap <- gplots::venn(list(Y = filter(merged_y, pval <= 0.1 & EffectSize > 0)$external_gene_name,
                             O = filter(merged_o, pval <= 0.1 & EffectSize > 0)$external_gene_name))

clipr::write_clip(attr(overlap, "intersections")$`Y:O`)


# GO analysis (young)
# input used at EnrichR <- clipr::write_clip(filter(merged, pval <= 0.05 & EffectSize > 0)[,2])
#Import in one go
files <- list.files(pattern = "*table.txt")
GOs <- c()
for (file in files){
  print(file)
  temp <- c()
  cat <- unlist(strsplit(x = file, fixed = TRUE, "_"))[2]
  temp <- read.table(file, stringsAsFactors = F, header = T, sep = "\t") %>%
    filter(Adjusted.P.value <= 0.05) %>%
    dplyr::select(Term, Overlap, Adjusted.P.value, Combined.Score) %>%
    mutate(type = rep(cat, times = length(.$Term)))
  
  GOs <- rbind.data.frame(GOs, temp)
}

#clean up
GO.cur <- GOs %>%
  as_tibble() %>%
  mutate(type = factor(type, levels = c("Biological", "Molecular", "Cellular"))) %>%
  mutate(log2 = -log(Adjusted.P.value, base = 2)) %>%
  filter(!grepl('Fc', Term)) %>%
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
  