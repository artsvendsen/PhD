library(tidyverse)
library(biomaRt)
library(viridis)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")


# Load BiomaRt
ensembl <- useMart("ensembl")
#if unresponsive
#ensembl <- useMart("ensembl", host="uswest.ensembl.org")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)

AS <- read.table("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/01_Aging_List_REAN.txt",
                 header = T, stringsAsFactors = F, sep = "\t") %>%
  filter(Freq_group >= 4) %>%
  filter(Gene != "")

anno <- getBM(
  attributes = c('external_gene_name', 'start_position', 'end_position'),
  filters    = 'external_gene_name',
  values     = AS$Gene,
  mart       = ensembl) %>%
  mutate(length = end_position - start_position) %>%
  dplyr::select(external_gene_name, length)
colnames(anno) <- c("Gene", "length")

df <- merge(x = AS, y = anno, by = "Gene")

ggplot(df, aes(x = length/1000, y = abs(Log2FC), color = Freq_group)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = lm) +
  scale_colour_viridis() +
  scale_x_log10() +
  labs(x = "gene lenght (kp)", y = "mean log2FC", colour = "AS consistency") +
  theme_pb() +
  theme(aspect.ratio = 1)

#subset meta and reanalysis
df2 <-data.frame(x = df$length, y = df$Log2FC)

#lm fit model for individual sets
fit <- lm(y ~ x, data = df2)
summary(fit)
