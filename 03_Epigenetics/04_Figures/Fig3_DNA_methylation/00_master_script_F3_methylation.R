##### Script containing all the necessary data (unless explicity stated) to generate
##### all the figures for Figure 4 of the Epigenentics Chapter, relating to DNA methylation
library(tidyverse)
library(GenomicFeatures)
library(ChIPseeker)
library(biomaRt)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/")

# Dmnt3b expression AS ----
read.csv("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/04_Aging_Signature_REAN_FC_NA.csv",
               header = T, stringsAsFactors = F) %>%
  filter(Genes == "Dnmt3b") %>%
  dplyr::select(-c(Avg, Consistency, Order)) %>%
  gather(study, value, -Genes) %>%
  filter(is.na(value) != TRUE) %>%
  ggplot(aes(x = Genes, y = value, color = study)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_point() +
  stat_summary(fun = mean,
               geom="errorbar",
               color = "red",
               aes(ymax = ..y.., ymin = ..y..),
               width = 0.3) +
  scale_y_continuous(limits = c(-2.5,2.5)) +
  labs(x = "", y = "Mean RNA Fold-Change") +
  theme_pb() +
  theme(aspect.ratio = 3)
  

# Colleting, overlapping and promoter annotation all CpGs from Beerman, Taiwo and Sun ----
## Load and collect promoters from Beerman and Sun -----
## BiomaRt
ensembl <- useMart("ensembl")
#if unresponsive
#ensembl <- useMart("ensembl", host="uswest.ensembl.org")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)

#Build a TxDb object from the actual annotation that is used to align samples (running locally):
#Assemble TxDb (takes a while)
# BiomartTxDb <- makeTxDbFromGFF(file = "/Users/Art/Drive/PhD/Scripts/GENCODE/gencode.vM13.annotation.gff3",
#                                       format = "gff3",
#                                       organism = "Mus musculus")
# saveDb(BiomartTxDb, file = "/Users/Art/Drive/PhD/Experiments/Epigenetics/BiomartTxDb_M13.sqlite")

#Lazy load
BiomartTxDb <- loadDb(file = "/Users/Art/Drive/PhD/Experiments/Epigenetics/BiomartTxDb_M13.sqlite")
seqlevelsStyle(BiomartTxDb) <- "UCSC"


#load peaks into GRrange format
files <- list(Beerman_hyper = "01_Beerman/01_raw_data/03_youngtoold_gain_mm10.bed",
              Beerman_hypo = "01_Beerman/01_raw_data/03_youngtoold_loss_mm10.bed",
              Sun_hyper = "03_Sun/01_raw_data/02_Sun_hypermethylated_mm10.bed",
              Sun_hypo = "03_Sun/01_raw_data/02_Sun_hypomethylated_mm10.bed")

peakAnnoList <- lapply(files, annotatePeak, TxDb=BiomartTxDb, tssRegion=c(-500, 500), verbose=FALSE)

# Get promoters
extract.promoter <- function(i){
  #Get gene annotations from reference annotation
  temp <- dplyr::as_tibble(peakAnnoList[[i]]@anno@elementMetadata@listData) %>%
    dplyr::select(geneChr, geneStart, geneEnd, annotation, geneId, distanceToTSS) %>%
    dplyr::mutate(geneId = substr(geneId, start = 1, stop = 18)) %>%
    dplyr::mutate(index = 1:length(geneChr)) %>%
    dplyr::filter(annotation == "Promoter")
  
  #Convert Ensembl IDs to Gene symbols
  anno <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name'),
    filters    = 'ensembl_gene_id',
    values     = temp$geneId,
    mart       = ensembl)
  colnames(anno) <- c("geneId", "Gene")
  
  #Collect real peaks back (ones that were in orgininal bed files, not gene annotations)
  temp2 <- as.data.frame(peakAnnoList[[i]]@anno@ranges[temp$index,]) %>%
    dplyr::mutate(peakStart = start - 1) %>% # start peaks get shifted 1bp downstream
    dplyr::mutate(peakEnd = end) %>%
    dplyr::select(-width)
  
  #Merge all parameters in a decent order
  temp3 <- cbind.data.frame(temp, temp2) %>%
    merge(., anno, by = "geneId") %>%
    dplyr::arrange(geneChr) %>%
    dplyr::mutate(geneChr = paste0("chr",geneChr)) %>%
    dplyr::select(geneChr, peakStart, peakEnd, distanceToTSS,Gene) %>%
    dplyr::select(Gene) %>%
    unique()
  
  return(temp3)
}
#test <- extract.feature(data = peakAnnoList, i = 2, promoter = T)
samples <- c("Beerman_hyper", "Beerman_hypo", "Sun_hyper", "Sun_hypo")
promoter.list <- list()
for (i in 1:4){
  print(samples[i])
  promoter.list[[samples[i]]] <- extract.promoter(i = i)
}

overlap.promoter <- gplots::venn(promoter.list, show.plot = FALSE, intersections=TRUE)
str(overlap.promoter)

## Import Taiwo  and AS data -----
## There were no information on Peak coordinates, so just importing gene symbols
taiwo_hyper <- read.table("02_Taiwo/01_Taiwo_hyper.txt", header = T, stringsAsFactors = F, sep = "\t")[,1] %>%
  str_to_title()
taiwo_hypo <- read.table("02_Taiwo/01_Taiwo_hypo.txt", header = T, stringsAsFactors = F, sep = "\t")[,1] %>%
  str_to_title()

# # Overlap Beerman, Taiwo, Sun promoter  -----
# # Taking genes which appears at least 2 out of 3 times in me3
# # hypermethylation in aging
# hyper.all <- gplots::venn(list(promoter.list[["Beerman_hyper"]], taiwo_hyper, promoter.list[["Sun_hyper"]]),
#                           show.plot = FALSE, intersections=TRUE)
# hyper.all <- c(
#   attr(hyper.all, "intersections")$'A:B',
#   attr(hyper.all, "intersections")$'A:C',
#   attr(hyper.all, "intersections")$'A:B:C')
# 
# ## hypomethylation in aging
# hypo.all <- gplots::venn(list(promoter.list[["Beerman_hypo"]], taiwo_hypo, promoter.list[["Sun_hypo"]]),
#                          show.plot = FALSE, intersections=TRUE)
# hypo.all <- c(
#   attr(hypo.all, "intersections")$'A:C',
#   attr(hypo.all, "intersections")$'B:C')
# 
# write.csv(file = "/Users/Art/Drive/PhD/Experiments/Epigenetics/04_Figures/Fig3_DNA_methylation/03_ALL_DNAmethylation_hypermethylated_promoter_Gene_symbols.csv",
#           hyper.all, row.names = F, quote = F)
# 
# write.csv(file = "/Users/Art/Drive/PhD/Experiments/Epigenetics/04_Figures/Fig3_DNA_methylation/03_ALL_DNAmethylation_hypomethylated_promoter_Gene_symbols.csv",
#           hypo.all, row.names = F, quote = F)




## Import AS gene symbols
AS <- read.table("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/01_Aging_List_REAN.txt",
                 header = T, stringsAsFactors = F, sep = "\t") %>%
  filter(Freq_group >= 4) %>%
  filter(Gene != "") %>%
  dplyr::select(Gene)

## Overlap Beerman, Taiwo, Sun promoter with AS genes -----
## Taking genes which appears at least 2 out of 3 times in me3 and are present in the AS
## hypermethylation in aging
hyper.all <- gplots::venn(list(promoter.list[["Beerman_hyper"]], taiwo_hyper, promoter.list[["Sun_hyper"]], AS),
                          show.plot = FALSE, intersections=TRUE)
hyper.all <- c(
  attr(hyper.all, "intersections")$'A:D',
  attr(hyper.all, "intersections")$'B:D',
  attr(hyper.all, "intersections")$'C:D',
  attr(hyper.all, "intersections")$'A:C:D')
## hypomethylation in aging
hypo.all <- gplots::venn(list(promoter.list[["Beerman_hypo"]], taiwo_hypo, promoter.list[["Sun_hypo"]], AS),
                         show.plot = FALSE, intersections=TRUE)
hypo.all <- c(
  attr(hypo.all, "intersections")$'A:D',
  attr(hypo.all, "intersections")$'B:D',
  attr(hypo.all, "intersections")$'C:D',
  attr(hypo.all, "intersections")$'A:C:D',
  attr(hypo.all, "intersections")$'B:C:D')

## Final list of genes which are hyper or hypomethylated in the AS (duplicates removed)
hyper.final <- attr(gplots::venn(list(hyper = hyper.all, hypo = hypo.all), 
                                 show.plot = FALSE, intersections=TRUE), "intersections")$hyper
hypo.final <- attr(gplots::venn(list(hyper = hyper.all, hypo = hypo.all),
                                show.plot = FALSE, intersections=TRUE), "intersections")$hypo
#write.csv(hypo.final, "hypo.final.csv", quote = F, row.names = F, col.names = F)
#write.csv(hyper.final, "hyper.final.csv", quote = F, row.names = F, col.names = F)
remove(hyper.all, hypo.all)

# 1.SF4A/F4A - AS genes and their methylation status ----
## Import AS with consistency and logFC and creat a extra methylation (status) col ----
AS.complete <- read.table("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/01_Aging_List_REAN.txt",
                          header = T, stringsAsFactors = F, sep = "\t") %>%
  filter(Freq_group >= 4) %>%
  filter(Gene != "")

#create methylation status column for AS.complete
status <- factor()
for (x in AS.complete$Gene){
  methylation <- c()
  if (x %in% hyper.final == TRUE){
    methylation <- "hyper" 
  }  else if (x %in% hypo.final == TRUE){
    methylation <- "hypo"
  } else {
    methylation <- "nop"
  }
  status <- rbind.data.frame(status, methylation)
}

#cbind status to AS.complete
AS.complete$status <- factor(status$X.nop.)
colnames(AS.complete) <- c("gene", "consistency", "fc", "methylation")

## final plot -----
AS.complete %>%
  ggplot(aes(x = fc, y = consistency, color = methylation)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 4, color = "grey") +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("#9e0142", "#3288bd", "grey")) +
  ggrepel::geom_text_repel(data = filter(AS.complete, methylation != "nop"), aes(x = fc, y = consistency, label = gene)) +
  labs(x = "mean LogFC", y = "AS consistency score") +
  scale_x_continuous(limits = c(-5,5)) +
  theme_pb() +
  theme(aspect.ratio = 1)

#Get start for the pie chart
AS.complete %>%
  mutate(direction = ifelse(fc > 0, "up", "do")) %>%
  group_by(methylation, direction) %>%
  count


# 2.SF4B - gene annotation of hyper and hypo beerman and sun ----

# Features number extraction
df <- c()
methylation <- c("hyper", "hypo", "hyper", "hypo")
group <- c("Beerman", "Beerman", "Sun", "Sun")
for (i in 1:4){
  temp <- peakAnnoList[[i]]@annoStat
  age <- rep(group[i], times = length(temp[,2]))
  mark <- rep(methylation[i], times = length(temp[,2]))
  df <- rbind.data.frame(df, cbind.data.frame(temp, age, mark))
}

# Genomic annotation plot for H3K4me3, H3K27me3 and H3K27me3 -----
df %>%
  mutate(age = factor(age, levels = c( "Beerman","Sun"))) %>%
  mutate(mark = factor(mark, levels = c("hyper", "hypo"))) %>%
  ggplot(aes(fill=Feature, y= Frequency, x = mark, label = round(Frequency, 2))) +
  geom_bar(position="stack", stat="identity") +
  labs(title = "Peak genomic annotation", x = "", y = "Frequency of annotation (%)", fill = "Anotation") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  facet_wrap(vars(age)) +
  theme_pb() +
  theme(aspect.ratio = 1)
