##### Script containing all the necessary data (unless explicity stated) to generate
##### all the figures for Figure 3 of the Epigenentics Chapter, relating to DNA methylation
library(tidyverse)
library(GenomicFeatures)
library(ChIPseeker)
library(biomaRt)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/")

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

## Genomic annotation hyper and hypo methylation sites for Sun and Beerman ----
#Features number extraction
df <- c()
samples <- c("Beerman_hyper", "Beerman_hypo", "Sun_hyper", "Sun_hypo")
for (i in 1:4){
  temp <- peakAnnoList[[i]]@annoStat
  sample <- rep(samples[i], times = length(temp[,2]))
  df <- rbind.data.frame(df, cbind.data.frame(temp, sample))
}

df %>%
  #mutate(sample = recode(sample, ATAC_broad_Y = "Young", ATAC_broad_O = "Aged", "ATAC_broad" = "Common")) %>%
  #mutate(sample = factor(sample, levels = c("Young", "Common", "Aged"))) %>%
  ggplot(aes(fill=Feature, y=Frequency, x = sample, label = round(Frequency,2))) +
  geom_bar(position="stack", stat="identity") +
  labs(title = "Peak genomic annotation", x = "", y = "Frequency of annotation (%)", fill = "Anotation") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_pb() +
  theme(aspect.ratio = 1)

# Get features
extract.features <- function(i){
  #Get gene annotations from reference annotation
  temp <- dplyr::as_tibble(peakAnnoList[[i]]@anno@elementMetadata@listData) %>%
    dplyr::select(geneChr, geneStart, geneEnd, annotation, geneId, distanceToTSS) %>%
    dplyr::mutate(geneId = substr(geneId, start = 1, stop = 18)) %>%
    dplyr::mutate(index = 1:length(geneChr))
    #dplyr::filter(annotation == "Promoter")
  
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
samples <- c("Beerman_hyper", "Beerman_hypo", "Sun_hyper", "Sun_hypo")
features.list <- list()
for (i in 1:4){
  print(samples[i])
  features.list[[samples[i]]] <- extract.features(i = i)
}

overlap.features <- gplots::venn(features.list)
str(overlap.features)

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
samples <- c("Beerman_hyper", "Beerman_hypo", "Sun_hyper", "Sun_hypo")
promoter.list <- list()
for (i in 1:4){
  print(samples[i])
  promoter.list[[samples[i]]] <- extract.promoter(i = i)
}

overlap.promoter <- gplots::venn(promoter.list)
str(overlap.promoter)

# Import AS gene symbols
AS <- read.table("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/01_Aging_List_REAN.txt",
                 header = T, stringsAsFactors = F, sep = "\t") %>%
  filter(Freq_group >= 4) %>%
  filter(Gene != "") %>%
  dplyr::select(Gene)

# Overlap Beerman, Sun and AS genes for all features
features.list[["AS"]] <- AS
overlap.features.AS <- gplots::venn(features.list)
str(overlap.features.AS)

# Overlap Beerman, Sun and AS genes for promoters only
promoter.list[["AS"]] <- AS
overlap.promoter.AS <- gplots::venn(promoter.list)
str(overlap.promoter.AS)




# Import Taiwo gene names
taiwo_hyper <- read.table("02_Taiwo/01_Taiwo_hyper.txt", header = T, stringsAsFactors = F, sep = "\t")[,1] %>%
  str_to_title()
taiwo_hypo <- read.table("02_Taiwo/01_Taiwo_hypo.txt", header = T, stringsAsFactors = F, sep = "\t")[,1] %>%
  str_to_title()

#overlap Beerman, Taiwo and Sun hypermethylated genes
hyper.all <- gplots::venn(list(features.list[["Beerman_hyper"]], taiwo_hyper, features.list[["Sun_hyper"]], AS))
hyper.all <- c(
  attr(hyper.all, "intersections")$'A:D',
               attr(hyper.all, "intersections")$'B:D',
               attr(hyper.all, "intersections")$'C:D',
               attr(hyper.all, "intersections")$'A:C:D')

#promoters
hyper.all <- gplots::venn(list(promoter.list[["Beerman_hyper"]], taiwo_hyper, promoter.list[["Sun_hyper"]], AS))
hyper.all <- c(
  attr(hyper.all, "intersections")$'A:D',
  attr(hyper.all, "intersections")$'B:D',
  attr(hyper.all, "intersections")$'C:D',
  attr(hyper.all, "intersections")$'A:C:D')

#overlap Beerman, Taiwo and Sun hypomethylated genes
hypo.all <- gplots::venn(list(features.list[["Beerman_hypo"]], taiwo_hypo, features.list[["Sun_hypo"]], AS))
hypo.all <- c(
  attr(hypo.all, "intersections")$'A:D',
              attr(hypo.all, "intersections")$'B:D',
              attr(hypo.all, "intersections")$'C:D',
              attr(hypo.all, "intersections")$'A:C:D',
              attr(hypo.all, "intersections")$'B:C:D')

#promoters
hypo.all <- gplots::venn(list(promoter.list[["Beerman_hypo"]], taiwo_hypo, promoter.list[["Sun_hypo"]], AS))
hypo.all <- c(
  attr(hypo.all, "intersections")$'A:D',
  attr(hypo.all, "intersections")$'B:D',
  attr(hypo.all, "intersections")$'C:D',
  attr(hypo.all, "intersections")$'A:C:D',
  attr(hypo.all, "intersections")$'B:C:D')


# FINAL hypo and hyper methylation genes (will be feeded to AS overlap) #Converge whatever you want to this point 
# and then it should take care of itself
hyper.final <- attr(gplots::venn(list(hyper = hyper.all, hypo = hypo.all)), "intersections")$hyper
hypo.final <- attr(gplots::venn(list(hyper = hyper.all, hypo = hypo.all)), "intersections")$hypo


# Import AS
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

#le plot
AS.complete %>%
  ggplot(aes(x = fc, y = consistency, color = methylation)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("red", "blue", "grey")) +
  ggrepel::geom_text_repel(data = filter(AS.complete, methylation != "nop"), aes(x = fc, y = consistency, label = gene)) +
  scale_x_continuous(limits = c(-5,5)) +
  theme_pb() +
  theme(aspect.ratio = 1)

#Overlap with ATAC-seq data
#get ATAC-seq promoters
#load peaks into GRrange format
files <- list(ATAC_broad = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad.bed",
              ATAC_broad_Y = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad_unique_Young.bed",
              ATAC_broad_O = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad_unique_Old.bed")

peakAnnoList <- lapply(files, annotatePeak, TxDb=BiomartTxDb, tssRegion=c(-500, 500), verbose=FALSE)

extract.promoters <- function(i){
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
samples <- c("ATAC_broad", "ATAC_broad_Y", "ATAC_broad_O")
ATAC.promoter.list <- list()
for (i in 1:3){
  print(samples[i])
  ATAC.promoter.list[[samples[i]]] <- extract.promoters(i = i)
}
ATAC.promoter.list
ATAC.AS <- gplots::venn(list(common = ATAC.promoter.list[[1]],
                      young = ATAC.promoter.list[[2]],
                      old = ATAC.promoter.list[[3]],
                      AS = AS))
ATAC.AS <- list(common = attr(ATAC.AS, "intersections")$'common:AS',
                young.common = attr(ATAC.AS, "intersections")$'common:young:AS',
                old = attr(ATAC.AS, "intersections")$'old:AS',
                old.common = attr(ATAC.AS, "intersections")$'common:old:AS')

#Check overlap of hyper and hypo methylated promoters with ATAC-seq promoters
ATAC.AS.hyper <-gplots::venn(list(common = ATAC.AS[[1]],
                                  young.common = ATAC.AS[[2]],
                                  old = ATAC.AS[[3]],
                                  old.common = ATAC.AS[[4]],
                                  hyper = hyper.final))
ATAC.AS.hyper <- list(hyper_commonATAC = attr(ATAC.AS.hyper, "intersections")$'common:hyper',
                      hyper_common_old_ATAC = attr(ATAC.AS.hyper, "intersections")$'old.common:hyper',
                        hyper_common_young_ATAC = attr(ATAC.AS.hyper, "intersections")$'young.common:hyper')

ATAC.AS.hypo <-gplots::venn(list(common = ATAC.AS[[1]],
                                  young.common = ATAC.AS[[2]],
                                  old = ATAC.AS[[3]],
                                  old.common = ATAC.AS[[4]],
                                  hypo = hypo.final))
ATAC.AS.hyper <- list(hypo_commonATAC = attr(ATAC.AS.hypo, "intersections")$'common:hypo',
                      hypo_common_old_ATAC = attr(ATAC.AS.hypo, "intersections")$'old.common:hypo')
