##### Script containing all the necessary data (unless explicity stated) to generate
##### all the figures for Figure 1 of the Epigenentics Chapter, relating to ChIP-seq data
library(tidyverse)
library(eulerr)
library(GenomicFeatures)
library(biomaRt)
library(ChIPseeker, verbose = F)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/04_Figures/Fig1_ChIP/")

# 1.F1A - Heatmap for Sara H3K4me3 and H3K27me3 samples from deeptools ----
##### code can be found at the workstation 
##### art@129.125.165.146
##### Documents/202107_Epigenetics_chapter/04_deeptools/03_TSS_enrichment/

# 2.F1B - Schematics of data analysis (illustrator only) ----
# 3.F1C - Euler diagram - Peak overlap H3K4me3 and H3K27me3 and H3K36me3 ----
# Number of overlaps
# 11187 01_H3K4me3.bed
# 1808 01_H3K4me3_unique_Old.bed
# 479 01_H3K4me3_unique_Young.bed
# 10245 02_H3K27me3.bed
# 4810 02_H3K27me3_unique_Old.bed
# 3189 02_H3K27me3_unique_Young.bed
# 19728 03_H3K36me3.bed
# 9745 03_H3K36me3_unique_Old.bed
# 22502 03_H3K36me3_unique_Young.bed

# Fits
ChIP <- euler(c(k4.y = 479, k4.a = 1808, "k4.y&k4.a" = 11187,
                k27.y = 3189, k27.a = 4810, "k27.y&k27.a" = 10245,
                k36.y = 22502, k36.a = 9745, "k36.y&k36.a" = 19728))
#Generate plot
plot(ChIP)

# 4.F1E - ChIPseeker annotation of H3K36me3 and SF1C H3K27me3 and H3K27me3 peaks ----
## Data Annotation -----
# Load BiomaRt
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
files <- list(H3K4me3 = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/02_overlap_ages/01_H3K4me3.bed",
              H3K4me3_Y = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/02_overlap_ages/01_H3K4me3_unique_Young.bed",
              H3K4me3_O = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/02_overlap_ages/01_H3K4me3_unique_Old.bed",
              H3K27me3 = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/02_overlap_ages/02_H3K27me3.bed",
              H3K27me3_Y = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/02_overlap_ages/02_H3K27me3_unique_Young.bed",
              H3K27me3_O = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/02_overlap_ages/02_H3K27me3_unique_Old.bed",
              H3K36me3 = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/02_overlap_ages/03_H3K36me3.bed",
              H3K36me3_Y = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/02_overlap_ages/03_H3K36me3_unique_Young.bed",
              H3K36me3_O = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/02_overlap_ages/03_H3K36me3_unique_Old.bed")

peakAnnoList <- lapply(files, annotatePeak, TxDb = BiomartTxDb, tssRegion=c(-500, 500), verbose = FALSE)
#plotAnnoBar(peakAnnoList)

# Features number extraction
df <- c()
samples <- c("H3K4me3", "H3K4me3_Y", "H3K4me3_O", "H3K27me", "H3K27me3_Y", "H3K27me3_O","H3K36me3", "H3K36me3_Y", "H3K36me3_O")
histone <- c(rep("H3K4me3", times =3), rep("H3K27me3", times =3), rep("H3K36me3", times =3))
group <- c(rep(c("common", "young", "aged"), times = 3))
for (i in 1:9){
  temp <- peakAnnoList[[i]]@annoStat
  age <- rep(group[i], times = length(temp[,2]))
  mark <- rep(histone[i], times = length(temp[,2]))
  df <- rbind.data.frame(df, cbind.data.frame(temp, age, mark))
}

# Genomic annotation plot for H3K4me3, H3K27me3 and H3K27me3 -----
df %>%
  mutate(age = factor(age, levels = c( "young","common", "aged"))) %>%
  mutate(mark = factor(mark, levels = c("H3K4me3", "H3K27me3", "H3K36me3"))) %>%
  #filter(age != "common") %>%
  ggplot(aes(fill=Feature, y= Frequency, x = age, label = round(Frequency, 1))) +
  geom_bar(position="stack", stat="identity") +
  labs(title = "Peak genomic annotation", x = "", y = "Frequency of annotation (%)", fill = "Anotation") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  facet_wrap(vars(mark)) +
  theme_pb() +
  theme(aspect.ratio = 1)




# 5.F1D - Bivalent peaks in ChIP-seq ----
## Import, annotate and extract promoter gene names from K4vs.K27 -----
# load peaks into GRrange format
bivalent.files <- list(bivalent.young = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/03_bivalent_peaks/01_bivalent_young.bed",
                       K4.young = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/03_bivalent_peaks/01_bivalent_young_unique_K4.bed",
                       K27.young = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/03_bivalent_peaks/01_bivalent_young_unique_K27.bed",
                       bivalent.old = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/03_bivalent_peaks/02_bivalent_old.bed",
                       K4.old = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/03_bivalent_peaks/02_bivalent_old_unique_K4.bed",
                       K27.old = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/03_bivalent_peaks/02_bivalent_old_unique_K27.bed")

# Annotate
bivalent.peakAnnoList <- lapply(bivalent.files, annotatePeak, TxDb = BiomartTxDb.custom, tssRegion=c(-500, 500), verbose = FALSE)

# Function which retrieves the true peaks from original bed files per annotation
extract.promoters <- function(i){
  #Get gene annotations from reference annotation
  temp <- dplyr::as_tibble(bivalent.peakAnnoList[[i]]@anno@elementMetadata@listData) %>%
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
  temp2 <- as.data.frame(bivalent.peakAnnoList[[i]]@anno@ranges[temp$index,]) %>%
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
samples <- names(bivalent.files)
bivalent.promoter.list <- list()
for (i in 1:6){
  print(samples[i])
  bivalent.promoter.list[[samples[i]]] <- extract.promoters(i = i)
}
## Bivalent promoter wich change with aging -----
bivalent.overlap.combinations <- gplots::venn(bivalent.promoter.list, show.plot = FALSE, intersections=TRUE)
bivalent.promoters.names <- list(attr(bivalent.overlap.combinations, "intersections")$'bivalent.young:bivalent.old', #bivalent that don't change
                                 attr(bivalent.overlap.combinations, "intersections")$'bivalent.young:K4.old', #bivalent in young which loses K27
                                 attr(bivalent.overlap.combinations, "intersections")$'bivalent.young:K27.old', #bivalent in young which loses K4
                                 attr(bivalent.overlap.combinations, "intersections")$'K4.young:bivalent.old', # K4 in young which is bivalent in old
                                 attr(bivalent.overlap.combinations, "intersections")$'K27.young:bivalent.old' # K27 in young which is bivalent in old
                                 )
str(bivalent.promoters.names)



# 6.F1F - Overlap of young/aged-specific and common promoters with Aging signature ----
AS <- read.table("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/01_Aging_List_REAN.txt",
                 header = T, stringsAsFactors = F, sep = "\t") %>%
  filter(Freq_group >= 4) %>%
  filter(Gene != "")

# Function which retrieves the true peaks from original bed files per annotation
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
samples <- c("H3K4me3", "H3K4me3_Y", "H3K4me3_O", "H3K27me3", "H3K27me3_Y", "H3K27me3_O","H3K36me3", "H3K36me3_Y", "H3K36me3_O")
promoter.list <- list()
for (i in 1:9){
  print(samples[i])
  promoter.list[[samples[i]]] <- as.vector(extract.promoters(i = i))
}

#creat df with all genes and differente hisotnes
names <- names(promoter.list)
histones.all.genes <- c()
for (name in names){
  print(name)
  Gene <- c()
  mark <- c()
  Gene <- promoter.list[[name]]$Gene
  mark <- rep(name, times = length(Gene))
  histones.all.genes <- rbind.data.frame(histones.all.genes, cbind.data.frame(Gene, mark))
}
filter(histones.all.genes, Gene == "Eed")
write.csv(file = "/Users/Art/Drive/PhD/Experiments/Epigenetics/04_Figures/Fig1_ChIP/zz_All_genes_merged_Gene_symbols.csv",
          hypo.all, row.names = F, quote = F)

# K4 overlap with AS -----
K4 <- gplots::venn(list(A =AS[,1], 
                        B = promoter.list[["H3K4me3"]][,1], 
                        C = promoter.list[["H3K4me3_Y"]][,1], 
                        D = promoter.list[["H3K4me3_O"]][,1]), show.plot = FALSE, intersections=TRUE)
str(K4)
plot(euler(c(A = 63, B = 9762, C = 284, D = 1173, E.1 = 1000,
             'A&B' = 126, 'A&C' = 4, 'A&D' = 22,'B&C' = 37, 'B&D' = 135, 'C&D' = 4,
             'A&B&C' = 3, 'A&B&D' = 3, 'B&C&D' = 1)))


#common
attr(K4, "intersections")$'A:B' %>%
  write.csv("zz_K4_common_specific.csv", quote = F, row.names = F)

#young
attr(K4, "intersections")$'A:C' %>%
  write.csv("zz_K4_young_specific.csv", quote = F, row.names = F)

#old 
attr(K4, "intersections")$'A:D' %>%
  write.csv("zz_K4_old_specific.csv", quote = F, row.names = F)


# K27 overlap with AS -----
K27 <- gplots::venn(list(A = AS[,1],
                         B = promoter.list[["H3K27me3"]][,1],
                         C = promoter.list[["H3K27me3_Y"]][,1],
                         D = promoter.list[["H3K27me3_O"]][,1]), show.plot = FALSE, intersections=TRUE)
str(K27)
plot(euler(c(A = 199, B = 2486, C = 826, D = 387, E.2 = 1000,
             'A&B' = 3, 'A&C' = 18,'B&C' = 71, 'B&D' = 14, 'C&D' = 4,
             'A&B&C' = 1, 'B&C&D' = 1)))

#common
c(attr(K27, "intersections")$'A:B',
  attr(K27, "intersections")$'A:B:C')%>%
  write.csv("zz_K27_common_specific.csv", quote = F, row.names = F)

#young
attr(K27, "intersections")$'A:C' %>%
  write.csv("zz_K27_young_specific.csv", quote = F, row.names = F)

# K36 overlap with AS -----
K36 <- gplots::venn(list(A = AS[,1],
                         B = promoter.list[["H3K36me3"]][,1],
                         C = promoter.list[["H3K36me3_Y"]][,1],
                         D = promoter.list[["H3K36me3_O"]][,1]), show.plot = FALSE, intersections=TRUE)
str(K36)
plot(euler(c(A = 124, B = 6728, C = 847, D = 730, E.3 = 1000,
             'A&B' = 66, 'A&C' = 2, 'A&D' = 21,'B&C' = 77, 'B&D' = 202, 'C&D' = 6,
             'A&B&C' = 2, 'A&B&D' = 6, 'B&C&D' = 8)))

#common
c(attr(K36, "intersections")$'A:B',
  attr(K36, "intersections")$'A:B:C',
  attr(K36, "intersections")$'A:B:D') %>%
  write.csv("zz_K36_common_specific.csv", quote = F, row.names = F)

#young
attr(K36, "intersections")$'A:C' %>%
  write.csv("zz_K36_young_specific.csv", quote = F, row.names = F)

#old 
attr(K36, "intersections")$'A:D' %>%
  write.csv("zz_K36_old_specific.csv", quote = F, row.names = F)



# Master Euler plot:
plot(euler(c(A.1 = 63, B.1 = 9762, C.1 = 284, D.1 = 284, E = 
             'A.1&B.1' = 126, 'A.1&C.1' = 4, 'A.1&D.1' = 22,'B.1&C.1' = 37, 'B.1&D.1' = 135, 'C.1&D.1' = 4,
             'A.1&B.1&C.1' = 3, 'A.1&B.1&D.1' = 3, 'B.1&C.1&D.1' = 1,
             
             A.2 = 199, B.2 = 2486, C.2 = 826, D.2 = 387,
             'A.2&B.2' = 3, 'A.2&C.2' = 18, 'A.2&D.2' = 22,'B.2&C.2' = 71, 'B.2&D.2' = 14, 'C.2&D.2' = 4,
             'A.2&B.2&C.2' = 1, 'B.2&C.2&D.2' = 1,
             
             A.3 = 124, B.3 = 6728, C.3 = 847, D.3 = 730,
             'A.3&B.3' = 66, 'A.3&C.3' = 2, 'A.3&D.3' = 21,'B.3&C.3' = 77, 'B.3&D.3' = 202, 'C.3&D.3' = 6,
             'A.3&B.3&C.3' = 2, 'A.3&B.3&D.3' = 6, 'B.3&C.3&D.3' = 8
)))





# 6.SF1A - Heatmap for Sara H3K4me3 and H3K27me3 samples from deeptools ----
##### code can be found at the workstation 
##### art@129.125.165.146
##### Documents/202107_Epigenetics_chapter/04_deeptools/03_TSS_enrichment/

# 7.SF1B - Number of called peaks Sara and Sun ChIP ----
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/")
data <- read.csv("00_ChIP_peaks_long.csv", header = T, stringsAsFactors = F) %>%
  mutate(study = recode(study, Sara = "Own", Sun = "Sun et al.;")) %>%
  mutate(mark = factor(mark, levels = c("H3K4me3", "H3K27me3", "H3K36me3"))) %>%
  mutate(age = factor(age, levels = c("young", "old"))) %>%
  ggplot() +
  geom_bar(aes(x = mark, y = peaks/1000, fill = age),
           position = "dodge", stat = "summary", fun = "mean") +
  geom_point(aes(x = mark, y = peaks/1000, group = age, shape = study, color = study), 
             position = position_dodge(width = 0.9), size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("grey", "red")) +
  scale_color_manual(values = c("#252525", "#969696")) +
  labs(x = "", y = "number of peaks (in K)", title = "ChIP-seqcalled peaks") +
  theme_pb() +
  theme(aspect.ratio = 1)
data

