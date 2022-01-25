# Using Ensembl-based genomic annotations for ChIPseeker
# the package works with UCSC-based annotations, so it's ready for TxDb format,
# one need to generate the complementary data for it to work with Ensembl data
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/")
library(GenomicFeatures, verbose = F)
library(ChIPseeker, verbose = F)
library(biomaRt)
library(tidyverse)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
# Load all mouse genomic data from Biomart and generate a TxDb object (it takes a while)
BiomartTxDb <- makeTxDbFromBiomart(biomart= "ensembl", dataset = "mmusculus_gene_ensembl")
#saveRDS(BiomartTxDb, file = "03_Annotation/BiomartTxDb.rds")

edb <- BiomartTxDb
#edb <- readRDS(file = "03_Annotation/BiomartTxDb.rds")

#Change chromossome names from Ensembl to UCSC (data is still intact, it's just for compability)
seqlevelsStyle(edb) <- "UCSC"


#load peaks into GRrange format
files <- list(ATAC_broad = "01_overlap_ages/01_ATAC_broad.bed",
              ATAC_broad_Y = "01_overlap_ages/01_ATAC_broad_unique_Young.bed",
              ATAC_broad_O = "01_overlap_ages/01_ATAC_broad_unique_Old.bed",
              ATAC_narrow = "01_overlap_ages/02_ATAC_narrow.bed",
              ATAC_narrow_Y = "01_overlap_ages/02_ATAC_narrow_unique_Young.bed",
              ATAC_narrow_O = "01_overlap_ages/02_ATAC_narrow_unique_Old.bed")

# peakAnnoList <- lapply(files, annotatePeak, TxDb=edb, tssRegion=c(-5000, 5000), verbose=FALSE)
# plotAnnoBar(peakAnnoList)

peakAnnoList <- lapply(files, annotatePeak, TxDb=edb, tssRegion=c(-500, 500), verbose=FALSE)
plotAnnoBar(peakAnnoList)

#Features number extraction
df <- c()
samples <- c("ATAC_broad", "ATAC_broad_Y", "ATAC_broad_O", "ATAC_narrow", "ATAC_narrow_Y", "ATAC_narrow_O")
for (i in 2:3){
  temp <- peakAnnoList[[i]]@annoStat
  sample <- rep(samples[i], times = length(temp[,2]))
  df <- rbind.data.frame(df, cbind.data.frame(temp, sample))
}

df %>%
  mutate(sample = recode(sample, ATAC_broad_Y = "Young", ATAC_broad_O = "Aged")) %>%
  mutate(sample = factor(sample, levels = c("Young", "Aged"))) %>%
  ggplot(aes(fill=Feature, y=Frequency, x = sample)) +
  geom_bar(position="stack", stat="identity") +
  labs(title = "Peak genomic annotation", x = "", y = "Frequency of annotation (%)", fill = "Anotation") +
  theme_pb() +
  theme(aspect.ratio = 1)
  
# Overlap Aging signature
AS <- read.csv("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/01_Aging_List_REAN.txt", stringsAsFactors = F, header = T, sep = "\t",) %>%
  filter(Freq_group >= 4) %>%
  filter(Gene != "")

# - select correct database from ensembl
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)

# Young-specific
df2 <- as.data.frame(peakAnnoList[[2]])
anno <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters    = 'ensembl_gene_id',
  values     = df2$geneId,
  mart       = ensembl)
colnames(anno) <- c("geneId", "Gene")
df2.anno <- merge(df2, anno, by = "geneId")
df3.young <- inner_join(AS, df2.anno, by = "Gene") %>%
  dplyr::select(Gene, annotation, Freq_group, Log2FC, seqnames, start, end, width, distanceToTSS) %>%
  mutate(annotation = gsub("\\s*\\([^\\)]+\\)","",as.character(annotation))) %>%
  filter(abs(distanceToTSS) <= 5000) %>%
  mutate(sample = "young")
#write.csv(df3, "02_Annotation/ChIPSeeker_Young_specific_annotated_peaks.csv", quote = F, row.names = F)

# Old-specific
df2 <- as.data.frame(peakAnnoList[[3]])
anno <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters    = 'ensembl_gene_id',
  values     = df2$geneId,
  mart       = ensembl)
colnames(anno) <- c("geneId", "Gene")
df2.anno <- merge(df2, anno, by = "geneId")
df3.aged <- inner_join(AS, df2.anno, by = "Gene") %>%
  dplyr::select(Gene, annotation, Freq_group, Log2FC, seqnames, start, end, width, distanceToTSS) %>%
  mutate(annotation = gsub("\\s*\\([^\\)]+\\)","",as.character(annotation))) %>%
  filter(abs(distanceToTSS) <= 5000) %>%
  mutate(sample = "aged")
#write.csv(df3, "02_Annotation/ChIPSeeker_Old_specific_annotated_peaks.csv", quote = F)

# Aging Signature vs. Annotated ATAC-seq
df3 <- rbind.data.frame(df3.young, df3.aged) %>%
  mutate(sample = factor(sample, levels = c("young", "aged"))) %>%
  mutate(annotation = factor(annotation, levels = c("Promoter", "3' UTR", "5' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic"))) %>%
  group_by(sample) %>%
  count(annotation) %>%
  ggplot(aes(x = sample, y = n, fill = annotation, label = n)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_y_continuous(limits = c(0,115)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  labs(title = "number of ATAC-seq peaks in AS genes", y = "number of peaks", x = "", fill = "Genomic region") +
  theme_pb() +
  theme(aspect.ratio = 1)
df3


#####################################################
#Extract peaks which are at promoters in samples:
# - select correct database from ensembl
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)

extract.promoters <- function(i){
  #Get gene annotations from reference annotation
  temp <- dplyr::as_tibble(peakAnnoList[[i]]@anno@elementMetadata@listData) %>%
    dplyr::select(geneChr, geneStart, geneEnd, annotation, geneId, distanceToTSS) %>%
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
  temp3 <- inner_join(temp, anno, by = "geneId") %>%
    cbind.data.frame(temp2) %>%
    dplyr::arrange(geneChr) %>%
    dplyr::mutate(geneChr = paste0("chr",geneChr)) %>%
    dplyr::select(geneChr, peakStart, peakEnd, distanceToTSS, Gene, geneId, geneStart, geneEnd)
  
  return(temp3)
}

promoters.common <- extract.promoters(i = 1)
promoters.young <- extract.promoters(i = 2)
promoters.aged <- extract.promoters(i = 3)
promoters.all <- data.frame(promoters(BiomartTxDb)@ranges)
