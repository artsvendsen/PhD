##### Extract peak coordinates from H3K4me3, H3K27me3 and H3K36me3 from core promoter regions
##### A. Svendsen Oct 2021
library(tidyverse)
library(GenomicFeatures)
library(biomaRt)
library(ChIPseeker, verbose = F)
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/05_ChIP_promoters/")

# Load BiomaRt
ensembl <- useMart("ensembl")
#if unresponsive
#ensembl <- useMart("ensembl", host="uswest.ensembl.org")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)

## Build a TxDb object from the actual annotation that is used to align samples (running locally):
BiomartTxDb.custom <- makeTxDbFromGFF(file = "/Users/Art/Drive/PhD/Scripts/GENCODE/gencode.vM13.annotation.gff3",
                                      format = "gff3",
                                      organism = "Mus musculus")

#Change chromossome names from Ensembl to UCSC (data is still intact, it's just for compability)
seqlevelsStyle(BiomartTxDb.custom) <- "UCSC"

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

peakAnnoList <- lapply(files, annotatePeak, TxDb = BiomartTxDb.custom, tssRegion=c(-500, 500), verbose = FALSE)

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
    # dplyr::select(Gene) %>%
    # unique()
  
  return(temp3)
}
samples <- c("H3K4me3", "H3K4me3_Y", "H3K4me3_O", "H3K27me", "H3K27me3_Y", "H3K27me3_O","H3K36me3", "H3K36me3_Y", "H3K36me3_O")
promoter.list <- list()
for (i in 1:9){
  print(samples[i])
  promoter.list[[samples[i]]] <- extract.promoters(i = i)
}

numbers <- c(rep("01_", times = 3), rep("02_", times = 3), rep("03_", times = 3))
for (i in 1:9){
  print(paste0(numbers[i], samples[i], "_core_promoters.bed"))
  write.table(promoter.list[[i]],
              paste0(numbers[i], samples[i], "_core_promoters.bed"),
              sep = "\t", row.names = F, col.names = F, quote = F)
}
