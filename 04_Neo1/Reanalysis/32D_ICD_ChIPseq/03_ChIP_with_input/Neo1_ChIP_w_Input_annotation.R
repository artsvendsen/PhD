library(tidyverse)
library(eulerr)
library(GenomicFeatures)
library(biomaRt)
library(ChIPseeker, verbose = F)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
setwd("/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/03_ChIP_with_input/")

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
files <- list(#Rep1_narrow = "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/03_ChIP_with_input/ChIP_Neo1_wcontrol_rep1_peaks.narrowPeak",
              #Rep2_narrow = "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/03_ChIP_with_input/ChIP_Neo1_wcontrol_rep2_peaks.narrowPeak",
              #Rep3_narrow = "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/03_ChIP_with_input/ChIP_Neo1_wcontrol_rep3_peaks.narrowPeak",
              Merge_narrow = "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/03_ChIP_with_input/ChIP_Neo1_merged_wcontrol_peaks.narrowPeak",
              #Rep1_broad = "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/03_ChIP_with_input/ChIP_Neo1_wcontrol_rep1_peaks.broadPeak",
              #Rep2_broad = "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/03_ChIP_with_input/ChIP_Neo1_wcontrol_rep2_peaks.broadPeak",
              #Rep3_broad = "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/03_ChIP_with_input/ChIP_Neo1_wcontrol_rep3_peaks.broadPeak",
              Merge_broad = "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/03_ChIP_with_input/ChIP_Neo1_merged_wcontrol_peaks.broadPeak")


peakAnnoList <- lapply(files, annotatePeak, TxDb=BiomartTxDb, tssRegion=c(-500, 500), verbose=FALSE)

#Features number extraction
df <- c()
samples <- c(#"Rep1_narrow", "Rep2_narrow", "Rep3_narrow", 
             "Merge_narrow", 
             #"Rep1_broad", "Rep2_broad", "Rep3_broad", 
             "Merge_broad")
for (i in 1:2){
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

#Get promots
extract.promoters <- function(i){
  #Get gene annotations from reference annotation
  temp <- dplyr::as_tibble(peakAnnoList[[i]]@anno@elementMetadata@listData) %>%
    dplyr::select(geneChr, geneStart, geneEnd, annotation, geneId, distanceToTSS) %>%
    dplyr::mutate(geneId = substr(geneId, start = 1, stop = 18)) %>%
    dplyr::mutate(index = 1:length(geneChr)) %>%
    dplyr::filter(distanceToTSS == 0) %>%
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
#samples <- names(bivalent.files)
promoter.list <- list()
for (i in 1:2){
  print(samples[i])
  promoter.list[[samples[i]]] <- extract.promoters(i = i)
}
str(promoter.list)

overlap <- unique(attr(gplots::venn(list(broad = promoter.list[[1]], narrow = promoter.list[[2]])), "intersections")$'broad:narrow')
narrow <- unique(attr(gplots::venn(list(Neo1 = promoter.list[[3]], Input = promoter.list[[4]])), "intersections")$Neo1)

clipr::write_clip(overlap)

write.table(x = promoter.list[[samples[i]]], file = "Gene_names_overlap_Broad_narrow_Promoter_ChIP_Neo1_w_Input.csv", sep = "\t", quote = F, col.names = F)



clipr::write_clip(promoter.list[[samples[3]]]$Gene)
