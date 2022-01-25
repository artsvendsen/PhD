##### Script containing all the necessary data (unless explicity stated) to generate
##### all the figures for Figure 2 of the Epigenentics Chapter, relating to ATAC data
library(tidyverse)
library(eulerr)
library(GenomicFeatures)
library(biomaRt)
library(ChIPseeker, verbose = F)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/04_Figures/Fig2_ATAC/")

# 1.F1A - Euler diagram with young and aged ATAC-seq peak overlap ----
#Overall overlap
# 116601 /Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad.bed
# 46294 /Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad_unique_Old.bed
# 38543 /Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad_unique_Young.bed

con <- c(young = 38543,
         aged = 46294,
         "young&aged" = 116601)
plot(euler(con), 
     labels = c("young", "aged"),
     fills = c("grey", "red"),
     line = 1,
     quantities = list(fontsize = 12))


# 2.F1B - number of genomic annotations  ----
## Import and data crunking ------
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
files <- list(ATAC_broad = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad.bed",
              ATAC_broad_Y = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad_unique_Young.bed",
              ATAC_broad_O = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad_unique_Old.bed")

peakAnnoList <- lapply(files, annotatePeak, TxDb=BiomartTxDb, tssRegion=c(-500, 500), verbose=FALSE)

## Genomic annotation for young, aged and common peaks ----
#Features number extraction
df <- c()
samples <- c("ATAC_broad", "ATAC_broad_Y", "ATAC_broad_O")
for (i in 1:3){
  temp <- peakAnnoList[[i]]@annoStat
  sample <- rep(samples[i], times = length(temp[,2]))
  df <- rbind.data.frame(df, cbind.data.frame(temp, sample))
}

df %>%
  mutate(sample = recode(sample, ATAC_broad_Y = "Young", ATAC_broad_O = "Aged", "ATAC_broad" = "Common")) %>%
  mutate(sample = factor(sample, levels = c("Young", "Common", "Aged"))) %>%
  ggplot(aes(fill=Feature, y=Frequency, x = sample, label = round(Frequency,2))) +
  geom_bar(position="stack", stat="identity") +
  labs(title = "Peak genomic annotation", x = "", y = "Frequency of annotation (%)", fill = "Anotation") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_pb() +
  theme(aspect.ratio = 1)


# 3.F1C/SF4 - Diffbind volcano ----
diffbind.vulc <- read.csv("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/07_Transcriptional_activation/ATAC-seq/01_diffbind/02_ATACseq_narrow_p0_01_DBP_diffbind.csv") %>%
  ggplot(aes(x = Fold, y = -log10(FDR), color = ifelse(Fold > 0, 'red', "grey"))) +
  geom_point(size = 7,shape = 21, alpha = 0.5) +
  scale_color_identity() +
  geom_vline(xintercept = c(-1,1),  color = "grey") +
  geom_hline(yintercept = 1, color = "grey") +
  labs(title = "Young vs. Aged accessible sites", y = "-log10(FDR)")+
  xlab(expression(atop("accessible sites",young %<->% aged))) +
  scale_y_continuous(limits = c(0,12)) +
  scale_x_continuous(limits = c(-2.5,2.5)) +
  theme_pb() +
  theme(aspect.ratio = 1,
        legend.position = "none")
diffbind.vulc

diffbind.DE <- read.csv("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/07_Transcriptional_activation/ATAC-seq/01_diffbind/02_ATACseq_narrow_p0_01_DBP_diffbind.csv") %>%
  mutate(age = ifelse(Fold > 0, "aged", "young")) %>%
  mutate(age = factor(age, levels = c("young", "aged"))) %>%
  group_by(age) %>%
  count() %>%
  ggplot(aes(x = age, y = n, fill = age)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label= n), vjust=-0) +
  scale_y_log10(limits = c(1,10000)) +
  scale_fill_manual(values = c("grey", "red")) +
  labs(title = "Young vs. Aged accessible sites", x = "", y = "number of DBPs") +
  theme_pb() +
  theme(aspect.ratio = 1)
diffbind.DE

# 4.F1D - ATAC-seq overlap with hematopoeitic enhancer overlap ----
# Done in prism

# 5.F1E - Age-specific Motif erichment ----
chrom <- read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/05_chromVAR/01_chromVAR_variability_out_narrow_merged_samples_DIRECTION.csv",
                  header = T, stringsAsFactors = F) %>%
  as_tibble() %>%
  arrange(p_value_adj) %>%
  dplyr::select(-c(X, old, bootstrap_lower_bound, bootstrap_upper_bound)) %>%
  mutate(name = str_remove(string = name, pattern = fixed("(var.2)"))) %>%
  dplyr::rename(old_zscore = `X2017_07_02_06_old_merged_rmdup.bam`) %>%
  dplyr::rename(young_zscore = `X2017_07_04_08_young_merged_rmdup.bam`) %>%
  assign("df", value = ., pos = 1) %>%
  ggplot(aes(x = old_zscore, y = -log10(p_value_adj), color = ifelse(direction == "old", 'red', 'grey'))) +
  geom_vline(xintercept = 0,  color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_point(size = 5, alpha = 0.8, shape = 21) +
  ylim(0, 200) +
  scale_color_identity() +
  geom_text_repel(data = head(filter(df, direction == "old"), 15), aes(x = old_zscore, y = -log10(p_value_adj), label = name), box.padding = 0.5) +
  geom_text_repel(data = head(filter(df, direction == "young"), 4), aes(x = old_zscore, y = -log10(p_value_adj), label = name), box.padding = 1) +
  labs(title = "ATAC-seq Motif enrichment analysis", x = "accessibility") +
  theme_pb() +
  theme(aspect.ratio = 1)
chrom

# 6.F1F - TF search and schematics ----
# Done in illustrator directly (no data analysis)

# 7.F1G/H - Aging Signature gene overlap with TF targets ----
#Aging Signature genes
AS <- read.table("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/01_Aging_List_REAN.txt",
                 header = T, stringsAsFactors = F, sep = "\t") %>%
  filter(Freq_group >= 4) %>%
  filter(Gene != "") %>%
  dplyr::select(Gene)

#Motif enriched genes
motif.enriched.genes <- read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/05_chromVAR/01_chromVAR_variability_out_broad_merged_samples.csv",
                                 header = T, stringsAsFactors = F) %>%
  as_tibble() %>%
  #dplyr::select(-c(X, old, bootstrap_lower_bound, bootstrap_upper_bound)) %>%
  mutate(Gene = str_remove(string = name, pattern = fixed("(var.2)"))) %>%
  dplyr::select(Gene)
motif.enriched.genes <- gsub(" ", "", unlist(strsplit(motif.enriched.genes$Gene, "::"))) %>%
  str_to_title()

## Overlap AS vs. Genes with enriched motifs ----
ASvsMotif <- gplots::venn(list(AS, motif.enriched.genes), show.plot = FALSE, intersections=TRUE)


## TF downstream targets
files <- list.files("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/05_chromVAR/03_Downstream_targets",
                      pattern = "_targets.mouse.tsv")
targets <-list()
for (file in files){
  target <- str_remove(file, pattern = "_targets.mouse.tsv")
  targets[[target]] <- read.table(paste0("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/05_chromVAR/03_Downstream_targets/",
                                         file), stringsAsFactors = F,header = F)[,2]
}

overlap.num <- c()
for (i in 1:4){
  overlap <- gplots::venn(list(AS, targets[[i]]), show.plot = FALSE, intersections=TRUE)
  overlap.num <- rbind.data.frame(overlap.num, lengths(attributes(overlap)$intersections))
}
overlap.num$tagets <- str_remove(files, pattern = "_targets.mouse.tsv")
colnames(overlap.num) <- c("AS", "all.targets", "overlap", "target")

## Overlap AS vs. downstream Genes with enriched motifs ----
plot(euler(c(AS = 213, Bcl6 = 16, Hnf4a = 40, Id2 = 16, Jun = 16,
             "AS&Bcl6" = 3, 
             "AS&Hnf4a" = 1,
             "AS&Jun" = 4
             )))

# 8.SF1 - Number of called peaks----
# 322627 /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_02_ATAC_5_3_old_CGAGGCTG_p0.01_peaks.broadPeak
# 242700 /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_04_ATAC_5_3_young_GCTACGCT_p0.01_peaks.broadPeak
# 468792 /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_06_ATAC_5_4_old_GTAGAGGA_p0.01_peaks.broadPeak
# 461860 /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_08_ATAC_5_4_young_AAGAGGCA_p0.01_peaks.broadPeak

df.peaks <- data.frame(sample = c("young", "young", "aged", "aged"),
                       rep = c("rep1", "rep2", "rep1", "rep2"),
                       peaks = c(242700, 461860, 322627 , 468792)) %>%
  mutate(sample = factor(sample, levels = c("young", "aged"))) %>%
  mutate(peaks = peaks/1000) %>%
  ggplot(aes(fill = sample)) +
  geom_bar(aes(x = sample, y = peaks), stat = "summary", fun = "mean") +
  geom_point(aes(x = sample, y = peaks)) +
  scale_y_continuous(limits = c(0, 500)) +
  scale_fill_manual(values = c("grey", "red")) +
  labs(title = "number of ATAC-seq peaks", x = "", y = "number of peaks (in K)") +
  theme_pb() +
  theme(aspect.ratio = 1,
        legend.position = "none")
df.peaks
# 9.SF2 - Pearson correlation of ATAC-seq samples ----
# done in deeptools
# 10.SF3 - Peak widths ----
## Promoters -----
#Get promoter widths
files <- list(
  young.specific.rep1 = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width/03_Young_specific_original_rep1.bed",
  young.specific.rep2 = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width/03_Young_specific_original_rep2.bed",
  old.specific.rep1 = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width/03_Old_specific_original_rep1.bed",
  old.specific.rep2 = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width/03_Old_specific_original_rep2.bed",
  young.common.rep1 = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width/03_Young_common_original_rep1.bed",
  young.common.rep2 = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width/03_Young_common_original_rep2.bed",
  old.common.rep1 = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width/03_Old_common_original_rep1.bed",
  old.common.rep2 = "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width/03_Old_common_original_rep2.bed")

peakAnnoList <- lapply(files, annotatePeak, TxDb = BiomartTxDb, tssRegion=c(-500, 500), verbose = FALSE)

extract.promoter.width <- function(i, segment, age, rep){
  #Get gene annotations from reference annotation
  temp <- dplyr::as_tibble(peakAnnoList[[i]]@anno@elementMetadata@listData) %>%
    dplyr::mutate(index = 1:length(geneChr)) %>%
    dplyr::filter(annotation == "Promoter")
  
  #Collect real peaks back (ones that were in orgininal bed files, not gene annotations)
  temp2 <- as.data.frame(peakAnnoList[[i]]@anno@ranges[temp$index,]) %>%
    dplyr::mutate(age = age) %>%
    dplyr::mutate(segment = segment) %>%
    dplyr::mutate(rep = rep) %>%
    dplyr::select(width, segment, age, rep)
  
  return(temp2)
}
ref <- data.frame(segment = c(rep("specific", times = 4), rep("common", times = 4)),
                  age = c(rep(c("young", "young", "aged", "aged"), times = 2)),
                  rep = c(rep(c(1,2), times = 4))
                  )
width.promoters <- c()
for (i in 1:8){
  print(i)
  temp <- c()
  temp <- extract.promoter.width(i = i, segment = ref[i,1], age = ref[i,2], rep = ref[i,3])
  width.promoters <- rbind.data.frame(width.promoters, temp)
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})
geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

width.promoters %>%
  mutate(segment = factor(segment, levels = c("common", "specific"))) %>%
  mutate(age = factor(age, levels = c("young", "aged"))) %>%
  mutate(rep = factor(rep, levels = c("1", "2"))) %>%
  ggplot(aes(x = segment, y = width, fill = age)) +
  geom_split_violin(trim = TRUE) +
  scale_fill_manual(values = c("grey", "red")) +
  geom_boxplot(width = 0.2, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0, aes(fill = age)) +
  scale_y_continuous(limits = c(0, 3000)) +
  labs(title = "Promoter Peak width", x = "", y = "Peak width (bp)") +
  #facet_wrap(vars(rep)) +
  theme_pb() +
  theme(aspect.ratio = 1)

summary(filter(width.promoters, segment == "common" & age == "young")$width)
summary(filter(width.promoters, segment == "common" & age == "aged")$width)

## Enhancers -----
# All peaks per sample
files <- list.files("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width/enhancers/", pattern = "01_*")
ref <- data.frame(segment = rep(c(rep("common", times = 2), rep("specific", times = 2)), times = 2),
                  age = c(rep("aged", times = 4), rep("young", times = 4)),
                  rep = c(rep(c(1,2), times = 4)))

width.enhancers <- c()
for (i in 1:8){
  print(i)
  width <- c()
  width <- read.table(paste0("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width/enhancers/",
                             files[i]), header = F, stringsAsFactors = F, sep = "\t")[,c(2,3)] %>%
    mutate(width = V3-V2)
  temp <- cbind.data.frame(width = width[,3], segment = ref[i,1], age = ref[i,2], rep = ref[i,3])
  width.enhancers <- rbind.data.frame(width.enhancers, temp)
}

width.enhancers %>%
  mutate(segment = factor(segment, levels = c("common", "specific"))) %>%
  mutate(age = factor(age, levels = c("young", "aged"))) %>%
  mutate(rep = factor(rep, levels = c("1", "2"))) %>%
  ggplot(aes(x = segment, y = width, fill = age)) +
  geom_split_violin(trim = TRUE) +
  scale_fill_manual(values = c("grey", "red")) +
  geom_boxplot(width = 0.2, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0, aes(fill = age)) +
  #scale_y_log10() +
  labs(title = "Enhancers Peak width", x = "", y = "Peak width (bp)") +
  #facet_wrap(vars(rep)) +
  theme_pb() +
  theme(aspect.ratio = 1)

summary(filter(width.enhancers, segment == "specific" & age == "young")$width)
summary(filter(width.enhancers, segment == "specific" & age == "aged")$width)

summary(filter(width.enhancers, segment == "common" & age == "young")$width)
summary(filter(width.enhancers, segment == "common" & age == "aged")$width)




# 10.Collect ATAC promoter names ----
# extract promoters
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
samples <- c("ATAC_broad", "ATAC_broad_Y", "ATAC_broad_O")
promoter.list <- list()
for (i in 1:3){
  print(samples[i])
  promoter.list[[samples[i]]] <- extract.promoter(i = i)
}

#str(promoter.list)
#temp <- unique(promoter.list[[3]])[,1]
# write.csv(file = "/Users/Art/Drive/PhD/Experiments/Epigenetics/03_overlap_ATAC_DNAmethylation/ATAC_promoters_unique_Gene_symbols_Old.csv",
#           temp, row.names = F, quote = F)


ATAC.AS <- gplots::venn(list(common = promoter.list[[1]], young = promoter.list[[2]], old = promoter.list[[3]], AS = AS))
str(ATAC.AS)

common.ATAC.promoters <- unique(c(attr(ATAC.AS, "intersections")$'common:old:AS',
                                  attr(ATAC.AS, "intersections")$'common:young:AS',  
                                  attr(ATAC.AS, "intersections")$'common:AS',
                                  attr(ATAC.AS, "intersections")$'common:young:old:AS'))
#write.csv(common.ATAC.promoters, "promoter_common_ATAC.csv", quote = F, row.names = F, col.names = F)

aged.ATAC.promomter <- attr(ATAC.AS, "intersections")$'old:AS'
#write.csv(aged.ATAC.promomter, "promoter_aged_ATAC.csv", quote = F, row.names = F, col.names = F)
