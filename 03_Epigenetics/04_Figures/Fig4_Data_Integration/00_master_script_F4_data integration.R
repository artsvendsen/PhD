##### Script containing all the necessary data (unless explicity stated) to generate
##### all the figures for Figure 4 of the Epigenentics Chapter, relating to integration of ChIP, ATAC and Methylation data
library(tidyverse)
library(GenomicFeatures)
library(ChIPseeker)
library(biomaRt)
library(dendextend)
library(ggrepel)
library(viridis)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/04_Figures/Fig4_Data_Integration/")

# 1.F1A - Stacked bar plots of chromatin state of ChIP age-specific peaks ----
ChIP.open.peaks <- read.table("/Users/Art/Drive/PhD/Experiments/Epigenetics/02_overlap_ATAC_ChIP/01_chromatin_accessibiity_in_ChiP_peaks/04_number_ChIP_open_sites.txt",
                         header = F, stringsAsFactors = T) %>%
  mutate(total = read.table("/Users/Art/Drive/PhD/Experiments/Epigenetics/02_overlap_ATAC_ChIP/01_chromatin_accessibiity_in_ChiP_peaks/05_number_total_ChIP_peaks.txt",
                            header = F, stringsAsFactors = T)[,1]) %>%
  mutate(sample = c("H3K4me3_Y", "H3K4me3_O", "H3K4me3", "H3K27me_Y", "H3K27me3_O", "H3K27me3","H3K36me3_Y", "H3K36me3_O", "H3K36me3")) %>%
  mutate(open.raw = V1) %>%
  mutate(open = open.raw/total *100) %>%
  mutate(closed = 100 - open) %>%
  dplyr::select(sample, open, closed) %>%
  mutate(mark = c(rep("H3K4me3", times = 3),
                  rep("H3K27me3", times = 3),
                  rep("H3K36me3", times = 3))) %>%
  mutate(mark = factor(mark, levels = c("H3K4me3", "H3K27me3", "H3K36me3"))) %>%
  gather(state, value, -c(sample, mark)) %>%
  mutate(sample = c(rep(c("Young", "Aged", "Common"), times = 6))) %>%
  mutate(sample = factor(sample, levels = c("Young", "Common", "Aged"))) %>%
  ggplot(aes(x = sample, y = value, fill = state, label = round(value, 2)))+
  geom_bar(stat = "identity", position = "stack") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  labs(x = "", y = "% of detected peaks", title = "chromatin state of age-related promoters", fill = "Chromatin state") +
  theme_pb() +
  facet_wrap(vars(mark)) +
  theme(aspect.ratio = 0.8)
ChIP.open.peaks



# 2.F1B - Stacked bar plots of chromatin state of methylation data ----
# import ATAC and DNA data
ATAC <- list(Common = read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/03_overlap_ATAC_DNAmethylation/ATAC_promoters_unique_Gene_symbols_Common.csv", header = T, stringsAsFactors = F),
             Young = read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/03_overlap_ATAC_DNAmethylation/ATAC_promoters_unique_Gene_symbols_Young.csv", header = T, stringsAsFactors = F),
             Aged = read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/03_overlap_ATAC_DNAmethylation/ATAC_promoters_unique_Gene_symbols_Old.csv", header = T, stringsAsFactors = F))
DNA <- list(hyper = read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/03_overlap_ATAC_DNAmethylation/DNAmethylation_hypermethylated_promoter_Gene_symbols.csv", header = T, stringsAsFactors = F),
            hypo = read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/03_overlap_ATAC_DNAmethylation/DNAmethylation_hypomethylated_promoter_Gene_symbols.csv", header = T, stringsAsFactors = F))

#Count the number of DNA promoters which have open chromatin in either yong, aged or are common
df <- c()
for( j in 1:2) {
  total.genes <- length(DNA[[j]]$x)
  for (i in 1:3){
    overlap <- length(attr(gplots::venn(list(A = DNA[[j]], B = ATAC[[i]]), show.plot = FALSE, intersections=TRUE), "intersections")$'A:B')
    dna <- names(DNA[j])
    atac <- names(ATAC[i])
    open <- overlap/total.genes * 100
    df <- rbind.data.frame(df, cbind.data.frame(dna, atac, overlap, open))
  }
}

df %>%
  mutate(close = 100 - open) %>%
  dplyr::select(-overlap) %>%
  gather(state, value, -c(dna, atac)) %>%
  mutate(dna = factor(dna, levels = c("hyper", "hypo"))) %>%
  mutate(atac = factor(atac, levels = c("Young", "Common", "Aged"))) %>%
  ggplot(aes(x = atac, y = value, fill = state, label = round(value, 2)))+
  geom_bar(stat = "identity", position = "stack") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  labs(x = "", y = "% of detected peaks", title = "chromatin state of age-related promoters", fill = "Chromatin state") +
  theme_pb() +
  facet_wrap(vars(dna)) +
  theme(aspect.ratio = 0.8)


# 3.F1C - Heatmap and volcano with different epigenentic regulation clusters ----
### import all data -----
# ChIP
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/04_Figures/")
K4.young <- read.csv("Fig1_ChIP/zz_K4_young_specific.csv", header = T, stringsAsFactors = F)
K4.old <- read.csv("Fig1_ChIP/zz_K4_old_specific.csv", header = T, stringsAsFactors = F)
K4.common <- read.csv("Fig1_ChIP/zz_K4_common_specific.csv", header = T, stringsAsFactors = F)

K27.young <- read.csv("Fig1_ChIP/zz_K27_young_specific.csv", header = T, stringsAsFactors = F)
K27.common <- read.csv("Fig1_ChIP/zz_K27_common_specific.csv", header = T, stringsAsFactors = F)

K36.young <- read.csv("Fig1_ChIP/zz_K36_young_specific.csv", header = T, stringsAsFactors = F)
K36.common <- read.csv("Fig1_ChIP/zz_K36_common_specific.csv", header = T, stringsAsFactors = F)
K36.old <- read.csv("Fig1_ChIP/zz_K36_old_specific.csv", header = T, stringsAsFactors = F)

#ATAC
ATAC.common <- read.csv("Fig2_ATAC/promoter_common_ATAC.csv", header = T, stringsAsFactors = F)
ATAC.old <- read.csv("Fig2_ATAC/promoter_aged_ATAC.csv", header = T, stringsAsFactors = F)

#DNA methylation
Hypo <- read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/hypo.final.csv", header = T, stringsAsFactors = F)
Hyper <- read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/hyper.final.csv", header = T, stringsAsFactors = F)

setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/04_Figures/Fig4_Data_Integration/")
#Aging Signature
AS <- read.table("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/01_Aging_List_REAN.txt",
                 header = T, stringsAsFactors = F, sep = "\t") %>%
  filter(Freq_group >= 4) %>%
  filter(Gene != "")

### Assemble all data in to a TRUE/FALSE matrix -----
RNA <- c()
K4 <- c()
K27 <- c()
K36 <- c()
ATAC <- c()
methylation <- c()
for (gene in AS$Gene){
  #RNA
  value <- c()
  young <- ifelse(filter(AS, Gene == gene)[,3] < 0, TRUE, FALSE)
  old <- ifelse(filter(AS, Gene == gene)[,3] > 0, TRUE, FALSE)
  value <- c(young, old)
  
  RNA <- rbind.data.frame(RNA, value)
  #Get presence in ChIP data
  # K4
  value <- c()
  young <- gene %in% K4.young$x
  common <- gene %in% K4.common$x
  old <- gene %in% K4.old$x
  value <- c(young, common, old)
  
  K4 <- rbind.data.frame(K4, value)
  
  # K27
  value <- c()
  young <- gene %in% K27.young$x
  common <- gene %in% K27.common$x
  value <- c(young, common)
  
  K27 <- rbind.data.frame(K27, value)
  
  # K36
  value <- c()
  young <- gene %in% K36.young$x
  common <- gene %in% K36.common$x
  old <- gene %in% K36.old$x
  value <- c(young, common, old)
  
  K36 <- rbind.data.frame(K36, value)
  
  #ATAC
  value <- c()
  common <- gene %in% ATAC.common$x
  old <- gene %in% ATAC.old$x
  value <- c(common, old)
  
  ATAC <- rbind.data.frame(ATAC, value)
  
  #Methylation
  value <- c()
  hypo <- gene %in% Hypo$x
  hyper <- gene %in% Hyper$x
  value <- c(hypo, hyper)
  
  methylation <- rbind.data.frame(methylation, value)
}


df <- cbind.data.frame(RNA,
                       K4,
                       K27,
                       K36,
                       ATAC,
                       methylation) %>%
  as.matrix()
colnames(df) <- c("RNA_DO","RNA_UP",
                  "K4_young", "K4_common", "K4_old",
                  "K27_young", "K27_common",
                  "K36_young", "K36_common", "K36_old",
                  "ATAC_common", "ATAC_old",
                  "DNA_hypo", "DNA_hyper")
rownames(df) <- AS$Gene
df <- df[,-c(4,7,9,11)]
head(df)
remove(K4.young, K4.old, K4.common, K27.common, K27.young, K36.old, K36.common, K36.young, ATAC.common, ATAC.old,Hyper, Hypo,
       K4,K27, K36, ATAC, methylation)

### Determine how many clusters there will be -----
# matrix with conversted values from TRUE/FALSE to 1/0
dat <- 1*df

out <- pheatmap::pheatmap(dat,
                          #cluster_cols = F,
                          border_color=FALSE,
                          silent = T)
elbow <- c()
for(h in seq(0,3, by =0.1)){
  #print(h)
  num.clusters <- max(cutree(out$tree_row, h=h))
  elbow <- rbind.data.frame(elbow, cbind.data.frame(h, num.clusters))
}
ggplot(elbow, aes(y = h, x = num.clusters)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 25) +
  #geom_hline(yintercept = 0.9) +
  theme(aspect.ratio = 1)

#Determine the clustering
height <- 0.9
plot(out$tree_row)
abline(h=height, col="red", lty=2, lwd=2)
max(cutree(out$tree_row, h=height))


### Heatmap with epigenetic clusters -----
#Add label for clusters (genes/fows)
my_gene_col <- data.frame(cluster = cutree(out$tree_row, h=height))
#Add label for different epigenetic levels (cols)
my_sample_col <- data.frame(level = colnames(dat))
row.names(my_sample_col) <- colnames(dat)

#Colors
my_color <- list(
  # level = c("#f768a1", #ATAC
  #           "#9C1B45", "#3287BD", #DNA
  #           "#FF8000", #K27
  #           "#810f7c", "#8856a7", #K36
  #           "#006d2", "#2ca25f", #K4
  #           "#ED3224", "#3D58A6" #RNA
  #           ),
  cluster = viridis(25)
)

out2 <- pheatmap::pheatmap(dat,
                           border_color=FALSE,
                           color = c("black", "#BDBDBD"),
                           cutree_rows = max(cutree(out$tree_row, h=height)),
                           annotation_col = my_sample_col, 
                           annotation_row = my_gene_col,
                           annotation_colors = my_color,
                           show_rownames = F,
                           cellheight=1, cellwidth = 40
)

### Volcano with epigenetic clusters -----
#Get genes and clusters
clusters <- data.frame(Gene = attr(sort(cutree(out$tree_row, h=height)), "names"),
                       cluster = sort(cutree(out$tree_row, h=height))) %>%
  merge(x = AS, y = ., by = "Gene") %>%
  arrange(-Freq_group, -Log2FC) %>%
  dplyr::select(Gene, cluster)

clusters %>%
  group_by(cluster) %>%
  count %>%
  ggplot(aes(x = factor(cluster),
             y = n,
             fill = factor(cluster),
             label = n)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  geom_text(size = 3, position = position_stack(vjust = 1)) +
  labs(x = "epigenetic group", y = "number of genes") +
  theme_pb() +
  theme(aspect.ratio = 1,
        legend.position = "none")
  
  
  

AS.complete <- as.data.frame(df) %>%
  mutate(Gene = rownames(df)) %>%
  merge(x = ., y = clusters, by = "Gene") %>%
  arrange(cluster)

selected.clusters <- c(1:3,5:9,11,12,17,19,20,24,26) #see 01_cluters_information.xlsx for selection
AS.epi.aging <- AS.complete %>%
  filter(cluster %in% selected.clusters) %>%
  dplyr::select(Gene, cluster) %>%
  merge(x = AS, y = ., by = "Gene")
  
# plot volcano
ggplot(AS.epi.aging, aes(x = Log2FC, y = Freq_group, color = factor(cluster))) +
  geom_point(data = AS, aes(x = Log2FC, y = Freq_group), color = "grey") +
  geom_vline(xintercept = c(0),  color = "grey") +
  geom_hline(yintercept = 4, color = "grey") +
  geom_point(alpha = 0.8) +
  scale_x_continuous(limits = c(-5,5)) +
  labs(x = "Fold-Change", y = "AS consistency score", color = "Epigenetic cluster") +
  theme_pb() +
  theme(aspect.ratio = 1)

### Generate tables with summary of all epigenetic information -----
#### Import data and generate tables -----
# ChIP
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/04_Figures/")
K4.young <- read.csv("Fig1_ChIP/zz_K4_young_specific.csv", header = T, stringsAsFactors = F)
K4.old <- read.csv("Fig1_ChIP/zz_K4_old_specific.csv", header = T, stringsAsFactors = F)
K4.common <- read.csv("Fig1_ChIP/zz_K4_common_specific.csv", header = T, stringsAsFactors = F)

K27.young <- read.csv("Fig1_ChIP/zz_K27_young_specific.csv", header = T, stringsAsFactors = F)
K27.common <- read.csv("Fig1_ChIP/zz_K27_common_specific.csv", header = T, stringsAsFactors = F)

K36.young <- read.csv("Fig1_ChIP/zz_K36_young_specific.csv", header = T, stringsAsFactors = F)
K36.common <- read.csv("Fig1_ChIP/zz_K36_common_specific.csv", header = T, stringsAsFactors = F)
K36.old <- read.csv("Fig1_ChIP/zz_K36_old_specific.csv", header = T, stringsAsFactors = F)

#ATAC
ATAC.common <- read.csv("Fig2_ATAC/promoter_common_ATAC.csv", header = T, stringsAsFactors = F)
ATAC.old <- read.csv("Fig2_ATAC/promoter_aged_ATAC.csv", header = T, stringsAsFactors = F)

#DNA methylation
Hypo <- read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/hypo.final.csv", header = T, stringsAsFactors = F)
Hyper <- read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/hyper.final.csv", header = T, stringsAsFactors = F)

#Aging Signature
AS <- read.table("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/01_Aging_List_REAN.txt",
                 header = T, stringsAsFactors = F, sep = "\t") %>%
  filter(Freq_group >= 4) %>%
  filter(Gene != "")
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/04_Figures/Fig4_Data_Integration/")
K4 <- c()
K27 <- c()
K36 <- c()
ATAC <- c()
methylation <- c()

for (gene in AS$Gene){
  #Get presence in ChIP data
  # K4
  young <- filter(K4.young, x == gene)
  common <- filter(K4.common, x == gene)
  old <- filter(K4.old, x == gene)
  value <- c()
  value <- ifelse(nrow(common) >= 1,
                  "common",
                  ifelse(nrow(young) >= 1,
                         "young",
                         ifelse(nrow(old) >= 1,
                                "old",
                                "nd")
                  )
  )
  K4 <- rbind.data.frame(K4, value)
  
  # K27
  young <- filter(K27.young, x == gene)
  common <- filter(K27.common, x == gene)
  value <- c()
  value <- ifelse(nrow(common) >= 1,
                  "common",
                  ifelse(nrow(young) >= 1,
                         "young",
                         "nd")
  )
  
  K27 <- rbind.data.frame(K27, value)
  
  # K36
  young <- filter(K36.young, x == gene)
  common <- filter(K36.common, x == gene)
  old <- filter(K36.old, x == gene)
  value <- c()
  value <- ifelse(nrow(common) >= 1,
                  "common",
                  ifelse(nrow(young) >= 1,
                         "young",
                         ifelse(nrow(old) >= 1,
                                "old",
                                "nd")
                  )
  )
  K36 <- rbind.data.frame(K36, value)
  
  #ATAC
  common <- filter(ATAC.common, x == gene)
  old <- filter(ATAC.old, x == gene)
  value <- c()
  value <- ifelse(nrow(common) >= 1,
                  "common",
                  ifelse(nrow(young) >= 1,
                         "old",
                         "nd")
  )
  ATAC <- rbind.data.frame(ATAC, value)
  
  #Methylation
  hypo <- filter(Hypo, x == gene)
  hyper <- filter(Hyper, x == gene)
  value <- c()
  value <- ifelse(nrow(hypo) >= 1,
                  "hypomethylated",
                  ifelse(nrow(hyper) >= 1,
                         "hypermethylated",
                         "nd")
  )
  methylation <- rbind.data.frame(methylation, value)
}

#Final table with all AS genes and their epigenentic information 
AS.epi <- AS %>%
  mutate(accessibility = ATAC[,1]) %>%
  mutate(H3K4m3 = K4[,1]) %>%
  mutate(H3K27m3 = K27[,1]) %>%
  mutate(H3K36m3 = K36[,1]) %>%
  mutate(DNA = methylation[,1]) %>%
  merge(x = ., y = clusters, by = "Gene") %>%
  arrange(-Freq_group, Log2FC) %>%
  as_tibble()

remove(K4.young, K4.old, K4.common, K27.common, K27.young, K36.old, K36.common, K36.young, ATAC.common, ATAC.old,Hyper, Hypo,
       K4,K27, K36, ATAC, methylation)

#write.csv(file = "04_Table_all_AS_genes_epigenetics.csv", AS.epi, row.names = F, quote = F)


#Final table with just genes with epigenetic regulation (up and down)
blacklist <- data.frame(Gene = c(filter(clusters, cluster == 4)$Gene,
               filter(clusters, cluster == 10)$Gene))

AS.epi.currated <- AS.epi %>%
  anti_join(x = ., y = blacklist, by = "Gene")

#write.csv(file = "04_Table_AS_genes_with_epigenetic_regulation.csv", AS.epi, row.names = F, quote = F)



# 3.F1D - Linear regression for all epigenentic data ----
## Linear regression -----
mytab<-AS.epi.currated
mytab[mytab=="nd"]<-"absent"
mylist<-colnames(mytab)

#lm() function consumes table easily, no adjustments required
fit=lm(Log2FC~0 + accessibility+ H3K4m3+ H3K27m3+
         H3K36m3+ DNA,mytab)
anova(fit)
summary(fit)

lm.fit <- data.frame(x = AS.epi.currated$Log2FC,
                     y = predict(fit))
## Linear regression plot -----
ggplot(lm.fit, aes(x,y)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_smooth(method="lm", formula= (y ~ x), alpha = 0.5) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(limits = c(-2.75,5)) +
  geom_label(
    aes(x = -2, y = 3), hjust = 0,
    label = paste("p-val < 2.2e-16",
                  "\nAdj R2 = 0.7504")) +
  labs(x = "mean log2FC", y = "Epigenetic regulation (lm fit)") +
  theme_pb() +
  theme(aspect.ratio = 1)


# 4.F1E - Euler diagram of epigenentic regulation of upregulated genes in aging ----
#Import data
# ChIP
K4 <- read.csv("Fig1_ChIP/zz_K4_old_specific.csv", header = T, stringsAsFactors = F)
K27 <- read.csv("Fig1_ChIP/zz_K27_young_specific.csv", header = T, stringsAsFactors = F)
K36 <- read.csv("Fig1_ChIP/zz_K36_old_specific.csv", header = T, stringsAsFactors = F)

#ATAC
ATAC.common <- read.csv("Fig2_ATAC/promoter_common_ATAC.csv", header = T, stringsAsFactors = F)
ATAC.old <- read.csv("Fig2_ATAC/promoter_aged_ATAC.csv", header = T, stringsAsFactors = F)
ATAC.common <- rbind(ATAC.common, ATAC.old)
#hypomethylated
hypo <- read.csv("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/hypo.final.csv", header = T, stringsAsFactors = F)

ATAC.num <- unique(ATAC.common)

#Filter ATAC open sites
K4.filtered  <- attr(gplots::venn(list(A = ATAC.common,
                                       B = K4)), "intersections")$'A:B'

K36.filtered  <- attr(gplots::venn(list(A = ATAC.common,
                                        B = K36)), "intersections")$'A:B'

hypo.filtered <- attr(gplots::venn(list(A = ATAC.common,
                                        B = hypo)), "intersections")$'A:B'

#Overlap
myV2 <- plotVenn(list(
  K4 = K4.filtered,
  K27 = K27[,1],
  K36 = K36.filtered,
  hypo = hypo.filtered))

myV2 <- plotVenn(nVennObj = myV2)
showSVG(nVennObj = myV2,
        opacity = 0.1,
        borderWidth = 3,
        labelRegions = F,
        setColors = c("#008837",
                      "#ff8000",
                      "#82027e",
                      "#3287bd")
)



read.csv("Fig4_Data_Integration/04_Table_AS_genes_with_epigenetic_regulation.csv",
         header = T, stringsAsFactors = T) %>%
  as_tibble() %>%
  dplyr::select(-c(Freq_group, Log2FC, accessibility)) %>%
  filter(DNA == "hypomethylated" | DNA == "hypermethylated") %>%
  gather(mark, value, -c(Gene,DNA)) %>%
  filter(value != "nd") %>%
  mutate(histone = paste0(mark,"_",value)) %>%
  mutate(histone = factor(histone, levels = c("H3K4m3_common", "H3K4m3_old", "H3K27m3_young", "H3K36m3_common", "H3K36m3_old"))) %>%
  mutate(mark = factor(mark, levels = c("H3K4m3", "H3K27m3", "H3K36m3"))) %>%
  dplyr::select(-c(value)) %>%
  group_by(DNA, histone) %>%
  count



ggplot(aes(x=histone, fill = DNA)) +
  geom_bar(position="stack") +
  scale_fill_manual(values = c("#9e0142", "#3288bd")) +
  labs(title = "ChIP vs. DNA methylation", x = "", y = "number of genes", fill = "DNA methylation") +
  #geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  facet_wrap(vars(mark),scales = "free") +
  theme_pb()
#theme(aspect.ratio = 1)