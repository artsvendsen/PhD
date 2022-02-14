#RNA-seq data analysis
#starting from .tab files generated from START aligner
#A.Svendsen/E.Zwart March 2018

library(dplyr)
library(edgeR)
library(biomaRt)
path <- "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis Chapter/Reanalysis/32D_ICD_RNAseq/"
setwd("/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis Chapter/Reanalysis/32D_ICD_RNAseq/")

#Biomart info
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)

#Labels
condition <- factor(rep(c("EV", "ICD", "NFLAG", "CFLAG"), times = 2))
condition <- relevel(condition, ref = "EV")
replicate <- factor(c(rep(1, times = 4), rep(2, times = 4)))

#Importing files
y <- readDGE(list.files("/Users/Art/Drive/PhD/Scripts/Svendsen/2018_Run1_32D_ICD/01_differential_expression/00_files/"),
             path = "/Users/Art/Drive/PhD/Scripts/Svendsen/2018_Run1_32D_ICD/01_differential_expression/00_files",
             columns = c(1,3), #first strand
             group = condition,
             labels = c("EV.1", "ICD.1", "NFLAG.1", "CFLAG.1", "EV.2", "ICD.2", "NFLAG.2", "CFLAG.2"),
             skip = 4) # skip the first 4 lines of each file
#plotMDS(y)
#Filtering
##Low CPM
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes = FALSE]
remove(keep)

#Design
design <- model.matrix(~0 + condition + replicate)
colnames(design) <- c("EV","CFLAG","ICD","NFLAG","Replicate")
design

#Normalization and dispension
y <- calcNormFactors(y,
                     lib.size = y$samples$lib.size,
                     method = "TMM",
                     p = 0.75)

y <- estimateDisp(y, design)
#plotBCV(y)

#General Linear Model
fit <- glmQLFit(y, design)
#lrt <- glmLRT(fit)
#summary(decideTestsDGE(lrt))
EVsICD <- glmLRT(fit, contrast = c(-1, 0, 1, 0, 0))

#Differential Expression
DE <- topTags(EVsICD,
              adjust.method = "BH",
              n = Inf)$table
DE.ICD <- DE[which(DE[,5] <= 0.05),]

#gene annotation
anno <- getBM(
  attributes = c('ensembl_gene_id_version','ensembl_gene_id', 'external_gene_name'),
  filters    = 'ensembl_gene_id',
  values     = substr(rownames(DE.ICD), start = 1, stop = 18),
  mart       = ensembl)

DE.ICD$ensembl_gene_id <- substr(rownames(DE.ICD), start = 1, stop = 18)
merged <- merge(anno, DE.ICD, by = "ensembl_gene_id")
row.names(merged) <- c()
merged <- merged %>%
  dplyr::select(ensembl_gene_id_version, external_gene_name, logFC, logCPM, FDR) %>%
  arrange(FDR)

#write.csv(merged, "01_RNAseq_DEG_20211210_32D_ICD_EV.csv", quote = F, row.names = F)
#write_clip(merged$external_gene_name)


#Volcano plot
overlap_ICD_AS <- c("Neo1", "Sell", "Vldlr", "Evc","Mef2c", "Klrb1b", "Gabra4", "Pclo", "Tbc1d8", "Amotl2")
overlap_ICD_AS <- filter(merged,  external_gene_name %in% overlap_ICD_AS)
  

library(ggplot2)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
ggplot(filter(merged, external_gene_name != "Neo1"), aes(x = logFC, y = -log10(FDR), color = ifelse(logFC > 0, '#3F7F02', "#0F7FFE"))) +
  geom_point(data = filter(DE, logCPM != 8.254274), size = 7,shape = 21, aes(x = logFC, y = -log10(FDR), color = "grey", alpha = 0.15))+
  geom_point(size = 7,shape = 21,  alpha = 0.5)+
  scale_color_identity() +
  geom_vline(xintercept = c(-1,1),  color = "grey") +
  geom_hline(yintercept = 10^-0.05, color = "grey") +
  labs(title = "RNA-seq EV vs. Neo1 ICD", y = "-log10(FDR)",  x = "log2FC")+
  #AS
  ggrepel::geom_text_repel(data = overlap_ICD_AS,
                           aes(x = logFC, y = -log10(FDR), label = external_gene_name, color = "#82027e")) +
  #GOIs
  ggrepel::geom_text_repel(data = filter(merged, logFC >= 0 & -log10(FDR) > 2.5),
                           aes(x = logFC, y = -log10(FDR), label = external_gene_name)) +
  ggrepel::geom_text_repel(data = filter(merged, logFC < -1),
                           aes(x = logFC, y = -log10(FDR), label = external_gene_name)) +
  scale_x_continuous(limits = c(-5,10)) +
  theme_pb() +
  theme(aspect.ratio = 1,
        legend.position = "none")

  