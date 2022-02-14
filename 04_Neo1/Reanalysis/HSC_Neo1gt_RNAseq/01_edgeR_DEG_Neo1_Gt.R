#RNA-seq data analysis
#starting from .tab files generated from START aligner
#A.Svendsen/E.Zwart March 2018

library(dplyr)
library(edgeR)
library(biomaRt)
path <- "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/HSC_Neo1gt_RNAseq/"
setwd(path)

#Biomart info
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)

# Young samples
#Labels
condition <- factor(c(rep("WT", times = 2), rep("KO", times = 2)))
replicate <- as.factor(c("1","2","1","2"))

#Importing files
y <- readDGE(list.files("/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/HSC_Neo1gt_RNAseq/01_readcounts/")[c(3,4,8,9)],
             path = "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/HSC_Neo1gt_RNAseq/01_readcounts/",
             columns = c(1,4),
             group = condition,
             labels = c("Young_WT_1", "Young_WT_2","Young_KO_1","Young_KO_2"),
             skip = 4) # skip the first 4 lines of each file
y$samples
plotMDS(y)
#Filtering
##Low CPM
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes = FALSE]
remove(keep)

#Design
design <- model.matrix(~0 + condition + replicate)
colnames(design) <- c("KO","WT", "rep2")
design

#Normalization and dispension
y <- calcNormFactors(y,
                     lib.size = y$samples$lib.size,
                     method = "TMM",
                     p = 0.75)

y <- estimateDisp(y, design)
plotBCV(y)

#General Linear Model
fit <- glmQLFit(y, design)
#lrt <- glmLRT(fit)
#summary(decideTestsDGE(lrt))
WTsGT <- glmLRT(fit, contrast = c(-1, 1, 0))

#Differential Expression
DE <- topTags(WTsGT,
              adjust.method = "BH",
              n = Inf)$table
DE.GT <- DE[which(DE[,5] <= 0.05),]

#gene annotation
anno <- getBM(
  attributes = c('ensembl_gene_id_version','ensembl_gene_id', 'external_gene_name'),
  filters    = 'ensembl_gene_id',
  values     = substr(rownames(DE.GT), start = 1, stop = 18),
  mart       = ensembl)

DE.GT$ensembl_gene_id <- substr(rownames(DE.GT), start = 1, stop = 18)
merged <- merge(anno, DE.GT, by = "ensembl_gene_id")
row.names(merged) <- c()
merged <- merged %>%
  dplyr::select(ensembl_gene_id_version, external_gene_name, logFC, logCPM, FDR) %>%
  arrange(FDR)

#write.csv(merged, "01_RNAseq_DEG_20211210_HSC_GT.csv", quote = F, row.names = F)
#write_clip(merged$external_gene_name)

#normalized reads
norm_reads_y <- cpm(y, normalized.lib.sizes = FALSE, log = FALSE)
norm_reads_y <- as.data.frame(norm_reads_y)
norm_reads_y$ensembl_gene_id <- substr(rownames(norm_reads_y), start = 1, stop = 18)
rownames(norm_reads_y) <- c()
merged$ensembl_gene_id <- substr(merged$ensembl_gene_id_version, start = 1, stop = 18)
norm_reads_Young <- merge(x = norm_reads_y, y = merged, by = "ensembl_gene_id") %>%
  dplyr::select(external_gene_name, Young_WT_1, Young_WT_2, Young_KO_1, Young_KO_2) %>%
  {. ->> temp} %>%
  data.matrix() %>%
  .[,-1]
rownames(norm_reads_Young) <- temp$external_gene_name
remove(temp)

#Volcano plot
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
ggplot(merged, aes(x = logFC, y = -log10(FDR), color = ifelse(logFC > 0, '#FC8008', "#0F7FFE"))) +
  geom_point(data = DE, size = 7,shape = 21, aes(x = logFC, y = -log10(FDR), color = "grey", alpha = 0.15))+
  geom_point(size = 7,shape = 21)+
  scale_color_identity() +
  geom_vline(xintercept = c(-1,1),  color = "grey") +
  geom_hline(yintercept = 10^-0.05, color = "grey") +
  labs(title = "RNA-seq EV vs. Neo1 ICD", y = "-log10(FDR)",  x = "log2FC")+
  ggrepel::geom_text_repel(data = filter(merged, logFC >= 2 & -log10(FDR) > 2.5),
                           aes(x = logFC, y = -log10(FDR), label = external_gene_name)) +
  ggrepel::geom_text_repel(data = filter(merged, logFC < -1),
                           aes(x = logFC, y = -log10(FDR), label = external_gene_name)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_x_continuous(limits = c(-5,10)) +
  theme_pb() +
  theme(aspect.ratio = 1,
        legend.position = "none")

# Young samples
#Labels
condition <- factor(c(rep("WT", times = 2), rep("KO", times = 3)))
replicate <- as.factor(c("1","2","1","2","3"))

#Importing files
y <- readDGE(list.files("/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/HSC_Neo1gt_RNAseq/01_readcounts/")[c(1,2,5,6,7)],
             path = "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/HSC_Neo1gt_RNAseq/01_readcounts/",
             columns = c(1,4),
             group = condition,
             labels = c("Old_WT_1", "Old_WT_2","Old_KO_1", "Old_KO_2", "Old_KO_3"),
             skip = 4) # skip the first 4 lines of each file
y$samples
plotMDS(y)
#Filtering
##Low CPM
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes = FALSE]
remove(keep)

#Design
design <- model.matrix(~0 + condition + replicate)
colnames(design) <- c("KO","WT", "rep2", "rep3")
design

#Normalization and dispension
y <- calcNormFactors(y,
                     lib.size = y$samples$lib.size,
                     method = "TMM",
                     p = 0.75)

y <- estimateDisp(y, design)
plotBCV(y)

#General Linear Model
fit <- glmQLFit(y, design)
#lrt <- glmLRT(fit)
#summary(decideTestsDGE(lrt))
WTsGT <- glmLRT(fit, contrast = c(-1, 1, 0, 0))

#Differential Expression
DE <- topTags(WTsGT,
              adjust.method = "BH",
              n = Inf)$table
DE.GT <- DE[which(DE[,5] <= 0.05),]

#gene annotation
anno <- getBM(
  attributes = c('ensembl_gene_id_version','ensembl_gene_id', 'external_gene_name'),
  filters    = 'ensembl_gene_id',
  values     = substr(rownames(DE.GT), start = 1, stop = 18),
  mart       = ensembl)

DE.GT$ensembl_gene_id <- substr(rownames(DE.GT), start = 1, stop = 18)
merged <- merge(anno, DE.GT, by = "ensembl_gene_id")
row.names(merged) <- c()
merged <- merged %>%
  dplyr::select(ensembl_gene_id_version, external_gene_name, logFC, logCPM, FDR) %>%
  arrange(FDR)

#write.csv(merged, "01_RNAseq_DEG_20211210_HSC_GT_old.csv", quote = F, row.names = F)
#write_clip(merged$external_gene_name)

#normalized reads
norm_reads_o <- cpm(y, normalized.lib.sizes = FALSE, log = FALSE)
norm_reads_o <- as.data.frame(norm_reads_o)
norm_reads_o$ensembl_gene_id <- substr(rownames(norm_reads_o), start = 1, stop = 18)
rownames(norm_reads_o) <- c()
merged$ensembl_gene_id <- substr(merged$ensembl_gene_id_version, start = 1, stop = 18)
norm_reads_Old <- merge(x = norm_reads_o, y = merged, by = "ensembl_gene_id") %>%
  dplyr::select(external_gene_name, Old_WT_1, Old_WT_2, Old_KO_1, Old_KO_2) %>%
  {. ->> temp} %>%
  data.matrix() %>%
  .[,-1]
rownames(norm_reads_Old) <- temp$external_gene_name
remove(temp)

#Volcano plot
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
ggplot(merged, aes(x = logFC, y = -log10(FDR), color = ifelse(logFC > 0, '#FC8008', "#0F7FFE"))) +
  geom_point(data = DE, size = 7,shape = 21, aes(x = logFC, y = -log10(FDR), color = "grey", alpha = 0.15))+
  geom_point(size = 7,shape = 21)+
  scale_color_identity() +
  geom_vline(xintercept = c(-1,1),  color = "grey") +
  geom_hline(yintercept = 10^-0.05, color = "grey") +
  labs(title = "RNA-seq EV vs. Neo1 ICD", y = "-log10(FDR)",  x = "log2FC")+
  ggrepel::geom_text_repel(data = filter(merged, logFC >= 2 & -log10(FDR) > 2.5),
                           aes(x = logFC, y = -log10(FDR), label = external_gene_name)) +
  ggrepel::geom_text_repel(data = filter(merged, logFC < -1),
                           aes(x = logFC, y = -log10(FDR), label = external_gene_name)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_x_continuous(limits = c(-5,10)) +
  theme_pb() +
  theme(aspect.ratio = 1,
        legend.position = "none")



# Heatmap of overlapping genes of young and aged
#merge young and aged
all_norm_reads <- merge(x = norm_reads_Young, y = norm_reads_Old, by = "row.names")
rownames(all_norm_reads) <- all_norm_reads$Row.names
all_norm_reads <- all_norm_reads %>%
  dplyr::select(-Row.names)

#labels
my_sample_col <- data.frame(age = c(rep("Young", times = 4), rep("Aged", times = 4)),
                            group = rep(c("WT", "WT", "Gt", "Gt"), times = 2))
my_sample_col$age <- factor(my_sample_col$age, levels = c("Young", "Aged"))
my_sample_col$group <- factor(my_sample_col$group, levels = c("WT", "Gt"))
row.names(my_sample_col) <- colnames(all_norm_reads)

#heatmap
pheatmap::pheatmap(all_norm_reads,
                   #annotation_row = AS,
                   annotation_col = my_sample_col,
                   cluster_rows = T,
                   cluster_cols = T,
                   # annotation_colors = my_colour,
                   scale = "row",
                   cutree_rows = 2,
                   cellheight=20, cellwidth = 20)
