#RNA-seq data analysis
#starting from .tab files generated from START aligner
#A.Svendsen/E.Zwart March 2018

library(dplyr)
library(edgeR)
library(biomaRt)
path <- "/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis Chapter/Reanalysis/32D_FL_RNAseq/"
setwd(path)

#Biomart info
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)

#Labels
condition <- factor(c(rep("EV", times = 2), rep("FL", times = 2)))
condition <- relevel(condition, ref = "EV")
replicate <- as.factor(c("1","2","1","2"))

#Importing files
y <- readDGE(list.files("/Users/Art/Drive/PhD/Scripts/Svendsen/201909_32DNeo1FL_RNAseq/02_STAR/"),
             path = "/Users/Art/Drive/PhD/Scripts/Svendsen/201909_32DNeo1FL_RNAseq/02_STAR/",
             columns = c(1,3), #first strand
             group = condition,
             labels = c("EV_1", "EV_3",
                        "FL_1", "FL_2"),
             skip = 4) # skip the first 4 lines of each file
#plotMDS(y)
#Filtering
##Low CPM
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes = FALSE]
remove(keep)

#Design
design <- model.matrix(~0 + condition + replicate)
colnames(design) <- c("EV","FL","Rep2")
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
EVsFL <- glmLRT(fit, contrast = c(-1, 1, 0))

#Differential Expression
DE <- topTags(EVsFL,
              adjust.method = "BH",
              n = Inf)$table
DE.FL <- DE[which(DE[,5] <= 0.05),]

#gene annotation
anno <- getBM(
  attributes = c('ensembl_gene_id_version','ensembl_gene_id', 'external_gene_name'),
  filters    = 'ensembl_gene_id',
  values     = substr(rownames(DE.FL), start = 1, stop = 18),
  mart       = ensembl)

DE.FL$ensembl_gene_id <- substr(rownames(DE.FL), start = 1, stop = 18)
merged <- merge(anno, DE.FL, by = "ensembl_gene_id")
row.names(merged) <- c()
merged <- merged %>%
  dplyr::select(ensembl_gene_id_version, external_gene_name, logFC, logCPM, FDR) %>%
  arrange(FDR)

#write.csv(merged, "01_RNAseq_DEG_20211210_32D_FL_EV.csv", quote = F, row.names = F)
#write_clip(merged$external_gene_name)


#Volcano plot
overlap_FL_AS <- c("Neo1", "Egr1", "Plek", "Btg2","Itgb3", "Mt1")
overlap_FL_AS <- filter(merged,  external_gene_name %in% overlap_FL_AS)


source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")
ggplot(filter(merged, external_gene_name != "Neo1"), aes(x = logFC, y = -log10(FDR), color = ifelse(logFC > 0, '#FC8008', "#0F7FFE"))) +
  geom_point(data = filter(DE, logCPM != 10.296344), size = 7,shape = 21, aes(x = logFC, y = -log10(FDR), color = "grey", alpha = 0.15))+
  geom_point(size = 7,shape = 21)+
  scale_color_identity() +
  geom_vline(xintercept = c(-1,1),  color = "grey") +
  geom_hline(yintercept = 10^-0.05, color = "grey") +
  labs(title = "RNA-seq EV vs. Neo1 ICD", y = "-log10(FDR)",  x = "log2FC")+
  #AS
  ggrepel::geom_text_repel(data = overlap_FL_AS,
                           aes(x = logFC, y = -log10(FDR), label = external_gene_name, color = "#82027e")) +
  #GOIs
  ggrepel::geom_text_repel(data = filter(merged, logFC >= 0 & -log10(FDR) > 3),
                           aes(x = logFC, y = -log10(FDR), label = external_gene_name)) +
  ggrepel::geom_text_repel(data = filter(merged, logFC < -1),
                           aes(x = logFC, y = -log10(FDR), label = external_gene_name)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_x_continuous(limits = c(-5,10)) +
  theme_pb() +
  theme(aspect.ratio = 1,
        legend.position = "none")

#heatmap
#normalized reads
norm_reads <- cpm(y, normalized.lib.sizes = FALSE, log = FALSE)
norm_reads <- as.data.frame(norm_reads)
norm_reads$ensembl_gene_id_version <- substr(rownames(norm_reads), start = 1, stop = 18)
rownames(norm_reads) <- c()
overlap_FL_AS$ensembl_gene_id_version <- substr(overlap_FL_AS$ensembl_gene_id_version, start = 1, stop = 18)
norm_reads_FL <- merge(x = norm_reads, y = overlap_FL_AS, by = "ensembl_gene_id_version") %>%
  dplyr::select(external_gene_name, EV_1, EV_3, FL_1, FL_2) %>%
  {. ->> temp} %>%
  data.matrix() %>%
  .[,-1]
rownames(norm_reads_FL) <- temp$external_gene_name
remove(temp)

#labels
AS <- read.csv("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/04_Aging_Signature_REAN_FC_NA.csv") %>%
  mutate(direction = ifelse(Avg > 0, "up", "do")) %>%
  filter(Genes %in% overlap_FL_AS$external_gene_name) %>%
  dplyr::select(Genes, direction)
AS <- AS[match(rownames(norm_reads_FL), AS$Genes),]
rownames(AS) <- AS$Genes
AS <- dplyr::select(AS, direction) %>%
  mutate(direction = factor(direction, levels = c("up", "do")))

my_sample_col <- data.frame(sample = factor(rep(c("EV", "EV", "FL", "FL"))))
my_sample_col$sample <- factor(my_sample_col$sample, levels = c("FL", "EV"))
row.names(my_sample_col) <- colnames(norm_reads_FL)

my_colour = list(
  AS = c(up = "#0F7FFE", do = "#0F7FFE"),
  my_sample_col = c(EV = "#0F7FFE", FL = "#FC8008"))

#heatmap
pheatmap::pheatmap(norm_reads_FL,
                   annotation_row = AS,
                   annotation_col = my_sample_col,
                   #cluster_rows = F,
                   cluster_cols = F,
                   annotation_colors = my_colour,
                   scale = "row",
                   cellheight=20, cellwidth = 20)

colnames(AS) <- c("Genes", "AS_FC")

FL <- overlap_FL_AS %>%
  dplyr::select(external_gene_name, logFC) %>%
  mutate(direction = ifelse(logFC > 0, "up", "do"))
colnames(FL) <- c("Genes", "FL_FC", "direction")

AS_FL <- inner_join(x = AS, y = FL, by = "Genes") %>%
     {. ->> temp} %>%
      data.matrix() %>%
  .[,-1]
rownames(AS_FL) <- temp$Genes
remove(temp)

pheatmap::pheatmap(AS_FL)
  