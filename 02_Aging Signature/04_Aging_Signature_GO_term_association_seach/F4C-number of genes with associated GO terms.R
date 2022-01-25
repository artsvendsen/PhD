##Alternative way to determine if proposed theories of HSC aging are also related to GO terms found in the AS
## Data was extracted from GO_groups.csv (generated with GO_group jupyter notebook script)
#A.Svendsen Feb 2020
library(tidyverse)
library(viridis)
library(ggrepel)
library(gridExtra)
library(patchwork)
setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")


##Plot of number of GENES with GO terms associated with each search keyword
ref.list <- read.table("00_sourceFiles/01_Aging_List_REAN.txt", header = T, stringsAsFactors = F, sep = "\t")
keywords <- read.csv("04_Aging_Signature_GO_term_association_seach/genes_AS.csv", header = T, stringsAsFactors =T)

df <- data.frame()
df.genes <- c()
for (category in levels(keywords$keyword)){
  genes <- filter(keywords, category == keyword)
  num.genes <- length(filter(keywords, category == keyword)[,1])
  df.genes <-rbind(df.genes, genes)
  df <- rbind(df, cbind.data.frame(category,num.genes))
}

AS <- filter(ref.list, Freq_group > 3)
colnames(AS) <- c("gene", "Freq_group", "Log2FC")
AS$gene <- as.factor(AS$gene)
df3 <- merge(x = df.genes, y = AS, by = "gene")

df$category <- factor(df$category, levels = c("membrane", "inflammation", "cell.cycle", "histone"))

ggplot(df, aes(x = category, y = num.genes, fill = category)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ff7f00","#e41a1c", "#377eb8", "#4daf4a", "#984ea3")) +
  geom_text(aes(label=num.genes)) +
  labs(x = "category", y = "number of genes in AS") +
  theme_pb() +
  theme(aspect.ratio = 1,
        #axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



#Genes and its respective consistency scores

scores <- data.frame(cat = c(), score = c())
for (cat in levels(keywords$keyword)){
  print(cat)
  gene.score <- c()
  cat.score <- c()
  genes.cat <- filter(keywords, keyword == cat)
  for (gene in genes.cat$gene){
    #print(gene)
    gene.score <- filter(ref.list, Gene == gene & Freq_group > 3)[,2]
    cat.score <- c(cat.score, gene.score)
  }
  scores <- rbind.data.frame(scores, 
                             cbind.data.frame(cat = rep(cat, times = length(cat.score)),
                                              score = cat.score
                                              #sum = sum(cat.score),
                                              #norm = sum(cat.score)/length(cat.score)
                             ))
}
scores$cat <- factor(scores$cat, c("membrane", "inflammation", "cell.cycle", "histone"))

p2 <- ggplot(scores, aes(x = cat, y = score)) + 
  geom_count(aes(color = ..n..)) +
  scale_size_area(max_size = 15) +
  scale_color_viridis() +
  scale_y_discrete(limits = 4:12) +
  labs(y = "consistency")+
  theme_pb() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        aspect.ratio = 1)
p2


mydata <- read.csv("04_Aging_Signature_GO_term_association_seach/Frequency of genes in each GO.csv", header = T, stringsAsFactors = F)

histrv<-hist(mydata$num.genes, breaks = 75)
df <- data.frame(breaks = histrv$breaks,
                 counts = c(histrv$counts,0))



ggplot(df, aes(x = breaks, y = counts)) + 
  geom_smooth(se = F) +
  scale_y_log10() +
  labs(y = "GO terms", x = "AS genes") +
  theme_pb() +
  theme(aspect.ratio = 1)
