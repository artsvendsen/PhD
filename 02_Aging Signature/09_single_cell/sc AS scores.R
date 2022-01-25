## PCA and Aging signature score for single cells from scRNA-seq analysis
##Data extracted from AgeS_RNAseq_SC.py ("generate freqs for pyramid")
#A. Svendsen Sept 2020

library(tidyverse)
library(patchwork)
setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/09_single_cell/")
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

#Import PCA and scores from files
#Grover
df1 <- read.table("01_PCA_Grover_AS_scores.txt",
               stringsAsFactors = F, header = T, sep = "\t") %>%
  mutate(age = factor(.$age, levels = c("yng", "old")))

#merged
grover.merge <- ggplot(df1, aes(x = PC1, y = PC2, color = score)) +
  geom_point(aes(shape = age), size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20)) +
  scale_x_continuous(limits = c(-30,30), breaks = c(-30,-15,0,15,30)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

grover.yng <- ggplot(filter(df1, age == "yng"), aes(x = PC1, y = PC2, color = score)) +
  geom_point(shape = 16, size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20)) +
  scale_x_continuous(limits = c(-30,30), breaks = c(-30,-15,0,15,30)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

grover.old <- ggplot(filter(df1, age == "old"), aes(x = PC1, y = PC2, color = score)) +
  geom_point(shape = 17, size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20)) +
  scale_x_continuous(limits = c(-30,30), breaks = c(-30,-15,0,15,30)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

#Kirshner
df2 <- read.table("01_PCA_Kirshner_AS_scores.txt",
                 stringsAsFactors = F, header = T, sep = "\t") %>%
  mutate(age = factor(.$age, levels = c("yng", "old")))

#merged
kir.merge <- ggplot(df2, aes(x = PC1, y = PC2, color = score)) +
  geom_point(aes(shape = age), size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20)) +
  scale_x_continuous(limits = c(-25,25)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

kir.yng <- ggplot(filter(df2, age == "yng"), aes(x = PC1, y = PC2, color = score)) +
  geom_point(shape = 16, size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20)) +
  scale_x_continuous(limits = c(-25,25)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

kir.old <- ggplot(filter(df2, age == "old"), aes(x = PC1, y = PC2, color = score)) +
  geom_point(shape = 17, size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20)) +
  scale_x_continuous(limits = c(-25,25)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

#kir.yng | kir.old | kir.merge

#Kowalczyk
df3 <- read.table("01_PCA_Kowalczyk_AS_scores.txt",
                 stringsAsFactors = F, header = T, sep = "\t") %>%
  mutate(age = factor(.$age, levels = c("yng", "old")))

#merged
kow.merge <- ggplot(df3, aes(x = PC1, y = PC2, color = score)) +
  geom_point(aes(shape = age), size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20), labels = c(-15,-7.5,0,7.5,15)) +
  scale_x_continuous(limits = c(-20,20), labels = c(-15,-7.5,0,7.5,15)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

kow.yng <- ggplot(filter(df3, age == "yng"), aes(x = PC1, y = PC2, color = score)) +
  geom_point(shape = 16, size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20), labels = c(-15,-7.5,0,7.5,15)) +
  scale_x_continuous(limits = c(-20,20), labels = c(-15,-7.5,0,7.5,15)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

kow.old <- ggplot(filter(df3, age == "old"), aes(x = PC1, y = PC2, color = score)) +
  geom_point(shape = 17, size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20), labels = c(-15,-7.5,0,7.5,15)) +
  scale_x_continuous(limits = c(-20,20), labels = c(-15,-7.5,0,7.5,15)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

#kow.yng | kow.old | kow.merge

#Mann
df4 <- read.table("01_PCA_Mann_AS_scores.txt",
                 stringsAsFactors = F, header = T, sep = "\t") %>%
  mutate(age = factor(.$age, levels = c("yng", "old")))

#merged
mann.merge <- ggplot(df4, aes(x = PC1, y = PC2, color = score)) +
  geom_point(aes(shape = age), size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20), labels = c(-20,-10,0,10,20)) +
  scale_x_continuous(limits = c(-20,20), labels = c(-20,-10,0,10,20)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

mann.yng <- ggplot(filter(df4, age == "yng"), aes(x = PC1, y = PC2, color = score)) +
  geom_point(shape = 16, size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20), labels = c(-20,-10,0,10,20)) +
  scale_x_continuous(limits = c(-20,20), labels = c(-20,-10,0,10,20)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

mann.old <- ggplot(filter(df4, age == "old"), aes(x = PC1, y = PC2, color = score)) +
  geom_point(shape = 17, size = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(-20,20), labels = c(-20,-10,0,10,20)) +
  scale_x_continuous(limits = c(-20,20), labels = c(-20,-10,0,10,20)) +
  scale_color_viridis_c(limits = c(0,1)) +
  theme_pb() +
  theme(aspect.ratio = 1)

(grover.yng | grover.old | grover.merge) /
(kir.yng | kir.old | kir.merge) /
(kow.yng | kow.old | kow.merge) /
(mann.yng | mann.old | mann.merge)

# Violin plots of sc AS scores
df1$study <- "Grover"
df2$study <- "Kirshner"
df3$study <- "Kowalczyk"
df4$study <- "Mann"
df <- rbind.data.frame(df1,df2,df3,df4)


ggplot(df, aes(x = age, y = score, fill = age)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NULL, fill = "white") +
  scale_fill_manual(values = c('grey', 'red')) +
  facet_grid(~ study) +
  theme_pb() +
  theme(aspect.ratio = 1)


#Significance per gene/cell in different studies
studies <- c("Grover", "Kirshner", "Kowalczyk", "Mann")
mydata <- data.frame()
for (study in studies){
  temp <- c()
  temp <- read.table(paste0("02_significance_",study,".txt"), header = T, stringsAsFactors = F)
  temp$study <- as.factor(study)
  mydata <- rbind.data.frame(mydata, temp)
}
mydata <- mydata %>% 
  mutate(sig = ifelse(pval <0.01 , .$pval, NA))

top <- mydata %>% 
  group_by(study) %>% 
  filter(!is.na(sig)) %>%
  dplyr::top_n(10, wt = -pval)

ggplot(mapping = aes(x,y, label = gene)) +
  geom_point(data = mydata, aes(color = ifelse(pval >= 0.01, 'grey', 'red')), size = 3)+
  scale_color_identity() +
  ggrepel::geom_text_repel(data = top,
                           show.legend = FALSE) +
  facet_grid(~ study, scales = "free") +
  theme_pb() +
  theme(aspect.ratio = 1)

