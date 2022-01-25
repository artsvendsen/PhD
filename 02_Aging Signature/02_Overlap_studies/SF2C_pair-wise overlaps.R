#Find the nubmer of overlapping genes in paired-wide fashion across reanl and reanalyzed datasets
#A. Svendsen
#Jul 2019 (reviewed Jan 2020)

setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/")
library(stringr)
library(tidyverse)
library(reshape2)
library(patchwork)
library(viridis)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

#Meta-analysis
meta_files <- list.files("meta/")
meta.studies <- c()
raw.meta.list <- list()
clean.meta.list <- list()
for (study in meta_files){
  study.name <- unlist(str_split(study, fixed("_")))[1]
  meta.studies <- rbind(meta.studies, study.name)
  raw.meta.list[[study.name]] <- read.table(paste0("meta/",study), header = F, stringsAsFactors = F, sep = "\t", skip = 2)[,2]
  
  for (gene.study in raw.meta.list[[study.name]]){
    temp <- gsub(" ", "", unlist(strsplit(gene.study, "///"))[1]) #Remove potential whitespaces and "///"
    clean.meta.list[[study.name]] <- rbind(clean.meta.list[[study.name]], temp)
  }
}
remove(raw.meta.list, temp)

#Meta-analysis number of paired comparisons
meta.num.overlaps <- data.frame()
meta.genes.overlaps <- data.frame()
for (a in meta.studies){
  for (b in meta.studies){
    temp <- length(Reduce(intersect,list(clean.meta.list[[a]], clean.meta.list[[b]])))
    temp2 <- Reduce(intersect,list(clean.meta.list[[a]], clean.meta.list[[b]]))
    meta.num.overlaps <- rbind(meta.num.overlaps, cbind(a,b,temp))
  }
}
colnames(meta.num.overlaps) <- c("intersect.a", "intersect.b", "num.overlap")

#Input Matrixes for heatmaps (number of overlaps and pearson correlation)
#number of overlaps
meta.num.overlaps.wide <- meta.num.overlaps %>%
  spread(key = intersect.a, value = num.overlap, convert = TRUE) %>%
  tibble::column_to_rownames(var = "intersect.b") %>%
  as.matrix()

#Correlation (Pearson)
meta.num.overlaps.cor <- round(cor(meta.num.overlaps.wide, method = "pearson"),2)

# Heatmap nubmer of overlapping genes
#http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

#Remove redundant boxes
cormat <- reorder_cormat(meta.num.overlaps.wide)
lower_tri <- get_upper_tri(cormat)
matrix.melt <- melt(lower_tri, na.rm = TRUE)

#Import fraction of intersects of Meta and ReAnl sets
intersect.fraction.meta <- filter(read.csv("07_fraction_intersected.csv", header = T, stringsAsFactors = F), analysis == "meta") %>%
  mutate(study = factor(study, levels = levels(matrix.melt$Var1)))

#Heatmap with number of overlapping studies + fraction of intersects
p1 <- ggplot(data = matrix.melt, aes(x = Var2, y = Var1, fill = log(value)))+
  scale_fill_viridis(alpha = 0.7) +
  geom_tile(color = "white")+
  geom_text(aes(x = Var2, y = Var1, label = value), color = "white", size = 3) +
  labs(y = " ", x = " ") +
  theme_pb()+
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = 'black'),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'),
        aspect.ratio = 1)

p1
p2 <- ggplot(data = intersect.fraction.meta, aes(x = study, y = study, fill = fraction_intersected)) +
  scale_fill_viridis(begin = 0, end = 1, alpha = 0.7) +
  geom_tile(color = "white")+
  geom_text(aes(x = study, y = study, label = round(fraction_intersected, 2)), color = "white", size = 3) +
  labs(y = " ", x = " ") +
  theme_pb()+
  theme(
    #legend.position = 'none',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = 'black'),
    legend.text = element_text(size = 15, color = 'black'),
    legend.title = element_text(size = 15, color = 'black'),
    aspect.ratio = 1)

p2
#p1 + p2



#Re-analysis
rean_files <- list.files("rean/", pattern = "*.csv")
rean.studies <- c()
raw.rean.list <- list()
clean.rean.list <- list()
for (study in rean_files){
  study.name <- unlist(str_split(study, fixed("_")))[1]
  rean.studies <- rbind(rean.studies, study.name)
  raw.rean.list[[study.name]] <- read.csv(paste0("rean/",study), header = F, stringsAsFactors = F, sep = "\t", skip = 2)[,2]
  
  for (gene.study in raw.rean.list[[study.name]]){
    temp <- gsub(" ", "", unlist(strsplit(gene.study, "///"))[1]) #Remove potential whitespaces and "///"
    clean.rean.list[[study.name]] <- rbind(clean.rean.list[[study.name]], temp)
  }
}
remove(raw.rean.list, temp)

#rean-analysis number of paired comparisons
rean.num.overlaps <- data.frame()
rean.genes.overlaps <- data.frame()
for (a in rean.studies){
  for (b in rean.studies){
    temp <- length(Reduce(intersect,list(clean.rean.list[[a]], clean.rean.list[[b]])))
    temp2 <- Reduce(intersect,list(clean.rean.list[[a]], clean.rean.list[[b]]))
    rean.num.overlaps <- rbind(rean.num.overlaps, cbind(a,b,temp))
  }
  #rean.genes.overlaps <- rbind(b,rean.genes.overlaps, temp2)
}
colnames(rean.num.overlaps) <- c("intersect.a", "intersect.b", "num.overlap")

#Input Matrixes for heatmaps (number of overlaps and pearson correlation)
#number of overlaps
rean.num.overlaps.wide <- rean.num.overlaps %>%
  spread(key = intersect.a, value = num.overlap, convert = TRUE) %>%
  tibble::column_to_rownames(var = "intersect.b") %>%
  as.matrix()

#Correlation (Pearson)
rean.num.overlaps.cor <- round(cor(rean.num.overlaps.wide, method = "pearson"),2)

# Heatmap nubmer of overlapping genes
#http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


#Remove redundant boxes
upper_tri <- get_upper_tri(cormat)
matrix.melt <- melt(upper_tri, na.rm = TRUE)

# lower_tri <- get_lower_tri(rean.num.overlaps.wide)
# matrix.melt <- melt(lower_tri, na.rm = TRUE)
#matrix.melt$Var1 <- factor(matrix.melt$Var1, levels = unique(matrix.melt$Var1[order(desc(matrix.melt$Var1))]))
#matrix.melt$Var2 <- factor(matrix.melt$Var2, levels = unique(matrix.melt$Var2[order(desc(matrix.melt$Var2))]))

#Import fraction of intersects of rean and ReAnl sets
intersect.fraction.rean <- filter(read.csv("07_fraction_intersected.csv", header = T, stringsAsFactors = F), analysis == "reanl") %>%
  mutate(study = factor(study, levels = levels(matrix.melt$Var1)))

#Heatmap with number of overlapping studies + fraction of intersects
p3 <- ggplot(data = matrix.melt, aes(x = Var2, y = Var1, fill = log(value)))+
  scale_fill_viridis(alpha = 0.7) +
  geom_tile(color = "white")+
  geom_text(aes(x = Var2, y = Var1, label = value), color = "white", size = 3) +
  #scale_fill_gradient2(low = "blue", high = "red", mid = "grey", midpoint = 0.5, space = "Lab") +
  labs(y = " ",
       x = " ") +
  theme_pb()+
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = 'black'),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'),
        aspect.ratio = 1)

p3
p4 <- ggplot(data = intersect.fraction.rean, aes(x = study, y = study, fill = fraction_intersected)) +
  scale_fill_viridis(begin = 0, end = 1, alpha = 0.7) +
  geom_tile(color = "white")+
  geom_text(aes(x = study, y = study, label = round(fraction_intersected, 2)), color = "white", size = 3) +
  labs(y = " ", x = " ") +
  theme_pb()+
  theme(
    #legend.position = 'none',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = 'black'),
    legend.text = element_text(size = 15, color = 'black'),
    legend.title = element_text(size = 15, color = 'black'),
    aspect.ratio = 1)

p4
