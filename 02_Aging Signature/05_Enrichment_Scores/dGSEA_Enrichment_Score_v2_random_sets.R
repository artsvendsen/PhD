x <- Sys.time()
##Calculation of discrete GSEA(dGSEA)
#A.Svendsen/E.Zwart Aug 2019
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(patchwork)
setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R") 

###!WARNING!###
#CHANGE DIRECTORY FILES IN LINES 8, 16,65,77

#Load reference source file
ref.list <- read.csv("00_sourceFiles/01_Aging_List_REAN.txt", header = T, stringsAsFactors = F, sep = "\t")
#Sort ref.list by Freq_group (high to low)
ref.list <- ref.list[order(-ref.list$Freq_group),]
colnames(ref.list) <- c("gene", "consistency", "log2FC")
#Add indexed as a seperate collumn
ref.list$index <- 1:length(ref.list[,1])

#Make control points for all genes in ref.list (long list)
control.points <- data.frame(index = c(),
                             consistency = c())
for (i in ref.list$index){
  if (i == 1){
    next #skipping this index 1 because otherwise there's a index 0 2 lines down.
  }else {
    if (ref.list[i, 2] != ref.list[i - 1, 2]){ #get index of the last gene in each consistency category
      control.points <- rbind(control.points, c(i-1, ref.list[i-1, 2])) #append index and consistency value to a new df (control.points)
    }
  }
}

#Select aging signature genes
aging.sig <- filter(ref.list, consistency >= 4)


#Add the index of the last gene in longlist
control.points <- rbind(control.points, c(nrow(ref.list), min(ref.list$consistency)))

colnames(control.points) <- c("index", "consistency")
control.points


##Calculations of Rank score; total Enrichment score; cumulative Enrichment score and control points present in each study of signature
rean_files <- list.files("00_sourceFiles/rean/", pattern = "*.csv")
studies <- c()
for (i in rean_files){
  temp <- unlist(str_split(i, fixed("_")))[1]
  studies <- c(studies, temp)
}
studies 
studies.list <- list()
short.list <- c()
temp2 <- c()
Scores <- c()
control.points.plot <- list()
floating.EnrichementScore <- data.frame()
#df <- c()
#debug <- c()
for (study in rean_files){
  #Read csv file of individual studyes and store in a list
  studies.list[[study]] <- read.csv(paste0("00_sourceFiles/rean/",
                                      study), header = T, stringsAsFactors = F, skip = 1, sep = "\t")[,2]
}
names(studies.list) <- studies

studies.list[["Aging_signature"]] <- filter(ref.list, consistency > 1)$gene
studies <- c(studies, "Aging_signature")
##
#Load study that will be compared to the studies in the aging signature// it will be appended into the list
##
studies <- c()
test.studies <- list.files("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/05_Enrichment_Scores/Random_sets/", pattern = "_genes.csv")
# test.studies <- c("WONG_EMBRYONIC_STEM_CELL_CORE", "Challen", "MARYANOVICH_DO", "Yamashita_UP_3h_TNFa", "von_Eyss", "von_Eyss-Top250")
for(set in test.studies){
  studies.list[[set]] <- read.csv(paste0("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/05_Enrichment_Scores/Random_sets/",
                                         set), header = T, stringsAsFactors = F)[,1]
  studies <- c(studies, set)  
}

##

# Perform calculation scores for source lists and test list
for (study in studies){
  print(study)
  #collect data from ref.list for individual genes in a study at a time and store in another list
  short.list <- c()
  for (gene.study in studies.list[[study]]) {
    formated.gene <- gsub(" ", "", unlist(strsplit(gene.study, "///"))[1]) #Remove potential whitespaces and "///"
    short.list <- rbind(short.list, filter(ref.list, gene == formated.gene)[,c(1,2,4)])
    
  }
  #debug <- rbind.data.frame(debug, cbind.data.frame(study, length(studies.list[[study]]), length(short.list[,1])))
  
  #Get all genes which are not in the ref.list (will be used for penalty)
  penalty.list <- anti_join(data.frame(gene = studies.list[[study]]), short.list, by = "gene")
  
  #Rank Score
  if (nrow(short.list) == 0) {
    RankScore <- 0
  } else {
  RankScore <- short.list %>%
    filter(consistency > 1) %>%
    select(consistency) %>%
    sum()
  }
  if (min(short.list$consistency) == 1) {
    penalty.1 <- short.list %>%
      filter(consistency == 1) %>%
      select(consistency) %>%
      sum() * -1
  } else {
    penalty.1 <- 0
  }
  
  penalty.NULL <- -1 * length(penalty.list$gene)
  
  RankScore.corrected <- RankScore + penalty.1 + penalty.NULL
  
  #Total Enrichment Score
  if (study == "Aging_signature"){
    short.list <- short.list[1:length(filter(ref.list, consistency >= 4)[,1]),]
  } else {
    short.list <- short.list
  }
  obeserved <- length(base::intersect(aging.sig$gene, short.list$gene))
  expected <- nrow(short.list) * (nrow(aging.sig) / nrow(ref.list))
  EnrichementScore <- (obeserved - expected) /  nrow(short.list) * 100
  temp2 <- cbind.data.frame(study, RankScore.corrected, EnrichementScore)
  Scores <- rbind(Scores, temp2)
  #df <- rbind.data.frame(df,cbind.data.frame(study, obeserved, expected, EnrichementScore))

  #Enrichment score for dGSEA plot
  for (i in control.points$index){
    found <- filter(short.list, index <= i)
    floating.obeserved <- nrow(found)
    floating.expected <- i * (nrow(short.list) / nrow(ref.list))
    floating.EnrichementScore <- rbind.data.frame(floating.EnrichementScore, cbind.data.frame(study, i, c(floating.obeserved - floating.expected)))
  }

  #Collect control points for individual study
  #control.points.plot <- rbind.data.frame(control.points.plot, cbind.data.frame(study, short.list))

  short.list <- c()
  temp2 <- c()
  RankScore <- c()
  penalty <- c()
}
#Formatting col names  
colnames(floating.EnrichementScore) <- c("study", "index", "float.Ench.Score")
colnames(Scores) <- c("Study", "Rank.Score", "Enrich.Score")
Scores
