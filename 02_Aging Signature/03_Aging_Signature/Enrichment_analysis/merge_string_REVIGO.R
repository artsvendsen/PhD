library(tidyverse)
setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
categories <- c("Process", "Component", "Function")
directions <- c("UP", "DO")


for (cat in categories){
  print(cat)
  for (dir in directions){
    print(dir)
    UP <- read.csv(paste0("03_Aging_Signature/Enrichment_analysis/enrichment.",cat,"_",dir,".tsv"), stringsAsFactors = F, header = T, sep = "\t")
    DO <- read.csv(paste0("03_Aging_Signature/Enrichment_analysis/enrichment.",cat,"_",dir,".tsv"), stringsAsFactors = F, header = T, sep = "\t")
    filtered <- read.csv(paste0("03_Aging_Signature/Enrichment_analysis/unique_",cat,"_",dir,"_filtered.csv"), stringsAsFactors = F, header = T)
    df <-c()
    for (GO in filtered[,1]){
      if (GO != ""){
        print(GO)
        df <- rbind(df, filter(original, X.term.ID == GO))
        write.csv(df, paste0("03_Aging_Signature/Enrichment_analysis/zz_unique_",cat,"_",dir,"_filtered_original_merged.csv")) 
      }
    }
  }
}

# UP <- read.csv("03_Aging_Signature/Enrichment_analysis/enrichment.Process_UP.tsv", stringsAsFactors = F, header = T, sep = "\t")
# DO <- read.csv("03_Aging_Signature/Enrichment_analysis/enrichment.Process_DO.tsv", stringsAsFactors = F, header = T, sep = "\t")
# ref <- read.csv("03_Aging_Signature/common_GOs_process.csv", stringsAsFactors = F)
# 
# df <- c()
# for (GO in ref$erm.ID){
#   temp.UP <- filter(UP, X.term.ID == GO)[,c(2,3,4,5)]
#   temp.DO <- filter(DO, X.term.ID == GO)[,c(3,4,5)]
#   df <- rbind(df, cbind(temp.UP, temp.DO))
# }
# write.csv(df, "merged_UP_DO_filtered.csv", quote = F, row.names = F)
