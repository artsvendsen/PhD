library(tidyverse)
library(pheatmap)
library(viridis)
studies <- c("Grover", "Kirshner", "Kowalczyk", "Mann")
AS <- read.csv("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/04_Aging_Signature_REAN_FC_NA.csv") %>%
  arrange(Order) %>%
  select(Genes)

for (study in studies){
  #Import sc data
  df <- read.csv(paste0("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/09_single_cell/normalized_counts_",
                        study,
                        ".csv"), header = T, stringsAsFactors = F)
  df_merge <- left_join(x = AS, y = df, by = "Genes")
  df_num <- as.matrix(df_merge[,-1])
  
  #Young/Old legend
  cat_df <- data.frame("age" = gsub("\\..*","",colnames(df[,-1]))) %>%
    mutate(age = factor(age, levels = c("young", "old")))
  rownames(cat_df) = colnames(df_num)
  my_color <- list(age = c(young = 'grey', old = 'red'))
  
  pheatmap(t(df_num), 
           main = paste0(study),
           annotation_row = cat_df,
           cellwidth = 2, cellheight = 2,
           show_rownames = F, show_colnames = F,
           color = viridis(10),
           annotation_colors = my_color,
           cluster_cols = F, 
           na_col = 'grey')
}

