##Extract gene names and directinality from each keyword category
library(tidyverse)
library(officer)
library(flextable)

keywords <- read.csv("04_Aging_Signature_GO_term_association_seach/genes_AS.csv", header = T, stringsAsFactors =T)
direction <- read.csv("00_sourceFiles/05_Aging_Signature_REAN_logFC_piramid.csv", stringsAsFactors = F, header = T)[,1:2]

df <- data.frame()
df.not.found <- c()
for (category in levels(keywords$keyword)){
  temp <- filter(keywords,  keyword == category)
  for (g in temp$gene){
    if( g %in% direction$Genes){
      dir <- filter(direction, Genes == g)[,2]
      df <- rbind(df, cbind.data.frame(g, dir, category))
    }else{
      df.not.found <- c(df.not.found,g)
    }
  }
}

for (g in direction$Genes){
  if (g %in% df$g == FALSE){
    dir <- filter(direction, Genes == g)[,2]
    category <- NA
    df <- rbind(df, cbind.data.frame(g, dir, category))
  }
}



for (cat in levels(df$category)){
  print(cat)
  df.cat <- df %>%
    filter(is.na(category))
  dat <- df.cat
  ft <- flextable(dat)
  ft <- color(ft, color = "#0000FF")
  ft <- color(ft, i = ~ dir > 0, color="#FF0000")
  ft <- void(ft, ~ dir + category)
  ft
  doc <- read_docx()
  doc <- body_add_flextable(doc, value = ft)
  fileout <- tempfile(fileext = ".docx")
  fileout <- paste0("not_found", "_genes.docx") # uncomment to write in your working directory
  print(doc, target = fileout)
}

