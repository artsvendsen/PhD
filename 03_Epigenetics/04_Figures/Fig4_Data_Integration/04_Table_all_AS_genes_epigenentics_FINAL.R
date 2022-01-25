library(dplyr)
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/04_Figures/Fig4_Data_Integration/")

df <- read.csv("04_Table_all_AS_genes_epigenetics.csv", header = T, stringsAsFactors = F) %>%
  as_tibble() %>%
  mutate(epi.regulation = case_when(grepl(x = accessibility, pattern = "(young|old)") |
                                  grepl(x = H3K4m3, pattern = "(young|old)") |
                                  grepl(x = H3K27m3, pattern = "(young|old)") |
                                  grepl(x = H3K36m3, pattern = "(young|old)") |
                                  grepl(x = DNA, pattern = "(hypomethylated|hypermethylated)") ~ "yes",
                                  TRUE ~ "no")) %>%
  select(Gene,Freq_group, Log2FC, epi.regulation, accessibility, H3K4m3, H3K27m3, H3K36m3, DNA, cluster)

df
write.csv(df, "04_Table_all_AS_genes_epigenetics_FINAL.csv", quote = F, row.names = F)
