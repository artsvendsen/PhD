library(tidyverse)
setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/05_Studies_overlap/zz_closest/")

#Beerman
hyper <- read.table("01_closest_io_hyper.bed", header = F, stringsAsFactors = F, sep = "\t")
hypo <- read.table("02_closest_io_hypo.bed", header = F, stringsAsFactors = F, sep = "\t")

df <- data.frame(sample = c(rep("hypermethylation", times = length(hyper$V10)), rep("hypomethylation", times = length(hypo$V7))),
                 distance = c(hyper$V10, hypo$V7))

ggplot(df, aes(x = distance, fill = sample)) +                       
  geom_histogram(position = "identity", alpha = 0.5, bins = 100) +
  labs(title = "closest distance Beerman vs. Sun") +
  scale_x_log10() +
  theme_bw()
