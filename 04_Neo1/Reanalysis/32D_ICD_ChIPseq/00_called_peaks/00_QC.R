library(ggplot2)
library(dplyr)
library(reshape2)
setwd("/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/00_called_peaks/")


peaks <- data.frame()
files <- list.files("broad/", pattern = "*.broadPeak")
n <- 1
for(file in files){
  print(file)
  print(n)
  temp <- read.table(file = paste0( "broad/", file), header = FALSE, stringsAsFactors = FALSE)
  temp2 <- data.frame(pval = temp$V9,
                      peak = factor(rep("broad", times = length(temp[,1]))),
                      sample = factor(rep(n, times = length(temp[,1]))))
  n <- n + 1
  peaks <- rbind(peaks, temp2)
}

files <- list.files("narrow/", pattern = "*.narrowPeak")
n <- 1
for(file in files){
  print(file)
  print(n)
  temp <- read.table(file = paste0("narrow/", file), header = FALSE, stringsAsFactors = FALSE)
  temp2 <- data.frame(pval = temp$V9,
                      peak = factor(rep("narrow", times = length(temp[,1]))),, 
                      sample = factor(rep(n, times = length(temp[,1]))))
  n <- n + 1
  peaks <- rbind(peaks, temp2)
}
str(peaks)

stats <- data.frame()
for (x in 1:25){
  cut.off <- filter(peaks, pval >= x) %>%
    group_by(peak) %>%
    tally()
  stats <- rbind.data.frame(stats, data.frame(peak = c("broad", "narrow"), 
                                              n =  c(cut.off$n[1], cut.off$n[2]),
                                              perct = c(cut.off$n[1]/470577, cut.off$n[2]/265356),
                                              cutoff = c(x, x))) 
}

stats %>%
  ggplot() +
  geom_vline(xintercept = c(3), color = "grey")+
  geom_line(aes(x = cutoff, y = n/1000000, color = peak)) +
  geom_line(aes(x = cutoff, y = perct, color = peak), linetype = "dashed") +
  scale_x_continuous(limits = c(0,25))+
  theme_bw()


#Write out filtered files
files <- list.files("broad/", pattern = "*.broadPeak")
for (file in files) {
  print(file)
  gapped.df <- read.table(file = paste0("broad/", file), header = FALSE, stringsAsFactors = FALSE)
  df.001 <- filter(gapped.df, V9 >= 4)
  write.table(df.001, paste0("00_filtered_peaks/", gsub("_peaks.broadPeak", "", file), "_pval0.0001.broadPeak"),
  sep = "\t",  row.names = FALSE, col.names = FALSE, quote = FALSE)
}

files <- list.files("narrow/", pattern = "*.narrowPeak")
for (file in files) {
  print(file)
  gapped.df <- read.table(file = paste0("narrow/", file), header = FALSE, stringsAsFactors = FALSE)
  df.001 <- filter(gapped.df, V9 >= 4)
  write.table(df.001, paste0("00_filtered_peaks/", gsub("_peaks.narrowPeak", "", file)[1], "_pval0.0001.narrowPeak"),
              sep = "\t",  row.names = FALSE, col.names = FALSE, quote = FALSE)
}
