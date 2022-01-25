#Differentially Bound Peaks (DBPs) for young and Old ATAC-seq data
#Analysis of broad peaks
#A. Svendsen Dez 2020

library(DiffBind,verbose = F)
library(dplyr)
setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/07_Transcriptional_activation/ATAC-seq/01_diffbind/")

#Broad Peaks
#Load File with metadata information
file <- read.csv("20201219_ATACseq_broad_p0_01.csv", header = T, stringsAsFactors = F)
print(file)

#Load reference data into a dba object
test <- dba(sampleSheet = file)

#differential peak analysis
test <- dba.count(test)
test <- dba.contrast(test, categories = DBA_CONDITION, minMembers = 2) #analyse groups by age (DBA_CONDITION)
test <- dba.analyze(test, method = DBA_DESEQ2, bTagwise = TRUE, bParallel = TRUE)

dba.plotPCA(test, label = DBA_ID)
plot(test)
df <- as.data.frame(dba.report(test))
write.csv(df, paste0("02_ATACseq_broad_p0_01", "_DBP_diffbind.csv"), quote = F, row.names = F)
df.bed <- df %>%
  select(seqnames, start, end) %>%
  arrange(seqnames, start)
write.table(df.bed, paste0("02_ATACseq_broad_p0_01", ".bed"), quote = F, row.names = F, col.names = F, sep = "\t")

#Narrow Peaks
#Load File with metadata information
file <- read.csv("20201219_ATACseq_narrow_p0_01.csv", header = T, stringsAsFactors = F)
print(file)

#Load reference data into a dba object
test <- dba(sampleSheet = file)

#differential peak analysis
test <- dba.count(test)
test <- dba.contrast(test, categories = DBA_CONDITION, minMembers = 2) #analyse groups by age (DBA_CONDITION)
test <- dba.analyze(test, method = DBA_DESEQ2, bTagwise = TRUE, bParallel = TRUE)

dba.plotPCA(test, label = DBA_ID)
plot(test)
df <- as.data.frame(dba.report(test))
write.csv(df, paste0("02_ATACseq_narrow_p0_01", "_DBP_diffbind.csv"), quote = F, row.names = F)
df.bed <- df %>%
  select(seqnames, start, end) %>%
  arrange(seqnames, start)
write.table(df.bed, paste0("02_ATACseq_narrow_p0_01", ".bed"), quote = F, row.names = F, col.names = F, sep = "\t")
