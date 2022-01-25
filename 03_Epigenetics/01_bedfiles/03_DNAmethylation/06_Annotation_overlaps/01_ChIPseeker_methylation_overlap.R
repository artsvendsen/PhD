# Using Ensembl-based genomic annotations for ChIPseeker
# the package works with UCSC-based annotations, so it's ready for TxDb format,
# one need to generate the complementary data for it to work with Ensembl data

setwd("/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/")

library(GenomicFeatures, verbose = F)
# Load all mouse genomic data from Biomart and generate a TxDb object (it takes a while)
BiomartTxDb <- makeTxDbFromBiomart(biomart= "ensembl", dataset = "mmusculus_gene_ensembl")
#saveRDS(BiomartTxDb, file = "03_Annotation/BiomartTxDb.rds")

edb <- BiomartTxDb
#edb <- readRDS(file = "03_Annotation/BiomartTxDb.rds")

#Change chromossome names from Ensembl to UCSC (data is still intact, it's just for compability)
seqlevelsStyle(edb) <- "UCSC"

library(ChIPseeker, verbose = F)
#load peaks into GRrange format
files <- list(hypermethylation = "05_Studies_overlap/01_hypermethylation.bed",
              hypermethylation_Beerman = "05_Studies_overlap/01_hypermethylation_unique_Beerman.bed",
              hypermethylation_Sun = "05_Studies_overlap/01_hypermethylation_unique_Sun.bed",
              hypomethylation = "05_Studies_overlap/02_hypomethylation.bed",
              hypomethylation_Beerman = "05_Studies_overlap/02_hypomethylation_unique_Beerman.bed",
              hypomethylation_Sun = "05_Studies_overlap/02_hypomethylation_unique_Sun.bed")

peakAnnoList <- lapply(files, annotatePeak, TxDb=edb, tssRegion=c(-5000, 5000), verbose=FALSE)
plotAnnoBar(peakAnnoList)

peakAnnoList <- lapply(files, annotatePeak, TxDb=edb, tssRegion=c(-500, 500), verbose=FALSE)
plotAnnoBar(peakAnnoList)

