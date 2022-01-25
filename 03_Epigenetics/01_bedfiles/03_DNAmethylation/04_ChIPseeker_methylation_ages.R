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
files <- list(Beerman_gain = "01_Beerman/01_raw_data/03_youngtoold_gain_mm10.bed",
              Beerman_loss = "01_Beerman/01_raw_data/03_youngtoold_loss_mm10.bed",
              Sun_hyper = "03_Sun/01_raw_data/02_Sun_hypermethylated_mm10.bed",
              Sun_hypo = "03_Sun/01_raw_data/02_Sun_hypomethylated_mm10.bed")

peakAnnoList <- lapply(files, annotatePeak, TxDb=edb, tssRegion=c(-5000, 5000), verbose=FALSE)
plotAnnoBar(peakAnnoList)

peakAnnoList <- lapply(files, annotatePeak, TxDb=edb, tssRegion=c(-500, 500), verbose=FALSE)
plotAnnoBar(peakAnnoList)

