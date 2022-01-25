#Number of ATAC-seq peaks (broad or narrow from MACS2 p 0.01) overlapping between young and old HSCs
#A.Svendsen Dez 2020
setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/07_Transcriptional_activation/ATAC-seq/")
library(eulerr)

#Broad peaks
fit1 <- euler(c("young" = 48925, 
                "old" = 59528,
                "young&old" = 123569))
plot(fit1,
     quantities = TRUE,
     #fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))

plot(fit1)

#Narrow peaks
fit2 <- euler(c("young" = 12299, 
                "old" = 17159,
                "young&old" = 81485))
plot(fit2,
     quantities = TRUE,
     #fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))

