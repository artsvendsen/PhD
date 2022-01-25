#!/usr/bash
enhancer_list="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/03_enhancer/amit_enhancer_list_mm10.bed"
ATACfolder="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/"

#Common
bedtools subtract -A -wa -a $enhancer_list -b $ATACfolder/01_ATAC_broad.bed > 01_ATAC_common_enhancers.bed

#Aged speicfic
bedtools subtract -A -wa -a 01_ATAC_common_enhancers.bed -b $ATACfolder/01_ATAC_broad_unique_Old.bed > 02_ATAC_aged_enhancers.bed

#Young speicfic
bedtools subtract -A -wa -a 02_ATAC_aged_enhancers.bed -b $ATACfolder/01_ATAC_broad_unique_Young.bed > 03_ATAC_young_enhancers.bed

wc -l $enhancer_list > number_peaks.txt
wc -l *.bed >> number_peaks.txt