#!/usr/bash
# Overlap of ATAC-seq peaks with promoter list for all hematopoeitic populations from Lara-Astiaso pub.
# List was lifted over from mm9 to mm10 by Eric Martin (UCSS) and emailed.
# A. Svendsen Aug 2021

env="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/04_Promoters/"
cd $env
promoter_list="EPC_promoter_list_sorted.bed"

#Broad
#All young peaks
bedtools intersect -wb -b $promoter_list -a "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/00_meged_replicates/01_ATAC_Young_broad.bed" > 01_promoter_overlap_merged_Young_broad.bed
#All old peaks
bedtools intersect -wb -b $promoter_list -a "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/00_meged_replicates/02_ATAC_Old_broad.bed" > 01_promoter_overlap_merged_Old_broad.bed
#Shared peaks between young and old
bedtools intersect -wb -b $promoter_list -a "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad.bed" > 02_promoter_overlap_common_broad.bed
#Young only
bedtools intersect -wb -b $promoter_list -a "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad_unique_Young.bed" > 02_promoter_overlap_unique_Young_broad.bed
#Old only
bedtools intersect -wb -b $promoter_list -a "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/01_ATAC_broad_unique_Old.bed" > 02_promoter_overlap_unique_Old_broad.bed


#Narrow
#All young peaks
bedtools intersect -wb -b $promoter_list -a "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/00_meged_replicates/03_ATAC_Young_narrow.bed" > 03_promoter_overlap_merged_Young_narrow.bed
#All old peaks
bedtools intersect -wb -b $promoter_list -a "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/00_meged_replicates/04_ATAC_Old_narrow.bed" > 03_promoter_overlap_merged_Old_narrow.bed
#Shared peaks between young and old
bedtools intersect -wb -b $promoter_list -a "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/02_ATAC_narrow.bed" > 04_promoter_overlap_common_narrow.bed
#Young only
bedtools intersect -wb -b $promoter_list -a "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/02_ATAC_narrow_unique_Young.bed" > 04_promoter_overlap_unique_Young_narrow.bed
#Old only
bedtools intersect -wb -b $promoter_list -a "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/02_ATAC_narrow_unique_Old.bed" > 04_promoter_overlap_unique_Old_narrow.bed

wc -l *.bed > num_peaks.txt