#!/usr/bash
# Overlap of Beerman and Sun merge regions with promoter list for all hematopoeitic populations from Lara-Astiaso pub.
# List was lifted over from mm9 to mm10 by Eric Martin (UCSS) and emailed.
# A. Svendsen Aug 2021

env="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/05_Studies_overlap"
cd $env
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/08_promoter"


bedtools intersect -wb -b "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/04_promoter/EPC_promoter_list_sorted.bed" -a 01_hypermethylation.bed > $out/01_promoter_hypermethylation_common.bed
bedtools intersect -wb -b "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/04_promoter/EPC_promoter_list_sorted.bed" -a 01_hypermethylation_unique_Beerman.bed > $out/01_promoter_hypermethylation_Beerman.bed
bedtools intersect -wb -b "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/04_promoter/EPC_promoter_list_sorted.bed" -a 01_hypermethylation_unique_Sun.bed > $out/01_promoter_hypermethylation_Sun.bed

bedtools intersect -wb -b "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/04_promoter/EPC_promoter_list_sorted.bed" -a 02_hypomethylation.bed > $out/02_promoter_hypomethylation_common.bed
bedtools intersect -wb -b "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/04_promoter/EPC_promoter_list_sorted.bed" -a 02_hypomethylation_unique_Beerman.bed > $out/02_promoter_hypomethylation_Beerman.bed
bedtools intersect -wb -b "/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/04_promoter/EPC_promoter_list_sorted.bed" -a 02_hypomethylation_unique_Sun.bed > $out/02_promoter_hypomethylation_Sun.bed

cd $out
wc -l *.bed > $out/num_peaks.txt