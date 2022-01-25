#!/usr/bash
# Coordinate comparison between young and old samples
# A. Svendsen Aug 2021

env="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation"
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/05_Studies_overlap/v2"
cd $env

# Gain of methylation marks with aging (Hypermethylation)
sampleName="01_hypermethylation"
fa="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/01_Beerman/01_raw_data/03_youngtoold_gain_mm10_sorted.bed"
fb="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/03_Sun/01_raw_data/02_Sun_hypermethylated_mm10.bed"

bedtools sort -i $fa > fa_sorted.bed
bedtools sort -i $fb > fb_sorted.bed

bedtools intersect -a fa_sorted.bed -b fb_sorted.bed > $out/$sampleName.bed

# bedtools intersect -v -a fa_sorted.bed -b $out/$sampleName.bed  > $out/$sampleName"_unique_Beerman.bed"
# bedtools intersect -v -a fb_sorted.bed -b $out/$sampleName.bed  > $out/$sampleName"_unique_Sun.bed"

bedtools subtract -A -a fa_sorted.bed -b $out/$sampleName.bed  > $out/$sampleName"_unique_Beerman.bed"
bedtools subtract -A -a fb_sorted.bed -b $out/$sampleName.bed  > $out/$sampleName"_unique_Sun.bed"

rm fa_sorted.bed fb_sorted.bed

# Gain of methylation marks with aging (Hypermethylation)
sampleName="02_hypomethylation"
fa="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/01_Beerman/01_raw_data/03_youngtoold_loss_mm10_sorted.bed"
fb="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/03_Sun/01_raw_data/02_Sun_hypomethylated_mm10.bed"

bedtools sort -i $fa > fa_sorted.bed
bedtools sort -i $fb > fb_sorted.bed

bedtools intersect -a fa_sorted.bed -b fb_sorted.bed > $out/$sampleName.bed

bedtools subtract -A -a fa_sorted.bed -b $out/$sampleName.bed  > $out/$sampleName"_unique_Beerman.bed"
bedtools subtract -A -a fb_sorted.bed -b $out/$sampleName.bed  > $out/$sampleName"_unique_Sun.bed" 

rm fa_sorted.bed fb_sorted.bed

cd $out
wc -l $out/*.bed > number_peaks_studies.txt