#!/usr/bash
# Coordinate comparison between young and old samples
# A. Svendsen Aug 2021

env="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation"
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/05_Studies_overlap/merged"
cd $out

# Gain of methylation marks with aging (Hypermethylation)
sampleName="01_hypermethylation_merge"
fa="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/01_Beerman/01_raw_data/03_youngtoold_gain_mm10_sorted.bed"
fb="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/03_Sun/01_raw_data/02_hypermethylated_with_aging_mm10_sorted.bed"

# isolating age-specific peaks
echo "Generating peak overlap..."
cat $fa $fb > cat_pre_sorted.bed
cat cat_pre_sorted.bed | cut -f 1-3 > cat_pre_sorted_trimmed.bed
bedtools sort -i cat_pre_sorted_trimmed.bed > cat_sorted.bed
bedtools merge -d 1000 -i cat_sorted.bed > reference.bed

bedtools sort -i  reference.bed > reference_sorted.bed
bedtools sort -i $fa > fa_sorted.bed
bedtools sort -i $fb > fb_sorted.bed

bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed  > $out/unique_fb.bed
bedtools intersect -v -a reference_sorted.bed -b fb_sorted.bed  > $out/unique_fa.bed

bedtools subtract -A -a reference_sorted.bed -b $out/unique_fa.bed > temp1.bed
bedtools subtract -A -a temp1.bed -b $out/unique_fb.bed > $out/$sampleName.bed

mv $out/unique_fa.bed $out/"$sampleName"_unique_Beerman.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_Sun.bed

rm cat_pre_sorted.bed cat_pre_sorted_trimmed.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

# Gain of methylation marks with aging (Hypermethylation)
sampleName="02_hypomethylation_merge"
fa="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/01_Beerman/01_raw_data/03_youngtoold_loss_mm10_sorted.bed"
fb="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/03_Sun/01_raw_data/02_hypomethylated_with_aging_mm10_sorted.bed"

# isolating age-specific peaks
echo "Generating peak overlap..."
cat $fa $fb > cat_pre_sorted.bed
cat cat_pre_sorted.bed | cut -f 1-3 > cat_pre_sorted_trimmed.bed
bedtools sort -i cat_pre_sorted_trimmed.bed > cat_sorted.bed
bedtools merge -d 1000 -i cat_sorted.bed > reference.bed

bedtools sort -i  reference.bed > reference_sorted.bed
bedtools sort -i $fa > fa_sorted.bed
bedtools sort -i $fb > fb_sorted.bed

bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed  > $out/unique_fb.bed
bedtools intersect -v -a reference_sorted.bed -b fb_sorted.bed  > $out/unique_fa.bed

bedtools subtract -A -a reference_sorted.bed -b $out/unique_fa.bed > temp1.bed
bedtools subtract -A -a temp1.bed -b $out/unique_fb.bed > $out/$sampleName.bed

mv $out/unique_fa.bed $out/"$sampleName"_unique_Beerman.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_Sun.bed

rm cat_pre_sorted.bed cat_pre_sorted_trimmed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

cd $out
wc -l $out/*.bed > number_peaks_studies.txt