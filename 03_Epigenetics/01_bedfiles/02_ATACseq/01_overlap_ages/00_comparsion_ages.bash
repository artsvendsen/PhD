#!/usr/bash
# Coordinate comparison between young and old samples
# A. Svendsen Aug 2021

env="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/00_meged_replicates"
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages"
cd $env

# Broad
sampleName="01_ATAC_broad"
fa="01_ATAC_Young_broad.bed"
fb="02_ATAC_Old_broad.bed"

#jaccard index
bedtools jaccard -a $fa -b $fb >> $out/jaccard_Young_Old_comparison.txt

# isolating age-specific peaks
echo "Generating peak overlap..."
cat $fa $fb > cat_pre_sorted.bed
bedtools sort -i cat_pre_sorted.bed > cat_sorted.bed
bedtools merge -i cat_sorted.bed > reference.bed

bedtools sort -i  reference.bed > reference_sorted.bed
bedtools sort -i $fa > fa_sorted.bed
bedtools sort -i $fb > fb_sorted.bed

bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed  > $out/unique_fb.bed
bedtools intersect -v -a reference_sorted.bed -b fb_sorted.bed  > $out/unique_fa.bed

bedtools subtract -A -a reference_sorted.bed -b $out/unique_fa.bed > temp1.bed
bedtools subtract -A -a temp1.bed -b $out/unique_fb.bed > $out/$sampleName.bed

mv $out/unique_fa.bed $out/"$sampleName"_unique_Young.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_Old.bed

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

# Narrow
sampleName="02_ATAC_narrow"
fa="03_ATAC_Young_narrow.bed"
fb="04_ATAC_Old_narrow.bed"
#jaccard index
bedtools jaccard -a $fa -b $fb >> $out/jaccard_Young_Old_comparison.txt

# isolating age-specific peaks
echo "Generating peak overlap..."
cat $fa $fb > cat_pre_sorted.bed
bedtools sort -i cat_pre_sorted.bed > cat_sorted.bed
bedtools merge -i cat_sorted.bed > reference.bed

bedtools sort -i  reference.bed > reference_sorted.bed
bedtools sort -i $fa > fa_sorted.bed
bedtools sort -i $fb > fb_sorted.bed

bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed  > $out/unique_fb.bed
bedtools intersect -v -a reference_sorted.bed -b fb_sorted.bed  > $out/unique_fa.bed

bedtools subtract -A -a reference_sorted.bed -b $out/unique_fa.bed > temp1.bed
bedtools subtract -A -a temp1.bed -b $out/unique_fb.bed > $out/$sampleName.bed

mv $out/unique_fa.bed $out/"$sampleName"_unique_Young.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_Old.bed

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

cd $out
wc -l $out/*.bed > number_peaks_ages.txt