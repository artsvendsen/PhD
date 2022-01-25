#!/usr/bash
# Coordinate comparison between H3K4me3 and H3K27me3 from Sara and Sun
# A. Svendsen Aug 2021

env="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/00_merged_replicates"
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/01_overlap_Sara_Sun"
cd $env

# H3K4me3
# Young
sampleName="01_H3K4me3_Young"
fa="01_H3K4me3_Young_Sara_max_merged_replicates.bed"
fb="05_H3K4me3_Young_Sun_max_merged_replicates.bed"

#jaccard index
bedtools jaccard -a $fa -b $fb > $out/jaccard_Sara_Sun_comparison.txt

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

mv $out/unique_fa.bed $out/"$sampleName"_unique_Sara.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_Sun.bed

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

# Old 
sampleName="02_H3K4me3_Old"
fa="02_H3K4me3_Old_Sara_max_merged_replicates.bed"
fb="06_H3K4me3_Old_Sun_max_merged_replicates.bed"

#jaccard index
bedtools jaccard -a $fa -b $fb >> $out/jaccard_Sara_Sun_comparison.txt

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

mv $out/unique_fa.bed $out/"$sampleName"_unique_Sara.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_Sun.bed

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

# H3K27me3
# Young
sampleName="03_H3K27me3_Young"
fa="03_H3K27me3_Young_Sara_max_merged_replicates.bed"
fb="07_H3K27me3_Young_Sun_max_merged_replicates.bed"

#jaccard index
bedtools jaccard -a $fa -b $fb >> $out/jaccard_Sara_Sun_comparison.txt

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

mv $out/unique_fa.bed $out/"$sampleName"_unique_Sara.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_Sun.bed

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

# Old
sampleName="04_H3K27me3_Old"
fa="04_H3K27me3_Old_Sara_max_merged_replicates.bed"
fb="08_H3K27me3_Old_Sun_max_merged_replicates.bed"

#jaccard index
bedtools jaccard -a $fa -b $fb >> $out/jaccard_Sara_Sun_comparison.txt

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

mv $out/unique_fa.bed $out/"$sampleName"_unique_Sara.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_Sun.bed

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed fa_sorted.bed fb_sorted.bed temp1.bed

cd $out
wc -l $out/*.bed > $out/num_peaks_studies.txt