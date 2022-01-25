#!/usr/bash
# Combine ATAC-seq peak called files
# The resulting peaks are wider than in individual replicates (the maximum area encompassing identified peak regions (“MAX”)
# see http://dx.doi.org/10.5936/csbj.201401002)
# A. Svendsen Aug 2021


env="/Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/"
cd $env
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/00_meged_replicates/"

## broad peaks
## young
sampleName="01_ATAC_Young_broad"
fa="2017_07_04_ATAC_5_3_young_GCTACGCT_p0.01_peaks.broadPeak"
fb="2017_07_08_ATAC_5_4_young_AAGAGGCA_p0.01_peaks.broadPeak"

#jaccard index
bedtools jaccard -a $fa -b $fb > $out/jaccard_Reps_comparison.txt

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

mv $out/unique_fa.bed $out/"$sampleName"_unique_rep1.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_rep2.bed

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

## old
sampleName="02_ATAC_Old_broad"
fa="2017_07_02_ATAC_5_3_old_CGAGGCTG_p0.01_peaks.broadPeak"
fb="2017_07_06_ATAC_5_4_old_GTAGAGGA_p0.01_peaks.broadPeak"

bedtools jaccard -a $fa -b $fb > $out/jaccard_Reps_comparison.txt

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

mv $out/unique_fa.bed $out/"$sampleName"_unique_rep1.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_rep2.bed

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

## Narrow
## Young
sampleName="03_ATAC_Young_narrow"
fa="2017_07_04_ATAC_5_3_young_GCTACGCT_p0.01_peaks.narrowPeak"
fb="2017_07_08_ATAC_5_4_young_AAGAGGCA_p0.01_peaks.narrowPeak"

bedtools jaccard -a $fa -b $fb > $out/jaccard_Reps_comparison.txt

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

mv $out/unique_fa.bed $out/"$sampleName"_unique_rep1.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_rep2.bed

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

## Old
sampleName="04_ATAC_Old_narrow"
fa="2017_07_02_ATAC_5_3_old_CGAGGCTG_p0.01_peaks.narrowPeak"
fb="2017_07_06_ATAC_5_4_old_GTAGAGGA_p0.01_peaks.narrowPeak"

bedtools jaccard -a $fa -b $fb > $out/jaccard_Reps_comparison.txt

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

mv $out/unique_fa.bed $out/"$sampleName"_unique_rep1.bed
mv $out/unique_fb.bed $out/"$sampleName"_unique_rep2.bed

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

cd $out
wc -l $out/*.bed > $out/num_peaks_reps.txt