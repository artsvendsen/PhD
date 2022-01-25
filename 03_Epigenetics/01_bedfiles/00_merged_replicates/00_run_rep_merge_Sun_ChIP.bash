#!/usr/bash
# Combine ChIP-seq peak called files; peaks that are detected at least 2 out of 3 times will be kept.
# The resulting peaks are wider than in individual replicates (the maximum area encompassing identified peak regions (“MAX”)
# see http://dx.doi.org/10.5936/csbj.201401002)
# A. Svendsen Aug 2021

## Sun's ChIP-seq replicates
env="/Users/Art/Drive/PhD/Scripts/Sun_et_al_2014_ChIPseq/04_peak_calling/H3K4me3/H3K4me3_narrow/"
cd $env
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq"

#K4 - Young samples
sampleName="05_H3K4me3_Young_Sun_max_merged_replicates.bed"
fa="SRR892978_peaks.narrowPeak"
fb="SRR892979_peaks.narrowPeak"

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
bedtools subtract -A -a temp1.bed -b $out/unique_fb.bed > $out/$sampleName

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed $out/unique_fb.bed $out/unique_fa.bed

#K4 - Old samples
sampleName="06_H3K4me3_Old_Sun_max_merged_replicates.bed"
fa="SRR892980_peaks.narrowPeak"
fb="SRR892981_peaks.narrowPeak"

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
bedtools subtract -A -a temp1.bed -b $out/unique_fb.bed > $out/$sampleName

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed $out/unique_fb.bed $out/unique_fa.bed

env="/Users/Art/Drive/PhD/Scripts/Sun_et_al_2014_ChIPseq/04_peak_calling/H3K27me3/standard/"
cd $env
#K27 - Young samples
sampleName="07_H3K27me3_Young_Sun_max_merged_replicates.bed"
fa="SRR892966_sorted_rmdup-W200-G600.scoreisland"
fb="SRR892967_sorted_rmdup-W200-G600.scoreisland"
fc="SRR892968_sorted_rmdup-W200-G600.scoreisland"

echo "Generating peak overlap..."
cat $fa $fb $fc > cat_pre_sorted.bed
bedtools sort -i cat_pre_sorted.bed > cat_sorted.bed
bedtools merge -i cat_sorted.bed > reference.bed

bedtools sort -i  reference.bed > reference_sorted.bed
bedtools sort -i $fa > fa_sorted.bed
bedtools sort -i $fb > fb_sorted.bed
bedtools sort -i $fc > fc_sorted.bed

bedtools intersect -v -a reference_sorted.bed -b fb_sorted.bed fc_sorted.bed > $out/unique_fa.bed
bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed fc_sorted.bed > $out/unique_fb.bed
bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed fb_sorted.bed > $out/unique_fc.bed

bedtools subtract -A -a reference_sorted.bed -b $out/unique_fa.bed > temp1.bed
bedtools subtract -A -a reference_sorted.bed -b $out/unique_fb.bed > temp2.bed
bedtools subtract -A -a temp2.bed -b $out/unique_fc.bed > $out/$sampleName

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed fa_sorted.bed fb_sorted.bed fc_sorted.bed temp1.bed temp2.bed $out/unique_fa.bed $out/unique_fb.bed $out/unique_fc.bed

#K27 - Old samples
sampleName="08_H3K27me3_Old_Sun_max_merged_replicates.bed"
fa="SRR892969_sorted_rmdup-W200-G600.scoreisland"
fb="SRR892970_sorted_rmdup-W200-G600.scoreisland"
fc="SRR892971_sorted_rmdup-W200-G600.scoreisland"

echo "Generating peak overlap..."
cat $fa $fb $fc > cat_pre_sorted.bed
bedtools sort -i cat_pre_sorted.bed > cat_sorted.bed
bedtools merge -i cat_sorted.bed > reference.bed

bedtools sort -i  reference.bed > reference_sorted.bed
bedtools sort -i $fa > fa_sorted.bed
bedtools sort -i $fb > fb_sorted.bed
bedtools sort -i $fc > fc_sorted.bed

bedtools intersect -v -a reference_sorted.bed -b fb_sorted.bed fc_sorted.bed > $out/unique_fa.bed
bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed fc_sorted.bed > $out/unique_fb.bed
bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed fb_sorted.bed > $out/unique_fc.bed

bedtools subtract -A -a reference_sorted.bed -b $out/unique_fa.bed > temp1.bed
bedtools subtract -A -a reference_sorted.bed -b $out/unique_fb.bed > temp2.bed
bedtools subtract -A -a temp2.bed -b $out/unique_fc.bed > $out/$sampleName

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed fa_sorted.bed fb_sorted.bed fc_sorted.bed temp1.bed temp2.bed $out/unique_fa.bed $out/unique_fb.bed $out/unique_fc.bed

env="/Users/Art/Drive/PhD/Scripts/Sun_et_al_2014_ChIPseq/04_peak_calling/H3K36me3/H3K36me3_standard/"
cd $env
#K36 - Young samples
sampleName="09_H3K36me3_Young_Sun_max_merged_replicates.bed"
fa="SRR892972_sorted_rmdup-W200-G600.scoreisland"
fb="SRR892973_sorted_rmdup-W200-G600.scoreisland"
fc="SRR892974_sorted_rmdup-W200-G600.scoreisland"

echo "Generating peak overlap..."
cat $fa $fb $fc > cat_pre_sorted.bed
bedtools sort -i cat_pre_sorted.bed > cat_sorted.bed
bedtools merge -i cat_sorted.bed > reference.bed

bedtools sort -i  reference.bed > reference_sorted.bed
bedtools sort -i $fa > fa_sorted.bed
bedtools sort -i $fb > fb_sorted.bed
bedtools sort -i $fc > fc_sorted.bed

bedtools intersect -v -a reference_sorted.bed -b fb_sorted.bed fc_sorted.bed > $out/unique_fa.bed
bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed fc_sorted.bed > $out/unique_fb.bed
bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed fb_sorted.bed > $out/unique_fc.bed

bedtools subtract -A -a reference_sorted.bed -b $out/unique_fa.bed > temp1.bed
bedtools subtract -A -a reference_sorted.bed -b $out/unique_fb.bed > temp2.bed
bedtools subtract -A -a temp2.bed -b $out/unique_fc.bed > $out/$sampleName

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed fa_sorted.bed fb_sorted.bed fc_sorted.bed temp1.bed temp2.bed $out/unique_fa.bed $out/unique_fb.bed $out/unique_fc.bed

#K36 - Old samples
sampleName="10_H3K36me3_Old_Sun_max_merged_replicates.bed"
fa="SRR892975_sorted_rmdup-W200-G600.scoreisland"
fb="SRR892976_sorted_rmdup-W200-G600.scoreisland"
fc="SRR892977_sorted_rmdup-W200-G600.scoreisland"

echo "Generating peak overlap..."
cat $fa $fb $fc > cat_pre_sorted.bed
bedtools sort -i cat_pre_sorted.bed > cat_sorted.bed
bedtools merge -i cat_sorted.bed > reference.bed

bedtools sort -i  reference.bed > reference_sorted.bed
bedtools sort -i $fa > fa_sorted.bed
bedtools sort -i $fb > fb_sorted.bed
bedtools sort -i $fc > fc_sorted.bed

bedtools intersect -v -a reference_sorted.bed -b fb_sorted.bed fc_sorted.bed > $out/unique_fa.bed
bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed fc_sorted.bed > $out/unique_fb.bed
bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed fb_sorted.bed > $out/unique_fc.bed

bedtools subtract -A -a reference_sorted.bed -b $out/unique_fa.bed > temp1.bed
bedtools subtract -A -a reference_sorted.bed -b $out/unique_fb.bed > temp2.bed
bedtools subtract -A -a temp2.bed -b $out/unique_fc.bed > $out/$sampleName

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed fa_sorted.bed fb_sorted.bed fc_sorted.bed temp1.bed temp2.bed $out/unique_fa.bed $out/unique_fb.bed $out/unique_fc.bed