#!/usr/bash
# Combine ChIP-seq peak called files; peaks that are detected at least 2 out of 3 times will be kept.
# The resulting peaks are wider than in individual replicates (the maximum area encompassing identified peak regions (“MAX”)
# see http://dx.doi.org/10.5936/csbj.201401002)
# A. Svendsen Aug 2021

## Sara's ChIP-seq replicates
env="/Volumes/AFSvendsen/Data_Analysis/20171005_Sara_ChIP-seq/20200708_peakcalling"
cd $env
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq"

#K4 - Young samples
sampleName="01_H3K4me3_Young_Sara_max_merged_replicates.bed"
fa="3_10_H3K4me3_Y__peaks.narrowPeak"
fb="15_16_H3K4me3_Y__peaks.narrowPeak"
fc="30_31_H3K4me3_Y__peaks.narrowPeak"

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

#K4 - Old samples
sampleName="02_H3K4me3_Old_Sara_max_merged_replicates.bed"
fa="20_H3K4me3_O__peaks.narrowPeak"
fb="21_H3K4me3_O__peaks.narrowPeak"
fc="28_29_H3K4me3_O__peaks.narrowPeak"

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

#K27 - Young samples
sampleName="03_H3K27me3_Young_Sara_max_merged_replicates.bed"
fa="1_7_8_K27_B6_Y_sorted.rmdup-W200-G600.scoreisland"
fb="13_14_K27_B6_Y_sorted.rmdup-W200-G600.scoreisland"
fc="26_27_K27_B6_Y_sorted.rmdup-W200-G600.scoreisland"
fd="36_K27_B6_Y_sorted.rmdup-W200-G600.scoreisland"

echo "Generating peak overlap..."
cat $fa $fb $fc $fd > cat_pre_sorted.bed
bedtools sort -i cat_pre_sorted.bed > cat_sorted.bed
bedtools merge -i cat_sorted.bed > reference.bed

bedtools sort -i  reference.bed > reference_sorted.bed
bedtools sort -i $fa > fa_sorted.bed
bedtools sort -i $fb > fb_sorted.bed
bedtools sort -i $fc > fc_sorted.bed
bedtools sort -i $fd > fd_sorted.bed

bedtools intersect -v -a reference_sorted.bed -b fb_sorted.bed fc_sorted.bed fd_sorted.bed > $out/unique_fa.bed
bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed fc_sorted.bed fd_sorted.bed > $out/unique_fb.bed
bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed fb_sorted.bed fd_sorted.bed > $out/unique_fc.bed
bedtools intersect -v -a reference_sorted.bed -b fa_sorted.bed fb_sorted.bed fc_sorted.bed > $out/unique_fd.bed

bedtools subtract -A -a reference_sorted.bed -b $out/unique_fa.bed > temp1.bed
bedtools subtract -A -a reference_sorted.bed -b $out/unique_fb.bed > temp2.bed
bedtools subtract -A -a reference_sorted.bed -b $out/unique_fc.bed > temp3.bed
bedtools subtract -A -a temp3.bed -b $out/unique_fd.bed > $out/$sampleName

rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed fa_sorted.bed fb_sorted.bed fc_sorted.bed fd_sorted.bed temp1.bed temp2.bed temp3.bed $out/unique_fa.bed $out/unique_fb.bed $out/unique_fc.bed $out/unique_fd.bed


#K27 - Old samples
sampleName="03_H3K27me3_Old_Sara_max_merged_replicates.bed"
fa="19_H3K27me3_O_sorted.rmdup-W200-G600.scoreisland"
fb="24_H3K27me3_O_sorted.rmdup-W200-G600.scoreisland"
fc="25_H3K27me3_O_sorted.rmdup-W200-G600.scoreisland"

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