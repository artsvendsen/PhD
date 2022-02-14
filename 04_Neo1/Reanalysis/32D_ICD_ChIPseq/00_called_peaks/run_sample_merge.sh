#!/usr/bash
# Combine ChIP-seq peak called files; peaks that are detected at least 2 out of 3 times will be kept.
# The resulting peaks are wider than in individual replicates (the maximum area encompassing identified peak regions (“MAX”)
# see http://dx.doi.org/10.5936/csbj.201401002)
# A. Svendsen Dez 2021

out="/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/00_called_peaks/"

#Neo1 Broad Peaks
env="/Users/Art/Drive/PhD/Experiments/Neogenin1/zz_Thesis_Chapter/Reanalysis/32D_ICD_ChIPseq/00_called_peaks/00_filtered_peaks/"
cd $env

sampleName="ChIP_Neo1_broad_max_merged_replicates.bed"
fa="ChIP_Neo1_1_S1_pval0.0001.broadPeak"
fb="ChIP_Neo1_2_S4_pval0.0001.broadPeak"
fc="ChIP_Neo1_3_S7_pval0.0001.broadPeak"

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


#Neo1 Narrow Peaks
sampleName="ChIP_Neo1_narrow_max_merged_replicates.bed"
fa="ChIP_Neo1_1_S1_pval0.0001.narrowPeak"
fb="ChIP_Neo1_2_S4_pval0.0001.narrowPeak"
fc="ChIP_Neo1_3_S7_pval0.0001.narrowPeak"

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