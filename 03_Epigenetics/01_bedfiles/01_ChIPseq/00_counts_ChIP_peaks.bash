#!/usr/bash
# Count number of called peaks of each replicate from Sara and Sun
# A. Svendsen Sept 2021

## Sara's ChIP-seq replicates
env="/Volumes/AFSvendsen/Data_Analysis/20171005_Sara_ChIP-seq/20200708_peakcalling"
out="/Users/Art/Drive/PhD/Experiments/Epigenetics"
cd $env

#K4
fa="3_10_H3K4me3_Y__peaks.narrowPeak"
fb="15_16_H3K4me3_Y__peaks.narrowPeak"
fc="30_31_H3K4me3_Y__peaks.narrowPeak"
wc -l $fa >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fb >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fc >> $out/00_ChIP_peak_number_all_samples.txt

fa="20_H3K4me3_O__peaks.narrowPeak"
fb="21_H3K4me3_O__peaks.narrowPeak"
fc="28_29_H3K4me3_O__peaks.narrowPeak"
wc -l $fa >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fb >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fc >> $out/00_ChIP_peak_number_all_samples.txt

#K27
fa="1_7_8_K27_B6_Y_sorted.rmdup-W200-G600.scoreisland"
fb="13_14_K27_B6_Y_sorted.rmdup-W200-G600.scoreisland"
fc="26_27_K27_B6_Y_sorted.rmdup-W200-G600.scoreisland"
fd="36_K27_B6_Y_sorted.rmdup-W200-G600.scoreisland"
wc -l $fa >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fb >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fc >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fd >> $out/00_ChIP_peak_number_all_samples.txt

fa="19_H3K27me3_O_sorted.rmdup-W200-G600.scoreisland"
fb="24_H3K27me3_O_sorted.rmdup-W200-G600.scoreisland"
fc="25_H3K27me3_O_sorted.rmdup-W200-G600.scoreisland"
wc -l $fa >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fb >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fc >> $out/00_ChIP_peak_number_all_samples.txt

## Sun's ChIP-seq replicates
env="/Users/Art/Drive/PhD/Scripts/Sun_et_al_2014_ChIPseq/04_peak_calling/H3K4me3/H3K4me3_narrow/"
cd $env

#K4
fa="SRR892978_peaks.narrowPeak"
fb="SRR892979_peaks.narrowPeak"
wc -l $fa >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fb >> $out/00_ChIP_peak_number_all_samples.txt

fa="SRR892980_peaks.narrowPeak"
fb="SRR892981_peaks.narrowPeak"
wc -l $fa >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fb >> $out/00_ChIP_peak_number_all_samples.txt

env="/Users/Art/Drive/PhD/Scripts/Sun_et_al_2014_ChIPseq/04_peak_calling/H3K27me3/standard/"
cd $env

#K27
fa="SRR892966_sorted_rmdup-W200-G600.scoreisland"
fb="SRR892967_sorted_rmdup-W200-G600.scoreisland"
fc="SRR892968_sorted_rmdup-W200-G600.scoreisland"
wc -l $fa >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fb >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fc >> $out/00_ChIP_peak_number_all_samples.txt

fa="SRR892969_sorted_rmdup-W200-G600.scoreisland"
fb="SRR892970_sorted_rmdup-W200-G600.scoreisland"
fc="SRR892971_sorted_rmdup-W200-G600.scoreisland"
wc -l $fa >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fb >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fc >> $out/00_ChIP_peak_number_all_samples.txt

env="/Users/Art/Drive/PhD/Scripts/Sun_et_al_2014_ChIPseq/04_peak_calling/H3K36me3/H3K36me3_standard/"
cd $env

#K36
fa="SRR892972_sorted_rmdup-W200-G600.scoreisland"
fb="SRR892973_sorted_rmdup-W200-G600.scoreisland"
fc="SRR892974_sorted_rmdup-W200-G600.scoreisland"
wc -l $fa >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fb >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fc >> $out/00_ChIP_peak_number_all_samples.txt

fa="SRR892975_sorted_rmdup-W200-G600.scoreisland"
fb="SRR892976_sorted_rmdup-W200-G600.scoreisland"
fc="SRR892977_sorted_rmdup-W200-G600.scoreisland"
wc -l $fa >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fb >> $out/00_ChIP_peak_number_all_samples.txt
wc -l $fc >> $out/00_ChIP_peak_number_all_samples.txt