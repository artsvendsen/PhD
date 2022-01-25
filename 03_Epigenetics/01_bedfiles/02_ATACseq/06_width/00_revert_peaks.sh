#!/usr/bash
original="/Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling"
reference="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages"
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width"

# get original peaks
# Young common 
bedtools intersect -wb -a $reference/01_ATAC_broad.bed -b $original/"2017_07_04_ATAC_5_3_young_GCTACGCT_p0.01_peaks.broadPeak" > $out/03_Young_common_original_rep1.bed
bedtools intersect -wb -a $reference/01_ATAC_broad.bed -b $original/"2017_07_08_ATAC_5_4_young_AAGAGGCA_p0.01_peaks.broadPeak" > $out/03_Young_common_original_rep2.bed

# Old common
bedtools intersect -wb -a $reference/01_ATAC_broad.bed -b $original/"2017_07_02_ATAC_5_3_old_CGAGGCTG_p0.01_peaks.broadPeak" > $out/03_Old_common_original_rep1.bed
bedtools intersect -wb -a $reference/01_ATAC_broad.bed -b $original/"2017_07_06_ATAC_5_4_old_GTAGAGGA_p0.01_peaks.broadPeak" > $out/03_Old_common_original_rep2.bed

#Young specifc
bedtools intersect -wb -a $reference/01_ATAC_broad_unique_Young.bed -b $original/"2017_07_04_ATAC_5_3_young_GCTACGCT_p0.01_peaks.broadPeak" > $out/03_Young_specific_original_rep1.bed
bedtools intersect -wb -a $reference/01_ATAC_broad_unique_Young.bed -b $original/"2017_07_08_ATAC_5_4_young_AAGAGGCA_p0.01_peaks.broadPeak" > $out/03_Young_specific_original_rep2.bed

#Old specifc 
bedtools intersect -wb -a $reference/01_ATAC_broad_unique_Old.bed -b $original/"2017_07_02_ATAC_5_3_old_CGAGGCTG_p0.01_peaks.broadPeak" > $out/03_Old_specific_original_rep1.bed
bedtools intersect -wb -a $reference/01_ATAC_broad_unique_Old.bed -b $original/"2017_07_06_ATAC_5_4_old_GTAGAGGA_p0.01_peaks.broadPeak" > $out/03_Old_specific_original_rep2.bed