#!/usr/bash

#Merge peaks between biolocial replicates
#test showed that in order of the files (a or b/b or a) does not interfere with the overlap result, so just use it as is.

#Broad peaks
#Old
#bedtools intersect -a /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_02_ATAC_5_3_old_CGAGGCTG_p0.01_peaks.broadPeak -b /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_06_ATAC_5_4_old_GTAGAGGA_p0.01_peaks.broadPeak | bedtools sort -i stdin > broad_overlap_old_0.01.bed
#Young
#bedtools intersect -a /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_04_ATAC_5_3_young_GCTACGCT_p0.01_peaks.broadPeak -b /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_08_ATAC_5_4_young_AAGAGGCA_p0.01_peaks.broadPeak | bedtools sort -i stdin > broad_overlap_young_0.01.bed
#Overlap Young and Old
#bedtools intersect -v -a broad_overlap_old_0.01.bed -b broad_overlap_young_0.01.bed > broad_unique_old_p0.01.bed
#bedtools intersect -v -b broad_overlap_old_0.01.bed -a broad_overlap_young_0.01.bed > broad_unique_young_p0.01.bed
#bedtools intersect -wa -a broad_overlap_old_0.01.bed -b broad_overlap_young_0.01.bed > broad_overlap_p0.01.bed

# wc -l broad_overlap_old_0.01.bed
# wc -l broad_overlap_young_0.01.bed
# wc -l broad_unique_old_p0.01.bed
# wc -l broad_unique_young_p0.01.bed
# wc -l broad_overlap_p0.01.bed

#Narrow peaks
#Old
bedtools intersect -a /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_02_ATAC_5_3_old_CGAGGCTG_p0.01_peaks.narrowPeak -b /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_06_ATAC_5_4_old_GTAGAGGA_p0.01_peaks.narrowPeak | bedtools sort -i stdin | bedtools merge -i stdin > narrow_overlap_old_0.01.bed
#Young
bedtools intersect -a /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_04_ATAC_5_3_young_GCTACGCT_p0.01_peaks.narrowPeak -b /Volumes/AFSvendsen/Data_Analysis/20171012_young_old_LTHSC/20200713_PeakCalling/2017_07_08_ATAC_5_4_young_AAGAGGCA_p0.01_peaks.narrowPeak | bedtools sort -i stdin | bedtools merge -i stdin > narrow_overlap_young_0.01.bed
#Overlap Young and Old
bedtools intersect -v -a narrow_overlap_old_0.01.bed -b narrow_overlap_young_0.01.bed > narrow_unique_old_p0.01.bed
bedtools intersect -v -b narrow_overlap_old_0.01.bed -a narrow_overlap_young_0.01.bed > narrow_unique_young_p0.01.bed
bedtools intersect -wa -a narrow_overlap_old_0.01.bed -b narrow_overlap_young_0.01.bed > narrow_overlap_p0.01.bed

wc -l narrow_overlap_old_0.01.bed
wc -l narrow_overlap_young_0.01.bed
wc -l narrow_unique_old_p0.01.bed
wc -l narrow_unique_young_p0.01.bed
wc -l narrow_overlap_p0.01.bed