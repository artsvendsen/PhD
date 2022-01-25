#!/usr/bash
# Check the chromatin status of the different histone marks in young, aged and common peaks. i.e "Are K4 peaks really open?"

chip="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/02_overlap_ages"
atac="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages"

# K4
# young-specific peaks
bedtools intersect -u -b $atac/01_ATAC_broad_unique_Young.bed -a $chip/01_H3K4me3_unique_Young.bed | bedtools sort -i stdin | bedtools merge -i stdin > 01_ATAC_ChIP_K04_young_overlap.bed
wc -l 01_ATAC_ChIP_K04_young_overlap.bed > 04_number_ChIP_open_sites.txt
# aged-specific peaks
bedtools intersect -u -b $atac/01_ATAC_broad_unique_Old.bed -a $chip/01_H3K4me3_unique_Old.bed | bedtools sort -i stdin | bedtools merge -i stdin > 01_ATAC_ChIP_K04_aged_overlap.bed
wc -l 01_ATAC_ChIP_K04_aged_overlap.bed >> 04_number_ChIP_open_sites.txt
# common peaks
bedtools intersect -u -b $atac/01_ATAC_broad.bed -a $chip/01_H3K4me3.bed | bedtools sort -i stdin | bedtools merge -i stdin > 01_ATAC_ChIP_K04_common_overlap.bed
wc -l 01_ATAC_ChIP_K04_common_overlap.bed >> 04_number_ChIP_open_sites.txt

# K27
# young-specific peaks
bedtools intersect -u -b $atac/01_ATAC_broad_unique_Young.bed -a $chip/02_H3K27me3_unique_Young.bed | bedtools sort -i stdin | bedtools merge -i stdin > 02_ATAC_ChIP_K27_young_overlap.bed
wc -l 02_ATAC_ChIP_K27_young_overlap.bed >> 04_number_ChIP_open_sites.txt
# aged-specific peaks
bedtools intersect -u -b $atac/01_ATAC_broad_unique_Old.bed -a $chip/02_H3K27me3_unique_Old.bed | bedtools sort -i stdin | bedtools merge -i stdin > 02_ATAC_ChIP_K27_aged_overlap.bed
wc -l 02_ATAC_ChIP_K27_aged_overlap.bed >> 04_number_ChIP_open_sites.txt
# common peaks
bedtools intersect -u -b $atac/01_ATAC_broad.bed -a $chip/02_H3K27me3.bed | bedtools sort -i stdin | bedtools merge -i stdin > 02_ATAC_ChIP_K27_common_overlap.bed
wc -l 02_ATAC_ChIP_K27_common_overlap.bed >> 04_number_ChIP_open_sites.txt

# K36
# young-specific peaks
bedtools intersect -u -b $atac/01_ATAC_broad_unique_Young.bed -a $chip/03_H3K36me3_unique_Young.bed | bedtools sort -i stdin | bedtools merge -i stdin > 03_ATAC_ChIP_K36_young_overlap.bed
wc -l 03_ATAC_ChIP_K36_young_overlap.bed >> 04_number_ChIP_open_sites.txt
# aged-specific peaks
bedtools intersect -u -b $atac/01_ATAC_broad_unique_Old.bed -a $chip/03_H3K36me3_unique_Old.bed | bedtools sort -i stdin | bedtools merge -i stdin > 03_ATAC_ChIP_K36_aged_overlap.bed
wc -l 03_ATAC_ChIP_K36_aged_overlap.bed >> 04_number_ChIP_open_sites.txt
# common peaks
bedtools intersect -u -b $atac/01_ATAC_broad.bed -a $chip/03_H3K36me3.bed | bedtools sort -i stdin | bedtools merge -i stdin > 03_ATAC_ChIP_K36_common_overlap.bed
wc -l 03_ATAC_ChIP_K36_common_overlap.bed >> 04_number_ChIP_open_sites.txt

wc -l $chip/01_H3K4me3_unique_Young.bed > 05_number_total_ChIP_peaks.txt
wc -l $chip/01_H3K4me3_unique_Old.bed >> 05_number_total_ChIP_peaks.txt
wc -l $chip/01_H3K4me3.bed >> 05_number_total_ChIP_peaks.txt
wc -l $chip/02_H3K27me3_unique_Young.bed >> 05_number_total_ChIP_peaks.txt
wc -l $chip/02_H3K27me3_unique_Old.bed >> 05_number_total_ChIP_peaks.txt
wc -l $chip/02_H3K27me3.bed >> 05_number_total_ChIP_peaks.txt
wc -l $chip/03_H3K36me3_unique_Young.bed >> 05_number_total_ChIP_peaks.txt
wc -l $chip/03_H3K36me3_unique_Old.bed >> 05_number_total_ChIP_peaks.txt
wc -l $chip/03_H3K36me3.bed >> 05_number_total_ChIP_peaks.txt
