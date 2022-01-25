#!/usr/bash
#Identify peaks which are overlaping with exon, intergenic, both or no annotated reagions
#A. Svendsen Dez 2020

bedtools sort -i 00_ATACseq_narrow_p0_01_DBP_diffbind.bed > 00_ATACseq_narrow_p0_01_DBP_diffbind_sorted.bed
#Run experimental peaks against annotation references (generated in ger_annotation.bash)
## overlaping peaks with exon annotation
bedtools intersect -b /Users/Art/Desktop/M13_annotation/m13_exon_formatted.bed -a 00_ATACseq_narrow_p0_01_DBP_diffbind_sorted.bed -wo > 03_overlap_narrow_exons_wo.bed
bedtools merge -i 03_overlap_narrow_exons_wo.bed > 03_overlap_narrow_exons_merged.bed

## non-overlapping peaks with exon annotation
bedtools subtract -a 00_ATACseq_narrow_p0_01_DBP_diffbind_sorted.bed -b 03_overlap_narrow_exons_wo.bed -A > 03_nonoverlap_narrow_exons.bed

## remaning peaks with intragenic annotation
bedtools intersect -a 03_nonoverlap_narrow_exons.bed -b /Users/Art/Desktop/M13_annotation/m13_intragenic.bed -wo > 03_overlap_narrow_intragenic_wo.bed
bedtools merge -i 03_overlap_narrow_intragenic_wo.bed > 03_overlap_narrow_intragenic_merged.bed
# intergenic peaks
bedtools subtract -a 03_nonoverlap_narrow_exons.bed -b 03_overlap_narrow_intragenic_wo.bed > 03_intergenic_narrow.bed

#Exons
bedtools intersect -b /Users/Art/Desktop/M13_annotation/m13_exon_formatted.bed -a 00_ATACseq_narrow_p0_01_DBP_diffbind_sorted.bed -wo > 04_overlap_narrow_exons_wo.bed

#Intragenic
bedtools intersect -b /Users/Art/Desktop/M13_annotation/m13_intragenic.bed -a 00_ATACseq_narrow_p0_01_DBP_diffbind_sorted.bed -wo > 04_overlap_narrow_intragenic_wo.bed

#Intergenic
bedtools subtract -b 04_overlap_narrow_exons_wo.bed -a /Users/Art/Drive/PhD/Experiments/Aging\ Signature/000_FINAL/07_Transcriptional_activation/ATAC-seq/02_DBPs_annotation/different_annotations/narrow/00_ATACseq_narrow_p0_01_DBP_diffbind_sorted.bed -A | bedtools subtract -b 04_overlap_narrow_exons_wo.bed -a stdin -A > 04_overlap_narrow_intrergenic_wo.bed

   #  4932 00_ATACseq_narrow_p0_01.bed
   #  4932 00_ATACseq_narrow_p0_01_DBP_diffbind.bed
   #  4932 00_ATACseq_narrow_p0_01_DBP_diffbind_sorted.bed
   #   573 03_intergenic_narrow.bed
   #  1218 03_nonoverlap_narrow_exons.bed
   #  3714 03_overlap_narrow_exons_merged.bed
   # 16456 03_overlap_narrow_exons_wo.bed
   #   645 03_overlap_narrow_intragenic_merged.bed
   #   672 03_overlap_narrow_intragenic_wo.bed
   # 16456 04_overlap_narrow_exons_wo.bed
   #  6411 04_overlap_narrow_intragenic_wo.bed
   #  1218 04_overlap_narrow_intrergenic_wo.bed
   # 62159 total