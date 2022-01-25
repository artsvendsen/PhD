#!/usr/bash
original="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width"
enhancer="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/03_enhancer/amit_enhancer_list_mm10.bed"
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/06_width/enhancers"

cd $original

bedtools intersect -wb -a $enhancer -b 03_Young_common_original_rep1.bed > $out/01_enhancer_Young_common_original_rep1.bed
bedtools intersect -wb -a $enhancer -b 03_Young_common_original_rep2.bed > $out/01_enhancer_Young_common_original_rep2.bed

bedtools intersect -wb -a $enhancer -b 03_Old_common_original_rep1.bed > $out/01_enhancer_Old_common_original_rep1.bed
bedtools intersect -wb -a $enhancer -b 03_Old_common_original_rep2.bed > $out/01_enhancer_Old_common_original_rep2.bed

bedtools intersect -wb -a $enhancer -b 03_Young_specific_original_rep1.bed > $out/01_enhancer_Young_specific_original_rep1.bed
bedtools intersect -wb -a $enhancer -b 03_Young_specific_original_rep2.bed > $out/01_enhancer_Young_specific_original_rep2.bed

bedtools intersect -wb -a $enhancer -b 03_Old_specific_original_rep1.bed > $out/01_enhancer_Old_specific_original_rep1.bed
bedtools intersect -wb -a $enhancer -b 03_Old_specific_original_rep1.bed > $out/01_enhancer_Old_specific_original_rep2.bed