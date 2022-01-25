#!/usr/bash
# Apply bedtools merge function in order to investigate if there are overlapping peaks
# A. Svendsen Aug 2021

fileFolder="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/"
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/02_ATACseq/01_overlap_ages/merged"

wc -l $fileFolder/*.bed > $out/number_peaks_original.txt

for file in $fileFolder/*.bed
do

	fileName=$(basename $file .bed)
	echo $fileName
	bedtools merge -i $file -d 2000 > $out/$fileName"_merged.bed"
done


wc -l $out/*.bed > $out/number_peaks_merged.txt
