#!/usr/bash
# Determining overlapping K4 and K27 peaks (bivalent) in young and aged sets
# A. Svendsen Sept 2021
overlap_peaks () {
	echo "#################################################"
	echo "New Sample Name is $1"
	echo "Input #1 (ChIP) is $2"
	echo "Input #2 (ATAC) is $3"
	echo "Output is $4"
	name_fa="K4"
	name_fb="K27"

	# isolating age-specific peaks
	echo "###########################"
	echo "1.merging..."
	cat $fa $fb > cat_pre_sorted.bed
	echo "2.sorting merged..."
	bedtools sort -i cat_pre_sorted.bed > cat_sorted.bed
	echo "3.building reference file..."
	bedtools merge -i cat_sorted.bed > reference.bed

	echo "4.sorting reference file..."
	bedtools sort -i  reference.bed > reference_sorted.bed
	echo "5.sorting $2 ..."
	bedtools sort -i $fa > fa_sorted.bed
	echo "6.sorting $3 ..."
	bedtools sort -i $fb > fb_sorted.bed

	echo "7.generating unique peaks in $3 ..."
	bedtools intersect -nonamecheck -v -a reference_sorted.bed -b fa_sorted.bed  > $out/unique_fb.bed
	echo "8.generating unique peaks in $2 ..."
	bedtools intersect -nonamecheck -v -a reference_sorted.bed -b fb_sorted.bed  > $out/unique_fa.bed

	echo "9.generating common peaks..."
	bedtools subtract -nonamecheck -A -a reference_sorted.bed -b $out/unique_fa.bed > temp1.bed
	bedtools subtract -nonamecheck -A -a temp1.bed -b $out/unique_fb.bed > $out/$sampleName.bed

	mv $out/unique_fa.bed $out/"$sampleName"_unique_$name_fa.bed
	mv $out/unique_fb.bed $out/"$sampleName"_unique_$name_fb.bed

	echo "10. cleaning up..."
	rm cat_pre_sorted.bed cat_sorted.bed reference.bed reference_sorted.bed  fa_sorted.bed fb_sorted.bed temp1.bed

	echo "done!"
	echo "#################################################"
}

#Young
path="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/01_overlap_Sara_Sun"
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/01_ChIPseq/03_bivalent_peaks"
sampleName="01_bivalent_young"
fa="$path/01_H3K4me3_Young.bed" 
fb="$path/03_H3K27me3_Young.bed"
overlap_peaks $sampleName $fa $fb $out

#Old
sampleName="02_bivalent_old"
fa="$path/02_H3K4me3_Old.bed" 
fb="$path/04_H3K27me3_Old.bed"
overlap_peaks $sampleName $fa $fb $out