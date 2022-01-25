# install bigwigtowig from conda:
# conda install -c bioconda ucsc-bigwigtowig

# Download Bigwig files from GSE (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41658)

bw="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/02_Taiwo/01_bw_files"
wig="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/02_Taiwo/02_wig_files"
bed="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/02_Taiwo/03_bed_files"

for file in $bw/*.bw
do
	fileName=$(basename $file .bw)
	echo $fileName
	## convert to intermediate wig format (from USCS kent tools)
	bigWigToWig $file $wig/$fileName.wig

	## convert wig to bed (from BEDOPS)
	wig2bed < $wig/$fileName.wig > $bed/$fileName.bed
done


