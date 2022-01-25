#!/usr/shell
env="/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/00_sourceFiles/Uniprot/"
genes="03_Aging_Signature_REAN_UniprotEntries.tab"


while read -r line
 do
    gene=`echo $line | cut -d " " -f 2`
    entry=`echo $line | cut -d " " -f 1`
    echo $gene
    echo $entry
    entry=${entry%$'\r'}
    #Get GOs associated with protein
    #curl "https://www.uniprot.org/uniprot/?sort=score&desc=&query=accession:$entry&format=tab&columns=id,genes(PREFERRED),go" >> Uniprot_query_GO.txt
    #Get locatilzation associated with protein
    curl "https://www.uniprot.org/uniprot/?sort=score&desc=&query=accession:$entry&format=tab&columns=id,genes(PREFERRED),comment%28SUBCELLULAR%20LOCATION%29" >> Uniprot_query_location.txt
	sleep 0.5
	echo " "

done < "$env/$genes"

#sed -n ’n;p’ Uniprot_query_GO.txt > Uniprot_query_GO_filtered.txt
sed -n ’n;p’ Uniprot_query_location.txt > Uniprot_query_location_filtered.txt