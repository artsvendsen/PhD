#!/usr/bash
env="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/05_Studies_overlap"
out="/Users/Art/Drive/PhD/Experiments/Epigenetics/01_bedfiles/03_DNAmethylation/05_Studies_overlap/zz_closest"

bedtools closest -io -d -a $env/01_hypermethylation_unique_Beerman.bed -b $env/01_hypermethylation_unique_Sun.bed > $out/01_closest_io_hyper.bed
bedtools closest -io -d -a $env/02_hypomethylation_unique_Beerman.bed -b $env/02_hypomethylation_unique_Sun.bed > $out/02_closest_io_hypo.bed