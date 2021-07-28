#!/bin/bash

split --suffix-length=3 -l 1000 -d /cluster/projects/bhklab/Data/biobank_merged_RNAseq/MuTect1/bed_hg19/S04380110_Covered.headless.bed
for f in x*; do mv "$f" "${f%$}.bed"; done
ls x* > ../bed_list
