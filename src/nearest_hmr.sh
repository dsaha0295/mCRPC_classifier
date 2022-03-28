#!/bin/bash
#Script to get nearest rHMR to TSS of select genes and output to txt file


genes=( CD8A CD8B CD274 )
gencode=/Users/ds/Desktop/projects/data/anno/gencode.v35.gene.sort.gtf
rhmr=/Users/ds/Desktop/projects/data/bs/rHMR_mCRPC_consensus.bed
output=/Users/ds/Desktop/projects/mCRPC_classifier/results/nearest_rhmr.txt


for g in "${genes[@]}"

do

grep -w "$g" "$gencode" | bedtools closest -a stdin -b "$rhmr" | awk -v var="$g" -v OFS="\t" '{print var, $21, $22, $23}' >> "$output" 

done

