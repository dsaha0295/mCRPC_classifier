#!/bin/bash


filename=$1

while read line;
do
grep -w "$line" gencode.v35.geneid_to_name.gtf | bedtools slop -i stdin -g hg38.chrom.sizes -b 1000000 | bedtools intersect -a mCRPC.rhmr.bed -b stdin -wa -wb | awk -v OFS="\t" -v c="$line" '{print $1, $2, $3, c}' >> DNAm.hmr.AR.sig.txt

done < "$filename"





