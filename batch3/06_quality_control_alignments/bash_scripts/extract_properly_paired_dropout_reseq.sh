#!/bin/bash

# This script extracts the pct of properly aligned reads from samtool's flagstat report

i=../output/I34772-L1_S19_L004.sorted.RG.dedup.bqsr.flagstat.tab
pct=$(grep "properly paired"  $i | cut -d'(' -f2 | cut -d'%' -f1)
sample_name=$(echo $i | cut -d'-' -f1 | cut -d'/' -f3)
echo ${sample_name}.dropout_reseqd $pct > tmp.txt #output line with dropout-reseqd sample
	
cat ../output/flagstat_pct_properly_paired_reads.tab ./tmp.txt > ./tmp2.txt # concatenate old list with dropout-reseqd sample

mv ./tmp2.txt ../output/flagstat_pct_properly_paired_reads.tab # move new list into the same name as old list

rm ./tmp.txt # remove temp files
rm ./tmp2.txt

