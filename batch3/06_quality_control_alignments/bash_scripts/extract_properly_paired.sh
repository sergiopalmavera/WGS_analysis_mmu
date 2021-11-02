#!/bin/bash

# This script extracts the pct of properly aligned reads from samtool's flagstat report

for i in ../output/*.flagstat.tab
do
	pct=$(grep "properly paired"  $i | cut -d'(' -f2 | cut -d'%' -f1)
	sample_name=$(echo $i | cut -d'-' -f1 | cut -d'/' -f3)
	echo $sample_name $pct
done > ../output/flagstat_pct_properly_paired_reads.tab
