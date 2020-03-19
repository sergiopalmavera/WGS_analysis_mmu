#!/bin/bash

for fl in ../output/I34772-L1_S19_L004.sorted.RG.dedup.bqsr.CollectWgsMetrics.*
do
	echo "# $fl"
	ln_pt1=$(grep -n "## METRICS CLASS" $fl | cut -d':' -f1)
	ln_pt2=$(grep -n "## HISTOGRAM" $fl | cut -d':' -f1)
	last_lin_pt1=$(( ln_pt2 -2 ))

	sed -n "$ln_pt1,${last_lin_pt1}p" $fl > ${fl/.txt/_pt1.txt}

	sed -n "$ln_pt2,$ p" $fl > ${fl/.txt/_pt2.txt}
	printf "\n"
done
