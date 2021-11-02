#!/bin/bash

# This script changes commas to dots in the first part of CollectWgsMetrics report

for fl in ../output/*_pt1.txt
do
	echo $fl
	sed 's/,/./g' $fl > ${fl/.txt/_points.txt}
	printf "\n"
done
