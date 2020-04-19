#!/bin/bash

# Remove direcotries
while read l
do
	nm=$(echo $l | awk '{print $1 "_" $2 "_" $3}')
	rm -r ../output/chrX/DB/$nm

	# Remove bed files
	rm ../output/chrX/intervals/${nm}.bed

	# Remove logs
	rm ../output/chrX/logs/ConsolidateGVCFs_chr${nm}.out
done < ../output/chrX/intervals/chrX_pt2.bed

