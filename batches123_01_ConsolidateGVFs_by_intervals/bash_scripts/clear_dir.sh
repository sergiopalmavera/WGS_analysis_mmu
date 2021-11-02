#!/bin/bash

# Print hostname and data for reference
hostname
date
chr=$1
from=$2
to=$3

echo "# Running Intervals for Chromosome $chr"

dir_intervals=../output/chr$chr/intervals
dir_logs=../output/chr$chr/logs
dir_db=../output/chr$chr/DB

# Define name of bed file with intervals for this subset
bed_ss=chr${chr}_from${from}to${to}.bed

# Remove all interval files
while read l
do
	echo "## Working on interval $l"
	# Prepare subdirectories and intervals
	nm=$(echo $l | awk '{print $1 "_" $2 "_" $3}') #prepare name for subdir and sub-bed file
	
	if [ -e $dir_db/$nm ]; then rm -r $dir_db/$nm; fi
	if [ -e $dir_intervals/${nm}.bed ]; then rm $dir_intervals/${nm}.bed; fi
	if [ -e $dir_logs/ConsolidateGVCFs_chr${nm}.out ]; then rm $dir_logs/ConsolidateGVCFs_chr${nm}.out; fi

done < $dir_intervals/$bed_ss

# Remove interval list
rm -r $dir_intervals/$bed_ss
