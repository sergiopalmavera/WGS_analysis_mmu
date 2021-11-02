#!/bin/bash

# Print hostname and data for reference
hostname
date
chr=X

echo "# Running Intervals for Chromosome $chr"

dir_intervals=../output/chr$chr/intervals
dir_logs=../output/chr$chr/logs
dir_db=../output/chr$chr/DB

sed -n 81,125p $dir_intervals/chr${chr}.bed > $dir_intervals/chr${chr}_pt3.bed

echo "## Number of intervals (also parallel jobs)"
wc -l $dir_intervals/chr${chr}_pt3.bed

while read l
do
	echo "## Working on interval $l"
	# Prepare subdirectories and intervals
	nm=$(echo $l | awk '{print $1 "_" $2 "_" $3}') #prepare name for subdir and sub-bed file
	echo $l > $dir_intervals/${nm}.bed # Export interval subset
	
	nohup ./ConsolidateGVCFs_template.sh $dir_intervals/${nm}.bed $dir_db/$nm &> $dir_logs/ConsolidateGVCFs_chr${nm}.out & # submit job
done < $dir_intervals/chr${chr}_pt3.bed

echo "## All jobs for chromosome submitted. One job per interval"
