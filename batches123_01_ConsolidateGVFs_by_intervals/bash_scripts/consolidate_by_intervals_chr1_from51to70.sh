#!/bin/bash

# Print hostname and data for reference
hostname
date
chr=1
from=51
to=70

echo "# Running Intervals for Chromosome $chr"

if [ ! -e ../output/chr$chr ] 
then 
	mkdir ../output/chr$chr ../output/chr$chr/DB ../output/chr$chr/intervals ../output/chr$chr/logs
fi


dir_intervals=../output/chr$chr/intervals
dir_logs=../output/chr$chr/logs
dir_db=../output/chr$chr/DB

bed_ss=chr${chr}_from${from}to${to}.bed

sed -n $from,${to}p $dir_intervals/chr${chr}.bed > $dir_intervals/$bed_ss

echo "## Number of intervals (also parallel jobs)"
wc -l $dir_intervals/$bed_ss

while read l
do
	echo "## Working on interval $l"
	# Prepare subdirectories and intervals
	nm=$(echo $l | awk '{print $1 "_" $2 "_" $3}') #prepare name for subdir and sub-bed file
	echo $l > $dir_intervals/${nm}.bed # Export interval subset
	
	nohup ./ConsolidateGVCFs_template.sh $dir_intervals/${nm}.bed $dir_db/$nm &> $dir_logs/ConsolidateGVCFs_chr${nm}.out & # submit job
	
	# Print job information oncl pid
	echo "## Job information:"
	ps -xjf | head -1
	ps -xjf | grep "$nm" | grep "bash" | head -1

done < $dir_intervals/$bed_ss

echo "## All jobs for chromosome submitted. One job per interval"
