#!/bin/bash

# Print hostname and data for reference
hostname
date
chr=2
from=41
to=59

echo "# Running Intervals for Chromosome $chr"

if [ ! -e ../output/chr$chr ] 
then 
	mkdir ../output/chr$chr ../output/chr$chr/DB ../output/chr$chr/intervals ../output/chr$chr/logs
fi

dir_intervals=../output/chr$chr/intervals
dir_logs=../output/chr$chr/logs
dir_db=../output/chr$chr/DB

BED=../../reference_genome_ensembl/intervals_unmasked_bound_by_indels # path to unmasked intervals

# Get intervals
# Do no create if already created for this chromosome
if [ ! -e $dir_intervals/chr${chr}.bed ]
then
	awk -v x=$chr '$1 == x {print}' $BED/intervals_bounded_by_indels_max5Mb.bed > $dir_intervals/chr${chr}.bed
fi

# Define name of bed file with intervals for this subset
bed_ss=chr${chr}_from${from}to${to}.bed

# Make subset of intervals for this script
# Do nothing if already exists (job already submitted probably) 
if [ -e $dir_intervals/$bed_ss ]
then
	printf "\n"
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "The bed file $dir_intervals/$bed_ss exists ..."	
	echo "This script was probably executed before ..."
	echo "Nothing will be done."
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	printf "\n"
else
	sed -n $from,${to}p $dir_intervals/chr${chr}.bed > $dir_intervals/$bed_ss #make intervals for this subset

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
fi

