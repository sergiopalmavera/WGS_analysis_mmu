#!/bin/bash

# Print hostname and data for reference
hostname
date
chr=8

echo "# Running Intervals for Chromosome $chr"

if [ -e ../output/chr$chr ] 
then 
	rm -r ../output/chr$chr; mkdir ../output/chr$chr ../output/chr$chr/DB
else
	mkdir ../output/chr$chr ../output/chr$chr/DB
fi

dir_chr=../output/chr$chr 
dir_chr_db=../output/chr$chr/DB

BED=../../reference_genome_ensembl/unmasked_intervals # path to unmasked intervals

# Get intervals
awk -v x=$chr '$1 == x {print}' $BED/intervals_unmasked_windowed_5Mb.bed > $dir_chr/chr${chr}.bed

echo "## Number of intervals (also parallel jobs)"
wc -l $dir_chr/chr${chr}.bed

while read l
do
	echo "## Working on interval $l"
	# Prepare subdirectories and intervals
	nm=$(echo $l | awk '{print $1 "_" $2 "_" $3}') #prepare name for subdir and sub-bed file
	echo $l > $dir_chr/${nm}.bed # Export interval subset
	
	nohup ./ConsolidateGVCFs_template.sh $dir_chr/${nm}.bed $dir_chr_db/$nm &> $dir_chr/ConsolidateGVCFs_chr${nm}.out & # submit job
done < $dir_chr/chr${chr}.bed

echo "## All jobs for chromosome submitted. One job per interval"
