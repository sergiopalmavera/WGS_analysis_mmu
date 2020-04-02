#!/bin/bash

# Print hostname and data for reference
hostname
date

chr=1

echo "# Running Intervals for Chromosome $chr"

mkdir ../output/chr$chr #Make a chromosome file
mkdir ./chr${chr}_logs # Make a subdir to store log files
mkdir ./intervals/$chr # Make a subdirectory for intervals

BED=/projekte/I2-SOS-FERT/GitHub/WGS_analysis_mmu/reference_genome_ensembl/unmasked_intervals # path to unmasked intervals

awk -v x=$chr '$1 == x {print}' $BED/intervals_unmasked.bed  > tmp_${chr}.bed # extract chromosome intervals

echo "## Total nr of intervals:"
wc -l tmp_${chr}.bed

while read l # loop over each interval for chromosome
do
	echo "## Working on interval $l"
	# Prepare subdirectories and intervals
	nm=$(echo $l | awk '{print $1 "_" $2 "_" $3}') #prepare name for subdir and sub-bed file
	echo $l > ./intervals/$chr/${nm}.bed # Export interval subset

	echo "## ConsolidateGVCFs"	
	nohup ./ConsolidateGVCFs_template.sh ./intervals/$chr/${nm}.bed ../output/chr$chr/$nm &> ./chr${chr}_logs/ConsolidateGVCFs_chr${nm}.out & # submit job
	printf "\n"
done < tmp_${chr}.bed

rm tmp_${chr}.bed

echo "## All jobs for chromosome submitted. One job per interval"
