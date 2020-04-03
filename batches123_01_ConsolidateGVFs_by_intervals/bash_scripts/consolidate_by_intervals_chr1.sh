#!/bin/bash

# Print hostname and data for reference
hostname
date

chr=1

echo "# Running Intervals for Chromosome $chr"

mkdir ../output/chr$chr #Make a chromosome file
mkdir ./chr${chr}_logs # Make a subdir to store log files
mkdir ./intervals/chr$chr # Make a subdirectory for intervals
dir_iv=./intervals/chr$chr # store directory for intervals for simplicity

BED=../../reference_genome_ensembl/unmasked_intervals # path to unmasked intervals
BEDTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/bedtools_version_2.29.2 #absolute path, modify accordingly

# Get intervals equal or larger than 5Mb
awk -v x=$chr '$1 == x {print}' $BED/intervals_unmasked.bed | awk '{print $0,  $3 - $2}' | awk -F' ' '$4 >= 5000000 {print $1 "\t" $2 "\t" $3}' > $dir_iv/chr${chr}_large_i.bed

# Export small intervals
awk -v x=$chr '$1 == x {print}' $BED/intervals_unmasked.bed | awk '{print $0,  $3 - $2}' | awk -F' ' '$4 < 5000000 {print $1 "\t" $2 "\t" $3}' > $dir_iv/chr${chr}_small_i.bed

# Split large intervals into intervals of max 5Mb
$BEDTOOLS/bedtools.static.binary makewindows -b $dir_iv/chr${chr}_large_i.bed -w 5000000 > $dir_iv/chr${chr}_large_i_split.bed

# Concatenate small intervals with large intervals split
cat $dir_iv/chr${chr}_small_i.bed $dir_iv/chr${chr}_large_i_split.bed > $dir_iv/intervals_chr${chr}.bed

echo "## Total nr of intervals:"
wc -l $dir_iv/intervals_chr${chr}.bed

while read l # loop over each interval for chromosome
do
	echo "## Working on interval $l"
	# Prepare subdirectories and intervals
	nm=$(echo $l | awk '{print $1 "_" $2 "_" $3}') #prepare name for subdir and sub-bed file
	echo $l > $dir_iv/${nm}.bed # Export interval subset

	echo "## ConsolidateGVCFs"	
	nohup ./ConsolidateGVCFs_template.sh $dir_iv/${nm}.bed ../output/chr$chr/$nm &> ./chr${chr}_logs/ConsolidateGVCFs_chr${nm}.out & # submit job
	printf "\n"
done < $dir_iv/intervals_chr${chr}.bed

echo "## All jobs for chromosome submitted. One job per interval"


