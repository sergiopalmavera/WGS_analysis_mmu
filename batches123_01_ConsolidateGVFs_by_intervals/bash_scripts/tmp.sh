#!/bin/bash

# Print hostname and data for reference
hostname
date
chr=$1
max_size=5000000
max_nr_jobs=50

echo "# Running Intervals for Chromosome $chr"

dir_iv=./TMP # store directory for intervals for simplicity

BED=../../reference_genome_ensembl/unmasked_intervals # path to unmasked intervals
BED_INDELS=../../reference_genome_ensembl/intervals_indels
BEDTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/bedtools_version_2.29.2 #absolute path, modify accordingly

# Get intervals equal or larger than max_size
awk -v x=$chr '$1 == x {print}' $BED/intervals_unmasked_excl_indels.bed | awk '{print $0,  $3 - $2}' | awk -v ms=$max_size -F' ' '$4 >= ms {print $1 "\t" $2 "\t" $3}' > $dir_iv/chr${chr}_large_i.bed

# Export small intervals
awk -v x=$chr '$1 == x {print}' $BED/intervals_unmasked_excl_indels.bed | awk '{print $0,  $3 - $2}' | awk -v ms=$max_size -F' ' '$4 < ms {print $1 "\t" $2 "\t" $3}' > $dir_iv/chr${chr}_small_i.bed

# Split large intervals into intervals of max_size
$BEDTOOLS/bedtools.static.binary makewindows -b $dir_iv/chr${chr}_large_i.bed -w $max_size > $dir_iv/chr${chr}_large_i_split.bed

# Concatenate small intervals with large intervals split
cat $dir_iv/chr${chr}_small_i.bed $dir_iv/chr${chr}_large_i_split.bed > $dir_iv/intervals_chr${chr}.bed

n_intervals=$(wc -l $dir_iv/intervals_chr${chr}.bed | awk '{print $1}')
echo "## Total nr of intervals excl indels: $n_intervals"

# Split intervals into list of max_nr_jobs
n_lines=$(( (n_intervals / max_nr_jobs) + 1 ))
echo "## Number of intervals per job: $n_lines (total jobs running in parallel = $max_nr_jobs)"

split -l $n_lines $dir_iv/intervals_chr${chr}.bed $dir_iv/intervals_chr${chr}_chunk --additional-suffix .bed --numeric-suffixes=1

# Check integrity of files
echo "## Sanity check: Number of intervals across jobs (should be same as nr of intervals = $n_intervals)"
cat $dir_iv/*_chunk* | wc -l



