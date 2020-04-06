#!/bin/bash


chr=$1

echo "# Running Intervals for Chromosome $chr"

BED=../../reference_genome_ensembl/unmasked_intervals # path to unmasked intervals
BEDTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/bedtools_version_2.29.2 #absolute path, modify accordingly

# Get intervals equal or larger than 5Mb
awk -v x=$chr '$1 == x {print}' $BED/intervals_unmasked.bed | awk '{print $0,  $3 - $2}' | awk -F' ' '$4 >= 5000000 {print $1 "\t" $2 "\t" $3}' > ./TMP/large_i.bed
echo "## Number of intervals >= 5Mb"
wc -l ./TMP/large_i.bed 

# Export small intervals
awk -v x=$chr '$1 == x {print}' $BED/intervals_unmasked.bed | awk '{print $0,  $3 - $2}' | awk -F' ' '$4 < 5000000 {print $1 "\t" $2 "\t" $3}' > ./TMP/small_i.bed
echo "## Number of intervals < 5Mb"
wc -l ./TMP/small_i.bed

# Split large intervals into intervals of max 5Mb
$BEDTOOLS/bedtools.static.binary makewindows -b ./TMP/large_i.bed -w 5000000 > ./TMP/large_i_split.bed
echo "## Number of large intervals split by max 5Mb"
wc -l ./TMP/large_i_split.bed

# Concatenate small intervals with large intervals split
cat ./TMP/small_i.bed ./TMP/large_i_split.bed > ./TMP/intervals.bed

echo "## Total nr of intervals:"
wc -l ./TMP/intervals.bed



