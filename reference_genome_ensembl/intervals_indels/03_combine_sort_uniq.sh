#!/bin/bash

#bedtools=/home/fb4/palma-vera/FBN_HOME/Tools/bedtools_version_2.29.2
bedops=~/FBN_HOME/Tools/bedops/bin/

echo "# Concatenate bed files ..."
cat ./mus_musculus_dels.bed ./mus_musculus_ins.bed > ./tmp.bed
printf "\n"

echo "# Make it lighter"
awk '{print $1 "\t" $2 "\t" $3 "\t"}' ./tmp.bed > ./tmp2.bed

echo "# Sort bed file ..."
$bedops/sort-bed-typical ./tmp2.bed > ./intervals_indels_sorted.bed
printf "\n"

echo "# Uniq ..."
uniq ./intervals_indels_sorted.bed > ./intervals_indels_sorted_uniq.bed


rm ./tmp.bed ./tmp2.bed

echo "# Done"

