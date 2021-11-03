#!/bin/bash

bedops=~/FBN_HOME/Tools/bedops/bin/

echo "# Concatenate bed files ..."
cat ../unmasked_intervals/intervals_unmasked_minus_indels.bed ../intervals_indels/intervals_indels_sorted_uniq.bed > ./tmp.bed
printf "\n"

echo "# Sort bed file ..."
$bedops/sort-bed-typical ./tmp.bed > ./intervals_sorted.bed
printf "\n"

echo "# Anyl duplications?"
wc -l ./intervals_sorted.bed
uniq ./intervals_sorted.bed | wc -l


rm ./tmp.bed
echo "# Done"

