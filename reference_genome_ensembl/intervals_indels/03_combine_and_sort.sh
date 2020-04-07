#!/bin/bash

bedtools=/home/fb4/palma-vera/FBN_HOME/Tools/bedtools_version_2.29.2

echo "# Concatenate bed files ..."
cat ./mus_musculus_dels.bed ./mus_musculus_ins.bed > ./tmp.bed
printf "\n"

echo "# Sort bed file ..."
$bedtools/bedtools.static.binary sort -i ./tmp.bed > ./mus_musculus_indels_sorted.bed
printf "\n"

rm ./tmp.bed

echo "# Done"

