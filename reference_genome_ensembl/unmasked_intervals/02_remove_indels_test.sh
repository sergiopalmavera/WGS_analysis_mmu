#!/bin/bash

BEDTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/bedtools_version_2.29.2

head ./intervals_unmasked.bed > ./TMP/A.bed
head ../intervals_indels/intervals_indels_sorted_uniq.bed > ./TMP/B.bed

cat ./TMP/A.bed
printf "\n"

cat ./TMP/B.bed
printf "\n"

$BEDTOOLS/bedtools.static.binary subtract -a ./TMP/A.bed -b ./TMP/B.bed 
