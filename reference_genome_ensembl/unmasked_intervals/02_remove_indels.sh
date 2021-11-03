#!/bin/bash

BEDTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/bedtools_version_2.29.2

A=./intervals_unmasked.bed
B=../intervals_indels/intervals_indels_sorted_uniq.bed

$BEDTOOLS/bedtools.static.binary subtract -a $A -b $B > ./intervals_unmasked_minus_indels.bed
