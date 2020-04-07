#!/bin/bash

BEDTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/bedtools_version_2.29.2

A=./intervals_unmasked.bed
B=../intervals_indels/mus_musculus_indels_sorted.bed

$BEDTOOLS/bedtools.static.binary subtract -a $A -b $B > ./intervals_unmasked_excl_indels.bed
