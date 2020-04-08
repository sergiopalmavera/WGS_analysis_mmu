#!/bin/bash

awk '{print $1 "\t" $2 "\t" $3}' mus_musculus_indels_sorted.bed > mus_musculus_indels_sorted_lite.bed
