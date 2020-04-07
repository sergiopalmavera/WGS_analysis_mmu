#!/bin/bash

awk '{print $2 "\t" $3 "\t" $3-$2}' mus_musculus_indels_sorted.bed | awk '{print $3}' | sort -r | uniq | head -10
