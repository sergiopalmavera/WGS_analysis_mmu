#!/bin/bash

chr=$1

BCFTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
OUT=../output

time $BCFTOOLS/bcftools concat $OUT/chr$1/*.vcf.gz -o $OUT/concat_chr${chr}.vcf.gz -Oz
