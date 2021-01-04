#!/bin/bash

VCFTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp

$VCFTOOLS/vcftools --hist-indel-len --gzvcf ../output/cohort_biallelicINDELs.vcf.gz --out ../output/cohort_biallelicINDELs
