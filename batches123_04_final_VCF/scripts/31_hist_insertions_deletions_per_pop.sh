#!/bin/bash

pop=$1

pop_vcf=cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.${pop}.vcf

VCFTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp

$VCFTOOLS/vcftools --hist-indel-len --vcf ../output/$pop_vcf --out ../output/${pop_vcf/.vcf/}

