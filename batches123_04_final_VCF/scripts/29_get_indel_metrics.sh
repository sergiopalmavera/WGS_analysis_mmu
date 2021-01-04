#!/bin/bash

VCFTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp

vcf=cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.vcf

$VCFTOOLS/vcftools --hist-indel-len --vcf ../output/$vcf --out ../output/${vcf/.vcf/}

