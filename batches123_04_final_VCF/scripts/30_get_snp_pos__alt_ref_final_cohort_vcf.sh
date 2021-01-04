#!/bin/bash

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

vcf=../output/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.vcf

echo "# get INDELs for population"
$bcftools/bcftools query -f '%CHROM %POS %REF %ALT\n' $vcf > ${vcf/.vcf/.INDELids2}
printf "\n\n"

echo "# Number of INDEL ids"
wc -l ${vcf/.vcf/.INDELids2}

