#!/bin/bash

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

vcf=../output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf

echo "# get SNPs for population"
$bcftools/bcftools query -f '%CHROM  %POS  %REF  %ALT\n' $vcf > ${vcf/.vcf/.SNPids2}
printf "\n\n"

echo "# Number of SNP ids"
wc -l ${vcf/.vcf/.SNPids2}

