#!/bin/bash

pop=$1

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

in_vcf=../output/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.vcf
pop_vcf=../output/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.allrecords.${pop}.vcf
samps=../../sample_info/vcf_samples_$pop

echo "# subset main vcf for population (include fixed ref alleles)"
$bcftools/bcftools view --samples-file $samps $in_vcf -o $pop_vcf
printf "\n\n"

echo "# Number of SNPs in pop:"
grep -v '^#' $pop_vcf | wc -l
printf "\n\n"
