#!/bin/bash

pop=$1

REF=../../reference_genome_ensembl
PICARD=/home/fb4/palma-vera/FBN_HOME/Tools/picard_2.18.11
bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
GATK=~/FBN_HOME/Tools/gatk-4.0.6.0

in_vcf=../output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf
pop_vcf=../output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.allrecords.${pop}.vcf
samps=../../sample_info/vcf_samples_$pop

echo "# subset main vcf for population (exclude ref alleles)"
$bcftools/bcftools view --samples-file $samps $in_vcf -o $pop_vcf
printf "\n\n"

echo "# Number of SNPs in pop:"
grep -v '^#' $pop_vcf | wc -l
printf "\n\n"
