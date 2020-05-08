#!/bin/bash

# this scripts changes genotypes to missing when a sample GQ is less than 20, DP less than 4 or DP > meanDP+4*sdDP

vcftools=~/FBN_HOME/Tools/vcftools_0.1.13/cpp
vars=../../batches123_03_VariantQualityScoreRecalibration/output
vcf_in=cohort_biallelicSNPs_VQSR90_PASS.vcf
vcf_out=${vcf_in/.vcf/_AddedMissingness}

$vcftools/vcftools --vcf $vars/$vcf_in --minGQ 20 --minDP 4 --maxDP 70 --out ../output/$vcf_out --recode-INFO-all
