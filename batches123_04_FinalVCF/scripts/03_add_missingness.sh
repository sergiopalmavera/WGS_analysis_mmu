#!/bin/bash

# this scripts changes genotypes to missing when a sample GQ is less than 20, DP less than 4 or DP > meanDP+4*sdDP

pct=95

vcftools=~/FBN_HOME/Tools/vcftools_0.1.13/cpp
vars=../../batches123_03_VariantQualityScoreRecalibration/output
vcf_in=cohort_biallelicSNPs_VQSR${pct}_PASS.vcf
vcf_out=${vcf_in/.vcf/_AddedMissingness}

$vcftools/vcftools --vcf $vars/$vcf_in --minGQ 20 --minDP 4 --maxDP 77 --recode --recode-INFO-all --out ../output/$vcf_out 
