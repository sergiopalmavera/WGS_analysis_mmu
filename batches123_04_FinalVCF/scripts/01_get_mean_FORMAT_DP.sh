#!/bin/bash

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

VAR_DIR=../../batches123_03_VariantQualityScoreRecalibration/output/

$bcftools/bcftools query -f '[%DP ]' $VAR_DIR/cohort_biallelicSNPs_VQSR90_PASS.vcf -o ../output/DP.tmp

#awk '{s+=$1}END{print "ave:",s/NR}' ../output/DP.tmp
#awk '{sum+=$1; sumsq+=$1*$1}END{print "sd", sqrt(sumsq/NR - (sum/NR)**2)}' ../output/DP.tmp
