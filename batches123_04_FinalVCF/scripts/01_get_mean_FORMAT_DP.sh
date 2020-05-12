#!/bin/bash

pct=99.9

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

VAR_DIR=../../batches123_03_VariantQualityScoreRecalibration/output/

$bcftools/bcftools query -f '[%DP ]\n' $VAR_DIR/cohort_biallelicSNPs_VQSR${pct}_PASS.vcf -o ../output/DP.tmp
