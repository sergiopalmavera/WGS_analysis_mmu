#!/bin/bash

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

$bcftools/bcftools view ../output/cohort_biallelicSNPs.vcf.gz -i 'F_PASS(GT!="mis") = 1' -Oz -o ../output/cohort_biallelicSNPs_noMissing.vcf.gz
