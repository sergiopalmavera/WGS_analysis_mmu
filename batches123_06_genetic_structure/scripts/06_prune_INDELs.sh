#!/bin/bash

export BCFTOOLS_PLUGINS=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/libexec/bcftools

in_vcf=cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.vcf
out_vcf=cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.ldpruned.vcf

~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin/bcftools +prune -l 0.2 -w 50kb ../../batches123_04_final_VCF/output/$in_vcf -o ../data/$out_vcf 
