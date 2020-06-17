#!/bin/bash

export BCFTOOLS_PLUGINS=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/libexec/bcftools

~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin/bcftools +prune -l 0.2 -w 50kb ../../batches123_04_final_VCF/output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf -o ../data/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ldpruned.vcf 
