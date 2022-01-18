#!/bin/bash

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9
vcf=../output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf
# https://darencard.net/blog/2017-01-13-missing-data-proportions-vcf/
paste \
<(bcftools query -f '[%SAMPLE\t]\n' $vcf | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' $vcf | awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2)
