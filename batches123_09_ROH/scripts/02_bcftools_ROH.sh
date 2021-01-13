#!/bin/bash

from_sample=$1
to_sample=$2

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
vcf_path=../../batches123_04_final_VCF/output
vcf_fl=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf
out_fl=samp${from_sample}_to_samp${to_sample}.roh

sed -n "$from_sample,${to_sample}p" ../../sample_info/vcf_samples.txt > ./TMP/samples_${from_sample}_${to_sample}.txt

echo "# Processing samples:"
cat ./TMP/samples_${from_sample}_${to_sample}.txt
printf "\n"
date
printf "\n"
time $bcftools/bcftools roh --AF-file ../output/${vcf_fl/.vcf/.AF.tab.gz} $vcf_path/$vcf_fl -o ../output/$out_fl -S ./TMP/samples_${from_sample}_${to_sample}.txt
