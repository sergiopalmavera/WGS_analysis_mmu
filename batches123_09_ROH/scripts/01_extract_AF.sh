#!/bin/bash

bcftools=~/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin

vcf_path=../../batches123_04_final_VCF/output

vcf_fl=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf

#$bcftools/bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%AF\n' $vcf_path/$vcf_fl | bgzip -c > ../output/${vcf_fl/.vcf/.AF.tab.gz} && tabix -s1 -b2 -e2 ../output/${vcf_fl/.vcf/.AF.tab.gz}
$bcftools/bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%AF\n' $vcf_path/$vcf_fl | bgzip -c > ../output/${vcf_fl/.vcf/.AF.tab.gz} && tabix -s1 -b2 -e2 ../output/${vcf_fl/.vcf/.AF.tab.gz}
