#!/bin/bash

echo "Final SNP VCF"
md5sum ../output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf

echo "Final INDEL VCF"
md5sum ../output/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.filtered.vcf
