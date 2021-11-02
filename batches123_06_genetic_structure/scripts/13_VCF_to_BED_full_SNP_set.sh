#!/bin/bash

fl_dir=../../batches123_04_final_VCF/output

plink --vcf $fl_dir/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf --make-bed --chr 1-19 --out ../data/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered


