#!/bin/bash

fl_dir=../data
plink --vcf $fl_dir/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ldpruned.vcf --make-bed --chr 1-19 --out $fl_dir/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ldpruned


