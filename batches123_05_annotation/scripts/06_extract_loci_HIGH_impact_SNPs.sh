#!/bin/bash

# Prepare intervlas to extract sites from final VCF with GATK's SelectVariants

awk '{print $1 ":" $2 "-" $2}' ../output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ann.HIGH.tab > ../output/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ann.HIGH.intervals
