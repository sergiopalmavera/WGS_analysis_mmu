#!/bin/bash

# This file renames csv, txt and html files that were not appended with the cohort vcf file.

mv snpEff_summary.csv ../output/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.ann.csv
mv snpEff_summary.genes.txt ../output/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.ann.genes.txt
mv snpEff_summary.html ../output/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.ann.html
