#!/bin/bash

GATK=~/FBN_HOME/Tools/gatk-4.0.6.0

$GATK/gatk VariantsToTable \
	-V ../output/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.vcf \
	-O ../output/cohort_biallelicINDELs_VQSR99_PASS_withmissingness.table \
	-F CHROM -F POS -GF GT

