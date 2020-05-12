#!/bin/bash

pct=95

GATK=~/FBN_HOME/Tools/gatk-4.0.6.0

$GATK/gatk VariantsToTable \
	-V ../output/cohort_biallelicSNPs_VQSR${pct}_PASS_AddedMissingness.recode.vcf \
	-O ../output/cohort_biallelicSNPs_VQSR${pct}_PASS_AddedMissingness.recode.table \
	-F CHROM -F POS -GF GT

