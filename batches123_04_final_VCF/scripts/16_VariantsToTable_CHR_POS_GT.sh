#!/bin/bash

GATK=~/FBN_HOME/Tools/gatk-4.0.6.0

vcf=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf

$GATK/gatk VariantsToTable \
	-V ../output/${vcf/.vcf/.HIGHimpact.vcf} \
	-O ../output/${vcf/.vcf/.HIGHimpact.table} \
	-F CHROM -F POS -GF GT



