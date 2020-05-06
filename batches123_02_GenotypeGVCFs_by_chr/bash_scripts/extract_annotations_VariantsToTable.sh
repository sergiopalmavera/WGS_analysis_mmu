#!/bin/bash

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0

# All biallelic SNPs
#$GATK/gatk VariantsToTable \
#	-V ../output/cohort_biallelicSNPs.vcf.gz \
#	-F CHROM -F POS -F FILTER -F QD -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F SOR \
#	-O ../output/cohort_biallelicSNPs.table

# Only biallelic SNPs without missingness
$GATK/gatk IndexFeatureFile -I ../output/cohort_biallelicSNPs_noMissing.vcf.gz

$GATK/gatk VariantsToTable \
	-V ../output/cohort_biallelicSNPs_noMissing.vcf.gz \
	-F CHROM -F POS -F FILTER -F QD -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F SOR \
	-O ../output/cohort_biallelicSNPs_noMissing.table
