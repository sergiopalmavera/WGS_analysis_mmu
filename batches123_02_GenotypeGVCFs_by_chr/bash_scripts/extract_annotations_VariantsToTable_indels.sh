#!/bin/bash

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0

vcf=cohort_biallelicINDELs.vcf.gz

# Only biallelic SNPs without missingness
$GATK/gatk IndexFeatureFile -I ../output/$vcf

$GATK/gatk VariantsToTable \
	-V ../output/$vcf \
	-F CHROM -F POS -F QD -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F SOR \
	-O ../output/${vcf/.vcf/.table}
