#!/bin/bash

REF=../../reference_genome_ensembl
PICARD=/home/fb4/palma-vera/FBN_HOME/Tools/picard_2.18.11
GATK=~/FBN_HOME/Tools/gatk-4.0.6.0

#echo "# Metrics for raw vcf"
#java -jar $PICARD/picard.jar CollectVariantCallingMetrics \
#	INPUT=../../batches123_02_GenotypeGVCFs_by_chr/output/cohort.vcf.gz \
#	OUTPUT=../metrics/cohort.table \
#	DBSNP=$REF/mus_musculus.vcf
#printf "\n"

#echo "# Metrics for vcf biSNPs, after VQSR (incl PASS and not PASS)"
#java -jar $PICARD/picard.jar CollectVariantCallingMetrics \
#	INPUT=../../batches123_03_VariantQualityScoreRecalibration/output/cohort_biallelicSNPs_VQSR95.vcf \
#	OUTPUT=../metrics/cohort_biallelicSNPs_VQSR95.table \
#	DBSNP=$REF/mus_musculus.vcf
#printf "\n"

vcf_fl=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf

$GATK/gatk IndexFeatureFile -F ../output/$vcf_fl

echo "# Metrics for vcf biSNPs, PASS, filtered by GQ and DP"
java -jar $PICARD/picard.jar CollectVariantCallingMetrics \
	INPUT=../output/$vcf_fl \
	OUTPUT=../output/${vcf_fl/.vcf/} \
	DBSNP=$REF/mus_musculus.vcf
