#!/bin/bash

GATK=~/FBN_HOME/Tools/gatk-4.0.6.0

vcf=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf

$GATK/gatk VariantsToTable \
	-V ../output/$vcf \
	-O ../output/${vcf/.vcf/_SNPs_in_RDD_genes_GT.table} \
	-L ../../00_dashboard/data/RDD_genes/snps_in_RDD_genes.intervals \
	-F CHROM -F POS -GF GT



