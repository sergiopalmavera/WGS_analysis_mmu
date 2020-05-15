#!/bin/bash

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=../../reference_genome_ensembl

$GATK/gatk IndexFeatureFile -F ../output/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.vcf

$GATK/gatk SelectVariants \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V ../output/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.vcf \
	-L ../output/keep_snps_NminPerGroup20.intervals \
	-O ../output/cohort_biallelicSNPs_VQSR95_PASS_AddedMissingness.recode.filtered.vcf 