#!/bin/bash

# Extract sites labeled as "HIGH" impact by SNPeff

GATK=~/FBN_HOME/Tools/gatk-4.0.6.0

REF=../../reference_genome_ensembl

intervals_path=../../batches123_05_annotation/output

vcf=cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.vcf

$GATK/gatk SelectVariants \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V ../output/$vcf \
	-L $intervals_path/cohort_biallelicSNPs_VQSR95_PASS_withmissingness.filtered.ann.HIGH.intervals \
	-O ../output/${vcf/.vcf/.HIGHimpact.vcf} 


