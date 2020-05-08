#!/bin/bash

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.0.6.0
REF=../../reference_genome_ensembl

$GATK/gatk SelectVariants \
	-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V ../output/cohort_biallelicSNPs_VQSR90.vcf \
	-O ../output/cohort_biallelicSNPs_VQSR90_PASS.vcf \
	--exclude-filtered
