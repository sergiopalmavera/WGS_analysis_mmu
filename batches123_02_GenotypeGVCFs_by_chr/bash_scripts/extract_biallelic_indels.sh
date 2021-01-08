#!/bin/bash

GATK=/home/fb4/palma-vera/FBN_HOME/Tools/gatk-4.1.5.0
in_vcf=cohort.vcf.gz
out_vcf=cohort_biallelicINDELs.vcf.gz

$GATK/gatk SelectVariants \
	-R ../../reference_genome_ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V ../output/$in_vcf \
	--select-type-to-include INDEL \
	--restrict-alleles-to BIALLELIC \
	-O ../output/$out_vcf
