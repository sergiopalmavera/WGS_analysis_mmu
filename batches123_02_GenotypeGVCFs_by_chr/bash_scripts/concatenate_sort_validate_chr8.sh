#!/bin/bash

chr=8

BCFTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/samtools_1.9/bcftools-1.9_installed/bin
GATK=~/FBN_HOME/Tools/gatk-4.1.5.0

out_nm=cohort_chr${chr}.vcf.gz

echo "# concatenate and sort ..."
time $BCFTOOLS/bcftools concat ../output/chr$chr/*.vcf.gz -Ou | $BCFTOOLS/bcftools sort -Oz -o ../output/$out_nm
printf "\n"

echo "# Indexing ..."
$GATK/gatk IndexFeatureFile -I ../output/$out_nm
printff "\n"

echo "# Validate vcf ..."
$GATK/gatk ValidateVariants \
	-R ../../reference_genome_ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa \
	-V ../output/$out_nm \
	--validation-type-to-exclude ALL
