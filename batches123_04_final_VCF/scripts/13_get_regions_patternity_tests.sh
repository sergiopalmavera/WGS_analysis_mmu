#!/bin/bash

# Extract regions of interest for patternity test

GATK=~/FBN_HOME/Tools/gatk-4.0.6.0
vcftools=/home/fb4/palma-vera/FBN_HOME/Tools/vcftools_0.1.13/cpp
REF=../../reference_genome_ensembl

for vcf in ../output/*.allrecords.*.vcf
do
	echo "# $vcf"

	echo "## Indexing"
	$GATK/gatk IndexFeatureFile -F $vcf
	
	echo "## Creating vcf subsets"
	$GATK/gatk SelectVariants \
		-R $REF/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		-V $vcf \
		-L 9:86562500-86565000 \
		-L 17:40906500-40908000 \
		-L 17:41101000-41103500 \
		-O ${vcf/.vcf/_paternity.vcf} 
	printf "\n"

	echo "## Get alles in regions of interest"
	$vcftools/vcftools --freq --vcf  ${vcf/.vcf/_paternity.vcf} --out ${vcf/.vcf/_paternity}
	printf "\n"
done

